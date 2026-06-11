#!/usr/bin/env python3
"""Threshold-sweep HGRS metric with ontology-stratified micro/macro summaries."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
from pathlib import Path
from typing import Any

from fpr_metric_common import (
    BINDING_TERMS,
    GO_RE,
    AnnotationMap,
    OntologyMap,
    ParentMap,
    evaluate_fast_mapping,
    parse_do_not_annotate,
    parse_existing_annotations,
    parse_go_parent_aspects,
    parse_go_parents,
    parse_reference_annotations,
    remove_redundancy_fast,
    set_metric_logging,
)


ONTOLOGIES = ("biological_process", "cellular_component", "molecular_function")
LOGGER = logging.getLogger("fpr_metric_sweep")
OUTPUT_FIELDS = [
    "filename",
    "ns",
    "tau",
    "pr_micro",
    "rc_micro",
    "f_micro",
    "pr_macro",
    "rc_macro",
    "f_macro",
    "E",
    "L",
    "M",
    "U",
    "Z",
    "n_predicted_proteins",
    "n_reference_proteins",
    "n_mapped_prediction_proteins",
]

ScoredPredictionMap = dict[str, dict[str, float]]
ScoredPredictionByOntology = dict[str, ScoredPredictionMap]

_WORKER_CONTEXT: dict[str, Any] = {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep score thresholds for the hierarchy-aware GO recovery metric "
            "and report ontology-stratified micro/macro summaries."
        )
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-p", "--prediction", help="single scored prediction TSV/CSV")
    input_group.add_argument("-d", "--prediction-dir", help="directory of scored prediction files; searched recursively")
    parser.add_argument("-t", "--existing", required=True, help="existing annotation TSV")
    parser.add_argument("-r", "--reference", required=True, help="new/reference annotation TSV")
    parser.add_argument("-n", "--do-not-annotate", required=True, help="GO do-not-annotate TSV")
    parser.add_argument("-g", "--go-parent", required=True, help="goparents lookup file")
    parser.add_argument("-o", "--out-file", required=True, help="threshold sweep summary TSV")
    parser.add_argument("--best-macro-output", help="best-threshold summary TSV selected by macro_f1")
    parser.add_argument("--best-micro-output", help="best-threshold summary TSV selected by micro_f1")
    parser.add_argument("--tau-step", type=Decimal, required=True, help="threshold step, for example 0.01")
    parser.add_argument("--tau-min", type=Decimal, default=Decimal("0"), help="minimum threshold, default 0")
    parser.add_argument("--tau-max", type=Decimal, default=Decimal("1"), help="maximum threshold, default 1")
    parser.add_argument("--predicted-format", choices=["auto", "tsv", "csv"], default="auto")
    parser.add_argument("--glob", default="*", help="file glob for recursive directory mode")
    parser.add_argument("--protein-col", type=int, default=0, help="0-based protein column in prediction file")
    parser.add_argument("--go-col", type=int, default=1, help="0-based GO column in prediction file")
    parser.add_argument("--score-col", type=int, default=2, help="0-based score column in prediction file")
    parser.add_argument(
        "--ontology",
        action="append",
        choices=ONTOLOGIES,
        help="ontology to evaluate; repeatable. Defaults to all three GO ontologies.",
    )
    parser.add_argument("--workers", type=int, default=1, help="process workers for tau/ontology jobs")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    return parser.parse_args()


def delimiter_from_format(path: Path, input_format: str) -> str:
    if input_format == "csv":
        return ","
    if input_format == "tsv":
        return "\t"
    if path.suffix.lower() == ".csv":
        return ","
    if path.suffix.lower() in {".tsv", ".tab"}:
        return "\t"
    with path.open(newline="") as handle:
        sample = handle.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error:
        return "\t"


def looks_like_header(row: list[str], protein_col: int, go_col: int, score_col: int) -> bool:
    if len(row) <= max(protein_col, go_col, score_col):
        return True
    if not GO_RE.fullmatch(row[go_col].strip()):
        return True
    try:
        float(row[score_col])
    except ValueError:
        return True
    return False


def parse_scored_predictions(
    path: str | Path,
    *,
    input_format: str,
    protein_col: int,
    go_col: int,
    score_col: int,
    existing: AnnotationMap,
    do_not_annotate: set[str],
    aspects: OntologyMap,
    ontologies: tuple[str, ...],
) -> ScoredPredictionByOntology:
    """Parse scored predictions once, applying fixed filters and splitting by ontology."""

    prediction_path = Path(path)
    delimiter = delimiter_from_format(prediction_path, input_format)
    ontology_set = set(ontologies)
    predictions: dict[str, defaultdict[str, dict[str, float]]] = {
        ontology: defaultdict(dict) for ontology in ontologies
    }
    total_rows = 0
    kept_rows = 0
    skipped_unknown_ontology = 0
    skipped_fixed_filter = 0
    duplicate_replacements = 0

    LOGGER.info("Parsing scored predictions from %s", prediction_path)

    with prediction_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header_checked = False
        skip_first = False
        for row_number, row in enumerate(reader, start=1):
            if not row or all(not field.strip() for field in row):
                continue
            if row[0].lstrip().startswith("#"):
                continue
            if not header_checked:
                skip_first = looks_like_header(row, protein_col, go_col, score_col)
                header_checked = True
                if skip_first:
                    continue
            if len(row) <= max(protein_col, go_col, score_col):
                raise ValueError(f"{prediction_path}: row {row_number} has too few columns")

            total_rows += 1
            protein = row[protein_col].strip()
            go_id = row[go_col].strip()
            if not protein or not go_id:
                skipped_fixed_filter += 1
                continue
            if go_id in BINDING_TERMS:
                skipped_fixed_filter += 1
                continue
            if go_id in existing.get(protein, set()):
                skipped_fixed_filter += 1
                continue
            if go_id in do_not_annotate:
                skipped_fixed_filter += 1
                continue
            ontology = aspects.get(go_id)
            if ontology not in ontology_set:
                skipped_unknown_ontology += 1
                continue
            try:
                score = float(row[score_col])
            except ValueError as exc:
                raise ValueError(f"{prediction_path}: row {row_number} has non-numeric score {row[score_col]!r}") from exc

            previous = predictions[ontology][protein].get(go_id)
            if previous is None or score > previous:
                if previous is not None:
                    duplicate_replacements += 1
                predictions[ontology][protein][go_id] = score
                kept_rows += 1

    parsed = {
        ontology: {protein: dict(terms) for protein, terms in ontology_predictions.items()}
        for ontology, ontology_predictions in predictions.items()
    }
    LOGGER.info(
        "Parsed %s: rows=%d kept_or_replaced=%d fixed_filter_skipped=%d unknown_ontology_skipped=%d duplicate_replacements=%d",
        prediction_path.name,
        total_rows,
        kept_rows,
        skipped_fixed_filter,
        skipped_unknown_ontology,
        duplicate_replacements,
    )
    for ontology in ontologies:
        n_proteins = len(parsed[ontology])
        n_terms = sum(len(terms) for terms in parsed[ontology].values())
        LOGGER.info("Prediction terms for %s/%s: proteins=%d terms=%d", prediction_path.name, ontology, n_proteins, n_terms)
    return parsed


def threshold_values(tau_min: Decimal, tau_max: Decimal, tau_step: Decimal) -> list[Decimal]:
    if tau_step <= 0:
        raise ValueError("--tau-step must be positive")
    if tau_min > tau_max:
        raise ValueError("--tau-min must be <= --tau-max")

    values: list[Decimal] = []
    tau = tau_max - tau_step # not including tau_max itself
    while tau >= tau_min:
        values.append(tau)
        tau -= tau_step
    if not values or values[-1] != tau_min:
        values.append(tau_min)
    return values


def filter_annotations_by_ontology(annotations: AnnotationMap, aspects: OntologyMap, ontology: str) -> AnnotationMap:
    result: AnnotationMap = {}
    for protein, terms in annotations.items():
        kept = {go_id for go_id in terms if aspects.get(go_id) == ontology}
        if kept:
            result[protein] = kept
    return result


def active_predictions_for_tau(
    predictions: ScoredPredictionMap,
    tau: float,
) -> AnnotationMap:
    active: AnnotationMap = {}
    for protein, scored_terms in predictions.items():
        kept = {go_id for go_id, score in scored_terms.items() if score >= tau}
        if kept:
            active[protein] = kept
    return active


def f1_score(precision: float, recall: float) -> float:
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0


def metric_counts(
    prediction_map: dict,
    prediction_map_count: dict,
    reference_nomap: dict[str, set[str]],
    reference: AnnotationMap,
) -> dict[str, Any]:
    per_protein: dict[str, dict[str, int]] = {}
    reference_proteins = set(reference)

    for protein in set(prediction_map) | reference_proteins:
        e_count = len(prediction_map.get(protein, {}).get("direct", {}))
        l_count = len(prediction_map_count.get(protein, set()))
        m_count = len(prediction_map.get(protein, {}).get("related", {}))
        u_count = len(prediction_map.get(protein, {}).get("unrelated", {}))
        z_count = len(reference_nomap.get(protein, set()))
        per_protein[protein] = {"E": e_count, "L": l_count, "M": m_count, "U": u_count, "Z": z_count}

    totals = {
        key: sum(values[key] for values in per_protein.values())
        for key in ("E", "L", "M", "U", "Z")
    }
    return {"per_protein": per_protein, **totals}


def summarize_counts(counts: dict[str, Any], prediction_map: dict, reference: AnnotationMap) -> dict[str, float | int]:
    e_count = int(counts["E"])
    l_count = int(counts["L"])
    m_count = int(counts["M"])
    u_count = int(counts["U"])
    z_count = int(counts["Z"])

    numerator = e_count + 0.75 * l_count + 0.5 * m_count
    micro_precision_denominator = e_count + l_count + m_count + u_count
    micro_recall_denominator = e_count + l_count + m_count + z_count
    micro_precision = numerator / micro_precision_denominator if micro_precision_denominator else 0.0
    micro_recall = numerator / micro_recall_denominator if micro_recall_denominator else 0.0

    precision_sum = 0.0
    recall_sum = 0.0
    for protein, values in counts["per_protein"].items():
        protein_numerator = values["E"] + 0.75 * values["L"] + 0.5 * values["M"]
        protein_precision_denominator = values["E"] + values["L"] + values["M"] + values["U"]
        protein_recall_denominator = values["E"] + values["L"] + values["M"] + values["Z"]
        if protein_precision_denominator:
            precision_sum += protein_numerator / protein_precision_denominator
        if protein_recall_denominator:
            recall_sum += protein_numerator / protein_recall_denominator

    n_mapped_prediction_proteins = len(prediction_map)
    n_reference_proteins = len(reference)
    macro_precision = precision_sum / n_mapped_prediction_proteins if n_mapped_prediction_proteins else 0.0
    macro_recall = recall_sum / n_reference_proteins if n_reference_proteins else 0.0

    return {
        "pr_micro": micro_precision,
        "rc_micro": micro_recall,
        "f_micro": f1_score(micro_precision, micro_recall),
        "pr_macro": macro_precision,
        "rc_macro": macro_recall,
        "f_macro": f1_score(macro_precision, macro_recall),
        "E": e_count,
        "L": l_count,
        "M": m_count,
        "U": u_count,
        "Z": z_count,
        "n_reference_proteins": n_reference_proteins,
        "n_mapped_prediction_proteins": n_mapped_prediction_proteins,
    }


def evaluate_tau_ontology(
    *,
    filename: str,
    tau_text: str,
    tau: float,
    ontology: str,
    predictions: ScoredPredictionMap,
    reference_by_ontology: dict[str, AnnotationMap],
    go_parents: ParentMap,
) -> dict[str, Any]:
    active = active_predictions_for_tau(predictions, tau)
    predicted = remove_redundancy_fast(active, go_parents, "predicted")
    reference = reference_by_ontology[ontology]
    prediction_map, prediction_map_count, reference_nomap = evaluate_fast_mapping(predicted, reference, go_parents)
    counts = metric_counts(prediction_map, prediction_map_count, reference_nomap, reference)
    summary = summarize_counts(counts, prediction_map, reference)
    return {
        "filename": filename,
        "tau": tau_text,
        "ns": ontology,
        "n_predicted_proteins": len(predicted),
        **summary,
    }


def init_worker(context: dict[str, Any]) -> None:
    global _WORKER_CONTEXT
    _WORKER_CONTEXT = context
    set_metric_logging(False)


def worker_evaluate(task: tuple[str, float, str]) -> dict[str, Any]:
    tau_text, tau, ontology = task
    return evaluate_tau_ontology(
        filename=_WORKER_CONTEXT["filename"],
        tau_text=tau_text,
        tau=tau,
        ontology=ontology,
        predictions=_WORKER_CONTEXT["predictions_by_ontology"][ontology],
        reference_by_ontology=_WORKER_CONTEXT["reference_by_ontology"],
        go_parents=_WORKER_CONTEXT["go_parents"],
    )


def format_float(value: Any) -> Any:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.3f}"
    return value


def best_rows_by_score(rows: list[dict[str, Any]], score_field: str) -> list[dict[str, Any]]:
    """Pick one best-threshold row per filename and ontology."""

    best: dict[tuple[str, str], dict[str, Any]] = {}
    for row in rows:
        key = (str(row["filename"]), str(row["ns"]))
        current = best.get(key)
        if current is None:
            best[key] = row
            continue
        if float(row[score_field]) > float(current[score_field]):
            best[key] = row
            continue
        if float(row[score_field]) == float(current[score_field]) and float(row["tau"]) > float(current["tau"]):
            best[key] = row
    return [best[key] for key in sorted(best)]


def write_rows(path: str | Path, rows: list[dict[str, Any]]) -> None:
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_float(row.get(field, "")) for field in OUTPUT_FIELDS})


def prediction_files(args: argparse.Namespace) -> list[Path]:
    if args.prediction:
        return [Path(args.prediction)]
    return sorted(
        path
        for path in Path(args.prediction_dir).rglob(args.glob)
        if path.is_file() and not path.name.startswith(".")
    )


def evaluate_prediction_file(
    prediction_file: Path,
    *,
    args: argparse.Namespace,
    existing: AnnotationMap,
    do_not_annotate: set[str],
    reference_by_ontology: dict[str, AnnotationMap],
    go_parents: ParentMap,
    aspects: OntologyMap,
    ontologies: tuple[str, ...],
    taus: list[Decimal],
) -> list[dict[str, Any]]:
    start_time = time.perf_counter()
    predictions_by_ontology = parse_scored_predictions(
        prediction_file,
        input_format=args.predicted_format,
        protein_col=args.protein_col,
        go_col=args.go_col,
        score_col=args.score_col,
        existing=existing,
        do_not_annotate=do_not_annotate,
        aspects=aspects,
        ontologies=ontologies,
    )
    tasks = [(f"{float(tau):.4f}", float(tau), ontology) for tau in taus for ontology in ontologies]
    context = {
        "filename": prediction_file.name,
        "predictions_by_ontology": predictions_by_ontology,
        "reference_by_ontology": reference_by_ontology,
        "go_parents": go_parents,
    }

    if args.workers > 1 and len(tasks) > 1:
        LOGGER.info("Running %d tau/ontology jobs for %s with %d workers", len(tasks), prediction_file.name, args.workers)
        with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker, initargs=(context,)) as executor:
            rows = list(executor.map(worker_evaluate, tasks))
        LOGGER.info("Completed %s in %.1f seconds", prediction_file.name, time.perf_counter() - start_time)
        return rows

    init_worker(context)
    LOGGER.info("Running %d tau/ontology jobs for %s serially", len(tasks), prediction_file.name)
    rows = [worker_evaluate(task) for task in tasks]
    LOGGER.info("Completed %s in %.1f seconds", prediction_file.name, time.perf_counter() - start_time)
    return rows


def main() -> int:
    start_time = time.perf_counter()
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(message)s")
    set_metric_logging(False)

    LOGGER.info("Loading GO parent lookup: %s", args.go_parent)
    go_parents = parse_go_parents(args.go_parent)
    aspects = parse_go_parent_aspects(args.go_parent)
    LOGGER.info("Loaded %d GO terms with parent entries and %d GO term aspects", len(go_parents), len(aspects))
    LOGGER.info("Loading annotation filters and reference annotations")
    do_not_annotate = parse_do_not_annotate(args.do_not_annotate)
    existing = parse_existing_annotations(args.existing)
    raw_reference = parse_reference_annotations(args.reference, existing, do_not_annotate)
    reference = remove_redundancy_fast(raw_reference, go_parents, "reference")

    ontologies = tuple(args.ontology) if args.ontology else ONTOLOGIES
    LOGGER.info("Evaluating ontologies: %s", ", ".join(ontologies))
    reference_by_ontology = {
        ontology: filter_annotations_by_ontology(reference, aspects, ontology)
        for ontology in ontologies
    }
    taus = threshold_values(args.tau_min, args.tau_max, args.tau_step)
    LOGGER.info("Evaluating %d thresholds from %s to %s with step %s", len(taus), args.tau_min, args.tau_max, args.tau_step)
    rows: list[dict[str, Any]] = []
    for prediction_file in prediction_files(args):
        LOGGER.info("Evaluating prediction file: %s", prediction_file)
        rows.extend(
            evaluate_prediction_file(
                prediction_file,
                args=args,
                existing=existing,
                do_not_annotate=do_not_annotate,
                reference_by_ontology=reference_by_ontology,
                go_parents=go_parents,
                aspects=aspects,
                ontologies=ontologies,
                taus=taus,
            )
        )

    write_rows(args.out_file, rows)
    LOGGER.info("Wrote full sweep output: %s (%d rows)", args.out_file, len(rows))
    if args.best_macro_output:
        best_macro_rows = best_rows_by_score(rows, "macro_f1")
        write_rows(args.best_macro_output, best_macro_rows)
        LOGGER.info("Wrote best macro output: %s (%d rows)", args.best_macro_output, len(best_macro_rows))
    if args.best_micro_output:
        best_micro_rows = best_rows_by_score(rows, "micro_f1")
        write_rows(args.best_micro_output, best_micro_rows)
        LOGGER.info("Wrote best micro output: %s (%d rows)", args.best_micro_output, len(best_micro_rows))
    LOGGER.info("Finished sweep in %.1f seconds", time.perf_counter() - start_time)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
