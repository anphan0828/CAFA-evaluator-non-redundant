#!/usr/bin/env python3
"""Threshold-sweep HGRS metric with ontology-stratified micro/macro summaries."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import time
from bisect import bisect_right
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
from pathlib import Path
from typing import Any, NamedTuple

from fpr_metric_common import (
    BINDING_TERMS,
    GO_RE,
    GO_ROOT_TERMS,
    AnnotationMap,
    OntologyMap,
    ParentMap,
    evaluate_fast_mapping,
    parse_do_not_annotate,
    parse_existing_annotations,
    parse_go_parent_aspects,
    parse_go_parents,
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


class ProteinMetric(NamedTuple):
    e_count: int
    l_count: int
    m_count: int
    u_count: int
    z_count: int
    precision: float
    recall: float


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
    parser.add_argument(
        "--terms-of-interest",
        help="optional file with one eligible GO ID per line; GO terms not listed are excluded",
    )
    parser.add_argument("-g", "--go-parent", required=True, help="goparents lookup file")
    parser.add_argument("-o", "--out-file", required=True, help="threshold sweep summary TSV")
    parser.add_argument("--best-macro-output", help="best-threshold summary TSV selected by f_macro")
    parser.add_argument("--best-micro-output", help="best-threshold summary TSV selected by f_micro")
    parser.add_argument("--map-output", help="single-file debug mapping output for the lowest-threshold rows")
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
    parser.add_argument("--workers", type=int, default=1, help="process workers for ontology-level incremental sweeps")
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


def parse_terms_of_interest(path: str | Path | None) -> set[str] | None:
    if not path:
        return None

    terms: set[str] = set()
    with Path(path).open() as handle:
        for line in handle:
            go_id = line.strip().split()[0] if line.strip() else ""
            if go_id:
                terms.add(go_id)
    LOGGER.info("Loaded terms-of-interest: %s (%d terms)", path, len(terms))
    return terms


def parse_reference_annotations_for_sweep(
    path: str | Path,
    existing: AnnotationMap,
    do_not_annotate: set[str],
    terms_of_interest: set[str] | None,
) -> AnnotationMap:
    reference: defaultdict[str, set[str]] = defaultdict(set)
    existing_skipped = 0
    binding_skipped = 0
    root_skipped = 0
    do_not_annotate_skipped = 0
    not_terms_of_interest_skipped = 0
    kept = 0

    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            protein, go_id = fields[0], fields[1]
            if go_id in existing.get(protein, set()):
                existing_skipped += 1
                continue
            if go_id in BINDING_TERMS:
                binding_skipped += 1
                continue
            if go_id in GO_ROOT_TERMS:
                root_skipped += 1
                continue
            if terms_of_interest is not None and go_id not in terms_of_interest:
                not_terms_of_interest_skipped += 1
                continue
            if go_id in do_not_annotate:
                do_not_annotate_skipped += 1
                continue
            reference[protein].add(go_id)
            kept += 1

    LOGGER.info(
        "Parsed reference annotations: kept=%d existing_skipped=%d binding_skipped=%d "
        "root_skipped=%d obsolete_or_not_terms_of_interest_removed=%d do_not_annotate_skipped=%d",
        kept,
        existing_skipped,
        binding_skipped,
        root_skipped,
        not_terms_of_interest_skipped,
        do_not_annotate_skipped,
    )
    if terms_of_interest is not None:
        LOGGER.info(
            "Reference terms removed as obsolete/outside terms-of-interest: %d",
            not_terms_of_interest_skipped,
        )
    return {protein: set(terms) for protein, terms in reference.items()}


def parse_scored_predictions(
    path: str | Path,
    *,
    input_format: str,
    protein_col: int,
    go_col: int,
    score_col: int,
    existing: AnnotationMap,
    do_not_annotate: set[str],
    terms_of_interest: set[str] | None,
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
    skipped_root = 0
    skipped_not_terms_of_interest = 0
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
            if go_id in GO_ROOT_TERMS:
                skipped_root += 1
                continue
            if terms_of_interest is not None and go_id not in terms_of_interest:
                skipped_not_terms_of_interest += 1
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
        "Parsed %s: rows=%d kept_or_replaced=%d fixed_filter_skipped=%d "
        "root_skipped=%d obsolete_or_not_terms_of_interest_removed=%d "
        "unknown_ontology_skipped=%d duplicate_replacements=%d",
        prediction_path.name,
        total_rows,
        kept_rows,
        skipped_fixed_filter,
        skipped_root,
        skipped_not_terms_of_interest,
        skipped_unknown_ontology,
        duplicate_replacements,
    )
    if terms_of_interest is not None:
        LOGGER.info(
            "Prediction terms removed as obsolete/outside terms-of-interest: %d",
            skipped_not_terms_of_interest,
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


def bucket_predictions_by_threshold(
    predictions: ScoredPredictionMap,
    taus: list[Decimal],
) -> list[defaultdict[str, list[str]]]:
    """Assign each prediction to the first descending threshold it passes."""

    ascending_thresholds = sorted((float(tau), index) for index, tau in enumerate(taus))
    threshold_values_ascending = [value for value, _ in ascending_thresholds]
    threshold_indices = [index for _, index in ascending_thresholds]
    buckets: list[defaultdict[str, list[str]]] = [defaultdict(list) for _ in taus]

    for protein, scored_terms in predictions.items():
        for go_id, score in scored_terms.items():
            if math.isnan(score):
                continue
            threshold_position = bisect_right(threshold_values_ascending, score) - 1
            if threshold_position < 0:
                continue
            buckets[threshold_indices[threshold_position]][protein].append(go_id)
    return buckets


def nonredundant_terms_for_protein(terms: set[str], go_parents: ParentMap) -> set[str]:
    redundant: set[str] = set()
    for go_id in terms:
        redundant.update(go_parents.get(go_id, set()) & terms)
    return terms - redundant


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


def metric_for_protein(
    predicted_terms: set[str],
    reference_terms: set[str],
    go_parents: ParentMap,
) -> ProteinMetric:
    """Calculate the mapping counts for one protein without materializing maps."""

    e_count = 0
    m_count = 0
    u_count = 0
    true_sets: set[frozenset[str]] = set()
    reference_map: set[str] = set()

    true_by_predicted: defaultdict[str, set[str]] = defaultdict(set)
    for go_r in reference_terms:
        for parent in go_parents.get(go_r, set()):
            if parent in predicted_terms:
                true_by_predicted[parent].add(go_r)

    for go_p in predicted_terms:
        if go_p in reference_terms:
            e_count += 1
            reference_map.add(go_p)
            continue

        true = true_by_predicted.get(go_p)
        if true:
            true_sets.add(frozenset(true))
            reference_map.update(true)
            continue

        related = go_parents.get(go_p, set()) & reference_terms
        if related:
            m_count += 1
            reference_map.update(related)
        else:
            u_count += 1

    l_count = len(true_sets)
    z_count = len(reference_terms - reference_map)
    numerator = e_count + 0.75 * l_count + 0.5 * m_count
    precision_denominator = e_count + l_count + m_count + u_count
    recall_denominator = e_count + l_count + m_count + z_count
    precision = numerator / precision_denominator if precision_denominator else 0.0
    recall = numerator / recall_denominator if recall_denominator else 0.0
    return ProteinMetric(
        e_count=e_count,
        l_count=l_count,
        m_count=m_count,
        u_count=u_count,
        z_count=z_count,
        precision=precision,
        recall=recall,
    )


def row_from_incremental_counts(
    *,
    filename: str,
    tau_text: str,
    ontology: str,
    n_predicted_proteins: int,
    n_reference_proteins: int,
    reference_term_count: int,
    protein_metrics: dict[str, ProteinMetric],
    totals: dict[str, float],
) -> dict[str, Any]:
    e_count = int(totals["E"])
    l_count = int(totals["L"])
    m_count = int(totals["M"])
    u_count = int(totals["U"])
    z_count = int(reference_term_count + totals["Z_DELTA"])

    numerator = e_count + 0.75 * l_count + 0.5 * m_count
    micro_precision_denominator = e_count + l_count + m_count + u_count
    micro_recall_denominator = e_count + l_count + m_count + z_count
    micro_precision = numerator / micro_precision_denominator if micro_precision_denominator else 0.0
    micro_recall = numerator / micro_recall_denominator if micro_recall_denominator else 0.0

    n_mapped_prediction_proteins = len(protein_metrics)
    macro_precision = totals["PRECISION_SUM"] / n_mapped_prediction_proteins if n_mapped_prediction_proteins else 0.0
    macro_recall = totals["RECALL_SUM"] / n_reference_proteins if n_reference_proteins else 0.0

    return {
        "filename": filename,
        "tau": tau_text,
        "ns": ontology,
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
        "n_predicted_proteins": n_predicted_proteins,
        "n_reference_proteins": n_reference_proteins,
        "n_mapped_prediction_proteins": n_mapped_prediction_proteins,
    }


def evaluate_prediction_snapshot(
    *,
    filename: str,
    tau_text: str,
    ontology: str,
    predictions_for_mapping: AnnotationMap,
    n_predicted_proteins: int,
    reference_by_ontology: dict[str, AnnotationMap],
    go_parents: ParentMap,
) -> dict[str, Any]:
    reference = reference_by_ontology[ontology]
    prediction_map, prediction_map_count, reference_nomap = evaluate_fast_mapping(predictions_for_mapping, reference, go_parents)
    counts = metric_counts(prediction_map, prediction_map_count, reference_nomap, reference)
    summary = summarize_counts(counts, prediction_map, reference)
    return {
        "filename": filename,
        "tau": tau_text,
        "ns": ontology,
        "n_predicted_proteins": n_predicted_proteins,
        **summary,
    }


def sweep_ontology(
    *,
    filename: str,
    ontology: str,
    predictions: ScoredPredictionMap,
    reference_by_ontology: dict[str, AnnotationMap],
    go_parents: ParentMap,
    taus: list[Decimal],
) -> list[dict[str, Any]]:
    LOGGER.debug("Preparing incremental threshold buckets for %s/%s", filename, ontology)
    buckets = bucket_predictions_by_threshold(predictions, taus)
    active_by_protein: defaultdict[str, set[str]] = defaultdict(set)
    nonredundant_by_protein: AnnotationMap = {}
    reference = reference_by_ontology[ontology]
    reference_proteins = set(reference)
    reference_term_count = sum(len(terms) for terms in reference.values())
    protein_metrics: dict[str, ProteinMetric] = {}
    totals = {"E": 0.0, "L": 0.0, "M": 0.0, "U": 0.0, "Z_DELTA": 0.0, "PRECISION_SUM": 0.0, "RECALL_SUM": 0.0}
    rows: list[dict[str, Any]] = []

    for tau_index, tau in enumerate(taus):
        dirty_proteins: set[str] = set()
        for protein, go_ids in buckets[tau_index].items():
            active_by_protein[protein].update(go_ids)
            dirty_proteins.add(protein)

        for protein in dirty_proteins:
            kept = nonredundant_terms_for_protein(active_by_protein[protein], go_parents)
            if kept:
                nonredundant_by_protein[protein] = kept
            else:
                nonredundant_by_protein.pop(protein, None)

            if protein not in reference_proteins:
                continue

            old_metric = protein_metrics.pop(protein, None)
            if old_metric is not None:
                totals["E"] -= old_metric.e_count
                totals["L"] -= old_metric.l_count
                totals["M"] -= old_metric.m_count
                totals["U"] -= old_metric.u_count
                totals["Z_DELTA"] -= old_metric.z_count - len(reference[protein])
                totals["PRECISION_SUM"] -= old_metric.precision
                totals["RECALL_SUM"] -= old_metric.recall

            if kept:
                new_metric = metric_for_protein(kept, reference[protein], go_parents)
                protein_metrics[protein] = new_metric
                totals["E"] += new_metric.e_count
                totals["L"] += new_metric.l_count
                totals["M"] += new_metric.m_count
                totals["U"] += new_metric.u_count
                totals["Z_DELTA"] += new_metric.z_count - len(reference[protein])
                totals["PRECISION_SUM"] += new_metric.precision
                totals["RECALL_SUM"] += new_metric.recall

        rows.append(
            row_from_incremental_counts(
                filename=filename,
                tau_text=f"{float(tau):.4f}",
                ontology=ontology,
                n_predicted_proteins=len(nonredundant_by_protein),
                n_reference_proteins=len(reference),
                reference_term_count=reference_term_count,
                protein_metrics=protein_metrics,
                totals=totals,
            )
        )
    return rows


def init_worker(context: dict[str, Any]) -> None:
    global _WORKER_CONTEXT
    _WORKER_CONTEXT = context
    set_metric_logging(False)


def worker_sweep_ontology(ontology: str) -> list[dict[str, Any]]:
    return sweep_ontology(
        filename=_WORKER_CONTEXT["filename"],
        ontology=ontology,
        predictions=_WORKER_CONTEXT["predictions_by_ontology"][ontology],
        reference_by_ontology=_WORKER_CONTEXT["reference_by_ontology"],
        go_parents=_WORKER_CONTEXT["go_parents"],
        taus=_WORKER_CONTEXT["taus"],
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


def lowest_threshold_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Pick the lowest-threshold row per filename and ontology."""

    selected: dict[tuple[str, str], dict[str, Any]] = {}
    for row in rows:
        key = (str(row["filename"]), str(row["ns"]))
        current = selected.get(key)
        if current is None or float(row["tau"]) < float(current["tau"]):
            selected[key] = row
    return [selected[key] for key in sorted(selected)]


def write_rows(path: str | Path, rows: list[dict[str, Any]]) -> None:
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_float(row.get(field, "")) for field in OUTPUT_FIELDS})


def nonredundant_predictions_at_tau(
    predictions: ScoredPredictionMap,
    tau: float,
    go_parents: ParentMap,
) -> AnnotationMap:
    active: AnnotationMap = {}
    for protein, scored_terms in predictions.items():
        kept = {go_id for go_id, score in scored_terms.items() if score >= tau}
        if kept:
            nonredundant = nonredundant_terms_for_protein(kept, go_parents)
            if nonredundant:
                active[protein] = nonredundant
    return active


def write_mapping_output(
    path: str | Path,
    *,
    filename: str,
    selected_rows: list[dict[str, Any]],
    predictions_by_ontology: ScoredPredictionByOntology,
    reference_by_ontology: dict[str, AnnotationMap],
    go_parents: ParentMap,
) -> None:
    fieldnames = ["filename", "ns", "tau", "id", "predicted", "map type", "reference"]
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in selected_rows:
            ontology = str(row["ns"])
            tau_text = str(row["tau"])
            tau = float(tau_text)
            reference = reference_by_ontology[ontology]
            predicted = nonredundant_predictions_at_tau(
                predictions_by_ontology[ontology],
                tau,
                go_parents,
            )
            predicted_for_mapping = {
                protein: terms
                for protein, terms in predicted.items()
                if protein in reference
            }
            prediction_map, _, _ = evaluate_fast_mapping(predicted_for_mapping, reference, go_parents)
            for protein in sorted(prediction_map):
                for map_type in sorted(prediction_map[protein]):
                    for go_id in sorted(prediction_map[protein][map_type]):
                        mapped = prediction_map[protein][map_type][go_id]
                        writer.writerow(
                            {
                                "filename": filename,
                                "ns": ontology,
                                "tau": tau_text,
                                "id": protein,
                                "predicted": go_id,
                                "map type": map_type,
                                "reference": str(mapped) if str(mapped).startswith("GO") else "",
                            }
                        )


def prediction_files(args: argparse.Namespace) -> list[Path]:
    if args.prediction:
        return [Path(args.prediction)]
    return sorted(
        path
        for path in Path(args.prediction_dir).rglob(args.glob)
        if path.is_file() and not path.name.startswith(".")
    )


def nonempty_prediction_ontologies(
    predictions_by_ontology: ScoredPredictionByOntology,
    ontologies: tuple[str, ...],
    filename: str,
) -> tuple[str, ...]:
    kept: list[str] = []
    for ontology in ontologies:
        predictions = predictions_by_ontology[ontology]
        n_proteins = len(predictions)
        n_terms = sum(len(terms) for terms in predictions.values())
        if n_proteins == 0 and n_terms == 0:
            LOGGER.info("Skipping %s/%s after filtering: proteins=0 terms=0", filename, ontology)
            continue
        kept.append(ontology)
    return tuple(kept)


def evaluate_prediction_file(
    prediction_file: Path,
    *,
    args: argparse.Namespace,
    existing: AnnotationMap,
    do_not_annotate: set[str],
    terms_of_interest: set[str] | None,
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
        terms_of_interest=terms_of_interest,
        aspects=aspects,
        ontologies=ontologies,
    )
    active_ontologies = nonempty_prediction_ontologies(predictions_by_ontology, ontologies, prediction_file.name)
    if not active_ontologies:
        LOGGER.info("Skipping %s after filtering: no ontologies with predictions", prediction_file.name)
        return []

    context = {
        "filename": prediction_file.name,
        "predictions_by_ontology": predictions_by_ontology,
        "reference_by_ontology": reference_by_ontology,
        "go_parents": go_parents,
        "taus": taus,
    }

    if args.workers > 1 and len(active_ontologies) > 1:
        worker_count = min(args.workers, len(active_ontologies))
        LOGGER.debug(
            "Running %d incremental ontology sweeps for %s with %d workers",
            len(active_ontologies),
            prediction_file.name,
            worker_count,
        )
        with ProcessPoolExecutor(max_workers=worker_count, initializer=init_worker, initargs=(context,)) as executor:
            row_groups = list(executor.map(worker_sweep_ontology, active_ontologies))
        rows = [row for group in row_groups for row in group]
    else:
        init_worker(context)
        LOGGER.debug(
            "Running %d incremental ontology sweeps for %s serially",
            len(active_ontologies),
            prediction_file.name,
        )
        rows = []
        for ontology in active_ontologies:
            rows.extend(worker_sweep_ontology(ontology))
        _WORKER_CONTEXT.clear()

    tau_order = {f"{float(tau):.4f}": index for index, tau in enumerate(taus)}
    ontology_order = {ontology: index for index, ontology in enumerate(ontologies)}
    rows.sort(key=lambda row: (tau_order[str(row["tau"])], ontology_order[str(row["ns"])]))
    if args.map_output:
        selected_rows = lowest_threshold_rows(rows)
        write_mapping_output(
            args.map_output,
            filename=prediction_file.name,
            selected_rows=selected_rows,
            predictions_by_ontology=predictions_by_ontology,
            reference_by_ontology=reference_by_ontology,
            go_parents=go_parents,
        )
        LOGGER.info(
            "Wrote debug mapping output for lowest threshold: %s (%d ontology rows)",
            args.map_output,
            len(selected_rows),
        )
    del predictions_by_ontology
    context.clear()
    LOGGER.info("Completed %s in %.1f seconds", prediction_file.name, time.perf_counter() - start_time)
    return rows


def main() -> int:
    start_time = time.perf_counter()
    args = parse_args()
    if args.map_output and not args.prediction:
        raise SystemExit("--map-output is only supported with single-file input (-p/--prediction)")
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(message)s")
    set_metric_logging(False)

    # LOGGER.info("Loading GO parent lookup: %s", args.go_parent)
    go_parents = parse_go_parents(args.go_parent)
    aspects = parse_go_parent_aspects(args.go_parent)
    LOGGER.info("Loaded %d GO terms with parent entries and %d GO term aspects", len(go_parents), len(aspects))
    LOGGER.debug("Loading annotation filters and reference annotations")
    do_not_annotate = parse_do_not_annotate(args.do_not_annotate)
    terms_of_interest = parse_terms_of_interest(args.terms_of_interest)
    existing = parse_existing_annotations(args.existing)
    raw_reference = parse_reference_annotations_for_sweep(
        args.reference,
        existing,
        do_not_annotate,
        terms_of_interest,
    )
    reference = remove_redundancy_fast(raw_reference, go_parents, "reference")

    ontologies = tuple(args.ontology) if args.ontology else ONTOLOGIES
    LOGGER.debug("Evaluating ontologies: %s", ", ".join(ontologies))
    reference_by_ontology = {
        ontology: filter_annotations_by_ontology(reference, aspects, ontology)
        for ontology in ontologies
    }
    taus = threshold_values(args.tau_min, args.tau_max, args.tau_step)
    LOGGER.info("Evaluating %d thresholds from %s to %s with step %s", len(taus), args.tau_min, args.tau_max, args.tau_step)
    rows: list[dict[str, Any]] = []
    for prediction_file in prediction_files(args):
        LOGGER.debug("Evaluating prediction file: %s", prediction_file)
        rows.extend(
            evaluate_prediction_file(
                prediction_file,
                args=args,
                existing=existing,
                do_not_annotate=do_not_annotate,
                terms_of_interest=terms_of_interest,
                reference_by_ontology=reference_by_ontology,
                go_parents=go_parents,
                aspects=aspects,
                ontologies=ontologies,
                taus=taus,
            )
        )

    write_rows(args.out_file, rows)
    LOGGER.debug("Wrote full sweep output: %s (%d rows)", args.out_file, len(rows))
    if args.best_macro_output:
        best_macro_rows = best_rows_by_score(rows, "f_macro")
        write_rows(args.best_macro_output, best_macro_rows)
        LOGGER.debug("Wrote best macro output: %s (%d rows)", args.best_macro_output, len(best_macro_rows))
    if args.best_micro_output:
        best_micro_rows = best_rows_by_score(rows, "f_micro")
        write_rows(args.best_micro_output, best_micro_rows)
        LOGGER.debug("Wrote best micro output: %s (%d rows)", args.best_micro_output, len(best_micro_rows))
    LOGGER.info("Finished sweep in %.1f seconds", time.perf_counter() - start_time)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
