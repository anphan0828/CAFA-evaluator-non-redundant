#!/usr/bin/env python3
"""Threshold-sweep HGRS metric with obsolete-term prediction rescue."""

from __future__ import annotations

import argparse
import csv
import logging
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
from pathlib import Path
from typing import NamedTuple

import fpr_metric_sweep as sweep
from fpr_metric_common import (
    BINDING_TERMS,
    GO_ROOT_TERMS,
    AnnotationMap,
    OntologyMap,
    ParentMap,
    extract_go_id,
    parse_do_not_annotate,
    parse_existing_annotations,
    remove_redundancy_fast,
    set_metric_logging,
)


LOGGER = logging.getLogger("fpr_metric_sweep_rescue")


class ObsoleteInfo(NamedTuple):
    obsolete: set[str]
    replaced_by: dict[str, set[str]]


class RescueEntry(NamedTuple):
    targets: set[str]
    reason: str


class FilteredRescueEntry(NamedTuple):
    targets: tuple[tuple[str, str], ...]
    reason: str
    fixed_filtered: int
    root_filtered: int
    not_toi_filtered: int
    do_not_filtered: int
    unknown_ontology_filtered: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep score thresholds for the hierarchy-aware GO recovery metric, "
            "with prediction rescue for obsolete GO terms."
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
        required=True,
        help="file with one eligible t1 GO ID per line; used to identify future-obsolete t0 terms",
    )
    parser.add_argument("-g", "--go-parent", required=True, help="t0 goparents lookup file")
    parser.add_argument("--t0-obo", required=True, help="t0 GO OBO file; used for already-obsolete replaced_by rescue")
    parser.add_argument("--t1-obo", required=True, help="t1 GO OBO file; used for future-obsolete replaced_by rescue")
    parser.add_argument("--out-dir", required=True, help="directory for rescue evaluation output TSV files")
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
        choices=sweep.ONTOLOGIES,
        help="ontology to evaluate; repeatable. Defaults to all three GO ontologies.",
    )
    parser.add_argument("--workers", type=int, default=1, help="process workers for ontology-level incremental sweeps")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    return parser.parse_args()


def parse_obo_obsolete_info(path: str | Path) -> ObsoleteInfo:
    obsolete: set[str] = set()
    replaced_by: dict[str, set[str]] = {}
    current_id: str | None = None
    current_obsolete = False
    current_replacements: set[str] = set()
    in_term = False

    def flush_term() -> None:
        if current_id and current_obsolete:
            obsolete.add(current_id)
            if current_replacements:
                replaced_by[current_id] = set(current_replacements)

    with Path(path).open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                flush_term()
                current_id = None
                current_obsolete = False
                current_replacements = set()
                in_term = True
                continue
            if line.startswith("[") and line != "[Term]":
                flush_term()
                current_id = None
                current_obsolete = False
                current_replacements = set()
                in_term = False
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                current_id = line.split("id: ", 1)[1].strip()
            elif line == "is_obsolete: true":
                current_obsolete = True
            elif line.startswith("replaced_by: "):
                go_id = extract_go_id(line)
                if go_id:
                    current_replacements.add(go_id)
        flush_term()

    LOGGER.debug(
        "Parsed OBO obsolete info from %s: obsolete=%d replaced_by=%d",
        path,
        len(obsolete),
        len(replaced_by),
    )
    return ObsoleteInfo(obsolete=obsolete, replaced_by=replaced_by)


def parse_go_parent_lookups(path: str | Path) -> tuple[ParentMap, OntologyMap, ParentMap]:
    parents: defaultdict[str, set[str]] = defaultdict(set)
    aspects: dict[str, str] = {}
    direct: defaultdict[str, set[str]] = defaultdict(set)

    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            child_id = extract_go_id(fields[0])
            parent_id = extract_go_id(fields[1])
            if not child_id or not parent_id:
                continue

            parents[child_id].add(parent_id)

            if len(fields) >= 4 and fields[3]:
                aspect = fields[3]
                aspects[child_id] = aspect
                aspects[parent_id] = aspect

            if len(fields) >= 3 and fields[2].strip() in {"is_a", "part_of"}:
                direct[child_id].add(parent_id)

    LOGGER.debug(
        "Loaded GO parent lookups from %s: parent_children=%d aspects=%d direct_parent_children=%d",
        path,
        len(parents),
        len(aspects),
        len(direct),
    )
    return (
        {child: set(parent_set) for child, parent_set in parents.items()},
        aspects,
        {child: set(parent_set) for child, parent_set in direct.items()},
    )


def build_rescue_maps(
    *,
    aspects: OntologyMap,
    direct_parents: ParentMap,
    terms_of_interest: set[str],
    t0_obsolete: ObsoleteInfo,
    t1_obsolete: ObsoleteInfo,
) -> tuple[dict[str, RescueEntry], dict[str, RescueEntry]]:
    future_obsolete: dict[str, RescueEntry] = {}
    t0_obsolete_rescue: dict[str, RescueEntry] = {}

    for go_id in aspects:
        if go_id in terms_of_interest:
            continue
        replacements = t1_obsolete.replaced_by.get(go_id)
        if replacements:
            future_obsolete[go_id] = RescueEntry(set(replacements), "future_obsolete_replaced_by")
            continue
        parents = direct_parents.get(go_id, set())
        if parents:
            future_obsolete[go_id] = RescueEntry(set(parents), "future_obsolete_parent")

    for go_id in t0_obsolete.obsolete:
        replacements = t0_obsolete.replaced_by.get(go_id)
        if replacements:
            t0_obsolete_rescue[go_id] = RescueEntry(set(replacements), "t0_obsolete_replaced_by")

    LOGGER.info(
        "Built rescue maps: future_obsolete=%d t0_obsolete_with_replaced_by=%d t0_obsolete_unrescued=%d",
        len(future_obsolete),
        len(t0_obsolete_rescue),
        len(t0_obsolete.obsolete - set(t0_obsolete_rescue)),
    )
    return future_obsolete, t0_obsolete_rescue


def build_eligible_term_ontology(
    *,
    aspects: OntologyMap,
    terms_of_interest: set[str],
    do_not_annotate: set[str],
    ontologies: tuple[str, ...],
) -> dict[str, str]:
    ontology_set = set(ontologies)
    eligible = {
        go_id: ontology
        for go_id, ontology in aspects.items()
        if ontology in ontology_set
        and go_id in terms_of_interest
        and go_id not in do_not_annotate
        and go_id not in BINDING_TERMS
        and go_id not in GO_ROOT_TERMS
    }
    LOGGER.info("Built eligible GO term lookup: %d terms", len(eligible))
    return eligible


def prefilter_rescue_map(
    rescue_map: dict[str, RescueEntry],
    *,
    eligible_term_ontology: dict[str, str],
    terms_of_interest: set[str],
    do_not_annotate: set[str],
) -> dict[str, FilteredRescueEntry]:
    filtered: dict[str, FilteredRescueEntry] = {}
    for go_id, entry in rescue_map.items():
        targets: list[tuple[str, str]] = []
        fixed_filtered = 0
        root_filtered = 0
        not_toi_filtered = 0
        do_not_filtered = 0
        unknown_ontology_filtered = 0

        for target_go_id in entry.targets:
            if target_go_id in BINDING_TERMS:
                fixed_filtered += 1
                continue
            if target_go_id in GO_ROOT_TERMS:
                root_filtered += 1
                continue
            if target_go_id not in terms_of_interest:
                not_toi_filtered += 1
                continue
            if target_go_id in do_not_annotate:
                do_not_filtered += 1
                continue
            ontology = eligible_term_ontology.get(target_go_id)
            if ontology is None:
                unknown_ontology_filtered += 1
                continue
            targets.append((target_go_id, ontology))

        filtered[go_id] = FilteredRescueEntry(
            targets=tuple(targets),
            reason=entry.reason,
            fixed_filtered=fixed_filtered,
            root_filtered=root_filtered,
            not_toi_filtered=not_toi_filtered,
            do_not_filtered=do_not_filtered,
            unknown_ontology_filtered=unknown_ontology_filtered,
        )

    n_with_targets = sum(1 for entry in filtered.values() if entry.targets)
    LOGGER.info(
        "Pre-filtered rescue map: source_terms=%d source_terms_with_targets=%d",
        len(filtered),
        n_with_targets,
    )
    return filtered


def insert_prediction(
    predictions: dict[str, defaultdict[str, dict[str, float]]],
    direct_predictions: dict[str, defaultdict[str, set[str]]],
    *,
    ontology: str,
    protein: str,
    go_id: str,
    score: float,
    is_direct: bool,
) -> str:
    protein_predictions = predictions[ontology][protein]
    if go_id not in protein_predictions:
        protein_predictions[go_id] = score
        if is_direct:
            direct_predictions[ontology][protein].add(go_id)
        return "stored"

    protein_direct_predictions = direct_predictions[ontology][protein]
    if go_id in protein_direct_predictions:
        return "ignored_existing_direct"

    if is_direct:
        protein_predictions[go_id] = score
        protein_direct_predictions.add(go_id)
        return "direct_replaced_rescued"

    return "ignored_existing_rescued"


def parse_scored_predictions_with_rescue(
    path: str | Path,
    *,
    input_format: str,
    protein_col: int,
    go_col: int,
    score_col: int,
    existing: AnnotationMap,
    do_not_annotate: set[str],
    terms_of_interest: set[str],
    reference_proteins: set[str],
    reference_proteins_by_ontology: dict[str, set[str]],
    eligible_term_ontology: dict[str, str],
    ontologies: tuple[str, ...],
    future_obsolete_rescue: dict[str, FilteredRescueEntry],
    t0_obsolete_rescue: dict[str, FilteredRescueEntry],
    t0_obsolete_terms: set[str],
) -> sweep.ScoredPredictionByOntology:
    prediction_path = Path(path)
    delimiter = sweep.delimiter_from_format(prediction_path, input_format)
    predictions: dict[str, defaultdict[str, dict[str, float]]] = {
        ontology: defaultdict(dict) for ontology in ontologies
    }
    direct_predictions: dict[str, defaultdict[str, set[str]]] = {
        ontology: defaultdict(set) for ontology in ontologies
    }

    total_rows = 0
    direct_candidates = 0
    rescued_candidates = 0
    stored = 0
    direct_replaced_rescued = 0
    duplicate_direct_ignored = 0
    duplicate_rescue_ignored = 0
    skipped_fixed_filter = 0
    skipped_root = 0
    skipped_no_reference_protein = 0
    skipped_no_reference_ontology = 0
    skipped_not_terms_of_interest = 0
    skipped_unknown_ontology = 0
    skipped_rescue_target_root = 0
    skipped_rescue_target_not_toi = 0
    skipped_rescue_target_do_not_annotate = 0
    skipped_rescue_target_unknown_ontology = 0
    rescue_future_replaced_by = 0
    rescue_future_parent = 0
    rescue_t0_replaced_by = 0
    rescue_t0_unrescued = 0
    rescue_no_surviving_targets = 0

    LOGGER.info("Parsing scored predictions with rescue from %s", prediction_path)

    with prediction_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header_checked = False
        for row_number, row in enumerate(reader, start=1):
            if not row or all(not field.strip() for field in row):
                continue
            if row[0].lstrip().startswith("#"):
                continue
            if not header_checked:
                header_checked = True
                if sweep.looks_like_header(row, protein_col, go_col, score_col):
                    continue
            if len(row) <= max(protein_col, go_col, score_col):
                raise ValueError(f"{prediction_path}: row {row_number} has too few columns")

            total_rows += 1
            protein = row[protein_col].strip()
            raw_go_id = row[go_col].strip()
            if not protein or not raw_go_id:
                skipped_fixed_filter += 1
                continue
            if protein not in reference_proteins:
                skipped_no_reference_protein += 1
                continue
            if raw_go_id in BINDING_TERMS:
                skipped_fixed_filter += 1
                continue
            if raw_go_id in GO_ROOT_TERMS:
                skipped_root += 1
                continue
            if raw_go_id in existing.get(protein, set()):
                skipped_fixed_filter += 1
                continue
            if raw_go_id in do_not_annotate:
                skipped_fixed_filter += 1
                continue

            is_direct = True
            entry: FilteredRescueEntry | None = None
            if raw_go_id in t0_obsolete_rescue:
                entry = t0_obsolete_rescue[raw_go_id]
                is_direct = False
                rescue_t0_replaced_by += 1
            elif raw_go_id in t0_obsolete_terms:
                rescue_t0_unrescued += 1
                continue
            elif raw_go_id in future_obsolete_rescue:
                entry = future_obsolete_rescue[raw_go_id]
                is_direct = False
                if entry.reason == "future_obsolete_replaced_by":
                    rescue_future_replaced_by += 1
                else:
                    rescue_future_parent += 1

            if entry is None:
                ontology = eligible_term_ontology.get(raw_go_id)
                if ontology is None:
                    if raw_go_id not in terms_of_interest:
                        skipped_not_terms_of_interest += 1
                    else:
                        skipped_unknown_ontology += 1
                    continue
                if protein not in reference_proteins_by_ontology[ontology]:
                    skipped_no_reference_ontology += 1
                    continue
                try:
                    score = float(row[score_col])
                except ValueError as exc:
                    raise ValueError(f"{prediction_path}: row {row_number} has non-numeric score {row[score_col]!r}") from exc
                direct_candidates += 1
                result = insert_prediction(
                    predictions,
                    direct_predictions,
                    ontology=ontology,
                    protein=protein,
                    go_id=raw_go_id,
                    score=score,
                    is_direct=True,
                )
                if result == "stored":
                    stored += 1
                elif result == "direct_replaced_rescued":
                    direct_replaced_rescued += 1
                elif result == "ignored_existing_direct":
                    duplicate_direct_ignored += 1
                else:
                    duplicate_rescue_ignored += 1
                continue

            skipped_fixed_filter += entry.fixed_filtered
            skipped_rescue_target_root += entry.root_filtered
            skipped_rescue_target_not_toi += entry.not_toi_filtered
            skipped_rescue_target_do_not_annotate += entry.do_not_filtered
            skipped_rescue_target_unknown_ontology += entry.unknown_ontology_filtered

            filtered_candidates: list[tuple[str, str]] = []
            for candidate_go_id, ontology in entry.targets:
                if protein not in reference_proteins_by_ontology[ontology]:
                    skipped_no_reference_ontology += 1
                    continue
                if candidate_go_id in existing.get(protein, set()):
                    skipped_fixed_filter += 1
                    continue
                filtered_candidates.append((candidate_go_id, ontology))

            if not filtered_candidates:
                if not is_direct:
                    rescue_no_surviving_targets += 1
                continue

            try:
                score = float(row[score_col])
            except ValueError as exc:
                raise ValueError(f"{prediction_path}: row {row_number} has non-numeric score {row[score_col]!r}") from exc

            for candidate_go_id, ontology in filtered_candidates:
                rescued_candidates += 1
                result = insert_prediction(
                    predictions,
                    direct_predictions,
                    ontology=ontology,
                    protein=protein,
                    go_id=candidate_go_id,
                    score=score,
                    is_direct=False,
                )
                if result == "stored":
                    stored += 1
                elif result == "direct_replaced_rescued":
                    direct_replaced_rescued += 1
                elif result == "ignored_existing_direct":
                    duplicate_direct_ignored += 1
                else:
                    duplicate_rescue_ignored += 1
                continue

    parsed = {ontology: dict(ontology_predictions) for ontology, ontology_predictions in predictions.items()}

    LOGGER.debug(
        "Parsed %s with rescue: rows=%d stored=%d direct_candidates=%d rescued_candidates=%d "
        "direct_replaced_rescued=%d duplicate_direct_ignored=%d duplicate_rescue_ignored=%d",
        prediction_path.name,
        total_rows,
        stored,
        direct_candidates,
        rescued_candidates,
        direct_replaced_rescued,
        duplicate_direct_ignored,
        duplicate_rescue_ignored,
    )
    LOGGER.debug(
        "Prediction rescue counts for %s: future_replaced_by=%d future_parent=%d "
        "t0_replaced_by=%d t0_unrescued=%d no_surviving_targets=%d",
        prediction_path.name,
        rescue_future_replaced_by,
        rescue_future_parent,
        rescue_t0_replaced_by,
        rescue_t0_unrescued,
        rescue_no_surviving_targets,
    )
    LOGGER.info(
        "Prediction filter counts for %s: fixed_filter_skipped=%d no_reference_protein_skipped=%d "
        "root_skipped=%d no_reference_ontology_skipped=%d direct_not_toi=%d "
        "direct_unknown_ontology=%d rescued_target_root=%d rescued_target_not_toi=%d "
        "rescued_target_do_not_annotate=%d rescued_target_unknown_ontology=%d",
        prediction_path.name,
        skipped_fixed_filter,
        skipped_no_reference_protein,
        skipped_root,
        skipped_no_reference_ontology,
        skipped_not_terms_of_interest,
        skipped_unknown_ontology,
        skipped_rescue_target_root,
        skipped_rescue_target_not_toi,
        skipped_rescue_target_do_not_annotate,
        skipped_rescue_target_unknown_ontology,
    )
    for ontology in ontologies:
        n_proteins = len(parsed[ontology])
        n_terms = sum(len(terms) for terms in parsed[ontology].values())
        LOGGER.info("Prediction terms for %s/%s: proteins=%d terms=%d", prediction_path.name, ontology, n_proteins, n_terms)
    return parsed


def evaluate_prediction_file(
    prediction_file: Path,
    *,
    args: argparse.Namespace,
    existing: AnnotationMap,
    do_not_annotate: set[str],
    terms_of_interest: set[str],
    reference_by_ontology: dict[str, AnnotationMap],
    reference_proteins: set[str],
    reference_proteins_by_ontology: dict[str, set[str]],
    go_parents: ParentMap,
    eligible_term_ontology: dict[str, str],
    ontologies: tuple[str, ...],
    taus: list[Decimal],
    future_obsolete_rescue: dict[str, FilteredRescueEntry],
    t0_obsolete_rescue: dict[str, FilteredRescueEntry],
    t0_obsolete_terms: set[str],
) -> list[dict[str, object]]:
    start_time = time.perf_counter()
    predictions_by_ontology = parse_scored_predictions_with_rescue(
        prediction_file,
        input_format=args.predicted_format,
        protein_col=args.protein_col,
        go_col=args.go_col,
        score_col=args.score_col,
        existing=existing,
        do_not_annotate=do_not_annotate,
        terms_of_interest=terms_of_interest,
        reference_proteins=reference_proteins,
        reference_proteins_by_ontology=reference_proteins_by_ontology,
        eligible_term_ontology=eligible_term_ontology,
        ontologies=ontologies,
        future_obsolete_rescue=future_obsolete_rescue,
        t0_obsolete_rescue=t0_obsolete_rescue,
        t0_obsolete_terms=t0_obsolete_terms,
    )
    active_ontologies = sweep.nonempty_prediction_ontologies(predictions_by_ontology, ontologies, prediction_file.name)
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
        with ProcessPoolExecutor(max_workers=worker_count, initializer=sweep.init_worker, initargs=(context,)) as executor:
            row_groups = list(executor.map(sweep.worker_sweep_ontology, active_ontologies))
        rows = [row for group in row_groups for row in group]
    else:
        sweep.init_worker(context)
        LOGGER.debug(
            "Running %d incremental ontology sweeps for %s serially",
            len(active_ontologies),
            prediction_file.name,
        )
        rows = []
        for ontology in active_ontologies:
            rows.extend(sweep.worker_sweep_ontology(ontology))
        sweep._WORKER_CONTEXT.clear()

    tau_order = {f"{float(tau):.4f}": index for index, tau in enumerate(taus)}
    ontology_order = {ontology: index for index, ontology in enumerate(ontologies)}
    rows.sort(key=lambda row: (tau_order[str(row["tau"])], ontology_order[str(row["ns"])]))
    if args.map_output:
        selected_rows = sweep.lowest_threshold_rows(rows)
        sweep.write_mapping_output(
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
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    go_parents, aspects, direct_parents = parse_go_parent_lookups(args.go_parent)

    do_not_annotate = parse_do_not_annotate(args.do_not_annotate)
    terms_of_interest = sweep.parse_terms_of_interest(args.terms_of_interest)
    if terms_of_interest is None:
        raise ValueError("--terms-of-interest is required for rescue mode")
    existing = parse_existing_annotations(args.existing)
    t0_obsolete = parse_obo_obsolete_info(args.t0_obo)
    t1_obsolete = parse_obo_obsolete_info(args.t1_obo)
    future_obsolete_rescue, t0_obsolete_rescue = build_rescue_maps(
        aspects=aspects,
        direct_parents=direct_parents,
        terms_of_interest=terms_of_interest,
        t0_obsolete=t0_obsolete,
        t1_obsolete=t1_obsolete,
    )

    raw_reference = sweep.parse_reference_annotations_for_sweep(
        args.reference,
        existing,
        do_not_annotate,
        terms_of_interest,
    )
    reference = remove_redundancy_fast(raw_reference, go_parents, "reference")

    ontologies = tuple(args.ontology) if args.ontology else sweep.ONTOLOGIES
    reference_by_ontology = {
        ontology: sweep.filter_annotations_by_ontology(reference, aspects, ontology)
        for ontology in ontologies
    }
    reference_proteins_by_ontology = {
        ontology: set(reference_by_ontology[ontology])
        for ontology in ontologies
    }
    reference_proteins = set().union(*reference_proteins_by_ontology.values())
    eligible_term_ontology = build_eligible_term_ontology(
        aspects=aspects,
        terms_of_interest=terms_of_interest,
        do_not_annotate=do_not_annotate,
        ontologies=ontologies,
    )
    future_obsolete_rescue_filtered = prefilter_rescue_map(
        future_obsolete_rescue,
        eligible_term_ontology=eligible_term_ontology,
        terms_of_interest=terms_of_interest,
        do_not_annotate=do_not_annotate,
    )
    t0_obsolete_rescue_filtered = prefilter_rescue_map(
        t0_obsolete_rescue,
        eligible_term_ontology=eligible_term_ontology,
        terms_of_interest=terms_of_interest,
        do_not_annotate=do_not_annotate,
    )
    taus = sweep.threshold_values(args.tau_min, args.tau_max, args.tau_step)
    LOGGER.debug("Evaluating %d thresholds from %s to %s with step %s", len(taus), args.tau_min, args.tau_max, args.tau_step)

    rows: list[dict[str, object]] = []
    for prediction_file in sweep.prediction_files(args):
        LOGGER.debug("Evaluating prediction file: %s", prediction_file)
        rows.extend(
            evaluate_prediction_file(
                prediction_file,
                args=args,
                existing=existing,
                do_not_annotate=do_not_annotate,
                terms_of_interest=terms_of_interest,
                reference_by_ontology=reference_by_ontology,
                reference_proteins=reference_proteins,
                reference_proteins_by_ontology=reference_proteins_by_ontology,
                go_parents=go_parents,
                eligible_term_ontology=eligible_term_ontology,
                ontologies=ontologies,
                taus=taus,
                future_obsolete_rescue=future_obsolete_rescue_filtered,
                t0_obsolete_rescue=t0_obsolete_rescue_filtered,
                t0_obsolete_terms=t0_obsolete.obsolete,
            )
        )

    all_output = out_dir / "evaluation_all.tsv"
    best_micro_output = out_dir / "evaluation_best_f_micro.tsv"
    best_macro_output = out_dir / "evaluation_best_f_macro.tsv"

    sweep.write_rows(all_output, rows)
    LOGGER.debug("Wrote full sweep output: %s (%d rows)", all_output, len(rows))

    best_micro_rows = sweep.best_rows_by_score(rows, "f_micro")
    sweep.write_rows(best_micro_output, best_micro_rows)
    LOGGER.debug("Wrote best micro output: %s (%d rows)", best_micro_output, len(best_micro_rows))

    best_macro_rows = sweep.best_rows_by_score(rows, "f_macro")
    sweep.write_rows(best_macro_output, best_macro_rows)
    LOGGER.debug("Wrote best macro output: %s (%d rows)", best_macro_output, len(best_macro_rows))
    LOGGER.debug("Finished rescue sweep in %.1f seconds", time.perf_counter() - start_time)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
