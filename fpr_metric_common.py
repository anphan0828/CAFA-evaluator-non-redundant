#!/usr/bin/env python3
"""Shared Python implementation for FPR GO annotation metrics."""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Callable, DefaultDict


GO_RE = re.compile(r"GO:\d+")
BINDING_TERMS = {"GO:0005488", "GO:0005515"}
GO_ROOT_TERMS = {"GO:0008150", "GO:0003674", "GO:0005575"}
LOG_ENABLED = True

AnnotationMap = dict[str, set[str]]
ParentMap = dict[str, set[str]]
OntologyMap = dict[str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Calculate FPR GO annotation metrics")
    parser.add_argument("-p", dest="prediction", required=True, help="predicted annotations TSV")
    parser.add_argument("-t", dest="existing", required=True, help="existing annotations TSV")
    parser.add_argument("-r", dest="reference", required=True, help="new/reference annotations TSV")
    parser.add_argument("-n", dest="do_not_annotate", required=True, help="GO do-not-annotate TSV")
    parser.add_argument("-g", dest="go_parent", required=True, help="goparents lookup file")
    parser.add_argument("-o", dest="out_file", required=True, help="mapping output TSV")
    parser.add_argument("-i", dest="in_file", help=argparse.SUPPRESS)
    parser.add_argument("-O", dest="unused_o", help=argparse.SUPPRESS)
    parser.add_argument("-e", dest="err_file", help=argparse.SUPPRESS)
    parser.add_argument("-v", dest="verbose", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("-V", dest="very_verbose", action="store_true", help=argparse.SUPPRESS)
    return parser.parse_args()


def log(message: str) -> None:
    if LOG_ENABLED:
        print(message, file=sys.stderr)


def set_metric_logging(enabled: bool) -> None:
    global LOG_ENABLED
    LOG_ENABLED = enabled


def extract_go_id(term: str) -> str | None:
    match = GO_RE.search(term)
    return match.group(0) if match else None


def parse_go_parents(path: str | Path) -> ParentMap:
    log("Working on GO child parent term lookup file.")
    parents: DefaultDict[str, set[str]] = defaultdict(set)
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            child_id = extract_go_id(fields[0])
            parent_id = extract_go_id(fields[1])
            if child_id and parent_id:
                parents[child_id].add(parent_id)
    return {child: set(parent_set) for child, parent_set in parents.items()}


def parse_go_parent_aspects(path: str | Path) -> OntologyMap:
    """Parse GO aspect labels from a goparents lookup file.

    The goparents file carries the GO aspect in column 4. The child and parent
    terms in each row are within the same aspect, so both can be assigned from
    the same row.
    """

    aspects: dict[str, str] = {}
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or not fields[3]:
                continue
            aspect = fields[3]
            child_id = extract_go_id(fields[0])
            parent_id = extract_go_id(fields[1])
            if child_id:
                aspects[child_id] = aspect
            if parent_id:
                aspects[parent_id] = aspect
    return aspects


def parse_do_not_annotate(path: str | Path) -> set[str]:
    log("Working on the do_not_annotate file.")
    terms: set[str] = set()
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0]:
                terms.add(fields[0])
    return terms


def parse_existing_annotations(path: str | Path) -> AnnotationMap:
    log("Working on the existing annotation file.")
    existing: DefaultDict[str, set[str]] = defaultdict(set)
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            protein, go_id = fields[0], fields[1]
            if go_id in BINDING_TERMS:
                continue
            existing[protein].add(go_id)
    return {protein: set(terms) for protein, terms in existing.items()}


def parse_reference_annotations(
    path: str | Path,
    existing: AnnotationMap,
    do_not_annotate: set[str],
) -> AnnotationMap:
    log("Working on the new annotation file (reference annotations).")
    new: DefaultDict[str, set[str]] = defaultdict(set)
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            protein, go_id = fields[0], fields[1]
            if go_id in existing.get(protein, set()):
                log(f"--{protein} {go_id} already in existing_annotation, removed from reference.")
                continue
            if go_id in BINDING_TERMS:
                continue
            if go_id in do_not_annotate:
                log(f"--{protein} {go_id} in the do_not_annotate list, removed from reference.")
                continue
            new[protein].add(go_id)
    return {protein: set(terms) for protein, terms in new.items()}


def parse_predicted_annotations(
    path: str | Path,
    existing: AnnotationMap,
    do_not_annotate: set[str],
) -> AnnotationMap:
    log("Working on the predicted annotation file.")
    predicted: DefaultDict[str, set[str]] = defaultdict(set)
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            protein, go_id = fields[0], fields[1]
            if go_id in BINDING_TERMS:
                continue
            # Keep parity with the current Perl script: its root-term filter
            # checks an undefined variable, so roots are not filtered here.
            if go_id in existing.get(protein, set()):
                log(f"--{protein} {go_id} already in existing_annotation, removed from predicted.")
                continue
            if go_id in do_not_annotate:
                log(f"--{protein} {go_id} in the do_not_annotate list, removed from predicted.")
                continue
            predicted[protein].add(go_id)
    return {protein: set(terms) for protein, terms in predicted.items()}


def remove_redundancy_pairwise(annotations: AnnotationMap, go_parents: ParentMap, label: str) -> AnnotationMap:
    """Perl-equivalent parent removal with pairwise per-protein comparisons."""

    result: AnnotationMap = {}
    for protein, terms in annotations.items():
        kept: set[str] = set()
        for go_id in terms:
            is_redundant = False
            for other in terms:
                if other == go_id:
                    continue
                if go_id in go_parents.get(other, set()):
                    is_redundant = True
                    log(f"--{protein} {go_id} is a parent GO term to {other}, removed from {label}.")
                    break
            if not is_redundant:
                kept.add(go_id)
        if kept:
            result[protein] = kept
    return result


def remove_redundancy_fast(annotations: AnnotationMap, go_parents: ParentMap, label: str) -> AnnotationMap:
    """Remove parent terms by intersecting each protein's terms with all ancestors."""

    result: AnnotationMap = {}
    for protein, terms in annotations.items():
        redundant: set[str] = set()
        for go_id in terms:
            redundant.update(go_parents.get(go_id, set()) & terms)
        kept = terms - redundant
        if kept:
            result[protein] = kept
        for go_id in redundant:
            log(f"--{protein} {go_id} is a parent GO term to a more specific term, removed from {label}.")
    return result


def evaluate_pairwise_mapping(predicted: AnnotationMap, reference: AnnotationMap, go_parents: ParentMap) -> tuple[dict, dict, dict]:
    """Perl-equivalent mapping using predicted x reference loops per protein."""

    log("Mapping annotations.")
    prediction_map: DefaultDict[str, DefaultDict[str, dict[str, str | int]]] = defaultdict(lambda: defaultdict(dict))
    prediction_map_count: DefaultDict[str, set[frozenset[str]]] = defaultdict(set)
    reference_map: DefaultDict[str, set[str]] = defaultdict(set)

    for protein, predicted_terms in predicted.items():
        if protein not in reference:
            continue
        for go_p in predicted_terms:
            direct: set[str] = set()
            true: set[str] = set()
            related: set[str] = set()
            for go_r in reference[protein]:
                if go_p == go_r:
                    direct.add(go_r)
                elif go_r in go_parents.get(go_p, set()):
                    related.add(go_r)
                elif go_p in go_parents.get(go_r, set()):
                    true.add(go_r)

            if direct:
                joined = ";".join(sorted(direct))
                prediction_map[protein]["direct"][go_p] = joined
                reference_map[protein].update(direct)
            elif true:
                joined = ";".join(sorted(true))
                prediction_map[protein]["true"][go_p] = joined
                prediction_map_count[protein].add(frozenset(true))
                reference_map[protein].update(true)
            elif related:
                joined = ";".join(sorted(related))
                prediction_map[protein]["related"][go_p] = joined
                reference_map[protein].update(related)
            else:
                prediction_map[protein]["unrelated"][go_p] = 1

    reference_nomap = build_reference_nomap(reference, reference_map)
    return prediction_map, prediction_map_count, reference_nomap


def evaluate_fast_mapping(predicted: AnnotationMap, reference: AnnotationMap, go_parents: ParentMap) -> tuple[dict, dict, dict]:
    """Mapping optimized to avoid predicted x reference Cartesian scans."""

    log("Mapping annotations.")
    prediction_map: DefaultDict[str, DefaultDict[str, dict[str, str | int]]] = defaultdict(lambda: defaultdict(dict))
    prediction_map_count: DefaultDict[str, set[frozenset[str]]] = defaultdict(set)
    reference_map: DefaultDict[str, set[str]] = defaultdict(set)

    for protein, predicted_terms in predicted.items():
        if protein not in reference:
            continue
        reference_terms = reference[protein]

        true_by_predicted: DefaultDict[str, set[str]] = defaultdict(set)
        for go_r in reference_terms:
            for parent in go_parents.get(go_r, set()):
                if parent in predicted_terms:
                    true_by_predicted[parent].add(go_r)

        for go_p in predicted_terms:
            if go_p in reference_terms:
                prediction_map[protein]["direct"][go_p] = go_p
                reference_map[protein].add(go_p)
                continue

            true = true_by_predicted.get(go_p, set())
            if true:
                joined = ";".join(sorted(true))
                prediction_map[protein]["true"][go_p] = joined
                prediction_map_count[protein].add(frozenset(true))
                reference_map[protein].update(true)
                continue

            related = go_parents.get(go_p, set()) & reference_terms
            if related:
                joined = ";".join(sorted(related))
                prediction_map[protein]["related"][go_p] = joined
                reference_map[protein].update(related)
            else:
                prediction_map[protein]["unrelated"][go_p] = 1

    reference_nomap = build_reference_nomap(reference, reference_map)
    return prediction_map, prediction_map_count, reference_nomap


def build_reference_nomap(reference: AnnotationMap, reference_map: dict[str, set[str]]) -> dict[str, set[str]]:
    reference_nomap: dict[str, set[str]] = {}
    for protein, terms in reference.items():
        unmapped = terms - reference_map.get(protein, set())
        if unmapped:
            reference_nomap[protein] = unmapped
    return reference_nomap


def compute_metrics(
    prediction_map: dict,
    prediction_map_count: dict,
    reference_nomap: dict[str, set[str]],
    reference: AnnotationMap,
) -> tuple[dict[str, float], dict[str, float], float, float, float]:
    log("Calculating precision, recall and F scores.")
    precision_by_protein: dict[str, float] = {}
    recall_by_protein: dict[str, float] = {}
    sum_precision = 0.0
    sum_recall = 0.0

    for protein in prediction_map:
        e_count = len(prediction_map[protein].get("direct", {}))
        l_count = len(prediction_map_count.get(protein, set()))
        m_count = len(prediction_map[protein].get("related", {}))
        unrelated_count = len(prediction_map[protein].get("unrelated", {}))
        z_count = len(reference_nomap.get(protein, set()))

        numerator = e_count + 0.75 * l_count + 0.5 * m_count
        precision_denominator = e_count + l_count + m_count + unrelated_count
        recall_denominator = e_count + l_count + m_count + z_count
        precision = numerator / precision_denominator if precision_denominator else 0.0
        recall = numerator / recall_denominator if recall_denominator else 0.0

        precision_by_protein[protein] = precision
        recall_by_protein[protein] = recall
        sum_precision += precision
        sum_recall += recall

    n_gene_mapped_predictions = len(prediction_map)
    n_gene_reference = len(reference)
    average_precision = sum_precision / n_gene_mapped_predictions if n_gene_mapped_predictions else 0.0
    average_recall = sum_recall / n_gene_reference if n_gene_reference else 0.0
    f_score = (
        2 * average_precision * average_recall / (average_precision + average_recall)
        if average_precision + average_recall
        else 0.0
    )
    return precision_by_protein, recall_by_protein, average_precision, average_recall, f_score


def write_mapping(path: str | Path, prediction_map: dict) -> None:
    with Path(path).open("w") as handle:
        handle.write("id\tpredicted\tmap type\treference\n")
        for protein in sorted(prediction_map):
            for map_type in sorted(prediction_map[protein]):
                for go_p in sorted(prediction_map[protein][map_type]):
                    go_r = prediction_map[protein][map_type][go_p]
                    go_r_text = str(go_r) if str(go_r).startswith("GO") else ""
                    handle.write(f"{protein}\t{go_p}\t{map_type}\t{go_r_text}\n")


def write_summary(
    precision_by_protein: dict[str, float],
    recall_by_protein: dict[str, float],
    average_precision: float,
    average_recall: float,
    f_score: float,
) -> None:
    print(f"Average precision\t{average_precision:.15g}")
    print(f"Average recall\t{average_recall:.15g}")
    print(f"F_score\t{f_score:.15g}")
    print("\n")
    print("Precision and Recall for individual proteins\n")
    print("id\tpredicion\trecall")
    for protein in sorted(precision_by_protein):
        print(f"{protein}\t{precision_by_protein[protein]:.15g}\t{recall_by_protein[protein]:.15g}")


def run_metric(
    *,
    remove_redundancy: Callable[[AnnotationMap, ParentMap, str], AnnotationMap],
    map_annotations: Callable[[AnnotationMap, AnnotationMap, ParentMap], tuple[dict, dict, dict]],
) -> None:
    args = parse_args()
    go_parents = parse_go_parents(args.go_parent)
    do_not_annotate = parse_do_not_annotate(args.do_not_annotate)
    existing = parse_existing_annotations(args.existing)
    new = parse_reference_annotations(args.reference, existing, do_not_annotate)
    reference = remove_redundancy(new, go_parents, "reference")
    raw_predicted = parse_predicted_annotations(args.prediction, existing, do_not_annotate)
    predicted = remove_redundancy(raw_predicted, go_parents, "predicted")
    prediction_map, prediction_map_count, reference_nomap = map_annotations(predicted, reference, go_parents)
    precision, recall, average_precision, average_recall, f_score = compute_metrics(
        prediction_map,
        prediction_map_count,
        reference_nomap,
        reference,
    )
    write_mapping(args.out_file, prediction_map)
    write_summary(precision, recall, average_precision, average_recall, f_score)
