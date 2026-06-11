#!/usr/bin/env python3
"""Perl-equivalent Python implementation of FPR_calculations.pl."""

from fpr_metric_common import evaluate_pairwise_mapping, remove_redundancy_pairwise, run_metric


if __name__ == "__main__":
    run_metric(remove_redundancy=remove_redundancy_pairwise, map_annotations=evaluate_pairwise_mapping)
