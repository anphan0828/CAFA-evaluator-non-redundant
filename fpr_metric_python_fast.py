#!/usr/bin/env python3
"""Optimized Python implementation of the FPR GO annotation metric."""

from fpr_metric_common import evaluate_fast_mapping, remove_redundancy_fast, run_metric


if __name__ == "__main__":
    run_metric(remove_redundancy=remove_redundancy_fast, map_annotations=evaluate_fast_mapping)
