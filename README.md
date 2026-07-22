# Alternative Metrics to Evaluate Functional Annotations

> This pipeline is adapted from the process that generated **Supplementary Table 6** of the PAN-GO annotation in the PAN-GO Nature paper:
> 📄 [https://www.nature.com/articles/s41586-025-08592-0](https://www.nature.com/articles/s41586-025-08592-0)

---

## Overview

The pipeline consists of three main steps:

1. **Generate a GO parent lookup file**
2. **Gather files required for the calculation**
3. **Compare GO annotations** between prediction and reference datasets, and calculate precision, recall, and F score

---

## Step 1: Generate the GO Parent File

This step creates a lookup file containing all parent GO terms.

### Instructions

**1. Download the GO ontology file:**

```
https://current.geneontology.org/ontology/go-basic.obo
```

**2. Run the script:**

```bash
perl findGOparent_OBO.pl -i go-basic.obo > goparents
```

### Notes

- Only `is_a` and `part_of` relationships are used
- All relationships remain within the same GO aspect

### Output Format — `goparents`

| Column | Description |
|--------|-------------|
| 1 | Child term |
| 2 | Parent term |
| 3 | Relationship (`is_a`, `part_of`, or blank for indirect) |
| 4 | GO aspect |

### Example

```
protein serine kinase activity(GO:0106310)    protein kinase activity(GO:0004672)    is_a    molecular_function
protein serine kinase activity(GO:0106310)    catalytic activity(GO:0003824)                 molecular_function
```

---

## Step 2: Gather Required Files

All annotation files are **tab-delimited**, with the protein identifier in column 1 and the GO identifier in column 2.

Prepare the following four files:

| File | Description |
|------|-------------|
| `predicted_annotations` | The annotations to be evaluated |
| `existing_annotations` | Experimental annotations on the date the predictions were made. Should include all parent terms via `is_a` or `part_of` within the same ontology aspect |
| `new_annotations` | Experimental annotations created **after** the prediction date. Should contain only the genomes present in the predictions |
| `do_not_annotate` list | Terms labelled `gocheck_do_not_annotate` in `go-basic.obo`. These are excluded from the calculation |


*Test files are provided in the Files/ folder.

---

## Step 3: Generate Metrics

### Script

```
FPR_calculations.pl
```

The script calculates **precision**, **recall**, and **F score**. Here is how it works:

---

### 3a. Prepare Predicted Annotations

1. Take the `predicted_annotations` file (tab-delimited). Please note that this program **doesn't take score into consideration**. If you want to use a score cutoff, please preprocess the file to include only the annotations that meet the cutoff.
2. Remove annotations matching the `existing_annotations` file
3. Remove annotations to binding (`GO:0005488`) and protein binding (`GO:0005515`)
4. Remove root terms: `GO:0008150`, `GO:0005575`, `GO:0003674`
5. Remove terms in the `do_not_annotate` list
6. Remove parent terms where more specific predicted terms exist for the same protein (non-redundant list)

---

### 3b. Prepare Reference Annotations

1. Take the `new_annotations` file
2. Remove annotations already present in `existing_annotations` (e.g., re-dated entries with no annotation change)
3. Remove annotations to binding (`GO:0005488`) and protein binding (`GO:0005515`)
4. Remove terms in the `do_not_annotate` list
5. Remove parent terms where more specific terms are present

---

### 3c. Mapping Types

For each protein with at least one predicted term after filtering, check if it has new experimental annotations and assign a mapping type:

| Type | Code | Description |
|------|------|-------------|
| `direct` | E | Predicted and reference GO terms are identical |
| `true` | L | Predicted term is **less specific** than the reference term (predicted is a parent of the reference) |
| `related` | M | Predicted term is **more specific** than the reference term (predicted is a child of the reference) |
| `unrelated` | U | Predicted and reference terms are not related |
| `no map` | — | Protein has no corresponding annotation in the reference file |

> **Note on `true` (L) counting:** The count reflects the number of most specific reference terms, not the number of predicted terms. For example, if 4 predicted terms all map to one reference child term, the count is 1.

Also, for each protein in `new_annotations`, count reference terms with no mapping to any predicted term (**Z**). This is why `new_annotations` should be restricted to predicted genomes only.

---

### 3d. Score Formulas

```
Precision (per protein) = [E + (0.75 × L) + (0.5 × M)] / [E + L + M + U]

Recall (per protein)    = [E + (0.75 × L) + (0.5 × M)] / [E + L + M + Z]

```

- **Average Precision** — averaged over all proteins with at least one predicted GO term **and** at least one experimental GO term
- **Average Recall** — averaged over all proteins with at least one experimental GO term
- **F Score** = (2 × Avg_Precision × Avg_Recall) / (Avg_Precision + Avg_Recall)
---

## Usage

```bash
perl scripts/FPR_calculations.pl \
  -p test_prediction \
  -t test_existing_annotations \
  -r test_new_annotations \
  -n GO_do_not_annotate_list \
  -g goparents \
  -o test.map > test.FPR
```

### Arguments

| Flag | Description |
|------|-------------|
| `-p` | Prediction file (tab-delimited) |
| `-t` | Existing annotation file (experimental annotations at time of prediction) |
| `-r` | New annotations made after predictions |
| `-n` | `GO_do_not_annotate` list |
| `-g` | `goparents` file |
| `-o` | Output mapping file (per protein) |
| `STDOUT` | Precision, recall, and F score summary |
| `STDERR` | Log file |

---

## Python Wrapper

The repository also includes `fpr_wrapper.py`, which runs both core Perl steps:
it creates `goparents` and the `gocheck_do_not_annotate` term list from a GO
OBO file, normalizes prediction input, and then calls `FPR_calculations.pl`.

The metric script already removes redundant parent terms from both predicted
and reference annotations before mapping, so the wrapper does not duplicate
that filtering logic.

Two Python metric implementations are also available and can be passed with
`--metric-perl-script`:

| Script | Behavior |
|--------|----------|
| `fpr_metric_python_equiv.py` | Python port of the current `FPR_calculations.pl` logic, preserving the same filtering and mapping behavior |
| `fpr_metric_python_fast.py` | Optimized Python implementation of the same metric using faster non-redundancy filtering and per-protein mapping |

Prediction input can be TSV or CSV. Plain TSV files are passed directly to the
metric script without Python parsing. CSV files, and TSV files with a header row, are
streamed line-by-line into a temporary two-column TSV. The wrapper auto-detects
the delimiter by default, or you can force it with `--predicted-format tsv` or
`--predicted-format csv`.

### Evaluate One Prediction File

```bash
python fpr_wrapper.py \
  -p Files/test_predicted \
  -i go-basic.obo \
  -t Files/test_existing_annotations \
  -r Files/test_new_annotations \
  --goparents-output goparents \
  --do-not-annotate-output GO_do_no_annotate_list \
  --map-output test.map \
  --fpr-output test.FPR \
  --log-level INFO \
  -o summary.tsv
```

The wrapper first runs `findGOparents_OBO.pl` on `go-basic.obo`. If
`--goparents-output` is omitted, the generated lookup file is kept only in a
temporary directory for that run. It also extracts OBO terms marked
`subset: gocheck_do_not_annotate` into the same two-column format as
`Files/GO_do_no_annotate_list`. Use `--do-not-annotate-output` to keep that
generated file, or `--do-not-annotate` to provide a precomputed list. It then
normalizes the prediction file to the TSV format expected by the metric script
and calls the selected Perl or Python metric implementation.

`summary.tsv` contains one row with `filename`, `average_precision`,
`average_recall`, and `F_score`. `--fpr-output` preserves the full metric STDOUT
report, including per-protein precision and recall.

Perl STDERR is suppressed by default to avoid very large per-annotation logs on
large runs. Use `--perl-stderr-log run.log` to keep it. The wrapper logs stage
progress and aggregate E/L/M/U/Z counts through Python's logging module; use
`--log-level DEBUG` for more detail.

To use the optimized Python metric script with the wrapper, add:

```bash
--metric-perl-script fpr_metric_python_fast.py
```

### Prepare GO Helper Artifacts Only

The wrapper can also generate `goparents` and the do-not-annotate list without
running a metric. This is useful before running `fpr_metric_sweep.py`, where
the ontology artifacts are supplied directly.

```bash
python fpr_wrapper.py \
  -i go-basic.obo \
  --goparents-output goparents \
  --do-not-annotate-output GO_do_no_annotate_list
```

In this artifact-only mode, prediction files, existing annotations, and new
annotations are not required. You can generate only one artifact by providing
only `--goparents-output` or only `--do-not-annotate-output`.

### Evaluate a Directory of Prediction Files

```bash
python fpr_wrapper.py \
  -d prediction_dir \
  -i go-basic.obo \
  -t Files/test_existing_annotations \
  -r Files/test_new_annotations \
  -o directory_summary.tsv
```

Directory mode recursively evaluates every file under `prediction_dir`.
It does not write per-protein mapping files or individual-protein summaries.
The output summary table has one row per input filename and aggregate metric
columns: `average_precision`, `average_recall`, and `F_score`.

Use `--glob` to restrict recursive input selection, for example
`--glob '*.csv'`.

---

## Threshold Sweep Metric

`fpr_metric_sweep.py` evaluates scored predictions over a user-defined
threshold grid and writes one row per prediction file, threshold, and GO
ontology. The output includes pooled micro metrics, protein-averaged macro
metrics, and E/L/M/U/Z counts. Generate `goparents` and
`GO_do_no_annotate_list` once with `fpr_wrapper.py` before running large
sweeps.

```bash
python fpr_metric_sweep.py \
  -p predictions.csv \
  -t Files/test_existing_annotations \
  -r Files/test_new_annotations \
  -n Files/GO_do_no_annotate_list \
  --terms-of-interest GO_terms_of_interest \
  -g goparents \
  --predicted-format csv \
  --tau-step 0.01 \
  --tau-min 0 \
  --tau-max 1 \
  --workers 4 \
  --log-level INFO \
  --best-macro-output hgrs_best_macro.tsv \
  --best-micro-output hgrs_best_micro.tsv \
  --map-output hgrs_lowest_tau.map \
  -o hgrs_sweep.tsv
```

Directory mode is also supported:

```bash
python fpr_metric_sweep.py \
  -d prediction_dir \
  --glob '*.csv' \
  -t Files/test_existing_annotations \
  -r Files/test_new_annotations \
  -n Files/GO_do_no_annotate_list \
  --terms-of-interest GO_terms_of_interest \
  -g goparents \
  --predicted-format csv \
  --tau-step 0.01 \
  --log-level INFO \
  --best-macro-output hgrs_best_macro.tsv \
  --best-micro-output hgrs_best_micro.tsv \
  -o hgrs_sweep.tsv
```

The full sweep file keeps one row per threshold and ontology. The optional
best-threshold files retain the same columns but keep only one row per input
file and ontology: `--best-macro-output` selects by highest `f_macro`, and
`--best-micro-output` selects by highest `f_micro`. Ties are resolved by the
higher threshold. Filenames are reported as basenames, and floating-point
values are written with three decimal places.

In single-file mode, `--map-output` writes a debug mapping file for the lowest
threshold rows after all filtering, thresholding, and redundancy removal. With
the default `--tau-min 0`, this maps all retained predictions regardless of
score. The mapping file includes `filename`, `ns`, `tau`, `id`, `predicted`,
`map type`, and `reference`.

The prediction file is expected to contain protein, GO term, and score columns
in columns 1, 2, and 3 by default. Use `--protein-col`, `--go-col`, and
`--score-col` for different 0-based column positions.

Use `--terms-of-interest` to provide a one-column, no-header file of eligible
GO IDs. When supplied, predictions and reference annotations outside that list
are excluded before redundancy filtering and threshold sweeping. Existing
annotations are not filtered by this list.

`pr_micro`, `rc_micro`, and `f_micro` are pooled-count micro-average metrics.
`pr_macro`, `rc_macro`, and `f_macro` average per-protein precision and recall
using the same denominators as the original metric, scoped to each ontology.

`--workers` uses process-based parallelism across ontology-level incremental
sweeps. More than 3 workers does not help when all three ontologies are
evaluated, because there are only three ontology sweep jobs. For very large
prediction files, process workers can exceed memory limits because each worker
needs access to parsed prediction state.
Use `--log-level INFO` to report parsing counts, per-ontology term counts, and
per-file runtime; use `--log-level WARNING` for quieter batch runs.

---

## Output Formats

### Mapping File (`-o`)

| Column | Description |
|--------|-------------|
| 1 | UniProt ID |
| 2 | Predicted GO term |
| 3 | Mapping type (`direct`, `true`, `related`, `unrelated`) |
| 4 | Mapped GO term from the new annotation file |

### FPR File (`STDOUT`)

```
Average Precision
Average Recall
F Score

UniProtID    precision    recall
...
```

---

## Reference

If you use this pipeline, please cite:

> PAN-GO: *Nature* (2025). [https://www.nature.com/articles/s41586-025-08592-0](https://www.nature.com/articles/s41586-025-08592-0)
