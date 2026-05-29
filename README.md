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
Precision (per protein) = [E + (0.75 × L) + (0.5 × M)] / [E + L + M + A]

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
