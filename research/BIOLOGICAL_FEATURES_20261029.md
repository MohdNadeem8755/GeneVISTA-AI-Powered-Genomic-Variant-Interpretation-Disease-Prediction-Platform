# Biological feature enrichment

Implementation: `src/biological_features.py`. Built artifact:
`data/processed/biological_features_20261029/features.sqlite`. Coverage and source
hashes are in the adjacent `report.json`; feature profiles are in `profiles.json`.

This is feature preparation, not retraining or deployment. No new accuracy is
claimed. The cached assays were independently measured, but their evidence may
also have contributed to ClinVar labels; independent generation does not establish
an independent evaluation benchmark.

## Features added or prepared

Completed build: **2,003,946 development rows**, comprising 1,633,384 training
and 370,562 validation records, with 35 enrichment fields. The default combined
training frame has 51 features. A 3,645-row validation loader check passed; all
five feature tests passed.

| Development coverage | Training | Validation |
| --- | ---: | ---: |
| ClinGen gene context | 751,819 | 170,436 |
| BRCA1 functional score | 464 | 135 |
| BRCA1 RNA score | 443 | 132 |
| TP53 SpliceAI | 502 | 118 |
| Matched functional assay (either gene) | 965 | 253 |

The build additionally counted feature availability on 247,503 calibration and
222,236 test records without reading labels or storing their feature rows. These
counts cover the full matched partitions, not the smaller previously reported
test extension. The report's per-score coverage counts are the quantitative
part of the overlap audit; actual upstream training intersections remain unknown.

| Evidence | Features and safeguards |
| --- | --- |
| BRCA1 MaveDB assay | Raw functional and RNA scores with verified transcript/coding-change matching. Negative values retained. No derived clinical criteria imported. |
| TP53 published assays | Kato, Kotler, Giacomelli, Kawaguchi and Funk assay fields where matched by the existing allele/protein checker. Missing interpretations remain missing. These are published assay interpretations, not new raw measurements. |
| Splicing | Published TP53 SpliceAI score, only for an exact supported allele/transcript match; broader genes remain unavailable. |
| Population frequency | Existing verified gnomAD v4.1 joint, population-maximum and South Asian AF, allele number, missingness flags and bounded negative-log frequency transforms. Missing is distinct from observed zero. No new population release was downloaded. |
| Inheritance | ClinGen gene-level AD/AR/XL/XLD/XLR/MT context indicators and other-mode count. Not a patient inheritance assignment. |
| Gene–disease evidence | Distinct disease count and counts by validity category; complete disease/inheritance/validity relations retained in a separate table. Multi-gene strings are not assigned to a single gene. |

ClinGen joins respect both classification date and source-file release date at
the requested cutoff. All other source availability dates have not been certified
for historical prediction times. This output is a static research snapshot;
do not use it in prospective temporal experiments without versioned as-of sources.

No clinical target, review status or PS3/BS3/PP3/BP4 code is selected by the builder.
Only training and validation features are stored. Calibration/test **feature
availability** is counted without reading labels; it is not used to select a model.

## Pretrained-score overlap audit

Exact per-variant upstream training/evaluation manifests are absent locally.
Therefore intersection sizes are **unknown**, not zero. The report quantifies
where each existing predictor is present in each partition; coverage does not
establish overlap. No predictor receives a certified-independent status.

| Predictor | Evidence reviewed | Decision |
| --- | --- | --- |
| AlphaMissense | Fine-tuned on human/primate population frequencies, with clinical/experimental benchmarking; not explicitly trained on those clinical labels. | Population/evaluation dependence remains possible; exclude from primary baseline, retain for sensitivity comparison. |
| CADD | Proxy-neutral human-derived versus simulated variants; model versions can incorporate other pretrained scores. | Exact component/version and validation overlap unresolved; sensitivity comparison only. |
| REVEL | Original study used HGMD disease variants and rare population variants plus multiple component predictors. | Shared clinical variants/component training are plausible; no exact exclusion possible without manifests. |
| SIFT converted rank score | Exact upstream version/training entities and ranking-reference overlap not verified from this cache. | Unresolved; excluded from primary profile. |
| PolyPhen-2 HVAR converted rank score | Exact upstream training manifest unavailable in this cache. | Unresolved; excluded from primary profile. |
| TP53 SpliceAI | Sequence-based splice prediction; exact pretrained-model overlap and source-table ascertainment not certified. | Retain explicitly as a predicted splice feature, and repeat comparisons without it. |

Sources: [AlphaMissense authors](https://deepmind.google/research/publications/21083/),
[CADD training documentation](https://cadd.bihealth.org/training),
[REVEL original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5065685/),
[SpliceAI repository](https://github.com/Illumina/SpliceAI),
[BRCA1 original assay](https://pmc.ncbi.nlm.nih.gov/articles/PMC6181777/).

To complete a strict overlap audit, obtain exact versioned upstream variant lists,
normalize them to GRCh38 alleles, intersect with validation/test keys and report
both overlap-excluded and full-cohort results. Protein-only identifiers additionally
require gene/transcript/protein mapping. Do not infer absence of overlap from
publication dates or unsupervised pretraining alone.

## Training adapter

```python
from src.biological_features import load_development

X, y, categorical_columns = load_development(
    "data/processed/annotation_v3/full_cohort.sqlite",
    "data/processed/biological_features_20261029/features.sqlite",
    "train", modulus=10, include_pretrained=False,
)
```

The default profile combines population, conservation and sequence features with
the new biological context. It excludes the five existing pathogenicity predictors.
The loader preserves existing partitions and checks annotation-source MD5. Fit
imputation/encoding on training folds only. Normalize missing categorical values
to an explicit token before CatBoost; preserve numerical NaNs. Use the same
sampling modulus for validation when running a pilot. `modulus=1` loads all
development records and needs more memory. Test label loading is rejected.

Required comparisons before deployment: baseline versus enriched; with/without
assays; with/without ClinGen context; with/without SpliceAI; then a separate
pretrained-predictor sensitivity comparison. Report gene-held-out and source-
independent performance, coverage, minority-class recall and balanced accuracy.
Assay coverage is restricted to two genes and cannot fix genome-wide performance
by itself. No blanket 95% guarantee follows from adding these features.

## Reproduce

```powershell
.venv\Scripts\python.exe -m src.biological_features --output data/processed/biological_features_NEW --cutoff 2026-10-29
.venv\Scripts\python.exe -m unittest tests.test_biological_features -v
```

Existing output directories are rejected. Original datasets and active models are
not modified. Keep source hashes and the availability report with any trained model.
