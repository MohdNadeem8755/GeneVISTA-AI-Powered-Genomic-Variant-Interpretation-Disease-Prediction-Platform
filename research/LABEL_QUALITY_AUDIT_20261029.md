# Label-quality audit — 2026-10-29

## Scope

Audited all **2,962,171 training/validation records** in the active annotation
cohort against the latest locally cached ClinVar variant summary (October 5,
2026). The original cohort uses September 2026 data. Scanned all 9,290,347
source rows, retaining only development IDs with supported GRCh38 germline/unknown
SNV coordinates. Comparisons require both VariationID and exact allele identity.
No labels, deployed models, or source databases were changed. No individual test
or calibration labels were read; existing cohort metadata includes historical
aggregate partition counts. The October snapshot is latest locally, not a claim
about the newest available online release.

## Findings

| Check | Count |
| --- | ---: |
| Exact matches in newer snapshot | 2,961,997 |
| Changed labels | 27,308 |
| Now marked conflicting | 11,737 |
| Last evaluated more than five years ago | 426,931 |
| Missing/unparseable evaluation dates | 12,235 |
| No exact newer match | 174 |
| Duplicate allele groups remaining in development | 0 |
| Duplicate VariationIDs remaining in development | 0 |
| Conflict-flagged records remaining in development | 0 |
| Conflict-flagged records already excluded by cohort construction | 1,500 |

Counts overlap and must not be summed as independent errors. The cohort conflict
flag combines differing labels and differing IDs. The zero duplicate count applies
to the already consolidated cohort; it does not establish that raw submissions
are independent. No duplicate exact observations were found in the newer-source
matches. Old evaluation dates indicate review priority, not proof of incorrect
labels; missing matches do not prove withdrawal.

The largest transitions are Benign → Benign/Likely benign (11,558), VUS →
conflicting (7,559), Likely benign → conflicting (2,564), and Pathogenic →
Pathogenic/Likely pathogenic (2,117). Combined labels are not necessarily
biological contradictions and must not be silently mapped to one five-class label.

## Review quality and expert-only evaluation

Among 2,003,946 annotation-matched development records:

| Review tier | Records |
| --- | ---: |
| Single submitter with criteria | 1,678,820 |
| Multiple submitters, no conflicts | 273,326 |
| Expert panel/practice guideline | 8,567 |
| Other | 43,233 |

The unchanged active model was evaluated separately on **1,558 expert-reviewed
validation records**, covering **0.4204%** of annotation-matched validation.
This uses the original snapshot's review status and labels. It is a reused
validation subset, not independent clinical validation or a new test result.

| Metric | Value |
| --- | ---: |
| Accuracy | 57.57% |
| Balanced accuracy / macro recall | 55.95% |
| Macro precision | 60.09% |
| Macro F1 | 0.5706 |

All five classes are represented: Benign 122, Likely benign 95, Likely pathogenic
336, Pathogenic 534, and VUS 471. This small, differently distributed subset is
not directly comparable to overall test accuracy. Filtering the evaluation to
expert-reviewed records does not produce 95% and sharply limits coverage.

## Recommended next dataset version

1. Review newly conflicting records first and quarantine them from exact five-class
   training until a documented resolution exists.
2. Reconcile changed labels against disease context and source evidence. Preserve
   old/new labels and snapshot dates. Do not treat all changes as old errors.
3. Keep combined labels separate; do not force them into an exact class.
4. Investigate 174 missing exact matches and missing dates before making changes.
5. Prioritize older evaluations for review without deleting them automatically.
6. Build a versioned high-review training cohort and retain a broader cohort for
   coverage analysis. Preserve locus partitions; reserve a genuinely later test
   before further tuning. Do not selectively remove difficult test examples.

These are proposed follow-up changes, not changes performed by this audit.

## Reproduce and inspect

```powershell
$env:OMP_NUM_THREADS='4'
.venv\Scripts\python.exe -m src.audit_label_quality --output reports/label_quality_NEW
```

Full counts, per-class metrics, confusion matrix, source SHA-256 and provenance:
`reports/label_quality_20261029/audit.json`.

Row-level review queue with reasons:
`reports/label_quality_20261029/review_queue.csv`.

Matching-stage database (development records only):
`reports/label_quality_20261029/latest_development.sqlite`.

The script rejects existing output directories. Inputs are opened read-only.
