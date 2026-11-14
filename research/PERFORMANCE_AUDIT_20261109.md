# Performance audit and versioned corrections — November 9, 2026

Scope: explain existing performance, inspect data/feature/training/reporting
flaws, correct confirmed code issues and present class-level performance.
No retraining, model replacement or new test-driven hyperparameter selection
was performed in this audit. Existing partition scores remain unchanged.

## Observed performance, not inferred causation

Validation contains 64,224 VUS records among 73,770 (87.06%). The ensemble
predicts VUS for 57,030 records (77.31%). Thus this model is not simply predicting
VUS for everything. Class weighting trades some majority recall for minority
recall, but produces many false Likely benign / Likely pathogenic predictions.

| Actual class | Actual count | Predicted count | Precision | Recall |
| --- | ---: | ---: | ---: | ---: |
| Benign | 1,039 | 1,140 | 64.04% | 70.26% |
| Likely benign | 3,819 | 7,473 | 29.24% | 57.21% |
| Likely pathogenic | 2,194 | 5,996 | 20.46% | 55.93% |
| Pathogenic | 2,494 | 2,131 | 66.12% | 56.50% |
| VUS | 64,224 | 57,030 | 95.70% | 84.98% |

5,144 actual VUS were predicted Likely benign; 4,011 were predicted Likely
pathogenic. Similar aggregate train/validation performance does not establish
good generalization, but does not show classic large-gap overfitting either.
Underfitting, label uncertainty and limited feature information cannot be
separated causally without controlled development experiments.

## Findings by pipeline stage

1. Target and dataset: this predicts recorded five-class labels, not patient
   disease. VUS denotes unresolved evidence, not a homogeneous biological state.
   Likely/definite distinctions require evidence strength not fully represented
   by the current biology features. The cohort includes annotation-matched
   GRCh38 SNVs only; results do not cover all four-million source records.
2. Labels: original duplicate/conflict exclusion exists. The later audit flagged
   27,308 changed labels and 11,737 newly conflicting records (overlapping sets),
   but did not resolve them. Changed labels are not automatically erroneous.
   A reviewed, versioned training-only cleaning policy is still needed; never
   relabel/drop evaluation examples to inflate a score or treat VUS as benign.
3. Coverage: full matched training cache has assay matches on 965/1,633,384
   records and TP53 SpliceAI on 502. This is a confirmed coverage limitation,
   not proof of its causal contribution to the score. Gene-level context is not
   a substitute for variant-level functional evidence. Missingness is preserved.
4. Feature representation: parser extracts nucleotide/amino substitutions, but
   the old numerical allowlist drops those identities and direct gene identity.
   Different substitutions can therefore collapse when other annotations match.
   Confirmed bug: substring-based event detection scans full names, so INS can
   trigger insertion and WFS1 frameshift flags. Its impact on measured performance
   has not been quantified. The legacy behavior is retained for saved-model
   reproducibility; a corrected input profile is now available for NEW fits.
5. Sampling: only ~20% of original training loci were used, including sampling
   minority classes. Actual fitting counts include 256,687 VUS and 4,012 Benign.
   No evidence establishes that 3–5 lakh is optimal. Compare development learning
   curves and training-only sampling policies; keep evaluation prevalence fixed.
6. Learning: old early stopping uses macro-F1; final selection prioritizes
   balanced accuracy. This is a tradeoff, not a mathematical error. It sacrifices
   precision/F1: selected validation macro-F1 0.5732 versus 0.6166 for the highest-
   accuracy single model. Early-selected 47-tree and 301-tree members may limit
   learning, but this does not prove premature stopping. Eight fits with one seed
   do not establish an optimal model. The corrected opt-in run aligns stopping
   and selection on macro-F1; per-class errors remain essential review criteria.
7. Evaluation: loci are disjoint, but genes can repeat. Validation has been reused
   and historical test scores inspected. Neither is an untouched prospective or
   unseen-gene benchmark. The later snapshot lacked sufficient Benign/Likely benign
   support. An independently frozen, adequately supported benchmark is still needed.
8. Reporting: overall accuracy was incorrectly easy to interpret as VUS recall.
   Dashboard now leads with raw/row-normalized confusion matrices, actual/predicted
   counts, per-class precision/recall/F1. Overall partition accuracy is at the end;
   balanced accuracy and macro-F1 are available under additional metrics.

## Implemented corrections and verification

- `src/variant_feature_revision.py`: new `variant_identity_v2` profile retains
  genomic/transcript substitutions, amino-acid identities, and gene context as
  categorical inputs; event detection only reads HGVS change expressions.
- `load_development(..., feature_profile='variant_identity_v2')` explicitly
  selects the new profile. Default is legacy so saved reports/inference remain
  consistent. Test loading remains prohibited. Unknown changes remain missing.
- `src.retune_biological_models --corrected --output <NEW_DIRECTORY>` uses that
  profile and selects by macro-F1, then balanced accuracy. Version/source hashes
  are recorded. New profile adds six categories (57 combined columns), and must
  be used consistently for later inference. Existing bundles have not changed.
- Regression tests verify different substitutions remain distinguishable,
  gene names cannot create events, actual HGVS events survive, labels cannot
  alter inputs, and the loader preserves legacy default behavior and test guards.
  A small synthetic integration fixture checks that revised inputs fit and return
  five normalized class probabilities; this is a software test, not training data
  or evidence of model quality. A loader SQLite connection leak found by the
  Windows fixture was corrected with explicit context-managed connection closing.
- Frontend tests reproduce saved confusion-matrix totals, precision and recall
  for every evaluated partition, including zero-support handling.
  Production frontend build and local asset/API HTTP checks passed.

These fixes are not evidence of accuracy improvement. Before a production change:
resolve a documented training-label policy; quantify feature coverage by class;
compare a single corrected model and the legacy baseline on fixed development
groups; perform feature/cleaning ablations separately; assess errors and variation
across seeds/groups; freeze the choice before a new adequately supported test.
Gene identity requires an unseen-gene evaluation if claiming generalization to
new genes. No extra synthetic data or forced class relabeling is required here.
