# A. Diagnostic analysis

The repository's `README.md` records 83.50% accuracy, 53.83% balanced accuracy,
and macro-F1 0.531 on an already inspected internal test of 59,112 variants.
It records recalls of 6.71% for Likely pathogenic and 6.65% for Benign. These are
historical results, not measurements from the new pipeline. Accuracy averages
over examples; balanced accuracy averages recall equally over classes. A large
gap demonstrates uneven class recall, but does not establish class imbalance as
the sole cause. Inspect prevalence, confusion pairs, annotation disagreements,
feature coverage, duplicate loci, and label-source quality before choosing a remedy.
See the [scikit-learn metric definitions](https://scikit-learn.org/stable/modules/model_evaluation.html).

The existing sequence/gene baseline has limited biological information. Suggested
features, all available at the prediction timestamp:

| Source | Features | Required audit |
| --- | --- | --- |
| Variant | Population frequency and ancestry coverage; consequence, splice, conservation and protein-change features | Missing frequency is not zero; check genome assembly and allele matching |
| Gene | Constraint, mechanism, inheritance, expression, domains | Gene holdout in addition to locus holdout |
| Phenotype | HPO terms, ontology ancestors, semantic similarity | Use independently available patient phenotypes; exclude annotations derived from the target classification |
| Evidence | ClinGen gene–disease validity and evidence count | Validity is not variant pathogenicity; retain disease and availability provenance |
| Time | Previous observed classification, days since change, previous evidence count and increment | Separate historical observation log; no shifting current target labels |
| Graph/text | Frozen graph embeddings, relation counts and biomedical text embeddings | Fit graph/text transforms on historical training material; remove target labels/edges and future citations |

Review status and submitter information can encode the process that produced the
target. Exclude them from the primary biological benchmark. Existing pathogenicity
scores may themselves have trained on overlapping ClinVar records; compare a
frequency/conservation-only ablation before claiming independent generalization.

# B. Recommended modeling strategy

The executable primary model is equal soft voting of a regularized RandomForest,
ExtraTrees, and logistic regression. A selected LightGBM, XGBoost or CatBoost
backend replaces ExtraTrees while keeping the other two learners. This avoids
requiring every optional library for a first run and supplies nonlinear and linear
models with complementary inductive biases. Select the backend before inspecting
the final test; compare alternatives only in development CV.

Optimize mean fold **balanced accuracy**, with mean fold **macro-F1** as the exact
tie-breaker. Accuracy is reported but never used for selection. A seeded randomized
parameter search using sklearn `ParameterSampler` evaluates depth, leaf size,
regularization, and class-weight strength. Each candidate gets fresh preprocessing
and estimators in each fold. The initial candidate has no class reweighting;
the remaining candidates test inverse-frequency sample weights raised to different
powers. The same folds are used for all candidates.

For mixed HPO/categorical inputs, weights are the default. Optional SMOTE,
Borderline-SMOTE and ADASYN run only inside `fit` on training folds and are rejected
for categorical/token schemas: interpolated HPO membership or gene IDs would be
biologically questionable. Numeric-only resampling remains an experiment, not an
assumption of improvement; ADASYN can fail when no suitable neighborhood exists.
Do not hide that failure or silently change methods. Missingness indicators are
also interpolated by numeric resampling, so prefer complete numeric features or
weights if their interpretation matters. See
[imbalanced-learn's leakage guidance](https://github.com/scikit-learn-contrib/imbalanced-learn/blob/master/doc/common_pitfalls.rst).

After choosing hyperparameters, refit on the training partition and optimize
per-class log-probability offsets on the separate decision-tuning partition.
Prediction is `argmax(log(p) + offset)`; multiclass classes remain mutually
exclusive. The final model is not refit afterward. Raw ensemble probabilities
are **uncalibrated model outputs**; offsets do not create calibrated probabilities.
If calibrated probabilities are required, reserve another independent chronological
calibration block before decision tuning, and assess log loss and reliability.

Graph-structured secondary option: relation-aware GNNs with historical edges,
inductive gene/variant holdouts and train-only embeddings. Text-heavy option:
frozen biomedical encoders plus the same tree/linear baseline, followed by
regularized fine-tuning if training support permits. Neural models may use focal
loss or weighted cross-entropy, tuned on development data. Class/feature-specific
experts should be combined by out-of-fold stacking, with no in-sample predictions
used to train the meta-learner.

# C. Complete Python code

Implementation: `src/robust_variant_training.py`. Legacy-data exporter:
`src/prepare_robust_pilot.py`. Regression tests: `tests/test_robust_variant_training.py`.

The trainer consumes CSV/CSV.gz and a JSON schema such as
`research/robust_variant_schema.json`. One row represents one prediction event,
with exactly one mutually exclusive target. Column names are arbitrary; only
the explicit feature lists are used. Numeric strings are coerced, nonfinite values
are marked missing, categories are one-hot encoded, and pipe-delimited HPO terms
are multi-hot encoded. Scaling and vocabulary fitting happen inside each fit.
Unknown categories at inference are ignored. IDs and timestamps are never direct
features. Numeric features are zero-filled **with separate missingness flags**.

Required schema fields:

| Field | Meaning |
| --- | --- |
| `target` | One clinical-significance class per row |
| `group` | Normalized assembly/chromosome/locus key; use a gene key for an unseen-gene experiment |
| `numeric`, `categorical`, `tokens` | Explicit, distinct feature allowlists |
| `time` | Optional prediction timestamp, mandatory for a temporal experiment |
| `label_time` | Required with `time`: when the target outcome became available |

With no `time`, splits use stratified group CV, including outer approximately
20% test and a further approximately 20% of development for decision tuning.
Group constraints make exact proportions impossible in some datasets. With
`time`, split by distinct timestamps (latest approximately 20% for test), use
expanding-window CV over unique timestamps, and keep tied dates together. Purge
overlapping groups and any training outcomes unavailable at the validation/test
cutoff. Chronology takes priority over stratification: rebalancing future rows
into the past is leakage. Every partition/fold must contain every class or the
run fails with an actionable error. Consult
[scikit-learn's grouped and temporal CV guidance](https://scikit-learn.org/stable/modules/cross_validation.html).

This is a strict unseen-locus temporal benchmark. A panel containing the same
loci at every snapshot may have no eligible training data after purging. For
forecasting changes of **known** variants, define a separate benchmark and
implement a documented entity-overlap policy; do not quietly relax these guards.

Optional `--history` accepts a separate CSV with the group column plus
`available_at`, `observed_label`, and cumulative `evidence_count`. It adds the last
strictly earlier observed class, previous count, previous count increment and
days since the last observed class change. Duplicate entity/timestamp observations
are rejected. History may include observations arising between test predictions
only when the intended deployment actually receives them then; this is sequential
as-of evaluation. For a batch forecast, freeze history at the batch cutoff.
The trainer never derives history by shifting the current target. For forecasting,
construct the future outcome and its availability time upstream; predicting the
current recorded label from its earlier version measures annotation persistence.

Outputs in a new, non-overwritable directory:

- `audit.json`: schema, seed, input/history hashes, versions and partition support.
- `partitions.csv`: source row index and partition, including purged rows.
- `development.json`: every candidate's fold and pooled reports, confusion matrices,
  accuracy, balanced accuracy, macro precision/recall/F1 and weighted/per-class metrics.
- `decision_tuning.json`: before/after reports, explicitly optimization-set estimates.
- `feature_importance.csv`: permutation importance of original feature columns
  on the tuning set, measured by balanced-accuracy decrease. This is descriptive
  and potentially optimistic after tuning; correlated features can dilute importance.
- `model.joblib`: preprocessing, ensemble, class order and selected offsets.
- `selection_frozen.json`: chosen parameters and offsets, saved before test scoring.
- `test_metrics.json`: held-out reports for the prespecified selected policy and
  the default argmax reference. Never select between them based on test results.

This implementation is single-label multiclass. HPO inputs are multi-valued, but
that does not make the target multi-label. Independent disease-label prediction
requires a separate masked multi-label loss, per-label support and PR metrics,
and per-label threshold tuning; do not feed serialized disease lists into this
classifier as if they were mutually exclusive classes.

# D. Expected realistic performance discussion

Executed pilot on 2026-10-29: 6,000 real records sampled from the legacy training
partition; 3,839 training, 960 decision-tuning and 1,201 held-out rows. Three
candidates (unweighted baseline plus two randomized configurations), three CV
folds and 40 trees per forest were used. This is an exploratory reused-source
experiment with limited sequence/gene features, not the richer temporal schema.

| Held-out metric | Default argmax | Development-selected offsets |
| --- | ---: | ---: |
| Accuracy | 82.18% | 78.85% |
| Balanced accuracy / macro recall | 49.80% | 53.67% |
| Macro precision | 51.94% | 49.35% |
| Macro F1 | 0.4707 | 0.4979 |
| Weighted F1 | 0.7991 | 0.7944 |

Selected-policy per-class results from the same held-out set:

| Class | Precision | Recall | F1 | Support |
| --- | ---: | ---: | ---: | ---: |
| Benign | 0.2424 | 0.1231 | 0.1633 | 65 |
| Likely benign | 0.8084 | 0.7345 | 0.7697 | 339 |
| Likely pathogenic | 0.1029 | 0.1892 | 0.1333 | 37 |
| Pathogenic | 0.3854 | 0.7255 | 0.5034 | 51 |
| Uncertain significance | 0.9282 | 0.9111 | 0.9196 | 709 |

Full reports and confusion matrices are in
`reports/robust_pilot_20261029/test_metrics.json`; all CV candidate reports are
in that directory's `development.json`. The selected configuration's mean CV
balanced accuracy was 0.4973 versus the unweighted candidate's 0.4155. These are
development estimates used in selection. This pilot does **not** establish an
improvement over the historical 59,112-row benchmark: the cohorts differ.

Neither 95% accuracy nor 90% balanced accuracy is guaranteed. No defensible
numerical forecast follows from the baseline alone. Improving minority recall
can lower overall accuracy. Benign versus Likely benign and Pathogenic versus
Likely pathogenic boundaries partly reflect evidence strength and curation
practice, not solely measurable molecular differences.

If targets remain unmet, investigate errors by class, gene, variant type, ancestry,
review quality and time. Obtain independent functional evidence and better minority
coverage; audit label disagreements rather than relabeling difficult cases to
inflate scores. Add richer biological features, test feature-group experts, then
consider multi-task or temporal graph models. Report confidence intervals from
group-level resampling and evaluate an untouched later snapshot or independent
curation source with audited overlap. Abstention may improve accuracy among
accepted predictions, but must always report coverage and must not replace the
full-cohort benchmark. Binary class merging changes the task and cannot establish
the original five-class target.

# E. How to run the code + checklist to avoid common errors

From the repository root, use the existing Python environment. Core sklearn mode
does not require optional dependencies. To enable every backend/resampler:

```powershell
.venv\Scripts\python.exe -m pip install -r research/requirements-robust-training.txt
$env:TEMP='D:\GeneVISTA\.runtime'
$env:TMP=$env:TEMP
.venv\Scripts\python.exe -m unittest tests.test_robust_variant_training -v
```

New temporal dataset, after preparing its schema and independently available
features/history:

```powershell
.venv\Scripts\python.exe -m src.robust_variant_training --data data/processed/events.csv.gz --schema research/robust_variant_schema.json --history data/processed/history.csv.gz --output reports/robust_temporal_NEW --backend catboost --trials 12 --folds 3 --evaluation-status new_prespecified_holdout
```

Reproduce a small exploratory pilot from the existing legacy training cache:

```powershell
.venv\Scripts\python.exe -m src.prepare_robust_pilot --output data/processed/robust_pilot_NEW --rows 6000
.venv\Scripts\python.exe -m src.robust_variant_training --data data/processed/robust_pilot_NEW/cohort.csv.gz --schema data/processed/robust_pilot_NEW/schema.json --output reports/robust_pilot_NEW --trials 2 --trees 40 --folds 3 --evaluation-status exploratory_reused_source
```

Inference from a trusted artifact, in Python launched from the repository root:

```python
import joblib
import pandas as pd
model = joblib.load("reports/robust_pilot_NEW/model.joblib")
new_rows = pd.read_csv("new_variants.csv")  # Same feature columns and units.
labels = model.predict(new_rows)
probabilities = model.predict_proba(new_rows)  # Uncalibrated, class order: model.classes_.
```

For a history-enabled model, first apply `add_history(new_rows, history,
model.schema)` using the same availability convention. Already engineered history
columns must not be added twice. New rows need only feature columns for prediction;
the target is not required. Never load untrusted joblib files.

Checklist:

- Fix the task, target horizon, cohort, class definitions and split protocol first.
- Normalize loci and audit duplicates, genome assembly, and allele identity.
- Inspect prevalence and per-class support; never optimize accuracy.
- Supply true prediction and label-availability times, not last-modified proxies.
- Join ClinGen/HPO/graph/text evidence as of the prediction time upstream.
- Exclude current target, review-derived proxies, future evidence and identifiers.
- Fit vocabularies, scaling, sampling and estimators only inside training folds.
- Keep hyperparameter selection, decision tuning and final testing separate.
- Retain complete reports and source provenance, including upstream score overlap.
- Use a new output directory; a new directory alone does not make reused data independent.
- Validate feature coverage and calibration before deployment; research outputs
  do not alter the application's deployed predictor automatically.

The tests use synthetic fixtures only to check code invariants. They are not
biomedical performance measurements. Optional library paths require execution in
an environment containing those libraries; package availability is not validation.

Verification completed on 2026-10-29: all seven regression tests passed, the
sklearn pilot completed end to end, its persisted model loaded and predicted in
a fresh Python process, and the CatBoost backend passed a fit/predict smoke test.
LightGBM, XGBoost and imbalanced-learn were absent and their optional execution
paths were not tested. The legacy cache exporter was run using
`D:\AncondaAPP\python.exe`, which has pyarrow; the trainer and tests used the
project `.venv`. All observed metrics above came from the real sklearn pilot.
