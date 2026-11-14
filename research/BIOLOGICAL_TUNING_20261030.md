# Biological-feature tuning — development only

Run: `reports/biological_tuning_20261030`. Script:
`src/tune_biological_models.py`.

Uses deterministic locus sampling (modulus 10) of existing disjoint training and
validation partitions: 163,973 training and 36,478 validation variants. No test
labels are loaded. Validation has been previously inspected and is used here for
selection, so these estimates cannot establish generalization or a 95% target.

Seven full-feature CatBoost configurations compare depths, learning rates, L2
regularization and inverse-class-frequency weight powers, including an unweighted
baseline. Each learner uses 400 iterations and seed 42. Two additional learners
use population/conservation/sequence or biological-context features. Four fixed
soft-voting mixtures combine the best full model with those feature specialists.

All runs use the same data and fixed argmax decision rule. Selection is
lexicographic balanced accuracy, then macro-F1. No accuracy-based selection,
threshold optimization, test inspection, or automatic production replacement.
This is randomized parameter search on a fixed development holdout, not nested
cross-validation. Primary features exclude pretrained pathogenicity predictors
whose evaluation overlap remains unresolved. Biological context includes TP53
SpliceAI, whose limitations are documented in the feature audit.

Outputs include `protocol.json`, `development.json`, `ensembles.json`,
`selection.json`, individual candidate artifacts and `selected.joblib`. Every
candidate report includes per-class precision/recall/F1, confusion matrix,
balanced accuracy, macro and weighted averages. `selection.json` is authoritative
for the winner when the run finishes.

## Completed results

All nine learners and four voting mixtures completed. The selection rule chose
the 2:1:1 soft-voting ensemble of the best full-feature learner, the
population/sequence specialist, and the biological-context specialist.

| Development result | Accuracy | Balanced accuracy | Macro F1 |
| --- | ---: | ---: | ---: |
| Selected 2:1:1 ensemble | 69.38% | 66.61% | 0.4834 |
| Highest-accuracy candidate | 91.24% | 54.93% | 0.6144 |
| Best balanced-accuracy single model | 68.74% | 66.59% | 0.4614 |

The ensemble was selected because its 66.61% balanced accuracy slightly exceeded
the best single model's 66.59% and the other ensemble mixtures. The tie-breaker
was not needed. The selected artifact is valid as a development candidate, but
its very small selection advantage may be noise and requires independent testing.

Neither the 95% accuracy target nor high balanced accuracy was achieved. The
unweighted model's 91.24% accuracy hides weak minority recall, reflected by its
54.93% balanced accuracy. No test data was evaluated and the deployed model was
not changed. The selected artifact reloaded successfully and returned finite,
normalized five-class probabilities; this validates serialization, not quality.

```powershell
$env:OMP_NUM_THREADS='4'
$env:OPENBLAS_NUM_THREADS='4'
.venv\Scripts\python.exe -m src.tune_biological_models --output reports/biological_tuning_NEW --modulus 10 --iterations 400
```

Load the resulting research artifact with `joblib.load` and use
`src.tune_biological_models.predict_bundle(bundle, feature_frame)` to return
probabilities in `bundle['classes']` order. It restores missing-category handling
and class alignment across members. Outputs are uncalibrated research estimates.

Next evaluation must freeze the candidate and use an untouched dataset with
compatible feature availability, source/version provenance and locus isolation.
The deployed model remains unchanged until such evaluation justifies replacement.
