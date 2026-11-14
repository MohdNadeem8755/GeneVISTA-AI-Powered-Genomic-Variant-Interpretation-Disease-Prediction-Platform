# Continue Step 4

COMPLETED 2026-10-31: Fresh run `reports/biological_retuning_fresh_20261031`
finished all eight fits and 84 ensembles. Do not restart it. Selected development
ensemble accuracy 81.51%, balanced accuracy 64.98%, macro F1 0.5732; highest
accuracy candidate 91.26%, balanced accuracy 55.39%. Production unchanged.
See `research/STEP4_FRESH_20261030.md` and the run's `selection.json`.

UPDATE 2026-10-30: User requested a fresh restart. The earlier run below is
complete. Current run is `reports/biological_retuning_fresh_20261031`; see
`research/STEP4_FRESH_20261030.md` for the current protocol and resume command.
Logs: `.runtime/step4-fresh-20261031.out.log` and `.runtime/step4-fresh-20261031.err.log`.
Do not restart the old run or run two writers for the new run.

## Earlier completed run

Current output: `reports/biological_tuning_20261030`.
Run session at dispatch: 67495. Check whether it is still running before launching
another process. If `selection.json` exists, training/selection finished; inspect
results and verify the selected artifact rather than rerunning.

If the process stopped, resume completed candidates:

```powershell
$env:OMP_NUM_THREADS='4'
$env:OPENBLAS_NUM_THREADS='4'
.venv\Scripts\python.exe -m src.tune_biological_models --output reports/biological_tuning_20261030 --modulus 10 --iterations 400 --resume
```

Each completed candidate is saved individually. An interrupted unfinished
candidate restarts; completed candidates are loaded, not retrained. Do not run
two writers simultaneously. A partial/corrupt candidate artifact needs inspection
before retrying. Preserve original data and settings.

Finish by checking selected-artifact reload/probabilities, reporting balanced
accuracy and macro-F1 with development-only caveats, and providing the local link
http://127.0.0.1:8000/ after verifying or starting the server. The user explicitly
requested a local browser link in the final answer. Do not replace the active
model based solely on development selection metrics.
