# Corrected feature comparison — November 9, 2026

User requested improving weaker-class predictions, especially precision of
Likely benign / Likely pathogenic and an approximately 80% precision aspiration
for Benign / Pathogenic. That aspiration is not a guaranteed result, nor grounds
to sacrifice recall silently or filter evaluation records.

Run: `reports/biological_corrected_20261109`.
Logs: `.runtime/corrected-20261109.out.log` and `.runtime/corrected-20261109.err.log`.
Command: `.venv\Scripts\python.exe -u -m src.retune_biological_models --corrected --output reports/biological_corrected_20261109`.
After interruption, verify there is no active writer, then add `--resume`.
`selection.json` indicates completion; candidate artifacts checkpoint each fit.

Use the exact previous development sample and locus-separated stopping split,
with original labels, to isolate the feature-profile/selection changes. New inputs
retain genomic/transcript substitutions, amino-acid identity and gene identity,
and correct gene-name false event flags. Stop and select on macro F1 (selection
tie-breaker balanced accuracy). Eight candidates and 84 pairwise votes are fixed
in advance. No new datasets, label edits, historical/final test access, or automatic
production replacement. Compare per-class precision, recall, F1 and confusion
matrices on identical validation support with `biological_retuning_fresh_20261031`.

This is reused development validation; any gains require independent confirmation.
Existing source coverage, label uncertainty and unseen-gene limitations remain.
The dashboard database counts are not training sample sizes: 4,243,470 retained
latest-snapshot variant IDs, 4,503,858 timeline IDs across six historical snapshots,
and 4,264,709 variant/gene-symbol pairs, not unique genes. They support search,
history and annotation lookup and must not be added together as independent rows.
