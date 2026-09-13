# Preparing the next annotation experiment

Run from `D:\GeneVISTA` using the local Python installation. Source downloads, cohort databases and reports stay on D:. These commands do not overwrite the production predictor or existing training scripts.

```powershell
$env:TEMP='D:\GeneVISTA\.runtime'
$env:TMP=$env:TEMP
$env:PYTHONDONTWRITEBYTECODE='1'
$env:OMP_NUM_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
& D:\AncondaAPP\python.exe -m src.annotation_pipeline prepare
```

Preparation scans the full cached ClinVar source and samples approximately 10% of genomic positions, retaining exact five-class germline/unknown-origin GRCh38 SNVs. The audit counts unsupported variants, duplicate source rows, ambiguous mappings and conflicting labels. Combined/unsupported labels are excluded, not relabeled. Different alternative alleles at a position stay in the same split. The new internal test excludes the original v1 sampled SNV loci and published example sites. The old models and their existing metrics remain unchanged. This first annotation experiment does not cover indels, structural variants or the complete ClinVar catalog.

After the download completes:

```powershell
& D:\AncondaAPP\python.exe -m src.annotation_pipeline merge --source 'D:\GeneVISTA\data\incoming\dbNSFP5.4a_grch38.gz' --checksum 'D:\GeneVISTA\data\incoming\dbNSFP5.4a_grch38.gz.md5'
& D:\AncondaAPP\python.exe -m src.annotation_pipeline compare
```

The merge verifies the compressed file's MD5 before import and reads gzip/BGZF to EOF, checking CRC. It matches chromosome, position, REF and ALT exactly, retaining only an explicit numeric score/frequency allowlist in SQLite. Raw clinical labels, IDs and review status never enter the comparison feature matrix. Missing frequency stays missing. Duplicate/transcript score entries use the maximum finite value; this exploratory aggregation is documented and is not transcript-specific inference. The version and verified checksum are recorded. An interrupted/invalid import stays unavailable to training and can be rerun from the start. The portable streaming fallback scans the entire 55 GB file and can take substantial time; it never decompresses the whole file onto disk. The `.tbi` can be kept alongside the download for a future indexed reader; this fallback does not require it.

Comparison runs a bounded per-class pilot of logistic regression and histogram gradient boosting using train and validation only. Calibration and test rows are not loaded. It reports accuracy, balanced accuracy and each class's precision/recall plus macro F1. Coordinate-ordered per-class caps are not a population-representative sample. Comparison does not automatically promote or publish a model. Before final training/evaluation, audit upstream score training overlap, feature coverage, transcript choices, labels and source licenses, then fit/calibrate on separate partitions and evaluate the locked test once. Additional score files cannot guarantee 99% accuracy.

## Deployment measurement

```powershell
& D:\AncondaAPP\python.exe deployment/profile_inference.py
& D:\AncondaAPP\python.exe deployment/profile_inference.py --database 'D:\GeneVISTA\data\app\genevista.sqlite' --output 'D:\GeneVISTA\reports\deployment\full_database_profile.json'
```

The profiler requires the local development package `psutil`; it is not a production dependency. Reports measure a fresh local process with both predictors and representative requests. Windows working-set measurements do not guarantee Render/Linux behavior or concurrency limits. A production decision also needs Linux/container testing and annotation lookup size after import. Keep raw dbNSFP out of the image and Git. New variants still need matching annotations or an explicit annotated-input workflow; the saved classifier alone cannot supply absent biological evidence.
