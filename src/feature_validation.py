"""Validation and baseline transforms for the temporal VUS feature notebook."""
from pathlib import Path
import re
import numpy as np
import pandas as pd

MISSING_VALUES = {'', '-', '.', 'na', 'n/a', 'nan', 'none', 'null'}
TARGETS = {'stayed_vus', 'became_benign', 'became_pathogenic'}


def load_label_manifest(label_folder, processed_folder):
    label_folder, processed_folder = Path(label_folder), Path(processed_folder)
    manifest = pd.read_csv(label_folder / 'label_manifest.csv', dtype='string', keep_default_na=False)
    required = {'file', 'split', 'feature_snapshot', 'outcome_snapshot', 'rows'} | TARGETS
    if manifest.empty or not required.issubset(manifest.columns):
        raise ValueError('Empty or incomplete label manifest')
    if not manifest['file'].is_unique or not manifest['feature_snapshot'].is_unique:
        raise ValueError('Duplicate label files or feature snapshots')
    for column in ('feature_snapshot', 'outcome_snapshot'):
        if not manifest[column].str.fullmatch(r'[0-9]{4}-(0[1-9]|1[0-2])').all():
            raise ValueError('Snapshot must have YYYY-MM format')
    if not manifest['feature_snapshot'].lt(manifest['outcome_snapshot']).all():
        raise ValueError('Outcomes must follow feature snapshots')
    order = manifest.sort_values('feature_snapshot').reset_index(drop=True)
    ranks = order['split'].map({'train': 0, 'validation': 1, 'test': 2})
    if ranks.isna().any() or not ranks.is_monotonic_increasing or not ranks.eq(0).any():
        raise ValueError('Invalid chronological train/validation/test splits')
    for rank in (1, 2):
        current = order.loc[ranks.eq(rank)]
        earlier = order.loc[ranks.lt(rank)]
        if not current.empty and not earlier.empty and earlier['outcome_snapshot'].max() > current['feature_snapshot'].min():
            raise ValueError('Earlier split outcomes extend beyond later split feature dates')
    for row in order.itertuples(index=False):
        if Path(row.file).name != row.file or '/' in row.file or '\\' in row.file:
            raise ValueError('Label filename must stay inside the label folder')
        counts = [getattr(row, name) for name in ('rows', *sorted(TARGETS))]
        if any(re.fullmatch(r'[0-9]+', value) is None for value in counts) or int(row.rows) <= 0:
            raise ValueError('Invalid manifest counts')
        if sum(int(getattr(row, name)) for name in TARGETS) != int(row.rows):
            raise ValueError('Manifest target counts do not add up')
        for path in (label_folder / row.file, processed_folder / f'clinvar_variants_{row.feature_snapshot}.csv.gz'):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
    return order


def load_labels(label_folder, row):
    labels = pd.read_csv(Path(label_folder) / row.file, dtype='string', keep_default_na=False)
    required = {'VariationID', 'feature_snapshot', 'outcome_snapshot', 'split', 'target'}
    if not required.issubset(labels.columns) or len(labels) != int(row.rows):
        raise ValueError('Label schema or row count does not match manifest')
    if not labels['VariationID'].str.fullmatch(r'[0-9]+').all() or not labels['VariationID'].is_unique:
        raise ValueError('Label IDs must be unique and nonmissing')
    for column in ('feature_snapshot', 'outcome_snapshot', 'split'):
        if not labels[column].eq(getattr(row, column)).all():
            raise ValueError(f'Label {column} differs from manifest')
    if not labels['target'].isin(TARGETS).all():
        raise ValueError('Unknown temporal target')
    counts = labels['target'].value_counts()
    if any(int(counts.get(target, 0)) != int(getattr(row, target)) for target in TARGETS):
        raise ValueError('Label target counts differ from manifest')
    return labels


def previous_snapshots(timeline_file, snapshots):
    columns = pd.read_csv(timeline_file, nrows=0, keep_default_na=False).columns.tolist()
    available = sorted(c for c in columns if re.fullmatch(r'[0-9]{4}-(0[1-9]|1[0-2])', c))
    if 'VariationID' not in columns or len(available) != len(columns) - 1:
        raise ValueError('Invalid timeline schema')
    mapping = {}
    for snapshot in snapshots:
        if snapshot not in available:
            raise ValueError(f'Snapshot {snapshot} missing from timeline')
        earlier = [value for value in available if value < snapshot]
        mapping[snapshot] = earlier[-1] if earlier else None
    return mapping


def build_baseline_features(rows, snapshot):
    required = {'VariationID', 'Type', 'ReviewStatus', 'NumberSubmitters', 'LastEvaluated'}
    if not required.issubset(rows.columns):
        raise ValueError(f'Missing source columns: {required - set(rows.columns)}')
    rows = rows.reset_index(drop=True)
    if rows['VariationID'].isna().any() or not rows['VariationID'].str.fullmatch(r'[0-9]+').all() or not rows['VariationID'].is_unique:
        raise ValueError('Source variant IDs must be unique, nonmissing numeric identifiers')
    features = pd.DataFrame(index=rows.index)

    features["VariationID"] = rows["VariationID"]
    features["feature_snapshot"] = snapshot

    # Preserve categories for an encoder fitted later on training data.
    for source, output in [
        ("Type", "variant_type"),
        ("ReviewStatus", "review_status")
    ]:
        values = rows[source].str.strip()
        features[output] = values.mask(values.isna() | values.str.lower().isin(MISSING_VALUES))

    submitters = pd.to_numeric(
        rows["NumberSubmitters"].str.strip(), errors="coerce"
    )

    valid_submitters = (
        submitters.notna()
        & np.isfinite(submitters)
        & submitters.ge(0)
        & submitters.mod(1).eq(0)
    ).fillna(False)

    features["number_submitters"] = (
        submitters.where(valid_submitters).astype("Float64")
    )
    features["submitter_count_missing"] = (
        ~valid_submitters
    ).astype("int8")

    # Use the start of the named snapshot month consistently.
    reference_date = pd.Timestamp(f"{snapshot}-01", tz="UTC")
    raw_dates = rows["LastEvaluated"].str.strip()

    date_lookup = {
        value: pd.to_datetime(value, errors="coerce", utc=True)
        for value in raw_dates.unique()
    }

    dates = pd.to_datetime(
        raw_dates.map(date_lookup),
        errors="coerce",
        utc=True
    )

    after_reference = dates.gt(reference_date)
    usable_dates = dates.notna() & ~after_reference

    features["days_since_evaluation"] = (
        (reference_date - dates)
        .dt.total_seconds()
        .div(86400)
        .where(usable_dates)
        .astype("Float64")
    )

    features["evaluation_date_missing"] = dates.isna().astype("int8")
    features["evaluation_date_after_reference"] = (
        after_reference.astype("int8")
    )

    assert features["days_since_evaluation"].dropna().ge(0).all()
    assert features["number_submitters"].dropna().ge(0).all()

    return features.reset_index(drop=True)