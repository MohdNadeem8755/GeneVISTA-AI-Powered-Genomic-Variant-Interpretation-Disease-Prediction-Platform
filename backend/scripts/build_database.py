r"""Build GeneVISTA's read-only app data from the prepared CSV files.

Run: python backend/scripts/build_database.py --project-root D:\GeneVISTA
Uses only the Python standard library. Does not train or load a model.
"""

import argparse
import csv
import gzip
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path


SNAPSHOTS = ["2022-01", "2023-01", "2024-01", "2025-01", "2026-01", "2026-09"]
STATES = {"absent_grch38", "present_excluded", "benign", "pathogenic", "vus", "conflicting", "other"}
LATEST = SNAPSHOTS[-1]
VARIANT_FIELDS = [
    "VariationID", "Name", "GeneSymbol", "GeneID", "HGNC_ID", "Type",
    "ClinicalSignificance", "clinical_state", "ReviewStatus", "NumberSubmitters", "snapshot"
]
STATE_COLUMNS = ["state_" + snapshot.replace("-", "_") for snapshot in SNAPSHOTS]


def positive_id(value):
    number = int(value)
    if number <= 0:
        raise ValueError(f"Invalid variant ID: {value!r}")
    return number


def csv_batches(path, required, batch_size=25000):
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in {path.name}: {sorted(missing)}")
        batch = []
        for row in reader:
            batch.append(row)
            if len(batch) == batch_size:
                yield batch
                batch = []
        if batch:
            yield batch


def build_database(project_root):
    processed = project_root / "data" / "processed"
    variant_file = processed / f"clinvar_variants_{LATEST}.csv.gz"
    timeline_file = processed / "clinvar_state_timelines.csv.gz"
    for source in [variant_file, timeline_file]:
        if not source.is_file():
            raise FileNotFoundError(source)

    output_folder = project_root / "data" / "app"
    output_folder.mkdir(parents=True, exist_ok=True)
    destination = output_folder / "genevista.sqlite"
    if destination.exists():
        raise FileExistsError(f"Database already exists; it was not changed: {destination}")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
    temporary = output_folder / f"genevista_{stamp}.partial.sqlite"
    connection = sqlite3.connect(temporary)

    try:
        connection.execute("PRAGMA foreign_keys = ON")
        connection.executescript("""
            CREATE TABLE variants (
                variation_id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                gene_symbol TEXT NOT NULL,
                ncbi_gene_id_raw TEXT NOT NULL,
                hgnc_id_raw TEXT NOT NULL,
                variant_type TEXT NOT NULL,
                clinical_significance TEXT NOT NULL,
                clinical_state TEXT NOT NULL,
                review_status TEXT NOT NULL,
                number_submitters_raw TEXT NOT NULL,
                snapshot TEXT NOT NULL
            );
            CREATE TABLE variant_gene_symbols (
                variation_id INTEGER NOT NULL REFERENCES variants(variation_id),
                gene_symbol TEXT NOT NULL,
                PRIMARY KEY (variation_id, gene_symbol)
            );
            CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        """)
        columns = ", ".join(f'"{name}" TEXT NOT NULL' for name in STATE_COLUMNS)
        connection.execute(f"CREATE TABLE timelines (variation_id INTEGER PRIMARY KEY, {columns})")

        variant_count = 0
        symbol_count = 0
        for batch in csv_batches(variant_file, VARIANT_FIELDS):
            records = []
            symbols = []
            for row in batch:
                variant_id = positive_id(row["VariationID"])
                if row["snapshot"] != LATEST or row["clinical_state"] not in STATES - {"absent_grch38", "present_excluded"}:
                    raise ValueError(f"Unexpected latest-snapshot state for {variant_id}")
                records.append((variant_id, *(row[field] for field in VARIANT_FIELDS[1:])))
                # These are source-reported symbols, not inferred ID mappings.
                names = {name.strip() for name in row["GeneSymbol"].split(";")}
                symbols.extend((variant_id, name) for name in sorted(names) if name and name.lower() not in {"-", "na", "nan"})

            with connection:
                connection.executemany("INSERT INTO variants VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", records)
                connection.executemany("INSERT INTO variant_gene_symbols VALUES (?, ?)", symbols)
            variant_count += len(records)
            symbol_count += len(symbols)
            if variant_count % 250000 == 0:
                print(f"Variant records loaded: {variant_count:,}", flush=True)

        print(f"All variant records loaded: {variant_count:,}", flush=True)
        timeline_count = 0
        for batch in csv_batches(timeline_file, ["VariationID"] + SNAPSHOTS):
            records = []
            for row in batch:
                values = [row[snapshot] for snapshot in SNAPSHOTS]
                if any(state not in STATES for state in values):
                    raise ValueError(f"Invalid timeline state for {row['VariationID']}")
                records.append((positive_id(row["VariationID"]), *values))
            with connection:
                connection.executemany("INSERT INTO timelines VALUES (?, ?, ?, ?, ?, ?, ?)", records)
            timeline_count += len(records)
            if timeline_count % 500000 == 0:
                print(f"Timelines loaded: {timeline_count:,}", flush=True)

        print("Building the gene-search index and checking the database...", flush=True)
        with connection:
            connection.execute("CREATE INDEX gene_symbol_lookup ON variant_gene_symbols (gene_symbol COLLATE NOCASE, variation_id)")

        if variant_count == 0 or timeline_count == 0:
            raise ValueError("An input dataset is empty.")
        if connection.execute("SELECT COUNT(*) FROM variants").fetchone()[0] != variant_count:
            raise ValueError("Variant count mismatch.")
        if connection.execute("SELECT COUNT(*) FROM timelines").fetchone()[0] != timeline_count:
            raise ValueError("Timeline count mismatch.")

        latest_column = STATE_COLUMNS[-1]
        mismatch = connection.execute(f"""
            SELECT COUNT(*) FROM variants v LEFT JOIN timelines t
            ON v.variation_id = t.variation_id
            WHERE t.variation_id IS NULL OR v.clinical_state != t.{latest_column}
        """).fetchone()[0]
        if mismatch:
            raise ValueError(f"{mismatch:,} variant records disagree with the latest timeline.")
        missing_variants = connection.execute(f"""
            SELECT COUNT(*) FROM timelines t LEFT JOIN variants v
            ON t.variation_id = v.variation_id
            WHERE t.{latest_column} NOT IN ('absent_grch38', 'present_excluded')
            AND v.variation_id IS NULL
        """).fetchone()[0]
        if missing_variants:
            raise ValueError("Retained latest-snapshot variants are missing from the lookup.")
        if connection.execute("PRAGMA foreign_key_check").fetchone() is not None:
            raise ValueError("Broken gene-symbol links.")

        metadata = {
            "schema_version": "1",
            "latest_snapshot": LATEST,
            "snapshots": json.dumps(SNAPSHOTS),
            "variant_count": str(variant_count),
            "timeline_count": str(timeline_count),
            "gene_symbol_link_count": str(symbol_count),
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "variant_source": str(variant_file),
            "timeline_source": str(timeline_file),
            "scope": "Latest retained variant details; all prepared GRCh38 timelines. No model predictions or HPO/ClinGen evidence tables yet."
        }
        with connection:
            connection.executemany("INSERT INTO metadata VALUES (?, ?)", metadata.items())
        if connection.execute("PRAGMA quick_check").fetchall() != [("ok",)]:
            raise ValueError("SQLite integrity check failed.")
        connection.close()

        # On Windows, rename fails rather than overwriting an existing file.
        if destination.exists():
            raise FileExistsError(destination)
        temporary.rename(destination)
    except Exception:
        connection.close()
        print(f"Build did not complete. Diagnostic partial file retained: {temporary}", flush=True)
        raise

    print("\nAPP DATABASE VERIFIED", flush=True)
    print(f"Variant details: {variant_count:,}")
    print(f"Timelines: {timeline_count:,}")
    print(f"Gene-symbol links: {symbol_count:,}")
    print(f"Database: {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(r"D:\GeneVISTA"))
    arguments = parser.parse_args()
    build_database(arguments.project_root.resolve())
