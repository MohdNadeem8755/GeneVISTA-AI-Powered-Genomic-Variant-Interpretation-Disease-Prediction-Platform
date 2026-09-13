"""Read-only, parameterized queries for the prepared GeneVISTA database."""

import json
import sqlite3
from contextlib import contextmanager
from pathlib import Path
from .evidence import disease_evidence

SNAPSHOTS = ["2022-01", "2023-01", "2024-01", "2025-01", "2026-01", "2026-09"]


class Database:
    def __init__(self, path):
        self.path = Path(path).resolve()

    @contextmanager
    def connect(self):
        connection = sqlite3.connect(self.path.as_uri() + "?mode=ro", uri=True, timeout=5)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA query_only = ON")
        try:
            yield connection
        finally:
            connection.close()

    def stats(self):
        with self.connect() as connection:
            metadata = dict(connection.execute("SELECT key, value FROM metadata"))
        if metadata.get("schema_version") != "1":
            raise ValueError("Unsupported database schema")
        return {
            "latest_snapshot": metadata["latest_snapshot"],
            "snapshots": json.loads(metadata["snapshots"]),
            "variant_count": int(metadata["variant_count"]),
            "timeline_count": int(metadata["timeline_count"]),
            "gene_symbol_link_count": int(metadata["gene_symbol_link_count"]),
            "scope": metadata["scope"],
        }

    def timeline(self, variant_id):
        with self.connect() as connection:
            row = connection.execute(
                "SELECT * FROM timelines WHERE variation_id = ?", (variant_id,)
            ).fetchone()
        if row is None:
            return None
        return {
            "variation_id": variant_id,
            "observations": [
                {"snapshot": snapshot, "state": row["state_" + snapshot.replace("-", "_")]}
                for snapshot in SNAPSHOTS
            ],
            "absence_definition": "absent_grch38 means absent from the prepared snapshot's GRCh38 records, not necessarily absent from ClinVar.",
        }

    def detail(self, variant_id):
        timeline = self.timeline(variant_id)
        if timeline is None:
            return None
        with self.connect() as connection:
            row = connection.execute(
                "SELECT * FROM variants WHERE variation_id = ?", (variant_id,)
            ).fetchone()
            genes = [item[0] for item in connection.execute(
                "SELECT gene_symbol FROM variant_gene_symbols WHERE variation_id = ?", (variant_id,)
            )]
        return {
            "variation_id": variant_id,
            "disease_evidence": disease_evidence(genes),
            "prediction_support": {"disease_probability": False, "clinical_class_probability": False,
                "reason": "The saved model predicts future VUS reclassification outcomes, not disease risk or five-class clinical classifications."},
            "latest_details": dict(row) if row else None,
            "latest_details_available": row is not None,
            "timeline": timeline["observations"],
            "absence_definition": timeline["absence_definition"],
        }

    def search(self, query, limit=20, after_id=0):
        query = query.strip()
        if not query or len(query) > 100:
            raise ValueError("Enter a variant ID or an exact gene symbol (1–100 characters).")
        if not 1 <= limit <= 100 or not 0 <= after_id <= 9223372036854775807:
            raise ValueError("Invalid pagination parameters.")
        if query.isascii() and query.isdecimal():
            variant_id = int(query)
            if not 1 <= variant_id <= 9223372036854775807:
                raise ValueError("Variant ID is outside the supported range.")
            result = self.detail(variant_id) if variant_id > after_id else None
            return {"query": query, "search_type": "variant_id", "items": [result] if result else [], "next_after_id": None}

        with self.connect() as connection:
            rows = connection.execute("""
                SELECT v.* FROM variant_gene_symbols g
                JOIN variants v ON v.variation_id = g.variation_id
                WHERE g.gene_symbol = ? COLLATE NOCASE
                  AND g.variation_id > ?
                GROUP BY v.variation_id
                ORDER BY v.variation_id LIMIT ?
            """, (query, after_id, limit + 1)).fetchall()
        more = len(rows) > limit
        rows = rows[:limit]
        return {
            "query": query,
            "search_type": "exact_gene_symbol",
            "items": [{"variation_id": row["variation_id"], "latest_details": dict(row), "latest_details_available": True} for row in rows],
            "next_after_id": rows[-1]["variation_id"] if more else None,
        }
