"""Build a compact, reproducible BRCA1 MaveDB evidence package from cached files."""
import argparse
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SCORE_SET = "urn:mavedb:00000097-0-2"
SCORES_URL = "https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000097-0-2/scores"
METADATA_URL = "https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000097-0-2"
NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_007294.3,NM_007294.4&rettype=gb&retmode=text"
SNV = re.compile(r"c\.([1-9]\d*)([ACGT])>([ACGT])$")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_transcripts(text):
    records = {}
    for block in text.split("//"):
        match = re.search(r"^VERSION\s+(NM_007294\.\d+)", block, re.M)
        cds = re.search(r"^     CDS\s+(\d+)\.\.(\d+)\s*$", block, re.M)
        if not match or not cds or "\nORIGIN" not in block:
            continue
        sequence = re.sub("[^acgt]", "", block.split("\nORIGIN", 1)[1]).upper()
        start, end = map(int, cds.groups())
        coding = sequence[start - 1:end]
        if len(coding) != 5592:
            raise ValueError("Unexpected BRCA1 CDS length.")
        records[match[1]] = coding
    if set(records) != {"NM_007294.3", "NM_007294.4"}:
        raise ValueError("Both BRCA1 transcript versions are required.")
    if records["NM_007294.3"] != records["NM_007294.4"]:
        raise ValueError("BRCA1 transcript coding sequences differ; mapping is not safe.")
    return {"accessions": ["NM_007294.3", "NM_007294.4"], "cds_identical": True,
            "cds_length": len(records["NM_007294.3"]),
            "cds_sha256": hashlib.sha256(records["NM_007294.3"].encode()).hexdigest(),
            "source_url": NCBI_URL}


def calibration(metadata):
    selected = next((item for item in metadata.get("scoreCalibrations", [])
                     if item.get("title", "").startswith("IGVF Coding Variant Focus Group")
                     and "Missense" in item.get("title", "")), None)
    if selected is None:
        raise ValueError("Expected IGVF missense calibration was not found.")
    ranges = {item["functionalClassification"]: item["range"]
              for item in selected["functionalClassifications"]}
    return {"urn": selected["urn"], "title": selected["title"],
            "non_functional_upper": ranges["abnormal"][1] or -1.328,
            "functional_upper": ranges["normal"][0] or -0.748}


def build(source_dir):
    metadata = json.loads((source_dir / "metadata.json").read_text(encoding="utf-8"))
    if metadata.get("urn") != SCORE_SET or metadata.get("numVariants") != 3893:
        raise ValueError("Unexpected MaveDB score-set identity or size.")
    if metadata.get("license", {}).get("shortName") != "CC0":
        raise ValueError("Expected CC0 MaveDB licence.")
    target = metadata.get("targetGenes", [{}])[0].get("targetAccession", {})
    if target.get("accession") != "NM_007294.3" or target.get("gene") != "BRCA1":
        raise ValueError("Unexpected BRCA1 target transcript.")
    proof = verify_transcripts((source_dir / "NM_007294_versions.gb").read_text(encoding="utf-8"))
    thresholds = calibration(metadata)
    scores = {}
    with (source_dir / "scores.csv").open(newline="", encoding="utf-8") as handle:
        for line, row in enumerate(csv.DictReader(handle), 2):
            change = row["hgvs_nt"].split(":", 1)[-1]
            if not SNV.fullmatch(change):
                continue
            if change in scores:
                raise ValueError(f"Duplicate MaveDB allele at row {line}.")
            score = float(row["score"])
            rna = None if row.get("score_rna") in {None, "", "NA"} else float(row["score_rna"])
            classification = __import__("src.brca1_mavedb", fromlist=["score_classification"]).score_classification(score, thresholds)
            scores[change] = {"protein_change": None if row.get("hgvs_pro") in {None, "", "NA"} else row["hgvs_pro"],
                              "score": score, "score_rna": rna, **classification}
    # The score set includes exonic coding SNVs and intronic splice SNVs. This
    # pilot imports only exact coding HGVS alleles; splice rows remain counted
    # in provenance rather than being silently treated as missing data.
    if len(scores) != 2803:
        raise ValueError("The score table does not contain the expected coding-SNV count.")
    package_metadata = {"schema_version": 1, "score_set": SCORE_SET,
        "source_transcript": "NM_007294.3", "gene": "BRCA1", "num_variants": len(scores),
        "source_score_rows": metadata["numVariants"], "excluded_splice_rows": metadata["numVariants"] - len(scores),
        "transcript_verification": proof, "calibration": thresholds,
        "sources": [{"label": "MaveDB BRCA1 SGE score set", "url": METADATA_URL},
                    {"label": "BRCA1 saturation genome editing study", "url": "https://doi.org/10.1038/s41586-018-0461-z"}],
        "files": {name: {"sha256": sha256(source_dir / name), "bytes": (source_dir / name).stat().st_size}
                  for name in ["metadata.json", "scores.csv", "NM_007294_versions.gb"]},
        "license": "CC0", "license_url": "https://creativecommons.org/publicdomain/zero/1.0/",
        "attribution": "MaveDB urn:mavedb:00000097-0-2; score table reduced to exact coding SNVs by GeneVISTA.",
        "created_utc": datetime.now(timezone.utc).isoformat()}

    return {"metadata": package_metadata, "scores": scores}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, default=ROOT / "data/incoming/brca1_mavedb")
    parser.add_argument("--output", type=Path, default=ROOT / "data/processed/brca1_mavedb_v1.json.gz")
    args = parser.parse_args()
    bundle = build(args.source_dir)
    payload = json.dumps(bundle, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(".tmp")
    temporary.write_bytes(gzip.compress(payload, mtime=0))
    temporary.replace(args.output)
    print(json.dumps({"score_set": bundle["metadata"]["score_set"], "variants": len(bundle["scores"]),
                      "artifact_bytes": args.output.stat().st_size}, indent=2))


if __name__ == "__main__":
    main()
