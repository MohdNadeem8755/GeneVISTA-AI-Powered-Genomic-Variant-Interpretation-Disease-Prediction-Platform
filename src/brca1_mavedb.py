"""Conservative BRCA1 MaveDB functional-score matching.

Scores are auxiliary functional evidence only. They are never converted into
patient risk, a five-class clinical label, or a model prediction.
"""
import math
import re

NAME = re.compile(r"(NM_007294\.\d+)\(BRCA1\):(c\.[1-9]\d*[ACGT]>[ACGT])(?: \(p\.([^)]*)\))?")


def parse_variant(name):
    match = NAME.fullmatch(name or "")
    if not match:
        return None
    transcript, coding_change, protein = match.groups()
    return {"transcript": transcript, "coding_change": coding_change, "protein_change": protein}


def score_classification(score, thresholds):
    if not isinstance(score, (int, float)) or not math.isfinite(score):
        return {"functional_class": "not_available", "criterion_candidate": None}
    if score < thresholds["non_functional_upper"]:
        return {"functional_class": "non-functional", "criterion_candidate": "PS3_Strong"}
    if score < thresholds["functional_upper"]:
        return {"functional_class": "intermediate", "criterion_candidate": None}
    return {"functional_class": "functional", "criterion_candidate": "BS3_Strong"}


def assess_variant(record, bundle):
    if not record or record.get("gene_symbol") != "BRCA1":
        return {"status": "out_of_scope", "reason": "This pilot covers BRCA1 SNVs only."}
    variant = parse_variant(record.get("name"))
    if not variant:
        return {"status": "unsupported_variant", "reason": "An exact BRCA1 coding SNV is required."}
    metadata = bundle.get("metadata", {})
    proof = metadata.get("transcript_verification", {})
    if (metadata.get("schema_version") != 1 or metadata.get("source_transcript") != "NM_007294.3"
            or proof.get("cds_identical") is not True or proof.get("accessions") != ["NM_007294.3", "NM_007294.4"]):
        return {"status": "incompatible_source", "reason": "MaveDB transcript verification needs rebuilding."}
    if variant["transcript"] not in {"NM_007294.3", "NM_007294.4"}:
        return {"status": "transcript_mismatch", "reason": "This transcript version is not covered by the verified mapping."}
    row = bundle.get("scores", {}).get(variant["coding_change"])
    if row is None:
        return {"status": "not_available", "reason": "This SNV is not covered by the selected MaveDB score set."}
    if variant.get("protein_change") and row.get("protein_change") not in {None, "NA", variant["protein_change"]}:
        return {"status": "mapping_conflict", "reason": "The recorded protein consequence differs from MaveDB."}
    return {"status": "available", "variant": variant, "score_set": metadata["score_set"],
            "score": row["score"], "rna_score": row.get("score_rna"),
            "functional_class": row["functional_class"],
            "criterion_candidate": row["criterion_candidate"],
            "calibration": metadata["calibration"], "review_required": True,
            "classification": None, "patient_disease_probability": None,
            "limitations": "Research functional assay evidence. This does not establish a clinical classification or patient disease risk.",
            "sources": metadata["sources"]}
