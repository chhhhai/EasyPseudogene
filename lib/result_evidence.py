"""Shared sequence-disruption evidence labels, independent of filter decisions.

These grades are descriptive rules, not calibrated biological confidence.
Miniprot identity is a PAF match-fraction proxy; GeneWise uses aligned query
amino-acid identity. The method and metric must accompany every grade.
"""
from pathlib import Path

from alignment_metrics import aligned_query_length, query_coverage, sequence_identity
from genewise_alignment import segments

EVIDENCE_FIELDS = ("evidence_query_id", "stop_count", "frameshift_count", "query_coverage", "alignment_identity",
                   "aligned_length", "identity_metric", "evidence_method", "evidence_level", "evidence_reason")
EVIDENCE_SEMANTICS = {
    "scope": "sequence disruption support; not a probability of true pseudogene status",
    "high": "detected stop/frameshift, query coverage >=0.8 and identity >=0.5",
    "medium": "detected stop/frameshift, query coverage >=0.2 and identity >=0.3",
    "low": "assessed alignment without a disruption or below the evidence thresholds",
    "unassessed": "no structural validation or incomplete alignment measurements",
    "filter_independence": "A high evidence grade can still be excluded by copy or other filters",
    "identity_metrics": "GeneWise: amino-acid identity; miniprot: PAF match-fraction proxy",
}


def grade(stops, frameshifts, coverage, identity, assessed):
    if not assessed or coverage is None or identity is None or not 0 <= coverage <= 1 or not 0 <= identity <= 1:
        return "unassessed", "structural_evidence_not_assessed_or_incomplete"
    if not stops + frameshifts:
        return "low", "no_detected_sequence_disruption"
    if coverage >= 0.8 and identity >= 0.5:
        return "high", "disruption_with_strong_alignment_support"
    if coverage >= 0.2 and identity >= 0.3:
        return "medium", "disruption_with_partial_alignment_support"
    return "low", "disruption_with_weak_alignment_support"


def evidence(stops=None, frameshifts=None, coverage=None, identity=None, aligned_length=None,
             method="not_assessed", assessed=False, metric=""):
    assessed = assessed and stops is not None and frameshifts is not None
    level, reason = grade(stops or 0, frameshifts or 0, coverage, identity, assessed)
    return dict(stop_count="" if stops is None else stops,
                frameshift_count="" if frameshifts is None else frameshifts,
                query_coverage="" if coverage is None else round(coverage, 4),
                alignment_identity="" if identity is None else round(identity, 4),
                aligned_length="" if aligned_length is None else aligned_length,
                identity_metric=metric, evidence_method=method, evidence_level=level, evidence_reason=reason)


def genewise_evidence(path, protein_length, stops, frameshifts, completed):
    path = Path(path)
    if not completed or not path.is_file():
        return evidence(method="GeneWise")
    parsed = list(segments(path.read_text()))
    if not parsed or not protein_length:
        return evidence(stops, frameshifts, method="GeneWise")
    query = "".join(q for q, _ in parsed)
    target = "".join(t for _, t in parsed)
    return evidence(stops, frameshifts, query_coverage(query, target, protein_length),
                    sequence_identity(query, target), aligned_query_length(query, target),
                    "GeneWise", True, "query_amino_acid_identity")
