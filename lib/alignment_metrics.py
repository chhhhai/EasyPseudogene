"""Quality measures on GeneWise protein-alignment columns.

GeneWise pretty alignments include spaces for introns and dashes for gaps.
Neither represents a query amino acid; counting them can yield coverage above
one and artificially high identity from matching spaces.
"""
from bisect import bisect_left


def query_positions(query, target):
    return [i for i, aa in enumerate(query[:len(target)]) if aa.isalpha()]


def aligned_query_length(query, target):
    return len(query_positions(query, target))


def query_coverage(query, target, protein_length):
    if not protein_length:
        return 0.0
    return min(aligned_query_length(query, target), protein_length) / protein_length


def sequence_identity(query, target, positions=None):
    if positions is None:
        positions = query_positions(query, target)
    if not positions:
        return 0.0
    return sum(query[i].upper() == target[i].upper() for i in positions) / len(positions)


def flank_identities(query, target, column, flank):
    positions = query_positions(query, target)
    offset = bisect_left(positions, column)
    after = offset + int(offset < len(positions) and positions[offset] == column)
    before_positions = positions[max(0, offset - flank):offset]
    after_positions = positions[after:after + flank]
    return (sequence_identity(query, target, before_positions),
            sequence_identity(query, target, after_positions))


def local_flank_identity(query, target, column, flank):
    """Identity of the combined query-amino-acid window around a disruption."""
    positions = query_positions(query, target)
    offset = bisect_left(positions, column)
    after = offset + int(offset < len(positions) and positions[offset] == column)
    local = positions[max(0, offset - flank):offset] + positions[after:after + flank]
    return sequence_identity(query, target, local)


def mutation_query_fraction(query, target, column):
    positions = query_positions(query, target)
    if not positions:
        return 0.0, 0, 0
    preceding = bisect_left(positions, column)
    covered = preceding + int(preceding < len(positions) and positions[preceding] == column)
    return covered / len(positions), covered, len(positions)
