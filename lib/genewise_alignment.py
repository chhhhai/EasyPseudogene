"""Read GeneWise pretty alignments, verifying X against the displayed DNA codon.

GeneWise displays coding-strand bases in three vertical rows on both strands.
An ambiguous codon translated as X is not evidence for a stop. Query labels can
be truncated by GeneWise, so use the row structure rather than label equality.
"""
import re

ROW = re.compile(r'^\S+\s+-?\d+\s+(\S.*)$')
STOPS = {'TAA', 'TAG', 'TGA'}
BASES = set('ACGTN')


def segments(block):
    lines = block.splitlines()
    for i, line in enumerate(lines[:-5]):
        match = ROW.match(line)
        if not match:
            continue
        offset = match.start(1)
        query = match.group(1).rstrip()
        dna_rows = lines[i + 3:i + 6]
        # GeneWise uses fixed columns. A long genomic label can fill the space
        # before its coordinate, so a whitespace-based row regex misses it.
        genomic_row = dna_rows[0]
        if not ROW.match(genomic_row):
            prefix = genomic_row[:offset].strip()
            if ':' not in prefix or not any(
                all(len(row) > offset + j and row[offset + j].upper() in BASES
                    for row in dna_rows)
                for j in range(len(query))
            ):
                continue
        target = lines[i + 2][offset:offset + len(query)].ljust(len(query))
        bases = [s.ljust(offset + len(query)) for s in dna_rows]
        verified = []
        for j, aa in enumerate(target):
            codon = ''.join(s[offset + j] for s in bases).upper()
            if aa in 'X*':
                verified.append('X' if codon in STOPS else '?')
            else:
                verified.append(aa)
        yield query, ''.join(verified)


def mutation_counts(block):
    alignments = list(segments(block))
    if alignments:
        target = ''.join(t for _, t in alignments)
        return target.count('X'), target.count('!')
    # Legacy compact outputs. Never infer a stop from an unverified plain X.
    return max(block.count('*'), block.count(':X[')), block.count('!')
