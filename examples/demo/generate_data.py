#!/usr/bin/env python3
"""Regenerate the bundled synthetic installation data; no external data needed."""
import hashlib
import json
import random
from pathlib import Path


def fasta(records):
    return ''.join('>' + name + '\n' + '\n'.join(sequence[i:i + 60]
                   for i in range(0, len(sequence), 60)) + '\n'
                   for name, sequence in records)


def main():
    output = Path(__file__).resolve().parent / 'data'
    output.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20261013)
    codons = dict(zip('ACDEFGHIKLMNPQRSTVWY',
                      ['GCT', 'TGC', 'GAT', 'GAA', 'TTC', 'GGT', 'CAT', 'ATT',
                       'AAA', 'CTG', 'ATG', 'AAT', 'CCT', 'CAA', 'CGT', 'TCT',
                       'ACT', 'GTT', 'TGG', 'TAT']))
    proteins = ['M' + ''.join(rng.choice(tuple(codons)) for _ in range(299))
                for _ in range(2)]
    coding = [''.join(codons[aa] for aa in protein) + 'TAA'
              for protein in proteins]
    # Internal stop at amino acid 140; remove one base in codon 216.
    disrupted = coding[0][:417] + 'TAA' + coding[0][420:]
    disrupted = disrupted[:646] + disrupted[647:]

    def contig(sequence):
        return (''.join(rng.choice('ACGT') for _ in range(600)) + sequence
                + ''.join(rng.choice('ACGT') for _ in range(600)))

    cross = [('cross_disrupted_A', contig(disrupted)),
             ('cross_intact_B', contig(coding[1]))]
    own_a, own_b = contig(coding[0]), contig(coding[1])
    self_records = [('self_gene_A', own_a), ('self_gene_B', own_b),
                    ('self_disrupted_A', contig(disrupted))]
    names = ['demo_A.t1', 'demo_B.t1']
    annotation = ['##gff-version 3']
    for letter in ['A', 'B']:
        chrom, gene, transcript = f'self_gene_{letter}', f'demo_{letter}', f'demo_{letter}.t1'
        annotation += [f'##sequence-region {chrom} 1 2103',
                       f'{chrom}\tdemo\tgene\t601\t1503\t.\t+\t.\tID={gene}',
                       f'{chrom}\tdemo\tmRNA\t601\t1503\t.\t+\t.\tID={transcript};Parent={gene}',
                       f'{chrom}\tdemo\tCDS\t601\t1500\t.\t+\t0\tID={transcript}.cds;Parent={transcript}']
    files = {
        'proteins.fa': fasta([(name + f' gene=demo_{letter}', sequence)
                              for name, letter, sequence in zip(names, 'AB', proteins)]),
        'cross_genome.fa': fasta(cross),
        'self_genome.fa': fasta(self_records),
        'self_annotation.gff3': '\n'.join(annotation) + '\n',
    }
    for name, content in files.items():
        (output / name).write_text(content)
    truth = dict(origin='Synthetic; generated with Python standard library',
                 seed=20261013, protein_length_aa=300,
                 mutations=dict(query='demo_A.t1', stop_codon_aa1=140,
                                deleted_base_in_coding_sequence0=646),
                 cross=dict(disrupted_query='demo_A.t1', disrupted_contig='cross_disrupted_A',
                            intact_query='demo_B.t1', intact_contig='cross_intact_B'),
                 self=dict(annotated_queries=names,
                           unannotated_disrupted_contig='self_disrupted_A',
                           note='An intact annotated copy of demo_A.t1 is present; standard copy filters can exclude its disrupted duplicate.'),
                 genome_bases=dict(cross=sum(len(s) for _, s in cross),
                                   self=sum(len(s) for _, s in self_records)),
                 files={name: hashlib.sha256(content.encode()).hexdigest()
                        for name, content in files.items()})
    (output / 'truth.json').write_text(json.dumps(truth, indent=2) + '\n')
    print(f'Wrote synthetic demo data to {output}')


if __name__ == '__main__':
    main()
