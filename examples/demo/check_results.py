#!/usr/bin/env python3
"""Verify the bundled synthetic demo; failures return a nonzero exit status."""
import argparse
import csv
import hashlib
import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlparse


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.targets = []

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ('href', 'src') and value:
                self.targets.append(value)


def check(output):
    data = Path(__file__).resolve().parent / 'data'
    truth = json.loads((data / 'truth.json').read_text())
    for name, expected in truth['files'].items():
        require(hashlib.sha256((data / name).read_bytes()).hexdigest() == expected,
                f'Demo data changed: {name}')
    checked = {}
    for mode in ('cross', 'self'):
        run = output / mode
        for rel in ('report.html', 'results/final_loci.tsv', 'results/final_loci.bed',
                    'results/final_loci.gff3', 'results/final_queries.txt',
                    'results/loci.tsv', 'results/queries.tsv', 'results/summary.json',
                    'qc/input_qc.json', 'metadata/run_manifest.json',
                    'metadata/reference_cache.json', 'logs/run.log'):
            require((run / rel).is_file(), f'Missing output: {run / rel}')
        cache = json.loads((run / 'metadata/reference_cache.json').read_text())
        require(cache.get('status') == 'complete', f'{mode}: run cache is not complete')
        require(cache['signature']['tools'].get('genewise'), f'{mode}: GeneWise was not resolved')
        require(cache['signature'].get('wise_config'), f'{mode}: configuration was not recorded')
        query_rows = rows(run / 'results/queries.tsv')
        require(len(query_rows) == 2 and {r['query_id'] for r in query_rows} == {'demo_A.t1', 'demo_B.t1'},
                f'{mode}: query IDs differ from bundled demo')
        final = rows(run / 'results/final_loci.tsv')
        expected_queries = {'demo_A.t1'} if mode == 'cross' else set()
        require({r['query_id'] for r in final} == expected_queries,
                f'{mode}: unexpected final candidate queries')
        require(len(final) == (1 if mode == 'cross' else 0),
                f'{mode}: expected final locus count differs')
        if mode == 'cross':
            require(int(final[0]['stop_count']) >= 1 and int(final[0]['frameshift_count']) >= 1,
                    'cross: engineered stop/frameshift was not recovered')
            intact = [r for r in rows(run / 'results/loci.tsv') if r['query_id'] == 'demo_B.t1']
            require(len(intact) == 1 and intact[0]['decision_reason'] == 'no_genewise_disruption',
                    'cross: intact negative control did not pass its expected check')
        else:
            selected = rows(run / 'results/loci.tsv')
            require(len(selected) == 1 and selected[0]['query_id'] == 'demo_A.t1'
                    and selected[0]['chrom'] == 'self_disrupted_A'
                    and selected[0]['decision_reason'] == 'other_full_length_homolog_filter'
                    and selected[0]['genewise_status'] == 'completed',
                    'self: annotation/copy-filter control did not pass its expected check')
        bed = [line.split('\t') for line in (run / 'results/final_loci.bed').read_text().splitlines() if line]
        gff = [line.split('\t') for line in (run / 'results/final_loci.gff3').read_text().splitlines()
               if line and not line.startswith('#')]
        require(len(bed) == len(gff) == len(final), f'{mode}: export counts disagree')
        for locus, interval, feature in zip(final, bed, gff):
            require(len(interval) >= 6 and len(feature) == 9, f'{mode}: malformed export row')
            require(interval[0] == feature[0] == locus['chrom']
                    and interval[5] == feature[6] == locus['strand']
                    and interval[3] == locus['locus_id'], f'{mode}: export identities disagree')
            attrs = dict(item.split('=', 1) for item in feature[8].split(';') if '=' in item)
            require(unquote(attrs.get('ID', '')) == locus['locus_id'], f'{mode}: GFF3 ID disagrees')
            require(int(interval[1]) == int(locus['start0']) and int(interval[2]) == int(locus['end0']),
                    f'{mode}: BED/TSV coordinates disagree')
            require(int(feature[3]) == int(locus['start0']) + 1 and int(feature[4]) == int(locus['end0']),
                    f'{mode}: GFF3/TSV coordinates disagree')
        parser = Links()
        parser.feed((run / 'report.html').read_text())
        for target in parser.targets:
            parsed = urlparse(target)
            if parsed.scheme or parsed.netloc or not parsed.path:
                continue
            require((run / unquote(parsed.path)).exists(), f'{mode}: broken report link {target}')
        summary = json.loads((run / 'results/summary.json').read_text())
        require(summary['final_loci'] == len(final), f'{mode}: summary count disagrees')
        checked[mode] = dict(final_loci=len(final), report=str(run / 'report.html'))
    result = dict(status='passed', scope='Synthetic installation check; not a biological accuracy benchmark',
                  runs=checked)
    viewer = output / 'cross/report.html'
    require(viewer.is_file(), 'Merged candidate report is missing')
    match = re.search(r'<script id="final-alignment-data" type="application/json">(.*?)</script>', viewer.read_text(), re.S)
    require(match is not None, 'Merged report has no alignment data')
    data_rows = json.loads(match.group(1))
    require(len(data_rows) == 1 and data_rows[0]['transcript_id'] == 'demo_A.t1',
            'Final viewer must contain only the final demo_A.t1 candidate')
    final_row = rows(output / 'cross/results/final_loci.tsv')[0]
    record = data_rows[0]
    require(record['stop_count'] == int(final_row['stop_count'])
            and record['frameshift_count'] == int(final_row['frameshift_count']),
            'Viewer and TSV disruption counts disagree')
    require(abs(record['identity'] / 100 - float(final_row['alignment_identity'])) < 0.0001,
            'Viewer and TSV identity disagree')
    require({m['type'] for m in record['mutations']} == {'stop', 'frameshift'},
            'Final viewer must highlight both engineered events')
    result['merged_report'] = str(viewer)
    (output / 'demo_check.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = check(args.output.resolve())
    except (ValueError, OSError, KeyError, TypeError, IndexError) as error:
        if args.output.is_dir():
            (args.output / 'demo_check.json').write_text(json.dumps(
                dict(status='failed', error=str(error)), indent=2) + '\n')
        raise SystemExit(f'DEMO FAIL: {error}')
    print('DEMO PASS: Cross/Self completed; expected candidates, exports and report links verified.')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
