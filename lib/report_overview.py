"""Compact GeneWise overview measurements for the unified offline report."""
import csv
import math
import re
from pathlib import Path

from run_layout import work_directory


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def overview_records(run_dir, finals, alignments):
    run_dir = Path(run_dir)
    path = run_dir / 'results/loci.tsv'
    if path.is_file():
        with path.open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
    else:
        rows = finals
    final_ids = {r['locus_id'] for r in finals}
    records = {}
    for row in rows:
        status = row.get('genewise_status')
        if status is not None and status != 'completed':
            continue
        if status is None and row.get('evidence_method') != 'GeneWise':
            continue
        locus = row['locus_id']
        identity = number(row.get('alignment_identity'))
        record = dict(locus_id=locus, final=locus in final_ids, score=None,
                      stops=number(row.get('stop_count', row.get('genewise_stops'))),
                      frameshifts=number(row.get('frameshift_count', row.get('genewise_frameshifts'))),
                      identity=identity * 100 if identity is not None else None,
                      length=number(row.get('aligned_length')))
        parsed = alignments.get(locus)
        if parsed:
            record.update(score=parsed['score'], stops=parsed['stop_count'],
                          frameshifts=parsed['frameshift_count'])
            if parsed['query_seq'] and parsed['target_seq']:
                record.update(identity=parsed['identity'], length=parsed['length'])
        records[locus] = record
    work = work_directory(run_dir)
    # Read scores once from combined outputs, without embedding all sequences.
    for path in (work / 'genewise/pseudogene-genewise.txt', run_dir / 'genewise_alignments.txt',
                 run_dir / 'genewise_candidate.txt'):
        if not path.is_file():
            continue
        locus = None
        with path.open() as handle:
            for line in handle:
                if line.strip() == '//':
                    locus = None
                query = re.search(r'Query protein:\s+(\S+)', line)
                if query:
                    locus = query.group(1)
                if locus in records:
                    score = re.search(r'Score\s+([0-9.]+)\s+bits', line)
                    if score:
                        records[locus]['score'] = number(score.group(1))
    for locus, record in records.items():
        if record['score'] is not None or Path(locus).name != locus or locus in ('.', '..'):
            continue
        path = work / 'genewise' / (locus + '.gw.txt')
        if path.is_file():
            with path.open() as handle:
                for line in handle:
                    score = re.search(r'Score\s+([0-9.]+)\s+bits', line)
                    if score:
                        record['score'] = number(score.group(1))
                        break
    return list(records.values())


def overview_section(records):
    if not records:
        return '<section id="overview"><h2>GeneWise overview</h2><p class="empty-state">No completed GeneWise measurements are available for this run. Missing measurements are not counted as zero.</p></section>'
    charts = [('score', 'GeneWise score distribution', 'Score (bits)', 'Alignments'),
              ('events', 'Disruption events', '', ''),
              ('mutations', 'Mutation counts per alignment', 'Events per alignment', 'Alignments'),
              ('identity', 'Sequence identity distribution', 'Amino-acid identity (%)', 'Alignments'),
              ('identity-score', 'Identity versus score', 'Score (bits)', 'Identity (%)'),
              ('length-identity', 'Length versus identity', 'Alignment length (aa)', 'Identity (%)'),
              ('length', 'Alignment length distribution', 'Alignment length (aa)', 'Alignments')]
    panels = ''.join(f'<div class="overview-chart"><h3>{title}</h3><svg id="overview-{key}" viewBox="0 0 400 270" role="img" aria-label="{title}"></svg><p class="chart-caption" id="overview-{key}-caption"></p></div>'
                     for key, title, _, _ in charts)
    labels = [('alignments', 'Alignment records'), ('stops', 'With stop codons'),
              ('frameshifts', 'With frameshifts'), ('both', 'With both mutations'), ('score', 'Average score (bits)')]
    cards = ''.join(f'<div><span>{label}</span><strong id="overview-stat-{key}">—</strong></div>' for key, label in labels)
    return f'''<section id="overview"><h2>GeneWise overview</h2>
<label for="overview-scope">Statistics scope</label> <select id="overview-scope"><option value="all">All completed GeneWise alignments</option><option value="final">Final candidate alignments</option></select>
<p id="overview-scope-note" aria-live="polite">All completed alignments include loci excluded by downstream filters. They are not the final pseudogene count. Candidate-table search and pagination do not change these statistics.</p>
<div class="overview-stats">{cards}</div><div class="overview-charts">{panels}</div>
<noscript><p>Enable JavaScript to display charts. Measurements remain available in the result tables.</p></noscript></section>'''
