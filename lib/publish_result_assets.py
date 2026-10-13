"""Write common GFF3 and an offline report with final-candidate alignments."""
import csv
import html
import json
from collections import Counter
from pathlib import Path
from urllib.parse import quote

from run_layout import work_directory
from visualize_genewise import format_alignment_html, parse_genewise_block
from report_overview import overview_records, overview_section


def final_alignments(run_dir, finals):
    """Match GeneWise evidence by locus ID, never by the parent query ID."""
    run_dir = Path(run_dir)
    wanted = {row['locus_id'] for row in finals}
    records = {}
    work = work_directory(run_dir)
    for locus in wanted:
        # Refuse path traversal in an externally supplied result table.
        if Path(locus).name != locus or locus in ('.', '..'):
            continue
        path = work / 'genewise' / (locus + '.gw.txt')
        if path.is_file():
            record = parse_genewise_block(path.read_text())
            if record['transcript_id'] == locus:
                records[locus] = record
    # Portable snapshots contain one raw candidate file; older runs may only
    # have the combined output. Stream it so rejected loci are not embedded.
    for path in (run_dir / 'genewise_candidate.txt', work / 'genewise/pseudogene-genewise.txt'):
        if wanted <= records.keys() or not path.is_file():
            continue
        with path.open() as handle:
            block = []
            for line in handle:
                block.append(line)
                if line.strip() == '//':
                    record = parse_genewise_block(''.join(block))
                    locus = record['transcript_id']
                    if locus in wanted and locus not in records:
                        records[locus] = record
                    block = []
            if block:
                record = parse_genewise_block(''.join(block))
                if record['transcript_id'] in wanted and record['transcript_id'] not in records:
                    records[record['transcript_id']] = record
    return records


def candidate_details(finals, alignments):
    esc = lambda value: html.escape(str(value), quote=True)

    def percent(value):
        try:
            return f'{float(value) * 100:.2f}%'
        except (ValueError, TypeError):
            return 'Not assessed'

    if not finals:
        return '<p class="empty-state">No final candidates passed the enabled filters.</p>'
    rows = []
    for index, row in enumerate(finals):
        locus = row['locus_id']
        record = alignments.get(locus)
        evidence = ''.join(f'<div><span>{esc(label)}</span><strong>{esc(value)}</strong></div>'
                           for label, value in [('Stop codons', row.get('stop_count', '')),
                                                ('Frameshifts', row.get('frameshift_count', '')),
                                                ('Query coverage', percent(row.get('query_coverage'))),
                                                ('Alignment identity', percent(row.get('alignment_identity')))])
        if record and record['query_seq'] and record['target_seq']:
            alignment = format_alignment_html(record['query_seq'], record['target_seq'], record['mutations'])
            events = html_table([('type', 'Event'), ('position', 'Alignment column'), ('description', 'Description')],
                                record['mutations'])
            body = f'<h3>GeneWise amino-acid alignment</h3><p class="legend"><b class="stop-key">* Stop codon</b> <b class="frame-key">! Frameshift</b> · Positions refer to alignment columns.</p><div class="alignment">{alignment}</div><h3>Disruption events</h3>{events}'
        else:
            body = '<p class="empty-state">GeneWise alignment is not available for this locus. Consult the recorded evidence method and result tables.</p>'
        search_text = ' '.join(str(row.get(k, '')) for k in ('locus_id', 'query_id', 'chrom', 'evidence_level')).lower()
        hidden = ' hidden' if index >= 20 else ''
        cells = ''.join('<td>' + value + '</td>' for value in [
            f'<button type="button" class="candidate-toggle" aria-expanded="false" aria-controls="candidate-detail-{index}" aria-label="Show alignment and evidence for {esc(locus)}">{esc(locus)}</button>',
            esc(row['query_id']),
            f"{esc(row['chrom'])}:{esc(row['start0'])}–{esc(row['end0'])} ({esc(row['strand'])})",
            f'<span class="grade">{esc(row.get("evidence_level", "unassessed"))}</span>',
            esc(row.get('stop_count', '')), esc(row.get('frameshift_count', '')),
            esc(percent(row.get('query_coverage'))), esc(percent(row.get('alignment_identity')))])
        rows.append(f'''<tr class="candidate-row" data-search="{esc(search_text)}"{hidden}>{cells}</tr>
<tr class="candidate-detail" id="candidate-detail-{index}" hidden><td colspan="8"><div class="candidate-body"><h3>{esc(locus)}</h3><p>Reference protein: <strong>{esc(row['query_id'])}</strong> · Evidence method: {esc(row.get('evidence_method', ''))}</p><div class="evidence-metrics">{evidence}</div><p class="decision">Filter decision: {esc(row.get('decision_reason', ''))}</p>{body}</div></td></tr>''')
    headings = ''.join('<th scope="col">' + label + '</th>' for label in
                       ['Candidate locus', 'Reference protein', 'Location', 'Evidence', 'Stops', 'Frameshifts', 'Coverage', 'Identity'])
    pages = (len(finals) + 19) // 20
    next_disabled = ' disabled' if pages == 1 else ''
    return f'''<label class="search-label" for="candidate-search">Find a candidate</label><input id="candidate-search" type="search" placeholder="Locus, protein, contig or evidence level">
<p id="search-status" aria-live="polite">Showing 1–{min(20, len(finals))} of {len(finals)} candidates · Page 1 of {pages}</p>
<div class="table-wrap"><table class="candidate-table"><thead><tr>{headings}</tr></thead><tbody>{''.join(rows)}<tr id="no-candidate-matches" hidden><td colspan="8">No candidates match your search.</td></tr></tbody></table></div>
<nav class="pagination" aria-label="Candidate pages"><button id="candidate-prev" type="button" disabled>Previous</button><span id="candidate-page">Page 1 of {pages}</span><button id="candidate-next" type="button"{next_disabled}>Next</button><span>20 candidates per page</span></nav>
<noscript><p>Enable JavaScript to search, change pages and expand alignments. Complete records are available in the downloadable result tables.</p></noscript>'''


def write_text(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.partial')
    temporary.write_text(content, encoding='utf-8')
    temporary.replace(path)


def table(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def gff3(rows):
    lines = ['##gff-version 3', '# Computational pseudogene candidates; evidence grades are not biological confidence probabilities.']
    seen = set()
    for row in rows:
        locus = row['locus_id']
        if not locus or locus in seen:
            raise ValueError('Empty or duplicate final locus ID')
        seen.add(locus)
        start0, end0 = int(row['start0']), int(row['end0'])
        if start0 < 0 or end0 <= start0 or row['strand'] not in ('+', '-', '.'):
            raise ValueError('Invalid final locus coordinates/strand: ' + locus)
        attrs = dict(ID=locus, Name=locus, reference_protein=row['query_id'])
        for key in ['evidence_query_id', 'evidence_level', 'evidence_method', 'stop_count',
                    'frameshift_count', 'query_coverage', 'alignment_identity', 'identity_metric',
                    'decision_reason']:
            if row.get(key) not in ('', None):
                attrs[key] = row[key]
        encoded = ';'.join(key + '=' + quote(str(value), safe='._:-') for key, value in attrs.items())
        lines.append('\t'.join([quote(row['chrom'], safe='._:-'), 'EasyPseudogene', 'pseudogene',
                               str(start0 + 1), str(end0), '.', row['strand'], '.', encoded]))
    return '\n'.join(lines) + '\n'


def html_table(fields, rows):
    header = ''.join('<th>' + html.escape(label) + '</th>' for _, label in fields)
    body = ''.join('<tr>' + ''.join('<td>' + html.escape(str(row.get(key, ''))) + '</td>'
                                  for key, _ in fields) + '</tr>' for row in rows)
    if not rows:
        body = f'<tr><td colspan="{len(fields)}">No records</td></tr>'
    return '<div class="table-wrap"><table><thead><tr>' + header + '</tr></thead><tbody>' + body + '</tbody></table></div>'


def render_report(run_dir, manifest, summary, finals):
    run_dir = Path(run_dir)
    esc = lambda value: html.escape(str(value), quote=True)
    mode = manifest.get('comparison_mode', '')
    workflow = manifest.get('workflow', 'standard')
    workflow_label = 'GeneWise workflow' if workflow == 'standard' else 'Catalog with structural evidence'
    note = ('Final candidates passed the standard GeneWise workflow, including the enabled family, copy and alignment-quality filters.'
            if workflow == 'standard' else
            'Final candidates require a stop or frameshift detected by miniprot structural validation and medium/high sequence evidence. The standard GeneWise family and copy filters were not applied.')
    if manifest.get('layout_preview'):
        note = 'Illustrative rows demonstrate table layout and pagination only. They are not analysis results.'
    total = summary.get('input_total_records', summary.get('input_queries', 0))
    cards = [('Input protein records', total), ('Accepted query proteins', summary.get('input_queries', 0)),
             ('Excluded input records', summary.get('input_excluded_records', 0)),
             ('Selected loci', summary.get('selected_loci', 0)), ('Final candidate loci', summary.get('final_loci', 0)),
             ('Final candidate queries', summary.get('final_queries', 0))]
    metrics = ''.join(f'<div class="card"><span>{esc(label)}</span><strong>{esc(value)}</strong></div>' for label, value in cards)
    counts = summary.get('evidence_levels', {})
    final_counts = Counter(r.get('evidence_level', 'unassessed') for r in finals)
    grades = [dict(level=level, all=counts.get(level, 0), final=final_counts.get(level, 0), rule=rule)
              for level, rule in [('high', 'Stop/frameshift; query coverage ≥0.8 and identity ≥0.5'),
                                  ('medium', 'Stop/frameshift; query coverage ≥0.2 and identity ≥0.3'),
                                  ('low', 'Assessed, with no disruption or below the evidence thresholds'),
                                  ('unassessed', 'Structural validation or required measurements are missing')]]
    levels = html_table([('level', 'Evidence level'), ('all', 'All selected loci'), ('final', 'Final candidates'), ('rule', 'Rule')], grades)
    stages_path = run_dir / 'qc/stage_counts.tsv'
    stages = html_table([('stage', 'Stage'), ('loci_or_inputs', 'Loci / inputs'), ('queries', 'Query proteins')],
                        table(stages_path) if stages_path.is_file() else [])
    alignments = final_alignments(run_dir, finals) if workflow == 'standard' else {}
    candidates = candidate_details(finals, alignments)
    candidate_json = json.dumps(list(alignments.values()), ensure_ascii=True).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    overview = overview_records(run_dir, finals, alignments) if workflow == 'standard' else []
    overview_html = overview_section(overview)
    overview_json = json.dumps(overview, ensure_ascii=True, allow_nan=False).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    overview_script = Path(__file__).with_name('report_overview.js').read_text()
    downloads = [('results/final_loci.tsv', 'Final candidates (TSV)'), ('results/final_loci.bed', 'Final candidates (BED)'),
                 ('results/final_loci.gff3', 'Final candidates (GFF3)'), ('results/final_queries.txt', 'Candidate query protein IDs'),
                 ('results/loci.tsv', 'All selected loci and filtering decisions'), ('results/queries.tsv', 'All accepted query proteins'),
                 ('results/summary.json', 'Summary statistics (JSON)'), ('qc/input_qc.json', 'Input quality control'),
                 ('qc/protein_exclusions.tsv', 'Excluded protein records'), ('qc/nonfinal_loci.tsv', 'Loci outside the final candidate set'),
                 ('metadata/evidence_status.tsv', 'Assessed and unassessed evidence'),
                 ('metadata/run_manifest.json', 'Run metadata'), ('metadata/effective_config.txt', 'Effective parameters'),
                 ('logs/run.log', 'Run log')]
    links = ''.join(f'<li><a href="{esc(path)}">{esc(label)}</a></li>' for path, label in downloads if (run_dir / path).is_file())
    qc_path = run_dir / 'qc/input_qc.json'
    qc = json.loads(qc_path.read_text()) if qc_path.is_file() else {}
    warnings = qc.get('annotation_warnings', [])
    qc_text = 'Input QC status: ' + esc(qc.get('status', 'not provided'))
    if warnings:
        qc_text += '. Annotation notes: ' + esc('; '.join(str(v) for v in warnings[:10]))
    scope = '<p class="scope-note">Synthetic installation demo · This example does not measure biological accuracy.</p>' if manifest.get('snapshot_only') else ''
    if manifest.get('layout_preview'):
        scope = '<p class="scope-note">Layout preview · 100 illustrative rows, not pseudogene discoveries.</p>'
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>EasyPseudogene — {esc(manifest.get('species', run_dir.name))}</title>
<style>body{{margin:0;background:#f5f7fa;color:#182333;font:15px/1.65 system-ui,sans-serif}}main{{max-width:1180px;margin:auto;padding:28px 20px}}h1,h2{{line-height:1.3}}.subtitle{{color:#506176}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px}}.card,section{{background:white;border:1px solid #dae1e9;border-radius:10px;padding:18px;margin-bottom:18px}}.card span{{display:block;color:#506176}}.card strong{{font-size:30px}}.table-wrap{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{text-align:left;padding:9px;border-bottom:1px solid #e3e8ef;white-space:nowrap}}th{{background:#f5f7fa}}a{{color:#195aa3}}ul{{padding-left:22px}}.note{{border-left:4px solid #6290b6;padding-left:14px}}footer{{color:#506176}}
header{{margin-bottom:25px}}h1{{font-size:32px;margin:8px 0}}h2{{font-size:22px}}.eyebrow{{color:#14766b;letter-spacing:2px;font-size:11px;font-weight:700}}.scope-note{{font-size:13px;color:#986316;background:#fff4da;padding:8px 12px;border-radius:6px}}nav{{display:flex;gap:20px;flex-wrap:wrap;margin-top:20px}}nav a{{text-decoration:none;font-weight:600}}.cards{{margin-bottom:20px}}.card{{margin-bottom:0;box-shadow:0 3px 10px #18233304}}.card:nth-child(5){{background:#0f6e65;color:white;border-color:#0f6e65}}.card:nth-child(5) span{{color:#d0f1e9}}.count-pill,.grade{{display:inline-block;border-radius:20px;background:#e0f3ec;color:#125d4d;padding:2px 10px;font-size:13px;vertical-align:middle}}.count-pill{{margin-left:8px}}.candidate-table{{font-size:13px}}.candidate-table th{{position:sticky;top:0}}.candidate-row:hover{{background:#f6faf9}}.candidate-toggle{{border:0;background:transparent;color:#195aa3;padding:4px 0;font:inherit;font-weight:600;cursor:pointer;text-align:left}}.candidate-toggle::before{{content:'▸';display:inline-block;width:18px}}.candidate-toggle[aria-expanded=true]::before{{content:'▾'}}.candidate-detail>td{{white-space:normal;background:#f8fafc;padding:0}}.candidate-body{{padding:16px 20px}}.candidate-body h3{{margin-top:0}}.pagination{{display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin-top:18px}}.pagination button{{border:1px solid #cad6e0;border-radius:6px;background:white;padding:8px 14px;color:#195aa3;font:inherit;cursor:pointer}}.pagination button:disabled{{color:#7e8995;background:#f5f7fa;cursor:default}}.pagination span{{font-size:13px;color:#506176}}button:focus-visible{{outline:2px solid #14766b;outline-offset:3px}}.evidence-metrics{{display:grid;grid-template-columns:repeat(auto-fit,minmax(120px,1fr));gap:10px}}.evidence-metrics>div{{background:#f5f7fa;border-radius:6px;padding:10px 14px}}.evidence-metrics span{{display:block;font-size:12px;color:#506176}}.evidence-metrics strong{{font-size:22px}}.decision{{font-size:12px;color:#506176}}.alignment{{overflow:auto;background:#f8fafc;border:1px solid #e1e7ed;border-radius:6px;padding:16px;font:13px/1.7 ui-monospace,SFMono-Regular,Consolas,monospace}}.alignment-block{{margin-bottom:12px}}.seq-line{{white-space:pre}}.seq-label{{display:inline-block;width:65px;color:#68788b}}.mutation{{font-weight:bold;border-radius:3px}}.stop,.stop-key{{background:#ffe0e5;color:#b41938}}.frameshift,.frame-key{{background:#ffe7c2;color:#945000}}.legend{{font-size:12px;color:#506176}}.legend b{{padding:3px 7px;border-radius:4px}}.search-label{{display:block;font-size:12px;font-weight:600;color:#506176;margin-bottom:5px}}input[type=search]{{box-sizing:border-box;width:100%;max-width:420px;padding:10px 12px;font:inherit;border:1px solid #cad6e0;border-radius:6px}}#search-status{{font-size:12px;color:#506176}}.diagnostics{{border-top:1px solid #e3e8ef;padding:12px 0}}.diagnostics summary{{cursor:pointer;font-weight:600}}.diagnostics .table-wrap{{margin-top:14px}}.downloads{{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:8px 20px}}.empty-state{{color:#506176;background:#f5f7fa;border-radius:6px;padding:16px}}html{{scroll-behavior:smooth}}summary:focus-visible,a:focus-visible,input:focus-visible{{outline:2px solid #14766b;outline-offset:3px}}
.overview-stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin:20px 0}}.overview-stats>div{{background:#f2f6fa;border-radius:8px;padding:14px}}.overview-stats span{{display:block;font-size:12px;color:#506176}}.overview-stats strong{{display:block;font-size:26px}}.overview-charts{{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,350px),1fr));gap:18px}}.overview-chart{{min-width:0;border:1px solid #dae1e9;border-radius:8px;padding:12px}}.overview-chart h3{{font-size:14px;margin:0 0 10px}}.overview-chart svg{{width:100%;height:auto;display:block}}.chart-caption{{font-size:12px;color:#506176;margin:8px 0 0}}#overview-scope{{font:inherit;border:1px solid #cad6e0;border-radius:6px;padding:7px;max-width:100%;background:white;color:#182333}}#overview-scope-note{{font-size:13px;color:#506176}}
</style></head>
<body><main><header><p class="eyebrow">EASYPSEUDOGENE / ANALYSIS RESULTS</p><h1>EasyPseudogene Analysis Report</h1>
<p class="subtitle">{esc(manifest.get('species', run_dir.name))} · {esc(mode.capitalize())} · {esc(workflow_label)} · Analysis version: {esc(manifest.get('analysis_version', ''))}</p>
{scope}<nav><a href="#overview">Overview charts</a><a href="#candidates">Final candidates</a><a href="#diagnostics">Filtering &amp; QC</a><a href="#downloads">Download results</a></nav></header>
<div class="cards">{metrics}</div>
{overview_html}
<section id="candidates"><h2>Final candidates <span class="count-pill">{len(finals)}</span></h2><p>Click a candidate ID to expand its alignment and disruption evidence. Showing 20 candidates per page; search covers all candidates. Only loci in final_loci.tsv appear here. Coordinates are 0-based, half-open.</p>{candidates}</section>
<section id="diagnostics"><h2>Filtering and quality control</h2><details class="diagnostics"><summary>Filtering stage counts</summary>{stages}</details><details class="diagnostics"><summary>Evidence level distribution</summary>{levels}</details><details class="diagnostics"><summary>Result interpretation and input QC</summary><p>{esc(note)}</p><p class="note">These computational candidates require biological review. High/medium/low grades describe sequence-disruption evidence, not the probability that a locus is a true pseudogene. Evidence grades and filtering decisions are recorded separately. GeneWise identity is amino-acid identity; miniprot identity is the PAF match proportion. These measurements are not directly equivalent. See evidence_status.tsv for assessed and unassessed evidence.</p><p>{qc_text}</p></details></section>
<section id="downloads"><h2>Results and quality-control files</h2><ul class="downloads">{links}</ul></section>
<footer>This report opens offline, including alignments and search. TSV/BED coordinates are 0-based, half-open; GFF3 coordinates are 1-based, inclusive.</footer></main>
<script id="final-alignment-data" type="application/json">{candidate_json}</script>
<script id="overview-data" type="application/json">{overview_json}</script>
<script>{overview_script}</script>
<script>
const search = document.getElementById('candidate-search');
if (search) {{
  const rows = Array.from(document.querySelectorAll('.candidate-row'));
  const previous = document.getElementById('candidate-prev');
  const next = document.getElementById('candidate-next');
  let page = 0;
  const pageSize = 20;
  const renderCandidates = () => {{
    const term = search.value.trim().toLowerCase();
    const matches = rows.filter(row => row.dataset.search.includes(term));
    const pages = Math.ceil(matches.length / pageSize);
    page = Math.max(0, Math.min(page, pages - 1));
    const visible = new Set(matches.slice(page * pageSize, (page + 1) * pageSize));
    rows.forEach(row => {{
      row.hidden = !visible.has(row);
      const button = row.querySelector('.candidate-toggle');
      const detail = document.getElementById(button.getAttribute('aria-controls'));
      detail.hidden = row.hidden || button.getAttribute('aria-expanded') !== 'true';
    }});
    const first = matches.length ? page * pageSize + 1 : 0;
    const last = Math.min((page + 1) * pageSize, matches.length);
    document.getElementById('search-status').textContent = `Showing ${{first}}–${{last}} of ${{matches.length}} candidates (${{rows.length}} total) · Page ${{pages ? page + 1 : 0}} of ${{pages}}`;
    document.getElementById('candidate-page').textContent = `Page ${{pages ? page + 1 : 0}} of ${{pages}}`;
    document.getElementById('no-candidate-matches').hidden = matches.length !== 0;
    previous.disabled = page === 0;
    next.disabled = !pages || page >= pages - 1;
  }};
  search.addEventListener('input', () => {{ page = 0; renderCandidates(); }});
  previous.addEventListener('click', () => {{ page--; renderCandidates(); }});
  next.addEventListener('click', () => {{ page++; renderCandidates(); }});
  rows.forEach(row => {{
    const button = row.querySelector('.candidate-toggle');
    button.addEventListener('click', () => {{
      button.setAttribute('aria-expanded', button.getAttribute('aria-expanded') === 'true' ? 'false' : 'true');
      renderCandidates();
    }});
  }});
  renderCandidates();
}}
</script></body></html>
'''


def publish_assets(run_dir):
    run_dir = Path(run_dir)
    finals = table(run_dir / 'results/final_loci.tsv')
    summary = json.loads((run_dir / 'results/summary.json').read_text())
    manifest_path = run_dir / 'metadata/run_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if summary['final_loci'] != len(finals):
        raise ValueError('Final locus table and summary counts differ')
    write_text(run_dir / 'results/final_loci.gff3', gff3(finals))
    work = work_directory(run_dir)
    manifest['layout_version'] = 4 if work != run_dir else 3
    manifest['files'].update(report='../report.html', final_gff3='../results/final_loci.gff3')
    if work != run_dir:
        manifest['files'].update(work_directory='../work/', reference_cache='reference_cache.json', logs='../logs/')
    if (run_dir / 'results/genewise_alignment.html').is_file():
        manifest['files']['genewise_alignment'] = '../results/genewise_alignment.html'
    write_text(manifest_path, json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    write_text(run_dir / 'report.html', render_report(run_dir, manifest, summary, finals))
