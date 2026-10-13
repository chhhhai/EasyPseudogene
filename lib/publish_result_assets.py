"""Write common GFF3 and a standalone HTML summary from published result tables."""
import csv
import html
import json
from collections import Counter
from pathlib import Path
from urllib.parse import quote

from run_layout import work_directory


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
        body = f'<tr><td colspan="{len(fields)}">无记录</td></tr>'
    return '<div class="table-wrap"><table><thead><tr>' + header + '</tr></thead><tbody>' + body + '</tbody></table></div>'


def render_report(run_dir, manifest, summary, finals):
    esc = lambda value: html.escape(str(value), quote=True)
    mode = manifest.get('comparison_mode', '')
    workflow = manifest.get('workflow', 'standard')
    workflow_label = '标准（GeneWise）' if workflow == 'standard' else '快速（结构证据筛选）'
    note = ('最终候选经过标准 GeneWise 筛选，包括已启用的家族、拷贝和比对质量过滤。'
            if workflow == 'standard' else
            '最终候选要求 miniprot 结构验证检测到终止或移码，并达到 medium/high 序列证据等级；未应用标准流程的 GeneWise 家族/拷贝过滤。')
    total = summary.get('input_total_records', summary.get('input_queries', 0))
    cards = [('输入蛋白记录', total), ('可用查询蛋白', summary.get('input_queries', 0)),
             ('排除的输入记录', summary.get('input_excluded_records', 0)),
             ('初筛位点', summary.get('selected_loci', 0)), ('最终候选位点', summary.get('final_loci', 0)),
             ('最终候选查询蛋白', summary.get('final_queries', 0))]
    metrics = ''.join(f'<div class="card"><span>{esc(label)}</span><strong>{esc(value)}</strong></div>' for label, value in cards)
    counts = summary.get('evidence_levels', {})
    final_counts = Counter(r.get('evidence_level', 'unassessed') for r in finals)
    grades = [dict(level=level, all=counts.get(level, 0), final=final_counts.get(level, 0), rule=rule)
              for level, rule in [('high', '有终止/移码；覆盖率 ≥0.8、相似度 ≥0.5'),
                                  ('medium', '有终止/移码；覆盖率 ≥0.2、相似度 ≥0.3'),
                                  ('low', '已评估，无破坏或低于证据阈值'),
                                  ('unassessed', '结构验证或必要测量缺失')]]
    levels = html_table([('level', '证据等级'), ('all', '全部位点'), ('final', '最终候选'), ('rule', '规则')], grades)
    stages_path = run_dir / 'qc/stage_counts.tsv'
    stages = html_table([('stage', '筛选阶段'), ('loci_or_inputs', '位点/输入数'), ('queries', '查询蛋白数')],
                        table(stages_path) if stages_path.is_file() else [])
    preview = html_table([('locus_id', '候选位点'), ('query_id', '参考蛋白'), ('chrom', '序列'),
                          ('start0', '起点（0-based）'), ('end0', '终点'), ('strand', '方向'),
                          ('evidence_level', '证据等级'), ('stop_count', '终止数'), ('frameshift_count', '移码数')], finals[:50])
    downloads = [('results/final_loci.tsv', '最终候选 TSV'), ('results/final_loci.bed', '最终候选 BED'),
                 ('results/final_loci.gff3', '最终候选 GFF3'), ('results/final_queries.txt', '候选参考蛋白列表'),
                 ('results/loci.tsv', '全部位点与筛选原因'), ('results/queries.tsv', '全部查询蛋白'),
                 ('results/summary.json', '统计 JSON'), ('qc/input_qc.json', '输入质控'),
                 ('qc/protein_exclusions.tsv', '排除的蛋白记录'), ('qc/nonfinal_loci.tsv', '未进入最终集的位点'),
                 ('metadata/evidence_status.tsv', '已评估及未评估的证据'),
                 ('metadata/run_manifest.json', '运行信息'), ('metadata/effective_config.txt', '有效参数'),
                 ('logs/run.log', '运行日志'), ('results/genewise_alignment.html', 'GeneWise 比对浏览')]
    links = ''.join(f'<li><a href="{esc(path)}">{esc(label)}</a></li>' for path, label in downloads if (run_dir / path).is_file())
    qc_path = run_dir / 'qc/input_qc.json'
    qc = json.loads(qc_path.read_text()) if qc_path.is_file() else {}
    warnings = qc.get('annotation_warnings', [])
    qc_text = '输入质控状态：' + esc(qc.get('status', '未提供'))
    if warnings:
        qc_text += '。注释提示：' + esc('; '.join(str(v) for v in warnings[:10]))
    return f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>EasyPseudogene — {esc(manifest.get('species', run_dir.name))}</title>
<style>body{{margin:0;background:#f5f7fa;color:#182333;font:15px/1.65 system-ui,sans-serif}}main{{max-width:1180px;margin:auto;padding:28px 20px}}h1,h2{{line-height:1.3}}.subtitle{{color:#506176}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px}}.card,section{{background:white;border:1px solid #dae1e9;border-radius:10px;padding:18px;margin-bottom:18px}}.card span{{display:block;color:#506176}}.card strong{{font-size:30px}}.table-wrap{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{text-align:left;padding:9px;border-bottom:1px solid #e3e8ef;white-space:nowrap}}th{{background:#f5f7fa}}a{{color:#195aa3}}ul{{padding-left:22px}}.note{{border-left:4px solid #6290b6;padding-left:14px}}footer{{color:#506176}}</style></head>
<body><main><h1>EasyPseudogene 分析报告</h1>
<p class="subtitle">{esc(manifest.get('species', run_dir.name))} · {esc(mode)} · {esc(workflow_label)} · {esc(manifest.get('analysis_version', ''))}</p>
<div class="cards">{metrics}</div>
<section><h2>结果说明</h2><p>{esc(note)}</p><p class="note">这些结果是需要生物学复核的计算候选。high/medium/low 描述序列破坏证据，不是真实假基因的概率；等级与筛选结果分别记录。GeneWise 相似度为氨基酸一致性，miniprot 相似度为 PAF 匹配比例，两者不能直接等同。各项证据是否已评估请查看 evidence_status.tsv。</p><p>{qc_text}</p></section>
<section><h2>证据等级分布</h2>{levels}</section>
<section><h2>筛选阶段统计</h2>{stages}</section>
<section><h2>最终候选预览</h2><p>显示前 {min(50, len(finals))} 个，共 {len(finals)} 个。TSV/BED 为 0-based 半开区间；GFF3 为 1-based 闭区间。</p>{preview}</section>
<section><h2>结果与质控文件</h2><ul>{links}</ul></section>
<footer>本报告可离线打开；完整候选、筛选原因及参数保存在旁边的结果文件中。</footer></main></body></html>
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
