#!/usr/bin/env python3
"""
Visualize GeneWise pseudogene identification results.
Generates an interactive HTML report with alignments, mutations, and statistics.

Usage:
    python3 visualize_genewise.py \
        --genewise pseudogene-genewise.txt \
        --output visualization.html \
        [--pseudogene-list high_confidence_pseudogene_ids.txt]
"""
import argparse
import re
import sys
from collections import defaultdict
from html import escape


def parse_genewise_block(block):
    """Parse a single GeneWise output block."""
    result = {
        'transcript_id': None,
        'query_seq': '',
        'target_seq': '',
        'mutations': [],  # List of (position, type, description)
        'score': None,
        'stop_count': 0,
        'frameshift_count': 0,
        'alignment_lines': []
    }
    
    # Extract transcript ID
    id_match = re.search(r'Query protein:\s+(\S+)', block)
    if id_match:
        result['transcript_id'] = id_match.group(1)
    
    # Extract score
    score_match = re.search(r'Score\s+([0-9.]+)\s+bits', block)
    if score_match:
        result['score'] = float(score_match.group(1))
    
    # Count mutations
    result['stop_count'] = len(re.findall(r'\*', block))
    stopx_count = len(re.findall(r':X\[', block))
    if result['stop_count'] == 0 and stopx_count > 0:
        result['stop_count'] = stopx_count
    result['frameshift_count'] = len(re.findall(r'!', block))
    
    # Extract alignment sequences
    lines = block.split('\n')
    qid_short = result['transcript_id'].split('.')[0] if result['transcript_id'] else None
    
    if qid_short:
        qseq_parts = []
        tseq_parts = []
        in_alignment = False
        
        for i, line in enumerate(lines):
            if line.startswith(qid_short):
                parts = line.split()
                if len(parts) >= 3:
                    qseg = parts[2]
                    # Target sequence is typically 2 lines below
                    if i + 2 < len(lines):
                        tline = lines[i + 2]
                        # Clean target sequence (keep amino acids, X, !, *)
                        tseg = ''.join([c for c in tline if c.isalpha() or c in ('X', '!', '*', '-')])
                        qseq_parts.append(qseg)
                        tseq_parts.append(tseg)
                        in_alignment = True
                    elif in_alignment:
                        # Continue collecting if we're in alignment region
                        tline = lines[i + 1] if i + 1 < len(lines) else ''
                        tseg = ''.join([c for c in tline if c.isalpha() or c in ('X', '!', '*', '-')])
                        if tseg:
                            tseq_parts.append(tseg)
        
        if qseq_parts:
            result['query_seq'] = ''.join(qseq_parts)
        if tseq_parts:
            result['target_seq'] = ''.join(tseq_parts)
    
    # Calculate identity if sequences are present
    if result['query_seq'] and result['target_seq']:
        matches = sum(1 for a, b in zip(result['query_seq'], result['target_seq']) if a == b)
        length = len(result['query_seq'])
        result['identity'] = (matches / length * 100) if length > 0 else 0
        result['length'] = length
    else:
        result['identity'] = 0
        result['length'] = 0
    
    # Find mutation positions
    if result['query_seq'] and result['target_seq']:
        # Find stop codons (*)
        for i, char in enumerate(result['target_seq']):
            if char == '*':
                pos = i + 1
                result['mutations'].append({
                    'position': pos,
                    'type': 'stop',
                    'char': '*',
                    'description': f'Stop codon at position {pos}'
                })
        
        # Find frameshifts (!)
        for i, char in enumerate(result['target_seq']):
            if char == '!':
                pos = i + 1
                result['mutations'].append({
                    'position': pos,
                    'type': 'frameshift',
                    'char': '!',
                    'description': f'Frameshift at position {pos}'
                })
        
        # Find :X[ pattern (stop at splice site) - approximate position
        if ':X[' in block:
            # Try to find approximate position in alignment
            result['mutations'].append({
                'position': None,
                'type': 'splice_stop',
                'char': 'X',
                'description': 'Stop codon at splice site'
            })
    
    return result


def parse_genewise_file(genewise_file):
    """Parse entire GeneWise output file."""
    try:
        with open(genewise_file, 'r') as fh:
            content = fh.read()
    except IOError as e:
        print(f"ERROR: Cannot read {genewise_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    blocks = content.split("//")
    results = []
    
    for block in blocks:
        if not block.strip():
            continue
        result = parse_genewise_block(block)
        if result['transcript_id']:
            results.append(result)
    
    return results


def filter_by_list(results, pseudogene_list_file):
    """Filter results to only include transcripts in the list."""
    if not pseudogene_list_file:
        return results
    
    try:
        with open(pseudogene_list_file, 'r') as fh:
            allowed_ids = {line.strip().split()[0] for line in fh if line.strip()}
    except IOError as e:
        print(f"WARNING: Cannot read {pseudogene_list_file}: {e}", file=sys.stderr)
        return results
    
    return [r for r in results if r['transcript_id'] in allowed_ids]


def format_alignment_html(query_seq, target_seq, mutations, max_width=80):
    """Format alignment with highlighted mutations for HTML."""
    html_parts = []
    
    # Group mutations by position for quick lookup
    mutation_map = {}
    for mut in mutations:
        if mut['position']:
            mutation_map[mut['position']] = mut
    
    # Split into chunks
    for start in range(0, len(query_seq), max_width):
        end = min(start + max_width, len(query_seq))
        q_chunk = query_seq[start:end]
        t_chunk = target_seq[start:end]
        
        # Format query sequence
        html_parts.append('<div class="alignment-block">')
        html_parts.append(f'<div class="seq-line"><span class="seq-label">Query:</span> <span class="seq-query">')
        
        for i, char in enumerate(q_chunk):
            pos = start + i + 1
            html_parts.append(escape(char))
        html_parts.append('</span></div>')
        
        # Format target sequence with mutations highlighted
        html_parts.append('<div class="seq-line"><span class="seq-label">Target:</span> <span class="seq-target">')
        
        for i, char in enumerate(t_chunk):
            pos = start + i + 1
            mut = mutation_map.get(pos)
            
            if mut:
                if mut['type'] == 'stop':
                    html_parts.append(f'<span class="mutation stop" title="{escape(mut["description"])}">{escape(char)}</span>')
                elif mut['type'] == 'frameshift':
                    html_parts.append(f'<span class="mutation frameshift" title="{escape(mut["description"])}">{escape(char)}</span>')
                else:
                    html_parts.append(escape(char))
            else:
                html_parts.append(escape(char))
        
        html_parts.append('</span></div>')
        html_parts.append('</div>')
    
    return ''.join(html_parts)


def generate_html_report(results, output_file):
    """Generate HTML visualization report."""
    import json
    
    # Calculate statistics
    total = len(results)
    with_stops = sum(1 for r in results if r['stop_count'] > 0)
    with_frameshifts = sum(1 for r in results if r['frameshift_count'] > 0)
    with_both = sum(1 for r in results if r['stop_count'] > 0 and r['frameshift_count'] > 0)
    
    # Prepare data for charts
    chart_data = {
        'scores': [r['score'] for r in results if r['score'] is not None],
        'stops': [r['stop_count'] for r in results],
        'frameshifts': [r['frameshift_count'] for r in results],
        'identities': [r.get('identity', 0) for r in results],
        'lengths': [r.get('length', 0) for r in results],
        'scatter': [{'x': r['score'], 'y': r.get('identity', 0), 'id': r['transcript_id']} for r in results if r['score'] is not None],
        'lenId': [{'x': r.get('length', 0), 'y': r.get('identity', 0), 'id': r['transcript_id']} for r in results]
    }
    
    avg_score = sum(chart_data['scores']) / len(chart_data['scores']) if chart_data['scores'] else None
    


    # Serialize all results to JSON for client-side rendering
    # Only keep necessary fields to reduce size
    json_results = []
    for r in results:
        json_results.append({
            'transcript_id': r['transcript_id'],
            'score': r['score'],
            'stop_count': r['stop_count'],
            'frameshift_count': r['frameshift_count'],
            'identity': r.get('identity', 0),
            'length': r.get('length', 0),
            'query_seq': r['query_seq'],
            'target_seq': r['target_seq'],
            'mutations': r['mutations']
        })
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GeneWise Pseudogene Visualization</title>
    <!-- Chart.js -->
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background-color: #f5f5f5;
            padding: 20px;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        h1 {{
            color: #2c3e50;
            margin-bottom: 10px;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        
        h2 {{
            color: #34495e;
            margin-top: 30px;
            margin-bottom: 15px;
            padding-left: 10px;
            border-left: 4px solid #3498db;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        
        .stat-card h3 {{
            font-size: 2.5em;
            margin-bottom: 5px;
        }}
        
        .stat-card p {{
            font-size: 0.9em;
            opacity: 0.9;
        }}

        .charts-container {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}

        .chart-box {{
            background: white;
            padding: 15px;
            border: 1px solid #ddd;
            border-radius: 8px;
            height: 300px;
        }}
        
        /* Table Styles */
        .controls {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
        }}

        .search-box {{
            padding: 10px;
            font-size: 1em;
            border: 2px solid #ddd;
            border-radius: 4px;
            width: 300px;
        }}

        .pseudogene-list {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.95em;
        }}

        .pseudogene-list th, .pseudogene-list td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}

        .pseudogene-list th {{
            background-color: #f8f9fa;
            font-weight: bold;
            color: #2c3e50;
            cursor: pointer;
            user-select: none;
        }}
        
        .pseudogene-list th:hover {{
            background-color: #e9ecef;
        }}

        .pseudogene-list tr:hover {{
            background-color: #f1f1f1;
        }}

        .btn-expand {{
            background: #3498db;
            color: white;
            border: none;
            padding: 5px 10px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 0.85em;
        }}

        .btn-expand:hover {{
            background: #2980b9;
        }}

        .alignment-row {{
            background-color: #fafafa;
            display: none;
        }}

        .alignment-row.show {{
            display: table-row;
        }}
        
        .alignment-cell {{
            padding: 20px !important;
        }}

        /* Pagination */
        .pagination {{
            display: flex;
            justify-content: center;
            align-items: center;
            gap: 15px;
            margin-top: 20px;
        }}

        .page-btn {{
            padding: 8px 15px;
            background: #fff;
            border: 1px solid #ddd;
            border-radius: 4px;
            cursor: pointer;
            transition: all 0.2s;
        }}

        .page-btn:hover:not(:disabled) {{
            background: #f1f1f1;
            border-color: #ccc;
        }}

        .page-btn:disabled {{
            opacity: 0.5;
            cursor: not-allowed;
        }}
        
        /* Badges */
        .badge {{
            padding: 4px 8px;
            border-radius: 12px;
            font-size: 0.8em;
            font-weight: bold;
            color: white;
            margin-right: 5px;
            display: inline-block;
        }}
        
        .badge-stop {{ background: #e74c3c; }}
        .badge-frameshift {{ background: #f39c12; }}
        .badge-score {{ background: #3498db; }}
        .badge-identity {{ background: #27ae60; }}

        /* Alignment Detail Styling */
        .alignment-detail h4 {{
            margin-bottom: 10px;
            color: #2c3e50;
        }}

        .alignment-block {{
            margin: 10px 0;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
            background: #fff;
            border: 1px solid #eee;
            padding: 10px;
            border-radius: 4px;
            overflow-x: auto;
        }}
        
        .seq-line {{
            margin: 2px 0;
            white-space: nowrap;
        }}
        
        .seq-label {{
            display: inline-block;
            width: 80px;
            font-weight: bold;
            color: #7f8c8d;
        }}
        
        .seq-query {{ color: #27ae60; }}
        .seq-target {{ color: #2980b9; }}
        
        .mutation {{
            font-weight: bold;
            padding: 0 2px;
            border-radius: 2px;
            cursor: help;
            color: white;
        }}
        
        .mutation.stop {{ background: #e74c3c; }}
        .mutation.frameshift {{ background: #f39c12; }}
        
        .mutations-list {{
            margin-bottom: 15px;
            padding: 10px;
            background: #fff3cd;
            border-left: 4px solid #f39c12;
            border-radius: 4px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 GeneWise Pseudogene Identification Results</h1>
        
        <div class="stats-grid">
            <div class="stat-card">
                <h3>{total}</h3>
                <p>Total Pseudogenes</p>
            </div>
            <div class="stat-card">
                <h3>{with_stops}</h3>
                <p>With Stop Codons</p>
            </div>
            <div class="stat-card">
                <h3>{with_frameshifts}</h3>
                <p>With Frameshifts</p>
            </div>
            <div class="stat-card">
                <h3>{with_both}</h3>
                <p>With Both Mutations</p>
            </div>
            {f'<div class="stat-card"><h3>{avg_score:.1f}</h3><p>Average Score</p></div>' if avg_score else ''}
        </div>
        
        <h2>Overview Charts</h2>
        <div class="charts-container">
            <div class="chart-box">
                <canvas id="scoreChart"></canvas>
            </div>
            <div class="chart-box">
                <canvas id="mechanismChart"></canvas>
            </div>
            <div class="chart-box">
                <canvas id="mutationChart"></canvas>
            </div>
            <div class="chart-box">
                <canvas id="identityChart"></canvas>
            </div>
            <div class="chart-box">
                <canvas id="scatterChart"></canvas>
            </div>
            <div class="chart-box">
                <canvas id="lenIdChart"></canvas>
            </div>
             <div class="chart-box">
                <canvas id="lengthChart"></canvas>
            </div>
        </div>

        <h2>Pseudogene Details</h2>
        
        <div class="controls">
            <input type="text" class="search-box" id="searchBox" placeholder="Search by transcript ID...">
            <div id="recordCount">Showing 0 results</div>
        </div>
        
        <table class="pseudogene-list">
            <thead>
                <tr>
                    <th width="50"></th>
                    <th width="200">Transcript ID</th>
                    <th>Score</th>
                    <th>Identity</th>
                    <th>Mutations</th>
                    <th>Length</th>
                </tr>
            </thead>
            <tbody id="tableBody">
                <!-- Rows injected by JS -->
            </tbody>
        </table>

        <div class="pagination">
            <button class="page-btn" id="prevBtn" onclick="changePage(-1)">Previous</button>
            <span id="pageInfo">Page 1</span>
            <button class="page-btn" id="nextBtn" onclick="changePage(1)">Next</button>
        </div>
    </div>
    
    <script>
        // Data passed from Python
        const chartData = {json.dumps(chart_data)};
        // Full Dataset
        const allPseudogenes = {json.dumps(json_results)};
        
        // --- State ---
        let currentPage = 1;
        const itemsPerPage = 50;
        let filteredData = [...allPseudogenes];

        // --- Initialization ---
        document.addEventListener('DOMContentLoaded', () => {{
            renderCharts();
            updateTable();
            
            // Search Listener
            document.getElementById('searchBox').addEventListener('input', (e) => {{
                const term = e.target.value.toLowerCase();
                filteredData = allPseudogenes.filter(item => 
                    item.transcript_id.toLowerCase().includes(term)
                );
                currentPage = 1;
                updateTable();
            }});
        }});

        // --- Table Rendering ---
        function updateTable() {{
            const tbody = document.getElementById('tableBody');
            tbody.innerHTML = '';
            
            const start = (currentPage - 1) * itemsPerPage;
            const end = start + itemsPerPage;
            const pageData = filteredData.slice(start, end);
            
            document.getElementById('recordCount').textContent = `Showing ${{filteredData.length}} results`;
            document.getElementById('pageInfo').textContent = `Page ${{currentPage}} of ${{Math.ceil(filteredData.length / itemsPerPage)}}`;
            document.getElementById('prevBtn').disabled = currentPage === 1;
            document.getElementById('nextBtn').disabled = end >= filteredData.length || filteredData.length === 0;

            if (pageData.length === 0) {{
                tbody.innerHTML = '<tr><td colspan="6" style="text-align:center; padding: 20px;">No results found</td></tr>';
                return;
            }}

            pageData.forEach((item, index) => {{
                const realIndex = start + index;
                // Main Row
                const tr = document.createElement('tr');
                tr.innerHTML = `
                    <td><button class="btn-expand" onclick="toggleDetails(${{realIndex}})">+</button></td>
                    <td style="font-family: monospace; font-weight: bold;">${{item.transcript_id}}</td>
                    <td><span class="badge badge-score">${{item.score ? item.score.toFixed(1) : 'N/A'}}</span></td>
                    <td><span class="badge badge-identity">${{item.identity.toFixed(1)}}%</span></td>
                    <td>
                        ${{item.stop_count > 0 ? `<span class="badge badge-stop">${{item.stop_count}} Stop</span>` : ''}}
                        ${{item.frameshift_count > 0 ? `<span class="badge badge-frameshift">${{item.frameshift_count}} Frame</span>` : ''}}
                        ${{item.stop_count === 0 && item.frameshift_count === 0 ? '<span style="color:#aaa">None</span>' : ''}}
                    </td>
                    <td>${{item.length}} aa</td>
                `;
                tbody.appendChild(tr);

                // Detail Row (Hidden)
                const detailTr = document.createElement('tr');
                detailTr.id = `detail-${{realIndex}}`;
                detailTr.className = 'alignment-row';
                detailTr.innerHTML = `
                    <td colspan="6" class="alignment-cell">
                        <div id="alignment-content-${{realIndex}}">Loading...</div>
                    </td>
                `;
                tbody.appendChild(detailTr);
            }});
        }}

        function changePage(delta) {{
            currentPage += delta;
            updateTable();
            // window.scrollTo({{ top: document.querySelector('.controls').offsetTop - 20, behavior: 'smooth' }});
        }}

        function toggleDetails(index) {{
            const row = document.getElementById(`detail-${{index}}`);
            const isHidden = !row.classList.contains('show');
            
            if (isHidden) {{
                row.classList.add('show');
                const container = document.getElementById(`alignment-content-${{index}}`);
                renderAlignment(container, filteredData[index]);
            }} else {{
                row.classList.remove('show');
            }}
        }}

        function renderAlignment(container, item) {{
            if (!item.query_seq || !item.target_seq) {{
                container.innerHTML = '<p>No alignment available</p>';
                return;
            }}

            let html = '<div class="alignment-detail">';
            
            // Mutation List
            if (item.mutations && item.mutations.length > 0) {{
                html += '<div class="mutations-list"><strong>Mutations:</strong><ul>';
                item.mutations.forEach(mut => {{
                     html += `<li>${{mut.description}}</li>`;
                }});
                html += '</ul></div>';
            }}

            html += '<h4>Sequence Alignment</h4>';
            
            // Alignment Blocks
            const maxWidth = 80;
            const querySeq = item.query_seq;
            const targetSeq = item.target_seq;
            
            // Index mutations for fast lookup
            const mutMap = {{}};
            item.mutations.forEach(m => {{ if(m.position) mutMap[m.position] = m; }});

            for (let start = 0; start < querySeq.length; start += maxWidth) {{
                const end = Math.min(start + maxWidth, querySeq.length);
                const qChunk = querySeq.slice(start, end);
                const tChunk = targetSeq.slice(start, end);

                html += '<div class="alignment-block">';
                
                // Query Line
                html += '<div class="seq-line"><span class="seq-label">Query:</span> <span class="seq-query">';
                html += qChunk.split('').map(c => `<span>${{c}}</span>`).join('');
                html += '</span></div>';

                // Target Line
                html += '<div class="seq-line"><span class="seq-label">Target:</span> <span class="seq-target">';
                
                tChunk.split('').forEach((char, i) => {{
                    const pos = start + i + 1;
                    const mut = mutMap[pos];
                    if (mut) {{
                        const className = mut.type === 'stop' ? 'stop' : (mut.type === 'frameshift' ? 'frameshift' : '');
                        html += `<span class="mutation ${{className}}" title="${{mut.description}}">${{char}}</span>`;
                    }} else {{
                        html += `<span>${{char}}</span>`;
                    }}
                }});
                
                html += '</span></div></div>';
            }}
            
            html += '</div>';
            container.innerHTML = html;
        }}

        // --- Helper Functions ---
        function createHistogram(data, binCount) {{
            if(!data || data.length === 0) return {{ bins: [], labels: [] }};
            const maxVal = Math.max(...data);
            const minVal = Math.min(...data);
            const range = maxVal - minVal;
            if (range === 0) return {{ bins: [data.length], labels: [`${{minVal}}`] }};

            const binSize = range / binCount;
            const bins = new Array(binCount).fill(0);
            const labels = [];
            
            for(let i=0; i<binCount; i++) {{
                const start = minVal + i * binSize;
                const end = start + binSize;
                labels.push(`${{start.toFixed(1)}}-${{end.toFixed(1)}}`);
            }}
            
            data.forEach(val => {{
                let binIndex = Math.floor((val - minVal) / binSize);
                if(binIndex >= binCount) binIndex = binCount - 1;
                if(binIndex < 0) binIndex = 0;
                bins[binIndex]++;
            }});
            
            return {{ bins, labels }};
        }}
        
        // --- Chart Renderers ---
        function renderCharts() {{
            // 1. Score Distribution
            const scoreHist = createHistogram(chartData.scores, 10);
            new Chart(document.getElementById('scoreChart'), {{
                type: 'bar',
                data: {{
                    labels: scoreHist.labels,
                    datasets: [{{
                        label: 'Score Distribution',
                        data: scoreHist.bins,
                        backgroundColor: 'rgba(52, 152, 219, 0.6)',
                        borderColor: 'rgba(52, 152, 219, 1)',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{ title: {{ display: true, text: 'Pseudogene Score Distribution' }} }}
                }}
            }});
            
            // 2. Mechanism (Pie Chart) - NEW
            // Calculate totals from stops and frameshifts arrays
            const totalStops = chartData.stops ? chartData.stops.reduce((a, b) => a + b, 0) : 0;
            const totalShifts = chartData.frameshifts ? chartData.frameshifts.reduce((a, b) => a + b, 0) : 0;
            
            new Chart(document.getElementById('mechanismChart'), {{
                type: 'doughnut',
                data: {{
                    labels: ['Stop Codons', 'Frameshifts'],
                    datasets: [{{
                        data: [totalStops, totalShifts],
                        backgroundColor: ['#e74c3c', '#f39c12'],
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{ 
                        title: {{ display: true, text: 'Mechanism of Inactivation (Total Count)' }},
                        subtitle: {{ display: true, text: `Stops: ${{totalStops}}, Frameshifts: ${{totalShifts}}` }}
                    }}
                }}
            }});
            
            // 3. Mutation Counts (Grouped Bar)
            function countFrequencies(arr) {{
                const counts = {{}};
                arr.forEach(x => counts[x] = (counts[x] || 0) + 1);
                const max = Math.max(0, ...Object.keys(counts).map(Number)); // Ensure at least 0
                const labels = Array.from({{length: max + 1}}, (_, i) => i);
                const data = labels.map(i => counts[i] || 0);
                return {{ labels, data }};
            }}
            
            const stopFreq = countFrequencies(chartData.stops);
            const shiftFreq = countFrequencies(chartData.frameshifts);
            const maxLen = Math.max(stopFreq.labels.length, shiftFreq.labels.length);
            const labels = Array.from({{length: maxLen}}, (_, i) => i);
            
            new Chart(document.getElementById('mutationChart'), {{
                type: 'bar',
                data: {{
                    labels: labels,
                    datasets: [
                        {{ label: 'Stop Codons', data: labels.map(i => stopFreq.data[i] || 0), backgroundColor: 'rgba(231, 76, 60, 0.6)' }},
                        {{ label: 'Frameshifts', data: labels.map(i => shiftFreq.data[i] || 0), backgroundColor: 'rgba(243, 156, 18, 0.6)' }}
                    ]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{ title: {{ display: true, text: 'Mutation Counts per Pseudogene' }} }},
                    scales: {{
                        x: {{ title: {{ display: true, text: 'Number of Mutations' }} }},
                        y: {{ title: {{ display: true, text: 'Count of Pseudogenes' }} }}
                    }}
                }}
            }});
            
            // 4. Identity
            const identityHist = createHistogram(chartData.identities, 10);
            new Chart(document.getElementById('identityChart'), {{
                type: 'bar',
                data: {{
                    labels: identityHist.labels,
                    datasets: [{{
                        label: 'Sequence Identity (%)',
                        data: identityHist.bins,
                        backgroundColor: 'rgba(39, 174, 96, 0.6)',
                        borderColor: 'rgba(39, 174, 96, 1)',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{ title: {{ display: true, text: 'Sequence Identity Distribution' }} }}
                }}
            }});
            
            // 5. Scatter: Identity vs Score
            new Chart(document.getElementById('scatterChart'), {{
                type: 'scatter',
                data: {{
                    datasets: [{{
                        label: 'Identity vs Score',
                        data: chartData.scatter,
                        backgroundColor: 'rgba(155, 89, 182, 0.6)'
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        title: {{ display: true, text: 'Identity vs. Score Correlation' }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    return context.raw.id + ': (' + context.raw.x + ', ' + context.raw.y.toFixed(1) + '%)';
                                }}
                            }}
                        }}
                    }},
                    scales: {{
                        x: {{ title: {{ display: true, text: 'Score' }} }},
                        y: {{ title: {{ display: true, text: 'Identity (%)' }} }}
                    }}
                }}
            }});

            // 6. Scatter: Length vs Identity - NEW
            new Chart(document.getElementById('lenIdChart'), {{
                type: 'scatter',
                data: {{
                    datasets: [{{
                        label: 'Length vs Identity',
                        data: chartData.lenId,
                        backgroundColor: 'rgba(22, 160, 133, 0.6)'
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        title: {{ display: true, text: 'Length vs Identity' }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    return context.raw.id + ': (' + context.raw.x + ' aa, ' + context.raw.y.toFixed(1) + '%)';
                                }}
                            }}
                        }}
                    }},
                    scales: {{
                        x: {{ title: {{ display: true, text: 'Alignment Length (aa)' }} }},
                        y: {{ title: {{ display: true, text: 'Identity (%)' }} }}
                    }}
                }}
            }});

             // 7. Length Distribution
            const lengthHist = createHistogram(chartData.lengths, 10);
            new Chart(document.getElementById('lengthChart'), {{
                type: 'bar',
                data: {{
                    labels: lengthHist.labels,
                    datasets: [{{
                        label: 'Alignment Length',
                        data: lengthHist.bins,
                        backgroundColor: 'rgba(52, 73, 94, 0.6)',
                        borderColor: 'rgba(52, 73, 94, 1)',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{ title: {{ display: true, text: 'Alignment Length Distribution' }} }}
                }}
            }});
        }}
    </script>
</body>
</html>
"""
    
    try:
        with open(output_file, 'w') as fh:
            fh.write(html)
        print(f"Visualization saved to: {output_file}", file=sys.stderr)
    except IOError as e:
        print(f"ERROR: Cannot write {output_file}: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser(
        description="Visualize GeneWise pseudogene identification results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Visualize all results
  python3 visualize_genewise.py \\
      --genewise output/species/genewise/pseudogene-genewise.txt \\
      --output visualization.html
  
  # Visualize only high-confidence pseudogenes
  python3 visualize_genewise.py \\
      --genewise output/species/genewise/pseudogene-genewise.txt \\
      --output visualization.html \\
      --pseudogene-list output/species/post_genewise_filter/high_conf_strict_pseudogene_ids.txt
        """
    )
    ap.add_argument('--genewise', required=True,
                     help='GeneWise output file (pseudogene-genewise.txt)')
    ap.add_argument('--output', required=True,
                     help='Output HTML file')
    ap.add_argument('--pseudogene-list',
                     help='Optional: Filter to only these transcript IDs (one per line)')
    
    args = ap.parse_args()
    
    print("Parsing GeneWise output...", file=sys.stderr)
    results = parse_genewise_file(args.genewise)
    print(f"Found {len(results)} pseudogene results", file=sys.stderr)
    
    if args.pseudogene_list:
        print(f"Filtering by pseudogene list...", file=sys.stderr)
        results = filter_by_list(results, args.pseudogene_list)
        print(f"After filtering: {len(results)} pseudogenes", file=sys.stderr)
    
    if not results:
        print("WARNING: No results to visualize", file=sys.stderr)
        sys.exit(1)
    
    print("Generating HTML report...", file=sys.stderr)
    generate_html_report(results, args.output)
    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
