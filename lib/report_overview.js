/* Offline GeneWise charts; all data are embedded in the report. */
(() => {
  const selector = document.getElementById('overview-scope');
  if (!selector) return;
  const records = JSON.parse(document.getElementById('overview-data').textContent);
  const ns = 'http://www.w3.org/2000/svg';
  const colors = { score: '#559dcc', stops: '#d44c60', frameshifts: '#dc962e', identity: '#409d79', scatter: '#9470b4', length: '#748b9e' };
  const finite = value => typeof value === 'number' && Number.isFinite(value);
  const fmt = value => value.toLocaleString('en-US', { maximumFractionDigits: 1 });
  function node(parent, tag, attrs = {}, label) {
    const el = document.createElementNS(ns, tag);
    Object.entries(attrs).forEach(([key, value]) => el.setAttribute(key, value));
    if (label !== undefined) el.textContent = label;
    parent.appendChild(el);
    return el;
  }
  function chart(key, description) {
    const svg = document.getElementById(`overview-${key}`);
    svg.replaceChildren();
    node(svg, 'title', {}, description);
    return svg;
  }
  function caption(key, text) { document.getElementById(`overview-${key}-caption`).textContent = text; }
  function empty(key) {
    const svg = chart(key, 'No measurements available');
    node(svg, 'text', { x: 200, y: 130, 'text-anchor': 'middle', fill: '#506176' }, 'No measurements available');
    caption(key, 'Missing measurements are excluded.');
  }
  function axes(svg, xmax, ymax, xlabel, ylabel, xmin = 0, ymin = 0, categories = null) {
    if (ylabel === 'Alignments') ymax = Math.max(4, Math.ceil(ymax / 4) * 4);
    const box = { left: 66, top: 20, width: 314, height: 190 };
    const x = value => box.left + (value - xmin) / (xmax - xmin) * box.width;
    const y = value => box.top + box.height - (value - ymin) / (ymax - ymin) * box.height;
    for (let i = 0; i <= 4; i++) {
      const xv = xmin + (xmax - xmin) * i / 4;
      const yv = ymin + (ymax - ymin) * i / 4;
      node(svg, 'line', { x1: box.left, x2: box.left + box.width, y1: y(yv), y2: y(yv), stroke: '#e3e8ef' });
      node(svg, 'text', { x: box.left - 7, y: y(yv) + 4, 'text-anchor': 'end', fill: '#506176', 'font-size': 11 }, fmt(yv));
      if (!categories) node(svg, 'text', { x: x(xv), y: 230, 'text-anchor': i === 0 ? 'start' : i === 4 ? 'end' : 'middle', fill: '#506176', 'font-size': 11 }, fmt(xv));
    }
    if (categories) {
      const stride = Math.max(1, Math.ceil(categories.length / 5));
      categories.forEach((label, index) => {
        if (index % stride === 0 || index === categories.length - 1) node(svg, 'text', { x: x(index + 0.5), y: 230, 'text-anchor': 'middle', fill: '#506176', 'font-size': 11 }, label);
      });
    }
    node(svg, 'path', { d: `M${box.left},${box.top}V210H380`, fill: 'none', stroke: '#93a1af' });
    node(svg, 'text', { x: 220, y: 258, 'text-anchor': 'middle', fill: '#182333', 'font-size': 12 }, xlabel);
    node(svg, 'text', { x: 14, y: 118, transform: 'rotate(-90 14 118)', 'text-anchor': 'middle', fill: '#182333', 'font-size': 12 }, ylabel);
    return { x, y, box };
  }
  function histogram(data, key, field, xlabel, color) {
    const values = data.map(record => record[field]).filter(finite);
    if (!values.length) return empty(key);
    let min = Infinity, max = -Infinity;
    values.forEach(value => { min = Math.min(min, value); max = Math.max(max, value); });
    if (min === max) { const pad = Math.max(1, Math.abs(min) * 0.02); min = Math.max(0, min - pad); max += pad; }
    if (field === 'identity') { min = Math.max(0, min); max = Math.min(100, max); }
    const bins = Array(10).fill(0);
    const width = (max - min) / bins.length;
    values.forEach(value => { bins[Math.min(9, Math.max(0, Math.floor((value - min) / width)))]++; });
    const svg = chart(key, `${values.length} measured alignments`);
    const a = axes(svg, max, Math.max(1, ...bins), xlabel, 'Alignments', min);
    bins.forEach((count, index) => {
      const bar = node(svg, 'rect', { x: a.x(min + index * width) + 1, y: a.y(count), width: a.box.width / 10 - 2, height: 210 - a.y(count), fill: color });
      node(bar, 'title', {}, `${fmt(min + index * width)}–${fmt(min + (index + 1) * width)}: ${count} alignments`);
    });
    caption(key, `${values.length} measured alignments; ${data.length - values.length} missing.`);
  }
  function donut(data) {
    const stops = data.reduce((sum, record) => sum + (finite(record.stops) ? record.stops : 0), 0);
    const frames = data.reduce((sum, record) => sum + (finite(record.frameshifts) ? record.frameshifts : 0), 0);
    const total = stops + frames;
    const svg = chart('events', `${stops} stop events and ${frames} frameshift events`);
    const radius = 70, circumference = 2 * Math.PI * radius;
    node(svg, 'circle', { cx: 200, cy: 110, r: radius, fill: 'none', stroke: '#e9eef3', 'stroke-width': 32 });
    if (total) {
      let offset = 0;
      [[stops, colors.stops, 'Stop codons'], [frames, colors.frameshifts, 'Frameshifts']].forEach(([count, color, label]) => {
        if (count) {
          const arc = node(svg, 'circle', { cx: 200, cy: 110, r: radius, fill: 'none', stroke: color, 'stroke-width': 32, 'stroke-dasharray': `${count / total * circumference} ${circumference}`, 'stroke-dashoffset': -offset, transform: 'rotate(-90 200 110)' });
          node(arc, 'title', {}, `${label}: ${count} events`);
        }
        offset += count / total * circumference;
      });
    }
    node(svg, 'text', { x: 200, y: 108, 'text-anchor': 'middle', fill: '#182333', 'font-size': 24 }, fmt(total));
    node(svg, 'text', { x: 200, y: 130, 'text-anchor': 'middle', fill: '#506176', 'font-size': 12 }, 'events');
    [[stops, colors.stops, 'Stops'], [frames, colors.frameshifts, 'Frameshifts']].forEach(([count, color, label], index) => {
      node(svg, 'rect', { x: 75, y: 213 + index * 22, width: 12, height: 12, fill: color });
      node(svg, 'text', { x: 95, y: 224 + index * 22, fill: '#182333', 'font-size': 12 }, `${label}: ${fmt(count)}`);
    });
    const missing = data.filter(record => !finite(record.stops) || !finite(record.frameshifts)).length;
    caption('events', `Counts are events, not loci. ${missing} alignment(s) have missing event measurements.`);
  }
  function mutations(data) {
    const measured = data.filter(record => finite(record.stops) || finite(record.frameshifts));
    if (!measured.length) return empty('mutations');
    const bins = Array.from({ length: 11 }, () => [0, 0]);
    measured.forEach(record => [record.stops, record.frameshifts].forEach((count, series) => {
      if (finite(count)) bins[Math.min(10, Math.max(0, Math.round(count)))][series]++;
    }));
    let last = 0;
    bins.forEach((counts, index) => { if (counts.some(count => count > 0)) last = index; });
    const size = Math.max(2, last + 1);
    const svg = chart('mutations', 'Distribution of stop and frameshift counts per alignment');
    const a = axes(svg, size, Math.max(1, ...bins.flat()), 'Events per alignment', 'Alignments', 0, 0,
      Array.from({ length: size }, (_, index) => index === 10 ? '10+' : String(index)));
    bins.slice(0, size).forEach((counts, index) => {
      counts.forEach((count, series) => {
        const bar = node(svg, 'rect', { x: a.x(index) + 2 + series * a.box.width / size * 0.42, y: a.y(count), width: a.box.width / size * 0.38, height: 210 - a.y(count), fill: series ? colors.frameshifts : colors.stops });
        node(bar, 'title', {}, `${series ? 'Frameshifts' : 'Stops'}, ${index === 10 ? '10+' : index} events: ${count} alignments`);
      });
    });
    caption('mutations', 'Red: stops · Orange: frameshifts. Missing measurements excluded; 10+ events grouped.');
  }
  function scatter(data, key, xfield, xlabel) {
    const points = data.filter(record => finite(record[xfield]) && finite(record.identity));
    if (!points.length) return empty(key);
    let xmax = 1;
    points.forEach(record => { xmax = Math.max(xmax, record[xfield]); });
    const svg = chart(key, `${points.length} measured alignments`);
    const a = axes(svg, xmax * 1.05, 100, xlabel, 'Identity (%)');
    const step = Math.max(1, Math.ceil(points.length / 2000));
    let shown = 0;
    points.forEach((record, index) => {
      if (index % step) return;
      const point = node(svg, 'circle', { cx: a.x(record[xfield]), cy: a.y(record.identity), r: 3, fill: key === 'identity-score' ? colors.scatter : colors.identity, opacity: 0.65 });
      node(point, 'title', {}, `${record.locus_id}: ${xlabel} ${fmt(record[xfield])}; identity ${fmt(record.identity)}%`);
      shown++;
    });
    caption(key, `${points.length} measured alignments; ${data.length - points.length} missing.` + (shown < points.length ? ` Showing ${shown} representative points.` : ''));
  }
  function render() {
    const data = records.filter(record => selector.value !== 'final' || record.final);
    const count = key => data.filter(record => finite(record[key]) && record[key] > 0).length;
    const scores = data.map(record => record.score).filter(finite);
    const values = { alignments: data.length, stops: count('stops'), frameshifts: count('frameshifts'), both: data.filter(record => record.stops > 0 && record.frameshifts > 0).length,
      score: scores.length ? fmt(scores.reduce((sum, value) => sum + value, 0) / scores.length) : 'Not assessed' };
    Object.entries(values).forEach(([key, value]) => { document.getElementById(`overview-stat-${key}`).textContent = value; });
    document.getElementById('overview-scope-note').textContent = `${data.length} measured alignment record(s). ` + (selector.value === 'final' ? 'Only final candidate loci are included.' : 'Includes loci excluded by downstream filters; this is not the final pseudogene count.') + ` Average score uses ${scores.length} available scores. Table search and pagination do not change this overview.`;
    histogram(data, 'score', 'score', 'Score (bits)', colors.score);
    donut(data);
    mutations(data);
    histogram(data, 'identity', 'identity', 'Amino-acid identity (%)', colors.identity);
    scatter(data, 'identity-score', 'score', 'Score (bits)');
    scatter(data, 'length-identity', 'length', 'Alignment length (aa)');
    histogram(data, 'length', 'length', 'Alignment length (aa)', colors.length);
  }
  selector.addEventListener('change', render);
  render();
})();
