import React, { useEffect, useRef } from 'react';
import { stateLabels } from './api.js';
import { DnaLogo } from './DnaAnimation.jsx';

export function Badge({ state }) {
  return <span className={`badge ${Object.hasOwn(stateLabels, state) ? state : ''}`}>{stateLabels[state] || 'Unavailable'}</span>;
}
export function Sidebar({ motion }) {
  return <aside className="sidebar">
    <a className="brand" href="#overview" aria-label="GeneVISTA home"><DnaLogo motion={motion} />Gene<span>VISTA</span></a>
    <div className="brand-caption">Genomic Variant Intelligence</div>
    <nav aria-label="Main navigation"><a href="#explorer">Variant search</a><a href="#predictions">Disease evidence</a><a href="#disease-ranking">Disease model</a><a href="#model-charts">Performance charts</a><a href="#molecular">Molecular studio</a><a href="#settings">Settings</a></nav>
  </aside>;
}
export function Stats({ data }) {
  const cards = [
    ['◈', 'Retained variant records', data?.variant_count, 'Latest processed snapshot'],
    ['↗', 'Variant timelines', data?.timeline_count, 'Across all prepared snapshots'],
    ['▦', 'Historical snapshots', data?.snapshots.length, 'January 2022 – September 2026'],
    ['⊙', 'Gene-symbol links', data?.gene_symbol_link_count, 'Latest processed snapshot'],
  ];
  return <section className="stats" aria-label="Dataset overview">{cards.map(([icon, title, value, note]) =>
    <article key={title}><span className="stat-icon" aria-hidden="true">{icon}</span><p>{title}</p><strong>{value?.toLocaleString() ?? '—'}</strong><small>{note}</small></article>
  )}</section>;
}
export function Results({ result, loading, onMore, onSelect }) {
  if (!result?.items.length) return null;
  return <section className="panel" aria-label="Search results">
    <div className="section-heading"><div><p className="eyebrow">SEARCH RESULTS</p><h2>{result.search_type === 'variant_id' ? `Variant ${result.query}` : `${result.query.toUpperCase()} variants`}</h2></div><span className="muted">{result.items.length.toLocaleString()} shown</span></div>
    <div className="table-scroll"><table><thead><tr><th>VARIANT ID</th><th>GENE</th><th>VARIANT NAME</th><th>GROUPED STATE</th><th><span className="sr-only">Action</span></th></tr></thead>
      <tbody>{result.items.map(item => <tr key={item.variation_id}>
        <td className="id">{item.variation_id}</td><td>{item.latest_details?.gene_symbol || '—'}</td>
        <td className="name">{item.latest_details?.name || 'Latest retained details unavailable; timeline available.'}</td>
        <td><Badge state={item.latest_details?.clinical_state || item.timeline?.at(-1)?.state} /></td>
        <td><button className="view-button" aria-label={`View variant ${item.variation_id}`} onClick={() => onSelect(item.variation_id)}>View →</button></td>
      </tr>)}</tbody></table></div>
    <div className="table-footer"><span className="muted">Gene search includes retained latest-snapshot records.</span>{result.next_after_id !== null && <button className="secondary" disabled={loading} onClick={onMore}>{loading ? 'Loading…' : 'Load more'}</button>}</div>
  </section>;
}
export function Detail({ selected, onClose, onRetry }) {
  const sectionRef = useRef(null);
  const headingRef = useRef(null);
  useEffect(() => {
    sectionRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    headingRef.current?.focus({ preventScroll: true });
  }, [selected?.id]);
  if (!selected) return null;
  const details = selected.data?.latest_details;
  const fields = details ? [
    ['Gene symbol', details.gene_symbol], ['Variant type', details.variant_type],
    ['Recorded classification', details.clinical_significance], ['Review status', details.review_status],
    ['Submitters', details.number_submitters_raw], ['Source snapshot', details.snapshot],
  ] : [['Latest snapshot state', stateLabels[selected.data?.timeline.at(-1)?.state]]];
  return <section className="panel details" ref={sectionRef} aria-label="Selected variant">
    <div className="section-heading"><div><p className="eyebrow">VARIANT PROFILE</p><h2 ref={headingRef} tabIndex={-1}>Variant {selected.id}</h2></div><button className="secondary" aria-label="Close variant details" onClick={onClose}>Close ×</button></div>
    {selected.loading && <p className="detail-status" role="status">Loading the variant and its history…</p>}
    {selected.error && <div className="detail-status" role="alert">{selected.error} <button className="secondary" onClick={onRetry}>Retry</button></div>}
    {selected.data && <div className="detail-content">
      <p className="variant-name">{details?.name || 'This variant has historical observations but no retained latest-snapshot details.'}</p>
      <dl className="detail-fields">{fields.map(([label, value]) => <div key={label}><dt>{label}</dt><dd>{value || 'Unavailable'}</dd></div>)}</dl>
      <div className="timeline-heading"><h3>Classification timeline</h3><span>Observed snapshot states</span></div>
      <div className="timeline" aria-label="Six snapshot observations">{selected.data.timeline.map(observation =>
        <div className="timeline-item" key={observation.snapshot}><time dateTime={observation.snapshot}>{observation.snapshot}</time><Badge state={observation.state} /></div>
      )}</div>
      <p className="footnote">Benign and pathogenic groups include their “likely” classifications. These are recorded classifications, not model predictions. An absent GRCh38 observation does not establish removal from ClinVar.</p>
    </div>}
  </section>;
}
export function About() {
  return <section id="about" className="about-grid">
    <article className="about-card"><span className="eyebrow">A LONGITUDINAL VIEW</span><h3>One variant. Multiple moments.</h3><p>Compare snapshots without treating an unavailable classification as a biological change.</p><div className="mini-timeline" aria-hidden="true"><span>2022</span><i /><span>2023</span><i /><span>2024</span><i /><span>2025</span><i /><span>2026</span></div></article>
    <article className="about-card"><span className="eyebrow">READ THE EVIDENCE IN CONTEXT</span><h3>Clear provenance, visible limits.</h3><p>Latest details include retained GRCh38 records. ClinGen curations describe gene–disease evidence. “Present, excluded” means a record failed preparation rules. Educational molecular illustrations do not establish variant-specific functional effects.</p></article>
  </section>;
}
