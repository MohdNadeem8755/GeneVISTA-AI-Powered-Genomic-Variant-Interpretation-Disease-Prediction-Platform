import React, { useEffect, useRef, useState } from 'react';
import { fetchApi, errorMessage } from './api.js';
import { Sidebar, Stats, Results, Detail, About } from './components.jsx';
import DnaAnimation from './DnaAnimation.jsx';
import EvidencePanel from './EvidencePanel.jsx';
import MolecularJourney from './MolecularJourney.jsx';
import RankingPanel, {ModelCharts} from './RankingPanel.jsx';

export default function App() {
  const [theme,setTheme] = useState(()=>{try{return localStorage.getItem('genevista-theme')==='light'?'light':'dark';}catch{return 'dark';}});
  useEffect(()=>{document.documentElement.dataset.theme=theme;try{localStorage.setItem('genevista-theme',theme);}catch{}},[theme]);
  const [motion, setMotion] = useState(() => !window.matchMedia('(prefers-reduced-motion: reduce)').matches);
  const [stats, setStats] = useState(null);
  const [statsError, setStatsError] = useState('');
  const [statsAttempt, setStatsAttempt] = useState(0);
  const [query, setQuery] = useState('');
  const [result, setResult] = useState(null);
  const [searching, setSearching] = useState(false);
  const [message, setMessage] = useState('');
  const [selected, setSelected] = useState(null);
  const searchController = useRef(null);
  const detailController = useRef(null);
  const inputRef = useRef(null);

  useEffect(() => {
    const controller = new AbortController();
    setStatsError('');
    fetchApi('/api/stats', controller.signal).then(setStats).catch(error => {
      if (!controller.signal.aborted) setStatsError(errorMessage(error));
    });
    return () => controller.abort();
  }, [statsAttempt]);
  useEffect(() => () => { searchController.current?.abort(); detailController.current?.abort(); }, []);

  function closeDetail() {
    detailController.current?.abort(); setSelected(null); inputRef.current?.focus();
  }
  async function openDetail(id) {
    detailController.current?.abort();
    const controller = new AbortController(); detailController.current = controller;
    setSelected({ id, loading: true });
    try {
      const data = await fetchApi('/api/variants/' + encodeURIComponent(id), controller.signal);
      if (!controller.signal.aborted) setSelected({ id, data, loading: false });
    } catch (error) {
      if (!controller.signal.aborted) setSelected({ id, error: errorMessage(error), loading: false });
    }
  }
  async function runSearch(value, more = false) {
    const cleaned = value.trim();
    if (!cleaned) {setMessage('Enter an exact gene symbol or numeric variant ID.'); return;}
    searchController.current?.abort();
    const controller = new AbortController(); searchController.current = controller;
    if (!more) {setResult(null); detailController.current?.abort(); setSelected(null); setQuery(cleaned);}
    setSearching(true); setMessage(more ? 'Loading more variants…' : `Searching for ${cleaned}…`);
    try {
      const params = new URLSearchParams({ q: cleaned, limit: '20', after_id: String(more ? result.next_after_id : 0) });
      const data = await fetchApi('/api/variants/search?' + params, controller.signal);
      if (controller.signal.aborted) return;
      setResult(previous => more ? { ...data, items: [...previous.items, ...data.items] } : data);
      setMessage(data.items.length || more ? '' : `No results for “${cleaned}”${stats?.dataset_mode === 'demo' ? ' in the demo sample' : ''}. Try BRCA1 or variant 37565.`);
      if (!more && data.search_type === 'variant_id' && data.items.length) openDetail(data.items[0].variation_id);
    } catch (error) {
      if (!controller.signal.aborted) setMessage(errorMessage(error));
    } finally {if (!controller.signal.aborted) setSearching(false);}
  }
  return <>
    <Sidebar motion={motion} />
    <main>
      <header className="topbar"><span>Workspace <span className="slash">/</span> <strong>Variant intelligence</strong></span><div className="topbar-status"><span className={`connection ${statsError ? 'offline' : stats ? 'online' : ''}`}>{statsError ? 'Data service unavailable' : stats ? 'Data service connected' : 'Connecting…'}</span><button className="motion-toggle" onClick={() => setMotion(value => !value)}>{motion ? 'Pause motion' : 'Play motion'}</button></div></header>
      <section id="overview" className="intro"><div><p className="eyebrow">GENEVISTA KNOWLEDGE BASE</p><h1>Variant evidence & classification history</h1><p className="intro-copy">Explore ClinVar records and their observed classifications across six snapshots.</p><details className="dataset-facts"><summary>View dataset facts</summary><p>GRCh38 records · January 2022 to September 2026. Search by exact gene symbol or numeric VariationID. Counts describe prepared records and historical coverage.</p></details></div><span className="release">{stats ? `LATEST SNAPSHOT · ${stats.latest_snapshot}` : 'CLINVAR SNAPSHOTS'}</span></section>
      {statsError && <div className="notice error" role="alert">{statsError} <button className="secondary" onClick={() => setStatsAttempt(value => value + 1)}>Retry connection</button></div>}
      {stats?.dataset_mode === 'demo' && <section className="demo-notice notice" aria-label="Demo dataset"><strong>Demo dataset · {stats.variant_count.toLocaleString()} variants</strong><p>This workspace contains a selected sample. Searches outside the sample may return no results. Counts describe the demo, while model evaluation uses the original historical cohorts.</p><details><summary>Included genes and sampling</summary><p>{stats.demo_genes.join(', ')}</p><p>{stats.scope}</p></details></section>}
      <Stats data={stats} />
      <section className="search-panel" id="explorer"><div><p className="eyebrow">EXPLORE THE DATA</p><h2>Start with a gene or variant.</h2></div>
        <form onSubmit={event => { event.preventDefault(); runSearch(query); }}>
          <label className="sr-only" htmlFor="query">Exact gene symbol or numeric variant ID</label>
          <div className="search-input"><span aria-hidden="true">⌕</span><input id="query" ref={inputRef} value={query} onChange={event => setQuery(event.target.value)} placeholder="Enter a gene symbol or VariationID" maxLength={100} required autoComplete="off" /><button disabled={searching} type="submit">{searching ? 'Searching…' : 'Search'} <span aria-hidden="true">→</span></button></div>
        </form>
        <div className="examples"><span>TRY AN EXAMPLE</span>{['BRCA1', 'TP53', '17660', '37565'].map(value => <button key={value} onClick={() => runSearch(value)}>{value === '17660' ? 'Variant 17660' : value === '37565' ? 'VUS 37565' : value}</button>)}<small>Exact gene symbols · Numeric variant IDs</small></div>
      </section>
      <nav className="section-tabs" aria-label="Evidence sections"><a href="#variant-results">Variant summaries</a><a href="#sequence">Sequence & substitutions</a><a href="#about">Data provenance</a><a href="https://www.ncbi.nlm.nih.gov/clinvar/" target="_blank" rel="noreferrer">ClinVar ↗</a></nav>
      <div id="variant-results" className="search-status" role="status" aria-live="polite">{message}</div>
      {!result && !searching && <section className="empty-results"><h2>Variant summaries</h2><p>Search for a gene or variant above to view classifications, review status, and snapshot history.</p><div className="empty-columns" aria-hidden="true"><span>Variant</span><span>Gene</span><span>Classification</span><span>History</span></div></section>}
      <Results result={result} loading={searching} onMore={() => runSearch(result.query, true)} onSelect={openDetail} />
      <Detail selected={selected} onClose={closeDetail} onRetry={() => openDetail(selected.id)} />
      <EvidencePanel selected={selected} />
      <RankingPanel selected={selected} />
      <ModelCharts />
      <MolecularJourney motion={motion} theme={theme} variant={selected?.data?.latest_details} />
      <div id="sequence"><DnaAnimation variant={selected?.data?.latest_details} motion={motion} onToggleMotion={() => setMotion(value => !value)} /></div>
      <About />
      <section id="settings" className="settings-module"><div><p className="eyebrow">SETTINGS</p><h2>Make this workspace yours</h2></div><label>Theme<select value={theme} onChange={event=>setTheme(event.target.value)}><option value="dark">Dark Mode</option><option value="light">Light Mode</option></select></label><label className="motion-setting"><input type="checkbox" checked={motion} onChange={event=>setMotion(event.target.checked)}/> Enable animations</label><p>Your theme is saved in this browser. The molecular lesson can also be explored one chapter at a time.</p></section>
      <footer>GeneVISTA <span>Variant evidence explorer · Research use</span><span>Data: prepared ClinVar snapshots</span></footer>
    </main>
  </>;
}
