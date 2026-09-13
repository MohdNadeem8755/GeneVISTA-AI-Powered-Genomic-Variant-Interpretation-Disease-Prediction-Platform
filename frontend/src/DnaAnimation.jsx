import React, { useEffect, useMemo, useState } from 'react';
import './dna.css';

const complement = { A: 'T', T: 'A', C: 'G', G: 'C' };
const baseNames = { A: 'Adenine', C: 'Cytosine', G: 'Guanine', T: 'Thymine' };
const context = ['A', 'C', 'T', 'G', 'A', 'T', 'C', 'G', 'T'];

function substitutionFromName(name = '') {
  const matches = [...name.matchAll(/(?:c|g|m|n)\.(?:-?\d+|\*\d+)(?:[+-]\d+)?([ACGT])>([ACGT])/g)];
  const changes = [...new Set(matches.map(match => `${match[1]}>${match[2]}`))];
  return changes.length === 1 ? changes[0].split('>') : null;
}

export function DnaLogo({ motion }) {
  return <span className={`brand-icon dna-brand ${motion ? '' : 'motion-paused'}`} aria-hidden="true">
    <svg viewBox="0 0 36 44" fill="none"><path d="M9 3C9 13 27 13 27 23S9 33 9 43" stroke="#70ded0" strokeWidth="2"/><path d="M27 3C27 13 9 13 9 23S27 33 27 43" stroke="#b6a0e4" strokeWidth="2"/>{[0, 1, 2, 3, 4, 5].map(index =>
      <g key={index} style={{ animationDelay: `${index * -0.22}s` }}>
        <path d={`M9 ${7 + index * 6}H27`} stroke="#6bbdad" strokeWidth="1.5" />
        <circle cx="9" cy={7 + index * 6} r="2.7" fill={index % 2 ? '#e6ba63' : '#70ded0'} />
        <circle cx="27" cy={7 + index * 6} r="2.7" fill={index % 2 ? '#b6a0e4' : '#d28e9c'} />
      </g>
    )}</svg>
  </span>;
}

export default function DnaAnimation({ variant, motion, onToggleMotion }) {
  const observed = useMemo(() => substitutionFromName(variant?.name), [variant?.name]);
  const [reference, setReference] = useState('A');
  const [alternate, setAlternate] = useState('G');
  const [changed, setChanged] = useState(false);
  useEffect(() => {
    setChanged(false);
    if (observed) { setReference(observed[0]); setAlternate(observed[1]); }
  }, [observed]);
  useEffect(() => {
    if (!motion) return;
    const timer = setInterval(() => setChanged(value => !value), 3500);
    return () => clearInterval(timer);
  }, [motion, reference, alternate]);
  const active = changed ? alternate : reference;
  const sequence = context.map((base, index) => index === 4 ? active : base);
  const isExample = !observed;
  function chooseReference(value) {
    setReference(value);
    if (value === alternate) setAlternate(value === 'G' ? 'A' : 'G');
    setChanged(false);
  }
  return <section className={`dna-panel ${motion ? '' : 'motion-paused'}`} aria-label="DNA substitution illustration">
    <div className="dna-copy"><p className="eyebrow">THE LANGUAGE OF DNA</p><h2>One letter can change<br /> the sequence.</h2><p>DNA pairs <strong>A with T</strong> and <strong>C with G</strong>. A substitution replaces one base; this diagram shows the corresponding pair in the resulting sequence.</p>
      <div className="base-legend">{Object.entries(baseNames).map(([base, name]) => <span key={base}><b className={`base-${base}`}>{base}</b>{name}</span>)}</div>
      <div className="dna-controls"><button className="secondary" onClick={onToggleMotion}>{motion ? 'Pause animations' : 'Play animations'}</button><button className="secondary" onClick={() => setChanged(value => !value)}>{changed ? 'Show reference' : 'Show substitution'}</button></div>
    </div>
    <div className="dna-visual">
      <div className="dna-visual-top"><span className="demo-tag">{isExample ? 'ILLUSTRATIVE EXAMPLE' : 'SUBSTITUTION FROM VARIANT NAME'}</span><span className={`change-state ${changed ? 'is-changed' : ''}`}>{changed ? 'Alternate sequence' : 'Reference sequence'}</span></div>
      <div className="substitution-title"><span className={`base-${reference}`}>{reference}</span><span className="change-arrow">→</span><span className={`base-${alternate}`}>{alternate}</span><small>Single-base substitution</small></div>
      {isExample && <div className="base-selectors"><label>Reference <select value={reference} onChange={event => chooseReference(event.target.value)}>{Object.keys(baseNames).map(base => <option key={base}>{base}</option>)}</select></label><label>Alternate <select value={alternate} onChange={event => {setAlternate(event.target.value); setChanged(false);}}>{Object.keys(baseNames).filter(base => base !== reference).map(base => <option key={base}>{base}</option>)}</select></label></div>}
      <div className="dna-strand" role="img" aria-label={`Schematic DNA strand. Central base ${reference} changes to ${alternate}. Currently showing ${active}, paired with ${complement[active]}. Surrounding bases are illustrative.`}>
        <div className="strand-label"><span>5′</span><span>3′</span></div>
        <div className="base-pairs">{sequence.map((base, index) => <div key={index} className={`base-pair ${index === 4 ? 'mutation-site' : ''} ${index === 4 && changed ? 'mutated' : ''}`}>
          <span className={`base-tile base-${base}`} key={`top-${base}`}>{base}</span>
          <span className={`pair-bonds ${base === 'C' || base === 'G' ? 'three-bonds' : ''}`} aria-hidden="true" />
          <span className={`base-tile base-${complement[base]}`} key={`bottom-${base}`}>{complement[base]}</span>
          {index === 4 && <span className="site-marker">variant site</span>}
        </div>)}</div>
        <div className="strand-label"><span>3′</span><span>5′</span></div>
      </div>
      <p className="dna-caption">{isExample ? (variant ? 'This variant name does not provide one unambiguous single-base substitution. Showing an example.' : 'Explore an example, or select a variant with a readable substitution.') : `Showing ${reference} → ${alternate} from the selected variant’s name.`} Surrounding bases are schematic, not a retrieved reference sequence. This animation does not indicate disease risk.</p>
    </div>
  </section>;
}
