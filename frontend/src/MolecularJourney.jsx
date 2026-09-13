import React, {useEffect, useRef, useState} from 'react';

const chapters = [
  ['DNA', 'A double-stranded information store', 'Complementary bases join the two strands: A pairs with T, and C pairs with G. The helix rotates here to reveal its structure.'],
  ['Gene region', 'A gene is a region of DNA', 'The highlighted region represents a protein-coding gene. Genes are part of DNA; they are not a separate material that DNA turns into.'],
  ['Variant', 'Compare a reference and alternate sequence', 'Variants are sequence differences. Select an example to explore a changed base or a deletion. A difference alone does not establish a harmful effect.'],
  ['RNA', 'Transcription carries the sequence forward', 'A protein-coding gene is transcribed into RNA. After processing, messenger RNA can carry coding information to a ribosome. RNA uses U in place of T.'],
  ['Protein', 'The ribosome reads three bases at a time', 'During translation, codons specify amino acids or a stop signal. A sequence change can preserve an amino acid, replace one, introduce a stop, or shift the reading frame.'],
  ['Function', 'Sequence is one part of the evidence', 'Some changes affect protein abundance, folding, binding, or activity; others have little effect. Functional impact requires evidence. This illustrative fold is not a computed structure or a disease prediction.']
];
const examples = {
  missense: {label:'Missense substitution',dna:'ATGGTATTTGGCTAA',rna:'AUG GUA UUU GGC UAA',protein:'Met · Val · Phe · Gly · Stop',effect:'GAA → GTA changes glutamate to valine in this example.'},
  synonymous: {label:'Synonymous substitution',dna:'ATGGAGTTTGGCTAA',rna:'AUG GAG UUU GGC UAA',protein:'Met · Glu · Phe · Gly · Stop',effect:'GAA → GAG still encodes glutamate. The amino-acid sequence is unchanged in this example.'},
  nonsense: {label:'Stop-gain substitution',dna:'ATGTAATTTGGCTAA',rna:'AUG UAA UUU GGC UAA',protein:'Met · Stop',effect:'GAA → TAA introduces an early stop. RNA surveillance may also change transcript abundance.'},
  frameshift: {label:'Single-base deletion',dna:'ATGAATTTGGCTAA',rna:'AUG AAU UUG GCU AA…',protein:'Met · Asn · Leu · Ala · …',effect:'Deleting one G shifts the downstream codon grouping. The final stop position is outside this short example.'}
};

function MolecularCanvas({stage,clock,kind,theme}) {
  const ref = useRef(null);
  useEffect(() => {
    const canvas=ref.current, ctx=canvas.getContext('2d');
    const width=canvas.clientWidth, height=canvas.clientHeight, ratio=Math.min(window.devicePixelRatio||1,2);
    canvas.width=width*ratio;canvas.height=height*ratio;ctx.scale(ratio,ratio);
    const sx=width/1000,sy=height/430;ctx.scale(sx,sy);
    const dark=theme==='dark'; const ink=dark?'#bdcfdd':'#34536c';
    const colors=['#3de1bd','#65a5ff','#f6c86c','#e88ccb'];
    function orb(x,y,r,color){const g=ctx.createRadialGradient(x-r*.3,y-r*.4,1,x,y,r);g.addColorStop(0,'#ffffff');g.addColorStop(.24,color);g.addColorStop(1,dark?'#112338':color);ctx.fillStyle=g;ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ctx.fill();}
    function line(x,y,x2,y2,color,w=2){ctx.strokeStyle=color;ctx.lineWidth=w;ctx.beginPath();ctx.moveTo(x,y);ctx.lineTo(x2,y2);ctx.stroke();}
    ctx.clearRect(0,0,1000,430);
    for(let i=0;i<45;i++){orb((i*127+clock*3)%1000,25+(i*73)%385,1.4,dark?'#284253':'#b9cfda');}
    ctx.font='14px Arial';ctx.fillStyle=ink;
    if(stage<3){
      ctx.save();ctx.translate(500,220);ctx.rotate(-.13);
      const nodes=[];
      for(let i=0;i<32;i++){let x=-410+i*26,phase=i*.42+clock*.85,y=Math.sin(phase)*78,z=Math.cos(phase);nodes.push({x,y,z,i});}
      nodes.forEach(({x,y,z,i},index)=>{let glow=stage===1&&i>10&&i<22;let color=glow?'#f6c86c':dark?'#335970':'#9bbaca';line(x,y,x,-y,color,glow?4:1.5);if(index){let prev=nodes[index-1];line(prev.x,prev.y,x,y,'#349db2',4);line(prev.x,-prev.y,x,-y,'#675fa9',4);}let mutation=stage===2&&i===16;orb(x,y,8+z*2,mutation?'#ff746e':colors[i%4]);orb(x,-y,8-z*2,colors[(i+2)%4]);if(mutation){ctx.strokeStyle='#ffb58b';ctx.lineWidth=2;ctx.beginPath();ctx.arc(x,y,20+Math.sin(clock*2)*3,0,Math.PI*2);ctx.stroke();}}
      );ctx.restore();ctx.fillStyle=ink;ctx.fillText(stage===1?'Highlighted region: illustrative protein-coding gene':stage===2?'Variant site highlighted · sequence comparison below':'Sugar–phosphate backbones and complementary base pairs',35,390);
    } else if(stage===3){
      ctx.strokeStyle=dark?'#314354':'#b0c4d2';ctx.lineWidth=2;ctx.beginPath();ctx.ellipse(210,210,160,165,0,0,Math.PI*2);ctx.stroke();ctx.fillText('NUCLEUS',160,35);
      for(let i=0;i<17;i++){let y=75+i*16;line(130+Math.sin(i*.5)*30,y,260-Math.sin(i*.5)*30,y,'#3c687a');orb(130+Math.sin(i*.5)*30,y,5,colors[i%4]);orb(260-Math.sin(i*.5)*30,y,5,colors[(i+2)%4]);}
      let progress=(clock*.15)%1;orb(300+progress*120,205,28,'#76b7d4');ctx.fillStyle=ink;ctx.fillText('RNA polymerase (schematic)',290,155);
      for(let i=0;i<24;i++){let x=350+i*22,y=235+Math.sin(i*.4+clock)*24;orb(x,y,7,'#f4b969');if(i<23)line(x,y,x+22,235+Math.sin((i+1)*.4+clock)*24,'#c99353',2);}
      ctx.fillStyle=ink;ctx.fillText('RNA transcript → processing → mRNA',430,330);
    } else if(stage===4){
      const offset=(clock*25)%480;
      for(let i=0;i<30;i++){let x=60+i*30;line(x,285,x+30,285,'#b99866',3);orb(x,285,8,colors[i%4]);}
      ctx.globalAlpha=.8;orb(220+offset,235,70,'#5ba0bb');orb(220+offset,300,42,'#897dbf');ctx.globalAlpha=1;
      const count=kind==='nonsense'?2:kind==='frameshift'?10:8;
      for(let i=0;i<count;i++){let x=220+offset+i*15,y=185-i*13+Math.sin(clock+i)*8; if(i)line(x-15,y+13,x,y,'#608c9e',3);orb(x,y,9,i===2&&kind==='missense'?'#ff8172':'#42cdb5');}
      ctx.fillStyle=ink;ctx.fillText('Ribosome reads mRNA codons',65,375);ctx.fillText('Emerging amino-acid chain',460,75);
    } else {
      const points=[];let count=kind==='nonsense'?20:85;
      for(let i=0;i<count;i++){let a=i*.31+clock*.3;let x=500+Math.sin(a)*145+Math.cos(i*.14)*80,y=210+Math.cos(a)*75+Math.sin(i*.2)*55,z=Math.sin(i*.45+clock*.3);points.push({x,y,z,i});}
      points.forEach((p,i)=>{if(i)line(points[i-1].x,points[i-1].y,p.x,p.y,'#388899',4);});points.sort((a,b)=>a.z-b.z).forEach(p=>orb(p.x,p.y,6+p.z*2,p.i===25&&kind==='missense'?'#ff8479':'#53c7c1'));
      ctx.fillStyle=ink;ctx.fillText('Illustrative molecular conformation · not a structure prediction',35,390);
    }
  },[stage,clock,kind,theme]);
  return <canvas ref={ref} className="molecular-canvas" role="img" aria-label={`Animated schematic: ${chapters[stage][0]}. ${chapters[stage][2]}`} />;
}

export default function MolecularJourney({motion,theme,variant}) {
  const [stage,setStage]=useState(0),[playing,setPlaying]=useState(false),[clock,setClock]=useState(0),[speed,setSpeed]=useState(1),[kind,setKind]=useState('missense');
  const [alternate,setAlternate]=useState(true);
  useEffect(()=>{if(!motion)setPlaying(false);},[motion]);
  useEffect(()=>{if(!playing||!motion)return;let frame,last=0,elapsed=0;const tick=now=>{if(!last)last=now;let dt=Math.min((now-last)/1000,.1)*speed;last=now;elapsed+=dt;setClock(value=>value+dt);if(elapsed>=9){elapsed=0;setStage(value=>(value+1)%chapters.length);}frame=requestAnimationFrame(tick);};frame=requestAnimationFrame(tick);return()=>cancelAnimationFrame(frame);},[playing,motion,speed]);
  const sample=examples[kind], dna=alternate?sample.dna:'ATGGAATTTGGCTAA';
  return <section className="molecular-module" id="molecular">
    <div className="science-heading"><div><p className="eyebrow">MOLECULAR STUDIO · INTERACTIVE LESSON</p><h2>From DNA sequence to protein function</h2></div><span className="scope-pill">Educational schematic</span></div>
    <div className="journey-stages" aria-label="Molecular chapters">{chapters.map(([name],i)=><button key={name} aria-pressed={stage===i} onClick={()=>{setStage(i);setPlaying(false);}}><small>0{i+1}</small>{name}</button>)}</div>
    <div className="molecular-stage"><div className="scene-caption"><span>CHAPTER 0{stage+1} / 06</span><h3>{chapters[stage][1]}</h3></div><MolecularCanvas stage={stage} clock={clock} kind={kind} theme={theme}/><div className="scene-controls"><button className="studio-primary" disabled={!motion} onClick={()=>setPlaying(v=>!v)}>{playing?'Pause lesson':'Play lesson'}</button><button onClick={()=>{setStage((stage+1)%6);setPlaying(false);}}>Next chapter →</button><label>Speed <select value={speed} onChange={e=>setSpeed(Number(e.target.value))}><option value="0.5">0.5×</option><option value="1">1×</option><option value="1.5">1.5×</option></select></label><span>{!motion?'Motion disabled in Settings':'Six chapters · ~54 seconds at 1×'}</span></div></div>
    <div className="lesson-notes"><p>{chapters[stage][2]}</p><div><label htmlFor="variant-example">Explore a coding-sequence example</label><select id="variant-example" value={kind} onChange={e=>{setKind(e.target.value);setStage(2);setPlaying(false);}}>{Object.entries(examples).map(([key,value])=><option key={key} value={key}>{value.label}</option>)}</select></div></div>
    <div className="sequence-lab"><div className="science-heading"><h3>Sequence → codons → amino acids</h3><button onClick={()=>setAlternate(v=>!v)}>{alternate?'Show reference':'Show alternate'}</button></div><div className="sequence-triptych"><div><small>{alternate?'ALTERNATE':'REFERENCE'} CODING DNA · 5′ → 3′</small><p className="nucleotide-line">{dna.split('').map((base,i)=><b key={i} className={`base-${base} ${alternate&&base!=='ATGGAATTTGGCTAA'[i]?'changed-base':''}`}>{base}</b>)}</p></div><div><small>mRNA · 5′ → 3′</small><p>{alternate?sample.rna:'AUG GAA UUU GGC UAA'}</p></div><div><small>AMINO-ACID SEQUENCE</small><p>{alternate?sample.protein:'Met · Glu · Phe · Gly · Stop'}</p></div></div><p>{alternate?sample.effect:'Reference example: a short, synthetic coding sequence. It does not represent a complete gene.'}</p></div>
    {variant&&<p className="selected-protein"><strong>Selected variant annotation:</strong> {variant.name}<br/>The lesson above uses synthetic sequences. No experimentally determined protein structure or variant-specific functional simulation is connected.</p>}
    <details className="evidence-method"><summary>Scientific context & learning resources</summary><p>This lesson follows a protein-coding example. Noncoding genes and variants can affect regulation or RNA processing. A gene is a DNA region, and a variant is a sequence difference, rather than a separate step in gene expression.</p><a href="https://www.genome.gov/genetics-glossary/Transcription" target="_blank" rel="noreferrer">NHGRI · Transcription ↗</a> · <a href="https://www.genome.gov/genetics-glossary/Translation" target="_blank" rel="noreferrer">NHGRI · Translation ↗</a></details>
  </section>;
}
