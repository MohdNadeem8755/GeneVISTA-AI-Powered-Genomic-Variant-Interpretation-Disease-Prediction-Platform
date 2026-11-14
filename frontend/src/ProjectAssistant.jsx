import React, {useEffect, useRef, useState} from 'react';
import {askAssistant, fetchApi, errorMessage} from './api';
import './ProjectAssistant.css';
import {matchingVoices, speechParts} from './assistantSpeech.mjs';

export default function ProjectAssistant() {
  const [open,setOpen]=useState(false);
  const [voices,setVoices]=useState([]);
  const [voiceId,setVoiceId]=useState(''), [rate,setRate]=useState(0.94);
  const [speaking,setSpeaking]=useState(false);
  const playback=useRef({token:0,timer:null,utterance:null});
  const input=useRef(null), bottom=useRef(null), launcher=useRef(null);
  const [question,setQuestion]=useState(''), [mode,setMode]=useState('guide');
  const [available,setAvailable]=useState({}), [messages,setMessages]=useState([]);
  const [busy,setBusy]=useState(false), [error,setError]=useState('');
  const [listening,setListening]=useState(false), [language,setLanguage]=useState('en-IN');
  const recognition=useRef(null);
  useEffect(()=>{
    const synth=window.speechSynthesis;
    const load=()=>setVoices(synth?.getVoices() || []);
    load();synth?.addEventListener('voiceschanged',load);
    return ()=>synth?.removeEventListener('voiceschanged',load);
  },[]);
  useEffect(()=>{if(open)input.current?.focus();},[open]);
  useEffect(()=>{if(open)bottom.current?.scrollIntoView({block:'nearest'});},[messages,busy,open]);
  function stopAudio(){playback.current.token++;clearTimeout(playback.current.timer);window.speechSynthesis?.cancel();playback.current.utterance=null;setSpeaking(false);}
  function close(){setOpen(false);recognition.current?.abort();stopAudio();launcher.current?.focus();}
  useEffect(()=>()=>{playback.current.token++;clearTimeout(playback.current.timer);window.speechSynthesis?.cancel();},[]);
  useEffect(()=>{fetchApi('/api/assistant/status').then(setAvailable).catch(()=>{});return ()=>{recognition.current?.abort();window.speechSynthesis?.cancel();};},[]);
  async function send(event){
    event.preventDefault();if(!question.trim() || busy)return;
    setBusy(true);setError('');const asked=question;setQuestion('');
    try{const reply=await askAssistant(asked,mode,language);setMessages(old=>[...old.slice(-9),{asked,...reply,language}]);}
    catch(e){setError(errorMessage(e));setQuestion(asked);}finally{setBusy(false);}
  }
  function listen(){
    if(listening){recognition.current?.stop();return;}
    const Speech=window.SpeechRecognition || window.webkitSpeechRecognition;
    if(!Speech){setError('Voice input unsupported. Please type.');return;}
    const mic=new Speech();recognition.current=mic;mic.lang=language;mic.interimResults=false;
    mic.onresult=e=>setQuestion(e.results[0][0].transcript);
    mic.onerror=e=>{setError('Voice input: '+e.error+'. Please type instead.');setListening(false);};
    mic.onend=()=>setListening(false);
    try{mic.start();setListening(true);setError('');}catch{setError('Cannot start microphone.');}
  }
  function speak(text, answerLanguage=language){
    if(!window.speechSynthesis){setError('Read aloud unsupported.');return;}
    const choices=matchingVoices(window.speechSynthesis.getVoices(),answerLanguage);
    const voice=choices.find(v=>v.voiceURI===voiceId) || choices[0];
    if(!voice){setError('No local voice available yet. Enable an English or Hindi voice in your system speech settings; text chat still works.');return;}
    setError('');
    stopAudio();setSpeaking(true);
    const token=playback.current.token, parts=speechParts(text,answerLanguage);
    function next(index){
      if(token!==playback.current.token)return;
      if(index>=parts.length){setSpeaking(false);playback.current.utterance=null;return;}
      const speech=new SpeechSynthesisUtterance(parts[index]);
      playback.current.utterance=speech;speech.voice=voice;speech.lang=voice.lang;speech.rate=rate;speech.pitch=1;
      speech.onend=()=>{if(token===playback.current.token)playback.current.timer=setTimeout(()=>next(index+1),180);};
      speech.onerror=()=>{if(token===playback.current.token){setSpeaking(false);setError('Audio could not play. Try another installed voice.');}};
      window.speechSynthesis.speak(speech);
    }
    next(0);
  }
  return <div className="gv-assistant">
    <button ref={launcher} className={'gv-launcher'+(open?' active':'')} aria-label={open?'Close assistant':'Open GeneVISTA assistant'} aria-expanded={open} aria-controls="project-assistant" onClick={()=>open?close():setOpen(true)}><svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.8" aria-hidden="true"><rect x="9" y="2" width="6" height="12" rx="3"/><path d="M5 10v2a7 7 0 0014 0v-2M12 19v3M8 22h8"/></svg><span className="gv-launch-label">Ask GeneVISTA</span></button>
    {open && <section id="project-assistant" className="gv-chat" role="region" aria-label="GeneVISTA assistant" onKeyDown={e=>{if(e.key==='Escape')close();}}>
    <header className="gv-chat-head"><div><span className="gv-status-dot"/> GeneVISTA <small>Your project companion</small></div><button onClick={close} aria-label="Close chat">×</button></header>
    <div className="gv-chat-scroll">
    <p className="eyebrow">GENEVISTA ASSISTANT</p><h2>Ask about the project</h2>
    <p>Ask about data, features, training and results. Each question is independent. Offline guide covers common topics; configured AI supports broader questions.</p>
    <label>Answer mode <select value={mode} onChange={e=>setMode(e.target.value)}><option value="guide">Offline topic guide</option><option value="local" disabled={!available.local}>Local AI {available.local?'':'(not configured)'}</option><option value="online" disabled={!available.online}>Online AI {available.online?'':'(not configured)'}</option></select></label>
    <label>Answer & voice language <select value={language} onChange={e=>{stopAudio();setVoiceId('');setLanguage(e.target.value);}}><option value="en-IN">English (India)</option><option value="hi-IN">Hindi (India)</option></select></label>
    <label>Speaking voice <select value={voiceId} onChange={e=>{stopAudio();setVoiceId(e.target.value);}}><option value="">Automatic · prefer Indian voice</option>{matchingVoices(voices,language).map(v=><option key={v.voiceURI} value={v.voiceURI}>{v.name} ({v.lang})</option>)}</select></label>
    <label>Speaking pace <select value={rate} onChange={e=>{stopAudio();setRate(Number(e.target.value));}}><option value={0.85}>Relaxed</option><option value={0.94}>Conversational</option><option value={1}>Normal</option></select></label>
    <button onClick={()=>speak(language==='hi-IN'?'नमस्ते! चलिए, GeneVISTA को समझते हैं। डीएनए में A की जगह T आ जाए, तो उसे एक बदलाव कहते हैं।':'Hello! Let’s explore GeneVISTA together. A change from A to T is a DNA substitution.')}>Test pronunciation</button>
    {speaking && <p role="status">Speaking…</p>}
    {mode==='online' && <p>Send shares your question and project report context with Google Gemini. Do not enter patient identifiers.</p>}
    <p>Microphone transcription may use your browser’s online service, even with offline answers. Read aloud uses installed local voices.</p>
    <div style={{display:'flex',gap:8,flexWrap:'wrap'}}>{['Explain GeneVISTA','Latest scores and precision','Explain dataset and split','Improve weak classes'].map(q=><button key={q} onClick={()=>setQuestion(q)}>{q}</button>)}</div>
    {error && <p role="alert">{error}</p>}
    <div aria-live="polite">{messages.map((message,i)=><article key={i} style={{borderTop:'1px solid #365061',marginTop:16,paddingTop:12}}><strong>You: {message.asked}</strong><p style={{whiteSpace:'pre-wrap'}}>{message.answer}</p><small>{message.mode} · Sources: {message.sources.join('; ') || 'No matching source'}</small><p><button onClick={()=>speak(message.answer,message.language)}>Read aloud</button></p></article>)}</div>
    {busy && <p role="status" className="gv-thinking">Thinking <span>•••</span></p>}<div ref={bottom}/></div>
    <form onSubmit={send} className="gv-composer"><label className="sr-only" htmlFor="assistant-question">Your question</label><textarea ref={input} id="assistant-question" placeholder="Ask about GeneVISTA…" value={question} onChange={e=>setQuestion(e.target.value)} maxLength={2000} rows={2} required/><div><button type="button" className={listening?'gv-listening':''} onClick={listen}>{listening?'● Stop mic':'🎙 Mic'}</button><button type="button" onClick={stopAudio}>Stop audio</button><button className="gv-send" disabled={busy || !question.trim()}>{busy?'…':'Send ↑'}</button></div></form>
  </section>}</div>;
}
