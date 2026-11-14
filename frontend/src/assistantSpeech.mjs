// Pronunciation changes affect audio only; scientific text remains unchanged.
export function spokenText(text, language) {
  const hindi = language.startsWith('hi');
  let result = text.replace(/\bGene\s*Vista\b/gi, hindi ? 'जीन विस्टा' : 'Gene Vista');
  if (hindi) {
    const terms = {DNA:'डी एन ए', RNA:'आर एन ए', VUS:'वी यू एस', CatBoost:'कैट बूस्ट',
      'Likely benign':'लाइकली बिनाइन', 'Likely pathogenic':'लाइकली पैथोजेनिक',
      Benign:'बिनाइन', Pathogenic:'पैथोजेनिक', precision:'प्रिसिजन', recall:'रिकॉल'};
    for (const [term, pronunciation] of Object.entries(terms)) {
      result = result.replace(new RegExp('\\b'+term+'\\b','gi'), pronunciation);
    }
    const bases = {A:'ए', T:'टी', C:'सी', G:'जी'};
    result = result.replace(/\b[ATCG]\b/g, base=>bases[base]);
  } else {
    // Only explicit substitutions: never change the English article “A”.
    result = result.replace(/\b([ATCG])\s*(?:→|>|to)\s*([ATCG])\b/g, (_,a,b)=>`base ${a}, to base ${b}`);
  }
  return result.replace(/→|>/g,hindi?' से ':' to ').replace(/\*\*|`/g,'');
}

export function matchingVoices(voices, language) {
  const normalized = value=>value.replace('_','-').toLowerCase();
  return voices.filter(v=>v.localService && normalized(v.lang).split('-')[0]===language.split('-')[0])
    .sort((a,b)=>Number(normalized(b.lang)===normalized(language))-Number(normalized(a.lang)===normalized(language)));
}

export function speechParts(text, language) {
  const spoken = spokenText(text,language);
  if (typeof Intl.Segmenter === 'function') {
    return [...new Intl.Segmenter(language,{granularity:'sentence'}).segment(spoken)]
      .map(item=>item.segment.trim()).filter(Boolean);
  }
  return spoken.split(/\n+/).map(part=>part.trim()).filter(Boolean);
}
