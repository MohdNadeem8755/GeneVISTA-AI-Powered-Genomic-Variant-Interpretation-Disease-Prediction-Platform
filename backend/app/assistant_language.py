PROJECT_EN = '''Let’s walk through GeneVISTA with an example. Suppose, at a particular position in DNA, A changes to T. That change is called a variant. The question is: what evidence do we have about its effect? The letters alone are not enough; the position, gene and biological context matter too.

GeneVISTA brings variant records and supporting evidence together. You can search a gene or a variant ID, inspect its recorded classification and history, and see a model prediction. It also ranks disease annotation matches; those rankings are not a person’s probability of developing a disease.

The model learns five labels: Benign, Likely benign, VUS, Likely pathogenic and Pathogenic. “Likely” describes lower certainty. VUS means the evidence is not sufficient to decide, so it must not be read as harmless.

Before training, the project filters records for the intended genome build and variant type, consolidates duplicate locations and excludes conflicting labels. That is why all the millions of database records are not automatically used for training. A snapshot simply means a copy of the data saved on a particular date.

Features describe the variant: its gene, DNA and protein changes, population frequency, conservation and available biological evidence. A functional assay measures an effect in an experiment; splice evidence concerns how RNA is joined. Missing measurements stay missing, rather than being treated as normal results.

CatBoost models and ensembles are compared on development data. Training fits the model, early stopping chooses when to stop, and validation helps choose a candidate. A genuinely untouched test is needed for a final independent assessment. The historical internal test has already been inspected.

On the dashboard, look at each class’s precision and recall, not only overall accuracy. Precision tells you how often a class prediction is right; recall tells you how many actual cases it finds. The corrected candidate improved precision for the two Likely classes but missed more of their actual cases. It is not deployed: individual predictions still use the active annotation model.

So the project’s purpose is to help explore and prioritise variant evidence for research. Its predictions support investigation; they do not establish a diagnosis.'''

PROJECT_HI = '''चलिए, एक आसान उदाहरण से समझते हैं। मान लीजिए डीएनए की किसी खास जगह पर A की जगह T आ गया। इस बदलाव को वेरिएंट कहते हैं। अब सवाल है कि इस बदलाव का असर क्या हो सकता है? सिर्फ A और T देखकर जवाब नहीं मिलता। बदलाव की जगह, जीन और उसके बारे में उपलब्ध सबूत भी चाहिए।

GeneVISTA इन्हीं वेरिएंट रिकॉर्ड और सबूतों को एक जगह दिखाता है। आप जीन या वेरिएंट आईडी खोज सकते हैं, उसकी दर्ज की गई श्रेणी और पुराना इतिहास देख सकते हैं, और मॉडल का अनुमान समझ सकते हैं। बीमारी से जुड़े जो मैच दिखते हैं, वे किसी इंसान को बीमारी होने की संभावना नहीं हैं।

मॉडल पाँच श्रेणियाँ सीखता है। Benign यानी उपलब्ध सबूत उस संदर्भ में बीमारी पैदा न करने का समर्थन करते हैं। Pathogenic यानी बीमारी पैदा करने का समर्थन करते हैं। Likely benign और Likely pathogenic में भरोसा थोड़ा कम होता है। VUS का मतलब है कि अभी फैसला करने लायक सबूत नहीं हैं। इसका मतलब सुरक्षित होना नहीं है।

ट्रेनिंग से पहले सही जीनोम संस्करण और चुने हुए बदलावों के रिकॉर्ड लिए जाते हैं। एक ही जगह के दोहराए रिकॉर्ड जोड़े जाते हैं और आपस में टकराते लेबल हटाए जाते हैं। इसलिए डेटाबेस के सारे लाखों रिकॉर्ड सीधे ट्रेनिंग में नहीं जाते। स्नैपशॉट बस किसी तारीख को सेव की गई डेटा की कॉपी है।

मॉडल को जीन, डीएनए और प्रोटीन का बदलाव, आबादी में उसकी आवृत्ति और उपलब्ध जैविक सबूत दिए जाते हैं। फंक्शनल ऐसे में प्रयोग करके असर मापा जाता है। स्प्लाइस का संबंध आरएनए के हिस्सों के जुड़ने से है। कोई जानकारी गायब हो तो उसे सामान्य नतीजा या शून्य नहीं मानते।

फिर CatBoost मॉडल और उनके समूह, यानी एन्सेम्बल, की तुलना होती है। ट्रेनिंग में मॉडल सीखता है। अर्ली स्टॉपिंग बताती है कि कब रोकना है। वैलिडेशन से मॉडल चुना जाता है। आखिरी निष्पक्ष जाँच के लिए ऐसा टेस्ट चाहिए जिसे चुनने या सुधारने के दौरान देखा न गया हो। हमारा पुराना आंतरिक टेस्ट पहले देखा जा चुका है।

परफॉर्मेंस समझने के लिए हर श्रेणी का प्रिसिजन और रिकॉल साथ देखिए। प्रिसिजन बताता है कि उस श्रेणी के दिए हुए अनुमान कितनी बार सही हैं। रिकॉल बताता है कि उसके असली मामलों में से कितने मिले। नए उम्मीदवार में दोनों Likely श्रेणियों का प्रिसिजन बढ़ा, लेकिन उनके कुछ और असली मामले छूटे। नया उम्मीदवार अभी लगाया नहीं गया है; व्यक्तिगत अनुमान पुराने सक्रिय मॉडल से आते हैं।

तो आसान शब्दों में, GeneVISTA बदलावों के सबूत समझने और शोध में प्राथमिकता तय करने में मदद करता है। यह किसी मरीज की बीमारी का पक्का फैसला नहीं करता।'''

HINDI = [PROJECT_HI,
 'Benign में बीमारी पैदा न करने के सबूत होते हैं, और Pathogenic में बीमारी पैदा करने के। Likely का मतलब है कि भरोसा कुछ कम है। VUS में अभी फैसला करने लायक सबूत नहीं हैं; इसे सुरक्षित मत समझिए।',
 'एक्युरेसी सभी रिकॉर्ड में सही जवाबों का हिस्सा है। बैलेंस्ड एक्युरेसी पाँचों श्रेणियों के रिकॉल का बराबर औसत है। प्रिसिजन पूछता है: दिए गए अनुमान कितनी बार सही थे? रिकॉल पूछता है: असली मामलों में से कितने मिले? दोनों साथ देखना ज़रूरी है।',
 'स्नैपशॉट किसी तारीख की सेव की हुई कॉपी है। सारे रिकॉर्ड ट्रेनिंग में नहीं जाते। सही जीनोम और बदलाव के प्रकार चुने जाते हैं, दोहराव जोड़ा जाता है और टकराते लेबल हटते हैं। मूल बँटवारा लगभग 65, 15, 10 और 10 प्रतिशत ट्रेनिंग, वैलिडेशन, कैलिब्रेशन और टेस्ट है; छँटाई और सैंपलिंग के बाद वास्तविक गिनती अलग है।',
 'ट्रेनिंग में मॉडल सीखता है। अर्ली स्टॉपिंग से रुकने का समय चुनते हैं। वैलिडेशन से उम्मीदवार चुनते हैं। आखिरी टेस्ट निष्पक्ष तभी है जब सुधार करते समय उसे न देखा हो। ट्रेनिंग और वैलिडेशन के पास-पास स्कोर अकेले अच्छे मॉडल का प्रमाण नहीं हैं।',
 'फंक्शनल ऐसे प्रयोग में वेरिएंट का असर मापता है। उसकी माप के तरीके में माइनस स्कोर सही हो सकता है। स्प्लाइस आरएनए के हिस्सों के जुड़ने से जुड़ा है। बहुत से वेरिएंट में ये सबूत नहीं हैं। गायब जानकारी को सामान्य या शून्य नहीं मानना चाहिए।',
 'सुधार के लिए बेहतर जैविक सबूत और साफ लेबल ज़रूरी हैं। नए उम्मीदवार में जीन की पहचान और डीएनए तथा प्रोटीन बदलाव जोड़े गए हैं। प्रिसिजन बढ़ने के साथ रिकॉल गिर सकता है, इसलिए दोनों जाँचिए। सिर्फ ज्यादा रिकॉर्ड लेने से सुधार की गारंटी नहीं है।']
