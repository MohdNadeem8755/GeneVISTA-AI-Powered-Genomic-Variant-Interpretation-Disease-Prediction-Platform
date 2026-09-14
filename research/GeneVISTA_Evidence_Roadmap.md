# GeneVISTA: an evidence-based route to more reliable variant interpretation

The recommended next design is a gene–disease-specific evidence engine, with machine learning supplying calibrated computational evidence and experimental assays supplying a separate source of biological evidence. Begin with a narrowly scoped TP53 germline/Li-Fraumeni pilot, and evaluate a BRCA2 DNA-binding-domain functional-evidence pilot separately. Retain the broad five-class classifier as a research comparator. No existing result demonstrates 99% overall five-class accuracy for GeneVISTA, and neither a new formula nor a larger download can guarantee it.

**1. The actual limitation in the current results**

The saved annotation model report contains 22,105 test records. Of these, 19,130 are labelled uncertain significance. An always-VUS predictor would therefore achieve 86.54% accuracy; the model's 90.69% is 4.14 percentage points higher. This does not make the improvement worthless, but it changes how the headline should be interpreted.

| Recorded class | Test records | Exact-class recall |
|---|---:|---:|
| Pathogenic | 740 | 54.73% |
| Likely pathogenic | 667 | 31.63% |
| VUS | 19,130 | 98.75% |
| Likely benign | 1,246 | 28.57% |
| Benign | 322 | 56.83% |

Balanced accuracy is 54.10% and macro F1 is 0.610. A missed exact class is not necessarily the opposite clinical assertion: for example, a Pathogenic record predicted VUS is an error on this task but is different from falsely declaring it Benign. A future evaluation needs the complete confusion matrix and clinically meaningful error categories rather than treating every error as equally consequential.

These numbers come from the existing model card, not a new evaluation. The expanded eight-gene experiment also failed to establish a specialist advantage: pooled accuracy was 81.16% versus 80.90% for specialists on the same 3,382 validation records. Partial class balancing improved pooled macro F1 to 0.636. Those figures describe a different cohort and are not comparable to the 90.69% test result. Sources: [local annotation model card](../models/annotation_prediction_v2/model_card.json), [specialist experiment](../deployment/SPECIALIST_RESULTS.md).

**2. The five-class target needs more than biological-effect scores**

ClinVar separates variant-level aggregates (VCV) from variant–condition records (RCV). A classification can differ with the condition, and ClinVar explicitly distinguishes these classifications from patient-specific interpretations.[1] GeneVISTA currently starts with a variant, predicts its recorded label, and separately ranks disease annotations. That architecture omits information needed for some disease-specific interpretations.

The proposed unit of analysis is an exact allele plus genome assembly, transcript version, gene, disease and inheritance mechanism. For example, a germline TP53 interpretation for Li-Fraumeni syndrome is a defined task; a universal statement that an allele causes any disease is not. Somatic cancer interpretation must remain separate.

Pathogenic, likely pathogenic, VUS, likely benign and benign summarize evidence and uncertainty, not five mutually exclusive biochemical mechanisms. A VUS can have insufficient or conflicting evidence; it should not be treated as a known harmless example. Consequently, training on functional predictions alone to reproduce all five labels partly asks the model to infer missing clinical evidence. This is an inference from the classification framework and the local feature audit, not a proven irreducible error rate.[2]

A useful auxiliary model may distinguish well-supported pathogenic from benign examples, leaving VUS unlabelled for that training objective. Its binary accuracy must never be reported as five-class accuracy. Similarly, ordinal regression that simply places VUS halfway between damaging and harmless effects is not a justified biological model.

**3. A formula exists, but evidence must justify its inputs**

A Bayesian evidence formulation uses a justified prior p0 and calibrated likelihood ratios LR for applicable evidence:

    posterior_odds = (p0 / (1 - p0)) × LR1 × LR2 × ... × LRk
    posterior_pathogenicity = posterior_odds / (1 + posterior_odds)

The multiplication assumes conditional independence or a valid joint treatment of dependent observations. It cannot be applied to arbitrary scores or percentages. The published ACMG Bayesian formulation provides approximate evidence multipliers of 2.08, 4.33, 18.7 and 350 for supporting through very-strong pathogenic evidence.[3]

As a mathematical illustration only, prior 0.10 and genuinely independent likelihood ratios 350 and 4.33 yield a posterior of approximately 99.41%. This does not demonstrate 99.41% model accuracy or a 99.41% chance of disease. The likelihood ratios, their independence, the prior, the disease mechanism and evidence applicability all require justification. The illustration is not a classification rule for any real variant.

The related points approach represents log-odds using positive and negative evidence points. Its original thresholds embed particular assumptions, including a prior near 0.10.[4] Production code should implement the applicable versioned expert-panel specification, including exceptions, rather than copy generic thresholds. For example, the TP53 CSpec page inspected for this report displays version 2.4 and treats a total of −1 as VUS by default, with a defined exception; this differs from a simple reading of the original generic points table.[5]

Do not multiply REVEL, AlphaMissense, CADD and a GeneVISTA model score as independent evidence. Their inputs overlap, and the GeneVISTA model itself includes population-frequency features. Counting its output and gnomAD again would duplicate information. ClinGen recommends calibrated computational evidence and careful treatment of dependencies; adding more agreeing predictors is not automatically additional evidence.[6]

**4. The best concrete evidence pilots**

TP53/Li-Fraumeni is a practical first implementation because the expert-panel publication supplies functional, computational and splice evidence spreadsheets, while CSpec supplies explicit gene–disease rules. The publication's 43-variant pilot reached clinically meaningful classifications for 93% of variants; that is classification yield, not 93% accuracy. Public example cases can verify faithful rule implementation but are not an independent test when they helped develop the rules.[7]

Pin the current TP53 CSpec release and inspect its supporting files before coding. Reconcile transcript versions: the publication discusses NM_000546.6, while the registry metadata currently lists NM_000546.5. Do not silently combine coordinates or protein changes across accessions. Current rules also constrain combining splice and protein-function evidence.[5]

A second promising pilot is BRCA2's C-terminal DNA-binding domain. A 2026 primary study combined two saturation-genome-editing datasets and reported 98.8% accuracy for an Integrated VarCall functional model. It processed 6,383 variants, but that is not the independent accuracy-test denominator. One reported comparison used 158 ClinVar missense standards whose classifications were not influenced by the SGE results, with 98.9% sensitivity and 100% specificity under specified functional categories. The abstract's 98.8% and that comparison are different reported endpoints; they should not be reconstructed into one five-class accuracy claim.[8]

This study is a research precedent for adding measured functional evidence in a limited biological setting. It does not establish that the same performance transfers to GeneVISTA, all genes, all variant types or patient disease risk. It is a strong reason to test functional-data integration before buying more computing power.

**5. What is already available and what to obtain**

A read-only inspection of the downloaded dbNSFP5.4a header found 508 columns. Current training imports 11 annotation features. Available but unused fields include ESM1b_score, popEVE_score, MANE, Ensembl_transcriptid, Interpro_domain, MPC_score, VEST4_score, MutPred2_score and BayesDel_noAF_score. These are additional predictors or annotation context, not independent laboratory experiments. Their presence was verified from the local file; per-gene coverage has not yet been quantified.

| Resource | Purpose | Acquisition and scope |
|---|---|---|
| Existing dbNSFP5.4a | Transcript-aware protein-effect features and functional-context predictors | Already in D:\GeneVISTA\data\incoming; selectively extract suitable columns. No new ESM download is needed for fields present here. |
| TP53 functional evidence, Table S3 | Assay measurements and preliminary PS3/BS3 assignments | Download the linked XLSX below; reconcile it with the current CSpec supporting files. |
| TP53 computational/splice evidence, Tables S2/S1 | Applicable computational and splice evidence | Download the linked XLSX files below; do not sum correlated criteria independently. |
| ClinGen CSpec and Evidence Repository | Rules, versioned assertions and supporting evidence | Search TP53 or BRCA2 and export the relevant records; preserve assertion and rule versions. |
| BRCA2 Integrated VarCall study | Measured functional evidence in the BRCA2 DBD | Use the paper's Supplementary Information/Data availability links, preserving raw assay provenance and control-set membership. |
| MaveDB | Published multiplexed assay measurements | Search by gene and assay, then export a relevant public score set; score meaning and assay applicability must be checked.[9] |
| ClinVar RCV/SCV information | Disease-specific assertions, submitter evidence and evaluation dates | Use selected records first; full RCV XML is available through NCBI's documented download route.[10] |
| Independent clinical cases | Clinical validation and patient-level disease modelling | Requires a suitable laboratory or research cohort with appropriate permissions; public annotation tables do not supply this endpoint. |

Direct TP53 files linked by the primary paper:

- [Functional evidence: Table S3 XLSX](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13073-025-01536-3/MediaObjects/13073_2025_1536_MOESM3_ESM.xlsx)
- [Computational evidence: Table S2 XLSX](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13073-025-01536-3/MediaObjects/13073_2025_1536_MOESM2_ESM.xlsx)
- [Splice evidence: Table S1 XLSX](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13073-025-01536-3/MediaObjects/13073_2025_1536_MOESM1_ESM.xlsx)

These links were verified from the publication; workbook contents have not yet been imported or audited locally. The current CSpec attachments take precedence if they update the published supplements. Record licences and permissions per dataset before redistribution; a public download is not a universal licence for clinical or commercial use.

A concrete public BRCA1 dataset is MaveDB score set urn:mavedb:00000097-0-2. Its metadata describes 3,893 variants, source transcript NM_007294.3, functional/RNA measurements and a CC0 licence. Metadata and documented public API access were checked; the score CSV has not yet been imported. Public reads require no API key.[9]

- [BRCA1 score-set metadata](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000097-0-2)
- [BRCA1 functional/RNA scores CSV](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000097-0-2/scores)
- [BRCA1 mapped variants](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000097-0-2/mapped-variants)

This is genuinely measured evidence beyond the unimported dbNSFP predictors. Its transcript, calibration versions and possible prior contribution to ClinVar labels must be reconciled. Normal assay function alone does not establish clinical benignity.
For the existing file, the next ML experiment should preserve the transcript–score relationship rather than take a maximum across unrelated transcripts. Check score direction, missingness and calibration; a more negative score may mean greater damage for some models. Exclude ClinVar classifications, assertion dates, clinical identifiers and other target-derived fields from predictive features. MANE offers standardized transcript correspondence, but clinically relevant alternative transcripts may still require consideration.[11]

**6. A practical experimental design**

First define the pilot's intended output: an evidence-supported TP53 germline classification for Li-Fraumeni syndrome, with missing evidence and review requirements visible. Keep the existing general classifier as a separately labelled research estimate. This changes the architecture while preserving the five-class goal.

Build an evidence ledger per variant–disease pair: source accession, exact allele, transcript, assay, evidence date, criterion, calibrated strength, dependency group and reviewer status. Record 'not evaluated' separately from evidence against pathogenicity. A missing family history must not become negative segregation evidence.

Validate applicable functional evidence using the disease mechanism and the specific assay's controls. Functional evidence does not automatically qualify as strong merely because a score exists.[12] Use a gene-specific rules engine to combine valid criteria and block prohibited combinations. Independently test allele mapping, rule boundaries, source-version mismatches, contradictory evidence and missing fields.

In parallel, conduct a bounded ML ablation using the current file: existing 11 features; transcript-aligned ESM1b/popEVE and context additions; and a supervised-predictor-free comparator. A model such as regularized gradient boosting is a reasonable comparator. A deeper neural network is not the key hypothesis until additional information shows incremental benefit. The experiment must distinguish gains from new evidence, extra training examples and changes in label selection.

Do not choose a feature or threshold using the already evaluated test set. Register a new protocol and lock all model, criterion, calibration and decision choices before a genuinely independent evaluation. A historical split is not automatically leakage-free if modern source predictors were trained on later labels. Source model training lists, assay controls, overlapping codons, related families and assertion dates require an overlap audit. If independence cannot be established, label the result retrospective research concordance.

A successful rule-reproduction test is evidence that the software follows the specification. It is not independent evidence that the specification or model is clinically accurate. Report these two evaluations separately.

**7. What a credible 99% claim would require**

Specify the denominator and error first: five-class exact agreement, sensitivity for pathogenic variants, predictive value of accepted calls, or patient diagnosis. These are different endpoints. Specify which genes and variant types are in scope and how uncertain or unsupported inputs are handled.

For a system allowed to abstain, report both selective error and coverage:

    selective_error = incorrect accepted predictions / accepted predictions
    coverage = accepted predictions / all eligible inputs

A model can have excellent accuracy on a small accepted subset while refusing most variants. That can be useful, but it is not 99% overall accuracy. A threshold chosen on validation data must be frozen and checked independently; choosing only successful test cases invalidates the claim.

Confidence bounds matter. With zero errors among n independent Bernoulli trials, an exact one-sided 95% upper error bound is 1 − 0.05^(1/n). At n=299 this is approximately 0.997%. This is a mathematical illustration under specific sampling assumptions, not a sufficient clinical study design.[13] Related variants, selected easy cases, subgroup differences and repeated tuning can invalidate that simple interpretation. Separate rare-class and clinical-harm endpoints need adequate sample sizes of their own.

Conformal prediction can be investigated for prediction sets or abstention. Its guarantees depend on the procedure and data assumptions; marginal set coverage does not equal class-specific clinical sensitivity or 99% correctness of singleton answers. Distribution shift remains a material issue. Selective classification and conformal risk control provide formal research tools for defining such goals, with different guarantees that must not be conflated.[15][16]

A paper reporting AUROC 0.99 has not necessarily reported 99% accuracy. For example, a recent dbNSFP/XGBoost study reports AUROC 0.991 on a pathogenic-versus-benign benchmark; it is neither GeneVISTA's five-class task nor evidence of transferable 99% accuracy.[14]

**8. Disease prediction is a separate endpoint**

Variant pathogenicity does not determine whether a particular person develops a disease. Personal risk depends on genotype, zygosity, inheritance, penetrance, age, family and clinical context, among other variables. ClinVar explicitly treats its classifications as variant-level rather than patient-specific.[1]

GeneVISTA can present disease associations and, with suitable evidence, prioritize hypotheses for review. A value such as '78% chance of disease' needs a model fitted and validated on linked patient outcomes, with an explicit disease and time horizon. Ranking scores must not be converted to risk percentages by normalization or by multiplying a pathogenicity score with a gene–disease association score.

**9. Recommended next deliverable**

Build a versioned TP53 evidence prototype and a transparent dataset manifest on D:, together with a BRCA2 functional-evidence feasibility check. Obtain the small public supplemental datasets first. Quantify how many exact variants have sufficient applicable evidence, and distinguish known expert assertions, rule-derived research classifications and unresolved inputs.

The first acceptance criteria should be correct mapping, reproducible evidence rules, explicit dependencies and error accounting. The next gate is external expert assessment on independent cases. Only measured evaluation can decide whether any carefully scoped output meets a 99% criterion. Broader deployment should follow that evidence, rather than a manually selected accuracy target.

**Sources**

1. NCBI. [Representation of classifications in ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/). Official documentation accessed September 2026.
2. Richards et al. (2015). [Standards and guidelines for the interpretation of sequence variants](https://www.nature.com/articles/gim201530). ACMG/AMP consensus guideline.
3. Tavtigian et al. (2018). [Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336098/).
4. Tavtigian et al. (2020). [Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines](https://pmc.ncbi.nlm.nih.gov/articles/PMC8011844/).
5. ClinGen TP53 VCEP. [Criteria Specification Registry, TP53, version 2.4 displayed at access](https://cspec.genome.network/cspec/ui/svi/svi/GN009). Accessed September 2026.
6. Pejaver et al. (2022). [Calibration of computational tools and ClinGen PP3/BP4 recommendations](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/).
7. Fortuno et al. (2025). [Updated TP53 expert-panel recommendations for Li-Fraumeni syndrome](https://link.springer.com/article/10.1186/s13073-025-01536-3). Includes Tables S1–S3.
8. Hu et al. (2026). [Combining multiplexed assays of variant effect for enhanced BRCA2 variant classification](https://www.nature.com/articles/s41467-026-71393-0). [Public full-text record](https://pmc.ncbi.nlm.nih.gov/articles/PMC13223280/).
9. [MaveDB](https://www.mavedb.org/), [public API instructions](https://www.mavedb.org/docs/mavedb/programmatic-access/api-quickstart.html), and [BRCA1 score-set metadata](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000097-0-2). Dataset metadata and access verified; inspect each score set's calibration provenance and licence.
10. NCBI. [Guide to ClinVar FTP files and e-utilities](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/). Includes variant–condition RCV XML downloads.
11. NCBI/EMBL-EBI. [Matched Annotation (MANE)](https://www.ncbi.nlm.nih.gov/refseq/MANE/).
12. ClinGen. [Recommendations for applying PS3/BS3 functional evidence](https://www.clinicalgenome.org/docs/recommendations-for-application-of-the-functional-evidence-ps3-bs3-criterion-using-the-acmg-amp-sequence-variant-interpretation/).
13. NIST. [Exact binomial confidence limits](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/exacbici.htm). The zero-error example above is a derived calculation, not a clinical recommendation.
14. MetaXVP authors (2026). [MetaXVP: an interpretable machine learning framework](https://www.nature.com/articles/s41598-026-50032-0). Primary report; binary benchmark metrics should not be read as five-class accuracy.

15. Geifman and El-Yaniv (2017). [Selective Classification for Deep Neural Networks](https://papers.neurips.cc/paper_files/paper/2017/hash/4a8423d5e91fda00bb7e46540e2b0cf1-Abstract.html).
16. Angelopoulos et al. [Conformal Risk Control](https://arxiv.org/abs/2208.02814). Research method; guarantees depend on the loss and sampling assumptions.
Local evidence inspected: models/annotation_prediction_v2/model_card.json; reports/annotation_v3/specialists_variable_features/comparison.json; src/annotation_pipeline.py; and the header of data/incoming/dbNSFP5.4a_grch38.gz. No new model was trained or promoted during this investigation.
