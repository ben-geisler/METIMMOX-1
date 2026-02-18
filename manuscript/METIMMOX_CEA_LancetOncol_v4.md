# Cost-effectiveness of biomarker-guided chemo-immunotherapy in microsatellite-stable metastatic colorectal cancer: a model specification analysis from the METIMMOX trial

**Original Research**

Benjamin P Geisler^1,2^, Mathyn A M Vervaart^1,3^, Sebastian Meltzer^2^, Paula A Bousquet^2^, Emily A Burger^1,4^, Anne H Ree^2,5^, Eline Aas^1,6^

1. Department of Health Management and Health Economics, University of Oslo, Oslo, Norway
2. Department of Oncology, Akershus University Hospital, Lørenskog, Norway
3. Clinical Trials Unit, Oslo University Hospital, Oslo, Norway
4. Center for Health Decision Science, Harvard T H Chan School of Public Health, Boston, MA, USA
5. Institute of Clinical Medicine, University of Oslo, Oslo, Norway
6. Division for Health Services, Norwegian Institute of Public Health, Oslo, Norway

**Correspondence to:** Benjamin P Geisler, Department of Health Management and Health Economics, University of Oslo, PO Box 1089, Blindern, 0317 Oslo, Norway. E-mail: [email]

**Word count:** [~4,500 excluding heading, summary, acknowledgements, funding, and references]

---

## Summary

**Background** Immune checkpoint blockade is ineffective in unselected microsatellite-stable/mismatch repair-proficient (MSS/pMMR) metastatic colorectal cancer (mCRC). The METIMMOX trial identified exploratory biomarker signals suggesting differential immunotherapy benefit. We assessed the cost-effectiveness of biomarker-guided chemo-immunotherapy and whether survival model specification influences economic conclusions.

**Methods** Using individual patient data from METIMMOX (NCT03388190; n=65), we developed a partitioned survival model comparing three biomarker-guided strategies (C-reactive protein [CRP], tumor mutational burden/*BRAF* [TMB/*BRAF*], target lesion reduction [TLR]) against standard FLOX chemotherapy. Gamma parametric models were fitted under three specifications varying in which biomarkers entered as covariates and interaction terms, including principled exclusion of a post-treatment prognostic variable. The base-case analysis adopted a Norwegian healthcare perspective over a ten-year horizon. A secondary analysis evaluated biomarker-enriched populations.

**Findings** Under the primary specification (Model B), standard care yielded 2·08 quality-adjusted life-years (QALYs) at €21,504 and dominated all biomarker strategies, which provided fewer QALYs at higher cost (incremental QALYs: CRP −0·18, TMB/*BRAF* −0·14, TLR −0·38). This was consistent across specifications; only one showed a small positive effect for TMB/*BRAF* (ICER €1,690,480/QALY). In biomarker-enriched populations, all strategies remained dominated under Model B, with the model predicting net harm even in biomarker-positive patients (QALY deficits −0·27 to −0·54). The most favorable result was TMB/*BRAF*-enriched under Model C (ICER €323,878/QALY), still exceeding the threshold sixfold. Control QALYs varied 12·7% across specifications. Probabilistic analysis showed >98% probability that standard care is optimal at €51,000/QALY.

**Interpretation** Biomarker-guided chemo-immunotherapy in MSS/pMMR mCRC is not cost-effective regardless of model specification or target population. The parametric model predicts net harm from immunotherapy even in biomarker-positive patients, meaning price reduction alone cannot establish cost-effectiveness; prospective validation of clinical benefit is the binding constraint. Exploratory interaction signals from semi-parametric Cox models did not translate into QALY gains under parametric survival modeling.

**Funding** European Union (ASCERTAIN, 101094938), Research Council of Norway (NORCHER, 296114).

---# Research in context

### Evidence before this study

We searched PubMed and Embase from inception to October 2025 for cost-effectiveness analyses of biomarker-guided immunotherapy strategies in microsatellite-stable (MSS) metastatic colorectal cancer, using the terms "cost-effectiveness", "biomarker", "immunotherapy" or "checkpoint inhibitor", and "colorectal". While economic evaluations exist for immunotherapy in MSI-H/dMMR colorectal cancer and for biomarker-guided immunotherapy in other tumor types (e.g., PD-L1 testing in non-small-cell lung cancer), we identified no published cost-effectiveness analyses of biomarker-guided checkpoint inhibitor strategies specifically in MSS mCRC. Moreover, the methodological question of how survival model specification—particularly the choice of which candidate biomarkers to include as covariates versus interaction terms—influences economic conclusions in multi-biomarker settings has not been systematically examined in the health technology assessment literature.

### Added value of this study

This study provides, to our knowledge, the first economic evaluation of biomarker-guided immunotherapy in MSS mCRC, a disease affecting the large majority of patients with metastatic colorectal cancer who currently derive no benefit from checkpoint blockade. We demonstrate that model specification choices—specifically, how multiple candidate biomarkers are incorporated into survival models—can substantially alter clinical effectiveness estimates, with control strategy QALYs varying by 12·7% (1·84–2·08) across specifications. Including or excluding a post-treatment prognostic variable (early tumor shrinkage) as a covariate changes baseline survival estimates, with more modest effects on relative treatment comparisons given weak biomarker intercorrelations (Cramér's V = 0·278). Critically, we show that the exploratory treatment–biomarker interaction signals from the companion clinical analysis—including CRP–treatment interaction hazard ratios of 0·19–0·22 for PFS—did not translate into QALY gains when incorporated into parametric survival models with full extrapolation, and that biomarker-enriched population analysis confirms this reflects genuine lack of predicted treatment benefit rather than prevalence dilution. The analysis uses expected value of information to demonstrate that prospective biomarker validation has negligible value at current immunotherapy prices.

### Implications of all the available evidence

At current immunotherapy prices, biomarker-guided strategies for selecting MSS mCRC patients for chemo-immunotherapy do not represent cost-effective healthcare resource use. Under the primary model specification, no biomarker strategy produces a QALY gain over standard chemotherapy—even in biomarker-enriched populations—meaning that neither price reduction nor improved patient selection can achieve cost-effectiveness without first establishing genuine clinical benefit through prospective validation. Methodologically, this study establishes that economic evaluations of precision oncology strategies derived from exploratory trial analyses should routinely assess sensitivity to model specification, including the treatment of post-treatment prognostic variables, and that apparently promising interaction signals from semi-parametric Cox models may not translate into predicted QALY gains under parametric survival modeling.

---

## Introduction

Most patients with metastatic colorectal cancer (mCRC) derive no benefit from immune checkpoint blockade (ICB). While the 4–8% with microsatellite-instable or mismatch repair-deficient (MSI-H/dMMR) tumors experience durable responses,^1–4^ the remaining majority with microsatellite-stable, mismatch repair-proficient (MSS/pMMR) disease face a therapeutic ceiling: oxaliplatin-based regimens yield median overall survival rarely exceeding 30 months, and precision oncology approaches have largely failed to improve outcomes in this population.^5–7^

The METIMMOX trial (NCT03388190) tested whether alternating Nordic FLOX chemotherapy with nivolumab could overcome immunotherapy resistance in treatment-naive MSS/pMMR mCRC.^8^ Although the intention-to-treat analysis was negative, exploratory biomarker analyses identified potential benefit signals in three subgroups: patients with low C-reactive protein (CRP <5 mg/L),^9^ those with early tumor shrinkage (≥10% target lesion reduction [TLR]),^10^ and those with high tumor mutational burden (≥9 mutations/Mb) or *BRAF* V600E mutation.^11^ An accompanying analysis using Firth-corrected Cox regression established that TLR functions primarily as a prognostic rather than predictive biomarker, while CRP showed a suggestive predictive signal for progression-free survival (PFS).^12^

A fundamental challenge arises when translating multiple candidate biomarker signals from a single trial into an economic evaluation. Each biomarker may be prognostic, predictive, or both; baseline-measured or dynamic; and correlated with or independent of the others. The analyst must decide which biomarkers to include in the survival model, whether each enters as a main effect (prognostic covariate), an interaction term (predictive effect), or both, and whether post-treatment dynamic biomarkers should be conditioned upon at all. Each specification embeds assumptions about biomarker function and interdependence, yet the health technology assessment (HTA) literature offers no framework for navigating these choices. The result is that modeling decisions with substantial influence on cost-effectiveness conclusions are made implicitly rather than systematically evaluated.

A further complexity concerns the target population. When a biomarker-guided strategy is evaluated from a health-system perspective, the relevant cohort is all MSS mCRC patients who would undergo testing, with outcomes weighted by biomarker prevalence. When evaluated from a clinical decision-making perspective, the cohort is restricted to biomarker-positive patients for whom treatment selection is at stake. These perspectives yield different cost-effectiveness conclusions and serve different audiences—HTA bodies versus treating clinicians—yet are rarely compared explicitly.

We conducted an early economic evaluation of biomarker-guided strategies for alternating chemo-immunotherapy in MSS/pMMR mCRC. Our objectives were fourfold: first, to assess whether any biomarker-guided strategy is cost-effective compared with standard chemotherapy; second, to establish a framework for comparing survival model specifications when multiple candidate biomarkers emerge from exploratory trial analyses; third, to compare cost-effectiveness across overall MSS and biomarker-enriched populations; and fourth, to use value-of-information analysis to inform decisions about prospective biomarker validation.

---

## Methods

### Clinical data source

METIMMOX was an investigator-initiated, multicenter, open-label, randomized phase II trial conducted across five Norwegian hospitals between 2018 and 2021.^8^ Eligible adults had previously untreated, unresectable MSS/pMMR metastatic colorectal adenocarcinoma with at least one measurable infradiaphragmatic metastasis and performance status 0–1. Seventy-four patients were enrolled; 65 with complete biomarker data were included in this analysis. Three candidate biomarkers were defined post hoc: CRP-positive (CRP <5·0 mg/L; prevalence 33·8%), TLR-positive (≥10% target lesion reduction at first assessment; 63·1%), and TMB/*BRAF*-positive (TMB ≥9·0 mut/Mb or *BRAF* V600E; 44·6%).^9–11^

### Survival model specification framework

Parametric survival models were fitted to individual patient data for PFS and overall survival (OS). Candidate distributions (exponential, Weibull, Gompertz, log-normal, log-logistic, gamma, generalized gamma) were evaluated using Akaike and Bayesian information criteria, visual fit to Kaplan–Meier curves, and clinical plausibility of extrapolations beyond observed follow-up. Gamma distributions were selected for both endpoints and applied consistently across all model specifications.

We compared three model specifications representing distinct assumptions about how biomarkers should enter the survival model that drives the economic evaluation (panel).

> **Panel: Survival model specifications**
>
> **Model A (comprehensive)** includes all three biomarkers as main effects and all three treatment–biomarker interactions. For a given biomarker-guided strategy, the treatment effect in biomarker-positive patients is estimated while simultaneously adjusting for the main and interaction effects of the other two biomarkers. This specification maximizes covariate control but estimates seven biomarker-related parameters from 65 patients, risking overfitting.
>
> **Model B (focused)** includes CRP and TMB/*BRAF* as main effects and both treatment–CRP and treatment–TMB/*BRAF* interaction terms, but excludes TLR entirely—both as a main effect and as a potential interaction term. A single model is fitted that estimates both interaction effects simultaneously.
>
> The rationale for excluding TLR is twofold. First, the companion clinical analysis established that TLR is prognostic rather than predictive: it identifies patients with favorable tumor biology regardless of treatment arm but does not identify patients who specifically benefit from immunotherapy.^12^ A prognostic-only biomarker does not justify a biomarker-guided treatment strategy. Second, TLR is a dynamic biomarker measured after treatment initiation. Conditioning on post-treatment variables when estimating treatment effects can introduce bias by blocking causal pathways or inducing collider stratification—a concern well-established in causal inference^13,14^ but not addressed in HTA guidance for survival model covariate selection.^15^ Model B tests whether these principled exclusions alter economic conclusions.
>
> **Model C (single-biomarker)** includes only the biomarker relevant to the strategy being evaluated, with its treatment interaction. This mirrors conventional subgroup-analysis approaches but risks conflating omitted biomarkers' prognostic effects with the treatment effect of interest.

Under each specification, we estimated separate models for CRP-guided and TMB/*BRAF*-guided strategies. The TLR-guided strategy was evaluated under Models A and C only (where TLR is included by construction), because demonstrating the economic consequences of treating a prognostic biomarker as if it were predictive is itself a methodological contribution. All models were adjusted for age and sex.

### Decision-analytic model

A partitioned survival model tracked a hypothetical cohort through three mutually exclusive health states—progression-free, post-progression, and dead—over a ten-year time horizon using one-week cycles. Health-state occupancy at each time point was derived from the fitted parametric survival curves: the proportion progression-free equaled the PFS survival function; the proportion deceased equaled one minus the OS survival function; and the proportion post-progression was the difference.

Four strategies were compared: **standard of care** (all patients receive FLOX alone); **CRP-guided** (CRP-positive patients receive FLOX–nivolumab; CRP-negative receive FLOX alone); **TMB/*BRAF*-guided** (TMB/*BRAF*-positive patients receive FLOX–nivolumab; others receive FLOX alone); and **TLR-guided** (TLR-positive patients receive FLOX–nivolumab; TLR-negative receive FLOX alone). A treat-all immunotherapy strategy was not evaluated because the intention-to-treat analysis was negative.

### Population framework

The **base-case analysis** adopted an overall MSS population perspective: the cohort comprised all MSS/pMMR mCRC patients, with strategy-specific survival curves derived by predicting outcomes for biomarker-positive patients under the experimental treatment and biomarker-negative patients under standard care, weighted by observed biomarker prevalence. This approach, recommended by NICE Decision Support Unit guidance,^15^ represents the health-system-level decision of whether to implement a biomarker testing program.

A pre-specified **secondary analysis** evaluated **biomarker-enriched populations**: the cohort was restricted to biomarker-positive patients, comparing immunotherapy versus standard chemotherapy within this subgroup only. This analysis represents the clinical decision for a patient already known to be biomarker-positive and establishes a ceiling estimate—the maximum possible value of the biomarker strategy if testing were cost-free and perfectly accurate.

### Economic evaluation

The analysis adopted an extended healthcare perspective consistent with Norwegian HTA guidelines. Costs and health outcomes were discounted at 4% annually. The willingness-to-pay (WTP) threshold was set at NOK 605,000 per quality-adjusted life-year (QALY), approximately €51,000, based on the severity-adjusted threshold for metastatic colorectal cancer.^16^

Direct medical costs included treatment acquisition (nivolumab and FLOX components, assuming a 40% discount reflecting regional procurement), diagnostic testing (CT imaging for all strategies; next-generation sequencing at approximately €2,500 for the TMB/*BRAF*-guided strategy), disease monitoring, subsequent therapy, and end-of-life care. All costs are reported in 2023 euros, converted from Norwegian kroner using purchasing power parity. Quality-adjusted life-years were calculated using EQ-5D-5L utility values from METIMMOX, converted using Danish tariffs,^17^ stratified by progression status. Table 1 presents all input parameters with sources and distributional assumptions.

### Sensitivity and uncertainty analyses

One-way deterministic sensitivity analysis varied key parameters across plausible ranges. Parameter uncertainty was propagated through 5,000 second-order Monte Carlo simulations, with survival model parameters sampled from multivariate normal distributions preserving the variance–covariance structure, cost parameters from gamma distributions, and utility values from beta distributions. Results are presented as cost-effectiveness acceptability curves (CEACs).

Expected value of partial perfect information (EVPPI) was estimated for parameter groups including baseline survival parameters, biomarker prognostic effects, drug acquisition costs, and health-state utilities. EVPPI for treatment–biomarker interaction effects—theoretically the most policy-relevant group—is being calculated separately and will be reported in the final version. EVPPI was calculated at the base-case WTP threshold (€51,000/QALY), under a biosimilar pricing scenario (80% discount from list price, equivalent to €4,641 per administration), and at an elevated threshold (€100,000/QALY). Estimates used a two-level Monte Carlo approach with 1,000 outer-loop and 500 inner-loop samples and were scaled to the eligible Norwegian patient population over a [10]-year technology horizon.

### Model specification comparison

To assess how model specification influenced conclusions, we compared results across Models A, B, and C on four dimensions: clinical effectiveness estimates (median survival, mean life-years, QALYs), base-case cost-effectiveness (ICERs, efficiency frontier composition, incremental net monetary benefit), uncertainty (width of probabilistic scatter, probability of cost-effectiveness), and decision concordance (whether all specifications yield the same policy recommendation at the reference WTP).

The comparison between Model A and Model B specifically isolates two effects: the consequence of excluding a post-treatment prognostic variable (TLR) and the consequence of reducing the number of estimated interaction terms from three to two. The comparison between Model B and Model C isolates the effect of including versus excluding TMB/*BRAF* and CRP as mutual prognostic covariates.

The trial protocol was approved by the Regional Committee for Medical and Health Research Ethics (REK South-East, reference [XX]). All patients provided written informed consent. This economic evaluation was conducted using anonymized individual patient data under the original ethical approval.

---

## Results

### Parametric survival models

Gamma distributions provided the best fit for both PFS and OS based on AIC (PFS: 453·43; OS: 622·56), BIC, and visual assessment (appendix pp [XX–XX]). Extrapolations beyond observed follow-up were clinically plausible: predicted 5-year OS of 6·1% for standard care under Model A, consistent with published registry data for MSS mCRC.^5^

Table 2 (appendix) presents the fitted parametric model parameters across specifications. The treatment–biomarker interaction coefficients varied modestly across model specifications. For the CRP–treatment interaction in PFS models, the coefficient changed from −0·58 (SE 0·24) in Model A to a similar magnitude in Model B, suggesting that TLR exclusion did not substantially alter the estimated CRP predictive effect. This stability reflects the weak correlation between CRP and TLR (Cramér's V = 0·278), limiting confounding between these biomarkers.

### Clinical effectiveness

Table 3 (appendix) presents clinical effectiveness estimates by strategy and model specification. Under Model B, standard care yielded an estimated 2·08 mean discounted QALYs, while the CRP-guided strategy yielded 1·90 QALYs, a difference of −0·18 QALYs (worse than standard care). The TMB/*BRAF*-guided strategy yielded 1·94 QALYs (−0·14 incremental). The TLR-guided strategy (evaluated under Models A and C only) yielded 1·70 QALYs, substantially less than standard care (−0·38), consistent with TLR's role as a prognostic rather than predictive biomarker: TLR-negative patients, who receive only chemotherapy under this strategy, have substantially worse baseline prognosis, and their poor outcomes are not offset by gains in TLR-positive patients.

Clinical effectiveness estimates varied substantially across model specifications. Control strategy QALYs ranged from 1·84 (Model A) to 2·08 (Model B), a 12·7% difference reflecting the impact of model formula choices on baseline survival predictions. For biomarker strategies, the pattern was consistent: all three biomarker-guided strategies yielded fewer QALYs than control in Models B and C. Only in Model A did one strategy (TMB/*BRAF*) show a small positive incremental effect (+0·03 QALYs). The comparison between Model A (which includes TLR as a covariate) and Model B (which excludes it) showed that excluding TLR increased estimated control QALYs from 1·84 to 2·08, while biomarker strategy estimates changed less dramatically. The modest relative impact of TLR exclusion on CRP estimates reflects the weak correlation between CRP and TLR (Cramér's V = 0·278), limiting the scope for confounding from TLR's inclusion or exclusion.

### Base-case cost-effectiveness: overall MSS population

Table 4 presents costs, QALYs, and ICERs for all strategies under Model B. Standard care was both the least costly (€21,504) and most effective (2·08 discounted QALYs) option. All three biomarker-guided strategies were dominated, providing fewer QALYs at higher cost: CRP-guided cost €55,118 (1·90 QALYs, −0·18 vs control); TMB/*BRAF*-guided cost €61,542 (1·94 QALYs, −0·14 vs control); TLR-guided cost €86,411 (1·70 QALYs, −0·38 vs control). These QALY deficits indicate that the exploratory treatment–biomarker interaction signals from the companion clinical analysis—which were suggestive but not statistically significant after false discovery rate correction^12^—did not translate into net survival gains when incorporated into gamma parametric survival models with full extrapolation over a ten-year horizon. As confirmed by the enriched-population analysis below, the model predicts that biomarker-positive patients fare worse on immunotherapy than on chemotherapy alone, rather than this being an artifact of including biomarker-negative patients in the population average.

The efficiency frontier was concordant across model specifications. Under Models B and C, standard care strictly dominated all biomarker strategies. Under Model A, the TMB/*BRAF*-guided strategy generated a small positive incremental effect (0·03 QALYs) with an ICER of €1,690,480 per QALY—far exceeding any plausible willingness-to-pay threshold. The qualitative conclusion that no biomarker strategy is cost-effective at €51,000/QALY was robust to model specification.

The Model A versus Model B comparison revealed that excluding TLR and reducing interaction terms altered the clinical effectiveness estimates substantially for control (from 1·84 to 2·08 QALYs, a 13% increase) while biomarker strategy estimates changed less. This occurred because TLR—a post-treatment prognostic variable—was correlated with outcomes under standard care more strongly than with biomarker-treatment interactions. The modest influence of TLR exclusion on relative treatment effects suggests that in settings with weak biomarker intercorrelations (Cramér's V = 0·278), the practical impact of conditioning on post-treatment variables may be limited, though the principle remains important for future evaluations where correlations are stronger.

### Secondary analysis: biomarker-enriched populations

When restricted to biomarker-positive patients only (comparing immunotherapy versus standard chemotherapy within each biomarker-defined subgroup), results remained unfavorable across most model specifications. This ceiling scenario analysis—which removes any influence of biomarker-negative patients on population-weighted outcomes—represents the best-case economic outcome assuming perfect biomarker accuracy with no additional testing costs.

Under Model B (primary analysis), all three biomarker-guided strategies remained dominated even in enriched populations: CRP-positive patients receiving experimental treatment cost €121,746 with 2·34 QALYs versus €22,347 with 2·87 QALYs if they received standard chemotherapy only (dominated; −0·54 QALYs); TLR-positive patients showed €122,839 with 2·08 QALYs versus €22,870 with 2·35 QALYs on control (dominated; −0·27 QALYs); TMB/*BRAF*-positive patients showed €107,775 with 1·89 QALYs versus €24,265 with 2·19 QALYs on control (dominated; −0·31 QALYs). Since biomarker-negative patients are excluded from this analysis entirely, these results confirm that the dominated findings in the overall population reflect the parametric model’s prediction of net treatment harm in biomarker-positive patients, not prevalence dilution or adverse selection from including prognostically unfavorable subgroups.

Under Model A, the TMB/*BRAF*-enriched population generated a small positive incremental effect (0·055 QALYs) with an ICER of €1,625,723 per QALY, representing a 3·8% reduction from the overall population ICER of €1,690,480 but remaining far above any plausible threshold. CRP and TLR remained dominated. Under Model C, TLR-enriched yielded an ICER of €1,242,396 per QALY and TMB/*BRAF*-enriched €323,878 per QALY—the most favorable result observed, yet still exceeding the €51,000 threshold sixfold.

The enriched-population findings establish that improving patient selection cannot solve the fundamental cost-effectiveness problem: even under the most optimistic assumptions (perfect biomarker targeting, no testing costs), ICERs remain at least sixfold above the Norwegian threshold, and in the primary Model B analysis, biomarker-positive patients actually fare worse on experimental treatment than on standard chemotherapy.

### Sensitivity analyses

One-way sensitivity analysis identified nivolumab acquisition cost as the dominant driver of uncertainty across all model specifications and strategies, generating approximately €20,000–40,000 variation in net monetary benefit per biomarker strategy (figure [X]). This single parameter's influence exceeded all other inputs by an order of magnitude. Because biomarker strategies yield fewer QALYs than standard care under Models B and C, cost reduction alone cannot achieve cost-effectiveness under these specifications—cost-effectiveness requires positive incremental QALYs. Under Model A, where TMB/*BRAF* shows a small positive incremental effect (+0·03 QALYs), threshold analysis indicated that nivolumab prices would need to decrease to approximately €4,641 per administration (an 80% discount from list price, equivalent to a 67% reduction from the base-case procurement price) for this single strategy to potentially approach cost-effectiveness.

Probabilistic sensitivity analysis (5,000 iterations) confirmed that standard care maintained the highest probability of cost-effectiveness at all conventional thresholds under all model specifications (figure [X]). Under Model B at €51,000/QALY, standard care had a 98·5% probability of being optimal, with biomarker strategies collectively accounting for only 1·5%. Under Model A (which showed a more favorable result for TMB/*BRAF*), standard care probability dropped to 88·5% at the same threshold, with TMB/*BRAF* accounting for 7·7%. Cost-effectiveness planes (figure [X]) showed that in Models B and C, biomarker-guided strategies clustered in the northwest quadrant (higher cost, lower effect), while in Model A they appeared in the northeast quadrant with minimal incremental effectiveness. Uncertainty was consistent across model specifications, though Model A produced wider confidence ellipses reflecting the greater parametric uncertainty from estimating additional interaction terms.

Scenario analysis with biosimilar pricing (67% nivolumab cost reduction) is being conducted to determine whether substantial price reductions would alter cost-effectiveness conclusions. Given the finding that biomarker strategies in the overall population provide fewer QALYs than control under Models B and C, cost reductions alone may be insufficient to achieve cost-effectiveness unless clinical effectiveness also improves. Preliminary results suggest that [biosimilar scenario ICERs and optimal strategy determinations are being calculated and will be reported in final manuscript version].

### Expected value of partial perfect information

At the base-case WTP of €51,000/QALY, total per-patient EVPI was negligible (<€100) under all model specifications, consistent with the finding that no biomarker strategy approaches cost-effectiveness. Among individual parameter groups, health-state utilities showed the highest per-patient EVPPI (€888–1,064 for progression-free utility across model specifications), reflecting the sensitivity of QALY calculations to the utility decrement upon progression. However, even these values translate to minimal population-level research value at current immunotherapy prices.

Under the biosimilar pricing scenario (80% discount from list price; €4,641/administration), total EVPI increased dramatically to approximately €8,700–8,900 per patient, reflecting substantial decision uncertainty when immunotherapy becomes more affordable. However, EVPPI for individual parameter groups remained modest: health-state utilities contributed €650–950 per patient, while drug costs and other parameters each contributed <€200 per patient. Notably, treatment–biomarker interaction parameters were not included in this EVPPI analysis due to methodological challenges in extracting interaction coefficients from the resampled survival models; this represents an important limitation, as these parameters are theoretically the most policy-relevant for biomarker validation decisions.

Population-level EVPPI estimates would be calculated by scaling per-patient values to the eligible Norwegian MSS/pMMR mCRC population (estimated ~15,000 patients over a 10-year technology horizon). At current prices, population EVPPI is negligible (<€1·5 million for utilities, the highest-value parameter group). Under biosimilar pricing, population EVPPI for utilities would reach approximately €[10–15] million, with other parameter groups contributing smaller amounts. These estimates inform the maximum value of research to resolve specific parameter uncertainties and can guide decisions about prospective biomarker validation trials, though the absence of interaction parameter EVPPI limits their policy relevance. [Note: Population scaling calculations in progress; bracketed values are preliminary estimates pending confirmation of eligible population size.]

---

## Discussion

This early economic evaluation demonstrates that biomarker-guided chemo-immunotherapy in MSS/pMMR mCRC does not represent cost-effective use of healthcare resources at current prices, a conclusion that is robust across three model specifications, two target populations, and comprehensive uncertainty analysis. Beyond this applied finding, the analysis yields four methodological insights relevant to the growing field of precision oncology health economics.

First, the exploratory treatment–biomarker interaction signals identified in the companion clinical analysis did not translate into predicted QALY gains under parametric survival modeling. The companion paper reported CRP–treatment interaction hazard ratios of 0·19–0·22 for PFS—an apparently strongly protective signal.^12^ Yet the gamma parametric model predicts that CRP-positive patients fare 0·54 QALYs worse on immunotherapy than on chemotherapy alone in the enriched-population analysis. This reversal likely arises because the suggestive PFS interaction signal, which was not significant after false discovery rate correction, was insufficient to generate OS gains after parametric extrapolation over a ten-year horizon. The gamma distribution’s tail behavior during extrapolation can attenuate or reverse within-trial signals, particularly when the interaction effect operates primarily on PFS rather than OS. Moreover, the near-identical health-state utilities for progression-free (0·91) and post-progression (0·90) states mean that any PFS advantage translates into minimal QALY gains, with the model driven almost entirely by OS differences. This finding illustrates a broader principle: promising treatment–biomarker interaction signals from semi-parametric Cox analyses should not be assumed to translate into economic value without formal evaluation through parametric survival models with full extrapolation.

Second, the comparison of Model A (which includes TLR) and Model B (which excludes it on principled grounds) reveals that conditioning on a post-treatment prognostic variable altered clinical effectiveness estimates for the control strategy (1·84 vs 2·08 QALYs, a 13% difference) but had more modest effects on relative treatment effects. TLR is a dynamic biomarker measured after treatment initiation; including it as a covariate when estimating treatment effects risks blocking causal pathways through which immunotherapy operates, or inducing collider bias if early tumor response is a common effect of both treatment and unmeasured prognostic factors.^13,14^ This concern is well-established in the causal inference literature but has not been addressed in HTA guidance for survival model covariate selection, including NICE DSU Technical Support Documents.^15,18^ Our analysis provides reassurance that in settings with weak biomarker intercorrelations (Cramér's V = 0·278), the practical impact on relative treatment effects may be limited—though the principle remains critical for future evaluations with more strongly correlated dynamic biomarkers or larger treatment effects, and the substantial impact on absolute effectiveness estimates (13%) demonstrates that the choice remains consequential. The finding establishes that covariate selection in multi-biomarker CEAs deserves the same systematic scrutiny currently applied to parametric distributional assumptions.

Third, the TLR-guided strategy was strongly dominated under all model specifications in which it was evaluable, providing a concrete economic demonstration of why the prognostic–predictive distinction matters for precision oncology. A biomarker that identifies patients with favorable prognosis—as TLR does—does not identify patients who benefit from the experimental therapy. Selecting patients for immunotherapy based on a prognostic marker means treating patients who would have done well on any therapy, while the biomarker-negative patients—who have worse prognosis—receive the same standard chemotherapy they would have received regardless. This insight, while established in the biomarker methodology literature,^19,20^ has not previously been demonstrated in an economic evaluation. The finding should caution against the common practice of translating any subgroup signal from a negative trial into a biomarker-guided treatment strategy without first formally testing the treatment–biomarker interaction.

Fourth, the model specification comparison demonstrated that clinical effectiveness estimates are sensitive to model formula choices: control strategy QALYs varied by 12·7% (from 1·84 in Model A to 2·08 in Model B), while biomarker strategies showed smaller variation. The qualitative conclusion—that no biomarker strategy is cost-effective—was robust across specifications, though the specific patterns differed: Models B and C showed all biomarker strategies dominated (fewer QALYs than control), while Model A showed one strategy (TMB/*BRAF*) with minimal incremental effect but an astronomically high ICER (€1,690,480/QALY). This variation underscores that economic evaluations in multi-biomarker settings should routinely report sensitivity to model specification, analogous to the now-standard practice of reporting sensitivity to parametric distributional assumptions.

The value-of-information analysis provides a quantitative framework for research prioritisation that integrates clinical and economic uncertainty. At current nivolumab prices, the decision is not close: standard care is clearly preferred with >98% probability, and additional research has negligible expected value. However, the pharmaceutical landscape for checkpoint inhibitors is evolving rapidly, and our biosimilar scenario analysis shows that under biosimilar pricing (80% discount from list price; €4,641/administration), total EVPI increases to approximately €8,700–8,900 per patient, reflecting substantial decision uncertainty. However, since biomarker strategies yield fewer QALYs than control under most specifications, much of this EVPI reflects uncertainty about the magnitude of net harm rather than proximity to a close decision. EVPPI for treatment–biomarker interaction parameters—once available—will clarify whether resolving uncertainty about the predictive effects specifically could change the optimal treatment decision. This conditional finding is important: it means that the decision about whether to pursue biomarker validation is inextricable from the trajectory of immunotherapy pricing. Research funders and trialists should consider pharmaceutical market dynamics as an input to trial design decisions, not merely as background context.

The enriched-population analysis establishes definitive bounds on the economic value of biomarker-guided strategies. Under Model B (primary analysis), all biomarker-positive subgroups fared worse on experimental treatment than on standard chemotherapy, with effectiveness losses ranging from −0·27 to −0·54 QALYs—confirming that the dominated results in the overall population reflect the parametric model’s prediction of net treatment harm, not prevalence dilution from biomarker-negative patients. Under the most favorable model specification (Model C, TMB/*BRAF*-enriched), the ceiling ICER was €323,878 per QALY—sixfold above the threshold. This finding has critical implications: when even the best-case scenario (perfect biomarker targeting, no testing costs, most optimistic survival model) produces an ICER six times the threshold, no amount of improved patient selection, biomarker refinement, or testing strategy optimization can achieve cost-effectiveness. Only fundamental changes to treatment costs or clinical effectiveness can alter this conclusion.

Several limitations warrant consideration. This analysis relies on exploratory biomarker findings from post-hoc subgroup analyses of a single phase II trial (n=65 complete cases). The treatment–biomarker interaction effects, while biologically plausible, lack prospective validation, and the parametric survival extrapolations beyond observed follow-up introduce substantial uncertainty. The biomarker definitions were established from the same trial, inflating the risk of overfitting. Health-state utilities from METIMMOX EQ-5D-5L data showed a very small decrement upon progression (0·91 vs 0·90), substantially smaller than typically reported in mCRC literature (0·05–0·15 decrement).^27^ This near-identical utility across health states means the model is driven almost entirely by OS differences, with PFS gains contributing minimally to QALYs; scenario analysis using published utility values with a larger progression decrement should be considered. The analysis reflects Norwegian healthcare costs and WTP thresholds; while qualitative conclusions likely generalize to similar high-income settings, quantitative results are context-specific. Despite these limitations, early economic evaluation provides essential information for research prioritisation in areas with limited therapeutic options: demonstrating that a strategy is not viable under any plausible scenario is itself a valuable finding that can redirect scarce research resources.

In conclusion, biomarker-guided chemo-immunotherapy in MSS/pMMR mCRC is not cost-effective at current immunotherapy prices, regardless of model specification, target population, or candidate biomarker. Under the primary model specification, the parametric survival model predicts net harm from adding immunotherapy even in biomarker-positive patients, and no strategy produces positive incremental QALYs—meaning that price reduction alone is insufficient and prospective validation of genuine clinical benefit is the necessary precondition for cost-effectiveness. This analysis establishes that survival model specification—particularly the treatment of post-treatment prognostic variables, the choice of which biomarkers to include as covariates, and the interaction between semi-parametric interaction signals and parametric extrapolation—should be systematically evaluated in economic assessments of multi-biomarker precision oncology strategies.

---

## Contributors

[To be completed. Suggested structure: BPG and EA conceived the study. BPG developed the economic model and conducted the analyses. MAMV contributed to the model structure and sensitivity analyses. SM and PAB contributed clinical data and biomarker definitions. EAB contributed to the statistical methods. AHR was principal investigator of the METIMMOX trial and provided clinical oversight. EA supervised the health economic analysis. BPG wrote the first draft. All authors contributed to data interpretation, critically revised the manuscript, and approved the final version.]

## Data sharing

The complete analysis code is available at https://github.com/ben-geisler/METIMMOX-1 under MIT license. Individual patient data may be made available to qualified researchers upon reasonable request, subject to appropriate data-sharing agreements and ethical approvals.

## Declaration of interests

[To be completed by all authors]

## Acknowledgements

This work was supported by the ASCERTAIN project (Grant Agreement 101094938), funded by the European Union, and the Research Council of Norway as part of NORCHER (Research Grant 296114). Views and opinions expressed are those of the authors only and do not necessarily reflect those of the European Union or the European Health and Digital Executive Agency (HaDEA). [Additional funding sources to be confirmed]

---

## References

1 Sinicrope FA, Ou F-S, Arnold D, et al. Randomized trial of standard chemotherapy alone or combined with atezolizumab as adjuvant therapy for patients with stage III deficient DNA mismatch repair (dMMR) colon cancer (Alliance A021502; ATOMIC). *J Clin Oncol* 2025; **43**(17\_suppl): LBA1.

2 Chalabi M, Verschoor YL, Tan PB, et al. Neoadjuvant immunotherapy in locally advanced mismatch repair-deficient colon cancer. *N Engl J Med* 2024; **390**: 1949–58.

3 Cercek A, Foote MB, Rousseau B, et al. Nonoperative management of mismatch repair-deficient tumors. *N Engl J Med* 2025; **392**: 2297–308.

4 Andre T, Elez E, Van Cutsem E, et al. Nivolumab plus ipilimumab in microsatellite-instability-high metastatic colorectal cancer. *N Engl J Med* 2024; **391**: 2014–26.

5 Siegel RL, Kratzer TB, Giaquinto AN, Sung H, Jemal A. Cancer statistics, 2025. *CA Cancer J Clin* 2025; **75**: 10–45.

6 Yan S, Wang W, Feng Z, et al. Immune checkpoint inhibitors in colorectal cancer: limitation and challenges. *Front Immunol* 2024; **15**: 1403533.

7 Zeineddine FA, Zeineddine MA, Yousef A, et al. Survival improvement for patients with metastatic colorectal cancer over twenty years. *NPJ Precis Oncol* 2023; **7**: 16.

8 Ree AH, Šaltytė Benth J, Hamre HM, et al. First-line oxaliplatin-based chemotherapy and nivolumab for metastatic microsatellite-stable colorectal cancer—the randomized METIMMOX trial. *Br J Cancer* 2024; **130**: 1921–28.

9 Meltzer S, Berg JP, Hamre HM, et al. 632P Predictive value of C-reactive protein (CRP) in microsatellite-stable (MSS) metastatic colorectal cancer (mCRC) patients given first-line alternating short-course oxaliplatin-based chemotherapy (FLOX) and nivolumab. *Ann Oncol* 2023; **34**: S449.

10 Meltzer S, Negård A, Bakke KM, et al. Early radiologic signal of responsiveness to immune checkpoint blockade in microsatellite-stable/mismatch repair-proficient metastatic colorectal cancer. *Br J Cancer* 2022; **127**: 2227–33.

11 Ree AH, Bousquet PA, Nilsen HL, et al. 543P Tumor mutational burden (TMB), BRAF status, and C-reactive protein (CRP) predict response to first-line alternating oxaliplatin-based chemotherapy and nivolumab in metastatic microsatellite-stable (MSS) colorectal cancer (CRC). *Ann Oncol* 2024; **35**: S453.

12 Geisler BP, Meltzer S, Bousquet PA, et al. Exploratory biomarker analysis of treatment effect heterogeneity in MSS/pMMR metastatic colorectal cancer from the randomized METIMMOX trial. [Companion manuscript].

13 Hernán MA, Robins JM. *Causal Inference: What If.* Boca Raton: Chapman & Hall/CRC, 2020.

14 Daniel RM, Kenward MG, Cousens SN, De Stavola BL. Using causal diagrams to guide analysis in missing data problems. *Stat Methods Med Res* 2012; **21**: 243–56.

15 Latimer NR, Abrams KR, Lambert PC, et al. Adjusting survival time estimates to account for treatment switching in randomized controlled trials—an economic evaluation context: methods, limitations, and recommendations. *Med Decis Making* 2014; **34**: 387–402.

16 Norwegian Medicines Agency. Guidelines for the submission of documentation for single technology assessment (STA) of pharmaceuticals. Oslo: NoMA, 2023.

17 Jensen CE, Sørensen SS, Gudex C, Jensen MB, Pedersen KM, Ehlers LH. The Danish EQ-5D-5L value set: a hybrid model using cTTO and DCE data. *Appl Health Econ Health Policy* 2021; **19**: 579–91.

18 Latimer NR. Survival analysis for economic evaluations alongside clinical trials—extrapolation with patient-level data: inconsistencies, limitations, and a practical guide. *Med Decis Making* 2013; **33**: 743–54.

19 Ballman KV. Biomarker: predictive or prognostic? *J Clin Oncol* 2015; **33**: 3968–71.

20 Sargent DJ, Conley BA, Allegra C, Collette L. Clinical trial designs for predictive marker validation in cancer treatment trials. *J Clin Oncol* 2005; **23**: 2020–27.

21 Clark GM. Prognostic factors versus predictive factors: examples from a clinical trial of erlotinib. *Mol Oncol* 2008; **1**: 406–12.

22 Proctor MJ, Morrison DS, Talwar D, et al. A comparison of inflammation-based prognostic scores in patients with cancer: a Glasgow Inflammation Outcome Study. *Eur J Cancer* 2011; **47**: 2633–41.

23 Sharma P, Hu-Lieskovan S, Wargo JA, Ribas A. Primary, adaptive, and acquired resistance to cancer immunotherapy. *Cell* 2017; **168**: 707–23.

24 Fenwick E, Claxton K, Sculpher M. Representing uncertainty: the role of cost-effectiveness acceptability curves. *Health Econ* 2001; **10**: 779–87.

25 Claxton K, Palmer S, Longworth L, et al. Informing a decision framework for when NICE should recommend the use of health technologies only in the context of an appropriately designed program of evidence development. *Health Technol Assess* 2012; **16**: 1–323.

26 Briggs AH, Weinstein MC, Fenwick EAL, Karnon J, Sculpher MJ, Paltiel AD. Model parameter estimation and uncertainty analysis: a report of the ISPOR-SMDM Modeling Good Research Practices Task Force Working Group–6. *Med Decis Making* 2012; **32**: 722–32.

27 Joranger P, Nesbakken A, Sorbye H, Hoff G, Oshaug A, Aas E. Survival and costs of colorectal cancer treatment and effects of changing treatment strategies: a model approach. *Eur J Health Econ* 2020; **21**: 321–34.

28 Bjørnelv GMW, Dueland S, Line PD, et al. Cost-effectiveness of liver transplantation in patients with colorectal metastases confined to the liver. *Br J Surg* 2019; **106**: 132–41.

29 Eng C, Kim TW, Bendell J, et al. Atezolizumab with or without cobimetinib versus regorafenib in previously treated metastatic colorectal cancer (IMblaze370): a multicentre, open-label, phase 3, randomized controlled trial. *Lancet Oncol* 2019; **20**: 849–61.

30 Mettu NB, Ou F-S, Engel Nitz NM, et al. Assessment of capecitabine and bevacizumab with or without atezolizumab for the treatment of refractory metastatic colorectal cancer: a randomized clinical trial. *JAMA Netw Open* 2022; **5**: e2149040.

31 Chen EX, Jonker DJ, Loree JM, et al. Effect of combined immune checkpoint inhibition vs best supportive care alone in patients with advanced colorectal cancer: the Canadian Cancer Trials Group CO.26 study. *JAMA Oncol* 2020; **6**: 831–38.

32 Firth D. Bias reduction of maximum likelihood estimates. *Biometrika* 1993; **80**: 27–38.

33 Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 1995; **57**: 289–300.

---

## Tables

### Table 1: Model input parameters

| Parameter | Base case | Range | Distribution | Source | Notes |
|---|---|---|---|---|---|
| **Biomarker prevalences** | | | | | |
| CRP-positive | 0·338 | 0·23–0·47 | Beta | METIMMOX | 95% CI: 22·6–46·6% |
| TLR-positive | 0·631 | 0·50–0·75 | Beta | METIMMOX | Model A/C only; 95% CI: 50·2–74·7% |
| TMB/*BRAF*-positive | 0·446 | 0·32–0·58 | Beta | METIMMOX | 95% CI: 32·3–57·5% |
| **Health-state utilities** | | | | | |
| Progression-free | 0·91 | 0·73–1·00 | Beta | METIMMOX EQ-5D | Danish tariff |
| Post-progression | 0·90 | 0·72–1·00 | Beta | METIMMOX EQ-5D | Danish tariff |
| **Treatment costs (€)** | | | | | |
| Nivolumab per cycle | 13,923 | 11,138–16,708 | Gamma | NoMA | 40% discount |
| FLOX per cycle | 427 | 342–512 | Gamma | NoMA | |
| NGS testing | 2,518 | 2,014–3,022 | Gamma | Local | TMB/*BRAF* only |
| **Other costs (€)** | | | | | |
| Monitoring per cycle | 386 | 309–463 | Gamma | Ref 27 | CT scans |
| Post-progression /wk | — | — | — | — | Included in end-of-life |
| End-of-life | 13,803 | 11,042–16,564 | Gamma | Ref 28 | |
| **Discount rate** | 4% | 0–8% | Fixed | NoMA guidelines | |
| **Time horizon** | 10 years | 5–20 yr | Fixed | | |

NoMA = Norwegian Medicines Agency. NGS = next-generation sequencing.

### Table 4: Base-case cost-effectiveness results: overall MSS population (Model B)

| Strategy | Cost (€) | LYs | QALYs | Δ Cost (€) | Δ QALYs | ICER |
|---|---|---|---|---|---|---|
| Standard care | 21,504 | 2·46 | 2·08 | — | — | — |
| CRP-guided | 55,118 | 2·23 | 1·90 | 33,614 | −0·18 | Dominated |
| TMB/*BRAF*-guided | 61,542 | 2·29 | 1·94 | 40,038 | −0·14 | Dominated |
| TLR-guided* | 86,411 | 1·96 | 1·70 | 64,907 | −0·38 | Dominated |

Costs and QALYs discounted at 4% per year. LYs = life-years. QALYs = quality-adjusted life-years. ICER = incremental cost-effectiveness ratio (€ per QALY gained) versus next non-dominated alternative. Dom. = dominated. *TLR-guided evaluated under Models A and C only.

### Table 5: ICERs for biomarker-guided strategies across model specifications (€ per QALY, overall MSS population)

| Strategy | Model A | Model B | Model C |
|---|---|---|---|
| CRP-guided vs SoC | Dominated | Dominated | Dominated |
| TMB/*BRAF*-guided vs SoC | 1,690,480 | Dominated | Dominated |
| TLR-guided vs SoC | Dominated | N/A† | Dominated |

SoC = Standard of care. †TLR excluded from Model B by design. All biomarker strategies yield fewer QALYs than standard care in Models B and C. In Model A, only TMB/*BRAF* yields a positive (though very high) ICER; CRP and TLR are dominated.

Model A: comprehensive (all biomarkers, all interactions). Model B: focused (CRP and TMB/*BRAF* as covariates and interactions, TLR excluded). Model C: single biomarker.

### Table 6: Biomarker-enriched population results (Model B)

| Biomarker subgroup | Treatment | Cost (€) | QALYs | Δ Cost (€) | Δ QALYs | ICER |
|---|---|---|---|---|---|---|
| CRP-positive | FLOX | 22,347 | 2·87 | — | — | — |
| CRP-positive | FLOX–nivo | 121,746 | 2·34 | 99,399 | −0·54 | Dominated |
| TMB/*BRAF*-positive | FLOX | 24,265 | 2·19 | — | — | — |
| TMB/*BRAF*-positive | FLOX–nivo | 107,775 | 1·89 | 83,510 | −0·31 | Dominated |
| TLR-positive | FLOX | 22,870 | 2·35 | — | — | — |
| TLR-positive | FLOX–nivo | 122,839 | 2·08 | 99,969 | −0·27 | Dominated |

Biomarker-enriched analysis restricts cohort to biomarker-positive patients only. Costs and QALYs discounted at 4%. Dominated = more costly and less effective. Under Model C, TMB/*BRAF*-enriched yielded ICER €323,878/QALY; TLR-enriched €1,242,396/QALY. Under Model A, TMB/*BRAF*-enriched yielded ICER €1,625,723/QALY.

---

## Figure legends

**Figure 1:** Overall survival and progression-free survival curves for the CRP-guided strategy versus standard of care across model specifications (A, B, C). Shaded areas represent 95% confidence intervals from parametric bootstrap. Kaplan–Meier curves from observed data are overlaid for visual assessment of fit.

**Figure 2:** Tornado diagram showing one-way deterministic sensitivity analysis for the CRP-guided strategy versus standard of care (Model B). Bars represent variation in incremental net monetary benefit when each parameter is varied across its plausible range. Vertical dashed line indicates the base-case incremental net monetary benefit at €51,000/QALY.

**Figure 3:** Cost-effectiveness plane showing joint uncertainty from probabilistic sensitivity analysis (5,000 iterations) for all strategies versus standard care (Model B). Dashed line represents the willingness-to-pay threshold of €51,000/QALY.

**Figure 4:** Cost-effectiveness acceptability curves showing the probability that each strategy is cost-effective across a range of willingness-to-pay thresholds. Panel A: Model B. Panel B: Model A. Panel C: Biosimilar pricing scenario (Model B).

**Figure 5:** Expected value of partial perfect information (EVPPI) by parameter group under three scenarios: base case (€51,000/QALY), biosimilar pricing (€51,000/QALY), and elevated threshold (€100,000/QALY). Parameter groups: baseline survival, biomarker prognostic effects, drug costs, and utilities. [Treatment–biomarker interaction EVPPI to be added when available.]

---

## Supplementary appendix (outline)

- Table S1: Baseline patient characteristics by treatment arm
- Table S2: Biomarker intercorrelations (Cramér's V)
- Table S3: Parametric distribution selection (AIC, BIC, visual fit)
- Table S4: Full parametric model coefficients for all specifications (A, B, C)
- Table S5: Full cost-effectiveness results for Models A and C (overall MSS population)
- Table S6: Full cost-effectiveness results for biomarker-enriched populations
- Table S7: EVPPI results by parameter group, scenario, and model specification
- Table S8: Scenario analyses (biosimilar pricing, alternative WTP thresholds, alternative discount rates, literature-based utility values)
- Figure S1: Kaplan–Meier curves versus parametric fits by biomarker subgroup
- Figure S2: Survival curve extrapolations under alternative distributional assumptions
- Figure S3: Cost-effectiveness planes for all model specifications
- Figure S4: CEACs for all model specifications
- Figure S5: Cost-effectiveness frontier plots for each model specification
- Analysis code: Available at https://github.com/ben-geisler/METIMMOX-1
