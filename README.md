# Evolution of the Incremental Preventive Effect of Fluoride Varnish

### Reproducible R/Stan code and Supplementary Material

This repository contains the reproducible statistical analysis and technical Supplement for the manuscript:

> **Evolution of the Incremental Preventive Effect of Fluoride Varnish: A Log-Ratio Meta-Analysis of Historical and Contemporary Trials**  
> Sara Antonijevic, Danielle Sitalo, and Brani Vidakovic  
> Department of Statistics, Texas A&M University

The repository is organized around one scientific question: **has the incremental preventive effect attributed to fluoride varnish changed as the background preventive environment has changed over time?** The analysis separates a historical nine-trial reference set from a contemporary extension and deliberately distinguishes the average historical effect from prediction in a new setting, historical-completeness sensitivity, temporal change, and explanatory study-level moderators.

The preferred continuous-outcome effect measure is the log ratio of treatment and control means,

$$
y_i = \log\left(\frac{\bar X_{Ti}}{\bar X_{Ci}}\right),
$$

with results translated to the clinically interpretable positive prevented-fraction scale

$$
\mathrm{PF}_{+,i}=1-\exp(y_i).
$$

Thus negative log ratios favor fluoride varnish, while positive values of $\mathrm{PF}_+$ represent prevention.

---

## Repository layout

```text
FluorideVarnish_GitHub/
├── README.md
├── .gitignore
├── MANIFEST_SHA256.txt
├── SUPPLEMENT_fluoridevarnish.pdf
    ├── 00_install_packages.R
    ├── 01_historical_classical_REML_HK.R
    ├── 02_historical_bayesmeta_prior_sensitivity.R
    ├── 03_historical_armlevel_Stan.R
    ├── historical_armlevel.stan
    ├── 04_Cochrane13_completeness.R
    ├── 05_temporal_meta_regression.R
    ├── calendar_time_meta_regression.stan
    ├── 06_explanatory_meta_regression.R
    ├── time_background_prevention.stan
    ├── time_background_burden.stan
    ├── 07_binary_outcomes_supporting_TEMPLATE.R
    └── data/
        ├── contemporary_continuous_variance_audit.csv
        └── binary_event_data_TEMPLATE.csv
```

The numbered R programs are intended to be read and run in order. Each program creates its own `output_*` directory under `code/`; generated output directories are ignored by Git because they can be regenerated from the source files.

---

## Analysis map

| Analysis | Main program | Stan model(s) | Role |
|---|---|---|---|
| Historical classical benchmark | `01_historical_classical_REML_HK.R` | — | REML, Hartung–Knapp inference, prediction, influence |
| Historical Bayesian NNHM | `02_historical_bayesmeta_prior_sensitivity.R` | — | Deterministic Bayesian normal–normal meta-analysis and prior sensitivity |
| Historical aggregate arm-level models | `03_historical_armlevel_Stan.R` | `historical_armlevel.stan` | Gaussian and standardized $t_4$ random effects, MCMC diagnostics, posterior prediction |
| Historical completeness | `04_Cochrane13_completeness.R` | — | Nine-study versus Cochrane-13 sensitivity on the legacy PF scale |
| Calendar-time analysis | `05_temporal_meta_regression.R` | `calendar_time_meta_regression.stan` | Final 14-study historical-to-contemporary meta-regression and sensitivity analysis |
| Explanatory moderators | `06_explanatory_meta_regression.R` | `time_background_prevention.stan`, `time_background_burden.stan` | Organized prevention and control-burden explanations of the calendar-time signal |
| Contemporary binary outcomes | `07_binary_outcomes_supporting_TEMPLATE.R` | — | Supporting risk-ratio analysis; template pending frozen verified binary inputs |

---

## Software requirements

The classical analyses require R and [`metafor`](https://cran.r-project.org/package=metafor). The deterministic Bayesian normal–normal analysis uses [`bayesmeta`](https://cran.r-project.org/package=bayesmeta). The MCMC analyses use [`cmdstanr`](https://mc-stan.org/cmdstanr/), [`posterior`](https://mc-stan.org/posterior/), CmdStan, and a working C++ toolchain.

The production Stan analyses reported in the Supplement used:

- R 4.5.3
- CmdStanR 0.9.0
- `posterior` 1.7.0
- CmdStan 2.40.0

Later compatible versions may work, but reproducible runs should archive `sessionInfo()`, CmdStan version information, random-number seeds, and the generated CSV summaries.

For package setup:

```r
setwd("code")
source("00_install_packages.R")
```

CmdStan itself must also be installed for the Stan analyses. See the CmdStanR documentation for platform-specific toolchain setup.

---

## Recommended run order

From R, set the working directory to `code/` and run:

```r
source("01_historical_classical_REML_HK.R")
source("02_historical_bayesmeta_prior_sensitivity.R")
source("03_historical_armlevel_Stan.R")
source("04_Cochrane13_completeness.R")
source("05_temporal_meta_regression.R")
source("06_explanatory_meta_regression.R")
```

Run the binary-outcome analysis only after the verified binary-event/design-adjusted input file has been frozen:

```r
source("07_binary_outcomes_supporting_TEMPLATE.R")
```

The calendar-time program has a `RUN_STAN` switch near the top. Set it to `FALSE` to run only the classical/meta-regression and sensitivity components when CmdStan is unavailable; set it to `TRUE` for the Bayesian calendar-time fit reported in the manuscript.

---

# Detailed program documentation

## `00_install_packages.R` — package setup

A small helper for installing/loading the main R packages required by the analysis. It is not part of the statistical analysis and does not generate manuscript results. CmdStan installation remains a separate system-level step because it depends on the local C++ toolchain.

---

## `01_historical_classical_REML_HK.R` — classical historical benchmark

This program reconstructs the historical nine-study analysis directly from frozen arm-level means, standard deviations, and sample sizes. It does **not** enter rounded log ratios copied from the manuscript.

### Statistical tasks

- constructs the first-order log ratio of means and its delta-method sampling variance;
- fits an intercept-only random-effects model by REML;
- uses Hartung–Knapp small-sample inference for the pooled mean;
- reports the Riley-style $t_{k-2}$ prediction interval for a new true study effect;
- transforms estimates and intervals from the log-ratio scale to the positive prevented-fraction scale;
- reports heterogeneity summaries including $\widehat\tau$, $Q$, $I^2$, and $H^2$;
- performs leave-one-study-out influence analysis;
- generates the historical effect plot and influence diagnostics.

The script explicitly protects the **primary first-order LRR definition** from a software-default second-order correction. The Supplement separately reports the targeted second-order sensitivity calculation.

### Expected benchmark

Approximately:

- pooled log ratio: `-0.5050`
- $\widehat\tau$: `0.2835`
- pooled $\mathrm{PF}_+$: `0.3965`
- Hartung–Knapp 95% CI for $\mathrm{PF}_+$: `[0.180, 0.556]`
- 95% prediction interval for $\mathrm{PF}_{+,\mathrm{new}}$: `[-0.265, 0.712]`

### Output directory

`code/output_01_historical_classical/`

Principal generated files include:

- `historical9_study_effects.csv`
- `historical9_REML_HK_summary.csv`
- `historical9_leave_one_out.csv`
- `historical9_REML_HK_forest_PF.pdf`
- `historical9_influence_diagnostics.pdf`
- `R_sessionInfo.txt`

---

## `02_historical_bayesmeta_prior_sensitivity.R` — Bayesian NNHM and heterogeneity-prior sensitivity

This program analyzes the same nine historical log ratios under the Bayesian normal–normal hierarchical model

$$
y_i\mid\theta_i\sim N(\theta_i,s_i^2),
\qquad
\theta_i\mid\mu,\tau\sim N(\mu,\tau^2),
$$

with a broad $\mu\sim N(0,1^2)$ prior.

### Heterogeneity priors

The script compares:

1. Half-Normal$(0,0.5^2)$ — primary prior;
2. Half-Normal$(0,0.25^2)$;
3. Half-$t_4$ with scale 0.5;
4. DuMouchel prior.

The calculation uses deterministic numerical integration through `bayesmeta`; it is **not** an MCMC analysis. Therefore $\widehat R$, ESS, divergences, and tree-depth diagnostics are not relevant to this program.

### Statistical outputs

- posterior pooled prevented fraction;
- posterior heterogeneity $\tau$;
- posterior predictive distribution for a new true study effect;
- $\Pr(\mathrm{PF}_{+,\mathrm{new}}>0\mid\mathcal D)$;
- $\Pr(\mathrm{PF}_{+,\mathrm{new}}>0.20\mid\mathcal D)$;
- study-specific posterior shrinkage;
- prior and posterior heterogeneity-density plots.

### Primary-prior benchmark

Approximately:

- pooled median $\mathrm{PF}_+$: `0.389`
- 95% CrI: `[0.215, 0.556]`
- median $\tau$: `0.302`
- predictive median $\mathrm{PF}_{+,\mathrm{new}}$: `0.386`
- 95% predictive interval: about `[-0.284, 0.732]`

### Output directory

`code/output_02_historical_bayesmeta/`

Principal outputs include the full/paper prior-sensitivity tables, prior calibration, study-specific shrinkage, posterior density plots, predictive density plots, a LaTeX table, and `R_sessionInfo.txt`.

---

## `03_historical_armlevel_Stan.R` + `historical_armlevel.stan` — aggregate arm-level Bayesian models

These files provide a robustness analysis that avoids first constructing a log ratio and plug-in sampling variance. Treatment and control sample means are modeled directly at the arm level.

Two random-effects specifications are fit from the same Stan program:

- Gaussian random effects;
- standardized Student-$t_4$ random effects.

The $t_4$ distribution is standardized to variance one, so $\tau$ remains interpretable as a between-study standard deviation under both models.

### Production MCMC settings

- 4 chains;
- 2000 warm-up iterations per chain;
- 4000 retained iterations per chain;
- `adapt_delta = 0.995`;
- `max_treedepth = 15`.

### What is checked

- rank-normalized split $\widehat R$;
- bulk and tail ESS;
- Monte Carlo standard errors;
- divergent transitions;
- maximum-tree-depth hits;
- E-BFMI;
- trace plots;
- posterior predictive arm-level checks;
- study-specific posterior effects;
- pooled and predictive prevented-fraction distributions.

### Expected pooled posterior medians

- Gaussian random effects: $\mathrm{PF}_+\approx0.433$;
- standardized $t_4$ random effects: $\mathrm{PF}_+\approx0.421$.

The main substantive distinction is in heterogeneity and prediction rather than in the historical pooled center.

### Output directory

`code/output_03_historical_Stan/`

The script writes model summaries, sampler diagnostics, CmdStan diagnostic reports, full parameter summaries, study-specific posterior effects, posterior predictive checks, trace plots, manuscript-ready LaTeX fragments, run configuration, and `R_sessionInfo.txt`.

Compiled Stan executables are deliberately not tracked in Git; they should be generated locally by CmdStanR from `historical_armlevel.stan`.

---

## `04_Cochrane13_completeness.R` — historical-completeness sensitivity

This analysis asks a narrower question than the primary historical model: **does the historical pooled center materially change when the original nine-study subset is expanded to the complete 13-trial permanent-surface evidence set represented in the 2013 Cochrane synthesis?**

It is intentionally carried out on the common legacy prevented-fraction scale rather than forcing four heterogeneous historical studies into the preferred modern log-ratio likelihood.

The additional studies are:

- Holm (1984);
- Borutta (1991);
- Hardman (2007);
- Gugwad (2011).

Their study-specific design/reconstruction issues are documented in the Supplement: shared controls, clustering, historical variance reconstruction, and a reconstructed PF value that is not representable as $1-\exp(\theta)$ for positive arm means.

### Methods

- REML random-effects pooling;
- ordinary and modified Hartung–Knapp inference;
- small-sample prediction intervals;
- optional Bayesian completeness comparison through `bayesmeta` when available.

### Expected benchmark

Approximately:

- nine-study Cochrane-scale pooled PF: `0.439`;
- Cochrane-13 pooled PF: `0.435`;
- $\tau_9$: `0.217`;
- $\tau_{13}$: `0.197`.

The central result is that completing the historical evidence set changes the pooled center by less than one percentage point.

### Output directory

`code/output_04_Cochrane13/`

Outputs include the reconstructed data, nine-versus-thirteen summary table, frequentist fits, optional Bayesian summary, forest/summary graphics, and session information.

---

## `05_temporal_meta_regression.R` + `calendar_time_meta_regression.stan` — final calendar-time production analysis

This is the principal historical-to-contemporary analysis. The primary continuous-outcome dataset contains 14 studies:

- the nine historical trials;
- Muñoz-Millán (2018);
- Latifi-Xhemajli (2019);
- McMahon (2020);
- Wang (2022);
- Zeng (2025).

Calendar time is centered at 2000 and scaled in decades:

$$
t_i=\frac{\mathrm{Year}_i-2000}{10}.
$$

The primary random-effects meta-regression is

$$
y_i=\beta_0+\beta_1t_i+u_i+\epsilon_i,
\qquad
u_i\sim N(0,\tau^2),
\qquad
\epsilon_i\sim N(0,s_i^2).
$$

A positive $\beta_1$ means that later treatment/control ratios move toward one; on the prevented-fraction scale this corresponds to attenuation of the **incremental** varnish effect.

### What the program does

- REML–Hartung–Knapp linear calendar-time meta-regression;
- historical-versus-contemporary era comparison;
- separate pooling of the five primary contemporary studies;
- maximum-likelihood AICc comparison of:
  - no calendar moderator,
  - linear calendar time,
  - historical-versus-contemporary era,
  - quadratic time,
  - three eras;
- leave-one-contemporary-study-out analysis;
- Wang cluster-size sensitivity;
- Jayasinghe quasi-Poisson variance sensitivity;
- Ghasemi ICC/design-effect sensitivity;
- optional Bayesian calendar-time Stan model;
- production tables and figures.

### Wang variance treatment

Wang is a cluster-randomized comparison. The production analysis uses the published ICC `0.20` and planning cluster size `m = 40`, giving

$$
\mathrm{DE}=1+(40-1)(0.20)=8.8.
$$

The ordinary individual-level LRR variance is multiplied by this design effect. This is described as an **approximate cluster adjustment**, because the article does not report a directly cluster-robust LRR standard error. Replacing the planning cluster size by the observed mean completed cluster size changes the fitted temporal slope only minimally.

### Jayasinghe and Ghasemi

These studies are sensitivity-only because a unique design-adjusted sampling variance cannot be recovered from the published continuous-outcome summaries:

- Jayasinghe: quasi-Poisson overdispersion grid;
- Ghasemi: ICC/design-effect grid for the four-school cluster trial.

### Expected primary result

Frequentist 14-study fit:

- $\widehat\beta_1\approx0.0986$ per decade;
- Hartung–Knapp 95% CI approximately `[-0.0445, 0.2417]`;
- $\widehat\tau\approx0.279$.

Bayesian calendar-time fit:

- posterior median $\beta_1\approx0.091$;
- 95% CrI approximately `[-0.036, 0.231]`;
- $\Pr(\beta_1>0\mid\mathcal D)\approx0.924$.

The AICc comparison is intentionally part of the interpretation: the intercept-only model has the smallest AICc, while linear time and the era model remain close. The temporal evidence is therefore treated as **suggestive/compatible with attenuation**, not as proof that calendar time is required by the data.

### Output directory

`code/output_05_temporal_regression/`

Principal outputs include scenario datasets, calendar-time summaries, era comparisons, contemporary pooling, AICc comparison, leave-one-modern-out results, Jayasinghe/Ghasemi variance grids, the variance audit, Bayesian results when enabled, and manuscript-ready text/tables/figures.

The human-readable contemporary variance audit is also tracked as:

`code/data/contemporary_continuous_variance_audit.csv`.

---

## `06_explanatory_meta_regression.R` + Stan models — what may explain attenuation?

This analysis uses the **same 14-study primary continuous-outcome dataset** as the calendar-time analysis and asks whether measured features of the preventive environment reduce the estimated calendar-time coefficient.

The organized-background-prevention indicator is set to one for:

- Tagliaferro;
- Muñoz-Millán;
- McMahon;
- Zeng.

The coding is deliberately narrow: routine community exposure to fluoridated toothpaste or water does not by itself constitute a protocol-specified organized control-arm intervention.

### Model sequence

1. time only;
2. time + organized background prevention;
3. time + organized prevention + control-arm burden.

The burden model is specifically designed to avoid naively regressing the observed response ratio on the observed control mean, because the control mean appears in the denominator of the response and is measured with error. `time_background_burden.stan` therefore uses an error-aware log-arm formulation that propagates uncertainty in the shared control-arm quantity.

### Expected posterior medians

Approximately:

- time only: $\beta_1=0.091$;
- time + prevention: $\beta_1=0.057$;
- time + prevention + burden: $\beta_1=0.053$.

The change from `0.091` to `0.057` is the main explanatory sensitivity result: measured organized prevention accounts for a non-negligible part of the fitted calendar-time signal. The analysis remains a small, observational, study-level meta-regression and **must not be interpreted as causal mediation**.

### Output directory

`code/output_06_explanatory_regression/`

Outputs include the frozen 14-study explanatory input table, posterior summaries for the three models, model diagnostic reports, and session information.

---

## `07_binary_outcomes_supporting_TEMPLATE.R` — supporting binary evidence

Several contemporary trials report child-level caries incidence rather than continuous DMF-type outcomes. These define a different estimand and are kept as a parallel evidence stream using the log risk ratio

$$
r_i=\log\left(\frac{p_{Ti}}{p_{Ci}}\right).
$$

The manuscript/Supplement describes five contemporary comparisons with crude risk ratios approximately:

- Muñoz-Millán: `0.81`;
- McMahon: `0.85`;
- Wang: `0.72`;
- Zeng: `1.08`;
- He: `0.81`.

The reported exploratory synthesis is approximately:

- pooled RR: `0.80`;
- $\widehat\tau_R$: `0.073`;
- approximate 95% interval: `[0.71, 0.89]`.

### Important reproducibility status

This file is intentionally marked **TEMPLATE**. The cleaned source archive did not contain the final frozen binary event-count table or verified design-adjusted log-risk-ratio standard errors for all studies. Those values should not be reconstructed by guesswork.

Before using this program, fill

`code/data/binary_event_data_TEMPLATE.csv`

with verified event counts or, where appropriate, directly reported design-adjusted logRR/SE values. For clustered trials, a valid design-adjusted effect and standard error are preferred to naive child-level binomial variances.

The binary analysis is supportive only and is not pooled with the continuous-outcome LRR analyses.

---

# Supplementary Material

The technical Supplement is in:

`/supplement/supplement01.tex`

It is a self-contained LaTeX document accompanying the manuscript. Its purpose is to provide the audit trail and methodological detail that would unnecessarily interrupt the main article: frozen inputs, exact effect construction, small-sample inference, prior calibration, MCMC diagnostics, prediction, study influence, historical reconstruction, modern variance auditing, temporal sensitivity, and explanatory models.

The repository copy has been adjusted so that its figure paths point to the output directories generated by the cleaned code. In particular, the Bayesian and Stan figures are expected under `code/output_02_historical_bayesmeta/` and `code/output_03_historical_Stan/`.

## Supplement sections and corresponding code

| Supplement section | Content | Reproducible source |
|---|---|---|
| **S1. Frequentist Benchmark for the Nine Historical Trials** | Frozen historical arm data; LRR construction; REML; Hartung–Knapp; Riley prediction; second-order LRR sensitivity; leave-one-out influence; funnel-plot interpretation | `01_historical_classical_REML_HK.R` for the primary benchmark/influence; second-order formula and numerical sensitivity are documented in the Supplement |
| **S2. Historical Completeness Sensitivity Analysis: the Cochrane-13 Evidence Set** | Cochrane-scale study data; Holm/Borutta/Hardman/Gugwad design issues; 9-vs-13 comparison; relation to the published synthesis | `04_Cochrane13_completeness.R` |
| **S3. Bayesian Normal–Normal Hierarchy and Prior Sensitivity** | Prior calibration; four heterogeneity priors; deterministic posterior computation; prediction; study-specific shrinkage | `02_historical_bayesmeta_prior_sensitivity.R` |
| **S4. Stan Production Analysis of the Aggregate Arm-Level Models** | Arm-level likelihood; Gaussian vs standardized $t_4$ random effects; NUTS diagnostics; PPCs; study-specific posterior effects | `03_historical_armlevel_Stan.R`, `historical_armlevel.stan` |
| **S5. Contemporary Continuous Outcomes and Temporal-Model Sensitivity** | Contemporary variance audit; Wang cluster adjustment; alternative temporal forms; AICc; primary 14-study fit; influence and Jayasinghe/Ghasemi sensitivities | `05_temporal_meta_regression.R`, `calendar_time_meta_regression.stan`, audit CSV |
| **S6. Contemporary Binary Outcomes as Supporting Evidence** | Parallel risk-ratio evidence stream | `07_binary_outcomes_supporting_TEMPLATE.R` after verified inputs are frozen |
| **S7. Explanatory Meta-Regression for the Calendar-Time Signal** | Organized prevention coding; time/prevention/burden model sequence; shared-denominator correction | `06_explanatory_meta_regression.R`, `time_background_prevention.stan`, `time_background_burden.stan` |
| **S8. Reproducibility and Final Analysis Record** | Software versions, seeds, frozen-data philosophy, relationship among evidence streams | Entire repository |

### Compiling the Supplement

The Supplement includes graphics produced by the Bayesian and Stan programs. A typical workflow is therefore:

```r
setwd("code")
source("02_historical_bayesmeta_prior_sensitivity.R")
source("03_historical_armlevel_Stan.R")
```

Then compile from the `supplement/` directory, for example with:

```bash
pdflatex supplement01.tex
pdflatex supplement01.tex
```

or with `latexmk` if available:

```bash
latexmk -pdf supplement01.tex
```

The LaTeX document uses standard packages including `geometry`, `amsmath`, `amssymb`, `booktabs`, `graphicx`, `listings`, `microtype`, and `hyperref`.

---

# Interpretation and evidence-stream separation

The repository intentionally keeps several analyses separate because they answer different statistical questions.

### 1. Historical nine-study LRR analysis

This is the preferred multiplicative continuous-outcome analysis and supplies the historical reference effect.

### 2. Cochrane-13 historical completeness

This is a legacy PF-scale sensitivity analysis. Its purpose is to show that the historical center is not an artifact of omitting four permanent-surface trials. It is not a substitute for the primary LRR likelihood.

### 3. Fourteen-study temporal analysis

This combines the historical nine with five primary contemporary continuous-outcome comparisons. It asks whether the treatment/control ratio has shifted over calendar time. Calendar time is interpreted as a marker of secular change, not as a causal mechanism.

### 4. Contemporary binary outcomes

These use risk ratios rather than ratios of continuous mean outcomes. They provide parallel supporting evidence and are not pooled with the continuous analysis.

### 5. Explanatory moderator models

These examine whether organized preventive background care and control-arm burden account for part of the temporal signal. They are observational study-level regressions, not causal mediation analyses.

---

# Core numerical benchmarks

These values are useful for checking a fresh run against the manuscript/Supplement.

| Analysis | Benchmark |
|---|---:|
| Historical REML pooled $\mathrm{PF}_+$ | 0.3965 |
| Historical HK 95% CI | [0.180, 0.556] |
| Historical $\widehat\tau$ | 0.2835 |
| Historical 95% prediction interval | [-0.265, 0.712] |
| Bayesian NNHM pooled median $\mathrm{PF}_+$ | 0.389 |
| Gaussian arm-level pooled median $\mathrm{PF}_+$ | 0.433 |
| Standardized $t_4$ arm-level pooled median $\mathrm{PF}_+$ | 0.421 |
| Cochrane-scale 9-study pooled PF | 0.439 |
| Cochrane-13 pooled PF | 0.435 |
| 14-study REML–HK calendar slope / decade | 0.0986 |
| 14-study HK interval for calendar slope | [-0.0445, 0.2417] |
| Bayesian calendar slope median | 0.091 |
| $\Pr(\beta_1>0\mid\mathcal D)$ | 0.924 |
| Time + organized prevention slope median | 0.057 |
| Time + prevention + burden slope median | 0.053 |

Small differences in last digits may occur with software-version changes or Monte Carlo variation, but material discrepancies should trigger a data/parameterization audit rather than being silently accepted.

---

# Reproducibility conventions

This repository follows several rules intended to make the analysis auditable:

- study-level effects are reconstructed from frozen input summaries whenever possible rather than copied from rounded manuscript tables;
- each evidence stream has a separate data object and output directory;
- the Cochrane-13 special historical reconstructions are not silently inserted into the primary LRR model;
- Wang is explicitly labeled as approximately cluster-adjusted;
- Jayasinghe and Ghasemi remain variance-sensitivity studies rather than being assigned invented standard errors;
- the binary analysis remains a template until its verified design-adjusted input file is frozen;
- Stan source is tracked, while compiled executables are not;
- scripts write `sessionInfo()` and/or sampler diagnostics to their output directories;
- final manuscript values should be generated from archived outputs rather than manually re-entered.

The production Stan seeds used for the historical arm-level analyses are 20260927 and 20260928 for the Gaussian and standardized $t_4$ fits, respectively.

---

# What is intentionally not in the repository

Generated result directories, compiled Stan executables, R session-history files, old nested ZIP archives, obsolete manuscript insertion fragments, and intermediate draft scripts are excluded. They are either reproducible from the tracked source or are superseded by the numbered production programs.

The `.gitignore` file prevents the major generated artifacts from being committed accidentally.

---

# Repository status

The historical continuous-outcome, Cochrane-13, 14-study temporal, and explanatory-regression workflows are represented by production source files. The contemporary binary workflow is intentionally marked as a template until the final verified binary extraction/design-adjusted variance file is frozen.

The Supplement should be treated as the detailed statistical audit trail for the manuscript; the code is the executable implementation of that audit trail.

---

## Correspondence

Brani Vidakovic  
Department of Statistics, Texas A&M University  
`brani@stat.tamu.edu`

For citation, please cite the associated manuscript. A repository DOI/citation can be added here if the repository is archived through Zenodo or a similar service.
