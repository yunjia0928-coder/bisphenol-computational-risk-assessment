# Bisphenol Computational Risk Assessment

Computational exposure and mixture risk assessment of bisphenols using R, ToxCast, high-throughput toxicokinetics (HTTK), _in vitro_-to-_in vivo_ extrapolation (IVIVE), and probabilistic modeling.

This repository presents the computational risk-assessment component of a peer-reviewed study on bisphenol A (BPA) and its structural analogs in plant-based foods. The study integrated dietary exposure data with high-throughput in vitro bioactivity and toxicokinetic modeling to derive human-equivalent points of departure and characterize potential mixture risks.

## Research Overview

Bisphenol A and its structural analogs are used in a wide range of consumer materials and may exhibit estrogenic activity. However, conventional animal toxicity data are limited for many BPA alternatives, making comprehensive risk assessment difficult.

The study investigated 16 bisphenols in 23 plant-based core foods in Taiwan and estimated dietary exposure across different age groups. To complement conventional toxicological evidence, ToxCast estrogen-receptor bioactivity data were integrated with high-throughput toxicokinetic modeling and IVIVE.

The computational workflow was used to:

1. integrate dietary exposure estimates across population groups;
2. identify estrogen-receptor-related ToxCast assays and bioactive concentrations;
3. simulate population toxicokinetics using the `httk` R package;
4. perform reverse dosimetry to convert _in vitro_ concentrations into oral equivalent doses;
5. derive human-equivalent points of departure (PODs);
6. compare _in vitro_–_in silico_ PODs with animal-derived PODs; and
7. characterize individual and mixture risks using margin-of-exposure (MOE) approaches and probabilistic modeling.

## Computational Workflow

```text
Bisphenol concentrations in foods
              ↓
Dietary exposure estimation
              ↓
ToxCast estrogen-receptor bioactivity
              ↓
Bioactive concentration (AC10)
              ↓
HTTK / population PBTK modeling in R
              ↓
IVIVE and reverse dosimetry
              ↓
Oral equivalent doses
              ↓
Human-equivalent points of departure
              ↓
Individual margin of exposure (MOE)
              ↓
Mixture margin of exposure (MOET)
              ↓
Probabilistic risk characterization
```

## HTTK and IVIVE Modeling

High-throughput toxicokinetic modeling was performed in R using the `httk` package.

The analysis incorporated chemical-specific toxicokinetic parameters and population variability to estimate steady-state plasma concentrations. Monte Carlo population simulations were used to characterize inter-individual toxicokinetic variability.

For reverse dosimetry, ToxCast AC10 values were converted to oral equivalent doses (OEDs) using the relationship:

```text
OED = AC10 × (1 mg/kg/day) / Css
```

where `Css` represents the simulated steady-state plasma concentration resulting from a unit oral dose.

The resulting oral equivalent doses were used to derive human-equivalent PODs for estrogen-receptor-mediated bioactivity.

## Mixture Risk Assessment

Risk characterization was conducted using the margin of exposure (MOE):

```text
MOE = POD / Dietary Intake
```

For bisphenols acting through the common estrogen-receptor pathway, individual MOEs were integrated into a mixture margin of exposure (MOET) using a reciprocal-sum approach.

Probabilistic simulations were used to account for variability in dietary exposure and toxicokinetics across population groups.

## Key Findings

The study found detectable concentrations of seven bisphenols in the analyzed plant-based foods, with BPA contributing the largest proportion of total detected bisphenol concentrations.

HTTK/IVIVE-derived points of departure were generally comparable to, or more conservative than, PODs derived from available animal studies.

At the population level, combined dietary exposure to the evaluated bisphenols was generally not considered a risk concern for estrogen-receptor pathway perturbation, although potential concerns could not be excluded for a small proportion of individuals under worst-case exposure scenarios.

## Repository Structure

```text
bisphenol-computational-risk-assessment/
│
├── README.md
├── .gitignore
└── httk_ivive_analysis.R
```

`httk_ivive_analysis.R` contains the R workflow used for the HTTK/IVIVE component of the published study.

## My Contributions

I contributed to the computational and quantitative components of the study, including formal analysis and visualization, as well as manuscript review and editing.

My contribution is reported in the published article under the CRediT authorship statement as:

- Formal analysis
- Visualization
- Writing – review & editing

## Related Publication

Chiu, C.-H., Sun, S.-H., Yao, Y.-J., Chuang, Y., Lee, Y.-T., & Lin, Y.-J. (2025). **Concentrations, composition profiles, and in vitro–in silico-based mixture risk assessment of bisphenol A and its analogs in plant-based foods.** *Environment International, 195*, 109229.

https://doi.org/10.1016/j.envint.2024.109229

## Key Methodological References

Pearce, R. G., Setzer, R. W., Strope, C. L., Sipes, N. S., & Wambaugh, J. F. (2017). **httk: R package for high-throughput toxicokinetics.** *Journal of Statistical Software, 79*(4), 1–26.

Lin, Y.-J., & Lin, Z. (2020). **In vitro-in silico-based probabilistic risk assessment of combined exposure to bisphenol A and its analogues by integrating ToxCast high-throughput in vitro assays with in vitro to in vivo extrapolation (IVIVE) via physiologically based pharmacokinetic (PBPK) modeling.** *Journal of Hazardous Materials, 399*, 122856.

## Data Availability

This repository focuses on the computational methodology and does not redistribute the underlying study dataset. Data availability follows the statement provided in the published article.

## Tools and Methods

- **R**
- **httk**
- **ToxCast / U.S. EPA CompTox**
- **High-throughput toxicokinetics (HTTK)**
- **Physiologically based toxicokinetic (PBTK) modeling**
- **In vitro-to-in vivo extrapolation (IVIVE)**
- **Monte Carlo simulation**
- **Probabilistic risk assessment**
- **Mixture risk assessment**
