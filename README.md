# Proteomic Signatures of BMI Loss in UK Biobank

Analysis code for a study identifying plasma proteins associated with subsequent
BMI decline in UK Biobank, deriving a proteomic BMI-loss score, and testing that
score against multi-organ phenotypes and incident disease.

---

## Repository structure

```
.
├── README.md
├── R/
│   ├── 0_baseline_characteristics.R
│   ├── 1_discover_protein_signatures.R
│   ├── 1.1_subgroup_analysis.R
│   ├── 1.2_sensitivity_analysis.R
│   ├── 2_derive_bmi_loss_score.R
│   ├── 3_organ_system_associations.R
│   ├── 4_phenome_wide_association.R
│   ├── 5_joint_exposure_analysis.R
│   └── 6_model_discrimination.R
└── config/
    └── paths.R
```

### Scripts

| Script | Purpose |
|---|---|
| `0_baseline_characteristics.R` | Table 1 across the four training/validation datasets |
| `1_discover_protein_signatures.R` | Build the analytic dataset (`raw_data.csv`); screen plasma proteins against BMI change under two train/validation designs |
| `1.1_subgroup_analysis.R` | Protein–BMI-change associations by age, BMI, and self-reported weight change |
| `1.2_sensitivity_analysis.R` | Exclude ideal lifestyle, medication users, and prevalent/early disease |
| `2_derive_bmi_loss_score.R` | LASSO model; derive and evaluate the proteomic BMI-loss score |
| `3_organ_system_associations.R` | Score vs brain MRI, cerebellum, blood pressure, biochemistry, body composition, cardiac MRI |
| `4_phenome_wide_association.R` | Cox models across 103 incident disease outcomes |
| `5_joint_exposure_analysis.R` | Self-reported weight loss × score tertile, joint and separate |
| `6_model_discrimination.R` | Time-dependent AUC for four nested prediction models |

Script numbers follow the order results appear in the manuscript, not the order
of execution. `1_discover_protein_signatures.R` builds `raw_data.csv` before
screening proteins, so it must run first; `0_baseline_characteristics.R` reads
that file and therefore runs second. See the execution order below.

---

## Pipeline

```
              source CSVs (BMI, Olink, covariates, disease)
                                 │
                 1_discover_protein_signatures.R
                                 │
                          raw_data.csv
              ┌──────────────────┼──────────────────┐
              │                  │                  │
    0_baseline_          1.1_subgroup      2_derive_bmi_loss_score.R
    characteristics      1.2_sensitivity            │
                                            BMIscore_predict.csv
                                       ┌───────┬────┴───┬────────┐
                                       │       │        │        │
                                       3       4        5        6
```

`raw_data.csv` and `BMIscore_predict.csv` are the two hub files. Any change
upstream propagates to every downstream result, so both should be regenerated
together rather than piecemeal.

---

## Running the analysis

Execution order differs from the numbering — `1_` must run before `0_`, and
`2_` must complete before scripts `3_` through `6_`.

```r
source("config/paths.R")

# Build the analytic dataset and screen proteins
source("R/1_discover_protein_signatures.R")

# Descriptive table (depends on raw_data.csv)
source("R/0_baseline_characteristics.R")

# Robustness of the protein–BMI-change associations
source("R/1.1_subgroup_analysis.R")
source("R/1.2_sensitivity_analysis.R")

# Derive the BMI-loss score
source("R/2_derive_bmi_loss_score.R")

# Downstream applications of the score
source("R/3_organ_system_associations.R")
source("R/4_phenome_wide_association.R")
source("R/5_joint_exposure_analysis.R")
source("R/6_model_discrimination.R")
```

### Configuration

All file paths live in `config/paths.R` rather than being hard-coded in each
script, so the pipeline can be relocated by editing a single file.

```r
# config/paths.R
PROJECT_ROOT <- "D:/科研数据/课题数据/1. UKB数据课题/"

path_bmi       <- file.path(PROJECT_ROOT, "UKB数据常用汇总/Obesity data/BMI.csv")
path_protein   <- file.path(PROJECT_ROOT, "UKB数据常用汇总/Omics data/protein_UKB_filled.csv")
path_covariate <- file.path(PROJECT_ROOT, "UKB数据常用汇总/covariates_impute.csv")
path_disease   <- file.path(PROJECT_ROOT, "All disease/UBK_disease.csv")

dir_results <- file.path(PROJECT_ROOT, "BMI Loss 课题/更新后结果")
dir_raw     <- file.path(PROJECT_ROOT, "BMI Loss 课题/原始数据")
```

---

## Requirements

R ≥ 4.1.

```r
install.packages(c(
  "dplyr", "survival", "glmnet", "ggplot2",
  "tableone", "flextable", "officer", "table1",
  "riskRegression", "prodlim", "pec"
))
```

| Package | Used by |
|---|---|
| `dplyr` | all scripts |
| `survival` | `4_`, `5_`, `6_` |
| `glmnet` | `2_` |
| `ggplot2`, `table1` | `2_` |
| `tableone`, `flextable`, `officer` | `0_` |
| `riskRegression`, `prodlim`, `pec` | `6_` |

---

## Data availability

Individual-level UK Biobank data cannot be redistributed. Access is through the
UK Biobank Access Management System (https://www.ukbiobank.ac.uk). This
repository contains analysis code only; no data files are included, and
`.gitignore` excludes `*.csv`, `*.docx`, and `*.pdf` to prevent accidental
commits of participant data.

---

## Method summary

**Analytic sample.** UK Biobank participants with Olink proteomic measurements at
baseline, BMI ≥ 18.5 kg/m² at recruitment, and positive follow-up time for the
composite cancer and non-cancer disease outcomes.

**Exposure derivation.** Plasma proteins were screened against BMI change under
two designs: an assessment-time split (train on `BMI_2 − BMI_0`, validate on
`BMI_1 − BMI_0` among those without a second follow-up) and a 60:40 random
split. Proteins significant in both training and validation entered a LASSO
model, and the resulting coefficients define the proteomic BMI-loss score.

**Downstream analyses.** The score was tested against brain and cardiac MRI
measures, blood pressure, clinical biochemistry, and bioelectrical impedance
body composition; against 103 incident disease outcomes in Cox models; and in
combination with self-reported weight change. Model discrimination was compared
across four nested Cox models using time-dependent AUC at 15 years.

**Multiple testing.** FDR correction is applied within each analysis family.

---

## Citation

*Manuscript in preparation. Citation details will be added on publication.*

## License

Code released under the MIT License. Use of UK Biobank data is governed
separately by the terms of the approved application.

## Contact

Please open an issue for questions about the code.
