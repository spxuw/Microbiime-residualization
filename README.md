# Masked by Age

**Uncovering age-independent and age-confounded microbial disease biomarkers using an age-residualization framework**

## Overview

Age is one of the strongest sources of biological and technical variation in human microbiome studies. Because many diseases are also age-associated, microbial features can appear disease-relevant simply because they track age. This creates a central problem for microbiome biomarker discovery: some taxa may represent disease-specific biology, whereas others may reflect age structure, age imbalance, or age-associated confounding.

**Masked by Age** provides an age-aware modeling framework to separate microbial disease signals into age-independent and age-confounded components. The framework asks a simple but important question:

> How much disease-predictive information remains in the microbiome after accounting for age?

Rather than treating age only as a covariate in a final model, this framework explicitly quantifies the incremental value of microbial features beyond an age-only baseline.

## Key idea

For each study or BioProject, we compare two predictive models for a binary disease phenotype:

1. **Age-only null model**  
   A baseline model using age alone to predict disease status.

2. **Full microbiome model**  
   A model using microbial taxa plus age to predict disease status.

The key metric is:

```text
Delta AUC = AUC(taxa + age model) - AUC(age-only model)
```

This value estimates the added disease-predictive information contributed by the microbiome beyond age.

- A **large positive Delta AUC** suggests that microbial profiles contain disease-associated signal not explained by age alone.
- A **small or near-zero Delta AUC** suggests that apparent microbial disease prediction may be largely age-confounded.
- A **negative Delta AUC** indicates that adding microbial features does not improve, and may destabilize, prediction beyond age.

## Conceptual interpretation

The framework distinguishes three broad biomarker scenarios:

### 1. Age-independent microbial disease biomarkers

These taxa remain predictive after accounting for age. They may represent disease-linked microbial shifts that are not simply explained by age differences between cases and controls.

### 2. Age-confounded microbial biomarkers

These taxa appear disease-associated because they are also strongly age-associated. Their predictive value may diminish once age is modeled explicitly.

### 3. Age-masked microbial biomarkers

Some disease-relevant microbial signals may be obscured by age-related variation. By modeling age explicitly, the framework helps reveal whether taxa provide additional signal beyond the dominant age effect.

## Workflow

The analysis is performed separately within each BioProject to avoid mixing study-specific effects.

### Step 1. Prepare input data

Each input CSV should contain:

| Column | Description |
|---|---|
| `Age` | Numeric age variable |
| `Phenotype` | Binary disease/control phenotype |
| `BioProject` | Study or BioProject identifier |
| Taxa columns | Microbial abundance features |

The current implementation focuses on BioProjects with binary phenotypes and sufficient sample size per class.

### Step 2. Run repeated stratified cross-validation

Within each BioProject, the framework uses repeated stratified K-fold cross-validation to preserve phenotype balance across folds.

### Step 3. Fit an age-only null model

The null model uses age alone:

```text
Phenotype ~ Age
```

This estimates how well age alone predicts disease status.

### Step 4. Fit full microbiome models

The full models use microbial taxa plus age:

```text
Phenotype ~ Taxa + Age
```

The current implementation supports:

- Random forest
- XGBoost

### Step 5. Estimate Delta AUC

For each fold, the framework computes:

```text
Delta AUC = AUC(full model) - AUC(age-only model)
```

Mean Delta AUC across repeated folds summarizes the age-adjusted incremental predictive value of the microbiome.

### Step 6. Assess significance by permutation

Phenotype labels are permuted to generate a null distribution of Delta AUC. One-sided permutation P values test whether observed Delta AUC exceeds what would be expected by chance.

### Step 7. Identify candidate biomarkers

Feature importance is averaged across cross-validation models rather than estimated from a single full-data model. This reduces overfitting and improves stability of candidate taxa prioritization.

## Outputs

The pipeline generates:

| Output | Description |
|---|---|
| `all_bioproject_metrics_binary.csv` | Per-BioProject AUROC, Delta AUC, permutation P values, and model summaries |
| `*_imp_phenotype_full_rf_top50.csv` | Top phenotype-associated taxa from random forest models |
| `*_imp_phenotype_full_xgb_top50.csv` | Top phenotype-associated taxa from XGBoost models |
| `*_imp_age_rf_top50.csv` | Top age-associated taxa from random forest age models |
| `*_imp_age_xgb_top50.csv` | Top age-associated taxa from XGBoost age models |
| `auroc.csv` | Summary AUROC results across datasets |

## Example usage

```r
results <- run_folder(
  folder = "../data/clean",
  k = 5,
  repeats = 50,
  seed = 345,
  min_n = 20,
  min_class = 10
)

write.csv(results, file = "../results/auroc.csv", row.names = FALSE)
```

## Methodological rationale

Many microbiome disease-classification studies report strong prediction accuracy without asking whether the signal is disease-specific or driven by demographic structure. This is especially problematic when disease cases are older than controls, or when control groups are not age-matched.

This framework treats age as a competing baseline predictor. By comparing microbiome models against an age-only null model, it quantifies whether microbial profiles add information beyond a major confounder.

In this sense, Delta AUC serves as an interpretable measure of **age-independent microbial disease information**.

## Recommended reporting

For each BioProject, report:

- sample size
- phenotype class balance
- age distribution by phenotype
- age-only AUROC
- full taxa-plus-age AUROC
- Delta AUC
- permutation P value
- top phenotype-associated taxa
- overlap between phenotype-associated and age-associated taxa

## Repository structure

```text
.
├── code/
│   └── run_age_residualization_models.R
├── data/
│   └── clean/
│       └── *.csv
├── results/
│   ├── auroc.csv
│   └── model_outputs_binary/
└── README.md
```

## Citation

If you use this framework, please cite:

**Masked by Age: uncovering age-independent and age-confounded microbial disease biomarkers using an age-residualization framework**

## Contact

For questions or collaborations, please contact:

**Xu-Wen Wang**  
Harvard Medical School / Brigham and Women's Hospital  
Channing Division of Network Medicine
