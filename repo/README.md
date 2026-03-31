ERP Study Repository

This repository contains the data and analysis scripts for the study:

From the app to the lexicon: Neurocognitive markers of word integration in a second language

The project investigates ERP responses (e.g., N400, LPC) in a semantic decision task, including both EEG and behavioral data analyses conducted in R.

📁 repo/
├── data/
│   ├── eeg_epochs/        # Preprocessed EEG data (.fif files; one per subject; N=17)
│   ├── eeg_data.csv       # Dataset used for ERP statistical analyses
│   ├── behavioral.csv     # Behavioral dataset (RTs, accuracy; N=24)
│
├── scripts/
│   ├── main_analysis.R    # Main statistical analyses (GLMM / LMM)
│
├── results/
│   ├── tables.pdf         # Summary tables reported in the manuscript
│   ├── tables.xlsx        # Same tables in editable format
│
├── README.md

⚙️ Requirements
EEG data

EEG data are provided as preprocessed epoch files in .fif format, compatible with MNE-Python.
Preprocessing steps (artifact rejection, ICA, filtering) are described in the manuscript.

R environment
All statistical analyses were conducted in R (version X.X.X).

Required packages:
install.packages(c(
  "tidyverse",
  "lme4",
  "lmerTest",
  "emmeans",
  "DHARMa",
  "glmmTMB"
))

▶️ How to reproduce the analyses
Download or clone this repository.
Open R or RStudio and set the working directory:
setwd("path/to/repository")
Run the main analysis script:
source("scripts/main_analysis.R")

This script:

Loads the ERP and behavioral datasets
Fits linear mixed-effects models
Generates the results reported in the manuscript

📊 Data description
EEG data
17 participants
Preprocessed epoch files (.fif)
One file per subject
ERP dataset (eeg_data.csv)

Contains trial-level ERP measures used in statistical analyses.

Key variables:

subject: participant ID
item: stimulus ID
word_type: experimental condition
relatedness: Related / Unrelated
mean_amplitude: ERP amplitude (µV)
time_window: e.g., N400 (300–500 ms), LPC (500–800 ms)

Behavioral dataset (behavioral.csv)

Includes:

Reaction times (RTs)
Accuracy
📈 Statistical analysis

ERP and behavioral data were analyzed using linear mixed-effects models.

Fixed effects: Word Type (or Language), Relatedness, and their interaction
Random effects: Subject and Item

Post hoc comparisons were computed using emmeans.

🔍 Reproducibility notes
The repository includes preprocessed EEG data only.
All analyses reported in the manuscript can be reproduced using the provided datasets and scripts.
Minor differences may occur depending on R and package versions.

📬 Contact

For questions or issues:

[Your name]
[Your email]