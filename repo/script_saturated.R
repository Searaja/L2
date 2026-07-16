# LIBRARIES

library(modelsummary)
library(lme4)
library(lmerTest)
library(dplyr)
library(glmmTMB)
library(emmeans)
library(ggplot2)
library(DHARMa)

library(broom.mixed)
library(flextable)
library(officer)
library(pagedown)

# TRUE = use pre-decided RE formulas (fast); FALSE = run full allFit/rePCA/LRT pipeline
LOCK_MODELS <- TRUE


# BEHAVIORAL DATA

df <- read.csv("data/behavioral.csv", header = TRUE)

df$subject           <- factor(df$file)
df$word_type         <- factor(df$categ)
df$relatedness       <- factor(df$corrAns)
df$accuracy          <- factor(df$correct)
df$item              <- factor(df$word_prime_merge)
levels(df$word_type)[levels(df$word_type) == 1] <- "L2-Remote"
levels(df$word_type)[levels(df$word_type) == 2] <- "L2-Recent"
levels(df$word_type)[levels(df$word_type) == 3] <- "L1"
levels(df$relatedness)[levels(df$relatedness) == "right"] <- "Related"
levels(df$relatedness)[levels(df$relatedness) == "left"]  <- "Unrelated"
df$word_type <- factor(df$word_type, levels = c("L1", "L2-Remote", "L2-Recent"))

# Planned contrasts — levels: L1, L2-Remote, L2-Recent
contrasts(df$word_type) <- cbind(
  Remote_vs_Recent = c(0,    0.5,  -0.5),
  L1_vs_L2         = c(2/3, -1/3,  -1/3)
)
# Levels: Unrelated, Related → β > 0 means Unrelated > Related
contrasts(df$relatedness) <- cbind(Unrelated_vs_Related = c(0.5, -0.5))

# =============================================================================
# RT OUTLIER VISUALIZATION
# Run these plots first, then set manual exclusions below before fitting models
# =============================================================================

print(
  ggplot(df, aes(x = subject, y = rt)) +
    geom_boxplot(outlier.colour = "red", outlier.size = 1) +
    facet_wrap(~ word_type, ncol = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7)) +
    labs(title = "RT by subject", y = "RT (s)")
)

print(
  ggplot(df, aes(x = rt, fill = relatedness)) +
    geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
    facet_wrap(~ word_type) +
    theme_minimal() +
    labs(title = "RT distribution by condition", x = "RT (ms)")
)

# --- MANUAL EXCLUSION — adjust after inspecting plots above ---
rt_min <- 0.2    # lower cutoff (s)
rt_max <- 3.0    # upper cutoff (s)
exclude_subjects_rt <- c()  # e.g. c("s01", "s05")

df_filt <- df %>%
  filter(rt >= rt_min, rt <= rt_max) %>%
  filter(!subject %in% exclude_subjects_rt) %>%
  droplevels()

contrasts(df_filt$word_type) <- cbind(
  Remote_vs_Recent = c(0,    0.5,  -0.5),
  L1_vs_L2         = c(2/3, -1/3,  -1/3)
)
contrasts(df_filt$relatedness) <- cbind(Unrelated_vs_Related = c(0.5, -0.5))


# =============================================================================
# MODEL 1: ACCURACY — Binomial GLMM (glmer)
# Reduced RE: (1 + word_type + relatedness | subject) + (1 | item)
# Rationale: rePCA(sat) rank 4/6; LRT base→red p<.001, red→sat p=.655
# Optimizer: Nelder_Mead (bobyqa singular; Nelder_Mead converged clean)
# =============================================================================

if (LOCK_MODELS) {
  model_Acc_red <- glmer(
    accuracy ~ word_type * relatedness +
      (1 + word_type + relatedness | subject) + (1 | item),
    data = df, family = binomial,
    control = glmerControl(optimizer = "Nelder_Mead"))
} else {
  model_Acc_base <- glmer(
    accuracy ~ word_type * relatedness +
      (1 | subject) + (1 | item),
    data = df, family = binomial)

  model_Acc_sat <- glmer(
    accuracy ~ word_type * relatedness +
      (1 + word_type * relatedness | subject) + (1 | item),
    data = df, family = binomial)

  fits_acc      <- allFit(model_Acc_sat)
  ss_acc        <- summary(fits_acc)
  print(ss_acc$which.OK)
  model_Acc_sat <- fits_acc[[which(ss_acc$which.OK)[2]]]

  VarCorr(model_Acc_sat)
  rePCA(model_Acc_sat)
  isSingular(model_Acc_sat)

  model_Acc_red <- glmer(
    accuracy ~ word_type * relatedness +
      (1 + word_type + relatedness | subject) + (1 | item),
    data = df, family = binomial)

  fits_acc_red  <- allFit(model_Acc_red)
  ss_acc_red    <- summary(fits_acc_red)
  print(ss_acc_red$which.OK)
  model_Acc_red <- fits_acc_red[[which(ss_acc_red$which.OK)[1]]]

  VarCorr(model_Acc_red)
  rePCA(model_Acc_red)
  isSingular(model_Acc_red)

  anova(model_Acc_base, model_Acc_red, model_Acc_sat)
}

summary(model_Acc_red)

sim_res_acc <- simulateResiduals(model_Acc_red)
dev.new(); plot(sim_res_acc)


# =============================================================================
# MODEL 2: RT — LMM on log-RT (lmer, replaces nlme::lme)
# Reduced RE: (1 + relatedness | subject) + (1 | item)
# Rationale: rePCA(sat) rank 5/6; word_type slopes r=.997 → dropped
#            LRT base→red p<.001; sat singular → red is reported model
# =============================================================================

if (LOCK_MODELS) {
  model_RT_red <- lmer(
    log(rt) ~ word_type * relatedness +
      (1 + relatedness | subject) + (1 | item),
    data = df_filt, REML = TRUE)
} else {
  model_RT_base <- lmer(
    log(rt) ~ word_type * relatedness +
      (1 | subject) + (1 | item),
    data = df_filt, REML = TRUE)

  model_RT_sat <- lmer(
    log(rt) ~ word_type * relatedness +
      (1 + word_type * relatedness | subject) + (1 | item),
    data = df_filt, REML = TRUE)

  fits_rt      <- allFit(model_RT_sat)
  ss_rt        <- summary(fits_rt)
  print(ss_rt$which.OK)
  model_RT_sat <- fits_rt[[which(ss_rt$which.OK)[1]]]

  VarCorr(model_RT_sat)
  rePCA(model_RT_sat)
  isSingular(model_RT_sat)

  model_RT_red <- lmer(
    log(rt) ~ word_type * relatedness +
      (1 + relatedness | subject) + (1 | item),
    data = df_filt, REML = TRUE)

  fits_rt_red  <- allFit(model_RT_red)
  ss_rt_red    <- summary(fits_rt_red)
  print(ss_rt_red$which.OK)
  model_RT_red <- fits_rt_red[[which(ss_rt_red$which.OK)[1]]]

  VarCorr(model_RT_red)
  rePCA(model_RT_red)
  isSingular(model_RT_red)

  anova(model_RT_base, model_RT_red, model_RT_sat)
}

summary(model_RT_red)

res_rt <- residuals(model_RT_red)
qqnorm(res_rt); qqline(res_rt)
plot(fitted(model_RT_red), res_rt); abline(h = 0)



# =============================================================================
# EEG DATA
# =============================================================================

data <- read.csv("data/epoch_amplitudes.csv", header = TRUE)

data$subject    <- as.factor(data$subject)
data$item       <- as.factor(data$item)
data$relatedness <- as.factor(data$relatedness)
data$word_type  <- as.factor(data$language)
data$channel    <- as.factor(data$channel)

levels(data$word_type)[levels(data$word_type) == "rem"] <- "L2-Remote"
levels(data$word_type)[levels(data$word_type) == "rec"] <- "L2-Recent"
levels(data$word_type)[levels(data$word_type) == "esp"] <- "L1"
levels(data$relatedness)[levels(data$relatedness) == 1] <- "Related"
levels(data$relatedness)[levels(data$relatedness) == 0] <- "Unrelated"

# --- MANUAL EXCLUSION cutoffs (applied per window×electrode below) ---
rt_min_eeg           <- 0.2 * 256  # lower RT cutoff (s * fs)
rt_max_eeg           <- 3.0 * 256  # upper RT cutoff (s * fs)
amp_cutoff           <- 15          # exclude |mean_amplitude| > this (µV)
exclude_subjects_eeg <- c()         # e.g. c("s01", "s05")

# =============================================================================
# MODELS 3–6: EEG AMPLITUDES — Gaussian GLMM (glmmTMB)
# Note: rePCA() is lme4-specific and not available for glmmTMB.
#       Use VarCorr() and inspect near-zero SDs manually.
#
# RE SIMPLIFICATION RESULTS BY WINDOW/ELECTRODE:
# N400 Cz (0.3, Cz): base model (1|subject) + (1|item)
#   - sat: near-perfect correlations throughout (r = .984, -.997, .987)
#   - (1 + word_type + relatedness | subject): L2-Remote/relatedness r = -.983
#   - (1 + word_type | subject): singular convergence (non-positive Hessian)
#   - (1 + relatedness | subject): Intercept/relatedness r = -1.000
#   - (1|subject) + (0 + relatedness|subject): relatedness SD ≈ 0 (.00032)
#   → reported model = base
# LPC  Cz (0.5, Cz): base model (1|subject) + (1|item)
#   - sat: Intercept/L2-Recent r = -1.000, L2-Recent/L2-Remote r = 1.000
#   - (1 + word_type + relatedness | subject): all correlations at ±1.000
#   → reported model = base
# N400 Pz (0.3, Pz): (1 + word_type || subject), no item
#   - sat: near-perfect correlations (r = -.996, .998, .981); item SD ≈ 0
#   - (1 + word_type + relatedness | subject): Intercept/relatedness r = .963; item SD ≈ 0
#   - (1 + word_type | subject): singular convergence; item SD ≈ 0
#   - (1 + word_type || subject): converged, clean VarCorr
#   → reported model = (1 + word_type || subject), no (1|item)
# LPC  Pz (0.5, Pz): base model (1|subject) + (1|item)
#   - sat: singular convergence
#   - most reduced converging model: (1 + word_type || subject) + (1|item)
#   - AIC/BIC of reduced worse than base → base preferred
#   → reported model = base
# =============================================================================

# Locked RE formulas — established from VarCorr/AIC diagnostics (see comments above)
eeg_locked_formulas <- list(
  N400_Cz = mean_amplitude ~ word_type * relatedness + (1 | subject) + (1 | item),
  LPC_Cz  = mean_amplitude ~ word_type * relatedness + (1 | subject) + (1 | item),
  N400_Pz = mean_amplitude ~ word_type * relatedness + (1 + word_type || subject),
  LPC_Pz  = mean_amplitude ~ word_type * relatedness + (1 | subject) + (1 | item)
)

eeg_configs <- list(
  N400_Cz = list(win = 0.3, ch = "Cz"),
  LPC_Cz  = list(win = 0.5, ch = "Cz"),
  N400_Pz = list(win = 0.3, ch = "Pz"),
  LPC_Pz  = list(win = 0.5, ch = "Pz")
)

eeg_models <- list()

for (nm in names(eeg_configs)) {
  cfg <- eeg_configs[[nm]]

  d <- data %>%
    filter(window_start == cfg$win, channel == cfg$ch) %>%
    filter(RT >= rt_min_eeg, RT <= rt_max_eeg) %>%
    filter(abs(mean_amplitude) <= amp_cutoff) %>%
    filter(!subject %in% exclude_subjects_eeg) %>%
    select(mean_amplitude, word_type, relatedness, subject, item) %>%
    na.omit() %>%
    droplevels()

  contrasts(d$word_type) <- cbind(
    Remote_vs_Recent = c(0,   -0.5,  0.5),
    L1_vs_L2         = c(2/3, -1/3, -1/3)
  )
  contrasts(d$relatedness) <- cbind(Unrelated_vs_Related = c(0.5, -0.5))

  if (LOCK_MODELS) {
    m_red <- glmmTMB(
      eeg_locked_formulas[[nm]],
      dispformula = ~ word_type * relatedness,
      data = d)
  } else {
    m_base <- glmmTMB(
      mean_amplitude ~ word_type * relatedness + (1 | subject) + (1 | item),
      dispformula = ~ word_type * relatedness, data = d)
    m_sat <- glmmTMB(
      mean_amplitude ~ word_type * relatedness +
        (1 + word_type * relatedness | subject) + (1 | item),
      dispformula = ~ word_type * relatedness, data = d)
    cat("\n---", nm, "VarCorr(sat) ---\n")
    print(VarCorr(m_sat))
    cat("\n---", nm, "LRT base vs sat ---\n")
    print(anova(m_base, m_sat))
    m_red <- m_base  # update to winning formula after inspecting output above
  }

  cat("\n===", nm, "===\n")
  print(summary(m_red))

  sim_res <- simulateResiduals(m_red)
  dev.new(); plot(sim_res, main = nm)

  res_e <- residuals(m_red)
  qqnorm(res_e, main = paste("QQ —", nm)); qqline(res_e)
  plot(fitted(m_red), res_e, main = paste("Fitted vs Resid —", nm)); abline(h = 0)

  eeg_models[[nm]] <- list(model = m_red, data = d)
}


# =============================================================================
# SUPPLEMENTARY TABLES
# Models: model_Acc_red (glmer), model_RT_red (lmer), eeg_models (glmmTMB × 4)
# Tables: fixed effects, random effects, Type II ANOVA only (no post-hoc)
# =============================================================================

fmt_p <- function(p) {
  case_when(
    p < .001 ~ "<.001",
    p < .01  ~ "<.01",
    p < .05  ~ "<.05",
    TRUE     ~ as.character(round(p, 3))
  )
}

make_ft <- function(df, caption = "") {
  flextable(df) %>%
    set_caption(caption) %>%
    bold(part = "header") %>%
    hline_top(part = "header", border = fp_border(width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border(width = 1)) %>%
    hline_bottom(part = "body",   border = fp_border(width = 1.5)) %>%
    autofit() %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "body")
}

# -----------------------------------------------------------------------------
# 1. model_Acc_red — Binomial GLMM (glmer)
# -----------------------------------------------------------------------------

tidy_acc <- tidy(model_Acc_red, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

tab_acc_fixed <- tidy_acc %>%
  transmute(
    `Predictor`    = term,
    `B`            = round(estimate, 3),
    `SE`           = round(std.error, 3),
    `95% CI lower` = round(conf.low, 3),
    `95% CI upper` = round(conf.high, 3),
    `z`            = round(statistic, 3),
    `p`            = fmt_p(p.value)
  )

tidy_acc_re <- tidy(model_Acc_red, effects = "ran_pars")
tab_acc_re <- tidy_acc_re %>%
  transmute(`Group` = group, `Term` = term, `SD` = round(estimate, 3))

# -----------------------------------------------------------------------------
# 2. model_RT_red — LMM (lmer, log-RT)
# -----------------------------------------------------------------------------

tidy_rt <- tidy(model_RT_red, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

tab_rt_fixed <- tidy_rt %>%
  transmute(
    `Predictor`    = term,
    `B`            = round(estimate, 3),
    `SE`           = round(std.error, 3),
    `95% CI lower` = round(conf.low, 3),
    `95% CI upper` = round(conf.high, 3),
    `gl`           = round(df, 1),
    `t`            = round(statistic, 3),
    `p`            = fmt_p(p.value)
  )

tidy_rt_re <- tidy(model_RT_red, effects = "ran_pars")
tab_rt_re <- tidy_rt_re %>%
  transmute(`Group` = group, `Term` = term, `SD` = round(estimate, 3))

# -----------------------------------------------------------------------------
# 3–6. eeg_models — Gaussian GLMM (glmmTMB), one per window×electrode
# -----------------------------------------------------------------------------

eeg_table_prefixes <- c(N400_Cz = "S3", LPC_Cz = "S4", N400_Pz = "S5", LPC_Pz = "S6")

eeg_table_list <- list()

for (nm in names(eeg_models)) {
  m      <- eeg_models[[nm]]$model
  label  <- gsub("_", " ", nm)
  prefix <- eeg_table_prefixes[nm]

  tidy_e  <- tidy(m, effects = "fixed",    component = "cond", conf.int = TRUE, conf.level = 0.95)
  tidy_re <- tidy(m, effects = "ran_pars", component = "cond")

  tab_fixed <- tidy_e %>%
    transmute(
      `Predictor`    = term,
      `B`            = round(estimate, 3),
      `SE`           = round(std.error, 3),
      `95% CI lower` = round(conf.low, 3),
      `95% CI upper` = round(conf.high, 3),
      `z`            = round(statistic, 3),
      `p`            = fmt_p(p.value)
    )

  tab_re <- tidy_re %>%
    transmute(`Group` = group, `Term` = term, `SD` = round(estimate, 3))

  eeg_table_list[[paste0(nm, "_fixed")]] <-
    make_ft(tab_fixed, sprintf("Table %sa. Fixed effects — EEG %s", prefix, label))
  eeg_table_list[[paste0(nm, "_re")]] <-
    make_ft(tab_re, sprintf("Table %sb. Random effects — EEG %s", prefix, label))
}

# -----------------------------------------------------------------------------
# EXPORT TO PDF (requires webshot2 + Chrome)
# -----------------------------------------------------------------------------

tables_pdf <- c(
  list(
    "Table S1a. Fixed effects — Accuracy"  = make_ft(tab_acc_fixed, "Table S1a. Fixed effects — Accuracy"),
    "Table S1b. Random effects — Accuracy" = make_ft(tab_acc_re,    "Table S1b. Random effects — Accuracy"),
    "Table S2a. Fixed effects — RT"        = make_ft(tab_rt_fixed,  "Table S2a. Fixed effects — RT"),
    "Table S2b. Random effects — RT"       = make_ft(tab_rt_re,     "Table S2b. Random effects — RT")
  ),
  eeg_table_list
)

do.call(save_as_html, c(tables_pdf, list(path = "supplementary_tables_saturated.html")))
chrome_print("supplementary_tables_saturated.html",
             output = "supplementary_tables_saturated.pdf")
message("✓ File saved: supplementary_tables_saturated.pdf")
