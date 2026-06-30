# LIBRARIES

library(modelsummary)
library(car)
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

k <- 5

df_filt <- df %>%
  group_by(corrAns, categ) %>%
  filter(abs(rt - mean(rt, na.rm = TRUE)) < k * sd(rt, na.rm = TRUE)) %>%
  ungroup()


# =============================================================================
# MODEL 1: ACCURACY — Binomial GLMM (glmer)
# =============================================================================

model_Acc_base <- glmer(
  accuracy ~ word_type * relatedness +
    (1 | subject) + (1 | item),
  data = df, family = binomial)

model_Acc_sat <- glmer(
  accuracy ~ word_type * relatedness +
    (1 + word_type * relatedness | subject) +
    (1 | item),
  data = df, family = binomial)

# Try all available optimizers and pick the best-converging one
fits_acc      <- allFit(model_Acc_sat)
ss_acc        <- summary(fits_acc)
print(ss_acc$which.OK)                          # which optimizers converged
model_Acc_sat <- fits_acc[[which(ss_acc$which.OK)[1]]]  # use first converging

# RE diagnostics on saturated
VarCorr(model_Acc_sat)
rePCA(model_Acc_sat)
isSingular(model_Acc_sat)

# Reduced model: drop interaction from RE (rePCA shows rank 4, not 6)
model_Acc_red <- glmer(
  accuracy ~ word_type * relatedness +
    (1 + word_type + relatedness | subject) +
    (1 | item),
  data = df, family = binomial)

fits_acc_red      <- allFit(model_Acc_red)
ss_acc_red        <- summary(fits_acc_red)
print(ss_acc_red$which.OK)
model_Acc_red     <- fits_acc_red[[which(ss_acc_red$which.OK)[1]]]

VarCorr(model_Acc_red)
rePCA(model_Acc_red)
isSingular(model_Acc_red)

# LRT chain: base -> reduced -> saturated
anova(model_Acc_base, model_Acc_red, model_Acc_sat)

summary(model_Acc_red)
Anova(model_Acc_red, type = "II")

sim_res_acc <- simulateResiduals(model_Acc_red)
plot(sim_res_acc)

# Interactions
p <- emmip(model_Acc_red, relatedness ~ word_type)
print(p + ggplot2::ggtitle("Accuracy — reduced RE"))
emmeans(model_Acc_red, pairwise ~ word_type,               type = "response")
emmeans(model_Acc_red, pairwise ~ word_type | relatedness, type = "response")
pairs(emmeans(model_Acc_red, ~ relatedness | word_type,    type = "response"))


# =============================================================================
# MODEL 2: RT — LMM on log-RT (lmer, replaces nlme::lme)
# Correct answers only, RT < 2.5 SD from mean — same filtering as script.R
# Note: heterogeneous variance (varIdent) is dropped; focus is RE structure
# =============================================================================

model_RT_base <- lmer(
  log(rt) ~ word_type * relatedness +
    (1 | subject) + (1 | item),
  data = df_filt, REML = TRUE)

model_RT_sat <- lmer(
  log(rt) ~ word_type * relatedness +
    (1 + word_type * relatedness | subject) +
    (1 | item),
  data = df_filt, REML = TRUE)

# Try all available optimizers and pick the best-converging one
fits_rt      <- allFit(model_RT_sat)
ss_rt        <- summary(fits_rt)
print(ss_rt$which.OK)
model_RT_sat <- fits_rt[[which(ss_rt$which.OK)[1]]]

# RE diagnostics on saturated
VarCorr(model_RT_sat)
rePCA(model_RT_sat)
isSingular(model_RT_sat)

# Reduced model: rePCA(sat) rank 5/6, then rank 3/4 on first reduction;
# word_type slopes correlated at r=.997 → drop word_type from RE
model_RT_red <- lmer(
  log(rt) ~ word_type * relatedness +
    (1 + relatedness | subject) +
    (1 | item),
  data = df_filt, REML = TRUE)

fits_rt_red      <- allFit(model_RT_red)
ss_rt_red        <- summary(fits_rt_red)
print(ss_rt_red$which.OK)
model_RT_red     <- fits_rt_red[[which(ss_rt_red$which.OK)[1]]]

VarCorr(model_RT_red)
rePCA(model_RT_red)
isSingular(model_RT_red)

# LRT chain: base -> reduced -> saturated
anova(model_RT_base, model_RT_red, model_RT_sat)

summary(model_RT_red)

res_rt <- residuals(model_RT_red)
qqnorm(res_rt); qqline(res_rt)
plot(fitted(model_RT_red), res_rt); abline(h = 0)

Anova(model_RT_red, type = "II")

p <- emmip(model_RT_red, relatedness ~ word_type)
print(p + ggplot2::ggtitle("RT — reduced RE"))
emmeans(model_RT_red, pairwise ~ word_type,               type = "response")
emmeans(model_RT_red, pairwise ~ word_type | relatedness, type = "response")
emmeans(model_RT_red, pairwise ~ relatedness,             type = "response")
emmeans(model_RT_red, pairwise ~ relatedness | word_type, type = "response")


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

# CHANGE window_start AND channel TO SEE RESULTS FOR DIFFERENT
# TIME WINDOWS AND ELECTRODE SITES:
# 0.3 FOR N400 (300-500 ms); 0.5 FOR LPC (500-800 ms)

data_filt <- data %>%
  filter(window_start == 0.5, channel == "Pz")

data_filt <- data_filt %>%
  group_by(subject, word_type, relatedness) %>%
  mutate(
    media_rt = mean(RT, na.rm = TRUE),
    sd_rt    = sd(RT, na.rm = TRUE)
  ) %>%
  filter(RT > media_rt - k * sd_rt,
         RT < media_rt + k * sd_rt) %>%
  ungroup()

data_filt <- data_filt %>%
  select(mean_amplitude, word_type, relatedness, subject, item) %>%
  na.omit()


# =============================================================================
# MODEL 3: EEG AMPLITUDES — Gaussian GLMM (glmmTMB)
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

model_eeg_base <- glmmTMB(
  mean_amplitude ~ word_type * relatedness +
    (1 | subject) + (1 | item),
  dispformula = ~ word_type * relatedness,
  data = data_filt)

model_eeg_sat <- glmmTMB(
  mean_amplitude ~ word_type * relatedness +
    (1 + word_type * relatedness | subject) +
    (1 | item),
  dispformula = ~ word_type * relatedness,
  data = data_filt)

# RE diagnostics on saturated (rePCA and allFit not available for glmmTMB)
VarCorr(model_eeg_sat)

# RE simplification: VarCorr(sat) showed near-perfect correlations throughout.
# Stepwise reduction (word_type+relatedness, word_type only, relatedness only,
# relatedness with ||) all produced singular convergence or SD≈0.
# Data supports intercepts only → reported model = base.
model_eeg_red <- model_eeg_base

# LRT: base vs saturated (documents over-parameterization of sat)
anova(model_eeg_base, model_eeg_sat)

summary(model_eeg_red)
Anova(model_eeg_red, type = "II")

sim_res_eeg <- simulateResiduals(model_eeg_red)
plot(sim_res_eeg)

emm_eeg <- emmeans(model_eeg_red, ~ word_type * relatedness)
plot(emm_eeg)

emmeans(model_eeg_red, pairwise ~ word_type,               adjust = "tukey")
emmeans(model_eeg_red, pairwise ~ relatedness,             adjust = "tukey")
emmeans(model_eeg_red, pairwise ~ word_type | relatedness, adjust = "tukey")
emmeans(model_eeg_red, pairwise ~ relatedness | word_type, adjust = "tukey")

res_eeg <- residuals(model_eeg_red)
qqnorm(res_eeg); qqline(res_eeg)
plot(fitted(model_eeg_red), res_eeg); abline(h = 0)


# =============================================================================
# SUPPLEMENTARY TABLES
# Models: model_Acc_red (glmer), model_RT_red (lmer), model_eeg_red (glmmTMB)
# =============================================================================

round_df <- function(df, digits = 3) {
  df %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
}

fmt_p <- function(p) {
  case_when(
    p < .001 ~ "<.001",
    p < .01  ~ "<.01",
    p < .05  ~ "<.05",
    TRUE     ~ as.character(round(p, 3))
  )
}

make_ft <- function(df, caption = "") {
  ft <- flextable(df) %>%
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
  ft
}

# -----------------------------------------------------------------------------
# 1. model_Acc_red — Binomial GLMM (glmer)
# -----------------------------------------------------------------------------

tidy_acc <- tidy(model_Acc_red, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

tab_acc_fixed <- tidy_acc %>%
  transmute(
    `Predictor`    = term,
    `B`            = estimate,
    `SE`           = std.error,
    `95% CI lower` = conf.low,
    `95% CI upper` = conf.high,
    `z`            = statistic,
    `p`            = fmt_p(p.value)
  ) %>%
  round_df(3)
tab_acc_fixed$p <- fmt_p(tidy_acc$p.value)

tidy_acc_re <- tidy(model_Acc_red, effects = "ran_pars")
tab_acc_re <- tidy_acc_re %>%
  transmute(`Group` = group, `Term` = term, `SD` = round(estimate, 3))

anova_acc <- Anova(model_Acc_red, type = "II") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Effect") %>%
  rename(`Chi²` = Chisq, `gl` = Df, `p` = `Pr(>Chisq)`) %>%
  mutate(`Chi²` = round(`Chi²`, 3), `p` = fmt_p(`p`))

# Contrasts
emm_acc_wt   <- emmeans(model_Acc_red, pairwise ~ word_type, type = "response")
tab_acc_c1   <- as.data.frame(summary(emm_acc_wt$contrasts)) %>%
  transmute(`Contrast` = contrast, `OR` = round(odds.ratio, 3),
            `SE` = round(SE, 3), `z` = round(z.ratio, 3),
            `p adjusted` = fmt_p(p.value))

emm_acc_wt_r <- emmeans(model_Acc_red, pairwise ~ word_type | relatedness, type = "response")
tab_acc_c2   <- as.data.frame(summary(emm_acc_wt_r$contrasts)) %>%
  transmute(`Contrast` = contrast, `Relatedness` = relatedness,
            `OR` = round(odds.ratio, 3), `SE` = round(SE, 3),
            `z` = round(z.ratio, 3), `p adjusted` = fmt_p(p.value))

emm_acc_r_wt <- emmeans(model_Acc_red, pairwise ~ relatedness | word_type, type = "response")
tab_acc_c3   <- as.data.frame(summary(emm_acc_r_wt$contrasts)) %>%
  transmute(`Contrast` = contrast, `Word type` = word_type,
            `OR` = round(odds.ratio, 3), `SE` = round(SE, 3),
            `z` = round(z.ratio, 3), `p adjusted` = fmt_p(p.value))

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

anova_rt <- Anova(model_RT_red, type = "II") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Effect") %>%
  rename(`Chi²` = Chisq, `df` = Df, `p_raw` = `Pr(>Chisq)`) %>%
  mutate(`Chi²` = round(`Chi²`, 3), `p` = fmt_p(p_raw)) %>%
  select(Effect, `Chi²`, df, p)

# Contrasts
emm_rt_wt   <- emmeans(model_RT_red, pairwise ~ word_type, type = "response")
tab_rt_c1   <- as.data.frame(summary(emm_rt_wt$contrasts)) %>%
  transmute(`Contrast` = contrast, `Ratio` = round(ratio, 3),
            `SE` = round(SE, 3), `gl` = round(df, 1),
            `t` = round(t.ratio, 3), `p adjusted` = fmt_p(p.value))

emm_rt_wt_r <- emmeans(model_RT_red, pairwise ~ word_type | relatedness, type = "response")
tab_rt_c2   <- as.data.frame(summary(emm_rt_wt_r$contrasts)) %>%
  transmute(`Contrast` = contrast, `Relatedness` = relatedness,
            `Ratio` = round(ratio, 3), `SE` = round(SE, 3),
            `gl` = round(df, 1), `t` = round(t.ratio, 3),
            `p adjusted` = fmt_p(p.value))

emm_rt_r    <- emmeans(model_RT_red, pairwise ~ relatedness, type = "response")
tab_rt_c3   <- as.data.frame(summary(emm_rt_r$contrasts)) %>%
  transmute(`Contrast` = contrast, `Ratio` = round(ratio, 3),
            `SE` = round(SE, 3), `gl` = round(df, 1),
            `t` = round(t.ratio, 3), `p adjusted` = fmt_p(p.value))

emm_rt_r_wt <- emmeans(model_RT_red, pairwise ~ relatedness | word_type, type = "response")
tab_rt_c4   <- as.data.frame(summary(emm_rt_r_wt$contrasts)) %>%
  transmute(`Contrast` = contrast, `Word type` = word_type,
            `Ratio` = round(ratio, 3), `SE` = round(SE, 3),
            `gl` = round(df, 1), `t` = round(t.ratio, 3),
            `p adjusted` = fmt_p(p.value))

# -----------------------------------------------------------------------------
# 3. model_eeg_red — Gaussian GLMM (glmmTMB)
# -----------------------------------------------------------------------------

tidy_eeg <- tidy(model_eeg_red, effects = "fixed", component = "cond",
                 conf.int = TRUE, conf.level = 0.95)

tab_eeg_fixed <- tidy_eeg %>%
  transmute(
    `Predictor`    = term,
    `B`            = round(estimate, 3),
    `SE`           = round(std.error, 3),
    `95% CI lower` = round(conf.low, 3),
    `95% CI upper` = round(conf.high, 3),
    `z`            = round(statistic, 3),
    `p`            = fmt_p(p.value)
  )

tidy_eeg_re <- tidy(model_eeg_red, effects = "ran_pars", component = "cond")
tab_eeg_re <- tidy_eeg_re %>%
  transmute(`Group` = group, `Term` = term, `SD` = round(estimate, 3))

anova_eeg <- Anova(model_eeg_red, type = "II") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Effect") %>%
  rename(`Chi²` = Chisq, `gl` = Df, `p` = `Pr(>Chisq)`) %>%
  mutate(`Chi²` = round(`Chi²`, 3), `p` = fmt_p(`p`))

# Contrasts
emm_eeg_wt   <- emmeans(model_eeg_red, pairwise ~ word_type,               adjust = "tukey")
tab_eeg_c1   <- as.data.frame(summary(emm_eeg_wt$contrasts)) %>%
  transmute(`Contrast` = contrast, `Estimate` = round(estimate, 3),
            `SE` = round(SE, 3), `gl` = round(df, 1),
            `z` = round(z.ratio, 3), `p adjusted` = fmt_p(p.value))

emm_eeg_r    <- emmeans(model_eeg_red, pairwise ~ relatedness,             adjust = "tukey")
tab_eeg_c2   <- as.data.frame(summary(emm_eeg_r$contrasts)) %>%
  transmute(`Contrast` = contrast, `Estimate` = round(estimate, 3),
            `SE` = round(SE, 3), `gl` = round(df, 1),
            `z` = round(z.ratio, 3), `p adjusted` = fmt_p(p.value))

emm_eeg_wt_r <- emmeans(model_eeg_red, pairwise ~ word_type | relatedness, adjust = "tukey")
tab_eeg_c3   <- as.data.frame(summary(emm_eeg_wt_r$contrasts)) %>%
  transmute(`Contrast` = contrast, `Relatedness` = relatedness,
            `Estimate` = round(estimate, 3), `SE` = round(SE, 3),
            `gl` = round(df, 1), `z` = round(z.ratio, 3),
            `p adjusted` = fmt_p(p.value))

emm_eeg_r_wt <- emmeans(model_eeg_red, pairwise ~ relatedness | word_type, adjust = "tukey")
tab_eeg_c4   <- as.data.frame(summary(emm_eeg_r_wt$contrasts)) %>%
  transmute(`Contrast` = contrast, `Word type` = word_type,
            `Estimate` = round(estimate, 3), `SE` = round(SE, 3),
            `gl` = round(df, 1), `z` = round(z.ratio, 3),
            `p adjusted` = fmt_p(p.value))

# -----------------------------------------------------------------------------
# EXPORT TO WORD
# -----------------------------------------------------------------------------

doc <- read_docx()

add_table_to_doc <- function(doc, df, titulo) {
  doc <- doc %>%
    body_add_par(titulo, style = "heading 2") %>%
    body_add_flextable(make_ft(df, titulo)) %>%
    body_add_par("", style = "Normal")
  doc
}

# --- Model 1: Accuracy ---
doc <- doc %>% body_add_par("Model 1: Accuracy (GLMM binomial — reduced RE)", style = "heading 1")
doc <- add_table_to_doc(doc, tab_acc_fixed, "Table S1a. Fixed effects — Accuracy")
doc <- add_table_to_doc(doc, tab_acc_re,    "Table S1b. Random effects — Accuracy")
doc <- add_table_to_doc(doc, anova_acc,     "Table S1c. Type II ANOVA (Wald χ²) — Accuracy")
doc <- add_table_to_doc(doc, tab_acc_c1,    "Table S1d. Contrasts: word_type (OR) — Accuracy")
doc <- add_table_to_doc(doc, tab_acc_c2,    "Table S1e. Contrasts: word_type | relatedness (OR) — Accuracy")
doc <- add_table_to_doc(doc, tab_acc_c3,    "Table S1f. Contrasts: relatedness | word_type (OR) — Accuracy")

# --- Model 2: RT ---
doc <- doc %>% body_add_par("Model 2: Reaction Times (LMM, log-RT — reduced RE)", style = "heading 1")
doc <- add_table_to_doc(doc, tab_rt_fixed, "Table S2a. Fixed effects — RT")
doc <- add_table_to_doc(doc, tab_rt_re,    "Table S2b. Random effects — RT")
doc <- add_table_to_doc(doc, anova_rt,     "Table S2c. Type II ANOVA (Wald χ²) — RT")
doc <- add_table_to_doc(doc, tab_rt_c1,    "Table S2d. Contrasts: word_type (ratio) — RT")
doc <- add_table_to_doc(doc, tab_rt_c2,    "Table S2e. Contrasts: word_type | relatedness (ratio) — RT")
doc <- add_table_to_doc(doc, tab_rt_c3,    "Table S2f. Contrasts: relatedness (ratio) — RT")
doc <- add_table_to_doc(doc, tab_rt_c4,    "Table S2g. Contrasts: relatedness | word_type (ratio) — RT")

# --- Model 3: EEG ---
doc <- doc %>% body_add_par("Model 3: EEG Amplitudes (glmmTMB — reduced RE)", style = "heading 1")
doc <- add_table_to_doc(doc, tab_eeg_fixed, "Table S3a. Fixed effects — EEG")
doc <- add_table_to_doc(doc, tab_eeg_re,    "Table S3b. Random effects — EEG")
doc <- add_table_to_doc(doc, anova_eeg,     "Table S3c. Type II ANOVA (Wald χ²) — EEG")
doc <- add_table_to_doc(doc, tab_eeg_c1,    "Table S3d. Contrasts: word_type — EEG")
doc <- add_table_to_doc(doc, tab_eeg_c2,    "Table S3e. Contrasts: relatedness — EEG")
doc <- add_table_to_doc(doc, tab_eeg_c3,    "Table S3f. Contrasts: word_type | relatedness — EEG")
doc <- add_table_to_doc(doc, tab_eeg_c4,    "Table S3g. Contrasts: relatedness | word_type — EEG")

print(doc, target = "supplementary_tables_saturated.docx")
message("✓ File saved: supplementary_tables_saturated.docx")
