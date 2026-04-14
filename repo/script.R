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
library(nlme)

library(broom.mixed)   # tidy() para modelos mixtos
library(flextable)
library(officer)


# BEHAVIORAL DATA

df <- read.csv("data/behavioral.csv", header = TRUE) 

df$subject           <- factor(df$file)
df$word_type           <- factor(df$categ)
df$relatedness <- factor(df$corrAns)
df$accuracy <- factor(df$correct)
df$item <- factor(df$word_prime_merge)
levels(df$word_type)[levels(df$word_type) == 1] <- "L2-Remote"
levels(df$word_type)[levels(df$word_type) == 2] <- "L2-Recent"
levels(df$word_type)[levels(df$word_type) == 3] <- "L1"
levels(df$relatedness)[levels(df$relatedness) == "right"] <- "Related"
levels(df$relatedness)[levels(df$relatedness) == "left"] <- "Unrelated"
df$word_type <- factor(df$word_type,
                       levels = c("L1", "L2-Remote", "L2-Recent"))

# Accuracy
# GLMM

model_Acc <- glmer(
  accuracy ~ word_type * relatedness +
    (1 | subject) + (1 | item), data = df, family = binomial)

summary(model_Acc)

drop1(model_Acc, test="Chisq")

sim_res <- simulateResiduals(model_Acc)
plot(sim_res)

Anova(model_Acc, type="II")

modelsummary(model_Acc, statistic = "conf.int")

# Interactions

p <- emmip(model_Acc, relatedness ~ word_type)
print(p + ggplot2::ggtitle("Accuracy"))
emmeans(model_Acc, pairwise ~ word_type , type = "response")
emmeans(model_Acc, pairwise ~ word_type | relatedness, type = "response")
pairs(emmeans(model_Acc, ~ relatedness | word_type, type = "response"))

# RT
# LMM - log-transformed
# Only correct answers and < 2.5 SD from mean

k <- 2.5

df_filt <- df %>%
  group_by(corrAns, categ) %>%
  filter(abs(rt - mean(rt, na.rm = TRUE)) < k * sd(rt, na.rm = TRUE)) %>%
  ungroup()

table(df_filt$word_type, df_filt$relatedness)

model_RT <- lme( log(rt) ~ word_type * relatedness 
      ,random = ~1 | subject/item, data = df_filt, weights = varIdent(form = ~1 | word_type))

summary(model_RT)

res <- residuals(model_RT)
qqnorm(res)
qqline(res)

plot(fitted(model_RT), res)
abline(h = 0)


Anova(model_RT, type="II")

p <- emmip(model_RT, relatedness ~ word_type)
print(p + ggplot2::ggtitle("RT"))
emmeans(model_RT, pairwise ~ word_type, type="response")
emmeans(model_RT, pairwise ~ word_type|relatedness, type="response")
emmeans(model_RT, pairwise ~ relatedness, type="response")
emmeans(model_RT, pairwise ~ relatedness|word_type, type="response")

# EEG DATA

data <- read.csv("data/epoch_amplitudes.csv", header = TRUE)

data$subject <- as.factor(data$subject)
data$item <- as.factor(data$item)
data$relatedness <- as.factor(data$relatedness)
data$word_type <- as.factor(data$language)
data$channel <- as.factor(data$channel)

levels(data$word_type)[levels(data$word_type) == "rem"] <- "L2-Remote"
levels(data$word_type)[levels(data$word_type) == "rec"] <- "L2-Recent"
levels(data$word_type)[levels(data$word_type) == "esp"] <- "L1"
levels(data$relatedness)[levels(data$relatedness) == 1] <- "Related"
levels(data$relatedness)[levels(data$relatedness) == 0] <- "Unrelated"


# CHANGE window_start AND channel TO SEE RESULTS FOR DIFFERENT
# TIME WINDOWS AND ELECTRODE SITES:
# 0.3 FOR N400 (300-500 ms); 0.5 FOR LPC (500-800 ms)
# LMM
# Only correct answers and RT < 2.5 SD from RT mean

data_filt <- data %>%
  filter(window_start == 0.5, channel == "Pz")

data_filt <- data_filt %>%
  group_by(subject,word_type,relatedness) %>%
  mutate(
    media_rt = mean(RT, na.rm = TRUE),
    sd_rt = sd(RT, na.rm = TRUE)
  ) %>%
  filter(RT > media_rt - k*sd_rt,     
         RT < media_rt + k*sd_rt) %>% 
  ungroup()


report <- data_filt %>%
  group_by(subject, word_type, relatedness) %>%
  summarise(n_trials = n(), .groups = "drop") %>%
  group_by(word_type, relatedness) %>%
  summarise(
    mediana_trials = median(n_trials),
    q1_trials = quantile(n_trials, 0.25),
    q3_trials = quantile(n_trials, 0.75),
    iqr_trials = IQR(n_trials),
    min_trials = min(n_trials),
    max_trials = max(n_trials),
    total_sujetos = n(),
    .groups = "drop"
  )

print(report)

#ggplot(data_filt, aes(x = subject, y = mean_amplitude)) +
#  geom_boxplot() +
#  theme_minimal()

data_filt <- data_filt %>%
  select(mean_amplitude , word_type , relatedness, subject, item) %>%
  na.omit()

model_eeg <- glmmTMB(
  mean_amplitude ~ word_type * relatedness +
    (1 | subject) + (1 | item),
  dispformula = ~ word_type*relatedness,
  data = data_filt
)

summary(model_eeg)

Anova(model_eeg, type="II")

hist(data_filt$mean_amplitude, breaks = 100)

sim_res <- simulateResiduals(model_eeg)
plot(sim_res)
testUniformity(sim_res)

emm <- emmeans(model_eeg, ~ word_type * relatedness)
plot(emm)

emmeans(model_eeg, pairwise ~ word_type, adjust = "tukey")
emmeans(model_eeg, pairwise ~ relatedness, adjust = "tukey")

emmeans(model_eeg, pairwise ~ word_type|relatedness, adjust = "tukey")
emmeans(model_eeg, pairwise ~ relatedness|word_type, adjust = "tukey")

res <- residuals(model_eeg)

qqnorm(res)
qqline(res)

plot(fitted(model_eeg), res)
abline(h = 0)
# =============================================================================
# SUPPLEMENTARY TABLES
# Models: model_Acc (glmer), model_RT (lme), model_eeg (glmmTMB)
# Requires: broom.mixed, flextable, officer, dplyr, lmerTest, car
# =============================================================================

library(lme4)
library(lmerTest)
library(nlme)
library(glmmTMB)
library(car)
library(broom.mixed)   # tidy() para modelos mixtos
library(flextable)
library(officer)
library(dplyr)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# Round all numeric columns in a data.frame
round_df <- function(df, digits = 3) {
  df %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
}

# Format p-values: <.001, <.01, <.05, or exact value
fmt_p <- function(p) {
  case_when(
    p < .001 ~ "<.001",
    p < .01  ~ "<.01",
    p < .05  ~ "<.05",
    TRUE     ~ as.character(round(p, 3))
  )
}

# Build a flextable with APA-style formatting
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
# 1. model_Acc — Binomial GLMM (glmer)
# -----------------------------------------------------------------------------

## 1a. Fixed effects
tidy_acc <- tidy(model_Acc, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

tab_acc_fixed <- tidy_acc %>%
  transmute(
    `Predictor`        = term,
    `B`                = estimate,
    `SE`               = std.error,
    `95% CI lower`       = conf.low,
    `95% CI upper`       = conf.high,
    `z`                = statistic,
    `p`                = fmt_p(p.value)
  ) %>%
  round_df(3) %>%
  mutate(`p` = fmt_p(as.numeric(gsub("<", "", `p`))))   # conserva el formato

# re-apply fmt_p sobre el valor original (no el redondeado)
tab_acc_fixed$p <- fmt_p(tidy_acc$p.value)

## 1b. Random effects
tidy_acc_re <- tidy(model_Acc, effects = "ran_pars")

tab_acc_re <- tidy_acc_re %>%
  transmute(
    `Group`   = group,
    `Term` = term,
    `SD`      = round(estimate, 3)
  )

## 1c. Type II ANOVA (Wald chi-square)
anova_acc <- Anova(model_Acc, type = "II") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Effect") %>%
  rename(`Chi²` = Chisq, `gl` = Df, `p` = `Pr(>Chisq)`) %>%
  mutate(
    `Chi²` = round(`Chi²`, 3),
    `p`    = fmt_p(`p`)
  )

# -----------------------------------------------------------------------------
# 2. model_RT — LMM with nlme::lme (log RT)
# -----------------------------------------------------------------------------

## 2a. Fixed effects
sum_rt <- summary(model_RT)
fe_rt  <- as.data.frame(sum_rt$tTable)

tab_rt_fixed <- fe_rt %>%
  tibble::rownames_to_column("Predictor") %>%
  transmute(
    `Predictor` = Predictor,
    `B`         = round(Value, 3),
    `SE`        = round(Std.Error, 3),
    `gl`        = round(DF, 0),
    `t`         = round(`t-value`, 3),
    `p`         = fmt_p(`p-value`)
  )

# IC 95% para lme (intervals())
ic_rt <- tryCatch({
  intervals(model_RT, which = "fixed")$fixed %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Predictor") %>%
    transmute(Predictor, `95% CI lower` = round(lower, 3), `95% CI upper` = round(upper, 3))
}, error = function(e) NULL)

if (!is.null(ic_rt)) {
  tab_rt_fixed <- tab_rt_fixed %>%
    left_join(ic_rt, by = "Predictor") %>%
    select(Predictor, B, SE, `95% CI lower`, `95% CI upper`, gl, t, p)
}

## 2b. Random effects (variances)
vc_rt <- VarCorr(model_RT)
tab_rt_re <- as.data.frame(vc_rt) %>%
  transmute(
    `Group`   = grp,
    `Term` = var1,
    `SD`      = round(sdcor, 3)
  ) %>%
  filter(!is.na(`Term`))

## 2c. Type II ANOVA
# car::Anova with nlme::lme returns Chisq/Df/Pr(>Chisq)
anova_rt <- Anova(model_RT, type = "II") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Effect") %>%
  rename(`Chi²` = Chisq, `df` = Df, `p_raw` = `Pr(>Chisq)`) %>%
  mutate(
    `Chi²` = round(`Chi²`, 3),
    `p`      = fmt_p(p_raw)
  ) %>%
  select(Effect, `Chi²`, df, p)

# -----------------------------------------------------------------------------
# 3. model_eeg — Gaussian GLMM with glmmTMB
# -----------------------------------------------------------------------------

## 3a. Fixed effects (conditional component)
tidy_eeg <- tidy(model_eeg, effects = "fixed", component = "cond",
                 conf.int = TRUE, conf.level = 0.95)

tab_eeg_fixed <- tidy_eeg %>%
  transmute(
    `Predictor`  = term,
    `B`          = round(estimate, 3),
    `SE`         = round(std.error, 3),
    `95% CI lower` = round(conf.low, 3),
    `95% CI upper` = round(conf.high, 3),
    `z`          = round(statistic, 3),
    `p`          = fmt_p(p.value)
  )

## 3b. Random effects
tidy_eeg_re <- tidy(model_eeg, effects = "ran_pars", component = "cond")

tab_eeg_re <- tidy_eeg_re %>%
  transmute(
    `Group`   = group,
    `Term` = term,
    `SD`      = round(estimate, 3)
  )

## 3c. Type II ANOVA
anova_eeg <- Anova(model_eeg, type = "II") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Effect") %>%
  rename(`Chi²` = Chisq, `gl` = Df, `p` = `Pr(>Chisq)`) %>%
  mutate(
    `Chi²` = round(`Chi²`, 3),
    `p`    = fmt_p(`p`)
  )

# -----------------------------------------------------------------------------
# HELPER FUNCTION: extract emmeans contrasts as a clean data.frame
# -----------------------------------------------------------------------------

# Takes a pairwise emmeans object and returns a formatted data.frame
fmt_emmeans_contrasts <- function(emm_pairs, stat_type = "t") {
  df <- as.data.frame(summary(emm_pairs$contrasts))
  
  # Common columns
  out <- df %>%
    rename_with(~ gsub("\\.", "_", .x)) %>%
    rename_with(tolower)
  
  # Detect statistic column (t.ratio or z.ratio)
  stat_col <- grep("ratio", names(out), value = TRUE)[1]
  
  out <- out %>%
    transmute(
      `Contrast`    = contrast,
      `Estimate`     = round(estimate, 3),
      `SE`           = round(std_error, 3),
      !!stat_type   := round(.data[[stat_col]], 3),
      `gl`           = if ("df" %in% names(out)) round(df, 1) else NA_real_,
      `p adjusted`   = fmt_p(p_value)
    )
  
  # Add odds/response ratio column if present (binomial/response scale models)
  if ("odds_ratio" %in% names(df)) {
    out <- out %>% mutate(`OR` = round(df$odds_ratio, 3), .after = `Estimate`)
  }
  if ("ratio" %in% names(df)) {
    out <- out %>% mutate(`Ratio` = round(df$ratio, 3), .after = `Estimate`)
  }
  
  out
}

# Version for conditional contrasts (~ A | B): adds moderator column
fmt_emmeans_contrasts_by <- function(emm_pairs, by_var, stat_type = "t") {
  df <- as.data.frame(summary(emm_pairs$contrasts))
  
  stat_col <- grep("ratio", names(df), value = TRUE)[1]
  p_col    <- grep("p.value|p_value", names(df), value = TRUE)[1]
  se_col   <- grep("SE|std.error|std_error", names(df), value = TRUE, ignore.case = TRUE)[1]
  
  out <- df %>%
    transmute(
      `Contrast`    = contrast,
      !!by_var      := .data[[by_var]],
      `Estimate`     = round(estimate, 3),
      `SE`           = round(.data[[se_col]], 3),
      !!stat_type   := round(.data[[stat_col]], 3),
      `gl`           = if ("df" %in% names(df)) round(df, 1) else NA_real_,
      `p adjusted`   = fmt_p(.data[[p_col]])
    )
  out
}

# -----------------------------------------------------------------------------
# CONTRASTS — Model 1: Accuracy (glmer, response scale = odds ratio)
# -----------------------------------------------------------------------------

# word_type (marginal)
emm_acc_wt    <- emmeans(model_Acc, pairwise ~ word_type, type = "response")
tab_acc_c1    <- as.data.frame(summary(emm_acc_wt$contrasts)) %>%
  transmute(
    `Contrast`  = contrast,
    `OR`         = round(odds.ratio, 3),
    `SE`         = round(SE, 3),
    `z`          = round(z.ratio, 3),
    `p adjusted` = fmt_p(p.value)
  )

# word_type | relatedness
emm_acc_wt_r  <- emmeans(model_Acc, pairwise ~ word_type | relatedness, type = "response")
tab_acc_c2    <- as.data.frame(summary(emm_acc_wt_r$contrasts)) %>%
  transmute(
    `Contrast`    = contrast,
    `Relatedness`  = relatedness,
    `OR`           = round(odds.ratio, 3),
    `SE`           = round(SE, 3),
    `z`            = round(z.ratio, 3),
    `p adjusted`   = fmt_p(p.value)
  )

# relatedness | word_type
emm_acc_r_wt  <- emmeans(model_Acc, pairwise ~ relatedness | word_type, type = "response")
tab_acc_c3    <- as.data.frame(summary(emm_acc_r_wt$contrasts)) %>%
  transmute(
    `Contrast`    = contrast,
    `Word type`    = word_type,
    `OR`           = round(odds.ratio, 3),
    `SE`           = round(SE, 3),
    `z`            = round(z.ratio, 3),
    `p adjusted`   = fmt_p(p.value)
  )

# -----------------------------------------------------------------------------
# CONTRASTS — Model 2: RT (lme, response scale = ratio of geometric means)
# -----------------------------------------------------------------------------

# word_type (marginal)
emm_rt_wt     <- emmeans(model_RT, pairwise ~ word_type, type = "response")
tab_rt_c1     <- as.data.frame(summary(emm_rt_wt$contrasts)) %>%
  transmute(
    `Contrast`  = contrast,
    `Ratio`      = round(ratio, 3),
    `SE`         = round(SE, 3),
    `gl`         = round(df, 1),
    `t`          = round(t.ratio, 3),
    `p adjusted` = fmt_p(p.value)
  )

# word_type | relatedness
emm_rt_wt_r   <- emmeans(model_RT, pairwise ~ word_type | relatedness, type = "response")
tab_rt_c2     <- as.data.frame(summary(emm_rt_wt_r$contrasts)) %>%
  transmute(
    `Contrast`    = contrast,
    `Relatedness`  = relatedness,
    `Ratio`        = round(ratio, 3),
    `SE`           = round(SE, 3),
    `gl`           = round(df, 1),
    `t`            = round(t.ratio, 3),
    `p adjusted`   = fmt_p(p.value)
  )

# relatedness (marginal)
emm_rt_r      <- emmeans(model_RT, pairwise ~ relatedness, type = "response")
tab_rt_c3     <- as.data.frame(summary(emm_rt_r$contrasts)) %>%
  transmute(
    `Contrast`  = contrast,
    `Ratio`      = round(ratio, 3),
    `SE`         = round(SE, 3),
    `gl`         = round(df, 1),
    `t`          = round(t.ratio, 3),
    `p adjusted` = fmt_p(p.value)
  )

# relatedness | word_type
emm_rt_r_wt   <- emmeans(model_RT, pairwise ~ relatedness | word_type, type = "response")
tab_rt_c4     <- as.data.frame(summary(emm_rt_r_wt$contrasts)) %>%
  transmute(
    `Contrast`    = contrast,
    `Word type`    = word_type,
    `Ratio`        = round(ratio, 3),
    `SE`           = round(SE, 3),
    `gl`           = round(df, 1),
    `t`            = round(t.ratio, 3),
    `p adjusted`   = fmt_p(p.value)
  )

# -----------------------------------------------------------------------------
# CONTRASTS — Model 3: EEG (glmmTMB, original scale = mean difference)
# -----------------------------------------------------------------------------

# word_type (marginal)
emm_eeg_wt    <- emmeans(model_eeg, pairwise ~ word_type, adjust = "tukey")
tab_eeg_c1    <- as.data.frame(summary(emm_eeg_wt$contrasts)) %>%
  transmute(
    `Contrast`  = contrast,
    `Estimate`   = round(estimate, 3),
    `SE`         = round(SE, 3),
    `gl`         = round(df, 1),
    `z`          = round(z.ratio, 3),
    `p adjusted` = fmt_p(p.value)
  )

# relatedness (marginal)
emm_eeg_r     <- emmeans(model_eeg, pairwise ~ relatedness, adjust = "tukey")
tab_eeg_c2    <- as.data.frame(summary(emm_eeg_r$contrasts)) %>%
  transmute(
    `Contrast`  = contrast,
    `Estimate`   = round(estimate, 3),
    `SE`         = round(SE, 3),
    `gl`         = round(df, 1),
    `z`          = round(z.ratio, 3),
    `p adjusted` = fmt_p(p.value)
  )

# word_type | relatedness
emm_eeg_wt_r  <- emmeans(model_eeg, pairwise ~ word_type | relatedness, adjust = "tukey")
tab_eeg_c3    <- as.data.frame(summary(emm_eeg_wt_r$contrasts)) %>%
  transmute(
    `Contrast`    = contrast,
    `Relatedness`  = relatedness,
    `Estimate`     = round(estimate, 3),
    `SE`           = round(SE, 3),
    `gl`           = round(df, 1),
    `z`            = round(z.ratio, 3),
    `p adjusted`   = fmt_p(p.value)
  )

# relatedness | word_type
emm_eeg_r_wt  <- emmeans(model_eeg, pairwise ~ relatedness | word_type, adjust = "tukey")
tab_eeg_c4    <- as.data.frame(summary(emm_eeg_r_wt$contrasts)) %>%
  transmute(
    `Contrast`    = contrast,
    `Word type`    = word_type,
    `Estimate`     = round(estimate, 3),
    `SE`           = round(SE, 3),
    `gl`           = round(df, 1),
    `z`            = round(z.ratio, 3),
    `p adjusted`   = fmt_p(p.value)
  )

# -----------------------------------------------------------------------------
# EXPORT TO WORD
# -----------------------------------------------------------------------------

doc <- read_docx()

# Function to add a table + heading to the doc
add_table_to_doc <- function(doc, df, titulo) {
  doc <- doc %>%
    body_add_par(titulo, style = "heading 2") %>%
    body_add_flextable(make_ft(df, titulo)) %>%
    body_add_par("", style = "Normal")   # blank line
  doc
}

# --- Model 1: Accuracy (GLMM binomial) ---
doc <- doc %>% body_add_par("Model 1: Accuracy (GLMM binomial)", style = "heading 1")
doc <- add_table_to_doc(doc, tab_acc_fixed, "Table S1a. Fixed effects — Accuracy")
doc <- add_table_to_doc(doc, tab_acc_re,    "Table S1b. Random effects — Accuracy")
doc <- add_table_to_doc(doc, anova_acc,     "Table S1c. Type II ANOVA (Wald χ²) — Accuracy")
doc <- add_table_to_doc(doc, tab_acc_c1,    "Table S1d. Contrasts: word_type (OR) — Accuracy")
doc <- add_table_to_doc(doc, tab_acc_c2,    "Table S1e. Contrasts: word_type | relatedness (OR) — Accuracy")
doc <- add_table_to_doc(doc, tab_acc_c3,    "Table S1f. Contrasts: relatedness | word_type (OR) — Accuracy")

# --- Model 2: Reaction times (LMM, log-RT) ---
doc <- doc %>% body_add_par("Model 2: Reaction Times (LMM, log-RT)", style = "heading 1")
doc <- add_table_to_doc(doc, tab_rt_fixed, "Table S2a. Fixed effects — RT")
doc <- add_table_to_doc(doc, tab_rt_re,   "Table S2b. Random effects — RT")
doc <- add_table_to_doc(doc, anova_rt,    "Table S2c. Type II ANOVA (Wald χ²) — RT")
doc <- add_table_to_doc(doc, tab_rt_c1,   "Table S2d. Contrasts: word_type (ratio) — RT")
doc <- add_table_to_doc(doc, tab_rt_c2,   "Table S2e. Contrasts: word_type | relatedness (ratio) — RT")
doc <- add_table_to_doc(doc, tab_rt_c3,   "Table S2f. Contrasts: relatedness (ratio) — RT")
doc <- add_table_to_doc(doc, tab_rt_c4,   "Table S2g. Contrasts: relatedness | word_type (ratio) — RT")

# --- Model 3: EEG amplitudes (glmmTMB) ---
doc <- doc %>% body_add_par("Model 3: EEG Amplitudes (glmmTMB)", style = "heading 1")
doc <- add_table_to_doc(doc, tab_eeg_fixed, "Table S3a. Fixed effects — EEG")
doc <- add_table_to_doc(doc, tab_eeg_re,   "Table S3b. Random effects — EEG")
doc <- add_table_to_doc(doc, anova_eeg,    "Table S3c. Type II ANOVA (Wald χ²) — EEG")
doc <- add_table_to_doc(doc, tab_eeg_c1,   "Table S3d. Contrasts: word_type — EEG")
doc <- add_table_to_doc(doc, tab_eeg_c2,   "Table S3e. Contrasts: relatedness — EEG")
doc <- add_table_to_doc(doc, tab_eeg_c3,   "Table S3f. Contrasts: word_type | relatedness — EEG")
doc <- add_table_to_doc(doc, tab_eeg_c4,   "Table S3g. Contrasts: relatedness | word_type — EEG")

# Save
print(doc, target = "supplementary_tables.docx")
message("✓ File saved: supplementary_tables.docx")