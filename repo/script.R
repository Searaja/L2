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
  filter(window_start == 0.3, channel == "Cz")

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
