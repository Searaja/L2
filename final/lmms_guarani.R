library(lme4)
library(lmerTest)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(performance)
library(emmeans)
library(ggplot2)
library(influence.ME)
library(brms)
library(car)
library(ggplot2)

datos <- read.csv("epoch_amplitudes_22_avg3.csv", header = TRUE)

datos$subject <- as.factor(datos$subject)
datos$item <- as.factor(datos$item)
datos$relatedness <- as.factor(datos$relatedness)
datos$language <- as.factor(datos$language)
datos$channel <- as.factor(datos$channel)

datos <- datos %>%
  filter(window_start == 0.3, channel == "Cz")

datos_filtrados <- datos %>%
  filter(subject!="CB", subject!="DK", subject!="EDA", subject!="MGM", subject!="RAC")

datos_filtrados <- datos_filtrados %>%
  group_by(subject,language,relatedness) %>%
  mutate(
    media_rt = mean(RT, na.rm = TRUE),
    sd_rt = sd(RT, na.rm = TRUE)
  ) %>%
  filter(RT > media_rt - 2.5*sd_rt,     # if SD filter changed here,
         RT < media_rt + 2.5*sd_rt) %>% # change in behavior too!
  ungroup()

resumen_por_sujeto <- datos_filtrados %>%
  group_by(subject, language, relatedness) %>%
  summarise(
    n_trials = n(), # Cuántos quedaron para este sujeto en esta condición
    .groups = "drop"
  )

# Ahora calculamos la MEDIA y el SD de esos conteos
reporte_final <- resumen_por_sujeto %>%
  group_by(language, relatedness) %>%
  summarise(
    media_trials = mean(n_trials),
    sd_trials = sd(n_trials),
    min_trials = min(n_trials), # Útil para ver el "peor" caso
    max_trials = max(n_trials),
    n_sujetos = n(),
    .groups = "drop"
  )

print(reporte_final)

library(dplyr)

reporte_robusto <- datos_filtrados %>%
  group_by(subject, language, relatedness) %>%
  summarise(n_trials = n(), .groups = "drop") %>%
  group_by(language, relatedness) %>%
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

print(reporte_robusto)
library(ggplot2)

ggplot(resumen_por_sujeto, aes(x = factor(subject), y = n_trials, fill = relatedness)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~language) +
  theme_minimal() +
  labs(title = "Trials restantes por sujeto y condición",
       x = "Sujeto", y = "Número de Trials") +
  theme(axis.text.x = element_text(angle = 45))

#datos_filtrados <- datos_filtrados %>%
#  group_by(subject) %>%
#  mutate(
#    media_ma = mean(mean_amplitude, na.rm = TRUE),
#    sd_ma = sd(mean_amplitude, na.rm = TRUE)
#  ) %>%
#  filter(mean_amplitude > media_ma - 5*sd_ma,
#         mean_amplitude < media_ma + 5*sd_ma) %>%
#  ungroup()

ggplot(datos_filtrados, aes(x = subject, y = mean_amplitude)) +
  geom_boxplot() +
  theme_minimal()

datos_filtrados <- datos_filtrados %>%
  select(mean_amplitude , language , relatedness, subject, item) %>%
  na.omit()

modelo <- glmmTMB(
  mean_amplitude ~ language * relatedness +
    (1 | subject) + (1 | item),
  dispformula = ~ language*relatedness,
  data = datos_filtrados
)

#infl <- influence(modelo, group = "subject")

#plot(infl, which = "cook")

summary(modelo)

Anova(modelo, type="III")

hist(datos_filtrados$mean_amplitude, breaks = 100)

sim_res <- simulateResiduals(modelo)
plot(sim_res)
testUniformity(sim_res)
#check_model(modelo)

emm <- emmeans(modelo, ~ language * relatedness)
emm

pairs(emm)

plot(emm)

# 1. Calcular las medias por grupo
emmeans(modelo, pairwise ~ language, adjust = "tukey")
emmeans(modelo, pairwise ~ relatedness, adjust = "tukey")

emmeans(modelo, pairwise ~ language|relatedness, adjust = "tukey")
emmeans(modelo, pairwise ~ relatedness|language, adjust = "tukey")

# 2. Obtienes las comparaciones de a pares
#comparaciones <- pairs(medias_ajustadas, adjust = "bonferroni")
#print(comparaciones)

res <- residuals(modelo)

qqnorm(res)
qqline(res)
# shapiro.test(res)

plot(fitted(modelo), res)
abline(h = 0)




