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

df <- read.csv("behavioral.csv", header = TRUE) 

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
# 1

df_filtrado <- df#  %>% filter(word_type != "L2-Recent")

# Acc
# Modelo lineal mixto generalizado con interacción

model_Acc <- glmer(
  accuracy ~ word_type * relatedness +
    (1 | subject) + (1 | item), data = df_filtrado, family = binomial)

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
pairs(emmeans(model_Acc, ~ relatedness | word_type))

# RT
# Modelo lineal mixto con interacción. Solo correctas y con SD 2.5

df_filtrado2 <- df_filtrado %>%
  group_by(corrAns, categ) %>%
  filter(abs(rt - mean(rt, na.rm = TRUE)) < 2.5 * sd(rt, na.rm = TRUE)) %>%
  ungroup()

model_RT <- lme( log(rt) ~ word_type * relatedness 
      ,random = ~1 | subject/item, data = df_filtrado2, weights = varIdent(form = ~1 | word_type))

summary(model_RT)

Anova(model_RT, type="II")

#qqnorm(resid(model_RT, type = "normalized"))
#qqline(resid(model_RT, type = "normalized"))

#p <- emmip(model_RT,  word_type ~ relatedness,lmerTest.limit = 5726,pbkrtest.limit = 5726)
#print(p + ggplot2::ggtitle("RT"))
emmeans(model_RT, pairwise ~ word_type, type="response")
emmeans(model_RT, pairwise ~ relatedness, type="response")


datos <- read.csv("epoch_amplitudes.csv", header = TRUE)

datos$subject <- as.factor(datos$subject)
datos$item <- as.factor(datos$item)
datos$relatedness <- as.factor(datos$relatedness)
datos$word_type <- as.factor(datos$language)
datos$channel <- as.factor(datos$channel)

datos <- datos %>%
  filter(window_start == 0.3, channel == "Cz")

datos_filtrados <- datos %>%
  filter(subject!="CB", subject!="DK", subject!="EDA", subject!="MGM", subject!="RAC")

datos_filtrados <- datos_filtrados %>%
  group_by(subject,word_type,relatedness) %>%
  mutate(
    media_rt = mean(RT, na.rm = TRUE),
    sd_rt = sd(RT, na.rm = TRUE)
  ) %>%
  filter(RT > media_rt - 2.5*sd_rt,     # if SD filter changed here,
         RT < media_rt + 2.5*sd_rt) %>% # change in behavior too!
  ungroup()

resumen_por_sujeto <- datos_filtrados %>%
  group_by(subject, word_type, relatedness) %>%
  summarise(
    n_trials = n(), # Cuántos quedaron para este sujeto en esta condición
    .groups = "drop"
  )

# Ahora calculamos la MEDIA y el SD de esos conteos
reporte_final <- resumen_por_sujeto %>%
  group_by(word_type, relatedness) %>%
  summarise(
    media_trials = mean(n_trials),
    sd_trials = sd(n_trials),
    min_trials = min(n_trials), # Útil para ver el "peor" caso
    max_trials = max(n_trials),
    n_sujetos = n(),
    .groups = "drop"
  )

print(reporte_final)

reporte_robusto <- datos_filtrados %>%
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

print(reporte_robusto)

ggplot(resumen_por_sujeto, aes(x = factor(subject), y = n_trials, fill = relatedness)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~word_type) +
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
  select(mean_amplitude , word_type , relatedness, subject, item) %>%
  na.omit()

modelo <- glmmTMB(
  mean_amplitude ~ word_type * relatedness +
    (1 | subject) + (1 | item),
  dispformula = ~ word_type*relatedness,
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

emm <- emmeans(modelo, ~ word_type * relatedness)
emm

pairs(emm)

plot(emm)

# 1. Calcular las medias por grupo
emmeans(modelo, pairwise ~ word_type, adjust = "tukey")
emmeans(modelo, pairwise ~ relatedness, adjust = "tukey")

emmeans(modelo, pairwise ~ word_type|relatedness, adjust = "tukey")
emmeans(modelo, pairwise ~ relatedness|word_type, adjust = "tukey")

# 2. Obtienes las comparaciones de a pares
#comparaciones <- pairs(medias_ajustadas, adjust = "bonferroni")
#print(comparaciones)

res <- residuals(modelo)

qqnorm(res)
qqline(res)
# shapiro.test(res)

plot(fitted(modelo), res)
abline(h = 0)
