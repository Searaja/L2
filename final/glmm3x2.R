# Read
df <- read.csv("df.csv", header = TRUE) 

library(lme4)
library(DHARMa)
library(lmerTest)
library(dplyr)
library(emmeans)
library(nlme)
library(modelsummary)
library(car)

# Convertir factores si aún no están como factor
df$subject           <- factor(df$file)
df$condition           <- factor(df$categ)
df$relatedness <- factor(df$corrAns)
df$accuracy <- factor(df$correct)
df$item <- factor(df$word_prime_merge)
levels(df$condition)[levels(df$condition) == 1] <- "L2-Remote"
levels(df$condition)[levels(df$condition) == 2] <- "L2-Recent"
levels(df$condition)[levels(df$condition) == 3] <- "L1"
levels(df$relatedness)[levels(df$relatedness) == "right"] <- "Related"
levels(df$relatedness)[levels(df$relatedness) == "left"] <- "Unrelated"
df$condition <- factor(df$condition,
                       levels = c("L1", "L2-Remote", "L2-Recent"))
# 1

df_filtrado <- df#  %>% filter(condition != "L2-Recent")

# Acc
# Modelo lineal mixto generalizado con interacción

model_Acc <- glmer(
  accuracy ~ condition * relatedness +
    (1 | subject) + (1 | item), data = df_filtrado, family = binomial)

summary(model_Acc)

drop1(model_Acc, test="Chisq")

sim_res <- simulateResiduals(model_Acc)
plot(sim_res)

Anova(model_Acc, type="II")

modelsummary(model_Acc, statistic = "conf.int")

# Interactions

p <- emmip(model_Acc, relatedness ~ condition)
print(p + ggplot2::ggtitle("Accuracy"))
emmeans(model_Acc, pairwise ~ condition , type = "response")
emmeans(model_Acc, pairwise ~ condition | relatedness, type = "response")
pairs(emmeans(model_Acc, ~ relatedness | condition))

# RT
# Modelo lineal mixto con interacción. Solo correctas y con SD 2.5

df_filtrado2 <- df_filtrado %>%
  group_by(corrAns, categ) %>%
  filter(abs(rt - mean(rt, na.rm = TRUE)) < 2.5 * sd(rt, na.rm = TRUE)) %>%
  ungroup()

model_RT <- lme( log(rt) ~ condition * relatedness 
      ,random = ~1 | subject/item, data = df_filtrado2, weights = varIdent(form = ~1 | condition))

summary(model_RT)

Anova(model_RT, type="II")

#qqnorm(resid(model_RT, type = "normalized"))
#qqline(resid(model_RT, type = "normalized"))

#p <- emmip(model_RT,  condition ~ relatedness,lmerTest.limit = 5726,pbkrtest.limit = 5726)
#print(p + ggplot2::ggtitle("RT"))
emmeans(model_RT, pairwise ~ condition, type="response")
emmeans(model_RT, pairwise ~ relatedness, type="response")
