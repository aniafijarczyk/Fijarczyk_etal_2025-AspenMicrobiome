rm(list=ls())

library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(chemometrics)
library(MASS)
library(DAAG)
library(phyloseq)
library(microbiome)
library(cowplot)
library(nlme)
library(MuMIn)
library(car)
library(ggeffects)
library(usdm) # vif
library(relaimpo) # calc.relimp(()

library(performance) # for Nakagawa conditional/marginal R2
library(partR2) # for part R2 values

library(betareg) # betareg()
library(lmtest) # lrtest()
library(StepBeta)



#=======================#
#--------DATA-----------#
#=======================#


df <- read.csv("15_diversity_vs_soil_tables.tsv", sep="\t", header = T)
dim(df) # n=42
head(df)



#====================================#
#       Checking relationships       #
#====================================#


# Dataset with 2 outliers removed
df.sub <- df %>% filter(P < 300) %>% filter(C_N < 70)
hist(df.sub$C_N)
hist(df.sub$P)
dim(df.sub) # n=40


# Selecting a suset of uncorrelated variables
env.sel <- c("C_N","P","K","pH_H2O","Na")


### Checking distributions of diversity

### simpson
head(df.sub)
p3A <- ggplot(df.sub) + aes(x = diversity_gini_simpson_its) +
  geom_histogram(bins=12, fill="skyblue")
  #geom_density(aes(color=Country))

p3B <- ggplot(df.sub) + aes(x = diversity_gini_simpson_its) +
  geom_density(aes(color=Country))

p3C <- ggplot(df.sub) + aes(x = log(diversity_gini_simpson_its)) +
  geom_density(aes(color=Country))

plot_grid(p3A, p3B, p3C)

### shannon
p4A <- ggplot(df.sub) + aes(x = diversity_shannon_its) +
  geom_histogram(bins=12, fill="skyblue")
#geom_density(aes(color=Country))

p4B <- ggplot(df.sub) + aes(x = diversity_shannon_its) +
  geom_density(aes(color=Country))

png("17_diversity_vs_soil_run_histograms_diversity.png", w = 1200, h=1000, res=150)
plot_grid(p3A, p3B, p4A, p4B, ncol=2)
dev.off()


# Conclusions
# Distribution of most soil parameters  differs slighlty between groups
# Within groups their distributions are close to normal
# simpson - values between 0 and 1, in reality between 0.7 and 1, skewed, not really differentiated between countries
# shannon - close to normal distribution values > 0, slightly different between countries




#====================================#
#              Models                #
#====================================#

# Simpson - betareg() because response is proportional (between 0 and 1),
# including Country as a fixed effect - there are differences between countries in soil params, but not so much in simpson
# 5 + country -> too many values
# Check effects for each param separately, 


hist(df.sub$P, breaks=20)
dim(df.sub)

#df.sub <- df.sub %>% filter(P <100)

breg.1 <- betareg(diversity_gini_simpson_its ~ pH_H2O, df.sub)
breg.2 <- betareg(diversity_gini_simpson_its ~ C_N, df.sub)
breg.3 <- betareg(diversity_gini_simpson_its ~ P, df.sub)
breg.4 <- betareg(diversity_gini_simpson_its ~ K, df.sub)
breg.5 <- betareg(diversity_gini_simpson_its ~ Na, df.sub)
breg.6 <- betareg(diversity_gini_simpson_its ~ Country, df.sub)

AIC(breg.1, breg.2, breg.3, breg.4, breg.5, breg.6)

summary(breg.1)
summary(breg.2)
summary(breg.3)
summary(breg.4)
summary(breg.5)
summary(breg.6)

# Order = P, pH, Country, K, Na, C_N
# I'm dropping C_N



#====================================#
#       Constructing the model       #
#====================================#


### Full model + selection
full.model <- betareg(diversity_gini_simpson_its ~ K*P*Na*pH_H2O*Country, df.sub)
# Stepwise model selection
reduced.model <- StepBeta(full.model, k = 2, dispersion = T)
summary(reduced.model)

mod <- betareg(diversity_gini_simpson_its ~ P + P:Country, df.sub)
summary(mod)

### Checking the model visually
ggplot(df.sub) + aes(x = P, y = diversity_gini_simpson_its) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(color = Country))

ggplot(df.sub[df.sub$P<75,]) + aes(x = P, y = diversity_gini_simpson_its) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(color = Country))

ggplot(df.sub[df.sub$P<75,]) + aes(x = P, y = diversity_gini_simpson_its) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm")


#====================================#
#     Model assumptions and fit      #
#====================================#

# betareg
# designed for modeling continuous response variables bounded between 0 and 1, excluding 0 and 1
# variance depends on the mean — smaller or larger means often show less variation

mod <- betareg(diversity_gini_simpson_its ~ P + P:Country, df.sub)
summary(mod)

# Singularity
performance::check_singularity(mod)

# Distribution of residuals
hist(residuals(mod), breaks=25)
# No clear pattern in residuals vs fitted values
# Mean–variance relationship reasonable - Beta variance depends on mean
plot(fitted(mod), residuals(mod))
plot(mod)
# point with high leverage
lev <- hatvalues(mod)
summary(lev)
names(lev[lev == max(lev)])
hl <- df.sub[which(hatvalues(mod) > 0.4), ]
hl

# Cook's distance
plot(form=sqrt(abs(resid(.)))~fitted(.),mod, 2)
cooks <- cooks.distance(mod)
idx <- which(cooks > 0.2)
df.sub[idx, ] # one lowest pH and high C_N

# Checking the model without those values
df.sub[df.sub$P>100,]
df.sub.cd <- df.sub %>% filter(!short.code %in% c("FLOR1.5"))
full.model.cd <- betareg(diversity_gini_simpson_its ~ K*P*Na*pH_H2O*Country, df.sub.cd)
reduced.model.cd <- StepBeta(full.model.cd, k = 2, dispersion = T)
summary(reduced.model.cd)
mod.cd <- betareg(diversity_gini_simpson_its ~ pH_H2O, df.sub.cd)
hist(residuals(mod.cd), breaks=25)
plot(fitted(mod.cd), residuals(mod.cd))
plot(mod.cd)
# Doesnt look like a good idea, there are other outliers now
  

# Colinearity
mod.add <- betareg(diversity_gini_simpson_its ~ P + C_N + pH_H2O + Country, data = df.sub)
check_collinearity(mod.add)
# !!! P has very high upper CI - potential colinearity dependent on data draws or model parameters


# Q-Q - doesn’t need to be normal
qqnorm(resid(mod)); qqline(resid(mod))

# Check assumptions - doesn’t need to be normal
shap <- shapiro.test(residuals(mod, type = "pearson"))
print(paste0("Shapiro p: ",shap$p.value))



#====================================#
#           Model summary            #
#====================================#

summary(mod)
coefficients(mod)
confint(mod, level = 0.95)
r2(mod)
summary(mod)$pseudo.r.squared



#====================================#
####            END               ####
#====================================#
