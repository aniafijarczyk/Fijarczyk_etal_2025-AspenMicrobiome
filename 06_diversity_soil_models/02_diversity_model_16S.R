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

df$region <- factor(df$region, levels = c("Boreal","Cold_temperate","Warm_temperate"))

ggplot(df) + aes(y = short.code, x =Year, color = region) +
  geom_point()

ggplot(df) + aes(y = short.code, x =MAT, color = region) +
  geom_point()


ggplot(df) + aes(y = short.code, x =diversity_gini_simpson, color = region, shape=as.factor(Year)) +
  geom_point()


ggplot(df) + aes(y = short.code, x =diversity_gini_simpson_its, color = region, shape=as.factor(Year)) +
  geom_point()

df$Site <- factor(df$Site, levels = c("AMOS","STFE", "ESSI","FORE", "STET","FLOR1","FP2","Santiago"))
ggplot(df) + aes(x = Site, y =diversity_gini_simpson_its, color = region) +
  geom_point(aes(shape=as.factor(Year))) +
  geom_boxplot(fill=NA)
ggplot(df) + aes(x = Site, y =diversity_shannon_its, color = region) +
  geom_point(aes(shape=as.factor(Year))) +
  geom_boxplot(fill=NA)

ggplot(df) + aes(x = Site, y =diversity_gini_simpson, color = region) +
  geom_point(aes(shape=as.factor(Year))) +
  geom_boxplot(fill=NA)
ggplot(df) + aes(x = Site, y =diversity_shannon, color = region) +
  geom_point(aes(shape=as.factor(Year))) +
  geom_boxplot(fill=NA)



#====================================#
#       Checking relationships       #
#====================================#


library(psych)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- corr.test(x,y, method = "spearman")["r"]
  txt <- format(c(r, 0.123456789), digits = 3, nsmall = 3)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor))
    cex <- 1/strwidth(txt)
  
  test <- corr.test(x,y, method = "spearman")
  Signif <- symnum(test$p, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", " "))
  
  text(0.5, 0.5, txt, cex = cex)
  text(.8, .8, Signif, cex = cex, col = 2)
}

colnames(df)
env.vars <- c("C_N","P","K","CEC","pH_H2O","Mg","Ca","Mn","Al","Fe","Na","Sable","Limon","Argile","MAT","PPT_WS","TD")

df.env <- df[env.vars]
head(df.env)

#png("16_diversity_vs_soil_run_correlations.png", w = 1500, h=1500, res=150)
pairs(df.env, lower.panel = panel.smooth, upper.panel = panel.cor)
#dev.off()


# Conclusion -> C_N and P have strong outliers that can drive some correlations
hist(df$C_N)
hist(df$P)

# Dataset with 2 outliers removed
df.sub <- df %>% filter(P < 300) %>% filter(C_N < 70)
hist(df.sub$C_N)
hist(df.sub$P)
dim(df.sub) # n=40

#png("16_diversity_vs_soil_run_correlations_filtered.png", w = 1500, h=1500, res=150)
df.env <- df.sub[env.vars]
pairs(df.env, lower.panel = panel.smooth, upper.panel = panel.cor)
#dev.off()

# Selecting a suset of uncorrelated variables
env.sel <- c("C_N","P","K","pH_H2O","Na")
df.env <- df.sub[env.sel]
pairs(df.env, lower.panel = panel.smooth, upper.panel = panel.cor)

#png("16_diversity_vs_soil_run_correlations_filtered_subset.png", w = 1500, h=1500, res=150)
df.env <- df.sub[env.sel]
pairs(df.env, lower.panel = panel.smooth, upper.panel = panel.cor)
#dev.off()

### Checking distributions of soil parameters
head(df.sub)
dg <- df.sub %>% dplyr::select(all_of(c(c("short.code","Site", "Country", "region"), env.sel))) %>% gather(key = "Parameter", value = "Value", all_of(env.sel))
dg$log_Value <- log(dg$Value)
head(dg)

p1 <- ggplot(dg) + aes(x = Value) +
  geom_histogram(bins=12, fill="skyblue") +
  #geom_density(aes(color=Country)) +
  facet_wrap(.~Parameter, scales="free", nrow=1)

p2 <- ggplot(dg) + aes(x = Value) +
  #geom_histogram(bins=15, fill="skyblue") +
  geom_density(aes(color=Country)) +
  facet_wrap(.~Parameter, scales="free", nrow=1)


#png("16_diversity_vs_soil_run_histograms.png", w = 1800, h=1200, res=150)
plot_grid(p1, p2, ncol=1)
#dev.off()


p1B <- ggplot(dg) + aes(x = log_Value) +
  geom_histogram(bins=12, fill="skyblue") +
  #geom_density(aes(color=Country)) +
  facet_wrap(.~Parameter, scales="free", nrow=1)

p2B <- ggplot(dg) + aes(x = log_Value) +
  #geom_histogram(bins=15, fill="skyblue") +
  geom_density(aes(color=Country)) +
  facet_wrap(.~Parameter, scales="free", nrow=1)


#png("16_diversity_vs_soil_run_histograms_logs.png", w = 1800, h=1200, res=150)
plot_grid(p1B, p2B, ncol=1)
#dev.off()





#====================================#
#  Checking diversity distribution   #
#====================================#

### simpson
head(df.sub)
p3A <- ggplot(df.sub) + aes(x = diversity_gini_simpson) +
  geom_histogram(bins=12, fill="skyblue")
  #geom_density(aes(color=Country))

p3B <- ggplot(df.sub) + aes(x = diversity_gini_simpson) +
  geom_density(aes(color=Country))

p3C <- ggplot(df.sub) + aes(x = log(diversity_gini_simpson)) +
  geom_density(aes(color=Country))

plot_grid(p3B, p3C)

### shannon
p4A <- ggplot(df.sub) + aes(x = diversity_shannon) +
  geom_histogram(bins=12, fill="skyblue")
#geom_density(aes(color=Country))

p4B <- ggplot(df.sub) + aes(x = diversity_shannon) +
  geom_density(aes(color=Country))

png("16_diversity_vs_soil_run_histograms_diversity.png", w = 1200, h=1000, res=150)
plot_grid(p3A, p3B, p4A, p4B, ncol=2)
dev.off()


# Conclusions
# Distribution of most soil parameters is differs slighlty between groups
# Within groups their distributions are close to normal
# simpson - values between 0 and 1, in reality between 0.98 and 1, skewed, more differentiated between countries
# shannon - close to normal distribution values > 0, similar for the two countries




#====================================#
#              Models                #
#====================================#

# Simpson - betareg() because response is proportional (between 0 and 1),
# including Country as a fixed effect because there are differences between countries in simpson and in soil params
# 5 + country -> too many value
# Check effects for each param separately, 


breg.1 <- betareg(diversity_gini_simpson ~ pH_H2O, df.sub)
breg.2 <- betareg(diversity_gini_simpson ~ C_N, df.sub)
breg.3 <- betareg(diversity_gini_simpson ~ P, df.sub)
breg.4 <- betareg(diversity_gini_simpson ~ K, df.sub)
breg.5 <- betareg(diversity_gini_simpson ~ Na, df.sub)
breg.6 <- betareg(diversity_gini_simpson ~ Country, df.sub)

AIC(breg.1, breg.2, breg.3, breg.4, breg.5, breg.6)

summary(breg.1)
summary(breg.2)
summary(breg.3)
summary(breg.4)
summary(breg.5)
summary(breg.6)


# Order = Country, K, C_N, Na, pH, P
# I'm dropping P, P is slightly correlated with K and has lowest effect

#====================================#
#       Constructing the model       #
#====================================#

### Full model + selection
full.model <- betareg(diversity_gini_simpson ~ K*C_N*Na*pH_H2O*Country, df.sub)
# Stepwise model selection
reduced.model <- StepBeta(full.model, k = 2, dispersion = T)
summary(reduced.model)

mod <- betareg(diversity_gini_simpson ~ Country + Na + C_N:Country + pH_H2O:Country + C_N:pH_H2O:Country, df.sub)
summary(mod)
coefficients(mod)

### Checking the model visually
ggplot(df.sub) + aes(x = Na + C_N + pH_H2O + C_N*pH_H2O, y = diversity_gini_simpson) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(group = Country))
ggplot(df.sub) + aes(x = Na, y = diversity_gini_simpson) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm")+
  theme(panel.background = element_rect(fill=NA))
ggplot(df.sub) + aes(x = C_N, y = diversity_gini_simpson) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(group = Country))+
  theme(panel.background = element_rect(fill=NA))
ggplot(df.sub) + aes(x = pH_H2O, y = diversity_gini_simpson) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(group = Country))+
  theme(panel.background = element_rect(fill=NA))
ggplot(df.sub) + aes(x = pH_H2O*C_N, y = diversity_gini_simpson) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(group = Country))+
  theme(panel.background = element_rect(fill=NA))


#====================================#
#     Model assumptions and fit      #
#====================================#

# betareg
# designed for modeling continuous response variables bounded between 0 and 1, excluding 0 and 1
# variance depends on the mean — smaller or larger means often show less variation


mod <- betareg(diversity_gini_simpson ~ Country + Na + C_N:Country + pH_H2O:Country + C_N:pH_H2O:Country, df.sub)
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

# checking if model on log values works better
full.model.log <- betareg(diversity_gini_simpson ~ log(K)*C_N*log(Na)*pH_H2O*Country, df.sub)
reduced.model.log <- StepBeta(full.model.log, k = 2, dispersion = T)
summary(reduced.model.log)
hist(residuals(reduced.model.log), breaks=25)
plot(fitted(reduced.model.log), residuals(reduced.model.log))
plot(reduced.model.log)
# Cook's distance
plot(form=sqrt(abs(resid(.)))~fitted(.),mod, 2)
cooks <- cooks.distance(mod)
idx <- which(cooks > 0.2)
df.sub[idx, ] # one lowest pH and high C_N
summary(df.sub$Na)
summary(df.sub$pH_H2O)
hist(df.sub$pH_H2O, breaks=12)
hist(df.sub$C_N, breaks=12)
hist(df.sub$Na, breaks=12)
summary(df.sub$C_N)
ggplot(df.sub) + aes(x = pH_H2O, y = C_N) +
  geom_point() + geom_smooth()
ggplot(df.sub) + aes(x = pH_H2O, y = C_N) +
  geom_point() + geom_smooth(method="lm")


# Checking the model without those values
#df.sub.cd <- df.sub %>% filter(!short.code %in% c("ESSI.001.S","FORE.003S"))
df.sub.cd <- df.sub %>% filter(!short.code %in% c("ESSI.001.S"))
full.model.cd <- betareg(diversity_gini_simpson ~ K*C_N*Na*pH_H2O*Country, df.sub.cd)
reduced.model.cd <- StepBeta(full.model.cd, k = 2, dispersion = T)
summary(reduced.model.cd)
mod.cd <- betareg(diversity_gini_simpson ~ Country + Na + C_N:Country + pH_H2O:Country + C_N:pH_H2O, df.sub.cd)
hist(residuals(mod.cd), breaks=25)
plot(fitted(mod.cd), residuals(mod.cd))
plot(mod.cd)
# Doesnt look like a good idea, there are other outliers now
  

# Colinearity
mod.add <- betareg(diversity_gini_simpson ~ Country + Na + C_N + pH_H2O, data = df.sub)
check_collinearity(mod.add)

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



### Pseudo R2 of the interaction terms
rmod <- betareg(diversity_gini_simpson ~ Country + Na + C_N:Country + pH_H2O:Country + C_N:pH_H2O:Country, df.sub)

rmod.country <- betareg(diversity_gini_simpson ~ Na + C_N:Country + pH_H2O:Country + C_N:pH_H2O:Country, df.sub)
rmod.na <- betareg(diversity_gini_simpson ~ Country + C_N:Country + pH_H2O:Country + C_N:pH_H2O:Country, df.sub)
rmod.cn <- betareg(diversity_gini_simpson ~ Country + Na + pH_H2O:Country + C_N:pH_H2O:Country, df.sub)
rmod.ph <- betareg(diversity_gini_simpson ~ Country + Na + C_N:Country + C_N:pH_H2O:Country, df.sub)
rmod.cnph <- betareg(diversity_gini_simpson ~ Country + Na + C_N:Country + pH_H2O:Country, df.sub)

r2_full <- summary(rmod)$pseudo.r.squared
r2_no_country <- summary(rmod.country)$pseudo.r.squared
r2_no_na <- summary(rmod.na)$pseudo.r.squared
r2_no_cn <- summary(rmod.cn)$pseudo.r.squared
r2_no_ph <- summary(rmod.ph)$pseudo.r.squared
r2_no_cnph <- summary(rmod.cnph)$pseudo.r.squared

r2_country <- r2_full - r2_no_country
r2_na <- r2_full - r2_no_na
r2_cn <- r2_full - r2_no_cn
r2_ph <- r2_full - r2_no_ph
r2_cnph <- r2_full - r2_no_cnph

r2_full
r2_country
r2_na
r2_cn
r2_ph
r2_cnph





head(df.sub)

# diversity_gini_simpson ~ Country + Na + C_N:Country + pH_H2O:Country + C_N:pH_H2O:Country

ggplot(df.sub) + aes(x = Na, y = diversity_gini_simpson) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(group = Country)) +
  facet_wrap(.~Country)

ggplot(df.sub) + aes(x = C_N*pH_H2O, y = diversity_gini_simpson) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(group = Country)) +
  facet_wrap(.~Country)

ggplot(df.sub) + aes(x = C_N, y = diversity_gini_simpson) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(group = Country)) +
  facet_wrap(.~Country)

ggplot(df.sub) + aes(x = pH_H2O, y = diversity_gini_simpson) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(group = Country)) +
  facet_wrap(.~Country)

ggplot(df.sub) + aes(x = C_N*pH_H2O, y = diversity_shannon) +
  geom_point(size=3, aes(color = Country)) +
  geom_smooth(method = "lm", aes(group = Country)) +
  facet_wrap(.~Country)

#====================================#
####            END               ####
#====================================#

