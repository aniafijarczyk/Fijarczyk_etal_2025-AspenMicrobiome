library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(cowplot)
library(nlme)
#library(performance) # for Nakagawa conditional/marginal R2
#library(partR2) # for part R2 values
#library(rcompanion) # to calculate pseudo R2 for gls models
library(emmeans) # to get Tukey test on gls model




#########################
#------FUNCTIONS--------#
#########################



fit_ANOVA <- function(dataset) {
  check.mod <- lm(stat ~ group1 / group2, dataset)
  shap <- shapiro.test(residuals(check.mod, type = "pearson"))
  print(paste0("Shapiro p: ",shap$p.value))
  hom.group1 <- bartlett.test(residuals(check.mod, type = "pearson") ~ dataset$group1)
  hom.group2 <- bartlett.test(residuals(check.mod, type = "pearson") ~ dataset$group2)
  print(paste0("Bartlett p GROUP1: ",hom.group1$p.value))
  print(paste0("Bartlett p GROUP2: ",hom.group2$p.value))
  # Getting results
  df.mod <- as.data.frame(anova(check.mod))
  TSS <- sum(df.mod$`Sum Sq`)
  df.mod$Variance <- df.mod$`Sum Sq`/TSS
  
  out.mod <- data.frame("shap" = c(shap$p.value),
                          "bart_GR1" = c(hom.group1$p.value),
                          "bart_GR2" = c(hom.group2$p.value),
                          "Df_GR1" = df.mod["group1","Df"],
                          "F_GR1" = df.mod["group1","F value"],
                          "R2_GR1" = df.mod["group1","Variance"],
                          "p_GR1" = df.mod["group1","Pr(>F)"],
                          "Df_GR2" = df.mod["group1:group2","Df"],
                          "F_GR2" = df.mod["group1:group2","F value"],
                          "R2_GR2" = df.mod["group1:group2","Variance"],
                          "p_GR2" = df.mod["group1:group2","Pr(>F)"])
  
  # Tukey
  reg.aov <- aov(stat ~ group1 / group2, data=dataset)
  tukey.test <- TukeyHSD(reg.aov)
  pairs <- rownames(tukey.test$group1)
  #tt = list()
  for (p in pairs) {
    #tt[[paste0("TT_",p)]] <- tukey.test$group1[`p`, 'p adj']
    out.mod[paste0("TT_",p)] <- tukey.test$group1[`p`, 'p adj']
  }

  out <- as.data.frame(t(out.mod))
  out$statistic <- rownames(out)
  rownames(out) <- c(1:length(out$V1))
  out.sort <- out[,c("statistic","V1")]
  return(out.sort)
}



fit_GLS <- function(dataframe) {
  mod.gls <- gls(stat ~ group1,
                 correlation = corCompSymm(form = ~ 1 | group1/group2),
                 weights = varIdent(form = ~ 1|group1),
                 dataframe)
  mod.gls.lev2 <- gls(stat ~ group2,
                      correlation = corCompSymm(form = ~ 1 | group1/group2),
                      weights = varIdent(form = ~ 1|group2),
                      dataframe)
  # Pseudo R2 for group1
  #R2.gls <- rcompanion::nagelkerke(mod.gls)
  #R2.gls.nag <- as.data.frame(R2.gls$Pseudo.R.squared.for.model.vs.null)$"Pseudo.R.squared"[3]
  # Pseudo R2 for group2
  #R2.gls.lev2 <- rcompanion::nagelkerke(mod.gls.lev2)
  #R2.gls.lev2.nag <- as.data.frame(R2.gls.lev2$Pseudo.R.squared.for.model.vs.null)$"Pseudo.R.squared"[3]
  
  df.mod.gls <- as.data.frame(anova(mod.gls))
  df.mod.gls.2 <- as.data.frame(anova(mod.gls.lev2))
  
  out.mod <- data.frame("Df_GR1" = df.mod.gls["group1","numDF"],
                        "F_GR1" = df.mod.gls["group1","F-value"],
                        "R2_GR1" = NA,
                        "p_GR1" = df.mod.gls["group1","p-value"],
                        "Df_GR2" = df.mod.gls.2["group2","numDF"],
                        "F_GR2" = df.mod.gls.2["group2","F-value"],
                        "R2_GR2" = NA,
                        "p_GR2" = df.mod.gls.2["group2","p-value"])
  
  # Tukey
  mod.gls.emm <- emmeans(mod.gls, "group1", data=dataframe)
  df.gls.emm <- data.frame(pairs(mod.gls.emm))
  pairs <- df.gls.emm$contrast
  for (p in pairs) {
    np <- gsub(" - ","-",p)
    out.mod[paste0("TT_",np)] <- df.gls.emm[df.gls.emm$contrast == p,'p.value']
  }
  
  out <- as.data.frame(t(out.mod))
  out$statistic <- rownames(out)
  rownames(out) <- c(1:length(out$V1))
  out.sort <- out[,c("statistic","V1")]
  return(out.sort)
}



#########################################
run_test <- function(dataset, metric_value, group_1, group_2) {
  
  ds <- dataset
  
  # Assign stat
  ds$stat <- ds[[metric_value]]
  ds$group1 <- ds[[group_1]]
  ds$group2 <- ds[[group_2]]
  check.mod <- lm(stat ~ group1 / group2, data=ds)
  bc <- boxcox(check.mod, lambda = seq(-10, 10), data = ds)
  best.lambda <- bc$x[which(bc$y==max(bc$y))]
  
  model.3.out <- fit_GLS(ds)
  return(model.3.out)
}

#FILENAME="test_file.tsv"
#df <- read.csv(FILENAME, sep="\t", header=T)
#head(df)

#run_nested_models(df, "chao1", "Country", "site")


run_nested_models <- function(dataset, metric_value, group_1, group_2) {
  
  ds <- dataset
  
  # Assign stat
  ds$stat <- ds[[metric_value]]
  ds$group1 <- ds[[group_1]]
  ds$group2 <- ds[[group_2]]
  
  # Check missing values
  st1 <- length(ds$stat)
  st2 <- sum(complete.cases(ds$stat))
  print(paste0("There are ",st1," total values and ",st1-st2," are missing"))
  
  # Fit model - check normality for selected groups
  print("Fit model 1")
  model.1.out <- fit_ANOVA(ds)
  model.1.out$MOD <- "ANOVA"
  model.1.out
  
  # Assert if model meets the assumptions
  print("Assert")
  ass1 <- model.1.out[model.1.out$statistic=="shap",'V1']>0.05
  ass2 <- model.1.out[model.1.out$statistic=="bart_GR1",'V1']>0.05
  if (ass1 == TRUE & ass2==TRUE) {
    model.assert <- "OK"
  } else {model.assert <- "-"}
  model.assert
  model.1.out$EVAL <- model.assert
  model.1.out
  
  # Value transformation
  print("Value transformation")
  check.mod <- lm(stat ~ group1 / group2, data=ds)
  bc <- boxcox(check.mod, lambda = seq(-10, 10), data = ds)
  best.lambda <- bc$x[which(bc$y==max(bc$y))]
  print("Create new dataframe")
  new_ds <- ds
  new_ds$stat <- ds$stat^best.lambda
  head(new_ds)
  
  # Fit model with transformed variable - check assumptions
  print("Fit model 2")
  model.2.out <- fit_ANOVA(new_ds)
  model.2.out$MOD <- "ANOVA_transformed"
  model.2.out
  
  # Assert if model meets the assumptions
  print("Assert 2")
  ass1 <- model.2.out[model.2.out$statistic=="shap",'V1']>0.05
  ass2 <- model.2.out[model.2.out$statistic=="bart_GR1",'V1']>0.05
  if (ass1 == TRUE & ass2==TRUE & model.assert == '-') {
    model.assert.2 <- "OK"
  } else {model.assert.2 <- "-"}
  model.2.out$EVAL <- model.assert.2
  model.2.out
  
  # Fit GLS model
  print("Fit GLS")
  model.3.out <- fit_GLS(new_ds)
  model.3.out$MOD <- "GLS"
  
  # Label model
  if (model.assert == "-" & model.assert.2 == "-") {
    model.assert.3 <- "OK"
  } else {model.assert.3 <- "-"}
  model.3.out$EVAL <- model.assert.3
  
  ### Combining models
  print("Combine")
  models <- rbind(model.1.out, model.2.out, model.3.out)
  return(models)

}


#df.models <- run_nested_models(df, "chao1", "Country", "site")
#head(df.models)

