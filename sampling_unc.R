library(ggplot2)
library(reshape)
library(dplyr)
library(caret)
library(ranger)
library(glmnet)
library(SuperLearner)
library("e1071")
library(tune)
library(metafor)
library(foreign)
library(ggsci)
library(plotrix)
library(leaps)
library(splines)
library(gam)
library("xgboost")
library("Matrix")
library(class)
library(mosaic)
set.seed(817)

# loading data ####


data = read.spss("gusto.sav", to.data.frame=TRUE)

summary(data)

factors = c("DAY30","SEX","KILLIP","A65", "SHO", "DIA", "HYP", "ANT","HRT", 
            "PMI", "HTN", "SMK", "LIP", "PAN", "FAM", "ST4", "TTR", "HIG")


data[factors] = lapply(data[factors], factor)

partA = data[data$ESAMP == 1, ]  #training data
partB = data[data$ESAMP == 2, ]  #test data




# standard error by regression  model####

se_ind = sample(1:nrow(partA), 10000)
sample_se = partA[se_ind , c("DAY30", "AGE","SHO", "ANT", "DIA",
                                   "HYP","HRT","TTR","SEX")]

## model
logmod_se = glm(DAY30 ~ ., family = "binomial", data = sample_se)
## predicted probabilities
preds_logmod = predict(logmod_se, newdata = sample_se,
                             type = "response")
##predictions on logit scale
preds_logmod_linear = predict(logmod_se, newdata = sample_se,
                                    type = "link")


##standard errors on probability and logit scale with predict()
pred_se_mod = predict(logmod_se, newdata = sample_se,  
                            se.fit = TRUE, type = "response")$se.fit

pred_se_mod_linear = predict(logmod_se, newdata = sample_se,  
                                   se.fit = TRUE, type = "link")$se.fit

## standard errors on probability scale by hand with delta method
pred_se_mod_hand = pred_se_mod_linear *(preds_logmod * (1 - preds_logmod))


ggplot(data = data.frame(cbind(pred_se_mod, pred_se_mod_hand)))+
  geom_point(x = pred_se_mod, y = pred_se_mod_hand)+
  scale_x_continuous(limits = c(0,0.05))+
  scale_y_continuous( limits = c(0,0.05))



#standard error by bootstrap resampling####

B=1000

#prob scale
preds_se = matrix(NA, nrow = nrow(sample_se), ncol = B)

for (i in 1:B){
  ind_b = sample(1:nrow(sample_se), nrow(sample_se), replace = TRUE)
  boot_resamp = sample_se[ind_b, ]
  boot_mod = glm(DAY30 ~ ., family = "binomial", data = boot_resamp)
  preds_se[ ,i] = predict(boot_mod, newdata = sample_se, 
                                type = "response") 
}

#linear scale
preds_se_linear = matrix(NA, nrow = nrow(sample_se), ncol = B)

for (i in 1:B){
  ind_b = sample(1:nrow(sample_se), nrow(sample_se), replace = TRUE)
  boot_resamp = sample_se[ind_b, ]
  boot_mod = glm(DAY30 ~ ., family = "binomial", data = boot_resamp)
  preds_se_linear[ ,i] = predict(boot_mod, newdata = sample_se, 
                          type = "link") 
}



#se and var on the probability scale
pred_se_boot = apply(preds_se, MARGIN = 1, FUN = sd)
pred_var_boot = apply(preds_se,  MARGIN = 1, FUN = var)

#se and var on the logit scale 
pred_se_boot_linear = apply(preds_se_linear, MARGIN = 1, FUN = sd)
pred_var_boot_linear = apply(preds_se_linear,  MARGIN = 1, FUN = var)




#plots on agreement between predict(se) and bootstrap####

##probability scale####

se_df = data.frame(cbind(pred_se_boot,
                            pred_se_mod))

names(se_df) = c("boot", "mod")

ggplot(data = se_df)+
  geom_density(aes(x=mod, color = "Model"), linewidth = 1 )+
  geom_density(aes(x=boot, color = "Bootstrap"), linewidth = 1)+
  scale_color_manual(name = "Estimation \nmethod", 
                     values = c("Model" = "darkgrey",
                                "Bootstrap" = "skyblue"))+
  xlab ("standard errors on probability scale")




ggplot(data = se_df)+
  geom_point(aes(x=mod, y = boot), alpha = 0.2)+
  geom_abline(slope = 1)+
  ylab("Standard errors estimated by bootstrap")+
  xlab("Standard errors estimated by the model")

ggsave("se_agr_probscale1.png")


##linear scale####
se_df_linear = data.frame(cbind(pred_se_boot_linear,
                                   pred_se_mod_linear))

names(se_df_linear) = c("boot", "mod")

ggplot(data = se_df_linear)+
  geom_density(aes(x=mod, color = "Model"), linewidth = 1 )+
  geom_density(aes(x=boot, color = "Bootstrap"), linewidth = 1)+
  scale_color_manual(name = "Estimation \nmethod", 
                     values = c("Model" = "darkgrey",
                                "Bootstrap" = "skyblue"))+
  xlab ("standard errors on linear scale")+
  theme(legend.position = "none")



ggplot(data = se_df_linear)+
  geom_point(aes(x = mod, y = boot), alpha = 0.2)+
  geom_abline(slope = 1)+
  ylab("Standard errors estimated by bootstrap")+
  xlab("Standard errors estimated by the model")







