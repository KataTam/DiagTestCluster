########################################################
### Script to calculate binary-scale diagnostic test ###
### characteristics in the context of clustering     ###
########################################################

# written by Katalin Tamasi between July 2021 and May 2022
# k.tamasi@umcg.nl

# References:
# Fleiss et al. (2003);
# Genders et al. (2012); 
# Hujoel, Moulton, & Loesche (1990); 
# Kirkwood & Sterne (2003);
# McDonald (2019);
# Mercaldo et al. (2007);
# Ronco & Biggeri  (1999);
# Williams (2000); and
# Ying et al. (2020) 

#############
### Input ###
#############

# data structure in wide format (i.e., one row per patient)
# id = patient identification number (identifies the clusters i =1…I)
# TP = number of true-positive observations
# FN = number of false-negative observations
# FP = number of false-positive observations
# TN = number of true-negative observations
# Test_type (optional) = indicator variable if multiple index tests are considered (eg, 0 = test A, 1 = test B)

##############
### Output ###
##############

# Results generated from input:

############################
## Patient-level analyses ##
############################

# Patient-level contingency table
# Prevalence-independent measures: Type_of_method, Sensitivity (%) Lower_CI Higher_CI, Specificity (%) Lower_CI Higher_CI
# Prevalence-dependent measures: Type_of_method, PPV (%) Lower_CI Higher_CI, NPV (%) Lower_CI Higher_CI
# Patient-level likelihood ratios
# 
############################
## Segment-level analyses ##
############################

# Segment-level contingency table
# Prevalence-independent measures: Type_of_method, Sensitivity (%) Lower_CI Higher_CI, Specificity (%) Lower_CI Higher_CI
# Prevalence-dependent measures Type_of_method, PPV (%) Lower_CI Higher_CI, NPV (%) Lower_CI Higher_CI
# Segment-level likelihood ratios (NB: CI's not adjusted for clustering!)
# Intracluster correlation coefficients
# 

####################
## Visualizations ##
####################

# Forest plots: Allow for comparison of outcomes obtained by each method.

##################
## Some remarks ##
##################

# -- The methods reported by Genders et al. (2012) are validated using their dataset and supplementary materials.
# They report slightly different CI's in the main text vs. in the supplementary materials, the latter of which
#  are exactly reproduced by the R script. CI's from the mixed effects logistic regression (called logistic
# random-effects model by Genders et al., 2012) are slightly different due to the default number of quadrature
# points when approximating the integral over the random effects structure (1 in R's glmer vs. 7 in STATA's xtmelogit).
# 
# -- When sourcing the script, use the "print.eval=TRUE" option so that the forest plots can be generated, e.g.: source("~yourpath/Advanced_screening_diagnostics_Tamasi_website.R", print.eval=TRUE)


#####################
### SCRIPT BEGINS ###
#####################

## Installing the required packages. This takes a couple of minutes but needs to be run only once. 
# Uncomment as necessary to install
# install.packages("dplyr")
# install.packages("data.table")
# install.packages("rstatix")
# install.packages("boot")
# install.packages("sandwich")
# install.packages("lmtest")
# install.packages("lme4")
# install.packages("geepack")
# install.packages("gee")
# install.packages("ClusterBootstrap")
# install.packages("readr")
# install.packages("stringr")
# install.packages("data.tree")
# install.packages("igraph")
# install.packages("DiagrammeR")
# install.packages("forestplot")
# install.packages("tidyr")
# install.packages("ggplot2")
# install.packages("beepr")
# install.packages("bootLR")
# install.packages("multcomp")

# loading packages
library(dplyr)
library(tidyr)
library(data.table)
library(rstatix)
library(boot)
library(sandwich)
library(lmtest)
library(lme4)
library(gee) # contains function gee: data needs to be sorted by id if not already. preferable to geeglm due to 
# giving a more informative output (working correlation matrix)
library(ClusterBootstrap)
library(readr)
library(stringr)
library(data.tree)
library(dplyr)
library(igraph)
library(DiagrammeR)
library(forestplot)
library(beepr)
library(bootLR)
library(multcomp)


# clearing workspace 
rm(list=ls())

# setting working directory, change path to yours
setwd("/Users/katalintamasi/Dropbox/Uni/Groningen/Teaching/Stats_literacy_journal_club/Session_5_Advanced_screening_diagnostics/CAD")

# reading in example dataset 
data <- read.table("CAD_data_Genders.csv", header=T, sep = ",")

# inspecting data for possible typos
#View(data)
#summary(data)


# defining test and disease on the patient level
data$Test <- ifelse(data$TP + data$FP != 0,1,0)
data$Inv_test <- (data$Test - 1) * -1
data$Disease <- ifelse(data$TP + data$FN != 0,1,0)
data$Inv_disease <- (data$Disease - 1) * -1

# 2x2 patient-level contingency table: can be used to draw table or natural frequency tree
#table(data$Disease, data$Test, useNA = "ifany" )

id_cont_table <- table(data$Disease, data$Test, useNA = "ifany")
rownames(id_cont_table) <- c("Non-diseased", "Diseased")
colnames(id_cont_table) <- c("Negative test", "Positive test")
# id_cont_table
id_sum_TP <- id_cont_table[2,2]
id_sum_FN <- id_cont_table[2,1]
id_sum_TN <- id_cont_table[1,1]
id_sum_FP <-id_cont_table[1,2]

# calculate patient-level LR test
set.seed(5434355)
id_bootlrtest <- BayesianLR.test(id_sum_TP, id_sum_FN+id_sum_TP, id_sum_TN, id_sum_TN+id_sum_FP, R = 10000)

id_blr <- matrix(c(id_bootlrtest$posLR, id_bootlrtest$negLR, id_bootlrtest$posLR.ci[1],  id_bootlrtest$negLR.ci[1], id_bootlrtest$posLR.ci[2],  id_bootlrtest$negLR.ci[2]), ncol = 3)

colnames(id_blr) <- c("Point estimate", "Lower_CI", "Higher_CI")
rownames(id_blr) <- c("Positive likelihood ratio", "Negative likeilhood ratio")

# segment-level contingency table
# sums
sum_TP <- as.integer(sum(data$TP))
sum_FN <- as.integer(sum(data$FN))
sum_TN <- as.integer(sum(data$TN))
sum_FP <- as.integer(sum(data$FP))

cont_table <- as.table(matrix(c(sum_TN, sum_FN, sum_FP, sum_TP), ncol = 2))
rownames(cont_table) <- c("Non-diseased", "Diseased")
colnames(cont_table) <- c("Negative test", "Positive test")
#cont_table

# calculate segment-level LR test
set.seed(5434355)
bootlrtest <- BayesianLR.test(sum_TP, sum_FN+sum_TP, sum_TN, sum_TN+sum_FP, R = 10000)

blr <- matrix(c(bootlrtest$posLR, bootlrtest$negLR, bootlrtest$posLR.ci[1],  bootlrtest$negLR.ci[1], bootlrtest$posLR.ci[2],  bootlrtest$negLR.ci[2]), ncol = 3)

colnames(blr) <- c("Point estimate", "Lower_CI", "Higher_CI")
rownames(blr) <- c("Positive likelihood ratio", "Negative likeilhood ratio")


# some calculations require data to be in long format (i.e., one row per segment)
# converting from wide to long
data.long <- gather(data, Status, Obs, TP:FP, factor_key=TRUE)

# sort
data.long <- data.long[order(data.long$id, data.long$Status), ]

# dropping rows with empty observations
data.long <- subset(data.long, data.long$Obs != 0)

# repeating rows, each according to the number of observations 
data.long <- as.data.frame(lapply(data.long, rep, data.long$Obs))

# creating Segment as an index variable
# use pipes ( %>% ) to do this in a single line of code
data.long <-data.long %>% group_by(id) %>% mutate(Segment=1:n())

# total number of observations
setDT(data.long)[,Observation:=seq(1,.N)]

#summary(data.long)

data.long$Test <- ifelse(data.long$Status == "TP" | data.long$Status == "FP",1,0)
data.long$Inv_test <- (data.long$Test - 1) * -1
data.long$Disease <- ifelse(data.long$Status == "TP" | data.long$Status == "FN",1,0)
data.long$Inv_disease <- (data.long$Disease - 1) * -1
#data.long$Test_type <- 0 #rep(c(0, 1), round(nrow(data.long)/2)) # alternating test types, e.g., 0 = test A, 1 = test B, need to be filled in by hand in the actual analysis

# Note from Genders et al., 2012 Supp.—In a paired study design, each observation is entered twice, once with test_type==0 
# and once with test_type==1.


# reordering columns
#colnames(data.long)
col_order <- c("Observation", "id", "Segment",
               "Test", "Inv_test", "Disease", "Inv_disease")
data.long <- data.long[, ..col_order]
#data.long
#summary(data.long)

# creating some necessary variables
data$sens <- data$TP / (data$TP + data$FN) # indidivual sensitivity (true positive rate)
data$spec <- data$TN / (data$FP + data$TN) # individual specificity (true negative rate)
data$D <- data$TP + data$FN # number of diseased observations
data$N <- data$FP + data$TN # number of non-diseased observations
data$j <- data$D + data$N # total no of observations

# sums
sum_TP <- as.integer(sum(data$TP))
sum_FN <- as.integer(sum(data$FN))
sum_TN <- as.integer(sum(data$TN))
sum_FP <- as.integer(sum(data$FP))
sum_D <- as.integer(sum(data$D))
sum_N <- as.integer(sum(data$N))        
sum_j <- as.integer(sum(data$j))

# weigthed sums
data$Nw <- data$j/sum_j * data$j # overall weighting
data$Nw_sens <- data$D/sum_D * data$D # weighting for sensitivity
data$Nw_spec <- data$N/sum_N * data$N # weighting for specificity

data$Test <- ifelse(data$TP + data$FP != 0,1,0)
data$Inv_test <- (data$Test - 1) * -1
data$Disease <- ifelse(data$TP + data$FN != 0,1,0)
data$Inv_disease <- (data$Disease - 1) * -1
# optional variable data$Test_type: need to be imputed in original data file



#################################
#### Patient-level analyses #####
#################################
# equivalent to aggregate analysis

#######################################
### Prevalence-independent measures ###
#######################################

###################################
## No adjustment for clustering  ##
###################################
# adjustment not possible if patients are independent -- if they are not, another level of clustering needs to be defined

# number of true positive on the patient level: 
# TP if anywhere within a patient disease is found using the reference test (Gender's example: conventional angiography) 
# and if anywhere within a patient the index test (CT angiography) is positive

table(data$Test, data$Disease)

id_no_TP <- sum(data$Test == 1 & data$Disease == 1) # number of patients with diseased observations who have at least one TP
id_no_FP <- sum(data$Test == 1 & data$Disease == 0) # 

# within patients with no diseased observations (id_sum_N), there is no false positive
id_no_TN <- sum(data$Test == 0 & data$Disease == 0) # number of patients with no diseased observations who have no FP
id_no_FN <- sum(data$Test == 0 & data$Disease == 1) # 

id_sum <- nlevels(as.factor(data$id)) # number of clusters (in this example, patients)
id_sum_D <- sum(data$D != 0) # number of patients with diseased observations
id_sum_N <- sum(data$N != 0) # number of patients with non-diseased observations 


#########################
## Binomial proportion ##
#########################

# sensitivity
id_sens <- id_no_TP/(id_no_TP+id_no_FN) 

# specificity
id_spec <- id_no_TN/(id_no_FP+id_no_TN) 

#######################################
# Binomial proportion, Wald intervals #
#######################################

# unadjusted SE (based on binomial distribution) and 95% CI (based on Wald normal approximation) 
# sensitivity
id_sens_SE_unadj <- sqrt((id_sens*(1-id_sens))/(id_sum_D))
id_sens_CI_unadj <- 1.96 * id_sens_SE_unadj
id_sens_lower_CI_unadj <- id_sens - id_sens_CI_unadj
id_sens_higher_CI_unadj <- id_sens + id_sens_CI_unadj    

# id_sens
# id_sens_lower_CI_unadj 
# id_sens_higher_CI_unadj

# specificity
id_spec_SE_unadj <- sqrt((id_spec*(1-id_spec))/(id_sum-id_sum_D))
id_spec_CI_unadj <- 1.96 * id_spec_SE_unadj
id_spec_lower_CI_unadj <- id_spec-id_spec_CI_unadj
id_spec_higher_CI_unadj <- id_spec+id_spec_CI_unadj

# id_spec
# id_spec_lower_CI_unadj
# id_spec_higher_CI_unadj


########################################
# Binomial proportion, exact intervals #
########################################

# sens
id_bt_sens <- binom_test(id_no_TP, id_sum_D)
id_sens_exact <- id_bt_sens$estimate
id_sens_lower_CI_exact <- id_bt_sens$conf.low
id_sens_higher_CI_exact <- id_bt_sens$conf.high

# id_sens_exact
# id_sens_lower_CI_exact
# id_sens_higher_CI_exact

# spec
id_bt_spec <- binom_test(id_no_TN, (id_sum-id_sum_D))
id_spec_exact <- id_bt_spec$estimate
id_spec_lower_CI_exact <- id_bt_spec$conf.low
id_spec_higher_CI_exact <- id_bt_spec$conf.high

# id_spec_exact
# id_spec_lower_CI_exact
# id_spec_higher_CI_exact

#######################
# Logistic regression #
#######################

# sensitivity
#glm, -1 suppresses intercept
glm1 <- glm(Test ~ Disease - 1, family = binomial(link="logit"), data=data)
#summary(glm1)

# probability scale (exp/(1+exp))
id_sens_exp <- exp(summary(glm1)$coefficients[,1]) / (1+(exp(summary(glm1)$coefficients[,1])))

id_sens_lower_CI_exp <- exp(summary(glm1)$coefficients[,1] + qnorm(0.025) *
                                     summary(glm1)$coefficients[,2]) / (1+exp(summary(glm1)$coefficients[,1] + qnorm(0.025) *
                                                                                               summary(glm1)$coefficients[,2]))

id_sens_higher_CI_exp <- exp(summary(glm1)$coefficients[,1] + qnorm(0.975) *
                                      summary(glm1)$coefficients[,2]) / (1+exp(summary(glm1)$coefficients[,1] + qnorm(0.975) *
                                                                                                summary(glm1)$coefficients[,2]))

# id_sens_exp
# id_sens_lower_CI_exp
# id_sens_higher_CI_exp

# specificity
#glm, -1 suppresses intercept
glm2 <- glm(Inv_test ~ Inv_disease - 1, family = binomial(link="logit"), data=data)
#summary(glm2)

# probability scale (exp/(1+exp))
id_spec_exp <- exp(summary(glm2)$coefficients[,1]) / (1+(exp(summary(glm2)$coefficients[,1])))

id_spec_lower_CI_exp <- exp(summary(glm2)$coefficients[,1] + qnorm(0.025) *
                                    summary(glm2)$coefficients[,2]) / (1+exp(summary(glm2)$coefficients[,1] + qnorm(0.025) *
                                                                                     summary(glm2)$coefficients[,2]))

id_spec_higher_CI_exp <- exp(summary(glm2)$coefficients[,1] + qnorm(0.975) *
                                     summary(glm2)$coefficients[,2]) / (1+exp(summary(glm2)$coefficients[,1] + qnorm(0.975) *
                                                                                      summary(glm2)$coefficients[,2]))

# id_spec_exp
# id_spec_lower_CI_exp
# id_spec_higher_CI_exp

#####################################
# Logistic regression, single model #
#####################################

# modeling sens and spec in a single model (e.g., Ronco & Biggeri, 1999)
glm1_single_model <- glm(Test ~ Disease, family=binomial(link = "logit"), data=data)  
#summary(glm1_single_model)

# calculating sensitivity from single model: exp(b1+b2) / (1+exp(b1+b2)), where b1 is the intercept and b2 is the coeff
id_sens_glm_single <- (exp(summary(glm1_single_model)$coefficients[[1]] + summary(glm1_single_model)$coefficients[[2]])) / 
        (1+(exp(summary(glm1_single_model)$coefficients[[1]] + summary(glm1_single_model)$coefficients[[2]])))

#names(fixef(glm1_single_model))[1]

# calculating 95% CI for sensitivity from single model: exp(L1+L2) / (1+exp(L1+L2)), where L1 and L2 are the lower bounds of 
# the intercept and coeff
# lincom to calculate 95% CI for sensitivity
#summary(glht(glm1_single_model, linfct = c("`(Intercept)` + Disease = 0")))

# Show the confidence interval
glm1_single_model.lh <- glht(glm1_single_model, linfct = c("`(Intercept)` + Disease = 0"))
a <- confint(glm1_single_model.lh)
# lower CI
#a$confint[[2]]
# higher CI
#a$confint[[3]]

id_sens_lower_CI_glm_single <- (exp(a$confint[[2]])) / 
        (1+exp(a$confint[[2]]))

id_sens_higher_CI_glm_single <- (exp(a$confint[[3]])) / 
        (1+exp(a$confint[[3]]))


# calculating specificity from single model: 1-exp(b1) / (1+exp(b1)), where b1 is the intercept
id_spec_glm_single <- 1-exp(summary(glm1_single_model)$coefficients[[1]]) / (1+exp(summary(glm1_single_model)$coefficients[[1]]))


id_spec_higher_CI_glm_single <- 1- (exp(summary(glm1_single_model)$coefficients[[1]] + qnorm(0.025) *
                                             summary(glm1_single_model)$coefficients[[3]])) / (1+exp(summary(glm1_single_model)$coefficients[[1]] + qnorm(0.025) *
                                                                                                             summary(glm1_single_model)$coefficients[[3]]))

id_spec_lower_CI_glm_single <- 1- (exp(summary(glm1_single_model)$coefficients[[1]] + qnorm(0.975) *
                                            summary(glm1_single_model)$coefficients[[3]])) / (1+exp(summary(glm1_single_model)$coefficients[[1]] + qnorm(0.975) *
                                                                                                            summary(glm1_single_model)$coefficients[[3]]))

# id_sens_exp
# id_sens_glm_single
# id_sens_lower_CI_exp
# id_sens_lower_CI_glm_single
# id_sens_higher_CI_exp
# id_sens_higher_CI_glm_single
# 
# id_spec_exp
# id_spec_glm_single
# id_spec_lower_CI_exp
# id_spec_lower_CI_glm_single
# id_spec_higher_CI_exp
# id_spec_higher_CI_glm_single

#####################################
### Prevalence-dependent measures ###
#####################################

#########################
## Binomial proportion ##
#########################

# predictive values for patient-level analyses -- they can be calculated given the available info: sens, spec + prev of disease
# calculations based on Ying et al., 2020
id_prev <- id_sum_D / id_sum # not always appropriate to estimate from study (sampling bias), so usually estimated from separate study

# PPV
id_PPV <- (id_sens * id_prev) / ((id_sens * id_prev) + (1-id_spec) * (1 - id_prev)) # patient-level positive predictive value

# NPV
id_NPV <- (id_spec * (1- id_prev)) / ((1-id_sens) * id_prev + id_spec * (1 - id_prev)) # patient-level negative predictive value

#####################################
# Binomial proportion, Delta method #
#####################################
# Mercaldo et al., 2007

# unadjusted standard error of predictive values
id_PPV_SE_unadj <- sqrt((((id_prev * (1 - id_spec)) * (1 - id_prev))^2 * ((id_sens * (1 - id_sens))/id_sum_D) + 
                                 (id_prev * id_sens * (1 - id_prev))^2 * ((id_spec * (1 - id_spec)) / (id_sum-id_sum_D))) / 
                                ((id_sens * id_prev + (1 - id_spec) * (1 - id_prev))^4))

id_PPV_CI_unadj <- 1.96*id_PPV_SE_unadj
id_PPV_lower_CI_unadj <- id_PPV-id_PPV_CI_unadj
id_PPV_higher_CI_unadj <- id_PPV+id_PPV_CI_unadj   

# id_PPV
# id_PPV_lower_CI_unadj 
# id_PPV_higher_CI_unadj


id_NPV_SE_unadj <- sqrt((((id_spec * (1 - id_prev)) * id_prev)^2 * ((id_sens * (1 - id_sens)) / id_sum_D) +
                                 ((1 - id_sens) * (1 - id_prev) * id_prev)^2 * ((id_spec * (1  - id_spec)) / (id_sum-id_sum_D))) /
                                (((1 - id_sens) * id_prev + id_spec * (1 - id_prev))^4))

id_NPV_CI_unadj <- 1.96*id_NPV_SE_unadj
id_NPV_lower_CI_unadj <- id_NPV-id_NPV_CI_unadj
id_NPV_higher_CI_unadj <- id_NPV+id_NPV_CI_unadj     


# id_NPV
# id_NPV_lower_CI_unadj 
# id_NPV_higher_CI_unadj
 
########################################
# Binomial proportion, exact intervals #
########################################

# PPV
id_bt_PPV <- binom_test(id_no_TP, (id_no_TP + id_no_FP))
id_PPV_exact <- id_bt_PPV$estimate
id_PPV_lower_CI_exact <- id_bt_PPV$conf.low
id_PPV_higher_CI_exact <- id_bt_PPV$conf.high

# id_PPV_exact
# id_PPV_lower_CI_exact
# id_PPV_higher_CI_exact

# NPV
id_bt_NPV <- binom_test(id_no_TN, (id_no_TN + id_no_FN))
id_NPV_exact <- id_bt_NPV$estimate
id_NPV_lower_CI_exact <- id_bt_NPV$conf.low
id_NPV_higher_CI_exact <- id_bt_NPV$conf.high

# id_NPV_exact
# id_NPV_lower_CI_exact
# id_NPV_higher_CI_exact

#######################
# Logistic regression #
#######################

# PPV
#glm, -1 suppresses intercept
glm3 <- glm(Disease ~ Test - 1, family = binomial(link="logit"), data=data)
#summary(glm3)

# probability scale (exp/(1+exp))
id_PPV_exp <- exp(summary(glm3)$coefficients[,1]) / (1+exp(summary(glm3)$coefficients[,1]))

id_PPV_lower_CI_exp <- (exp(summary(glm3)$coefficients[,1] + qnorm(0.025) *
                                    summary(glm3)$coefficients[,2]) / (1+exp(summary(glm3)$coefficients[,1] + qnorm(0.025) *
                                                                                     summary(glm3)$coefficients[,2])))

id_PPV_higher_CI_exp <- (exp(summary(glm3)$coefficients[,1] + qnorm(0.975) *
                                     summary(glm3)$coefficients[,2]) / (1+exp(summary(glm3)$coefficients[,1] + qnorm(0.975) *
                                                                                      summary(glm3)$coefficients[,2])))

# id_PPV_exp
# id_PPV_lower_CI_exp
# id_PPV_higher_CI_exp

# NPV
#glm, -1 suppresses intercept
glm4 <- glm(Inv_disease ~ Inv_test - 1, family = binomial(link="logit"), data=data)
#summary(glm4)

# probability scale (exp/(1+exp))
id_NPV_exp <- exp(summary(glm4)$coefficients[,1]) / (1+exp(summary(glm4)$coefficients[,1]))

id_NPV_lower_CI_exp <- (exp(summary(glm4)$coefficients[,1] + qnorm(0.025) *
                                    summary(glm4)$coefficients[,2]) / (1+exp(summary(glm4)$coefficients[,1] + qnorm(0.025) *
                                                                                     summary(glm4)$coefficients[,2])))

id_NPV_higher_CI_exp <- (exp(summary(glm4)$coefficients[,1] + qnorm(0.975) *
                                     summary(glm4)$coefficients[,2]) / (1+exp(summary(glm4)$coefficients[,1] + qnorm(0.975) *
                                                                                      summary(glm4)$coefficients[,2])))

# id_NPV_exp
# id_NPV_lower_CI_exp
# id_NPV_higher_CI_exp

#####################################
# Logistic regression, single model #
#####################################

# PPV
glm3_single_model <- glm(Disease ~ Test, family = binomial(link="logit"), data=data)
#summary(glm3_single_model)

# calculating PPV from single model: exp(b1+b2) / (1+exp(b1+b2)), where b1 is the intercept and b2 is the coeff
id_PPV_glm_single <- exp(glm3_single_model$coefficients[1] + glm3_single_model$coefficients[2]) / 
        (1+ exp(glm3_single_model$coefficients[1] + glm3_single_model$coefficients[2]))
        

# SE(b1+b2) = sqrt(Var(b1)+Var(b2)+2Cov(b1, b2))

#Var(b1)
#vcov(summary(glm3_single_model))[1]

#Var(b2)
#vcov(summary(glm3_single_model))[4]

#Cov(b1,b2)
#vcov(summary(glm3_single_model))[2]

SE_b1b2 <- sqrt(vcov(summary(glm3_single_model))[1]+vcov(summary(glm3_single_model))[4]+2*vcov(summary(glm3_single_model))[2])
#SE_b1b2

id_PPV_lower_CI_glm_single <- (exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.025) *
                                          SE_b1b2)) / (1+(exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.025) *
                                                                      SE_b1b2)))

id_PPV_higher_CI_glm_single <- (exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.975) *
                                           SE_b1b2)) / (1+(exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.975) *
                                                                       SE_b1b2)))

# id_PPV_exp
# id_PPV_glm_single
# id_PPV_lower_CI_exp
# id_PPV_lower_CI_glm_single
# id_PPV_higher_CI_exp
# id_PPV_higher_CI_glm_single

# calculating NPV from single model: 1-exp(b1) / (1+exp(b1)), where b1 is the intercept
id_NPV_glm_single <- 1 - exp(glm3_single_model$coefficients[1]) / (1+exp(glm3_single_model$coefficients[1]))

# 95%CI 
id_NPV_higher_CI_glm_single <- 1 - exp((glm3_single_model$coefficients[1]) + qnorm(0.025) *
                                           sqrt(vcov(summary(glm3_single_model))[1])) / (1 + exp((glm3_single_model$coefficients[1]) + qnorm(0.025) *
                                                                                                        sqrt(vcov(summary(glm3_single_model))[1])))

id_NPV_lower_CI_glm_single <- 1 - exp((glm3_single_model$coefficients[1]) + qnorm(0.975) *
                                              sqrt(vcov(summary(glm3_single_model))[1])) / (1 + exp((glm3_single_model$coefficients[1]) + qnorm(0.975) *
                                                                                                            sqrt(vcov(summary(glm3_single_model))[1])))

# id_NPV_exp
# id_NPV_glm_single
# id_NPV_lower_CI_exp
# id_NPV_lower_CI_glm_single
# id_NPV_higher_CI_exp
# id_NPV_higher_CI_glm_single


##########################
## Logit transformation ##
##########################

# PPV
id_PPV_logit <- log(id_PPV/(1-id_PPV))
id_PPV_SE_logit <- sqrt(((1-id_sens)/id_sens)*(1/id_sum_D)+(id_spec/(1-id_spec))*(1/id_sum_N))

id_PPV_CI_logit <- 1.96 * id_PPV_SE_logit
id_PPV_lower_CI_logit <- id_PPV_logit - id_PPV_CI_logit
id_PPV_higher_CI_logit <- id_PPV_logit + id_PPV_CI_logit      

# exponentiating the CI values to backtransform
id_PPV_lower_CI_logit_exp <- exp(id_PPV_lower_CI_logit)/((1+exp(id_PPV_lower_CI_logit)))
id_PPV_higher_CI_logit_exp <- exp(id_PPV_higher_CI_logit)/((1+exp(id_PPV_higher_CI_logit)))

# id_PPV
# id_PPV_lower_CI_logit_exp
# id_PPV_higher_CI_logit_exp

# NPV
id_NPV_logit <- log(id_NPV/(1-id_NPV))
id_NPV_SE_logit <- sqrt((id_sens/(1-id_sens))*(1/id_sum_D)+((1-id_spec)/id_spec)*(1/id_sum_N))
id_NPV_CI_logit <- 1.96*id_NPV_SE_logit
id_NPV_lower_CI_logit <- id_NPV_logit - id_NPV_CI_logit
id_NPV_higher_CI_logit <- id_NPV_logit + id_NPV_CI_logit     

# exponentiating the CI values to backtransform
id_NPV_lower_CI_logit_exp <- exp(id_NPV_lower_CI_logit)/((1+exp(id_NPV_lower_CI_logit)))
id_NPV_higher_CI_logit_exp <- exp(id_NPV_higher_CI_logit)/((1+exp(id_NPV_higher_CI_logit)))

# id_NPV
# id_NPV_lower_CI_logit_exp
# id_NPV_higher_CI_logit_exp

#########################
# Continuity adjustment #
#########################
# suggested when sample sizes are small and PV's are close to 0 or 1

alpha <- 0.05 # risk of type I error
k <- qnorm(1-alpha/2) # (1-alpha/2)th quantile of the inverse normal distribution: with an alpha of 5%, it is approx. 1.96

# first calculate the continuity-corrected versions of sens and spec
id_sens_cont <- (id_sum_D * id_sens + k^2/2)/(id_sum_D+k^2)
id_spec_cont <- (id_sum_N * id_spec + k^2/2)/(id_sum_N+k^2)

# then calculate the continuity-corrected versions of PV's:
id_PPV_cont <- (id_sens_cont * id_prev) / ((id_sens_cont * id_prev) + (1-id_spec_cont) * (1 - id_prev)) # positive predictive value
id_NPV_cont <- (id_spec_cont * (1- id_prev)) / ((1-id_sens_cont) * id_prev + id_spec_cont * (1 - id_prev)) # negative predictive value

# calculating continuity-corrected CI's with delta method 
# PPV
id_PPV_SE_cont <- sqrt((((id_prev * (1 - id_spec_cont)) * (1 - id_prev))^2 * ((id_sens_cont * (1 - id_sens_cont))/id_sum_D) + 
                                (id_prev * id_sens_cont * (1 - id_prev))^2 * ((id_spec_cont * (1 - id_spec_cont)) / id_sum_N)) / 
                               ((id_sens_cont * id_prev + (1 - id_spec_cont) * (1 - id_prev))^4))

id_PPV_CI_cont <- 1.96 * id_PPV_SE_cont
id_PPV_lower_CI_cont <- id_PPV_cont - id_PPV_CI_cont
id_PPV_higher_CI_cont <- id_PPV_cont + id_PPV_CI_cont      

# id_PPV_cont
# id_PPV_lower_CI_cont
# id_PPV_higher_CI_cont

# NPV
id_NPV_SE_cont <- sqrt((((id_spec_cont * (1 - id_prev)) * id_prev)^2 * ((id_sens_cont * (1 - id_sens_cont)) / id_sum_D) +
                                ((1 - id_sens_cont) * (1 - id_prev) * id_prev)^2 * ((id_spec_cont * (1  - id_spec_cont)) / id_sum_N)) /
                               (((1 - id_sens_cont) * id_prev + id_spec_cont * (1 - id_prev))^4))

id_NPV_CI_cont <- 1.96 * id_NPV_SE_cont
id_NPV_lower_CI_cont <- id_NPV_cont - id_NPV_CI_cont
id_NPV_higher_CI_cont <- id_NPV_cont + id_NPV_CI_cont     

# id_NPV_cont
# id_NPV_lower_CI_cont
# id_NPV_higher_CI_cont

##################################################
# Logit transformation and continuity correction #
##################################################

# PPV
id_PPV_cont_logit <- log(id_PPV_cont/(1-id_PPV_cont))
id_PPV_SE_cont_logit <- sqrt(((1-id_sens_cont)/id_sens_cont)*(1/id_sum_D)+(id_spec_cont/(1-id_spec_cont))*(1/id_sum_N))
id_PPV_CI_cont_logit <- 1.96 * id_PPV_SE_cont_logit
id_PPV_lower_CI_cont_logit <- id_PPV_cont_logit - id_PPV_CI_cont_logit
id_PPV_higher_CI_cont_logit <- id_PPV_cont_logit + id_PPV_CI_cont_logit      

# exponentiating the CI values to backtransform
id_PPV_lower_CI_cont_logit_exp <- exp(id_PPV_lower_CI_cont_logit)/((1+exp(id_PPV_lower_CI_cont_logit)))
id_PPV_higher_CI_cont_logit_exp <- exp(id_PPV_higher_CI_cont_logit)/((1+exp(id_PPV_higher_CI_cont_logit)))

# id_PPV_cont
# id_PPV_lower_CI_cont_logit_exp
# id_PPV_higher_CI_cont_logit_exp

# NPV
id_NPV_cont_logit <- log(id_NPV_cont/(1-id_NPV_cont))
id_NPV_SE_cont_logit <- sqrt((id_sens_cont/(1-id_sens_cont))*(1/id_sum_D)+((1-id_spec_cont)/id_spec_cont)*(1/id_sum_N))
id_NPV_CI_cont_logit <- 1.96 * id_NPV_SE_cont_logit
id_NPV_lower_CI_cont_logit <- id_NPV_cont_logit - id_NPV_CI_cont_logit
id_NPV_higher_CI_cont_logit <- id_NPV_cont_logit + id_NPV_CI_cont_logit     

# exponentiating the CI values to backtransform
id_NPV_lower_CI_cont_logit_exp <- exp(id_NPV_lower_CI_cont_logit)/((1+exp(id_NPV_lower_CI_cont_logit)))
id_NPV_higher_CI_cont_logit_exp <- exp(id_NPV_higher_CI_cont_logit)/((1+exp(id_NPV_higher_CI_cont_logit)))

# id_NPV_cont
# id_NPV_lower_CI_cont_logit_exp
# id_NPV_higher_CI_cont_logit_exp



##################################################
##### Specimen/segment/sample-level analyses #####
##################################################

#######################################
### Prevalence-independent measures ###
#######################################

##################################
## No adjustment for clustering ##
##################################

#########################
## Binomial proportion ##
#########################

# sensitivity
sens <- sum_TP/(sum_TP+sum_FN) # global sensitivity

# specificity
spec <- sum_TN/(sum_FP+sum_TN) # global specificity

#######################################
# Binomial proportion, Wald intervals #
#######################################

# sensitivity
sens_SE_unadj <- sqrt((sens*(1-sens))/(sum_TP+sum_FN))
sens_CI_unadj <- 1.96*sens_SE_unadj
sens_lower_CI_unadj <- sens-sens_CI_unadj
sens_higher_CI_unadj <- sens+sens_CI_unadj   

# sens
# sens_lower_CI_unadj
# sens_higher_CI_unadj

# specificity
spec_SE_unadj <- sqrt((spec*(1-spec))/(sum_FP+sum_TN))
spec_CI_unadj <- 1.96*spec_SE_unadj
spec_lower_CI_unadj <- spec-spec_CI_unadj
spec_higher_CI_unadj <- spec+spec_CI_unadj       

# spec
# spec_lower_CI_unadj
# spec_higher_CI_unadj

##########################################################
# Binomial proportion, exact (Clopper-Pearson) intervals #
##########################################################
# as it is based on the exact binomial distribution, it is sometimes too conservative

# sensitivity
bt_sens <- binom_test(sum_TP, sum_D)
sens_exact <- bt_sens$estimate
sens_lower_CI_exact <- bt_sens$conf.low
sens_higher_CI_exact <- bt_sens$conf.high

# sens_exact
# sens_lower_CI_exact
# sens_higher_CI_exact

# specificity
bt_spec <- binom_test(sum_TN, sum_N)
spec_exact <- bt_spec$estimate
spec_lower_CI_exact <- bt_spec$conf.low
spec_higher_CI_exact <- bt_spec$conf.high

# spec_exact
# spec_lower_CI_exact
# spec_higher_CI_exact

#########################
## Logistic regression ##
#########################
# long format is required

# sensitivity
# fitting glm, -1 suppresses intercept
glm1 <- glm(Test ~ Disease - 1, family = binomial(link="logit"), data=data.long)
#summary(glm1)

# converting estimate and CI's to probability scale (exp/(1+exp))
sens_exp <- exp(summary(glm1)$coefficients[,1])  / (1+exp(summary(glm1)$coefficients[,1] ))

sens_lower_CI_exp <- exp(summary(glm1)$coefficients["Disease",1] + qnorm(0.025) *
                                  summary(glm1)$coefficients["Disease",2]) / (1+exp(summary(glm1)$coefficients["Disease",1] + qnorm(0.025) *
                                                                                             summary(glm1)$coefficients["Disease",2]))

sens_higher_CI_exp <- exp(summary(glm1)$coefficients["Disease",1] + qnorm(0.975) *
                                   summary(glm1)$coefficients["Disease",2]) / (1+exp(summary(glm1)$coefficients["Disease",1] + qnorm(0.975) *
                                                                                              summary(glm1)$coefficients["Disease",2]))
# sens_exp
# sens_lower_CI_exp
# sens_higher_CI_exp

# specificity
# fitting glm
glm2 <- glm(Inv_test ~ Inv_disease - 1, family = binomial(link="logit"), data=data.long)
#summary(glm2)

# converting estimate and CI's to probability scale (exp/(1+exp))
spec_exp <- (exp(summary(glm2)$coefficients["Inv_disease",1] + qnorm(0.5) *
                         summary(glm2)$coefficients["Inv_disease",2]) / (1+exp(summary(glm2)$coefficients["Inv_disease",1] + qnorm(0.5) *
                                                                                        summary(glm2)$coefficients["Inv_disease",2])))

spec_lower_CI_exp <- (exp(summary(glm2)$coefficients["Inv_disease",1] + qnorm(0.025) *
                                  summary(glm2)$coefficients["Inv_disease",2]) / (1+exp(summary(glm2)$coefficients["Inv_disease",1] + qnorm(0.025) *
                                                                                                 summary(glm2)$coefficients["Inv_disease",2])))

spec_higher_CI_exp <- (exp(summary(glm2)$coefficients["Inv_disease",1] + qnorm(0.975) *
                                   summary(glm2)$coefficients["Inv_disease",2]) / (1+exp(summary(glm2)$coefficients["Inv_disease",1] + qnorm(0.975) *
                                                                                                  summary(glm2)$coefficients["Inv_disease",2])))


#####################################
# Logistic regression, single model #
#####################################

# modeling sens and spec in a single model (e.g., Ronco & Biggeri, 1999)
glm1_single_model <- glm(Test ~ Disease, family=binomial(link = "logit"), data=data.long)  
#summary(glm1_single_model)


# calculating sensitivity from single model: exp(b1+b2) / (1+exp(b1+b2)), where b1 is the intercept and b2 is the coeff
sens_glm_single <- (exp(summary(glm1_single_model)$coefficients[[1]] + summary(glm1_single_model)$coefficients[[2]])) / 
        (1+(exp(summary(glm1_single_model)$coefficients[[1]] + summary(glm1_single_model)$coefficients[[2]])))

#names(fixef(glm1_single_model))[1]

# calculating 95% CI for sensitivity from single model: exp(L1+L2) / (1+exp(L1+L2)), where L1 and L2 are the lower bounds of 
# the intercept and coeff
# lincom to calculate 95% CI for sensitivity
#summary(glht(glm1_single_model, linfct = c("`(Intercept)` + Disease = 0")))

# Show the confidence interval
glm1_single_model.lh <- glht(glm1_single_model, linfct = c("`(Intercept)` + Disease = 0"))
a <- confint(glm1_single_model.lh)
# lower CI
#a$confint[[2]]
# higher CI
#a$confint[[3]]

sens_lower_CI_glm_single <- (exp(a$confint[[2]])) / 
        (1+exp(a$confint[[2]]))

sens_higher_CI_glm_single <- (exp(a$confint[[3]])) / 
        (1+exp(a$confint[[3]]))


# calculating specificity from single model: 1-exp(b1) / (1+exp(b1)), where b1 is the intercept
spec_glm_single <- 1-exp(summary(glm1_single_model)$coefficients[[1]]) / (1+exp(summary(glm1_single_model)$coefficients[[1]]))


spec_higher_CI_glm_single <- 1- (exp(summary(glm1_single_model)$coefficients[[1]] + qnorm(0.025) *
                                               summary(glm1_single_model)$coefficients[[3]])) / (1+exp(summary(glm1_single_model)$coefficients[[1]] + qnorm(0.025) *
                                                                                                                 summary(glm1_single_model)$coefficients[[3]]))

spec_lower_CI_glm_single <- 1- (exp(summary(glm1_single_model)$coefficients[[1]] + qnorm(0.975) *
                                              summary(glm1_single_model)$coefficients[[3]])) / (1+exp(summary(glm1_single_model)$coefficients[[1]] + qnorm(0.975) *
                                                                                                                summary(glm1_single_model)$coefficients[[3]]))

# sens_exp
# sens_glm_single
# sens_lower_CI_exp
# sens_lower_CI_glm_single
# sens_higher_CI_exp
# sens_higher_CI_glm_single
# 
# 
# spec_exp
# spec_glm_single
# spec_lower_CI_exp
# spec_lower_CI_glm_single
# spec_higher_CI_exp
# spec_higher_CI_glm_single


###############################
## Adjustment for clustering ##
###############################

############################
# Preliminary calculations #
############################

# calculations for ratio estimation adjustment
mean_cluster_size <- sum_j / id_sum # mean cluster size
j_Dmean <- sum_D / id_sum_D # mean cluster size, diseased patients
j_Nmean <- sum_N / id_sum_N # mean cluster size, non-diseased patients

sum_Nw <- sum(data$Nw) # weighted mean cluster size
sum_Nw_sens <- sum(data$Nw_sens) # weighted mean cluster size, diseased patients
sum_Nw_spec <- sum(data$Nw_spec) # weighted mean cluster size, non-diseased patients

# natural estimator calculations for ICC (separately for sens and spec) from Genders et al., 2012 <- Fleiss et al., 2003, p. 443
data$num_ICC_sens <- data$TP * (data$TP - 1) -2 * sens * (data$D - 1) * data$TP + data$D * (data$D - 1) * sens * sens 
data$den_ICC_sens <- data$D * (data$D - 1) * sens * (1- sens)

sum_num_ICC_sens <- sum(data$num_ICC_sens)
sum_den_ICC_sens <- sum(data$den_ICC_sens)

ICC_sens <- sum_num_ICC_sens / sum_den_ICC_sens # intraclass correlation for sens
#ICC_sens

data$num_ICC_spec <- data$TN * (data$TN - 1) -2 * spec * (data$N - 1) * data$TN + data$N * (data$N - 1) * spec * spec 
data$den_ICC_spec <- data$N * (data$N - 1) * spec * (1- spec)

sum_num_ICC_spec <- sum(data$num_ICC_spec)
sum_den_ICC_spec <- sum(data$den_ICC_spec)

ICC_spec <- sum_num_ICC_spec / sum_den_ICC_spec # intraclass correlation for spec
ICC_spec

VIF_sens <- 1 + (sum_Nw_sens - 1) * ICC_sens # variance inflation factor for sensitivity
VIF_spec <- 1 + (sum_Nw_spec - 1) * ICC_spec # variance inflation factor for specificity

data$D_sens <- data$D * data$sens 
data$N_spec <- data$N * data$spec  
data$D_j_Dmean_2 <- (data$D/j_Dmean)^2
data$N_j_Nmean_2 <- (data$N/j_Nmean)^2
data$ind_sens_sens_2 <- (data$sens - sens)^2
data$ind_spec_spec_2 <- (data$spec - spec)^2
data$rt <- data$D_j_Dmean_2 * data$ind_sens_sens_2
data$su <- data$N_j_Nmean_2 * data$ind_spec_spec_2

sum_rt <- sum(data$rt, na.rm=T)
sum_su <- sum(data$su, na.rm=T)

##############################
# Ratio estimator adjustment #
##############################

# sensitivity
sens_SE_readj <- sqrt((1/(id_sum_D*(id_sum_D-1)))*(sum_rt))
sens_CI_readj <- 1.96 * sens_SE_readj
sens_lower_CI_readj <- sens - sens_CI_readj
sens_higher_CI_readj <- sens + sens_CI_readj       

# specificity
spec_SE_readj <- sqrt((1/(id_sum_N*(id_sum_N-1)))*(sum_su))
spec_CI_readj <- 1.96 * spec_SE_readj
spec_lower_CI_readj <- spec - spec_CI_readj
spec_higher_CI_readj <- spec + spec_CI_readj       

########################################
# Variance inflation factor adjustment #
########################################

# sensitivity
sens_SE_vifadj <- sens_SE_unadj * sqrt(VIF_sens)
sens_CI_vifadj <- 1.96 * sens_SE_vifadj
sens_lower_CI_vifadj <- sens - sens_CI_vifadj
sens_higher_CI_vifadj <- sens + sens_CI_vifadj       

# specificity
spec_SE_vifadj <- spec_SE_unadj * sqrt(VIF_spec)
spec_CI_vifadj <- 1.96 * spec_SE_vifadj
spec_lower_CI_vifadj <- spec - spec_CI_vifadj
spec_higher_CI_vifadj <- spec + spec_CI_vifadj       

######################################
# Logistic regression with robust SE #
######################################
# sandwich estimator, c.f., Williams (2000), implementation aided by McDonald (2019)

# sens
# general sandwich estimator: no regard for the type of cluster (in stata: vce (robust))
# when we have no information about how the data is clustered together, but we know they are clustered
sens_lower_CI_robust_logit <- coefci(glm1, vcov. = vcovHC(glm1, type = 'HC1'))[1] 
sens_higher_CI_robust_logit <- coefci(glm1, vcov. = vcovHC(glm1, type = 'HC1'))[2]

# exponentiating to backtransform
sens_lower_CI_robust <- exp(sens_lower_CI_robust_logit) / (1+exp(sens_lower_CI_robust_logit)) 
sens_higher_CI_robust <- exp(sens_higher_CI_robust_logit) / (1+exp(sens_higher_CI_robust_logit))

# we need a sandwich estimator that clusters on id (in stata: vce (cluster id))
sens_lower_CI_robust_cl_logit <- coefci(glm1, vcov. = vcovCL(glm1, cluster = ~id, type = 'HC1'))[1] 
sens_higher_CI_robust_cl_logit <- coefci(glm1,  vcov. = vcovCL(glm1, cluster = ~id, type = 'HC1'))[2]

# exponentiating to backtransform
sens_lower_CI_robust_cl <- exp(sens_lower_CI_robust_cl_logit) / (1+exp(sens_lower_CI_robust_cl_logit)) 
sens_higher_CI_robust_cl <- exp(sens_higher_CI_robust_cl_logit) / (1+exp(sens_higher_CI_robust_cl_logit))

# spec
# general sandwich estimator: no regard for the type of cluster (in stata: vce (robust))
# when we have no information about how the data is clustered together, but we know they are clustered
spec_lower_CI_robust_logit <- coefci(glm2, vcov. = vcovHC(glm2, type = 'HC1'))[1] 
spec_higher_CI_robust_logit <- coefci(glm2, vcov. = vcovHC(glm2, type = 'HC1'))[2]

# exponentiating to backtransform
spec_lower_CI_robust <- exp(spec_lower_CI_robust_logit) / (1+exp(spec_lower_CI_robust_logit)) # 
spec_higher_CI_robust <- exp(spec_higher_CI_robust_logit) / (1+exp(spec_higher_CI_robust_logit))

# we need a sandwich estimator that clusters on id (in stata: vce (cluster id))
spec_lower_CI_robust_cl_logit <- coefci(glm2, vcov. = vcovCL(glm2, cluster = ~id, type = 'HC1'))[1] 
spec_higher_CI_robust_cl_logit <- coefci(glm2,  vcov. = vcovCL(glm2, cluster = ~id, type = 'HC1'))[2]

# exponentiating to backtransform
spec_lower_CI_robust_cl <- exp(spec_lower_CI_robust_cl_logit) / (1+exp(spec_lower_CI_robust_cl_logit)) 
spec_higher_CI_robust_cl <- exp(spec_higher_CI_robust_cl_logit) / (1+exp(spec_higher_CI_robust_cl_logit))

# spec_lower_CI_robust
# spec_lower_CI_robust_cl
# spec_higher_CI_robust
# spec_higher_CI_robust_cl

####################################################
# Logistic regression with robust SE, single model #
####################################################
# sandwich estimator, c.f., Williams (2000), implementation aided by McDonald (2019)

# sens
# general sandwich estimator: no regard for the type of cluster (in stata: vce (robust))
# when we have no information about how the data is clustered together, but we know they are clustered

# SE(b1+b2) = sqrt(Var(b1)+Var(b2)+2Cov(b1+b2))

# adjusted var(b1)
#vcovHC(glm1_single_model, type = 'HC1')[1]

# adjusted var(b2)
#vcovHC(glm1_single_model, type = 'HC1')[4]

# adjusted Cov(b1 ,b2)
#vcovHC(glm1_single_model, type = 'HC1')[2]

SE_b1b2 <- sqrt(vcovHC(glm1_single_model, type = 'HC1')[1]+vcovHC(glm1_single_model, type = 'HC1')[4]+2*vcovHC(glm1_single_model, type = 'HC1')[2])

sens_lower_CI_robust_single <- (exp((glm1_single_model$coefficients[1]+glm1_single_model$coefficients[2]) + qnorm(0.025) *
                                        SE_b1b2)) / (1+(exp((glm1_single_model$coefficients[1]+glm1_single_model$coefficients[2]) + qnorm(0.025) *
                                                                    SE_b1b2)))

sens_higher_CI_robust_single <- (exp((glm1_single_model$coefficients[1]+glm1_single_model$coefficients[2]) + qnorm(0.975) *
                                         SE_b1b2)) / (1+(exp((glm1_single_model$coefficients[1]+glm1_single_model$coefficients[2]) + qnorm(0.975) *
                                                                     SE_b1b2)))
# sens_lower_CI_robust
# sens_lower_CI_robust_single
# sens_higher_CI_robust
# sens_higher_CI_robust_single

# spec 95%CI 
spec_higher_CI_robust_single <- 1 - exp((glm1_single_model$coefficients[1]) + qnorm(0.025) *
                                                sqrt(vcovHC(glm1_single_model, type = 'HC1')[1])) / (1 + exp((glm1_single_model$coefficients[1]) + qnorm(0.025) *
                                                                                                               sqrt(vcovHC(glm1_single_model, type = 'HC1')[1])))

spec_lower_CI_robust_single <- 1 - exp((glm1_single_model$coefficients[1]) + qnorm(0.975) *
                                               sqrt(vcovHC(glm1_single_model, type = 'HC1')[1])) / (1 + exp((glm1_single_model$coefficients[1]) + qnorm(0.975) *
                                                                                                                    sqrt(vcovHC(glm1_single_model, type = 'HC1')[1])))

# spec_lower_CI_robust
# spec_lower_CI_robust_single
# spec_higher_CI_robust
# spec_higher_CI_robust_single

# we need a sandwich estimator that clusters on id (in stata: vce (cluster id))
# SE(b1+b2) = sqrt(Var(b1)+Var(b2)+2Cov(b1+b2))
# SE_b1b2

# adjusted var(b1)
#vcovHC(glm1_single_model, cluster= ~id, type = 'HC1')[1]

# adjusted var(b2)
#vcovHC(glm1_single_model, cluster= ~id, type = 'HC1')[4]

# adjusted Cov(b1 ,b2)
#vcovHC(glm1_single_model, cluster= ~id, type = 'HC1')[2]

# SE_b1b2
SE_b1b2 <- sqrt(vcovHC(glm1_single_model, cluster= ~id, type = 'HC1')[1]+vcovHC(glm1_single_model, cluster= ~id, type = 'HC1')[4]+2*vcovHC(glm1_single_model, cluster= ~id, type = 'HC1')[2])
#SE_b1b2

sens_lower_CI_robust_cl_single <- (exp((glm1_single_model$coefficients[1]+glm1_single_model$coefficients[2]) + qnorm(0.025) *
                                            SE_b1b2)) / (1+(exp((glm1_single_model$coefficients[1]+glm1_single_model$coefficients[2]) + qnorm(0.025) *
                                                                        SE_b1b2)))

sens_higher_CI_robust_cl_single <- (exp((glm1_single_model$coefficients[1]+glm1_single_model$coefficients[2]) + qnorm(0.975) *
                                             SE_b1b2)) / (1+(exp((glm1_single_model$coefficients[1]+glm1_single_model$coefficients[2]) + qnorm(0.975) *
                                                                         SE_b1b2)))

# sens_lower_CI_robust_cl
# sens_lower_CI_robust_cl_single
# sens_higher_CI_robust_cl
# sens_higher_CI_robust_cl_single

# spec 95%CI robust cluster
spec_higher_CI_robust_cl_single <- 1 - exp((glm1_single_model$coefficients[1]) + qnorm(0.025) *
                                                sqrt(vcovHC(glm1_single_model, cluster= ~id, type = 'HC1')[1])) / (1 + exp((glm1_single_model$coefficients[1]) + qnorm(0.025) *
                                                                                                                     sqrt(vcovHC(glm1_single_model, cluster= ~id, type = 'HC1')[1])))

spec_lower_CI_robust_cl_single <- 1 - exp((glm1_single_model$coefficients[1]) + qnorm(0.975) *
                                               sqrt(vcovHC(glm1_single_model, cluster= ~id, type = 'HC1')[1])) / (1 + exp((glm1_single_model$coefficients[1]) + qnorm(0.975) *
                                                                                                                    sqrt(vcovHC(glm1_single_model, cluster= ~id, type = 'HC1')[1])))

# spec_lower_CI_robust_cl
# spec_lower_CI_robust_cl_single
# spec_higher_CI_robust_cl
# spec_higher_CI_robust_cl_single


#####################################
# Mixed effects logistic regression #
#####################################
# id in random effect structure as random intercept

# sensitivity
# NB: syntax for no intercept is different in glmer compared to glm
# nAGQ means 'number of adaptive Gauss-Hermite quadrature points', and sets how glmer will integrate out the random 
# effects when fitting the mixed model. When nAGQ is greater than 1, then adaptive quadrature is used with nAGQ points. When nAGQ =
#1, the Laplace approximation is used, and when nAGQ = 0, the integral is 'ignored'.
# 'glmer' uses the Laplacian approximation as default, corresponding to adaptive Gauss-Hermite approximation 
# with only 1 point, while 'xtmelogit' uses 7 points.
glmer1 <-  (glmer(Test ~ 1 + (1 | id), family=binomial(link = "logit"), nAGQ = 1, data=data.long, subset=Disease==1)) 
#summary(glmer1)

#(vc <- VarCorr(glmer1))  ## default print method: standard dev and corr
## both variance and std.dev.
#print(vc,comp=c("Variance","Std.Dev."),digits=2)
## variance only
#print(vc,comp=c("Variance"))
#as.data.frame(vc)
#as.data.frame(vc,order="lower.tri")
#icc(glmer1)

# converting estimate and CI's to probability scale (exp/(1+exp))
sens_glmer <-  ((exp(summary(glmer1)$coefficients[,1] + qnorm(0.5) *
                                        summary(glmer1)$coefficients[,2]) / (1+exp(summary(glmer1)$coefficients[,1] + qnorm(0.5) *
                                                                                                     summary(glmer1)$coefficients[,2]))))

sens_lower_CI_glmer <-  ((exp(summary(glmer1)$coefficients[,1] + qnorm(0.025) *
                                        summary(glmer1)$coefficients[,2]) / (1+exp(summary(glmer1)$coefficients[,1] + qnorm(0.025) *
                                                                                                    summary(glmer1)$coefficients[,2]))))

sens_higher_CI_glmer <-  ((exp(summary(glmer1)$coefficients[,1] + qnorm(0.975) *
                                         summary(glmer1)$coefficients[,2]) / (1+exp(summary(glmer1)$coefficients[,1] + qnorm(0.975) *
                                                                                                     summary(glmer1)$coefficients[,2]))))

# sens_glmer
# sens_lower_CI_glmer
# sens_higher_CI_glmer


# specificity: 
#glmer2 <-  glmer(Inv_test ~ 1 + (1 | id), family=binomial(link = "logit"), nAGQ = 1, data=data.long, subset=Inv_disease==1) 
# convergence warning, probably due to low number of patients in one condition: changing optimizer
glmer2 <-  glmer(Inv_test ~ 1 + (1 | id), family=binomial(link = "logit"), nAGQ = 1, data=data.long, subset=Inv_disease==1, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
#summary(glmer2)

# converting estimate and CI's to probability scale (exp/(1+exp))
spec_glmer <-  ((exp(summary(glmer2)$coefficients[,1] + qnorm(0.5) *
                           summary(glmer2)$coefficients[,2]) / (1+exp(summary(glmer2)$coefficients[,1] + qnorm(0.5) *
                                                                                           summary(glmer2)$coefficients[,2]))))

spec_lower_CI_glmer <-  ((exp(summary(glmer2)$coefficients[,1] + qnorm(0.025) *
                                    summary(glmer2)$coefficients[,2]) / (1+exp(summary(glmer2)$coefficients[,1] + qnorm(0.025) *
                                                                                                    summary(glmer2)$coefficients[,2]))))

spec_higher_CI_glmer <-  ((exp(summary(glmer2)$coefficients[,1] + qnorm(0.975) *
                                     summary(glmer2)$coefficients[,2]) / (1+exp(summary(glmer2)$coefficients[,1] + qnorm(0.975) *
                                                                                                     summary(glmer2)$coefficients[,2]))))

# spec_glmer
# spec_lower_CI_glmer
# spec_higher_CI_glmer

###################################################
# Mixed effects logistic regression, single model #
###################################################

# modeling sens and spec in a single model (e.g., Ronco & Biggeri, 1999)
glmer1_single_model <-  (glmer(Test ~ Disease + (1 | id), family=binomial(link = "logit"), nAGQ = 1, data=data.long))  
#summary(glmer1_single_model)

#icc(glmer1_single_model)
# intercept coef
#summary(glmer1_single_model)$coef[[1]]

# disease coef
#summary(glmer1_single_model)$coef[[2]]

# intercept SE
#summary(glmer1_single_model)$coef[[3]]

# disease SE
#summary(glmer1_single_model)$coef[[4]]

#exp(summary(glmer1_single_model)$coef[[1]])

# calculating sensitivity from single model: exp(b1+b2) / (1+exp(b1+b2)), where b1 is the intercept and b2 is the coeff
sens_glmer_single <-  ((exp(summary(glmer1_single_model)$coefficients[[1]] + summary(glmer1_single_model)$coefficients[[2]])) / 
        (1+(exp(summary(glmer1_single_model)$coefficients[[1]] + summary(glmer1_single_model)$coefficients[[2]]))))

#names(fixef(glmer1_single_model))[1]

# calculating 95% CI for sensitivity from single model: exp(L1+L2) / (1+exp(L1+L2)), where L1 and L2 are the lower bounds of 
# the intercept and coeff
# lincom to calculate 95% CI for sensitivity
#summary(glht(glmer1_single_model, linfct = c("`(Intercept)` + Disease = 0")))

# Show the confidence interval
glmer1_single_model.lh <-  (glht(glmer1_single_model, linfct = c("`(Intercept)` + Disease = 0")))
a <-  (confint(glmer1_single_model.lh))
# lower CI
#a$confint[[2]]
# higher CI
#a$confint[[3]]

sens_lower_CI_glmer_single <-  ((exp(a$confint[[2]])) / 
        (1+exp(a$confint[[2]])))

sens_higher_CI_glmer_single <-  ((exp(a$confint[[3]])) / 
        (1+exp(a$confint[[3]])))

# another method
# SE(b1+b2) = sqrt(Var(b1)+Var(b2)+2Cov(b1+b2))
# SE_b1b2 <- sqrt(vcov(glmer1_single_model)[1]+vcov(glmer1_single_model)[4]+2*(vcov(glmer1_single_model)[2]))
# SE_b1b2

# calculating specificity from single model: 1-exp(b1) / (1+exp(b1)), where b1 is the intercept
spec_glmer_single <-  (1-exp(summary(glmer1_single_model)$coefficients[[1]]) / (1+exp(summary(glmer1_single_model)$coefficients[[1]])))


spec_higher_CI_glmer_single <-  (1- (exp(summary(glmer1_single_model)$coefficients[[1]] + qnorm(0.025) *
                                         summary(glmer1_single_model)$coefficients[[3]])) / (1+exp(summary(glmer1_single_model)$coefficients[[1]] + qnorm(0.025) *
                                                                                                    summary(glmer1_single_model)$coefficients[[3]])))

spec_lower_CI_glmer_single <-  (1- (exp(summary(glmer1_single_model)$coefficients[[1]] + qnorm(0.975) *
                                        summary(glmer1_single_model)$coefficients[[3]])) / (1+exp(summary(glmer1_single_model)$coefficients[[1]] + qnorm(0.975) *
                                                                                                   summary(glmer1_single_model)$coefficients[[3]])))

# calculating higher 95% CI of the spec from CI of b1, the intercept from exp confint: similar to model, but using profile CI's
#spec_higher_CI_glmer_single <- 1-exp(confint(glmer1_single_model))[[2]] / (1+exp(confint(glmer1_single_model))[[2]])


# calculating lower 95% CI of the spec from CI of b1, the intercept from exp confint
#spec_lower_CI_glmer_single <- 1-exp(confint(glmer1_single_model))[[5]] / (1+exp(confint(glmer1_single_model))[[5]]) 

# sens_glmer
# sens_glmer_single
# 
# sens_lower_CI_glmer
# sens_lower_CI_glmer_single
# sens_lower_CI_exact
# 
# sens_higher_CI_glmer
# sens_higher_CI_glmer_single
# sens_higher_CI_exact
# 
# spec_glmer
# spec_glmer_single
# 
# spec_lower_CI_glmer
# spec_lower_CI_glmer_single
# spec_lower_CI_exact
# 
# spec_higher_CI_glmer
# spec_higher_CI_glmer_single
# spec_higher_CI_exact


# exp(confint(glmer1_single_model))[[5]]
# 
# 
# estimates <- round(coef(glm1), 3)
# 
# results_tab <- tidy(glm1) %>%
#         mutate_if(is.numeric, funs(round(.,3))) %>%
#         var_mapping(term) %>%
#         dplyr::rename(Estimate = term, Beta = estimate, SE = std.error, z = statistic, p = p.value)
# 
# 
# #A function to return odds ratios and confidence intervals based on normal approximations.
# 
# OR_SE <- function(glm_object, variable){
#         
#         d <- exp(coef(glm_object)[variable])
#         v <- vcov(glm_object)[variable, variable]
#         
#         output_df <- as_tibble(list(Estimate = variable,
#                                     Beta     = coef(glm_object)[variable],
#                                     OR       = d,
#                                     `OR SE`  = sqrt(d*v*d))) %>%
#                 mutate(`Lower CI` = OR - 1.96*`OR SE`,
#                        `Upper CI` = OR + 1.96*`OR SE`)
#         
#         output_df
#         
# }
# #Map over all estimates and reduce to tibble.
# 
# or_tibble <- map(names(coef(glm1)), function(i) OR_SE(glm1, i)) %>%
#         reduce(bind_rows)
# kable(or_tibble %>% var_mapping(Estimate), align = c("l", rep("c", 5)))
# 
# 
# library(knitr)
# kable(results_tab, align = c("l",rep("c",4)))
# 
# 
# spec_CI_single_  <- 1-exp(summary(glmer1_single_model)$coefficients[[3]]^2) / (1+exp(summary(glmer1_single_model)$coefficients[[3]]^2))        
# vcov(glmer1_single_model)[1]
# vcov(glmer1_single_model)[2]
# vcov(glmer1_single_model)[3]
# vcov(glmer1_single_model)[4]
# 
# 
# summary(glmer1_single_model)
# #spec_lower_CI_single <- spec_single - spec_CI_single
#         


##########################################
# Generalized estimating equations (GEE) #
##########################################

# test result (+/-) is modeled as the outcome variable, disease status (from a reference test,
# gold standard) as the predictor, with a logit link

# NB: data needs to be sorted by id first (if not already)
# sens 
# in order to ignore the error message "Cgee: error: logistic model for probability has fitted value very close to 1.
# estimates diverging; iteration terminated."
# we wrap statement in try so that the error message will be displayed but execution will continue:

# produce NA when no FN patients are available
# this way, the gee1 model runs, otherwise it produces fatal error
# depending on your data structure, you might need to extend the conditional replacement to capture no or small number of subjects per cell
# in the patient-level contingency table with the other gee models as well

#if (id_no_FN != 0) {
#        sens_gee <- NA
#        sens_lower_CI_gee <- NA
#        sens_higher_CI_gee <- NA
#} else {
        gee1 <-  (gee(Test ~ 1, data = data.long, subset=Disease==1, id=id, family=binomial(link = "logit"), corstr="exchangeable"))
        #alternative correlation matrix
        #gee1 <- gee(Test ~ 1, data = data.long, subset=Disease==1, id=id, family=binomial(link = "logit"), corstr="unstructured")
        #summary(gee1)

        # converting estimate and CI's to probability scale (exp/(1+exp)) using robust SE's (column 4)
        sens_gee <-  (exp(summary(gee1)$coefficients[,1]) / (1+exp(summary(gee1)$coefficients[,1])))
        
        sens_lower_CI_gee <-  (exp(summary(gee1)$coefficients[,1] + qnorm(0.025) *
                                             summary(gee1)$coefficients[,4]) / (1+exp(summary(gee1)$coefficients[,1] + qnorm(0.025) *
                                                                                              summary(gee1)$coefficients[,4])))
        
        sens_higher_CI_gee <-  (exp(summary(gee1)$coefficients[,1] + qnorm(0.975) *
                                              summary(gee1)$coefficients[,4]) / (1+exp(summary(gee1)$coefficients[,1] + qnorm(0.975) *
                                                                                               summary(gee1)$coefficients[,4])))
        
#}

# sens_gee
# sens_lower_CI_gee
# sens_higher_CI_gee

# spec
gee2 <-  (gee(Inv_test ~ 1, data = data.long, subset=Inv_disease==1, id=id, family=binomial(link = "logit"), corstr="exchangeable"))
#gee2 <- gee(Inv_test ~ 1, data = data.long, subset=Inv_disease==1, id=id, family=binomial(link = "logit"), corstr="unstructured")
#summary(gee2)

# converting estimate and CI's to probability scale (exp/(1+exp)) using robust SE's (column 4)
spec_gee <-  (exp(summary(gee2)$coefficients[,1]) / (1+exp(summary(gee2)$coefficients[,1])))

spec_lower_CI_gee <-  (exp(summary(gee2)$coefficients[,1] + qnorm(0.025) *
                                  summary(gee2)$coefficients[,4]) / (1+exp(summary(gee2)$coefficients[,1] + qnorm(0.025) *
                                                                                   summary(gee2)$coefficients[,4])))

spec_higher_CI_gee <-  (exp(summary(gee2)$coefficients[,1] + qnorm(0.975) *
                                   summary(gee2)$coefficients[,4]) / (1+exp(summary(gee2)$coefficients[,1] + qnorm(0.975) *
                                                                                    summary(gee2)$coefficients[,4])))

# spec_gee
# spec_lower_CI_gee
# spec_higher_CI_gee

########################################################
# Generalized estimating equations (GEE), single model #
########################################################

gee1_single_model <-  (gee(Test ~ Disease, data = data.long, id=id, family=binomial(link = "logit"), corstr="exchangeable"))
#gee1_single_model <- gee(Test ~ Disease, data = data.long, id=id, family=binomial(link = "logit"), corstr="unstructured")
#summary(gee1_single_model)

#coef(summary(gee1_single_model))

# calculating sensitivity from single model: exp(b1+b2) / (1+exp(b1+b2)), where b1 is the intercept and b2 is the coeff
sens_gee_single <-  ((exp(summary(gee1_single_model)$coefficients[[1]] + summary(gee1_single_model)$coefficients[[2]])) / 
        (1+(exp(summary(gee1_single_model)$coefficients[[1]] + summary(gee1_single_model)$coefficients[[2]]))))


# gee is not compatible neither with coefci nor with glht, using a different method to calculate 95% CI

# calculating 95% CI for sensitivity from single model: exp(L1+L2) / (1+exp(L1+L2)), where L1 and L2 are the lower bounds of 

# handcoding linear combinations to obtain the 95% CI of the combinations of the intercept and coef
# SE_b1b2 <- sqrt(vcov(glmer1_single_model)[1]+vcov(glmer1_single_model)[4]+2*(vcov(glmer1_single_model)[2]))
# SE_b1b2

#Var(b1)
#gee1_single_model$robust.variance[1]

#Var(b2)
#gee1_single_model$robust.variance[4]

#Cov(b1,b2)
#gee1_single_model$robust.variance[2]

SE_b1b2 <-  (sqrt(gee1_single_model$robust.variance[1]+gee1_single_model$robust.variance[4]+2*gee1_single_model$robust.variance[2]))
#SE_b1b2

sens_lower_CI_gee_single <-  ((exp((summary(gee1_single_model)$coefficients[[1]] + summary(gee1_single_model)$coefficients[[2]]) + qnorm(0.025) *
                                                  SE_b1b2)) / (1+(exp((summary(gee1_single_model)$coefficients[[1]] + summary(gee1_single_model)$coefficients[[2]]) + qnorm(0.025) *
                                                                              SE_b1b2))))

sens_higher_CI_gee_single <-  ((exp((summary(gee1_single_model)$coefficients[[1]] + summary(gee1_single_model)$coefficients[[2]]) + qnorm(0.975) *
                                           SE_b1b2)) / (1+(exp((summary(gee1_single_model)$coefficients[[1]] + summary(gee1_single_model)$coefficients[[2]]) + qnorm(0.975) *
                                                                       SE_b1b2))))


# calculating specificity from single model: 1-exp(b1) / (1+exp(b1)), where b1 is the intercept
spec_gee_single <-  (1-exp(summary(gee1_single_model)$coefficients[[1]]) / (1+exp(summary(gee1_single_model)$coefficients[[1]])))


spec_higher_CI_gee_single <-  (1- (exp(summary(gee1_single_model)$coefficients[[1]] + qnorm(0.025) *
                                               summary(gee1_single_model)$coefficients[[7]])) / (1+exp(summary(gee1_single_model)$coefficients[[1]] + qnorm(0.025) *
                                                                                                                 summary(gee1_single_model)$coefficients[[7]])))

spec_lower_CI_gee_single <-  (1- (exp(summary(gee1_single_model)$coefficients[[1]] + qnorm(0.975) *
                                              summary(gee1_single_model)$coefficients[[7]])) / (1+exp(summary(gee1_single_model)$coefficients[[1]] + qnorm(0.975) *
                                                                                                                summary(gee1_single_model)$coefficients[[7]])))

# sens_gee
# sens_gee_single
# sens_lower_CI_gee
# sens_lower_CI_gee_single
# sens_higher_CI_gee
# sens_higher_CI_gee_single

# spec_gee
# spec_gee_single
# spec_lower_CI_gee
# spec_lower_CI_gee_single
# spec_higher_CI_gee
# spec_higher_CI_gee_single

#####################
# Cluster bootstrap #
#####################

# sens
# you may decrease bootstrapping iterations (B) if processing takes too long
clusboot1 <- clusbootglm(Test ~ Disease - 1, clusterid=id, data=data.long, family=binomial, B = 10000, confint.level = 0.95)
#summary(clusboot1)
sens_boot <- exp(clusboot1$boot.coefs)  / (1+exp(clusboot1$boot.coefs))
sens_lower_CI_boot <- exp(clusboot1$percentile.interval[1])  / (1+exp(clusboot1$percentile.interval[1]))
sens_higher_CI_boot <- exp(clusboot1$percentile.interval[2])  / (1+exp(clusboot1$percentile.interval[2]))

# sens_boot
# sens_lower_CI_boot
# sens_higher_CI_boot

# spec
clusboot2 <- clusbootglm(Inv_test ~ Inv_disease - 1, clusterid=id, data=data.long, family=binomial, B = 10000, confint.level = 0.95)

spec_boot <- exp(clusboot2$boot.coefs)  / (1+exp(clusboot2$boot.coefs))
spec_lower_CI_boot <- exp(clusboot2$percentile.interval[1])  / (1+exp(clusboot2$percentile.interval[1]))
spec_higher_CI_boot <- exp(clusboot2$percentile.interval[2])  / (1+exp(clusboot2$percentile.interval[2]))

# spec_boot
# spec_lower_CI_boot
# spec_higher_CI_boot


###################################
# Cluster bootstrap, single model #
###################################

# sens
clusboot1_single_model <- clusbootglm(Test ~ Disease, clusterid=id, data=data.long, family=binomial, B = 10000, confint.level = 0.95)

# calculating sensitivity from single model: exp(b1+b2) / (1+exp(b1+b2)), where b1 is the intercept and b2 is the coeff
sens_boot_single <- exp(clusboot1_single_model$boot.coefs[1]+clusboot1_single_model$boot.coefs[2])  / (1+exp(clusboot1_single_model$boot.coefs[1]+clusboot1_single_model$boot.coefs[2]))

# SE(b1+b2) = sqrt(Var(b1)+Var(b2)+2Cov(b1, b2))
# SE_b1b2 <- sqrt(vcov(glmer1_single_model)[1]+vcov(glmer1_single_model)[4]+2*(vcov(glmer1_single_model)[2]))
# SE_b1b2

#Var(b1)
#cov(clusboot1_single_model$coefficients)[1]

#Var(b2)
#cov(clusboot1_single_model$coefficients)[4]

#Cov(b1,b2)
#cov(clusboot1_single_model$coefficients)[2]

SE_b1b2 <- sqrt(cov(clusboot1_single_model$coefficients)[1]+cov(clusboot1_single_model$coefficients)[4]+2*cov(clusboot1_single_model$coefficients)[2])
#SE_b1b2

sens_lower_CI_boot_single <- (exp((clusboot1_single_model$boot.coefs[1]+clusboot1_single_model$boot.coefs[2]) + qnorm(0.025) *
                                             SE_b1b2)) / (1+(exp((clusboot1_single_model$boot.coefs[1]+clusboot1_single_model$boot.coefs[2]) + qnorm(0.025) *
                                                                      SE_b1b2)))

sens_higher_CI_boot_single <- (exp((clusboot1_single_model$boot.coefs[1]+clusboot1_single_model$boot.coefs[2]) + qnorm(0.975) *
                                          SE_b1b2)) / (1+(exp((clusboot1_single_model$boot.coefs[1]+clusboot1_single_model$boot.coefs[2]) + qnorm(0.975) *
                                                                      SE_b1b2)))


# sens_boot
# sens_boot_single
# sens_lower_CI_boot
# sens_lower_CI_boot_single
# sens_higher_CI_boot
# sens_higher_CI_boot_single

# calculating specificity from single model: 1-exp(b1) / (1+exp(b1)), where b1 is the intercept
spec_boot_single <- 1 - exp(clusboot1_single_model$boot.coefs[1]) / (1+exp(clusboot1_single_model$boot.coefs[1]))

spec_higher_CI_boot_single <- 1 - exp(clusboot1_single_model$percentile.interval[1])  / (1 + exp(clusboot1_single_model$percentile.interval[1]))
spec_lower_CI_boot_single <- 1 - exp(clusboot1_single_model$percentile.interval[3])  / (1 + exp(clusboot1_single_model$percentile.interval[3]))

# spec_boot
# spec_boot_single
# spec_lower_CI_boot
# spec_lower_CI_boot_single
# spec_higher_CI_boot
# spec_higher_CI_boot_single

#####################################
### Prevalence-dependent measures ###
#####################################

# predictive values can be calculated given sens, spec + prevalence of disease
# calculations based on Ying et al., 2020 and Mercaldo et al., 2009
prev <- sum_D / sum_j # not always appropriate to estimate from study (sampling bias, deliberate oversampling), so usually estimated from separate study


# positive predictive value
PPV <- (sens * prev) / ((sens * prev) + (1-spec) * (1 - prev)) 

# negative predictive value
NPV <- (spec * (1- prev)) / ((1-sens) * prev + spec * (1 - prev)) 

#########################
## Binomial proportion ##
#########################
# NB: not applicable when prevalence cannot be calculated from the study

# PPV
bt_PPV <- binom_test(sum_TP, (sum_TP + sum_FP))
PPV_exact <- bt_PPV$estimate

# NPV 
bt_NPV <- binom_test(sum_TN, (sum_TN + sum_FN))
NPV_exact <- bt_NPV$estimate

##################################
## No adjustment for clustering ##
##################################

#####################################
# Binomial proportion, Delta method #
#####################################
# Mercaldo et al., 2007


# PPV 
PPV_SE_unadj <- sqrt((((prev * (1 - spec)) * (1 - prev))^2 * ((sens * (1 - sens))/sum_D) + 
                              (prev * sens * (1 - prev))^2 * ((spec * (1 - spec)) / sum_N)) / 
                             ((sens * prev + (1 - spec) * (1 - prev))^4))

PPV_CI_unadj <- 1.96*PPV_SE_unadj
PPV_lower_CI_unadj <- PPV-PPV_CI_unadj
PPV_higher_CI_unadj <- PPV+PPV_CI_unadj   

# PPV
# PPV_higher_CI_unadj
# PPV_lower_CI_unadj


# NPV
NPV_SE_unadj <- sqrt((((spec * (1 - prev)) * prev)^2 * ((sens * (1 - sens)) / sum_D) +
                              ((1 - sens) * (1 - prev) * prev)^2 * ((spec * (1  - spec)) / sum_N)) /
                             (((1 - sens) * prev + spec * (1 - prev))^4))

NPV_CI_unadj <- 1.96*NPV_SE_unadj
NPV_lower_CI_unadj <- NPV-NPV_CI_unadj
NPV_higher_CI_unadj <- NPV+NPV_CI_unadj     

# NPV
# NPV_higher_CI_unadj
# NPV_lower_CI_unadj

##########################################################
# Binomial proportion, exact (Clopper-Pearson) intervals #
##########################################################
# NB: exact variances for predictive values are misleading if prevalence cannot be determined from the data

# PPV
PPV_lower_CI_exact <- bt_PPV$conf.low
PPV_higher_CI_exact <- bt_PPV$conf.high

# PPV_exact
# PPV_lower_CI_exact
# PPV_higher_CI_exact


# NPV
NPV_lower_CI_exact <- bt_NPV$conf.low
NPV_higher_CI_exact <- bt_NPV$conf.high

# NPV_exact
# NPV_lower_CI_exact
# NPV_higher_CI_exact

# ratio estimator and VIF adjustment are not applicable to PPV/NPV CI's because PV's are contingent on the prevalence as well. 
# Mercaldo looks at standard (delta method) predictive values and their CI's vs. continuity corrected estimates and CI's vs. 
# their logistic transformations. Other CI's exist for binomial proportions but Wilson's rely on the assumption of not having
# more than two roots and Jeffrey's on the beta distribution. 

#########################
## Logistic regression ##
#########################
# NB: not applicable when prevalence cannot be calculated from the study 

# long format is required 
# PPV
#glm, -1 suppresses intercept
glm3 <- glm(Disease ~ Test - 1, family = binomial(link="logit"), data=data.long)
#summary(glm3)

# probability scale (exp/(1+exp))
PPV_exp <- exp(summary(glm3)$coefficients[,1]) / (1+exp(summary(glm3)$coefficients[,1]))


PPV_lower_CI_exp <- exp(summary(glm3)$coefficients[,1] + qnorm(0.025) *
                        summary(glm3)$coefficients[,2]) / (1+exp(summary(glm3)$coefficients[,1] + qnorm(0.025) *
                                                                                summary(glm3)$coefficients[,2]))


PPV_higher_CI_exp <- (exp(summary(glm3)$coefficients[,1] + qnorm(0.975) *
                        summary(glm3)$coefficients[,2]) / (1+exp(summary(glm3)$coefficients[,1] + qnorm(0.975) *
                                                                                summary(glm3)$coefficients[,2])))

# PPV_exp
# PPV_lower_CI_exp
# PPV_higher_CI_exp

# NPV
#glm, -1 suppresses intercept
glm4 <- glm(Inv_disease ~ Inv_test - 1, family = binomial(link="logit"), data=data.long)
#summary(glm4)

# probability scale (exp/(1+exp))
NPV_exp <- exp(summary(glm4)$coefficients[,1]) / (1+exp(summary(glm4)$coefficients[,1]))


NPV_lower_CI_exp <- exp(summary(glm4)$coefficients[,1] + qnorm(0.025) *
                                summary(glm4)$coefficients[,2]) / (1+exp(summary(glm4)$coefficients[,1] + qnorm(0.025) *
                                                                                 summary(glm4)$coefficients[,2]))


NPV_higher_CI_exp <- (exp(summary(glm4)$coefficients[,1] + qnorm(0.975) *
                                  summary(glm4)$coefficients[,2]) / (1+exp(summary(glm4)$coefficients[,1] + qnorm(0.975) *
                                                                                   summary(glm4)$coefficients[,2])))

# NPV_exp
# NPV_lower_CI_exp
# NPV_higher_CI_exp

#####################################
# Logistic regression, single model #
#####################################

# PPV
glm3_single_model <- glm(Disease ~ Test, family = binomial(link="logit"), data=data.long)
#summary(glm3_single_model)

# calculating PPV from single model: exp(b1+b2) / (1+exp(b1+b2)), where b1 is the intercept and b2 is the coeff
PPV_glm_single <- exp(glm3_single_model$coefficients[1] + glm3_single_model$coefficients[2]) / 
        (1+ exp(glm3_single_model$coefficients[1] + glm3_single_model$coefficients[2]))



# SE(b1+b2) = sqrt(Var(b1)+Var(b2)+2Cov(b1, b2))

#Var(b1)
#vcov(summary(glm3_single_model))[1]

#Var(b2)
#vcov(summary(glm3_single_model))[4]

#Cov(b1,b2)
#vcov(summary(glm3_single_model))[2]

SE_b1b2 <- sqrt(vcov(summary(glm3_single_model))[1]+vcov(summary(glm3_single_model))[4]+2*vcov(summary(glm3_single_model))[2])
SE_b1b2

PPV_lower_CI_glm_single <- (exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.025) *
                                           SE_b1b2)) / (1+(exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.025) *
                                                                       SE_b1b2)))

PPV_higher_CI_glm_single <- (exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.975) *
                                            SE_b1b2)) / (1+(exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.975) *
                                                                        SE_b1b2)))

# PPV_exp
# PPV_glm_single
# PPV_lower_CI_exp
# PPV_lower_CI_glm_single
# PPV_higher_CI_exp
# PPV_higher_CI_glm_single

# calculating NPV from single model: 1-exp(b1) / (1+exp(b1)), where b1 is the intercept
NPV_glm_single <- 1 - exp(glm3_single_model$coefficients[1]) / (1+exp(glm3_single_model$coefficients[1]))

# 95%CI 
NPV_higher_CI_glm_single <- 1 - exp((glm3_single_model$coefficients[1]) + qnorm(0.025) *
                                               sqrt(vcov(summary(glm3_single_model))[1])) / (1 + exp((glm3_single_model$coefficients[1]) + qnorm(0.025) *
                                                                                                             sqrt(vcov(summary(glm3_single_model))[1])))

NPV_lower_CI_glm_single <- 1 - exp((glm3_single_model$coefficients[1]) + qnorm(0.975) *
                                              sqrt(vcov(summary(glm3_single_model))[1])) / (1 + exp((glm3_single_model$coefficients[1]) + qnorm(0.975) *
                                                                                                            sqrt(vcov(summary(glm3_single_model))[1])))
# NPV_exp
# NPV_glm_single
# NPV_lower_CI_exp
# NPV_lower_CI_glm_single
# NPV_higher_CI_exp
# NPV_higher_CI_glm_single


########################
# Logit transformation #
########################

# When predictive values are close to 0 or 1, their CI's can be out of range. Logit transformation ensures that the CI's always 
# fall between the possible values

# Calculations on the logistic scale
# PPV
PPV_logit <- log(PPV/(1-PPV))
PPV_SE_logit <- sqrt(((1-sens)/sens)*(1/sum_D)+(spec/(1-spec))*(1/sum_N))

PPV_CI_logit <- 1.96*PPV_SE_logit
PPV_lower_CI_logit <- PPV_logit-PPV_CI_logit
PPV_higher_CI_logit <- PPV_logit+PPV_CI_logit      

# exponentiating the CI values to backtransform
PPV_lower_CI_logit_exp <- exp(PPV_lower_CI_logit)/((1+exp(PPV_lower_CI_logit)))
PPV_higher_CI_logit_exp <- exp(PPV_higher_CI_logit)/((1+exp(PPV_higher_CI_logit)))

# PPV
# PPV_lower_CI_logit_exp
# PPV_higher_CI_logit_exp

# NPV
NPV_logit <- log(NPV/(1-NPV))
NPV_SE_logit <- sqrt((sens/(1-sens))*(1/sum_D)+((1-spec)/spec)*(1/sum_N))
NPV_CI_logit <- 1.96*NPV_SE_logit
NPV_lower_CI_logit <- NPV_logit-NPV_CI_logit
NPV_higher_CI_logit <- NPV_logit+NPV_CI_logit     

# exponentiating the CI values to backtransform
NPV_lower_CI_logit_exp <- exp(NPV_lower_CI_logit)/((1+exp(NPV_lower_CI_logit)))
NPV_higher_CI_logit_exp <- exp(NPV_higher_CI_logit)/((1+exp(NPV_higher_CI_logit)))

# NPV
# NPV_lower_CI_logit_exp
# NPV_higher_CI_logit_exp

#########################
# Continuity adjustment #
#########################
# suggested when sample sizes are small and PV's are close to 0 or 1

# first calculate the continuity-corrected versions of sens and spec
sens_cont <- (sum_D * sens + k^2/2)/(sum_D+k^2)
spec_cont <- (sum_N * spec + k^2/2)/(sum_N+k^2)

# then calculate the continuity-corrected versions of PV's:
# PPV
PPV_cont <- (sens_cont * prev) / ((sens_cont * prev) + (1-spec_cont) * (1 - prev)) # positive predictive value

# calculating continuity-corrected CI's with delta method 
PPV_SE_cont <- sqrt((((prev * (1 - spec_cont)) * (1 - prev))^2 * ((sens_cont * (1 - sens_cont))/sum_D) + 
                              (prev * sens_cont * (1 - prev))^2 * ((spec_cont * (1 - spec_cont)) / sum_N)) / 
                             ((sens_cont * prev + (1 - spec_cont) * (1 - prev))^4))

PPV_CI_cont <- 1.96 * PPV_SE_cont
PPV_lower_CI_cont <- PPV_cont - PPV_CI_cont
PPV_higher_CI_cont <- PPV_cont + PPV_CI_cont       

# PPV_cont
# PPV_lower_CI_cont
# PPV_higher_CI_cont

# NPV
NPV_cont <- (spec_cont * (1- prev)) / ((1-sens_cont) * prev + spec_cont * (1 - prev)) # negative predictive value

NPV_SE_cont <- sqrt((((spec_cont * (1 - prev)) * prev)^2 * ((sens_cont * (1 - sens_cont)) / sum_D) +
                              ((1 - sens_cont) * (1 - prev) * prev)^2 * ((spec_cont * (1  - spec_cont)) / sum_N)) /
                             (((1 - sens_cont) * prev + spec_cont * (1 - prev))^4))

NPV_CI_cont <- 1.96 * NPV_SE_cont
NPV_lower_CI_cont <- NPV_cont - NPV_CI_cont
NPV_higher_CI_cont <- NPV_cont + NPV_CI_cont     

# NPV_cont
# NPV_lower_CI_cont
# NPV_higher_CI_cont


##################################################
# Logit transformation and continuity correction #
##################################################

# PPV
PPV_cont_logit <- log(PPV_cont/(1-PPV_cont))
PPV_SE_cont_logit <- sqrt(((1-sens_cont)/sens_cont)*(1/sum_D)+(spec_cont/(1-spec_cont))*(1/sum_N))

PPV_CI_cont_logit <- 1.96 * PPV_SE_cont_logit
PPV_lower_CI_cont_logit <- PPV_cont_logit - PPV_CI_cont_logit
PPV_higher_CI_cont_logit <- PPV_cont_logit + PPV_CI_cont_logit      

# exponentiating the CI values to backtransform
PPV_lower_CI_cont_logit_exp <- exp(PPV_lower_CI_cont_logit)/((1+exp(PPV_lower_CI_cont_logit)))
PPV_higher_CI_cont_logit_exp <- exp(PPV_higher_CI_cont_logit)/((1+exp(PPV_higher_CI_cont_logit)))

# PPV_cont
# PPV_lower_CI_cont_logit_exp
# PPV_higher_CI_cont_logit_exp


# NPV
NPV_cont_logit <- log(NPV_cont/(1-NPV_cont))
NPV_SE_cont_logit <- sqrt((sens_cont/(1-sens_cont))*(1/sum_D)+((1-spec_cont)/spec_cont)*(1/sum_N))
NPV_CI_cont_logit <- 1.96*NPV_SE_cont_logit
NPV_lower_CI_cont_logit <- NPV_cont_logit - NPV_CI_cont_logit
NPV_higher_CI_cont_logit <- NPV_cont_logit + NPV_CI_cont_logit     

# exponentiating the CI values to backtransform
NPV_lower_CI_cont_logit_exp <- exp(NPV_lower_CI_cont_logit)/((1+exp(NPV_lower_CI_cont_logit)))
NPV_higher_CI_cont_logit_exp <- exp(NPV_higher_CI_cont_logit)/((1+exp(NPV_higher_CI_cont_logit)))

# NPV_cont
# NPV_lower_CI_cont_logit_exp
# NPV_higher_CI_cont_logit_exp


######################################
# Logistic regression with robust SE #
######################################
# sandwich estimator, c.f., Williams (2000), implementation aided by McDonald (2019)

# PPV
# general sandwich estimator: no regard for the type of cluster (in stata: vce (robust))
# when we have no information about how the data is clustered together, but we know they are clustered, it might be an option
PPV_lower_CI_robust_logit <- coefci(glm3, vcov. = vcovHC(glm3, type = 'HC1'))[1] 
PPV_higher_CI_robust_logit <- coefci(glm3, vcov. = vcovHC(glm3, type = 'HC1'))[2]

# exponentiating to backtransform
PPV_lower_CI_robust <- exp(PPV_lower_CI_robust_logit) / (1+exp(PPV_lower_CI_robust_logit)) # does not completely match stata output 
PPV_higher_CI_robust <- exp(PPV_higher_CI_robust_logit) / (1+exp(PPV_higher_CI_robust_logit))

# we need a sandwich estimator that clusters on id (in stata: vce (cluster id))
PPV_lower_CI_robust_cl_logit <- coefci(glm3, vcov. = vcovCL(glm3, cluster = ~id, type = 'HC1'))[1] 
PPV_higher_CI_robust_cl_logit <- coefci(glm3,  vcov. = vcovCL(glm3, cluster = ~id, type = 'HC1'))[2]

# exponentiating to backtransform
PPV_lower_CI_robust_cl <- exp(PPV_lower_CI_robust_cl_logit) / (1+exp(PPV_lower_CI_robust_cl_logit)) 
PPV_higher_CI_robust_cl <- exp(PPV_higher_CI_robust_cl_logit) / (1+exp(PPV_higher_CI_robust_cl_logit))

# PPV_lower_CI_robust
# PPV_higher_CI_robust
# PPV_lower_CI_robust_cl
# PPV_higher_CI_robust_cl

# NPV
# general sandwich estimator: no regard for the type of cluster (in stata: vce (robust))
# when we have no information about how the data is clustered together, but we know they are clustered, it might be an option
NPV_lower_CI_robust_logit <- coefci(glm4, vcov. = vcovHC(glm4, type = 'HC1'))[1] 
NPV_higher_CI_robust_logit <- coefci(glm4, vcov. = vcovHC(glm4, type = 'HC1'))[2]

# exponentiating to backtransform
NPV_lower_CI_robust <- exp(NPV_lower_CI_robust_logit) / (1+exp(NPV_lower_CI_robust_logit)) # does not completely match stata output 
NPV_higher_CI_robust <- exp(NPV_higher_CI_robust_logit) / (1+exp(NPV_higher_CI_robust_logit))

# we need a sandwich estimator that clusters on id (in stata: vce (cluster id))
NPV_lower_CI_robust_cl_logit <- coefci(glm4, vcov. = vcovCL(glm4, cluster = ~id, type = 'HC1'))[1] 
NPV_higher_CI_robust_cl_logit <- coefci(glm4,  vcov. = vcovCL(glm4, cluster = ~id, type = 'HC1'))[2]

# exponentiating to backtransform
NPV_lower_CI_robust_cl <- exp(NPV_lower_CI_robust_cl_logit) / (1+exp(NPV_lower_CI_robust_cl_logit)) 
NPV_higher_CI_robust_cl <- exp(NPV_higher_CI_robust_cl_logit) / (1+exp(NPV_higher_CI_robust_cl_logit))

# NPV_lower_CI_robust
# NPV_higher_CI_robust
# NPV_lower_CI_robust_cl
# NPV_higher_CI_robust_cl

####################################################
# Logistic regression with robust SE, single model #
####################################################
# sandwich estimator, c.f., Williams (2000), implementation aided by McDonald (2019)

# PPV
# general sandwich estimator: no regard for the type of cluster (in stata: vce (robust))
# when we have no information about how the data is clustered together, but we know they are clustered, it might be an option

# SE(b1+b2) = sqrt(Var(b1)+Var(b2)+2Cov(b1+b2))
# SE_b1b2

# adjusted var(b1)
vcovHC(glm3_single_model, type = 'HC1')[1]

# adjusted var(b2)
vcovHC(glm3_single_model, type = 'HC1')[4]

# adjusted Cov(b1 ,b2)
vcovHC(glm3_single_model, type = 'HC1')[2]

# SE_b1b2
SE_b1b2 <- sqrt(vcovHC(glm3_single_model, type = 'HC1')[1]+vcovHC(glm3_single_model, type = 'HC1')[4]+2*vcovHC(glm3_single_model, type = 'HC1')[2])
SE_b1b2

PPV_lower_CI_robust_single <- (exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.025) *
                                            SE_b1b2)) / (1+(exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.025) *
                                                                        SE_b1b2)))

PPV_higher_CI_robust_single <- (exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.975) *
                                             SE_b1b2)) / (1+(exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.975) *
                                                                         SE_b1b2)))

# PPV_lower_CI_robust
# PPV_lower_CI_robust_single
# PPV_higher_CI_robust
# PPV_higher_CI_robust_single


# NPV 95%CI 
NPV_higher_CI_robust_single <- 1 - exp((glm3_single_model$coefficients[1]) + qnorm(0.025) *
                                                sqrt(vcovHC(glm3_single_model, type = 'HC1')[1])) / (1 + exp((glm3_single_model$coefficients[1]) + qnorm(0.025) *
                                                                                                                     sqrt(vcovHC(glm3_single_model, type = 'HC1')[1])))

NPV_lower_CI_robust_single <- 1 - exp((glm3_single_model$coefficients[1]) + qnorm(0.975) *
                                               sqrt(vcovHC(glm3_single_model, type = 'HC1')[1])) / (1 + exp((glm3_single_model$coefficients[1]) + qnorm(0.975) *
                                                                                                                    sqrt(vcovHC(glm3_single_model, type = 'HC1')[1])))

# NPV_lower_CI_robust
# NPV_lower_CI_robust_single
# NPV_higher_CI_robust
# NPV_higher_CI_robust_single


# we need a sandwich estimator that clusters on id (in stata: vce (cluster id))
# SE(b1+b2) = sqrt(Var(b1)+Var(b2)+2Cov(b1+b2))
# SE_b1b2

# adjusted var(b1)
vcovHC(glm3_single_model, cluster= ~id, type = 'HC1')[1]

# adjusted var(b2)
vcovHC(glm3_single_model, cluster= ~id, type = 'HC1')[4]

# adjusted Cov(b1 ,b2)
vcovHC(glm3_single_model, cluster= ~id, type = 'HC1')[2]

# SE_b1b2
SE_b1b2 <- sqrt(vcovHC(glm3_single_model, cluster= ~id, type = 'HC1')[1]+vcovHC(glm3_single_model, cluster= ~id, type = 'HC1')[4]+2*vcovHC(glm3_single_model, cluster= ~id, type = 'HC1')[2])
SE_b1b2

PPV_lower_CI_robust_cl_single <- (exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.025) *
                                               SE_b1b2)) / (1+(exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.025) *
                                                                           SE_b1b2)))

PPV_higher_CI_robust_cl_single <- (exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.975) *
                                                SE_b1b2)) / (1+(exp((glm3_single_model$coefficients[1]+glm3_single_model$coefficients[2]) + qnorm(0.975) *
                                                                            SE_b1b2)))

# PPV_lower_CI_robust_cl
# PPV_lower_CI_robust_cl_single
# PPV_higher_CI_robust_cl
# PPV_higher_CI_robust_cl_single


# NPV 95%CI robust cluster
NPV_higher_CI_robust_cl_single <- 1 - exp((glm3_single_model$coefficients[1]) + qnorm(0.025) *
                                                   sqrt(vcovHC(glm3_single_model, cluster= ~id, type = 'HC1')[1])) / (1 + exp((glm3_single_model$coefficients[1]) + qnorm(0.025) *
                                                                                                                                      sqrt(vcovHC(glm3_single_model, cluster= ~id, type = 'HC1')[1])))

NPV_lower_CI_robust_cl_single <- 1 - exp((glm3_single_model$coefficients[1]) + qnorm(0.975) *
                                                  sqrt(vcovHC(glm3_single_model, cluster= ~id, type = 'HC1')[1])) / (1 + exp((glm3_single_model$coefficients[1]) + qnorm(0.975) *
                                                                                                                                     sqrt(vcovHC(glm3_single_model, cluster= ~id, type = 'HC1')[1])))

# NPV_lower_CI_robust_cl
# NPV_lower_CI_robust_cl_single
# NPV_higher_CI_robust_cl
# NPV_higher_CI_robust_cl_single

#####################################
# Mixed effects logistic regression #
#####################################
# id in random effect structure as random intercept

# PPV
# NB: syntax for no intercept is different in glmer compared to glm
# nAGQ means 'number of adaptive Gauss-Hermite quadrature points', and sets how glmer will integrate out the random 
# effects when fitting the mixed model. When nAGQ is greater than 1, then adaptive quadrature is used with nAGQ points. When nAGQ =
#1, the Laplace approximation is used, and when nAGQ = 0, the integral is 'ignored'.
glmer3 <-  (glmer(Disease ~ 1 + (1 | id), family=binomial(link = "logit"), nAGQ = 1, data=data.long, subset=Test==1))  
#summary(glmer3)


# converting estimate and CI's to probability scale (exp/(1+exp))
PPV_glmer <-  ((exp(summary(glmer3)$coefficients[,1] + qnorm(0.5) *
                           summary(glmer3)$coefficients[,2]) / (1+exp(summary(glmer3)$coefficients[,1] + qnorm(0.5) *
                                                                              summary(glmer3)$coefficients[,2]))))

PPV_lower_CI_glmer <-  ((exp(summary(glmer3)$coefficients[,1] + qnorm(0.025) *
                                    summary(glmer3)$coefficients[,2]) / (1+exp(summary(glmer3)$coefficients[,1] + qnorm(0.025) *
                                                                                       summary(glmer3)$coefficients[,2]))))

PPV_higher_CI_glmer <-  ((exp(summary(glmer3)$coefficients[,1] + qnorm(0.975) *
                                     summary(glmer3)$coefficients[,2]) / (1+exp(summary(glmer3)$coefficients[,1] + qnorm(0.975) *
                                                                                        summary(glmer3)$coefficients[,2]))))

# PPV_glmer
# PPV_lower_CI_glmer
# PPV_higher_CI_glmer

# 'glmer' uses the Laplacian approximation as default, corresponding to adaptive Gauss-Hermite approximation 
# with only 1 point, while 'xtmelogit' uses 7 points.

# NPV
glmer4 <-  (glmer(Inv_disease ~ 1 + (1 | id), family=binomial(link = "logit"), nAGQ = 1, data=data.long, subset=Inv_test==1))  
#summary(glmer4)

# converting estimate and CI's to probability scale (exp/(1+exp))
NPV_glmer <-  ((exp(summary(glmer4)$coefficients[,1] + qnorm(0.5) *
                           summary(glmer4)$coefficients[,2]) / (1+exp(summary(glmer4)$coefficients[,1] + qnorm(0.5) *
                                                                              summary(glmer4)$coefficients[,2]))))

NPV_lower_CI_glmer <-  ((exp(summary(glmer4)$coefficients[,1] + qnorm(0.025) *
                                    summary(glmer4)$coefficients[,2]) / (1+exp(summary(glmer4)$coefficients[,1] + qnorm(0.025) *
                                                                                       summary(glmer4)$coefficients[,2]))))

NPV_higher_CI_glmer <-  ((exp(summary(glmer4)$coefficients[,1] + qnorm(0.975) *
                                     summary(glmer4)$coefficients[,2]) / (1+exp(summary(glmer4)$coefficients[,1] + qnorm(0.975) *
                                                                                        summary(glmer4)$coefficients[,2]))))

# NPV_glmer
# NPV_lower_CI_glmer
# NPV_higher_CI_glmer

###################################################
# Mixed effects logistic regression, single model #
###################################################

# glmer single model
glmer3_single_model <-  (glmer(Disease ~ Test + (1 | id), nAGQ = 1, family=binomial(link = "logit"), data=data.long))  
#summary(glmer3_single_model)

# calculating PPV from single model: exp(b1+b2) / (1+exp(b1+b2)), where b1 is the intercept and b2 is the coeff
PPV_glmer_single <-  (exp(summary(glmer3_single_model)$coef[1] + summary(glmer3_single_model)$coef[2]) / 
        (1+ exp(summary(glmer3_single_model)$coef[1] + summary(glmer3_single_model)$coef[2])))



# SE(b1+b2) = sqrt(Var(b1)+Var(b2)+2Cov(b1, b2))

#Var(b1)
#vcov(summary(glmer3_single_model))[1]

#Var(b2)
#vcov(summary(glmer3_single_model))[4]

#Cov(b1,b2)
#vcov(summary(glmer3_single_model))[2]

SE_b1b2 <-  (sqrt(vcov(summary(glmer3_single_model))[1]+vcov(summary(glmer3_single_model))[4]+2*vcov(summary(glmer3_single_model))[2]))
#SE_b1b2

PPV_lower_CI_glmer_single <-  ((exp((summary(glmer3_single_model)$coef[1]+summary(glmer3_single_model)$coef[2]) + qnorm(0.025) *
                                        SE_b1b2)) / (1+(exp((summary(glmer3_single_model)$coef[1]+summary(glmer3_single_model)$coef[2]) + qnorm(0.025) *
                                                                    SE_b1b2))))

PPV_higher_CI_glmer_single <-  ((exp((summary(glmer3_single_model)$coef[1]+summary(glmer3_single_model)$coef[2]) + qnorm(0.975) *
                                         SE_b1b2)) / (1+(exp((summary(glmer3_single_model)$coef[1]+summary(glmer3_single_model)$coef[2]) + qnorm(0.975) *
                                                                     SE_b1b2))))

# PPV_glmer
# PPV_glmer_single
# PPV_lower_CI_glmer
# PPV_lower_CI_glmer_single
# PPV_higher_CI_glmer
# PPV_higher_CI_glmer_single

# calculating NPV from single model: 1-exp(b1) / (1+exp(b1)), where b1 is the intercept
NPV_glmer_single <-  (1 - exp(summary(glmer3_single_model)$coef[1]) / (1+exp(summary(glmer3_single_model)$coef[1])))

# 95%CI 
NPV_higher_CI_glmer_single <-  (1 - exp((summary(glmer3_single_model)$coef[1]) + qnorm(0.025) *
                                            sqrt(vcov(summary(glmer3_single_model))[1])) / (1 + exp((summary(glmer3_single_model)$coef[1]) + qnorm(0.025) *
                                                                                                          sqrt(vcov(summary(glmer3_single_model))[1]))))

NPV_lower_CI_glmer_single <-  (1 - exp((summary(glmer3_single_model)$coef[1]) + qnorm(0.975) *
                                           sqrt(vcov(summary(glmer3_single_model))[1])) / (1 + exp((summary(glmer3_single_model)$coef[1]) + qnorm(0.975) *
                                                                                                         sqrt(vcov(summary(glmer3_single_model))[1]))))
# NPV_glmer
# NPV_glmer_single
# NPV_lower_CI_glmer
# NPV_lower_CI_glmer_single
# NPV_higher_CI_glmer
# NPV_higher_CI_glmer_single



##########################################
# Generalized estimating equations (GEE) #
##########################################

# disease (+ vs. -) is modeled as the outcome variable on a subset of id's with positive tests, with a logit link

# NB: data needs to be sorted by id first (if not already)
gee3 <-  (gee(Disease ~ 1, data = data.long, subset=Test==1, id=id, family=binomial(link = "logit"), corstr="exchangeable"))
#gee3 <- gee(Disease ~ 1, data = data.long, subset=Test==1, id=id, family=binomial(link = "logit"), corstr="unstructured")
#gee3 <- gee(Disease ~ 1, data = data.long, subset=Test==1, id=id, family=binomial(link = "logit"), corstr="independence")
#summary(gee3) 
# too sparse data, Warning message:
#In gee(Disease ~ 1, data = data.long, subset = Test == 1, id = id,  :
#               Working correlation estimate not positive definite
       

# converting estimate and CI's to probability scale (exp/(1+exp)) using robust SE's (column 4)
PPV_gee <-  (exp(summary(gee3)$coefficients[,1]) / (1+exp(summary(gee3)$coefficients[,1])))

PPV_lower_CI_gee <-  (exp(summary(gee3)$coefficients[,1] + qnorm(0.025) *
                                 summary(gee3)$coefficients[,4]) / (1+exp(summary(gee3)$coefficients[,1] + qnorm(0.025) *
                                                                                  summary(gee3)$coefficients[,4])))

PPV_higher_CI_gee <-  (exp(summary(gee3)$coefficients[,1] + qnorm(0.975) *
                                  summary(gee3)$coefficients[,4]) / (1+exp(summary(gee3)$coefficients[,1] + qnorm(0.975) *
                                                                                   summary(gee3)$coefficients[,4])))

# PPV_gee
# PPV_lower_CI_gee
# PPV_higher_CI_gee

# NPV
# inv_disease (+ vs. -) is modeled as the outcome variable on a subset of id's with negative tests, with a logit link

gee4 <-  (gee(Inv_disease ~ 1, data = data.long, subset=Inv_test==1, id=id, family=binomial(link = "logit"), corstr="exchangeable"))
#gee4 <- gee(Inv_disease ~ 1, data = data.long, subset=Inv_test==1, id=id, family=binomial(link = "logit"), corstr="unstructured")
#summary(gee4)

# converting estimate and CI's to probability scale (exp/(1+exp)) using robust SE's (column 4)
NPV_gee <-  (exp(summary(gee4)$coefficients[,1]) / (1+exp(summary(gee4)$coefficients[,1])))

NPV_lower_CI_gee <-  (exp(summary(gee4)$coefficients[,1] + qnorm(0.025) *
                                 summary(gee4)$coefficients[,4]) / (1+exp(summary(gee4)$coefficients[,1] + qnorm(0.025) *
                                                                                  summary(gee4)$coefficients[,4])))

NPV_higher_CI_gee <-  (exp(summary(gee4)$coefficients[,1] + qnorm(0.975) *
                                  summary(gee4)$coefficients[,4]) / (1+exp(summary(gee4)$coefficients[,1] + qnorm(0.975) *
                                                                                   summary(gee4)$coefficients[,4])))

# NPV_gee
# NPV_lower_CI_gee
# NPV_higher_CI_gee

########################################################
# Generalized estimating equations (GEE), single model #
########################################################

# disease (+ vs. -) is modeled as the outcome variable on a subset of id's with positive tests, with a logit link
# NB: data needs to be sorted by id first (if not already)
# PPV 
gee3_single_model <-  (gee(Disease ~ Test, data = data.long, id=id, family=binomial(link = "logit"), corstr="exchangeable"))
#gee3_single_model <- gee(Disease ~ Test, data = data.long, id=id, family=binomial(link = "logit"), corstr="unstructured")
#summary(gee3_single_model)

#coef(summary(gee3_single_model))

# calculating PPV from single model: exp(b1+b2) / (1+exp(b1+b2)), where b1 is the intercept and b2 is the coeff
PPV_gee_single <-  ((exp(summary(gee3_single_model)$coefficients[[1]] + summary(gee3_single_model)$coefficients[[2]])) / 
        (1+(exp(summary(gee3_single_model)$coefficients[[1]] + summary(gee3_single_model)$coefficients[[2]]))))


# handcoding linear combinations to obtain the 95% CI of the combinations of the intercept and coef
# SE_b1b2 <- sqrt(vcov(glmer1_single_model)[1]+vcov(glmer1_single_model)[4]+2*(vcov(glmer1_single_model)[2]))
# SE_b1b2

#Var(b1)
#gee3_single_model$robust.variance[1]

#Var(b2)
#gee3_single_model$robust.variance[4]

#Cov(b1,b2)
#gee3_single_model$robust.variance[2]

SE_b1b2 <-  (sqrt(gee3_single_model$robust.variance[1]+gee3_single_model$robust.variance[4]+2*gee3_single_model$robust.variance[2]))
#SE_b1b2

PPV_lower_CI_gee_single <-  ((exp((summary(gee3_single_model)$coefficients[[1]] + summary(gee3_single_model)$coefficients[[2]]) + qnorm(0.025) *
                                         SE_b1b2)) / (1+(exp((summary(gee3_single_model)$coefficients[[1]] + summary(gee3_single_model)$coefficients[[2]]) + qnorm(0.025) *
                                                                     SE_b1b2))))

PPV_higher_CI_gee_single <-  ((exp((summary(gee3_single_model)$coefficients[[1]] + summary(gee3_single_model)$coefficients[[2]]) + qnorm(0.975) *
                                          SE_b1b2)) / (1+(exp((summary(gee3_single_model)$coefficients[[1]] + summary(gee3_single_model)$coefficients[[2]]) + qnorm(0.975) *
                                                                      SE_b1b2))))


# calculating NPV from single model: 1-exp(b1) / (1+exp(b1)), where b1 is the intercept
NPV_gee_single <-  (1-exp(summary(gee3_single_model)$coefficients[[1]]) / (1+exp(summary(gee3_single_model)$coefficients[[1]])))


NPV_higher_CI_gee_single <-  (1- (exp(summary(gee3_single_model)$coefficients[[1]] + qnorm(0.025) *
                                             summary(gee3_single_model)$coefficients[[7]])) / (1+exp(summary(gee3_single_model)$coefficients[[1]] + qnorm(0.025) *
                                                                                                             summary(gee3_single_model)$coefficients[[7]])))

NPV_lower_CI_gee_single <-  (1- (exp(summary(gee3_single_model)$coefficients[[1]] + qnorm(0.975) *
                                            summary(gee3_single_model)$coefficients[[7]])) / (1+exp(summary(gee3_single_model)$coefficients[[1]] + qnorm(0.975) *
                                                                                                            summary(gee3_single_model)$coefficients[[7]])))

# PPV_gee
# PPV_gee_single
# PPV_lower_CI_gee
# PPV_lower_CI_gee_single
# PPV_higher_CI_gee
# PPV_higher_CI_gee_single
# 
# NPV_gee
# NPV_gee_single
# NPV_lower_CI_gee
# NPV_lower_CI_gee_single
# NPV_higher_CI_gee
# NPV_higher_CI_gee_single

#####################
# Cluster bootstrap #
#####################

# taking the output of the cluster bootstrap sens and spec estimates sens_boot and spec_boot as input 

# PPV

# calculated an individual PPV and NPV for each of those coefs and took the median, the 2.5th and the 
# 97.5th percentiles of the resulting distribution

# sorting clusboot1$coefficients and clusboot1$coefficients


# converting the sens and spec coeffs to probabilities
clusboot1$sens <- exp(clusboot1$coefficients) / (1+exp(clusboot1$coefficients))
clusboot2$spec <- exp(clusboot2$coefficients) / (1+exp(clusboot2$coefficients))

PPV_boots <- (clusboot1$sens * prev) / ((clusboot1$sens * prev) + (1 - clusboot2$spec) * (1 - prev))

PPV_boot <- median(PPV_boots)
PPV_lower_CI_boot <- quantile(PPV_boots, probs = 0.025, na.rm = T)
PPV_higher_CI_boot <- quantile(PPV_boots, probs = 0.975, na.rm = T)

#hist(PPV_boots)
#abline(v=c(PPV_boot, PPV_lower_CI_boot, PPV_higher_CI_boot), col="red")

# ggplot() + geom_histogram(aes(PPV_boots), color= "white") +
#         xlab("Histogram of PPV estimate based on 10000 samples") +
#         geom_vline(aes(xintercept = c(PPV_boot, PPV_lower_CI_boot, PPV_higher_CI_boot)), colour="red") +
#         theme_bw()

# PPV_boot
# PPV_lower_CI_boot
# PPV_higher_CI_boot

# NPV
NPV_boots <- (clusboot2$spec * (1-prev)) / ((1-clusboot1$sens) * prev + clusboot2$spec * (1 - prev))
#hist(NPV_boots, brakes=5)
#hist(NPV_boots, breaks = seq(min(NPV_boots), max(NPV_boots), length.out = 11))

NPV_boot <- median(NPV_boots)
NPV_lower_CI_boot <- quantile(NPV_boots, probs = 0.025, na.rm = T)
NPV_higher_CI_boot <- quantile(NPV_boots, probs = 0.975, na.rm = T)

# NPV_boot
# NPV_lower_CI_boot
# NPV_higher_CI_boot


###################################
# Cluster bootstrap, single model #
###################################

# PPV
clusboot3_single_model <- clusbootglm(Disease ~ Test, clusterid=id, data=data.long, family=binomial, B = 10000, confint.level = 0.95)

# calculating PPV from single model: exp(b1+b2) / (1+exp(b1+b2)), where b1 is the intercept and b2 is the coeff
PPV_boot_single <- exp(clusboot3_single_model$boot.coefs[1]+clusboot3_single_model$boot.coefs[2])  / (1+exp(clusboot3_single_model$boot.coefs[1]+clusboot3_single_model$boot.coefs[2]))

# SE(b1+b2) = sqrt(Var(b1)+Var(b2)+2Cov(b1, b2))
# SE_b1b2 <- sqrt(vcov(glmer1_single_model)[1]+vcov(glmer1_single_model)[4]+2*(vcov(glmer1_single_model)[2]))
# SE_b1b2

#Var(b1)
#cov(clusboot3_single_model$coefficients)[1]

#Var(b2)
#cov(clusboot3_single_model$coefficients)[4]

#Cov(b1,b2)
#cov(clusboot3_single_model$coefficients)[2]

SE_b1b2 <- sqrt(cov(clusboot3_single_model$coefficients)[1]+cov(clusboot3_single_model$coefficients)[4]+2*cov(clusboot3_single_model$coefficients)[2])
#SE_b1b2

PPV_lower_CI_boot_single <- (exp((clusboot3_single_model$boot.coefs[1]+clusboot3_single_model$boot.coefs[2]) + qnorm(0.025) *
                                          SE_b1b2)) / (1+(exp((clusboot3_single_model$boot.coefs[1]+clusboot3_single_model$boot.coefs[2]) + qnorm(0.025) *
                                                                      SE_b1b2)))

PPV_higher_CI_boot_single <- (exp((clusboot3_single_model$boot.coefs[1]+clusboot3_single_model$boot.coefs[2]) + qnorm(0.975) *
                                           SE_b1b2)) / (1+(exp((clusboot3_single_model$boot.coefs[1]+clusboot3_single_model$boot.coefs[2]) + qnorm(0.975) *
                                                                       SE_b1b2)))


# PPV_boot
# PPV_boot_single
# PPV_lower_CI_boot
# PPV_lower_CI_boot_single
# PPV_higher_CI_boot
# PPV_higher_CI_boot_single

# calculating NPV from single model: 1-exp(b1) / (1+exp(b1)), where b1 is the intercept
NPV_boot_single <- 1 - exp(clusboot3_single_model$boot.coefs[1]) / (1+exp(clusboot3_single_model$boot.coefs[1]))

NPV_higher_CI_boot_single <- 1 - exp(clusboot3_single_model$percentile.interval[1])  / (1 + exp(clusboot3_single_model$percentile.interval[1]))
NPV_lower_CI_boot_single <- 1 - exp(clusboot3_single_model$percentile.interval[3])  / (1 + exp(clusboot3_single_model$percentile.interval[3]))

# NPV_boot
# NPV_boot_single
# NPV_lower_CI_boot
# NPV_lower_CI_boot_single
# NPV_higher_CI_boot
# NPV_higher_CI_boot_single


#################
#### Outputs ####
#################


# outputting results from patient-level analyses, prevalence-indep measures: sens spec, binomial proportion (exact), binomial proportion (Wald), logistic regression
results_id_sensspec <- data.frame(Type_of_method =c('Binomial (Exact)', 'Binomial (Exact)', 
                                                    'Binomial (Wald)', 'Binomial (Wald)', 
                                                    'Logistic regression', 'Logistic regression', 
                                                    'Logistic regression, single m.', 'Logistic regression, single m.'),
                                  Type_of_diagnostic = rep(c('Sensitivity', 'Specificity'),4),
                                  Estimate = (c(id_sens_exact, id_spec_exact, 
                                                id_sens, id_spec, 
                                                id_sens_exp, id_spec_exp, 
                                                id_sens_glm_single, id_spec_glm_single)),
                                  Lower_CI = c(id_sens_lower_CI_exact, id_spec_lower_CI_exact,
                                               id_sens_lower_CI_unadj, id_spec_lower_CI_unadj, 
                                               id_sens_lower_CI_exp, id_spec_lower_CI_exp, 
                                               id_sens_lower_CI_glm_single, id_spec_lower_CI_glm_single),
                                  Higher_CI = c(id_sens_higher_CI_exact, id_spec_higher_CI_exact, 
                                                id_sens_higher_CI_unadj, id_spec_higher_CI_unadj, 
                                                id_sens_higher_CI_exp, id_spec_higher_CI_exp, 
                                                id_sens_higher_CI_glm_single, id_spec_higher_CI_glm_single))

# outputting results from patient-level analyses, prevalence-dep measures: predictive values
results_id_PV <- data.frame(Type_of_method =c('Binomial (Exact)', 'Binomial (Exact)', 
                                              'Binomial (Wald)', 'Binomial (Wald)', 
                                              'Logit transform', 'Logit transform', 
                                              'Cont-corrected', 'Cont-corrected', 
                                              'Logit + Cont', 'Logit + Cont', 
                                              'Logistic regression', 'Logistic regression', 
                                              'Logistic regression, single m.', 'Logistic regression, single m.'),
                            Type_of_diagnostic = rep(c('PPV', 'NPV'),7),
                            Estimate = (c(id_PPV_exact, id_NPV_exact, 
                                          id_PPV, id_NPV, 
                                          id_PPV, id_NPV, 
                                          id_PPV_cont, id_NPV_cont, 
                                          id_PPV_cont, id_NPV_cont, 
                                          id_PPV_exp, id_NPV_exp, 
                                          id_PPV_glm_single, id_NPV_glm_single)),
                            Lower_CI = c(id_PPV_lower_CI_exact, id_NPV_lower_CI_exact, 
                                         id_PPV_lower_CI_unadj, id_NPV_lower_CI_unadj, 
                                         id_PPV_lower_CI_logit_exp, id_NPV_lower_CI_logit_exp, 
                                         id_PPV_lower_CI_cont, id_NPV_lower_CI_cont, 
                                         id_PPV_lower_CI_cont_logit_exp, id_NPV_lower_CI_cont_logit_exp, 
                                         id_PPV_lower_CI_exp, id_NPV_lower_CI_exp, 
                                         id_PPV_lower_CI_glm_single, id_NPV_lower_CI_glm_single),
                            Higher_CI = c(id_PPV_higher_CI_exact, id_NPV_higher_CI_exact, 
                                          id_PPV_higher_CI_unadj, id_NPV_higher_CI_unadj, 
                                          id_PPV_higher_CI_logit_exp, id_NPV_higher_CI_logit_exp, 
                                          id_PPV_higher_CI_cont, id_NPV_higher_CI_cont, 
                                          id_PPV_higher_CI_cont_logit_exp, id_NPV_higher_CI_cont_logit_exp, 
                                          id_PPV_higher_CI_exp, id_NPV_higher_CI_exp, 
                                          id_PPV_higher_CI_glm_single, id_NPV_higher_CI_glm_single))

# outputting results from methods with no adjustment and method with adjustment 
results_sensspec <- data.frame(Type_of_method =c('Binomial (Exact)', 'Binomial (Exact)', 
                                                 'Binomial (Wald)', 'Binomial (Wald)', 
                                                 'Logistic regression', 'Logistic regression', 
                                                 'Logistic regression, single m.', 'Logistic regression, single m.', 
                                                 'Ratio estimator', 'Ratio estimator',
                                                 'VIF', 'VIF', 
                                                 'Logistic regression robust SE (gen)', 'Logistic regression robust SE (gen)', 
                                                 'Logistic regression robust SE (gen), single m.', 'Logistic regression robust SE (gen), single m.',
                                                 'Logistic regression robust SE (id)', 'Logistic regression robust SE (id)', 
                                                 'Logistic regression robust SE (id), single m.', 'Logistic regression robust SE (id), single m.', 
                                                 'Mixed effects logistic regression', 'Mixed effects logistic regression', 
                                                 'Mixed effects logistic regression, single m.', 'Mixed effects logistic regression, single m.',
                                                 "GEE", 'GEE', 
                                                 "GEE, single m.", 'GEE, single m.',
                                                 "Cluster bootstrap", "Cluster bootstrap", 
                                                 "Cluster bootstrap, single m.", "Cluster bootstrap, single m."),    
                               Type_of_diagnostic = rep(c('Sensitivity', 'Specificity'),16),
                               Estimate = (c(sens_exact, spec_exact, 
                                             sens, spec, 
                                             sens_exp, spec_exp, 
                                             sens_glm_single, spec_glm_single, 
                                             sens, spec, 
                                             sens, spec, 
                                             sens_exp, spec_exp, 
                                             sens_exp, spec_exp,
                                             sens_exp, spec_exp,
                                             sens_exp, spec_exp, 
                                             sens_glmer, spec_glmer, 
                                             sens_glmer_single, spec_glmer_single, 
                                             sens_gee, spec_gee, 
                                             sens_gee_single, spec_gee_single, 
                                             sens_boot, spec_boot, 
                                             sens_boot_single, spec_boot_single)),
                               Lower_CI = c(sens_lower_CI_exact, spec_lower_CI_exact, 
                                            sens_lower_CI_unadj, spec_lower_CI_unadj, 
                                            sens_lower_CI_exp, spec_lower_CI_exp, 
                                            sens_lower_CI_glm_single, spec_lower_CI_glm_single, 
                                            sens_lower_CI_readj, spec_lower_CI_readj, 
                                            sens_lower_CI_vifadj, spec_lower_CI_vifadj, 
                                            sens_lower_CI_robust, spec_lower_CI_robust, 
                                            sens_lower_CI_robust_single, spec_lower_CI_robust_single, 
                                            sens_lower_CI_robust_cl, spec_lower_CI_robust_cl,
                                            sens_lower_CI_robust_cl_single, spec_lower_CI_robust_cl_single, 
                                            sens_lower_CI_glmer, spec_lower_CI_glmer, 
                                            sens_lower_CI_glmer_single, spec_lower_CI_glmer_single, 
                                            sens_lower_CI_gee, spec_lower_CI_gee, 
                                            sens_lower_CI_gee_single, spec_lower_CI_gee_single, 
                                            sens_lower_CI_boot, spec_lower_CI_boot, 
                                            sens_lower_CI_boot_single, spec_lower_CI_boot_single),
                               Higher_CI = c(sens_higher_CI_exact, spec_higher_CI_exact, 
                                             sens_higher_CI_unadj, spec_higher_CI_unadj, 
                                             sens_higher_CI_exp, spec_higher_CI_exp, 
                                             sens_higher_CI_glm_single, spec_higher_CI_glm_single, 
                                             sens_higher_CI_readj, spec_higher_CI_readj, 
                                             sens_higher_CI_vifadj, spec_higher_CI_vifadj, 
                                             sens_higher_CI_robust, spec_higher_CI_robust, 
                                             sens_higher_CI_robust_single, spec_higher_CI_robust_single, 
                                             sens_higher_CI_robust_cl, spec_higher_CI_robust_cl,
                                             sens_higher_CI_robust_cl_single, spec_higher_CI_robust_cl_single,
                                             sens_higher_CI_glmer, spec_higher_CI_glmer, 
                                             sens_higher_CI_glmer_single, spec_higher_CI_glmer_single, 
                                             sens_higher_CI_gee, spec_higher_CI_gee, 
                                             sens_higher_CI_gee_single, spec_higher_CI_gee_single, 
                                             sens_higher_CI_boot, spec_higher_CI_boot, 
                                             sens_higher_CI_boot_single, spec_higher_CI_boot_single))


# outputting results from Method 1 (no adjustment) and Method 2 (adjustment with ratio estimator and variance inflation factor) 
results_PV <- data.frame(Type_of_method =c( 'Binomial (Exact)', 'Binomial (Exact)',
                                           'Binomial (Wald)', 'Binomial (Wald)',
                                           'Logit transform', 'Logit transform', 
                                           'Cont-corrected', 'Cont-corrected', 
                                           'Logit + Cont', 'Logit + Cont', 
                                           'Logistic regression', 'Logistic regression', 
                                           'Logistic regression, single m.', 'Logistic regression, single m.', 
                                           'Logistic regression robust SE (gen)', 'Logistic regression robust SE (gen)', 
                                           'Logistic regression robust SE (gen), single m.', 'Logistic regression robust SE (gen), single m.', 
                                           'Logistic regression robust SE (id)', 'Logistic regression robust SE (id)', 
                                           'Logistic regression robust SE (id), single m.', 'Logistic regression robust SE (id), single m.', 
                                           'Mixed effects logistic regression', 'Mixed effects logistic regression',
                                           'Mixed effects logistic regression, single m.', 'Mixed effects logistic regression, single m.', 
                                           "GEE", 'GEE', 
                                           "GEE, single m.", 'GEE, single m.', 
                                           "Cluster bootstrap", "Cluster bootstrap", 
                                           "Cluster bootstrap, single m.", "Cluster bootstrap, single m."),
                         Type_of_diagnostic = rep(c('PPV', 'NPV'),17),
                         Estimate = (c(PPV_exact, NPV_exact, 
                                       PPV, NPV, 
                                       PPV, NPV, 
                                       PPV_cont, NPV_cont, 
                                       PPV_cont, NPV_cont, 
                                       PPV_exp, NPV_exp, 
                                       PPV_exp, NPV_exp, 
                                       PPV_exp, NPV_exp, 
                                       PPV_exp, NPV_exp, 
                                       PPV_exp, NPV_exp, 
                                       PPV_exp, NPV_exp, 
                                       PPV_glmer, NPV_glmer, 
                                       PPV_glmer_single, NPV_glmer_single,
                                       PPV_gee, NPV_gee, 
                                       PPV_gee_single, NPV_gee_single,
                                       PPV_boot, NPV_boot,
                                       PPV_boot_single, NPV_boot_single)),
                         Lower_CI = c(PPV_lower_CI_exact, NPV_lower_CI_exact, # exact
                                      PPV_lower_CI_unadj, NPV_lower_CI_unadj, # wald
                                      PPV_lower_CI_logit_exp, NPV_lower_CI_logit_exp, # logit transform
                                      PPV_lower_CI_cont, NPV_lower_CI_cont, # cont-corrected
                                      PPV_lower_CI_cont_logit_exp, NPV_lower_CI_cont_logit_exp,  # logit+cont
                                      PPV_lower_CI_exp, NPV_lower_CI_exp, # logreg
                                      PPV_lower_CI_glm_single, NPV_lower_CI_glm_single, # logreg single
                                      PPV_lower_CI_robust, NPV_lower_CI_robust, # robust
                                      PPV_lower_CI_robust_single, NPV_lower_CI_robust_single, # robust single
                                      PPV_lower_CI_robust_cl, NPV_lower_CI_robust_cl, # robust id 
                                      PPV_lower_CI_robust_cl_single, NPV_lower_CI_robust_cl_single, # robust id single
                                      PPV_lower_CI_glmer, NPV_lower_CI_glmer, # glmer
                                      PPV_lower_CI_glmer_single, NPV_lower_CI_glmer_single, # glmer single
                                      PPV_lower_CI_gee, NPV_lower_CI_gee, 
                                      PPV_lower_CI_gee_single, NPV_lower_CI_gee_single,
                                      PPV_lower_CI_boot, NPV_lower_CI_boot, 
                                      PPV_lower_CI_boot_single, NPV_lower_CI_boot_single),
                         Higher_CI = c(PPV_higher_CI_exact, NPV_higher_CI_exact, 
                                       PPV_higher_CI_unadj, NPV_higher_CI_unadj, 
                                       PPV_higher_CI_logit_exp, NPV_higher_CI_logit_exp, 
                                       PPV_higher_CI_cont, NPV_higher_CI_cont, 
                                       PPV_higher_CI_cont_logit_exp, NPV_higher_CI_cont_logit_exp, 
                                       PPV_higher_CI_exp, NPV_higher_CI_exp, 
                                       PPV_higher_CI_glm_single, NPV_higher_CI_glm_single,
                                       PPV_higher_CI_robust, NPV_higher_CI_robust, 
                                       PPV_higher_CI_robust_single, NPV_higher_CI_robust_single, 
                                       PPV_higher_CI_robust_cl, NPV_higher_CI_robust_cl, 
                                       PPV_higher_CI_robust_cl_single, NPV_higher_CI_robust_cl_single, 
                                       PPV_higher_CI_glmer, NPV_higher_CI_glmer, 
                                       PPV_higher_CI_glmer_single, NPV_higher_CI_glmer_single, 
                                       PPV_higher_CI_gee, NPV_higher_CI_gee, 
                                       PPV_higher_CI_gee_single, NPV_higher_CI_gee_single,
                                       PPV_higher_CI_boot, NPV_higher_CI_boot, 
                                       PPV_higher_CI_boot_single, NPV_higher_CI_boot_single))



# convert estimates and CI's to percentage
results_id_sensspec$Estimate <- results_id_sensspec$Estimate*100
results_id_sensspec$Lower_CI <- results_id_sensspec$Lower_CI*100
results_id_sensspec$Higher_CI <- results_id_sensspec$Higher_CI*100

results_id_PV$Estimate <- results_id_PV$Estimate*100
results_id_PV$Lower_CI <- results_id_PV$Lower_CI*100
results_id_PV$Higher_CI <- results_id_PV$Higher_CI*100

results_sensspec$Estimate <- results_sensspec$Estimate*100
results_sensspec$Lower_CI <- results_sensspec$Lower_CI*100
results_sensspec$Higher_CI <- results_sensspec$Higher_CI*100

results_PV$Estimate <- results_PV$Estimate*100
results_PV$Lower_CI <- results_PV$Lower_CI*100
results_PV$Higher_CI <- results_PV$Higher_CI*100

# separate results tables for each diagnostic
results_id_sens <- results_id_sensspec[seq(1, nrow(results_id_sensspec), 2), ]
results_id_spec <- results_id_sensspec[seq(2, nrow(results_id_sensspec), 2), ]
results_id_PPV <- results_id_PV[seq(1, nrow(results_id_PV), 2), ]
results_id_NPV <- results_id_PV[seq(2, nrow(results_id_PV), 2), ]

results_sens <- results_sensspec[seq(1, nrow(results_sensspec), 2), ]
results_spec <- results_sensspec[seq(2, nrow(results_sensspec), 2), ]
results_PPV <- results_PV[seq(1, nrow(results_PV), 2), ]
results_NPV <- results_PV[seq(2, nrow(results_PV), 2), ]


# # storing ICC info
ICC_info <- as.data.frame(cbind(ICC_sens, ICC_spec))
colnames(ICC_info) <- c("ICC_within_diseased_population", "ICC_within_nondiseased_population")

# 
# cbinding sens and spec results
names(results_id_sens)[3] <- "Sensitivity (%)"
names(results_id_spec)[3] <- "Specificity (%)"
results_id_sensspec_c <- cbind(as.data.frame(results_id_sens[,-2]), as.data.frame(results_id_spec[,3:length(results_id_spec)]))

names(results_sens)[3] <- "Sensitivity (%)"
names(results_spec)[3] <- "Specificity (%)"
results_sensspec_c <- cbind(as.data.frame(results_sens[,-2]), as.data.frame(results_spec[,3:length(results_spec)]))

# cbinding PPV and NPV results
names(results_id_PPV)[3] <- "PPV (%)"
names(results_id_NPV)[3] <- "NPV (%)"
results_id_PV_c <- cbind(as.data.frame(results_id_PPV[,-2]), as.data.frame(results_id_NPV[,3:length(results_id_NPV)]))

names(results_PPV)[3] <- "PPV (%)"
names(results_NPV)[3] <- "NPV (%)"
results_PV_c <- cbind(as.data.frame(results_PPV[,-2]), as.data.frame(results_NPV[,3:length(results_NPV)]))


# resetting row numbers
rownames(results_id_sensspec_c) = seq(length=nrow(results_id_sensspec_c))
rownames(results_sensspec_c) = seq(length=nrow(results_sensspec_c))
rownames(results_id_PV_c) = seq(length=nrow(results_id_PV_c))
rownames(results_PV_c) = seq(length=nrow(results_PV_c))


# Start a sink file with a CSV extension
sink("Results.csv")

# Write the first dataframe, with a title and final line separator 
cat('Patient-level analyses')
cat('\n')
cat('\n')
cat('Patient-level contingency table')
cat('\n')
cat('____________________________')
cat('\n')
write.csv(id_cont_table)
cat('____________________________')
cat('\n')
cat('\n')
cat('Prevalence-independent measures')
cat('\n')
write.csv(results_id_sensspec_c)
cat('\n')
cat('____________________________')
cat('\n')
cat('\n')
cat('Prevalence-dependent measures')
cat('\n')
write.csv(results_id_PV_c)
cat('____________________________')
cat('\n')
cat('\n')
cat('Patient-level likelihood ratios')
cat('\n')
write.csv(id_blr)
cat('95% CIs computed via BCa bootstrapping.') 
cat('\n')
cat('____________________________')
cat('\n')
cat('\n')


# Write the 2nd dataframe to the same sink
cat('Segment-level analyses')
cat('\n')
cat('\n')
cat('Segment-level contingency table')
cat('\n')
cat('____________________________')
cat('\n')
write.csv(cont_table)
cat('____________________________')
cat('\n')
cat('\n')
cat('Prevalence-independent measures')
cat('\n')
write.csv(results_sensspec_c)
cat('\n')
cat('____________________________')
cat('\n')
cat('\n')
cat('Prevalence-dependent measures')
cat('\n')
write.csv(results_PV_c)
cat('____________________________')
cat('\n')
cat('\n')
cat('Segment-level likelihood ratios')
cat('\n')
cat('NB: CIs not adjusted for clustering!')
cat('\n')
write.csv(blr)
cat('95% CIs computed via BCa bootstrapping.') 
cat('\n')
cat('____________________________')
cat('\n')
cat('\n')

cat('Intracluster correlation coefficients')
cat('\n')
write.csv(ICC_info)
cat('____________________________')
cat('\n')

# Close the sink
sink()


########################
#### Visualizations ####
########################

######################
#### Forest plots ####
######################

# renaming the point estimate
results_id_sens$Sensitivity <- results_id_sens$Estimate
results_id_spec$Specificity <- results_id_spec$Estimate
results_id_PPV$PPV <- results_id_PPV$Estimate
results_id_NPV$NPV <- results_id_NPV$Estimate
results_sens$Sensitivity <- results_sens$Estimate
results_spec$Specificity <- results_spec$Estimate
results_PPV$PPV <- results_PPV$Estimate
results_NPV$NPV <- results_NPV$Estimate

        
# initialize plot
png("1_id_sensspec_forestplot.png", units = "in", width = 9, height = 10, res = 300)
# patient-level sens and spec
 forestplot(mean=cbind(results_id_sens$Sensitivity, results_id_spec$Specificity), 
           lower=cbind(results_id_sens$Lower_CI, results_id_spec$Lower_CI), 
           upper=cbind(results_id_sens$Higher_CI, results_id_spec$Higher_CI), 
           labeltext=list(paste(results_id_sens$Type_of_method)),
           legend=c("Sensitivity", "Specificity"),
           zero = 100,
           #zero = c(True_id_sens*100, True_id_spec*100),
           legend.pos=list(x=0.8,y=.4),
           legend.gp = gpar(col="#AAAAAA"), 
           legend.r=unit(.1, "snpc"),
           #clip=c(-.2, .2), 
           xticks=c(30, 40, 50, 60, 70, 80, 90, 100),
           boxsize=0.1,
           ci.vertices=T,
           col = fpColors(
                   box = c("darkred", "blue"),
                   line = c("darkred", "blue"), 
                   zero="black"),
           pch=13,
           grid=T,
           lwd.ci=3,
           # Set the different functions
           fn.ci_norm=c("fpDrawNormalCI", "fpDrawCircleCI"),
           xlab="%",
           title = "Patient-level prevalence-independent measures",
           new_page=TRUE)
dev.off()


# patient-level PPV and NPV
png("2_id_PV_forestplot.png", units = "in", width = 9, height = 10, res = 300)
forestplot(mean=cbind(results_id_PPV$PPV, results_id_NPV$NPV), 
           lower=cbind(results_id_PPV$Lower_CI, results_id_NPV$Lower_CI), 
           upper=cbind(results_id_PPV$Higher_CI, results_id_NPV$Higher_CI), 
           labeltext=list(paste(results_id_PPV$Type_of_method)),
           legend=c("PPV", "NPV"), 
           zero = 100,
           #zero = c(True_id_PPV*100, True_id_NPV*100),
           legend.pos=list(x=0.8,y=.4),
           legend.gp = gpar(col="#AAAAAA"), 
           legend.r=unit(.1, "snpc"),
           #clip=c(-.2, .2), 
           xticks=c(30, 40, 50, 60, 70, 80, 90, 100),
           boxsize=0.15,
           ci.vertices=T,
           col = fpColors(
                   box = c("purple", "darkgreen"),
                   line = c("purple", "darkgreen"), 
                   zero="black"),
           pch=13,
           grid=T,
           lwd.ci=3,
           # Set the different functions
           fn.ci_norm=c("fpDrawDiamondCI", "fpDrawPointCI"),
           xlab="%",
           title = "Patient-level prevalence-dependent measures",
           new_page=TRUE)
dev.off()


# segment-level sens and spec
png("3_sensspec_forestplot.png", units = "in", width = 9, height = 10, res = 300)
forestplot(mean=cbind(results_sens$Sensitivity, results_spec$Specificity), 
           lower=cbind(results_sens$Lower_CI, results_spec$Lower_CI), 
           upper=cbind(results_sens$Higher_CI, results_spec$Higher_CI), 
           labeltext=list(paste(results_sens$Type_of_method)),
           legend=c("Sensitivity", "Specificity"), 
           #zero = c(True_sample_sens*100, True_sample_spec*100),
           zero = 100,
           legend.pos=list(x=0.8,y=.4),
           legend.gp = gpar(col="#AAAAAA"), 
           legend.r=unit(.1, "snpc"),
           #clip=c(-.2, .2), 
           xticks=c(30, 40, 50, 60, 70, 80, 90, 100),
           boxsize=0.15,
           ci.vertices=T,
           col = fpColors(
                   box = c("darkred", "blue"),
                   line = c("darkred", "blue"), 
                   zero=c("black")),
           pch=13,
           grid=T,
           lwd.ci=3,
           # Set the different functions
           fn.ci_norm=c("fpDrawNormalCI", "fpDrawCircleCI"),
           xlab="%",
           title = "Segment-level prevalence-independent measures",
           new_page=TRUE)
dev.off()

# segment-level PPV and NPV
png("4_PV_forestplot.png", units = "in", width = 9, height = 10, res = 300)
forestplot(mean=cbind(results_PPV$PPV, results_NPV$NPV), 
           lower=cbind(results_PPV$Lower_CI, results_NPV$Lower_CI), 
           upper=cbind(results_PPV$Higher_CI, results_NPV$Higher_CI), 
           labeltext=list(paste(results_PPV$Type_of_method)),
           legend=c("PPV", "NPV"), 
           #zero = c(True_sample_PPV*100, True_sample_NPV*100),
           zero = 100,
           legend.pos=list(x=0.8,y=.4),
           legend.gp = gpar(col="#AAAAAA"), 
           legend.r=unit(.1, "snpc"),
           #clip=c(-.2, .2), 
           xticks=c(30, 40, 50, 60, 70, 80, 90, 100),
           boxsize=0.15,
           ci.vertices=T,
           col = fpColors(
                   box = c("purple", "darkgreen"),
                   line = c("purple", "darkgreen"), 
                   zero="black"),
           pch=13,
           grid=T,
           lwd.ci=3,
           # Set the different functions
           fn.ci_norm=c("fpDrawDiamondCI", "fpDrawPointCI"),
           xlab="%",
           title = "Segment-level prevalence-dependent measures",
           new_page=TRUE)
dev.off()

# emitting beep sound to alert you when script finished running
beep(1)

###################
### SCRIPT ENDS ###
###################

