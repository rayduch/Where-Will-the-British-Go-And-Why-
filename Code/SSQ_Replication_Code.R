################################################################################
##                                                                            ##
##                    §Conjoint Analysis Replication Code                     ##
##                                  SSQ Paper                                 ##
##                                                                            ##
##  Raymond Duch, Denise Laroze, Constantin Reinprecht, and Thomas Robinson   ##
##                                                                            ##
##                                11th April 2018                             ##
##                                                                            ##
################################################################################

library(openxlsx)
library(xtable)
library(gridExtra)
library(coefplot)
library(dummies)
library(nnet)
library(stargazer)
theme_set(theme_bw())
setwd("/Users/tomrobinson/OneDrive/CESS/Trump Conjoint/Data/SSQ Paper/")
#setwd("/Users/Denise Laroze Prehn/Dropbox/CESS-Santiago/archive/Trump/Data")

# Include code to replace standard errors with robust S.Es
robustse.f <- function(model, cluster, df_correction) {
  ## Huber-White heteroskedasticity-robust standard error calculation and table generation code for lm and glm models in R.
  ##Written by Joshua Gubler ~  http://joshuagubler.com.  Note that the second half of this function is just a wrapper for the excellent "multiwaycov" package here: https://cran.r-project.org/web/packages/multiwayvcov/multiwayvcov.pdf .  Love the authors of that package...
  ##Last updated: 3 July 2017  
  
  #model = the model you estimated, now calculated with robust standard errors
  #cluster = the name of the variable on which you will cluster. Put a tilda in front of it (e.g. ~ state).  If you don't put in a cluster, you will simply get huber-white robust errors.
  #df_correction: If you do not want the number of levels in your cluster variable to count against your degrees of freedom (like the xt- options in Stata), then type "F".  Otherwise, type "T" and these levels will be counted against your degrees of freedom
  
  require(sandwich)
  require(lmtest)
  require(multiwayvcov)
  if(missing(cluster)) {
    name <- deparse(substitute(model))
    modelname <- paste(name,"rob",sep=".")
    model$se <- coeftest(model, vcov=vcovHC(model,"HC1"))[,2]
    model$vcovHC <- vcovHC(model,"HC1")
    assign(modelname,model,envir = .GlobalEnv)
    coeftest(model, vcov=vcovHC(model,"HC1"))
  } else {
    name <- deparse(substitute(model))
    vcovCL <- cluster.vcov(model, cluster, df_correction = df_correction)
    model$vcovCL <- vcovCL
    modelname <- paste(name,"clustrob",sep=".")
    model$se <- coeftest(model, vcovCL)[,2]
    assign(modelname,model,envir = .GlobalEnv)
    #coeftest(model, vcovCL)
  }
}

#### REPLICATION CODE ####

conjoint1 <- read.csv("conjoint1.csv")
conjoint2 <- read.csv("conjoint2.csv")
conjoint3 <- read.csv("conjoint3.csv")

#### Analysis - Full Sample ####
conjoint1 <- within(conjoint1, econ <- relevel(econ, ref = "Annual GDP Growth of 4%"))
conjoint1 <- within(conjoint1, service <- relevel(service, ref = "Average international ranking of service salaries: 70th Percentile"))
conjoint1 <- within(conjoint1, immigration <- relevel(immigration, ref = "Change in visa processing centres"))
conjoint1 <- within(conjoint1, education <- relevel(education, ref = "Average international ranking of universities: 60th Percentile"))

conjoint2 <- within(conjoint2, econ <- relevel(econ, ref = "Annual GDP Growth of 4%"))
conjoint2 <- within(conjoint2, service <- relevel(service, ref = "Average international ranking of service salaries: 70th Percentile"))
conjoint2 <- within(conjoint2, immigration <- relevel(immigration, ref = "Change in visa processing centres"))
conjoint2 <- within(conjoint2, education <- relevel(education, ref = "Average international ranking of universities: 60th Percentile"))

conjoint3 <- within(conjoint3, econ <- relevel(econ, ref = "Annual GDP Growth of 4%"))
conjoint3 <- within(conjoint3, service <- relevel(service, ref = "Average international ranking of service salaries: 70th Percentile"))
conjoint3 <- within(conjoint3, education <- relevel(education, ref = "Average international ranking of universities: 60th Percentile"))

# Logit models
model1 <- glm(destination ~ social + econ + service + immigration + education + likely,family=binomial(link='logit'),data=conjoint1)
model1_controls <- glm(destination ~ social + econ + service + immigration + education + age + gender + ideology + likely,family=binomial(link='logit'),data=conjoint1)

model2 <- glm(destination ~ social + econ + service + immigration + education+ likely,family=binomial(link='logit'),data=conjoint2)
model2_controls <- glm(destination ~ social + econ + service + immigration + education + age + gender + ideology + likely,family=binomial(link='logit'),data=conjoint2)

model3 <- glm(destination ~ social + econ + service + country + education + likely,family=binomial(link='logit'),data=conjoint3)
model3_controls <- glm(destination ~ social + econ + service + country + education + age + gender + ideology + likely,family=binomial(link='logit'),data=conjoint3)
model3_controls_same1 <- glm(destination ~ social + econ + service + country + education + age + gender + ideology + likely,family=binomial(link='logit'),data=conjoint3[conjoint3$c.same == 1,])
model3_controls_same0 <- glm(destination ~ social + econ + service + country + education + age + gender + ideology + likely,family=binomial(link='logit'),data=conjoint3[conjoint3$c.same == 0,])


# Cluster-robust standard errors
robustse.f(model1, ~id, F)
robustse.f(model1_controls, ~id, F)
robustse.f(model2, ~id, F)
robustse.f(model2_controls, ~id, F)
robustse.f(model3, ~id, F)
robustse.f(model3_controls, ~id, F)
robustse.f(model3_controls_same0, ~id, F)
robustse.f(model3_controls_same1, ~id, F)


# Main table - LaTeX
stargazer(model1.clustrob,model2.clustrob, model3.clustrob,
          covariate.labels = c("Generous family allowance",
                               "No minimum wage or income support",
                               "GDP 2 percent",
                               "GDP 6 percent",
                               "Service salaries 50th pc",
                               "Service salaries 90th pc",
                               "Deportation of all illegal immigrants",
                               "Point-system visa",
                               "Muslim Ban",
                               "Canada",
                               "U.S.A.",
                               "University Ranking 40th pc",
                               "University Ranking 90th pc",
                               "Likelihood of emigrating"),
          dep.var.caption = "Model",
          dep.var.labels.include = FALSE,
          se = list(model1.clustrob$se,model2.clustrob$se, model3.clustrob$se),
          no.space = TRUE)

#LaTeX with controls
stargazer(model1_controls.clustrob,model2_controls.clustrob,model3_controls.clustrob,
          covariate.labels = c("Generous family allowance",
                               "No minimum wage or income support",
                               "GDP 2 percent",
                               "GDP 6 percent",
                               "Service salaries 50th pc",
                               "Service salaries 90th pc",
                               "Deportation of all illegal immigrants",
                               "Point-system visa",
                               "Muslim Ban",
                               "Canada",
                               "U.S.A.",
                               "University Ranking 40th pc",
                               "University Ranking 90th pc",
                               "Age",
                               "Gender: Male",
                               "Gender: Other",
                               "Ideological self-placement",
                               "Likelihood of emigrating"),
          dep.var.caption = "Model",
          dep.var.labels.include = FALSE,
          se=list(model1_controls.clustrob$se,model2_controls.clustrob$se,model3_controls.clustrob$se),
          no.space = TRUE)

# LaTeX with controls + same country breakout
stargazer(model3_controls_same0.clustrob,model3_controls_same1.clustrob,
          covariate.labels = c("Generous family allowance",
                               "No minimum wage or income support",
                               "GDP 2 percent",
                               "GDP 6 percent",
                               "Service salaries 50th pc",
                               "Service salaries 90th pc",
                               "Canada",
                               "U.S.A.",
                               "University Ranking 40th pc",
                               "University Ranking 90th pc",
                               "Age",
                               "Gender: Male",
                               "Gender: Other",
                               "Ideological self-placement",
                               "Likelihood of emigrating"),
          dep.var.caption = "Model",
          dep.var.labels.include = FALSE,
          se=list(model3_controls_same0.clustrob$se,model3_controls_same1.clustrob$se),
          no.space = TRUE)

## Coef plots

#Model 1
plotdf1 <- data.frame(estimate = coeftest(model1.clustrob, model1.clustrob$vcovCL)[,1],SE = coeftest(model1.clustrob, model1.clustrob$vcovCL)[,2])
plotdf1 <- plotdf1[-1,]
plotdf1 <- rbind(plotdf1[1,],
                 c(0,0),
                 plotdf1[2:11,])
plotdf1 <- rbind(plotdf1[1:4,],
                 c(0,0),
                 plotdf1[5:12,])
plotdf1 <- rbind(plotdf1[1:7,],
                 c(0,0),
                 plotdf1[8:13,])
plotdf1 <- rbind(plotdf1[1:10,],
                 c(0,0),
                 plotdf1[11:14,])
plotdf1 <- rbind(plotdf1[1:13,],
                 c(0,0),
                 plotdf1[14:15,])
plotdf1$coef <- c("Generous family allowance",
                  "Basic minimum wage",
                  "No minimum wage or income support",
                  "GDP 2 percent",
                  "GDP 4 percent",
                  "GDP 6 percent",
                  "Service salaries 50th pc",
                  "Service salaries 70th pc",
                  "Service salaries 90th pc",
                  "Point-system visa",
                  "Change in visa processing centres",
                  "Muslim Ban",
                  "University Ranking 40th pc",
                  "University Ranking 60th pc",
                  "University Ranking 90th pc",
                  "Likelihood of emigrating")

plotdf1$coef <- factor(plotdf1$coef, levels = c("Likelihood of emigrating",
                                                "University Ranking 90th pc",
                                                "University Ranking 60th pc",
                                                "University Ranking 40th pc",
                                                "Muslim Ban",
                                                "Change in visa processing centres",
                                                "Point-system visa",
                                                "Service salaries 90th pc",
                                                "Service salaries 70th pc",
                                                "Service salaries 50th pc",
                                                "GDP 6 percent",
                                                "GDP 4 percent",
                                                "GDP 2 percent",
                                                "No minimum wage or income support",
                                                "Basic minimum wage",
                                                "Generous family allowance",
                                                "Intercept"))


plotdf1$LCI <- plotdf1$estimate-1.96*plotdf1$SE
plotdf1$UCI <- plotdf1$estimate+1.96*plotdf1$SE



ggplot(plotdf1, aes(x=coef)) +
  geom_point(y=plotdf1$estimate) +
  geom_linerange(max=plotdf1$UCI, min=plotdf1$LCI) +
  geom_hline(yintercept=0, linetype = "dashed") +
  ylim(-1,1) +
  labs(x="",y="") +
  coord_flip()

ggsave(paste0("conjoint1.png" ,""), width = 15, height = 10, units = c("cm"), dpi = 300)

#Model 2
plotdf2 <- data.frame(estimate = coeftest(model2.clustrob, model2.clustrob$vcovCL)[,1],SE = coeftest(model2.clustrob, model2.clustrob$vcovCL)[,2])
plotdf2 <- plotdf2[-1,]
plotdf2 <- rbind(plotdf2[1,],
                 c(0,0),
                 plotdf2[2:11,])
plotdf2 <- rbind(plotdf2[1:4,],
                 c(0,0),
                 plotdf2[5:12,])
plotdf2 <- rbind(plotdf2[1:7,],
                 c(0,0),
                 plotdf2[8:13,])
plotdf2 <- rbind(plotdf2[1:10,],
                 c(0,0),
                 plotdf2[11:14,])
plotdf2 <- rbind(plotdf2[1:13,],
                 c(0,0),
                 plotdf2[14:15,])
plotdf2$coef <- c("Generous family allowance",
                  "Basic minimum wage",
                  "No minimum wage or income support",
                  "GDP 2 percent",
                  "GDP 4 percent",
                  "GDP 6 percent",
                  "Service salaries 50th pc",
                  "Service salaries 70th pc",
                  "Service salaries 90th pc",
                  "Deportation of all illegal immigrants",
                  "Change in visa processing centres",
                  "Point-system visa",
                  "University Ranking 40th pc",
                  "University Ranking 60th pc",
                  "University Ranking 90th pc",
                  "Likelihood of emigrating")

plotdf2$LCI <- plotdf2$estimate-1.96*plotdf2$SE
plotdf2$UCI <- plotdf2$estimate+1.96*plotdf2$SE

plotdf2$coef <- factor(plotdf2$coef, levels = c("Likelihood of emigrating",
                                                "University Ranking 90th pc",
                                                "University Ranking 60th pc",
                                                "University Ranking 40th pc",
                                                "Deportation of all illegal immigrants",
                                                "Change in visa processing centres",
                                                "Point-system visa",
                                                "Service salaries 90th pc",
                                                "Service salaries 70th pc",
                                                "Service salaries 50th pc",
                                                "GDP 6 percent",
                                                "GDP 4 percent",
                                                "GDP 2 percent",
                                                "No minimum wage or income support",
                                                "Basic minimum wage",
                                                "Generous family allowance",
                                                "Intercept"))

ggplot(plotdf2, aes(x=coef)) +
  geom_point(y=plotdf2$estimate) +
  geom_linerange(max=plotdf2$UCI, min=plotdf2$LCI) +
  geom_hline(yintercept=0, linetype = "dashed") +
  ylim(-1,1) +
  labs(x="",y="") +
  coord_flip()

ggsave(paste0("conjoint2.png" ,""), width = 15, height = 10, units = c("cm"), dpi = 300)

#Model 3
plotdf3 <- data.frame(estimate = coeftest(model3.clustrob, model3.clustrob$vcovCL)[,1],SE = coeftest(model3.clustrob, model3.clustrob$vcovCL)[,2])
plotdf3 <- plotdf3[-1,]
plotdf3 <- rbind(plotdf3[1,],
                 c(0,0),
                 plotdf3[2:11,])
plotdf3 <- rbind(plotdf3[1:4,],
                 c(0,0),
                 plotdf3[5:12,])
plotdf3 <- rbind(plotdf3[1:7,],
                 c(0,0),
                 plotdf3[8:13,])
plotdf3 <- rbind(plotdf3[1:10,],
                 c(0,0),
                 plotdf3[11:14,])
plotdf3 <- rbind(plotdf3[1:13,],
                 c(0,0),
                 plotdf3[14:15,])
plotdf3$coef <- c("Generous family allowance",
                  "Basic minimum wage",
                  "No minimum wage or income support",
                  "GDP 2 percent",
                  "GDP 4 percent",
                  "GDP 6 percent",
                  "Service salaries 50th pc",
                  "Service salaries 70th pc",
                  "Service salaries 90th pc",
                  "Canada",
                  "Australia",
                  "U.S.A.",
                  "University Ranking 40th pc",
                  "University Ranking 60th pc",
                  "University Ranking 90th pc",
                  "Likelihood of emigrating")

plotdf3$LCI <- plotdf3$estimate-1.96*plotdf3$SE
plotdf3$UCI <- plotdf3$estimate+1.96*plotdf3$SE

plotdf3$coef <- factor(plotdf3$coef, levels = c("Likelihood of emigrating",
                                                "University Ranking 90th pc",
                                                "University Ranking 60th pc",
                                                "University Ranking 40th pc",
                                                "U.S.A.",
                                                "Australia",
                                                "Canada",
                                                "Service salaries 90th pc",
                                                "Service salaries 70th pc",
                                                "Service salaries 50th pc",
                                                "GDP 6 percent",
                                                "GDP 4 percent",
                                                "GDP 2 percent",
                                                "No minimum wage or income support",
                                                "Basic minimum wage",
                                                "Generous family allowance",
                                                "Intercept"))

ggplot(plotdf3, aes(x=coef)) +
  geom_point(y=plotdf3$estimate) +
  geom_linerange(max=plotdf3$UCI, min=plotdf3$LCI) +
  geom_hline(yintercept=0, linetype = "dashed") +
  ylim(-1,1) +
  labs(x="",y="") +
  coord_flip()

ggsave(paste0("conjoint3.png" ,""), width = 15, height = 10, units = c("cm"), dpi = 300)

# Combine into one plot

plotdf1 <- cbind(plotdf1,treatment = "Treatment: 1")
plotdf2 <- cbind(plotdf2,treatment = "Treatment: 2")
plotdf3 <- cbind(plotdf3,treatment = "Treatment: 3")

plotdf <- rbind(plotdf1,plotdf2,plotdf3)
plotdf$coef <- factor(plotdf$coef, levels = c("Likelihood of emigrating",
                                              "University Ranking 90th pc",
                                              "University Ranking 60th pc",
                                              "University Ranking 40th pc",
                                              "U.S.A.",
                                              "Australia",
                                              "Canada",
                                              "Deportation of all illegal immigrants",
                                              "Muslim Ban",
                                              "Change in visa processing centres",
                                              "Point-system visa",
                                              "Service salaries 90th pc",
                                              "Service salaries 70th pc",
                                              "Service salaries 50th pc",
                                              "GDP 6 percent",
                                              "GDP 4 percent",
                                              "GDP 2 percent",
                                              "No minimum wage or income support",
                                              "Basic minimum wage",
                                              "Generous family allowance",
                                              "Intercept"))

ggplot(data = plotdf, aes(x=coef)) +
  facet_grid(.~treatment) +
  geom_point(aes(y=estimate)) +
  geom_linerange(aes(max=UCI, min=LCI)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  labs(x="",y="") +
  coord_flip()

ggsave(paste0("conjoint_combined.png"), width = 15, height = 20, units = c("cm"), dpi = 300)


#### Left-Right Divide (+ Gender Breakout) Analysis ####

## Generate models
# Set reference categories
conjoint1 <- within(conjoint1, econ <- relevel(econ, ref = "Annual GDP Growth of 2%"))
conjoint1 <- within(conjoint1, service <- relevel(service, ref = "Average international ranking of service salaries: 50th Percentile"))
conjoint1 <- within(conjoint1, education <- relevel(education, ref = "Average international ranking of universities: 40th Percentile"))

conjoint2 <- within(conjoint2, econ <- relevel(econ, ref = "Annual GDP Growth of 2%"))
conjoint2 <- within(conjoint2, service <- relevel(service, ref = "Average international ranking of service salaries: 50th Percentile"))
conjoint2 <- within(conjoint2, education <- relevel(education, ref = "Average international ranking of universities: 40th Percentile"))

conjoint3 <- within(conjoint3, econ <- relevel(econ, ref = "Annual GDP Growth of 2%"))
conjoint3 <- within(conjoint3, service <- relevel(service, ref = "Average international ranking of service salaries: 50th Percentile"))
conjoint3 <- within(conjoint3, education <- relevel(education, ref = "Average international ranking of universities: 40th Percentile"))

# Logit models
model1_left <- glm(destination ~ social + econ + service + immigration + education + likely,family=binomial(link='logit'),data=conjoint1[conjoint1$ideology < 5,])
model1_right <- glm(destination ~ social + econ + service + immigration + education+ likely,family=binomial(link='logit'),data=conjoint1[conjoint1$ideology > 5,])
model1_centre <- glm(destination ~ social + econ + service + immigration + education + likely,family=binomial(link='logit'),data=conjoint1[conjoint1$ideology == 5,])
model1_male <- glm(destination ~ social + econ + service + immigration + education + likely,family=binomial(link='logit'),data=conjoint1[conjoint1$gender == "Male",])
model1_female <- glm(destination ~ social + econ + service + immigration + education + likely,family=binomial(link='logit'),data=conjoint1[conjoint1$gender == "Female",])

model2_left <- glm(destination ~ social + econ + service + immigration + education + likely,family=binomial(link='logit'),data=conjoint2[conjoint2$ideology < 5,])
model2_right <- glm(destination ~ social + econ + service + immigration + education + likely,family=binomial(link='logit'),data=conjoint2[conjoint2$ideology > 5,])
model2_centre <- glm(destination ~ social + econ + service + immigration + education + likely,family=binomial(link='logit'),data=conjoint2[conjoint2$ideology == 5,])
model2_male <- glm(destination ~ social + econ + service + immigration + education + likely,family=binomial(link='logit'),data=conjoint2[conjoint2$gender == "Male",])
model2_female <- glm(destination ~ social + econ + service + immigration + education + likely,family=binomial(link='logit'),data=conjoint2[conjoint2$gender == "Female",])

model3_left <- glm(destination ~ social + econ + service + country + education + likely,family=binomial(link='logit'),data=conjoint3[conjoint3$ideology < 5,])
model3_right <- glm(destination ~ social + econ + service + country + education + likely,family=binomial(link='logit'),data=conjoint3[conjoint3$ideology > 5,])
model3_centre <- glm(destination ~ social + econ + service + country + education + likely,family=binomial(link='logit'),data=conjoint3[conjoint3$ideology == 5,])
model3_male <- glm(destination ~ social + econ + service + country + education + likely,family=binomial(link='logit'),data=conjoint3[conjoint3$gender == "Male",])
model3_female <- glm(destination ~ social + econ + service + country + education + likely,family=binomial(link='logit'),data=conjoint3[conjoint3$gender == "Female",])

# Cluster-robust standard errors
robustse.f(model1_left, ~id, F)
robustse.f(model1_centre, ~id, F)
robustse.f(model1_right, ~id, F)
robustse.f(model1_male, ~id, F)
robustse.f(model1_female, ~id, F)
robustse.f(model2_left, ~id, F)
robustse.f(model2_centre, ~id, F)
robustse.f(model2_right, ~id, F)
robustse.f(model2_male, ~id, F)
robustse.f(model2_female, ~id, F)
robustse.f(model3_left, ~id, F)
robustse.f(model3_centre, ~id, F)
robustse.f(model3_right, ~id, F)
robustse.f(model3_male, ~id, F)
robustse.f(model3_female, ~id, F)

## LaTeX Tables

#LaTeX Treatment 1 with breakouts
stargazer(model1_male.clustrob,model1_female.clustrob,model1_left.clustrob,model1_centre.clustrob,model1_right.clustrob,
          covariate.labels = c("Generous family allowance",
                               "No minimum wage or income support",
                               "GDP 4 percent",
                               "GDP 6 percent",
                               "Service salaries 70th pc",
                               "Service salaries 90th pc",
                               "Point-system visa",
                               "Muslim Ban",
                               "University Ranking 60th pc",
                               "University Ranking 90th pc",
                               "Likelihood of emigrating"),
          dep.var.caption = "Model",
          dep.var.labels.include = FALSE,
          se=list(model1_male.clustrob$se,model1_female.clustrob$se,model1_left.clustrob$se,model1_centre.clustrob$se,model1_right.clustrob$se),
          no.space = TRUE)

#LaTeX Treatment 2 with breakouts
stargazer(model2_male.clustrob,model2_female.clustrob,model2_left.clustrob,model2_centre.clustrob,model2_right.clustrob,
          covariate.labels = c("Generous family allowance",
                               "No minimum wage or income support",
                               "GDP 4 percent",
                               "GDP 6 percent",
                               "Service salaries 70th pc",
                               "Service salaries 90th pc",
                               "Deportation of all illegal immigrants",
                               "Point-system visa",
                               "University Ranking 60th pc",
                               "University Ranking 90th pc",
                               "Likelihood of emigrating"),
          dep.var.caption = "Model",
          dep.var.labels.include = FALSE,
          se = list(model2_male.clustrob$se,model2_female.clustrob$se,model2_left.clustrob$se,model2_centre.clustrob$se,model2_right.clustrob$se),
          no.space = TRUE)

#LaTeX Treatment 3 with breakouts
stargazer(model3_male.clustrob,model3_female.clustrob,model3_left.clustrob,model3_centre.clustrob,model3_right.clustrob,
          covariate.labels = c("Generous family allowance",
                               "No minimum wage or income support",
                               "GDP 4 percent",
                               "GDP 6 percent",
                               "Service salaries 70th pc",
                               "Service salaries 90th pc",
                               "Canada",
                               "U.S.A.",
                               "University Ranking 60th pc",
                               "University Ranking 90th pc",
                               "Likelihood of emigrating"),
          dep.var.caption = "Model",
          dep.var.labels.include = FALSE,
          se = list(model3_male.clustrob$se,model3_female.clustrob$se,model3_left.clustrob$se,model3_centre.clustrob$se,model3_right.clustrob$se),
          no.space = TRUE)

## Coef plot
# Left group
plotdf_l <- data.frame(estimate = coeftest(model1_left.clustrob, model1_left.clustrob$vcovCL)[,1],SE = coeftest(model1_left.clustrob, model1_left.clustrob$vcovCL)[,2],group = "Left",treatment = 1)
plotdf_l <- rbind(plotdf_l, data.frame(estimate = coeftest(model2_left.clustrob, model2_left.clustrob$vcovCL)[,1],SE = coeftest(model2_left.clustrob, model2_left.clustrob$vcovCL)[,2],group = "Left",treatment = 2))
plotdf_l <- rbind(plotdf_l, data.frame(estimate = coeftest(model3_left.clustrob, model3_left.clustrob$vcovCL)[,1],SE = coeftest(model3_left.clustrob, model3_left.clustrob$vcovCL)[,2],group = "Left",treatment = 3))

plotdf_l <- rbind(plotdf_l[1:2,],c(0,0,"Left",1),plotdf_l[3:36,])
plotdf_l <- rbind(plotdf_l[1:4,],c(0,0,"Left",1),plotdf_l[5:37,])
plotdf_l <- rbind(plotdf_l[1:7,],c(0,0,"Left",1),plotdf_l[8:38,])
plotdf_l <- rbind(plotdf_l[1:11,],c(0,0,"Left",1),plotdf_l[12:39,])
plotdf_l <- rbind(plotdf_l[1:13,],c(0,0,"Left",1),plotdf_l[14:40,])
plotdf_l <- rbind(plotdf_l[1:19,],c(0,0,"Left",2),plotdf_l[20:41,])
plotdf_l <- rbind(plotdf_l[1:21,],c(0,0,"Left",2),plotdf_l[22:42,])
plotdf_l <- rbind(plotdf_l[1:24,],c(0,0,"Left",2),plotdf_l[25:43,])
plotdf_l <- rbind(plotdf_l[1:28,],c(0,0,"Left",2),plotdf_l[29:44,])
plotdf_l <- rbind(plotdf_l[1:30,],c(0,0,"Left",2),plotdf_l[31:45,])
plotdf_l <- rbind(plotdf_l[1:36,],c(0,0,"Left",3),plotdf_l[37:46,])
plotdf_l <- rbind(plotdf_l[1:38,],c(0,0,"Left",3),plotdf_l[39:47,])
plotdf_l <- rbind(plotdf_l[1:41,],c(0,0,"Left",3),plotdf_l[42:48,])
plotdf_l <- rbind(plotdf_l[1:45,],c(0,0,"Left",3),plotdf_l[46:49,])
plotdf_l <- rbind(plotdf_l[1:47,],c(0,0,"Left",3),plotdf_l[48:50,])

# Right group
plotdf_r <- data.frame(estimate = coeftest(model1_right.clustrob, model1_right.clustrob$vcovCL)[,1],SE = coeftest(model1_right.clustrob, model1_right.clustrob$vcovCL)[,2],group = "Right",treatment = 1)
plotdf_r <- rbind(plotdf_r, data.frame(estimate = coeftest(model2_right.clustrob, model2_right.clustrob$vcovCL)[,1],SE = coeftest(model2_right.clustrob, model2_right.clustrob$vcovCL)[,2],group = "Right",treatment = 2))
plotdf_r <- rbind(plotdf_r, data.frame(estimate = coeftest(model3_right.clustrob, model3_right.clustrob$vcovCL)[,1],SE = coeftest(model3_right.clustrob, model3_right.clustrob$vcovCL)[,2],group = "Right",treatment = 3))

plotdf_r <- rbind(plotdf_r[1:2,],c(0,0,"Right",1),plotdf_r[3:36,])
plotdf_r <- rbind(plotdf_r[1:4,],c(0,0,"Right",1),plotdf_r[5:37,])
plotdf_r <- rbind(plotdf_r[1:7,],c(0,0,"Right",1),plotdf_r[8:38,])
plotdf_r <- rbind(plotdf_r[1:11,],c(0,0,"Right",1),plotdf_r[12:39,])
plotdf_r <- rbind(plotdf_r[1:13,],c(0,0,"Right",1),plotdf_r[14:40,])
plotdf_r <- rbind(plotdf_r[1:19,],c(0,0,"Right",2),plotdf_r[20:41,])
plotdf_r <- rbind(plotdf_r[1:21,],c(0,0,"Right",2),plotdf_r[22:42,])
plotdf_r <- rbind(plotdf_r[1:24,],c(0,0,"Right",2),plotdf_r[25:43,])
plotdf_r <- rbind(plotdf_r[1:28,],c(0,0,"Right",2),plotdf_r[29:44,])
plotdf_r <- rbind(plotdf_r[1:30,],c(0,0,"Right",2),plotdf_r[31:45,])
plotdf_r <- rbind(plotdf_r[1:36,],c(0,0,"Right",3),plotdf_r[37:46,])
plotdf_r <- rbind(plotdf_r[1:38,],c(0,0,"Right",3),plotdf_r[39:47,])
plotdf_r <- rbind(plotdf_r[1:41,],c(0,0,"Right",3),plotdf_r[42:48,])
plotdf_r <- rbind(plotdf_r[1:45,],c(0,0,"Right",3),plotdf_r[46:49,])
plotdf_r <- rbind(plotdf_r[1:47,],c(0,0,"Right",3),plotdf_r[48:50,])

# Centre group
plotdf_c <- data.frame(estimate = coeftest(model1_centre.clustrob, model1_centre.clustrob$vcovCL)[,1],SE = coeftest(model1_centre.clustrob, model1_centre.clustrob$vcovCL)[,2],group = "Centre",treatment = 1)
plotdf_c <- rbind(plotdf_c, data.frame(estimate = coeftest(model2_centre.clustrob, model2_centre.clustrob$vcovCL)[,1],SE = coeftest(model2_centre.clustrob, model2_centre.clustrob$vcovCL)[,2],group = "Centre",treatment = 2))
plotdf_c <- rbind(plotdf_c, data.frame(estimate = coeftest(model3_centre.clustrob, model3_centre.clustrob$vcovCL)[,1],SE = coeftest(model3_centre.clustrob, model3_centre.clustrob$vcovCL)[,2],group = "Centre",treatment = 3))

plotdf_c <- rbind(plotdf_c[1:2,],c(0,0,"Centre",1),plotdf_c[3:36,])
plotdf_c <- rbind(plotdf_c[1:4,],c(0,0,"Centre",1),plotdf_c[5:37,])
plotdf_c <- rbind(plotdf_c[1:7,],c(0,0,"Centre",1),plotdf_c[8:38,])
plotdf_c <- rbind(plotdf_c[1:11,],c(0,0,"Centre",1),plotdf_c[12:39,])
plotdf_c <- rbind(plotdf_c[1:13,],c(0,0,"Centre",1),plotdf_c[14:40,])
plotdf_c <- rbind(plotdf_c[1:19,],c(0,0,"Centre",2),plotdf_c[20:41,])
plotdf_c <- rbind(plotdf_c[1:21,],c(0,0,"Centre",2),plotdf_c[22:42,])
plotdf_c <- rbind(plotdf_c[1:24,],c(0,0,"Centre",2),plotdf_c[25:43,])
plotdf_c <- rbind(plotdf_c[1:28,],c(0,0,"Centre",2),plotdf_c[29:44,])
plotdf_c <- rbind(plotdf_c[1:30,],c(0,0,"Centre",2),plotdf_c[31:45,])
plotdf_c <- rbind(plotdf_c[1:36,],c(0,0,"Centre",3),plotdf_c[37:46,])
plotdf_c <- rbind(plotdf_c[1:38,],c(0,0,"Centre",3),plotdf_c[39:47,])
plotdf_c <- rbind(plotdf_c[1:41,],c(0,0,"Centre",3),plotdf_c[42:48,])
plotdf_c <- rbind(plotdf_c[1:45,],c(0,0,"Centre",3),plotdf_c[46:49,])
plotdf_c <- rbind(plotdf_c[1:47,],c(0,0,"Centre",3),plotdf_c[48:50,])


plotdf <- rbind(plotdf_l,plotdf_r,plotdf_c)

plotdf$coef <- c("Intercept",
                 "Generous family allowance",
                 "Basic minimum wage",
                 "No minimum wage or income support",
                 "GDP 2 percent",
                 "GDP 4 percent",
                 "GDP 6 percent",
                 "Service salaries 50th pc",
                 "Service salaries 70th pc",
                 "Service salaries 90th pc",
                 "Point-system visa",
                 "Change in visa processing centres",
                 "Muslim Ban",
                 "University Ranking 40th pc",
                 "University Ranking 60th pc",
                 "University Ranking 90th pc",
                 "Likelihood of emigrating",
                 "Intercept",
                 "Generous family allowance",
                 "Basic minimum wage",
                 "No minimum wage or income support",
                 "GDP 2 percent",
                 "GDP 4 percent",
                 "GDP 6 percent",
                 "Service salaries 50th pc",
                 "Service salaries 70th pc",
                 "Service salaries 90th pc",
                 "Deportation of all illegal immigrants",
                 "Change in visa processing centres",
                 "Point-system visa",
                 "University Ranking 40th pc",
                 "University Ranking 60th pc",
                 "University Ranking 90th pc",
                 "Likelihood of emigrating",
                 "Intercept",
                 "Generous family allowance",
                 "Basic minimum wage",
                 "No minimum wage or income support",
                 "GDP 2 percent",
                 "GDP 4 percent",
                 "GDP 6 percent",
                 "Service salaries 50th pc",
                 "Service salaries 70th pc",
                 "Service salaries 90th pc",
                 "Canada",
                 "Australia",
                 "U.S.A.",
                 "University Ranking 40th pc",
                 "University Ranking 60th pc",
                 "University Ranking 90th pc",
                 "Likelihood of emigrating")

plotdf$coef <- factor(plotdf$coef, levels = c("Likelihood of emigrating",
                                              "University Ranking 90th pc",
                                              "University Ranking 60th pc",
                                              "University Ranking 40th pc",
                                              "U.S.A.",
                                              "Australia",
                                              "Canada",
                                              "Muslim Ban",
                                              "Deportation of all illegal immigrants",
                                              "Change in visa processing centres",
                                              "Point-system visa",
                                              "Service salaries 90th pc",
                                              "Service salaries 70th pc",
                                              "Service salaries 50th pc",
                                              "GDP 6 percent",
                                              "GDP 4 percent",
                                              "GDP 2 percent",
                                              "No minimum wage or income support",
                                              "Basic minimum wage",
                                              "Generous family allowance",
                                              "Intercept"))

plotdf <- plotdf[plotdf$coef != "Intercept",]
plotdf$estimate <- as.numeric(plotdf$estimate)
plotdf$SE <- as.numeric(plotdf$SE)
plotdf$LCI <- plotdf$estimate-1.96*plotdf$SE
plotdf$UCI <- plotdf$estimate+1.96*plotdf$SE

plotdf$treatment <- factor(paste0("Treatment: ",plotdf$treatment))
plotdf$group <- factor(plotdf$group, levels=c("Left","Centre","Right"))

ggplot(data = plotdf, aes(x=coef, group=group)) +
  facet_grid(.~treatment) +
  geom_point(aes(y=estimate, colour = factor(group)),position=position_dodge(width = 0.7)) +
  geom_linerange(aes(max=UCI, min=LCI,colour = factor(group)),position=position_dodge(width = 0.7)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  scale_colour_manual(values=c("red","green","blue")) +
  labs(x="",y="") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(title="")) +
  coord_flip()

ggsave(paste0("conjoint_ideology.png"), width = 15, height = 20, units = c("cm"), dpi = 300)


#### Balance tests #####
conjoint3$c.same <- NULL
names(conjoint3)[names(conjoint3) == 'country'] <- 'immigration'
conjoint <- rbind(conjoint1,conjoint2,conjoint3)

conjoint <- within(conjoint, econ <- relevel(econ, ref = "Annual GDP Growth of 4%"))
conjoint <- within(conjoint, service <- relevel(service, ref = "Average international ranking of service salaries: 70th Percentile"))
conjoint <- within(conjoint, education <- relevel(education, ref = "Average international ranking of universities: 60th Percentile"))


bal_social <- multinom(social ~ age + gender + likely + ideology, data = conjoint)
bal_econ <- multinom(econ ~ age + gender + likely + ideology, data = conjoint)
bal_service <- multinom(service ~ age + gender + likely + ideology, data = conjoint)
bal_education <- multinom(education ~ age + gender + likely + ideology, data = conjoint)
bal_immigration <- multinom(immigration ~ age + gender + likely + ideology, data = conjoint)

stargazer(bal_social,
          covariate.labels = c("Age","Gender: Male","Gender: Other","Likelihood of emigrating"))
stargazer(bal_econ,
          covariate.labels = c("Age","Gender: Male","Gender: Other","Likelihood of emigrating"))
stargazer(bal_service,
          covariate.labels = c("Age","Gender: Male","Gender: Other","Likelihood of emigrating"))
stargazer(bal_education,
          covariate.labels = c("Age","Gender: Male","Gender: Other","Likelihood of emigrating"))
stargazer(bal_immigration,
          covariate.labels = c("Age","Gender: Male","Gender: Other","Likelihood of emigrating"))

#### Descriptive Stats ####
desc <- data.frame(Variable = as.character(),
                   Mean = as.double(),
                   SD = as.double(),
                   Min. = as.numeric(),
                   Max. = as.numeric(),
                   N = as.numeric(),
                   stringsAsFactors = FALSE)


age <- conjoint1$age[1:196]
aus <- conjoint1$aus[1:196]
can <- conjoint1$can[1:196]
us <- conjoint1$us[1:196]
interest <- conjoint1$interest[1:196]
likely <- conjoint1$likely[1:196]
ideology <- conjoint1$ideology[1:196]
gender <- conjoint1$gender[1:196]

tab1a <- list(age, aus, can,us, interest, likely, ideology)

for (i in tab1a) {
  row <- c("",
           round(mean(i, na.rm=TRUE),2),
           round(sd(i, na.rm=TRUE),2),
           min(i, na.rm=TRUE),
           max(i, na.rm=TRUE),
           length(i[!is.na(i)])
  )
  desc[nrow(desc)+1,] <- row
}
desc$Variable <- c("Age","Australia","Canada","U.S.A.","Interest","Likely","Ideology")

# Add gender
gender <- dummy(gender)
desc[nrow(desc)+1,] <- c("Female", round(mean(gender[,1]),2),round(sd(gender[,1]),2),0,1,length(gender[,1]))


#### Conjoint descriptives####
conjoint1 <- droplevels(conjoint1)
conjoint2 <- droplevels(conjoint2)
conjoint3 <- droplevels(conjoint3)

cj <- rbind(conjoint1, conjoint2, conjoint3)

cj <- cj[,-c(6:16)]

cj <- dummy.data.frame(cj)

dummyMeanscj <- round(as.double(colMeans(cj)),2)
dummyN <- as.numeric(colSums(cj))
dummySD <- round(as.double(apply(cj, 2, sd)),2)

var <- c("Basic hourly minimum wage",
         "Generous guaranteed monthly family allowance",
         "No state minimum wage or income support",
         "Annual GDP Growth of 2%",
         "Annual GDP Growth of 4%",
         "Annual GDP Growth of 6%",
         "Service salaries: 50th Percentile",
         "Service salaries: 70th Percentile",
         "Service salaries: 90th Percentile",
         "Change in visa processing centres", 
         "Implementation of point-system",
         "Restriction on Muslim immigration/tourist visas",
         "Deportation of all illegal immigrants",
         "Australia",
         "Canada",
         "U.S.A.",
         "Ranking of universities: 40th Percentile",
         "Ranking of universities: 60th Percentile",
         "Ranking of universities: 90th Percentile")
tab2 <- cbind(var,dummyMeanscj,dummySD,0,1,dummyN)
colnames(tab2) <- colnames(desc)

desc<-rbind(desc, tab2)
print(xtable(desc, digits=2),include.rownames=FALSE, digits=2)

# Check favourability score differences

t.test(aus, can)
t.test(aus, us)
t.test(us, can)

#### Descriptive Stats Graphs ####

df <- data.frame(age = conjoint1$age[1:196],gender = as.factor(as.character(conjoint1$gender[1:196])), conjoint1$ideology[1:196])

#Age
ggplot(df, aes(x=age)) +
  geom_density(fill="red") +
  labs(y="Density", x = "Age")

ggsave(paste0("age.png" ,""), width = 15, height = 10, units = c("cm"), dpi = 300)

#Gender
ggplot(df, aes(x=gender)) +
  geom_bar(fill="red",aes(y=(..count..)*100/sum(..count..))) +
  ylim(0,100) +
  labs(y="Percent", x="Gender")

ggsave(paste0("gender.png" ,""), width = 15, height = 10, units = c("cm"), dpi = 300)

#Ideology
ggplot(df, aes(x=ideology)) +
  geom_density(fill="red") +
  labs(y="Density", x = "Ideology") +
  xlim(0,10)

ggsave(paste0("ideology.png" ,""), width = 15, height = 10, units = c("cm"), dpi = 300)