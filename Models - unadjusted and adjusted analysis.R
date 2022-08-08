###########################################################
#' GEE and logistic unadjusted and adjusted models for: 
#'      "Clinical and radiologic factors associated with 
#'      detection of Mycobacterium tuberculosis in 
#'      children under 5 years old using invasive and 
#'      noninvasive sample collection techniques - Kenya"
#'
#' Author: Jonathan Smith, PhD, MPH
#'         Global Tuberculosis Branch
#'         Division of Global HIV and Tuberculosis
#'         Centers for Disease Control and Prevention
###########################################################
# rm(list = ls())

# load libraries
library(gee)
library(dplyr)

#' READ IN DATA
#' Change filepath to where data and code are 
#' saved from github. Default to working directory
flpth <- getwd()
geedata <- read.csv(paste0(flpth, "/TOTO_gee_model_data.csv"))

# Convert all covariates into R factors
idx <- 1:ncol(geedata)
geedata[idx] <- lapply(geedata[idx], as.factor)


########################################################
## Global setup for both unadjusted and adjusted models
########################################################

#' Set `all` to:
#'  TRUE for all symptomatic children (primary analysis) 
#'  FALSE for only those with a TB diagnosis (secondary analysis)
all <- TRUE

#' Function to obtain Wald confidence intervals for GEE analysis
#' @param model_object GEE model
#' @param ci Desired confidence interval level (default 95% CI)
#' 
gee_summary <- function(model_object, ci = 0.95) {
  a <- coef(summary(model_object))
  mult <- qnorm((1 + ci) / 2)
  restab <- with(
    as.data.frame(a),
    cbind(
      est = Estimate,
      lwr =  Estimate - mult * `Robust S.E.`,
      upr = Estimate + mult * `Robust S.E.`
    )
  )
  rownames(restab) <- rownames(a)
  return(data.frame(restab))
}


sample_types <- names(table(geedata$sample))
var_names <- c("sex","hiv_cat","agegroup","wfa","ie_final","he_final", "cough_di", "fever_di", "lethargy_di",
               "sympcount","cxr","qft_result","tst_result")

###################################################
## Unadjusted models
###################################################
gee_unadjusted <- list()
for(i in sample_types){
  for(j in var_names){
    if (i != "inspt"){
      if(all){
        gee_unadjusted[[i]][[j]] <- round(exp(gee_summary(gee(result ~ get(j), data = geedata[geedata$sample == i,],
                                                              id = id, family = binomial, corstr = "exchangeable"))[-1,]),2)
      } else {
        gee_unadjusted[[i]][[j]] <- round(exp(gee_summary(gee(result ~ get(j), data = geedata[geedata$sample == i & geedata$casedef %in% c(1:2),],
                                                              id = id, family = binomial, corstr = "exchangeable"))[-1,]),2)
      }
    }
    else {
      if (j == "hiv_cat"){
        gee_unadjusted[["inspt"]][[j]] <- rbind(cbind(est = NA, lwr = NA, upr = NA), cbind(est = NA, lwr = NA, upr = NA)) 
      }
      else {
        if(all){
          gee_unadjusted[["inspt"]][[j]] <- round(exp(gee_summary(gee(result ~ get(j), data = geedata[geedata$sample == i,],
                                                                      id = id, family = binomial, corstr = "exchangeable"))[-1,]),2)
        } else {
          gee_unadjusted[["inspt"]][[j]] <- round(exp(gee_summary(gee(result ~ get(j), data = geedata[geedata$sample == i & geedata$casedef %in% c(1:2),],
                                                                      id = id, family = binomial, corstr = "exchangeable"))[-1,]),2)
        }}
    }}
}

# Arrange results
unadj_gee_list <- lapply(gee_unadjusted, function(x) do.call(rbind, x))
unadj_gee_table <- do.call(data.frame, unadj_gee_list)

#' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#'  Standard Logistic regression for any positive specimen
#' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
anysamp <- geedata %>% group_by(id) %>% mutate(anypos = any(result == 1)*1)

if(all){
  anysamp <- unique(anysamp[,-c(2:5)])
} else {
  anysamp <- unique(anysamp[anysamp$casedef %in% 1:2,-c(2:5)])
}

logit_unadjusted <- list()
for (j in var_names){
  logit_unadjusted[[j]] <- cbind(exp(coef(glm(anypos ~ get(j), data = anysamp, family = "binomial")))[-1],
                                 exp(confint(glm(anypos ~ get(j), data = anysamp, family = "binomial")))[-1,1],
                                 exp(confint(glm(anypos ~ get(j), data = anysamp, family = "binomial")))[-1,2])
  rownames(logit_unadjusted[[j]]) <- sprintf(paste0(j,"%s"),seq(1:nrow(logit_unadjusted[[j]])))
}
unadj_logit_table <- do.call(rbind, logit_unadjusted)
rownames(unadj_logit_table) <- rownames(unadj_gee_table)
colnames(unadj_logit_table) <- c("any.est","any.lwr","any.upr")


##################################################
## Adjusted models
##################################################

#' Specify best multivariable model determined by quantitative and 
#' qualitative information (see supplemental materials)

adj_models <- list(sex = as.formula("result ~ sex + ie_final + he_final + lethargy_di + cxr"),
                   hiv_cat = as.formula("result ~ hiv_cat"),
                   agegroup = as.formula("result ~ agegroup"),
                   wfa = as.formula("result ~ wfa"),
                   ie_final = as.formula("result ~ ie_final + he_final + lethargy_di + cxr"),
                   he_final = as.formula("result ~ he_final + ie_final"),
                   cough_di = as.formula("result ~ cough_di"),
                   fever_di = as.formula("result ~ fever_di"),
                   lethargy_di = as.formula("result ~ lethargy_di"),
                   sympcount = as.formula("result ~ sympcount"),
                   cxr = as.formula("result ~ cxr + ie_final"),
                   qft_result = as.formula("result ~ qft_result + he_final + lethargy_di + cxr"),
                   tst_result = as.formula("result ~ tst_result + he_final + lethargy_di + cxr"))

gee_adjusted <- list()
for(i in sample_types){
  for(j in 1:length(adj_models)){
    if( i != "inspt"){
      if(all){
        gee_adjusted[[i]][[var_names[j]]] <- round(exp(gee_summary(gee(adj_models[[j]], data = geedata[geedata$sample %in% i,],
                                                                       id = id, family = binomial, corstr = "exchangeable"))[-1,]),2)
      } else {
        gee_adjusted[[i]][[var_names[j]]] <- round(exp(gee_summary(gee(adj_models[[j]], data = geedata[geedata$sample %in% i & geedata$casedef %in% c(1:2),],
                                                                       id = id, family = binomial, corstr = "exchangeable"))[-1,]),2)
      }
    } else {
      if (j == 2){
        gee_adjusted[["inspt"]][[var_names[j]]] <- rbind(cbind(est = NA, lwr = NA, upr = NA),
                                                         cbind(est = NA, lwr = NA, upr = NA)) 
      } else {
        if(all){
          gee_adjusted[[i]][[var_names[j]]] <- round(exp(gee_summary(gee(adj_models[[j]], data = geedata[geedata$sample %in% i,],
                                                                         id = id, family = binomial, corstr = "exchangeable"))[-1,]),2)
        } else {
          gee_adjusted[[i]][[var_names[j]]] <- round(exp(gee_summary(gee(adj_models[[j]], data = geedata[geedata$sample %in% i & geedata$casedef %in% c(1:2),],
                                                                         id = id, family = binomial, corstr = "exchangeable"))[-1,]),2)
        }
      }
    } 
  }
}

adj_gee_list <- lapply(gee_adjusted, function(x) do.call(rbind, x))
adj_gee_table <- do.call(data.frame, adj_gee_list)

adj_gee_table <- adj_gee_table[-c(2:5, 12:14, 16, 23, 25:27, 29:31),]

# Classical models
adj_models_logit <- list(sex = as.formula("anypos ~ sex + ie_final + he_final + lethargy_di + cxr"),
                         hiv_status = as.formula("anypos ~ hiv_cat"),
                         agegroup = as.formula("anypos ~ agegroup"),
                         wfa = as.formula("anypos ~ wfa"),
                         ie_final = as.formula("anypos ~ ie_final + he_final + lethargy_di + cxr"),
                         he_final = as.formula("anypos ~ he_final + ie_final"),
                         cough_di = as.formula("anypos ~ cough_di"),
                         fever_di = as.formula("anypos ~ fever_di"),
                         lethargy_di = as.formula("anypos ~ lethargy_di"),
                         sympcount = as.formula("anypos ~ sympcount"),
                         cxr = as.formula("anypos ~ cxr + hiv_cat + ie_final"),
                         qft_result = as.formula("anypos ~ qft_result + he_final + lethargy_di + cxr"),
                         tst_result = as.formula("anypos ~ tst_result + he_final + lethargy_di + cxr"))

logit_adjusted <- list()
for (j in 1:length(adj_models_logit)){
  logit_adjusted[[var_names[j]]] <- cbind(exp(coef(glm(adj_models_logit[[j]], data = anysamp, family = "binomial")))[-1],
                                          exp(confint(glm(adj_models_logit[[j]], data = anysamp, family = "binomial")))[-1,1],
                                          exp(confint(glm(adj_models_logit[[j]], data = anysamp, family = "binomial")))[-1,2])
  rownames(logit_adjusted[[var_names[j]]]) <- sprintf(paste0(var_names[j],"%s"),seq(1:nrow(logit_adjusted[[var_names[j]]])))
}
adj_logit_table <- do.call(rbind, logit_adjusted)
adj_logit_table <- round(adj_logit_table[c("sex1","hiv_cat1","hiv_cat2", "agegroup1", "agegroup2",
                                           "wfa1","ie_final1","he_final1","cough_di1","fever_di1",
                                           "lethargy_di1", "sympcount1","sympcount2","cxr1",
                                           "qft_result1","tst_result1"),],2)

#rownames(adj_logit_table) <- rownames(adj_gee_table)
colnames(adj_logit_table) <- c("any.est","any.lwr","any.upr")
