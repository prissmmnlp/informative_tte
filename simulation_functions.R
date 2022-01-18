library(survival)
library(survminer)
library(dplyr)
library(SurvCorr)
library(tranSurv)
library(foreach)
library(doParallel)
library(survRM2)
library(tibble)
library(readr)
library(caret)



logit <- function(x){log(x/(1-x))}

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}




create_patient_mort_hazard_vector <- function(max_time, std_dev, already_genotyped, genotyping_mort_multiplier, pop_mort_baseline_hazard, pop_mort_gamma, biomarker, biomarker_mult, biomarker_gamma, random_delta){
  
  updowns = ((c(1:max_time)/max_time) - 0.5)*pop_mort_gamma
  logit_hazard_vector = rnorm(1,mean=logit(pop_mort_baseline_hazard) + already_genotyped*genotyping_mort_multiplier, sd=std_dev) + cumsum(ifelse(rbinom(max_time, 1, 0.5 + updowns)==0, -1, 1) * random_delta)
  # port mort gamma goes from -1 (decreasing hazard) to +1 (increasing hazard)
  
  output_logits = logit_hazard_vector + biomarker*biomarker_mult + biomarker*biomarker_mult*biomarker_gamma
  output_probs = sigmoid(output_logits)
  return(output_probs)
  
}


create_progression_hazard_vector <- function(themean, std_dev, mort_hazard_vector){
  return(sigmoid(logit(mort_hazard_vector) + rnorm(1,mean=themean,sd=std_dev)))
}





create_patient_testing_hazard_vector <- function(max_time, already_genotyped, multiplier, increasing_multiplier, patient_mort_hazard_vector, pop_genomic_baseline_hazard, 
                                                 pop_genomic_gamma, genomic_sd, random_genomic_delta){
  
  if(already_genotyped==1){
    return(rep(1.0, max_time))
  }
  # derive independent testing prob
  updowns = ((c(1:max_time)/max_time) - 0.5)*pop_genomic_gamma
  logit_hazard_vector = rnorm(1,mean=logit(pop_genomic_baseline_hazard), sd=genomic_sd) + cumsum(ifelse(rbinom(max_time, 1, 0.5 + updowns)==0, -1, 1) * random_genomic_delta)
  
  
  if(increasing_multiplier==0){
    output_vector = sigmoid(logit_hazard_vector + multiplier*logit(patient_mort_hazard_vector))
  } else{
    times = c(1:max_time)
    # make the relationship between testing risk and mortality risk vary with time (minimal at first, sqrt(time) thereafter)
    output_vector = sigmoid(logit_hazard_vector + multiplier*sqrt(times)*logit(patient_mort_hazard_vector))
  }
  
  
  return(output_vector)
}


create_patient_right_censoring_hazard_vector <- function(right_censor_hazard = 0.05, max_time){
  hazard = right_censor_hazard
  return(rep(hazard, max_time))
}

create_patient_observation_vector <- function(observed_proba = 0.1, mort_hazard_multiplier=0, mort_hazard_vector){
  observed_logit = logit(observed_proba) + mort_hazard_multiplier*logit(mort_hazard_vector)
  observed_prob = sigmoid(observed_logit)
  # always observe last month for now
  observed_prob[length(mort_hazard_vector)] = 1
  return(observed_prob)
}


simulate_single_cohort <- function(num_patients=1000, len_follow_up=60, biomarker_effect = 2, biomarker_prevalence = 0.3,
                                   biomarker_gamma = 0.0,
                                   pre_genotyping_proportion = 0.0,
                                   pop_mort_baseline_hazard = 0.02,
                                   pop_mort_gamma = 0,
                                   pop_genomic_baseline_hazard = 0.25,
                                   pop_genomic_gamma = 0,
                                   genomic_sd = 1,
                                   random_genomic_delta = 0.5,
                                   patient_mort_hazard_std_dev = 2, patient_mort_hazard_random_delta = 0.5,
                                   patient_testing_hazard_multiplier = 2, 
                                   patient_testing_hazard_multiplier_increasing = 0,
                                   progression_hazard_logit_mean = 1,
                                   progression_hazard_std_dev = 1,
                                   right_censor_hazard = 0.05,
                                   informed_right_censoring = 0,
                                   observation_prob = 0.1,
                                   informative_observation_multiplier = 0){
  
  
  already_genotyped = rbinom(num_patients, 1, pre_genotyping_proportion)
  

  biomarker = rbinom(num_patients, 1, biomarker_prevalence)
  biomarker_mult = biomarker_effect
  
  cohort = vector("list", num_patients)
  for(pt in 1:num_patients) {
    results = list(10)
    
    # risks
    # mortality hazard
    results[[1]] = create_patient_mort_hazard_vector(len_follow_up, already_genotyped[pt], patient_testing_hazard_multiplier, std_dev = patient_mort_hazard_std_dev, pop_mort_baseline_hazard, pop_mort_gamma, biomarker[pt], biomarker_mult, biomarker_gamma, random_delta = patient_mort_hazard_random_delta)
    
    # testing hazard
    results[[2]] = create_patient_testing_hazard_vector(len_follow_up, already_genotyped[pt], multiplier = patient_testing_hazard_multiplier, patient_testing_hazard_multiplier_increasing, results[[1]], pop_genomic_baseline_hazard, pop_genomic_gamma, genomic_sd, random_genomic_delta)
    
    # progression hazard
    results[[3]] = create_progression_hazard_vector(progression_hazard_logit_mean, progression_hazard_std_dev, results[[1]])
    
    # right censoring hazard
    results[[4]] = create_patient_right_censoring_hazard_vector(right_censor_hazard, len_follow_up)
    
    # observation hazard
    results[[5]] = create_patient_observation_vector(observation_prob, informative_observation_multiplier, results[[1]])
    
    # events
    
    for(i in 1:5){
      results[[i+5]] = rbinom(len_follow_up, 1, results[[i]])
      #results[[i+5]] = vapply(results[[i]], function(x){rbinom(1,1,x)}, FUN.VALUE = numeric(1))
      
    }
    
    # observation times
    results[[11]] = which(results[[10]] == 1)
    
    
    # results[[12]] survival times
    
    # results[[13]] death indicator
    # 14 testing time
    # 15 testing indicator
    # 16 progression times
    # 17 progression indicator
    # 18 right censor times
    # 19 right censor indicator
    
    j = 1
    for(i in c(1,3,5,7)){
      
      results[[i+11]] = ifelse(sum(results[[j+5]]) > 0, which.max(results[[j+5]]), NA)
      results[[i+12]] = ifelse(is.na(results[[i+11]]), 0, 1)
      results[[i+11]] = ifelse(results[[i+12]] == 0, len_follow_up, results[[i+11]])
      j = j + 1
      
    }
    
    #20 mortality hazard at time of testing
    results[[20]] = results[[1]][results[[14]]]
    
    cohort[[pt]] = results
  }
  
  # get vector of survival times across cohort
  survival_times <- vapply(cohort, function(x) x[[12]], FUN.VALUE = numeric(1))
  died <- vapply(cohort, function(x) x[[13]], FUN.VALUE = numeric(1))
  testing_times <- vapply(cohort, function(x) x[[14]], FUN.VALUE = numeric(1))
  tested <- vapply(cohort, function(x) x[[15]], FUN.VALUE = numeric(1))
  
  # force testing times to zero for patients previously tested 
  testing_times <- ifelse(already_genotyped == 1, 0, testing_times)
  tested <- ifelse(already_genotyped == 1, 1, tested)
  
  progression_times <- vapply(cohort, function(x) x[[16]], FUN.VALUE = numeric(1))
  progressed <- vapply(cohort, function(x) x[[17]], FUN.VALUE = numeric(1))
  right_censor_times <- vapply(cohort, function(x) x[[18]], FUN.VALUE = numeric(1))
  right_censored <- vapply(cohort, function(x) x[[19]], FUN.VALUE = numeric(1))
  mortality_risk_at_testing <- vapply(cohort, function(x) x[[20]], FUN.VALUE = numeric(1))
  
  patient_obs_times <- lapply(cohort, function(x) x[[11]])
  
  
  
  true_cohort = data.frame(cbind(biomarker, survival_times, died, testing_times, tested, progression_times, progressed, right_censor_times, right_censored, mortality_risk_at_testing)) %>% mutate(pfs_times = pmin(progression_times, survival_times)) %>% mutate(pfs_event = ifelse(died==0 & pfs_times == len_follow_up, 0, 1))
  
  
  
  true_cohort <- true_cohort %>% mutate(observed = ifelse(tested == 1 & (survival_times >= testing_times) & (right_censor_times >= testing_times), TRUE, FALSE))
  observed_obs_times <- patient_obs_times[true_cohort$observed]
  observed <- true_cohort %>% filter(observed==TRUE)
  
  # deal with informative right censoring
  if(informed_right_censoring==0){
    observed <- observed %>% mutate(progression_times = ifelse(right_censored == 1 & right_censor_times < progression_times, right_censor_times, progression_times)) %>% mutate(progressed = ifelse(right_censored == 1 & progression_times == right_censor_times, 0, progressed)) %>% mutate(survival_times = ifelse(right_censored == 1 & right_censor_times < survival_times, right_censor_times, survival_times)) %>% mutate(died = ifelse(right_censored == 1 & survival_times == right_censor_times, 0, died))
  } else {
    observed <- observed %>% mutate(progression_times = ifelse(right_censored == 1 & right_censor_times < progression_times, right_censor_times, progression_times)) %>% mutate(progressed = ifelse(right_censored == 1 & progression_times == right_censor_times, 0, progressed)) %>% mutate(survival_times = ifelse(died==0 & right_censored == 1 & right_censor_times < survival_times, right_censor_times, survival_times))
    
  }
  
  
  # adjust progression time based on observation time
  for(i in which(observed$progressed ==1)){
    obs_times = observed_obs_times[[i]]
    observed[i,]$progression_times = obs_times[obs_times >= observed[i,]$progression_times][1]
  }
  
  
  
  observed <- observed %>% mutate(pfs_times = pmin(survival_times, progression_times)) %>% mutate(pfs_event = ifelse(progressed==1 | died == 1, 1, 0))
  
  
  return(list(true_cohort=true_cohort, observed_cohort=observed))
  
}



eval_differences_naive <- function(true_cohort, observed_cohort){
  true_os_biomarker_results <- tryCatch({
    summary(coxph(Surv(survival_times, died) ~ biomarker, data = true_cohort))$coefficients
  }, error=function(cond) {return(NA)})
  
  true_os_biomarker_coef = true_os_biomarker_results[1]
  true_os_biomarker_lower = true_os_biomarker_results[1] - 1.96 * true_os_biomarker_results[3]
  true_os_biomarker_upper = true_os_biomarker_results[1] + 1.96 * true_os_biomarker_results[3]
  os_biomarker_difference_exists = ifelse(true_os_biomarker_results[5] < 0.05, 1, 0)
  
  
  observed_os_biomarker_results <- tryCatch({
    summary(coxph(Surv(survival_times, died) ~ biomarker, data = observed_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  observed_os_biomarker_coef = observed_os_biomarker_results[1]
  observed_os_biomarker_lower = observed_os_biomarker_results[1] - 1.96 * observed_os_biomarker_results[3]
  observed_os_biomarker_upper = observed_os_biomarker_results[1] + 1.96 * observed_os_biomarker_results[3]
  
  os_biomarker_difference_detected = ifelse(observed_os_biomarker_results[5] < 0.05, 1, 0)
  os_biomarker_type1_error = ifelse(os_biomarker_difference_detected ==1 & os_biomarker_difference_exists == 0, 1, 0)
  os_biomarker_type2_error = ifelse(os_biomarker_difference_detected ==0 & os_biomarker_difference_exists == 1, 1, 0)
  
  
  os_biomarker_coef_diff = observed_os_biomarker_coef - true_os_biomarker_coef
  
  
  
  true_pfs_biomarker_results <- tryCatch({
    summary(coxph(Surv(pfs_times, pfs_event) ~ biomarker, data = true_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  true_pfs_biomarker_coef = true_pfs_biomarker_results[1]
  true_pfs_biomarker_lower = true_pfs_biomarker_results[1] - 1.96 * true_pfs_biomarker_results[3]
  true_pfs_biomarker_upper = true_pfs_biomarker_results[1] + 1.96 * true_pfs_biomarker_results[3]
  pfs_biomarker_difference_exists = ifelse(true_pfs_biomarker_results[5] < 0.05, 1, 0)
  
  observed_pfs_biomarker_results <- tryCatch({
    summary(coxph(Surv(pfs_times, pfs_event) ~ biomarker, data = observed_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  observed_pfs_biomarker_coef = observed_pfs_biomarker_results[1]
  observed_pfs_biomarker_lower = observed_pfs_biomarker_results[1] - 1.96 * observed_pfs_biomarker_results[3]
  observed_pfs_biomarker_upper = observed_pfs_biomarker_results[1] + 1.96 * observed_pfs_biomarker_results[3]
  
  pfs_biomarker_difference_detected = ifelse(observed_pfs_biomarker_results[5] < 0.05, 1, 0)
  pfs_biomarker_type1_error = ifelse(pfs_biomarker_difference_detected ==1 & pfs_biomarker_difference_exists == 0, 1, 0)
  pfs_biomarker_type2_error = ifelse(pfs_biomarker_difference_detected ==0 & pfs_biomarker_difference_exists == 1, 1, 0)
  
  pfs_biomarker_coef_diff = observed_pfs_biomarker_coef - true_pfs_biomarker_coef
  
  
  
  prop_observed = nrow(observed_cohort) / nrow(true_cohort)
  
  os_true_median <- tryCatch({
    os_true_survfit <- survfit(Surv(survival_times, died) ~ 1, data = true_cohort)
    unname(summary(os_true_survfit)$table['median'])
  }, error = function(cond) {return(NA)})
  
  
  os_observed_median <- tryCatch({
    os_observed_survfit <- survfit(Surv(survival_times, died) ~ 1, data = observed_cohort)
    unname(summary(os_observed_survfit)$table['median'])  
  }, error = function(cond) {return(NA)})
  
  # observed minus true
  os_observed_minus_true = os_observed_median - os_true_median
  
  
  pfs_true_median <- tryCatch({
    pfs_true_survfit <- survfit(Surv(pfs_times, pfs_event) ~ 1, data = true_cohort)
    unname(summary(pfs_true_survfit)$table['median'])
  }, error = function(cond) {return(NA)})
  
  
  pfs_observed_median <- tryCatch({
    pfs_observed_survfit <- survfit(Surv(pfs_times, pfs_event) ~ 1, data = observed_cohort)
    unname(summary(pfs_observed_survfit)$table['median'])  
  }, error = function(cond) {return(NA)})
  
  # observed minus true
  os_observed_minus_true = os_observed_median - os_true_median  
  pfs_observed_minus_true = pfs_observed_median - pfs_true_median
  
  
  return(list(prop_observed=prop_observed,
              os_true_median=os_true_median,
              os_observed_median=os_observed_median,
              os_observed_minus_true=os_observed_minus_true,
              pfs_true_median=pfs_true_median,
              pfs_observed_median=pfs_observed_median,
              pfs_observed_minus_true=pfs_observed_minus_true,
              true_os_biomarker_coef=true_os_biomarker_coef,
              observed_os_biomarker_coef=observed_os_biomarker_coef,
              observed_os_biomarker_lower=observed_os_biomarker_lower,
              observed_os_biomarker_upper=observed_os_biomarker_upper,
              os_biomarker_difference_exists=os_biomarker_difference_exists,
              os_biomarker_difference_detected=os_biomarker_difference_detected,
              os_biomarker_type1_error=os_biomarker_type1_error,
              os_biomarker_type2_error=os_biomarker_type2_error,
              os_biomarker_coef_diff=os_biomarker_coef_diff,
              true_pfs_biomarker_coef=true_pfs_biomarker_coef,
              observed_pfs_biomarker_coef=observed_pfs_biomarker_coef,
              observed_pfs_biomarker_lower=observed_pfs_biomarker_lower,
              observed_pfs_biomarker_upper=observed_pfs_biomarker_upper,
              pfs_biomarker_difference_exists=pfs_biomarker_difference_exists,
              pfs_biomarker_difference_detected=pfs_biomarker_difference_detected,
              pfs_biomarker_type1_error=pfs_biomarker_type1_error,
              pfs_biomarker_type2_error=pfs_biomarker_type2_error,
              pfs_biomarker_coef_diff=pfs_biomarker_coef_diff))
}




eval_differences_left_trunc <- function(true_cohort, observed_cohort){
  
  
  
  true_os_biomarker_results <- tryCatch({
    summary(coxph(Surv(survival_times, died) ~ biomarker, data = true_cohort))$coefficients
  }, error=function(cond) {return(NA)})
  
  true_os_biomarker_coef = true_os_biomarker_results[1]
  true_os_biomarker_lower = true_os_biomarker_results[1] - 1.96 * true_os_biomarker_results[3]
  true_os_biomarker_upper = true_os_biomarker_results[1] + 1.96 * true_os_biomarker_results[3]
  os_biomarker_difference_exists = ifelse(true_os_biomarker_results[5] < 0.05, 1, 0)
  
  
  observed_os_biomarker_results <- tryCatch({
    summary(coxph(Surv(testing_times, survival_times, died) ~ biomarker, data = observed_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  observed_os_biomarker_coef = observed_os_biomarker_results[1]
  observed_os_biomarker_lower = observed_os_biomarker_results[1] - 1.96 * observed_os_biomarker_results[3]
  observed_os_biomarker_upper = observed_os_biomarker_results[1] + 1.96 * observed_os_biomarker_results[3]
  
  os_biomarker_difference_detected = ifelse(observed_os_biomarker_results[5] < 0.05, 1, 0)
  os_biomarker_type1_error = ifelse(os_biomarker_difference_detected ==1 & os_biomarker_difference_exists == 0, 1, 0)
  os_biomarker_type2_error = ifelse(os_biomarker_difference_detected ==0 & os_biomarker_difference_exists == 1, 1, 0)
  
  
  os_biomarker_coef_diff = observed_os_biomarker_coef - true_os_biomarker_coef
  
  
  
  true_pfs_biomarker_results <- tryCatch({
    summary(coxph(Surv(pfs_times, pfs_event) ~ biomarker, data = true_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  true_pfs_biomarker_coef = true_pfs_biomarker_results[1]
  true_pfs_biomarker_lower = true_pfs_biomarker_results[1] - 1.96 * true_pfs_biomarker_results[3]
  true_pfs_biomarker_upper = true_pfs_biomarker_results[1] + 1.96 * true_pfs_biomarker_results[3]
  pfs_biomarker_difference_exists = ifelse(true_pfs_biomarker_results[5] < 0.05, 1, 0)
  
  observed_pfs_biomarker_results <- tryCatch({
    summary(coxph(Surv(testing_times, pfs_times, pfs_event) ~ biomarker, data = observed_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  observed_pfs_biomarker_coef = observed_pfs_biomarker_results[1]
  observed_pfs_biomarker_lower = observed_pfs_biomarker_results[1] - 1.96 * observed_pfs_biomarker_results[3]
  observed_pfs_biomarker_upper = observed_pfs_biomarker_results[1] + 1.96 * observed_pfs_biomarker_results[3]
  
  pfs_biomarker_difference_detected = ifelse(observed_pfs_biomarker_results[5] < 0.05, 1, 0)
  pfs_biomarker_type1_error = ifelse(pfs_biomarker_difference_detected ==1 & pfs_biomarker_difference_exists == 0, 1, 0)
  pfs_biomarker_type2_error = ifelse(pfs_biomarker_difference_detected ==0 & pfs_biomarker_difference_exists == 1, 1, 0)
  
  pfs_biomarker_coef_diff = observed_pfs_biomarker_coef - true_pfs_biomarker_coef
  
  
  
  prop_observed = nrow(observed_cohort) / nrow(true_cohort)
  
  os_true_median <- tryCatch({
    os_true_survfit <- survfit(Surv(survival_times, died) ~ 1, data = true_cohort)
    unname(summary(os_true_survfit)$table['median'])
  }, error = function(cond) {return(NA)})
  
  
  os_observed_median <- tryCatch({
    os_observed_survfit <- survfit(Surv(testing_times, survival_times, died) ~ 1, data = observed_cohort)
    unname(summary(os_observed_survfit)$table['median'])  
  }, error = function(cond) {return(NA)})
  
  # observed minus true
  os_observed_minus_true = os_observed_median - os_true_median
  
  
  pfs_true_median <- tryCatch({
    pfs_true_survfit <- survfit(Surv(pfs_times, pfs_event) ~ 1, data = true_cohort)
    unname(summary(pfs_true_survfit)$table['median'])
  }, error = function(cond) {return(NA)})
  
  
  pfs_observed_median <- tryCatch({
    pfs_observed_survfit <- survfit(Surv(testing_times, pfs_times, pfs_event) ~ 1, data = observed_cohort)
    unname(summary(pfs_observed_survfit)$table['median'])  
  }, error = function(cond) {return(NA)})
  
  # observed minus true
  os_observed_minus_true = os_observed_median - os_true_median  
  pfs_observed_minus_true = pfs_observed_median - pfs_true_median
  
  
  return(list(prop_observed=prop_observed,
              os_true_median=os_true_median,
              os_observed_median=os_observed_median,
              os_observed_minus_true=os_observed_minus_true,
              pfs_true_median=pfs_true_median,
              pfs_observed_median=pfs_observed_median,
              pfs_observed_minus_true=pfs_observed_minus_true,
              true_os_biomarker_coef=true_os_biomarker_coef,
              observed_os_biomarker_coef=observed_os_biomarker_coef,
              observed_os_biomarker_lower=observed_os_biomarker_lower,
              observed_os_biomarker_upper=observed_os_biomarker_upper,
              os_biomarker_difference_exists=os_biomarker_difference_exists,
              os_biomarker_difference_detected=os_biomarker_difference_detected,
              os_biomarker_type1_error=os_biomarker_type1_error,
              os_biomarker_type2_error=os_biomarker_type2_error,
              os_biomarker_coef_diff=os_biomarker_coef_diff,
              true_pfs_biomarker_coef=true_pfs_biomarker_coef,
              observed_pfs_biomarker_coef=observed_pfs_biomarker_coef,
              observed_pfs_biomarker_lower=observed_pfs_biomarker_lower,
              observed_pfs_biomarker_upper=observed_pfs_biomarker_upper,
              pfs_biomarker_difference_exists=pfs_biomarker_difference_exists,
              pfs_biomarker_difference_detected=pfs_biomarker_difference_detected,
              pfs_biomarker_type1_error=pfs_biomarker_type1_error,
              pfs_biomarker_type2_error=pfs_biomarker_type2_error,
              pfs_biomarker_coef_diff=pfs_biomarker_coef_diff))
  
  
}


eval_differences_left_trunc_timetotest_covariate <- function(true_cohort, observed_cohort){
  
  
  
  true_os_biomarker_results <- tryCatch({
    summary(coxph(Surv(survival_times, died) ~ biomarker, data = true_cohort))$coefficients
  }, error=function(cond) {return(NA)})
  
  true_os_biomarker_coef = true_os_biomarker_results[1]
  true_os_biomarker_lower = true_os_biomarker_results[1] - 1.96 * true_os_biomarker_results[3]
  true_os_biomarker_upper = true_os_biomarker_results[1] + 1.96 * true_os_biomarker_results[3]
  os_biomarker_difference_exists = ifelse(true_os_biomarker_results[5] < 0.05, 1, 0)
  
  
  observed_os_biomarker_results <- tryCatch({
    summary(coxph(Surv(testing_times, survival_times, died) ~ biomarker + testing_times, data = observed_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  observed_os_biomarker_coef = observed_os_biomarker_results[1]
  observed_os_biomarker_lower = observed_os_biomarker_results[1] - 1.96 * observed_os_biomarker_results[3]
  observed_os_biomarker_upper = observed_os_biomarker_results[1] + 1.96 * observed_os_biomarker_results[3]
  
  os_biomarker_difference_detected = ifelse(observed_os_biomarker_results[5] < 0.05, 1, 0)
  os_biomarker_type1_error = ifelse(os_biomarker_difference_detected ==1 & os_biomarker_difference_exists == 0, 1, 0)
  os_biomarker_type2_error = ifelse(os_biomarker_difference_detected ==0 & os_biomarker_difference_exists == 1, 1, 0)
  
  
  os_biomarker_coef_diff = observed_os_biomarker_coef - true_os_biomarker_coef
  
  
  
  true_pfs_biomarker_results <- tryCatch({
    summary(coxph(Surv(pfs_times, pfs_event) ~ biomarker, data = true_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  true_pfs_biomarker_coef = true_pfs_biomarker_results[1]
  true_pfs_biomarker_lower = true_pfs_biomarker_results[1] - 1.96 * true_pfs_biomarker_results[3]
  true_pfs_biomarker_upper = true_pfs_biomarker_results[1] + 1.96 * true_pfs_biomarker_results[3]
  pfs_biomarker_difference_exists = ifelse(true_pfs_biomarker_results[5] < 0.05, 1, 0)
  
  observed_pfs_biomarker_results <- tryCatch({
    summary(coxph(Surv(testing_times, pfs_times, pfs_event) ~ biomarker + testing_times, data = observed_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  observed_pfs_biomarker_coef = observed_pfs_biomarker_results[1]
  observed_pfs_biomarker_lower = observed_pfs_biomarker_results[1] - 1.96 * observed_pfs_biomarker_results[3]
  observed_pfs_biomarker_upper = observed_pfs_biomarker_results[1] + 1.96 * observed_pfs_biomarker_results[3]
  
  pfs_biomarker_difference_detected = ifelse(observed_pfs_biomarker_results[5] < 0.05, 1, 0)
  pfs_biomarker_type1_error = ifelse(pfs_biomarker_difference_detected ==1 & pfs_biomarker_difference_exists == 0, 1, 0)
  pfs_biomarker_type2_error = ifelse(pfs_biomarker_difference_detected ==0 & pfs_biomarker_difference_exists == 1, 1, 0)
  
  pfs_biomarker_coef_diff = observed_pfs_biomarker_coef - true_pfs_biomarker_coef
  
  
  
  prop_observed = nrow(observed_cohort) / nrow(true_cohort)
  
  os_true_median <- tryCatch({
    os_true_survfit <- survfit(Surv(survival_times, died) ~ 1, data = true_cohort)
    unname(summary(os_true_survfit)$table['median'])
  }, error = function(cond) {return(NA)})
  
  
  os_observed_median <- tryCatch({
    os_observed_survfit <- survfit(Surv(testing_times, survival_times, died) ~ 1, data = observed_cohort)
    unname(summary(os_observed_survfit)$table['median'])  
  }, error = function(cond) {return(NA)})
  
  # observed minus true
  os_observed_minus_true = os_observed_median - os_true_median
  
  
  pfs_true_median <- tryCatch({
    pfs_true_survfit <- survfit(Surv(pfs_times, pfs_event) ~ 1, data = true_cohort)
    unname(summary(pfs_true_survfit)$table['median'])
  }, error = function(cond) {return(NA)})
  
  
  pfs_observed_median <- tryCatch({
    pfs_observed_survfit <- survfit(Surv(testing_times, pfs_times, pfs_event) ~ 1, data = observed_cohort)
    unname(summary(pfs_observed_survfit)$table['median'])  
  }, error = function(cond) {return(NA)})
  
  # observed minus true
  os_observed_minus_true = os_observed_median - os_true_median  
  pfs_observed_minus_true = pfs_observed_median - pfs_true_median
  
  
  return(list(prop_observed=prop_observed,
              os_true_median=os_true_median,
              os_observed_median=os_observed_median,
              os_observed_minus_true=os_observed_minus_true,
              pfs_true_median=pfs_true_median,
              pfs_observed_median=pfs_observed_median,
              pfs_observed_minus_true=pfs_observed_minus_true,
              true_os_biomarker_coef=true_os_biomarker_coef,
              observed_os_biomarker_coef=observed_os_biomarker_coef,
              observed_os_biomarker_lower=observed_os_biomarker_lower,
              observed_os_biomarker_upper=observed_os_biomarker_upper,
              os_biomarker_difference_exists=os_biomarker_difference_exists,
              os_biomarker_difference_detected=os_biomarker_difference_detected,
              os_biomarker_type1_error=os_biomarker_type1_error,
              os_biomarker_type2_error=os_biomarker_type2_error,
              os_biomarker_coef_diff=os_biomarker_coef_diff,
              true_pfs_biomarker_coef=true_pfs_biomarker_coef,
              observed_pfs_biomarker_coef=observed_pfs_biomarker_coef,
              observed_pfs_biomarker_lower=observed_pfs_biomarker_lower,
              observed_pfs_biomarker_upper=observed_pfs_biomarker_upper,
              pfs_biomarker_difference_exists=pfs_biomarker_difference_exists,
              pfs_biomarker_difference_detected=pfs_biomarker_difference_detected,
              pfs_biomarker_type1_error=pfs_biomarker_type1_error,
              pfs_biomarker_type2_error=pfs_biomarker_type2_error,
              pfs_biomarker_coef_diff=pfs_biomarker_coef_diff))
  
  
}


eval_differences_left_trunc_timetotest_andrisk_covariate <- function(true_cohort, observed_cohort){
  
  
  
  true_os_biomarker_results <- tryCatch({
    summary(coxph(Surv(survival_times, died) ~ biomarker, data = true_cohort))$coefficients
  }, error=function(cond) {return(NA)})
  
  true_os_biomarker_coef = true_os_biomarker_results[1]
  true_os_biomarker_lower = true_os_biomarker_results[1] - 1.96 * true_os_biomarker_results[3]
  true_os_biomarker_upper = true_os_biomarker_results[1] + 1.96 * true_os_biomarker_results[3]
  os_biomarker_difference_exists = ifelse(true_os_biomarker_results[5] < 0.05, 1, 0)
  
  
  observed_os_biomarker_results <- tryCatch({
    summary(coxph(Surv(testing_times, survival_times, died) ~ biomarker + testing_times + mortality_risk_at_testing, data = observed_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  observed_os_biomarker_coef = observed_os_biomarker_results[1]
  observed_os_biomarker_lower = observed_os_biomarker_results[1] - 1.96 * observed_os_biomarker_results[3]
  observed_os_biomarker_upper = observed_os_biomarker_results[1] + 1.96 * observed_os_biomarker_results[3]
  
  os_biomarker_difference_detected = ifelse(observed_os_biomarker_results[5] < 0.05, 1, 0)
  os_biomarker_type1_error = ifelse(os_biomarker_difference_detected ==1 & os_biomarker_difference_exists == 0, 1, 0)
  os_biomarker_type2_error = ifelse(os_biomarker_difference_detected ==0 & os_biomarker_difference_exists == 1, 1, 0)
  
  
  os_biomarker_coef_diff = observed_os_biomarker_coef - true_os_biomarker_coef
  
  
  
  true_pfs_biomarker_results <- tryCatch({
    summary(coxph(Surv(pfs_times, pfs_event) ~ biomarker, data = true_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  true_pfs_biomarker_coef = true_pfs_biomarker_results[1]
  true_pfs_biomarker_lower = true_pfs_biomarker_results[1] - 1.96 * true_pfs_biomarker_results[3]
  true_pfs_biomarker_upper = true_pfs_biomarker_results[1] + 1.96 * true_pfs_biomarker_results[3]
  pfs_biomarker_difference_exists = ifelse(true_pfs_biomarker_results[5] < 0.05, 1, 0)
  
  observed_pfs_biomarker_results <- tryCatch({
    summary(coxph(Surv(testing_times, pfs_times, pfs_event) ~ biomarker + testing_times, data = observed_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  observed_pfs_biomarker_coef = observed_pfs_biomarker_results[1]
  observed_pfs_biomarker_lower = observed_pfs_biomarker_results[1] - 1.96 * observed_pfs_biomarker_results[3]
  observed_pfs_biomarker_upper = observed_pfs_biomarker_results[1] + 1.96 * observed_pfs_biomarker_results[3]
  
  pfs_biomarker_difference_detected = ifelse(observed_pfs_biomarker_results[5] < 0.05, 1, 0)
  pfs_biomarker_type1_error = ifelse(pfs_biomarker_difference_detected ==1 & pfs_biomarker_difference_exists == 0, 1, 0)
  pfs_biomarker_type2_error = ifelse(pfs_biomarker_difference_detected ==0 & pfs_biomarker_difference_exists == 1, 1, 0)
  
  pfs_biomarker_coef_diff = observed_pfs_biomarker_coef - true_pfs_biomarker_coef
  
  
  
  prop_observed = nrow(observed_cohort) / nrow(true_cohort)
  
  os_true_median <- tryCatch({
    os_true_survfit <- survfit(Surv(survival_times, died) ~ 1, data = true_cohort)
    unname(summary(os_true_survfit)$table['median'])
  }, error = function(cond) {return(NA)})
  
  
  os_observed_median <- tryCatch({
    os_observed_survfit <- survfit(Surv(testing_times, survival_times, died) ~ 1, data = observed_cohort)
    unname(summary(os_observed_survfit)$table['median'])  
  }, error = function(cond) {return(NA)})
  
  # observed minus true
  os_observed_minus_true = os_observed_median - os_true_median
  
  
  pfs_true_median <- tryCatch({
    pfs_true_survfit <- survfit(Surv(pfs_times, pfs_event) ~ 1, data = true_cohort)
    unname(summary(pfs_true_survfit)$table['median'])
  }, error = function(cond) {return(NA)})
  
  
  pfs_observed_median <- tryCatch({
    pfs_observed_survfit <- survfit(Surv(testing_times, pfs_times, pfs_event) ~ 1, data = observed_cohort)
    unname(summary(pfs_observed_survfit)$table['median'])  
  }, error = function(cond) {return(NA)})
  
  # observed minus true
  os_observed_minus_true = os_observed_median - os_true_median  
  pfs_observed_minus_true = pfs_observed_median - pfs_true_median
  
  
  return(list(prop_observed=prop_observed,
              os_true_median=os_true_median,
              os_observed_median=os_observed_median,
              os_observed_minus_true=os_observed_minus_true,
              pfs_true_median=pfs_true_median,
              pfs_observed_median=pfs_observed_median,
              pfs_observed_minus_true=pfs_observed_minus_true,
              true_os_biomarker_coef=true_os_biomarker_coef,
              observed_os_biomarker_coef=observed_os_biomarker_coef,
              observed_os_biomarker_lower=observed_os_biomarker_lower,
              observed_os_biomarker_upper=observed_os_biomarker_upper,
              os_biomarker_difference_exists=os_biomarker_difference_exists,
              os_biomarker_difference_detected=os_biomarker_difference_detected,
              os_biomarker_type1_error=os_biomarker_type1_error,
              os_biomarker_type2_error=os_biomarker_type2_error,
              os_biomarker_coef_diff=os_biomarker_coef_diff,
              true_pfs_biomarker_coef=true_pfs_biomarker_coef,
              observed_pfs_biomarker_coef=observed_pfs_biomarker_coef,
              observed_pfs_biomarker_lower=observed_pfs_biomarker_lower,
              observed_pfs_biomarker_upper=observed_pfs_biomarker_upper,
              pfs_biomarker_difference_exists=pfs_biomarker_difference_exists,
              pfs_biomarker_difference_detected=pfs_biomarker_difference_detected,
              pfs_biomarker_type1_error=pfs_biomarker_type1_error,
              pfs_biomarker_type2_error=pfs_biomarker_type2_error,
              pfs_biomarker_coef_diff=pfs_biomarker_coef_diff))
  
  
}

eval_differences_post_testing_only <- function(true_cohort, observed_cohort){
  
  
  observed_cohort <- observed_cohort %>% filter(testing_times == 0)
  
  true_os_biomarker_results <- tryCatch({
    summary(coxph(Surv(survival_times, died) ~ biomarker, data = true_cohort))$coefficients
  }, error=function(cond) {return(NA)})
  
  true_os_biomarker_coef = true_os_biomarker_results[1]
  true_os_biomarker_lower = true_os_biomarker_results[1] - 1.96 * true_os_biomarker_results[3]
  true_os_biomarker_upper = true_os_biomarker_results[1] + 1.96 * true_os_biomarker_results[3]
  os_biomarker_difference_exists = ifelse(true_os_biomarker_results[5] < 0.05, 1, 0)
  
  
  observed_os_biomarker_results <- tryCatch({
    summary(coxph(Surv(survival_times, died) ~ biomarker, data = observed_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  observed_os_biomarker_coef = observed_os_biomarker_results[1]
  observed_os_biomarker_lower = observed_os_biomarker_results[1] - 1.96 * observed_os_biomarker_results[3]
  observed_os_biomarker_upper = observed_os_biomarker_results[1] + 1.96 * observed_os_biomarker_results[3]
  
  os_biomarker_difference_detected = ifelse(observed_os_biomarker_results[5] < 0.05, 1, 0)
  os_biomarker_type1_error = ifelse(os_biomarker_difference_detected ==1 & os_biomarker_difference_exists == 0, 1, 0)
  os_biomarker_type2_error = ifelse(os_biomarker_difference_detected ==0 & os_biomarker_difference_exists == 1, 1, 0)
  
  
  os_biomarker_coef_diff = observed_os_biomarker_coef - true_os_biomarker_coef
  
  
  
  true_pfs_biomarker_results <- tryCatch({
    summary(coxph(Surv(pfs_times, pfs_event) ~ biomarker, data = true_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  true_pfs_biomarker_coef = true_pfs_biomarker_results[1]
  true_pfs_biomarker_lower = true_pfs_biomarker_results[1] - 1.96 * true_pfs_biomarker_results[3]
  true_pfs_biomarker_upper = true_pfs_biomarker_results[1] + 1.96 * true_pfs_biomarker_results[3]
  pfs_biomarker_difference_exists = ifelse(true_pfs_biomarker_results[5] < 0.05, 1, 0)
  
  observed_pfs_biomarker_results <- tryCatch({
    summary(coxph(Surv(pfs_times, pfs_event) ~ biomarker, data = observed_cohort))$coefficients
  }, error = function(cond) {return(NA)})
  
  observed_pfs_biomarker_coef = observed_pfs_biomarker_results[1]
  observed_pfs_biomarker_lower = observed_pfs_biomarker_results[1] - 1.96 * observed_pfs_biomarker_results[3]
  observed_pfs_biomarker_upper = observed_pfs_biomarker_results[1] + 1.96 * observed_pfs_biomarker_results[3]
  
  pfs_biomarker_difference_detected = ifelse(observed_pfs_biomarker_results[5] < 0.05, 1, 0)
  pfs_biomarker_type1_error = ifelse(pfs_biomarker_difference_detected ==1 & pfs_biomarker_difference_exists == 0, 1, 0)
  pfs_biomarker_type2_error = ifelse(pfs_biomarker_difference_detected ==0 & pfs_biomarker_difference_exists == 1, 1, 0)
  
  pfs_biomarker_coef_diff = observed_pfs_biomarker_coef - true_pfs_biomarker_coef
  
  
  
  prop_observed = nrow(observed_cohort) / nrow(true_cohort)
  
  os_true_median <- tryCatch({
    os_true_survfit <- survfit(Surv(survival_times, died) ~ 1, data = true_cohort)
    unname(summary(os_true_survfit)$table['median'])
  }, error = function(cond) {return(NA)})
  
  
  os_observed_median <- tryCatch({
    os_observed_survfit <- survfit(Surv(survival_times, died) ~ 1, data = observed_cohort)
    unname(summary(os_observed_survfit)$table['median'])  
  }, error = function(cond) {return(NA)})
  
  # observed minus true
  os_observed_minus_true = os_observed_median - os_true_median
  
  
  pfs_true_median <- tryCatch({
    pfs_true_survfit <- survfit(Surv(pfs_times, pfs_event) ~ 1, data = true_cohort)
    unname(summary(pfs_true_survfit)$table['median'])
  }, error = function(cond) {return(NA)})
  
  
  pfs_observed_median <- tryCatch({
    pfs_observed_survfit <- survfit(Surv(pfs_times, pfs_event) ~ 1, data = observed_cohort)
    unname(summary(pfs_observed_survfit)$table['median'])  
  }, error = function(cond) {return(NA)})
  
  # observed minus true
  os_observed_minus_true = os_observed_median - os_true_median  
  pfs_observed_minus_true = pfs_observed_median - pfs_true_median
  
  
  return(list(prop_observed=prop_observed,
              os_true_median=os_true_median,
              os_observed_median=os_observed_median,
              os_observed_minus_true=os_observed_minus_true,
              pfs_true_median=pfs_true_median,
              pfs_observed_median=pfs_observed_median,
              pfs_observed_minus_true=pfs_observed_minus_true,
              true_os_biomarker_coef=true_os_biomarker_coef,
              observed_os_biomarker_coef=observed_os_biomarker_coef,
              observed_os_biomarker_lower=observed_os_biomarker_lower,
              observed_os_biomarker_upper=observed_os_biomarker_upper,
              os_biomarker_difference_exists=os_biomarker_difference_exists,
              os_biomarker_difference_detected=os_biomarker_difference_detected,
              os_biomarker_type1_error=os_biomarker_type1_error,
              os_biomarker_type2_error=os_biomarker_type2_error,
              os_biomarker_coef_diff=os_biomarker_coef_diff,
              true_pfs_biomarker_coef=true_pfs_biomarker_coef,
              observed_pfs_biomarker_coef=observed_pfs_biomarker_coef,
              observed_pfs_biomarker_lower=observed_pfs_biomarker_lower,
              observed_pfs_biomarker_upper=observed_pfs_biomarker_upper,
              pfs_biomarker_difference_exists=pfs_biomarker_difference_exists,
              pfs_biomarker_difference_detected=pfs_biomarker_difference_detected,
              pfs_biomarker_type1_error=pfs_biomarker_type1_error,
              pfs_biomarker_type2_error=pfs_biomarker_type2_error,
              pfs_biomarker_coef_diff=pfs_biomarker_coef_diff))
  
  
}

