#!/usr/bin/env Rscript
args = commandArgs()

setwd(".")
packages <- c("devtools","here","ranger","MatchIt","parallel","boot","caret")
userLib <- "~/R/R_LIBS_USER"
.libPaths(userLib)

ifelse(!dir.exists(userLib), dir.create(userLib), FALSE)

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN') 
  }
}

for (package in packages) {
  library(package, character.only=T, lib.loc=.libPaths())
}

devtools::install_github("hadley/tidyverse")
devtools::install_github("hadley/ggplot2")

for (package in c("tidyverse","ggplot2")) {
  library(package, character.only=T)
}

start_time <- Sys.time()

expit <- function(x){1/(1+exp(-x))}
itr=2000
nSim <- 1:2

print(paste("Numer of simulations:",nSim))

set.seed(123)
#' Load data that will be used for plasmode:

eager_bl <- read_csv(here("data","eager_base_imputed.csv")) %>% 
  mutate(meanAP = (BPS + 2*BPD)/3,
         liveBirth = as.numeric(outcome=="live birth"),
         time_try_pregnant=if_else(as.numeric(time_try_pregnant)>10,10,as.numeric(time_try_pregnant)),
         eligibility = as.numeric(eligibility=="new")) %>% 
  select(age,time_try_pregnant,BMI,meanAP,treatment,eligibility,liveBirth)
eager_bl


#' Fit PS model and change parameters to get simulation model
mod_ps <- glm(treatment ~ BMI + meanAP + eligibility + age + time_try_pregnant,data=eager_bl,family=binomial("logit"))

coef_ps <- coef(mod_ps)
coef_ps[1] <- coef_ps[1] - 4
coef_ps[6] <- coef_ps[6] + .5
coef_ps

#' Fit outcome model and change parameters to get simulation model
mod_mu <- glm(liveBirth ~ treatment + BMI + meanAP + eligibility + age + time_try_pregnant,data=eager_bl,family=binomial("logit"))
summary(mod_mu)

coef_mu <- coef(mod_mu)
coef_mu[2] <- log(1.75) # true conditional treatment effect on OR scale
coef_mu[1] <- coef_mu[1] - 1.25
coef_mu

#' Get true marginal treatment effect on OR scale using g computation and 1e6 observations
index <- sample(1:nrow(eager_bl),1e7,replace=T)
sim_cov <- eager_bl[index,] %>% select(BMI,meanAP,eligibility,age,time_try_pregnant) %>% as.matrix()
sim_treat <- rbinom(nrow(sim_cov),1,expit(cbind(1,sim_cov)%*%coef_ps))
sim_y <- rbinom(nrow(sim_cov),1,expit(cbind(1,sim_treat,sim_cov)%*%coef_mu))

print(paste("Mean of the exposure:",mean(sim_treat)))
print(paste("Mean of the outcome:",mean(sim_y)))

big_dat <- tibble(sim_y,sim_treat,data.frame(sim_cov))
names(big_dat) <- c("outcome","treatment","BMI","meanAP","eligibility","age","time_try_pregnant")

mod_true <- glm(outcome ~ treatment + BMI + meanAP + eligibility + age + time_try_pregnant, data=big_dat, family=binomial("logit"))
mu1 <- mean(predict(mod_true,newdata=transform(big_dat,treatment=1),type="response"))
mu0 <- mean(predict(mod_true,newdata=transform(big_dat,treatment=0),type="response"))

true_or <- log((mu1/(1-mu1))/(mu0/(1-mu0)))
true_or

# sample data
index <- sample(1:nrow(eager_bl),nrow(eager_bl),replace=T)
sim_cov <- eager_bl[index,] %>% select(BMI,meanAP,eligibility,age,time_try_pregnant) %>% as.matrix()
sim_treat <- rbinom(nrow(sim_cov),1,expit(cbind(1,sim_cov)%*%coef_ps))
sim_y <- rbinom(nrow(sim_cov),1,expit(cbind(1,sim_treat,sim_cov)%*%coef_mu))

e <- tibble(sim_y,sim_treat,data.frame(sim_cov))
names(e) <- c("outcome","treatment","BMI","meanAP","eligibility","age","time_try_pregnant")

#' create simulation function that implements g computation with a regression estimator, random forest, and matching procedure
nCore <- detectCores()

print(paste("Numbre of cores:",nCore))

g_boot <- function(x){
  set.seed(x)
  print(paste("Iteration:",x))
  # simulation
  index <- sample(1:nrow(eager_bl),nrow(eager_bl),replace=T)
  sim_cov <- eager_bl[index,] %>% select(BMI,meanAP,eligibility,age,time_try_pregnant) %>% as.matrix()
  sim_treat <- rbinom(nrow(sim_cov),1,expit(cbind(1,sim_cov)%*%coef_ps))
  sim_y <- rbinom(nrow(sim_cov),1,expit(cbind(1,sim_treat,sim_cov)%*%coef_mu))
  
  e <- tibble(sim_y,sim_treat,data.frame(sim_cov))
  names(e) <- c("outcome","treatment","BMI","meanAP","eligibility","age","time_try_pregnant")
  
  # g Computation (bootstrap valid)
  mod_gComp <- glm(outcome ~ treatment + BMI + meanAP + eligibility + age + time_try_pregnant, data=e, family=binomial("logit"))
  mu1_gComp <- mean(predict(mod_true,newdata=transform(e,treatment=1),type="response"))
  mu0_gComp <- mean(predict(mod_true,newdata=transform(e,treatment=0),type="response"))
  
  gComp_or <- log((mu1_gComp/(1-mu1_gComp))/(mu0_gComp/(1-mu0_gComp)))
  
  # PS matching (bootstrap not valid)
  mod_match <- matchit(treatment ~ BMI + meanAP + eligibility + age + time_try_pregnant, method = "nearest", data = e)
  mu1_match <- mean(match.data(mod_match, "treat")$outcome)
  mu0_match <- mean(match.data(mod_match, "control")$outcome)
  
  match_or <- log((mu1_match/(1-mu1_match))/(mu0_match/(1-mu0_match)))
  
  # Random Forest (bootstrap not valid)
  tgrid <- expand.grid(
    .mtry = 1:5,
    .splitrule = "gini",
    .min.node.size = seq(5,20,by=5)
  )
  e <- e %>% mutate(outcome1=if_else(outcome==1,"event","nonevent"))
  mod_rf0 <- train(outcome1 ~ BMI + meanAP + eligibility + age + time_try_pregnant, 
                   data = subset(e,treatment==0),
                   method = "ranger",
                   trControl = trainControl(method="cv", number = 5, verboseIter = F, classProbs = T),
                   tuneGrid = tgrid,
                   num.trees = 250)
  mod_rf1 <- train(outcome1 ~ BMI + meanAP + eligibility + age + time_try_pregnant, 
                   data = subset(e,treatment==1),
                   method = "ranger",
                   trControl = trainControl(method="cv", number = 5, verboseIter = F, classProbs = T),
                   tuneGrid = tgrid,
                   num.trees = 250)

  # mod_rf0 <- ranger(factor(outcome) ~ BMI + meanAP + eligibility + age + time_try_pregnant, 
  #                   data = subset(e,treatment==0), min.node.size = 40,
  #                   num.trees = 500,mtry=3,splitrule="gini")
  # mod_rf1 <- ranger(factor(outcome) ~ BMI + meanAP + eligibility + age + time_try_pregnant, 
  #                   data = subset(e,treatment==1), min.node.size = 40,
  #                   num.trees = 500,mtry=3,splitrule="gini")
  
  # mu1_rf <- mean(as.numeric(predict(mod_rf1,data=e)$pred=="1"))
  # mu0_rf <- mean(as.numeric(predict(mod_rf0,data=e)$pred=="0"))
  
  mu1_rf <- mean(predict(mod_rf1,data=e,type="prob")[,1])
  mu0_rf <- mean(predict(mod_rf0,data=e,type="prob")[,1])
  rf_or <- log((mu1_rf/(1-mu1_rf))/(mu0_rf/(1-mu0_rf)))

  # bootstrap function
  sim_boot_fun <- function(data,j,
                           .mtry0=mod_rf0$bestTune$mtry,
                           .min.node.size0=mod_rf0$bestTune$min.node.size,
                           .mtry1=mod_rf1$bestTune$mtry,
                           .min.node.size1=mod_rf1$bestTune$min.node.size){
    temp <- data[j,]

    # g Computation (bootstrap valid)
    mod_gComp <- glm(outcome ~ treatment + BMI + meanAP + eligibility + age + time_try_pregnant, data=temp, family=binomial("logit"))
    mu1_gComp <- mean(predict(mod_gComp,newdata=transform(temp,treatment=1),type="response"))
    mu0_gComp <- mean(predict(mod_gComp,newdata=transform(temp,treatment=0),type="response"))
    
    est_gComp <- log((mu1_gComp/(1-mu1_gComp))/(mu0_gComp/(1-mu0_gComp)))
    
    # PS matching (bootstrap not valid)
    mod_match <- matchit(treatment ~ BMI + meanAP + eligibility + age + time_try_pregnant, method = "nearest", data = temp)
    mu1_match <- mean(match.data(mod_match, "treat")$outcome)
    mu0_match <- mean(match.data(mod_match, "control")$outcome)
    
    est_match <- log((mu1_match/(1-mu1_match))/(mu0_match/(1-mu0_match)))
    
    # Random Forest (bootstrap not valid)
    mod_rf0 <- ranger(factor(outcome) ~ BMI + meanAP + eligibility + age + time_try_pregnant, 
                      data = subset(temp,treatment==0),probability=T, min.node.size = .min.node.size0,
                      num.trees = 250,mtry=.mtry0,splitrule="gini")
    mod_rf1 <- ranger(factor(outcome) ~ BMI + meanAP + eligibility + age + time_try_pregnant, 
                      data = subset(temp,treatment==1),probability=T, min.node.size = .min.node.size1,
                      num.trees = 250,mtry=.mtry1,splitrule="gini")
    
    mu1_rf <- mean(predict(mod_rf1,data=e,type="response")$pred[,2])
    mu0_rf <- mean(predict(mod_rf0,data=e,type="response")$pred[,2])
    est_rf <- log((mu1_rf/(1-mu1_rf))/(mu0_rf/(1-mu0_rf)))
    
    est <- c(est_gComp,est_match,est_rf)
    
    return(est)
  }
  #bootstrapping
  boot_temp <- boot(data=e,statistic = sim_boot_fun, R=itr,parallel = "multicore", ncpus=nCore)
  ci_gComp <- boot.ci(boot.out = boot_temp, index = 1, type = c("norm","perc","bca")) ###CAUTION, index =1
  ci_match <- boot.ci(boot.out = boot_temp, index = 2, type = c("norm","perc","bca")) ###CAUTION, index =2
  ci_rf <- boot.ci(boot.out = boot_temp, index = 3, type = c("norm","perc","bca")) ###CAUTION, index =3
  #save output
  gComp_res <- tibble(MC=x, 
                      estimate=gComp_or,
                      t0=boot_temp$t0[1],
                      t=mean(boot_temp$t[,1]),
                      se=sd(boot_temp$t[,1]),
                      wald_lcl=ci_gComp$normal[2],
                      wald_ucl=ci_gComp$normal[3],
                      perc_lcl=ci_gComp$percent[4],
                      perc_ucl=ci_gComp$percent[5],
                      bca_lcl =ci_gComp$bca[4],
                      bca_ucl =ci_gComp$bca[5])
  match_res <- tibble(MC=x, 
                      estimate=match_or,
                      t0=boot_temp$t0[2],
                      t=mean(boot_temp$t[,2]),
                      se=sd(boot_temp$t[,2]),
                      wald_lcl=ci_match$normal[2],
                      wald_ucl=ci_match$normal[3],
                      perc_lcl=ci_match$percent[4],
                      perc_ucl=ci_match$percent[5],
                      bca_lcl =ci_match$bca[4],
                      bca_ucl =ci_match$bca[5])
  rf_res <-    tibble(MC=x,
                      estimate=rf_or,
                      t0=boot_temp$t0[3],
                      t=mean(boot_temp$t[,3]),
                      se=sd(boot_temp$t[,3]),
                      wald_lcl=ci_rf$normal[2],
                      wald_ucl=ci_rf$normal[3],
                      perc_lcl=ci_rf$percent[4],
                      perc_ucl=ci_rf$percent[5],
                      bca_lcl =ci_rf$bca[4],
                      bca_ucl =ci_rf$bca[5])
  
  est1 <- rbind(gComp_res,match_res,rf_res)
  
  est2 <- tibble(MC=x,
                 t_gComp=boot_temp$t[,1],
                 t_Match=boot_temp$t[,2],
                 t_rF=boot_temp$t[,3])

  if(x==1){
    write_csv(est2,here(paste0("output/bs_distributions",max(nSim),".csv")),append = F,col_names = T)
  } else{
    write_csv(est2,here(paste0("output/bs_distributions",max(nSim),".csv")),append = T,col_names = F)
  }
  
  return(est1)
    
}

system.time(res0 <- lapply(nSim,function(a) g_boot(a)))
res <-  bind_rows(res0)
nn <- c("gComp","Match","rF")
res$type <- rep(nn,max(nSim))
res <- res %>% mutate(wald2_lcl = t0 - 1.96*se,
                      wald2_ucl = t0 + 1.96*se)

print(res)

print(max(nSim))

print(paste0("output/sim_res",max(nSim),".csv"))

write_csv(res,here(paste0("output/sim_res",max(nSim),".csv")),append = F,col_names = T)

run_time <- Sys.time() - start_time

print(paste("Run Time:",run_time))

