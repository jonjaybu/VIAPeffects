##################################################################################################
############################## OUTCOME MODELS ####################################################
##################################################################################
library(tidyverse)
library(splines)

## note: code below was adapted from Murray, Caniglia & Petito, Causal Survival Analysis workshop code, available at https://github.com/eleanormurray/CausalSurvivalAnalysisWorkshop
## Jay et al are responsible for any errors

## pdf is the panel dataframe without clones; pdf_c1 is the panel dataframe, with clones, for trial 1; pdf_c2 is the same for trial 2

################################### Trial 1 ###############################################

## Create interaction terms between exposure (any_tx_x) and time, to allow distinct trajectories
pdf_c1$atx_time <- pdf_c1$any_tx_x * pdf_c1$time
pdf_c1$atx_time2 <- pdf_c1$any_tx_x * pdf_c1$time2

## Fit logistic regression w/ IP weights
### note: in tests (not shown), quadratic time terms showed best AIC vs. linear, 4th degree poly, splines (BS/NS)
adj_plr_fit <- glm(outcome ~ any_tx_x + #treatment indicator 
                     time + time2 + #time and quadratic time terms
                     atx_time + atx_time2 + #interacted w/ treatment
                     black + latinx + age + cismale + #demographics
                     factor(inj_type=="GSW") + #whether index injury was gunshot (vs. stab)
                     admit + #whether admitted to hospital for index injury
                     ISS_none + #whether trauma activation for index injury
                     MH + #whether mental health disorder dx pre-injury
                     factor(inj_year)+  #injury year indicators
                     ice, #neighborhood race-poverty index of concentration at the extremes
                   data=pdf_c1, 
                   family=binomial(), 
                   weights=sw_all) #stabilized IP weights combining protocol adherence, admin censoring, death censoring 

##### Create simulated data where everyone adheres and doesn't adhere
## baseline values
baseline <- pdf[pdf$time==0,]

#treateds
treateds1 <- baseline[rep(1:nrow(baseline), each=length(time_periods)),]
treateds1$time <- rep(time_periods, times=nrow(baseline))
treateds1$time2 <- treateds1$time * treateds1$time
treateds1$any_tx_x <- 1 #make everyone treated
treateds1$atx_time <- treateds1$any_tx_x * treateds1$time
treateds1$atx_time2 <- treateds1$any_tx_x * treateds1$time2
treateds1$p <- 1- predict(adj_plr_fit, newdata = treateds1, type="response") #predicted survival density at each time point
treateds1$s <- ave(treateds1$p, treateds1$clientID, FUN=cumprod) #calculate survival
treateds1$ci <- 1- treateds1$s #calculate cumulative incidence

#controls
controls1 <- treateds1 #same baseline characteristics
controls1$any_tx_x <- 0 #make everyone a control
controls1$atx_time <- controls1$any_tx_x * controls1$time
controls1$atx_time2 <- controls1$any_tx_x * controls1$time2
controls1$p <- 1- predict(adj_plr_fit, newdata = controls1, type="response")
controls1$s <- ave(controls1$p, controls1$clientID, FUN=cumprod) 
controls1$ci <- 1- controls1$s 

# Step 3. Calculate standardized survival at each time
# Create concatenated dataset, only keep ci, any_tx_x, and time
both1 <- rbind(treateds1, controls1)
both1 <- both1[,c('ci', 'any_tx_x', 'time')]

# Calculate the mean cuminc at each visit within each treated group
results1 <- aggregate(ci ~ time + any_tx_x, FUN=mean, data=both1)


################################### Trial 2 ###############################################
### note: 4th-degree natural spline showed best AIC vs. linear, quadratic, 4th degree poly, lower/higher degree splines (BS/NS)
adj_plr_fit <- glm(outcome ~ high_dose_x * ns(time, df=4) + #includes interaction using *
                     black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                   data=pdf_c2, family=binomial(), weights=sw_all)

##### Create simulated data where everyone adheres and doesn't adhere 
baseline <- pdf[pdf$time==0,]

#treateds
treateds2 <- baseline[rep(1:nrow(baseline), each=length(time_periods)),]
treateds2$time <- rep(time_periods, times=nrow(baseline))
treateds2$time2 <- treateds2$time * treateds2$time
treateds2$high_dose_x <- 1
treateds2$hdx_time <- treateds2$high_dose_x * treateds2$time
treateds2$hdx_time2 <- treateds2$high_dose_x * treateds2$time2
treateds2$p <- 1- predict(adj_plr_fit, newdata = treateds2, type="response")
treateds2$s <- ave(treateds2$p, treateds2$clientID, FUN=cumprod) 
treateds2$ci <- 1- treateds2$s  

#controls
controls2 <- treateds2
controls2$high_dose_x <- 0
controls2$hdx_time <- controls2$high_dose_x * controls2$time
controls2$hdx_time2 <- controls2$high_dose_x * controls2$time2
controls2$p <- 1- predict(adj_plr_fit, newdata = controls2, type="response")
controls2$s <- ave(controls2$p, controls2$clientID, FUN=cumprod) 
controls2$ci <- 1- controls2$s 

# Step 3. Calculate standardized survival at each time
both2 <- rbind(treateds2, controls2)
both2 <- both2[,c('ci', 'high_dose_x', 'time')]

# Calculate the mean cuminc at each visit within each treated group
results2 <- aggregate(ci ~ time + high_dose_x, FUN=mean, data=both2)



