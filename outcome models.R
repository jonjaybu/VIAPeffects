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

#adherers
adherers1 <- baseline[rep(1:nrow(baseline), each=length(time_periods)),]
adherers1$time <- rep(time_periods, times=nrow(baseline))
adherers1$time2 <- adherers1$time * adherers1$time
adherers1$any_tx_x <- 1 #make everyone an adherer
adherers1$atx_time <- adherers1$any_tx_x * adherers1$time
adherers1$atx_time2 <- adherers1$any_tx_x * adherers1$time2
adherers1$p <- 1- predict(adj_plr_fit, newdata = adherers1, type="response") #predicted survival density at each time point
adherers1$s <- ave(adherers1$p, adherers1$clientID, FUN=cumprod) #calculate survival
adherers1$ci <- 1- adherers1$s #calculate cumulative incidence

#nonadherers
nonadherers1 <- adherers1 #same baseline characteristics
nonadherers1$any_tx_x <- 0 #make everyone a nonadherer
nonadherers1$atx_time <- nonadherers1$any_tx_x * nonadherers1$time
nonadherers1$atx_time2 <- nonadherers1$any_tx_x * nonadherers1$time2
nonadherers1$p <- 1- predict(adj_plr_fit, newdata = nonadherers1, type="response")
nonadherers1$s <- ave(nonadherers1$p, nonadherers1$clientID, FUN=cumprod) 
nonadherers1$ci <- 1- nonadherers1$s 

# Step 3. Calculate standardized survival at each time
# Create concatenated dataset, only keep ci, any_tx_x, and time
both1 <- rbind(adherers1, nonadherers1)
both1 <- both1[,c('ci', 'any_tx_x', 'time')]

# Calculate the mean cuminc at each visit within each adherer group
results1 <- aggregate(ci ~ time + any_tx_x, FUN=mean, data=both1)

# Fix labels, for visualization
results1$adhrf <- factor(results1$any_tx_x, labels = c("Nonadherers", "Adherers"))


################################### Trial 2 ###############################################
### note: 4th-degree natural spline showed best AIC vs. linear, quadratic, 4th degree poly, lower/higher degree splines (BS/NS)
adj_plr_fit <- glm(outcome ~ high_dose_x * ns(time, df=4) + #includes interaction using *
                     black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                   data=pdf_c2, family=binomial(), weights=sw_all)

##### Create simulated data where everyone adheres and doesn't adhere 
baseline <- pdf[pdf$time==0,]

#adherers
adherers2 <- baseline[rep(1:nrow(baseline), each=length(time_periods)),]
adherers2$time <- rep(time_periods, times=nrow(baseline))
adherers2$time2 <- adherers2$time * adherers2$time
adherers2$high_dose_x <- 1
adherers2$hdx_time <- adherers2$high_dose_x * adherers2$time
adherers2$hdx_time2 <- adherers2$high_dose_x * adherers2$time2
adherers2$p <- 1- predict(adj_plr_fit, newdata = adherers2, type="response")
adherers2$s <- ave(adherers2$p, adherers2$clientID, FUN=cumprod) 
adherers2$ci <- 1- adherers2$s  

#nonadherers
nonadherers2 <- adherers2
nonadherers2$high_dose_x <- 0
nonadherers2$hdx_time <- nonadherers2$high_dose_x * nonadherers2$time
nonadherers2$hdx_time2 <- nonadherers2$high_dose_x * nonadherers2$time2
nonadherers2$p <- 1- predict(adj_plr_fit, newdata = nonadherers2, type="response")
nonadherers2$s <- ave(nonadherers2$p, nonadherers2$clientID, FUN=cumprod) 
nonadherers2$ci <- 1- nonadherers2$s 

# Step 3. Calculate standardized survival at each time
both2 <- rbind(adherers2, nonadherers2)
both2 <- both2[,c('ci', 'high_dose_x', 'time')]

# Calculate the mean cuminc at each visit within each adherer group
results2 <- aggregate(ci ~ time + high_dose_x, FUN=mean, data=both2)

# Fix labels, for visualization
results2$adhrf <- factor(results2$high_dose_x, labels = c("Nonadherers", "Adherers"))


