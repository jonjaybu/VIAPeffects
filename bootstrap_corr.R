##################################################################################################
############################## BOOTSTRAP ####################################################
##################################################################################

n_sample <- length(unique(c1strap$clientID))

######################### ONLY TRIAL 1 BOOTSTRAP SHOWN BELOW -- SAME APPROACH FOR TRIAL 2 

strap1 <- function(pdf_c1, j){
  ## randomly resample individuals
  resample <- sample(unique(pdf_c1$clientID), n_sample, replace = T) 
  c1strap <- data.frame(clientID = resample)
  c1strap <- left_join(c1strap, pdf_c1, relationship = "many-to-many")
  
  ############# REPEAT WEIGHTS, MODELING PROCEDURE WITHIN EACH SAMPLE  
  ####### IPW weight calculation
  ## protocol censoring
  d_Cproto_mod2 <- glm(not_Cproto ~ time + time2 + 
                         black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                       data=c1strap, family=binomial())
  n_Cproto_mod2 <- glm(not_Cproto ~ time + time2,
                       data=c1strap, family=binomial())
  
  c1strap$risk_proto <- 0
  c1strap$risk_proto[c1strap$time<=4] <- 1
  
  c1strap$d_Cproto_mod2 <- predict(d_Cproto_mod2, type="response")
  c1strap$d_Cproto_mod2[c1strap$risk_proto==0] <- 1
  
  c1strap$n_Cproto_mod2 <- predict(n_Cproto_mod2, type="response")
  c1strap$n_Cproto_mod2[c1strap$risk_proto==0] <- 1
  
  c1strap$w_proto <- c1strap$n_Cproto_mod2/c1strap$d_Cproto_mod2  ##stabilized weights
  c1strap <- c1strap %>% group_by(clientID) %>% mutate(sw_proto = cumprod(w_proto))
  
  c1strap$sw_proto[c1strap$sw_proto < quantile(c1strap$sw_proto, 0.01)] <- quantile(c1strap$sw_proto, 0.01)
  c1strap$sw_proto[c1strap$sw_proto > quantile(c1strap$sw_proto, 0.99)] <- quantile(c1strap$sw_proto, 0.99)
  
  ### administrative censoring
  d_Cadmin_mod2 <- glm(not_Cadmin ~ time + time2 + 
                         black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                       data=c1strap, family=binomial())
  n_Cadmin_mod2 <- glm(not_Cadmin ~ time + time2,
                       data=c1strap, family=binomial())
  
  c1strap$d_Cadmin_mod2 <- predict(d_Cadmin_mod2, type="response")
  c1strap$d_Cadmin_mod2[c1strap$time<=52] <- 1 #outside risk window
  c1strap$n_Cadmin_mod2 <- predict(n_Cadmin_mod2, type="response")
  c1strap$n_Cadmin_mod2[c1strap$time<=52] <- 1 #outside risk window
  
  c1strap$w_admin <- c1strap$n_Cadmin_mod2/c1strap$d_Cadmin_mod2 
  c1strap <- c1strap %>% group_by(clientID) %>% mutate(sw_admin = cumprod(w_admin))
  
  c1strap$sw_admin[c1strap$sw_admin < quantile(c1strap$sw_admin, 0.01)] <- quantile(c1strap$sw_admin, 0.01)
  c1strap$sw_admin[c1strap$sw_admin > quantile(c1strap$sw_admin, 0.99)] <- quantile(c1strap$sw_admin, 0.99)
  
  ### death censoring
  d_Cdeath_mod2 <- glm(not_Cdeath ~ time + time2 + 
                         black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                       data=c1strap, family=binomial())
  n_Cdeath_mod2 <- glm(not_Cdeath ~ time + time2,
                       data=c1strap, family=binomial())
  
  c1strap$d_Cdeath_mod2 <- predict(d_Cdeath_mod2, type="response")
  c1strap$n_Cdeath_mod2 <- predict(n_Cdeath_mod2, type="response")
  c1strap$w_death <- c1strap$n_Cdeath_mod2/c1strap$d_Cdeath_mod2 
  c1strap <- c1strap %>% group_by(clientID) %>% mutate(sw_death = cumprod(w_death))
  
  c1strap$sw_death[c1strap$sw_death < quantile(c1strap$sw_death, 0.01)] <- quantile(c1strap$sw_death, 0.01)
  c1strap$sw_death[c1strap$sw_death > quantile(c1strap$sw_death, 0.99)] <- quantile(c1strap$sw_death, 0.99)
  
  ### combine all weights
  c1strap <- c1strap %>% mutate(sw_all = sw_proto*sw_admin*sw_death)
  
  ########### Weighted Cumulative Incidence Curves -------------------------------
  ##### Estimate weighted outcome regression with interactions
  adj_plr_fit <- glm(outcome ~ any_tx_x + 
                       time + time2 + atx_time + atx_time2 + 
                       black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                     data=c1strap, family=binomial(), weights=sw_all)
  
  ##### Create simulated data where everyone adheres and doesn't adhere
  baseline <- pdf[pdf$time==0,]
  
  #treateds
  treateds1 <- baseline[rep(1:nrow(baseline), each=length(time_periods)),]
  treateds1$time <- rep(time_periods, times=nrow(baseline))
  treateds1$time2 <- treateds1$time * treateds1$time
  treateds1$any_tx_x <- 1
  treateds1$atx_time <- treateds1$any_tx_x * treateds1$time
  treateds1$atx_time2 <- treateds1$any_tx_x * treateds1$time2
  treateds1$p <- 1- predict(adj_plr_fit, newdata = treateds1, type="response")
  treateds1$s <- ave(treateds1$p, treateds1$clientID, FUN=cumprod) 
  treateds1$ci <- 1- treateds1$s  
  
  #controls
  controls1 <- treateds1
  controls1$any_tx_x <- 0
  controls1$atx_time <- controls1$any_tx_x * controls1$time
  controls1$atx_time2 <- controls1$any_tx_x * controls1$time2
  controls1$p <- 1- predict(adj_plr_fit, newdata = controls1, type="response")
  controls1$s <- ave(controls1$p, controls1$clientID, FUN=cumprod) 
  controls1$ci <- 1- controls1$s 
  
  # Step 3. Calculate standardized survival at each time
  # Create concatenated dataset, only keep s, any_tx_x, and time
  both1 <- rbind(treateds1, controls1)
  both1 <- both1[,c('ci', 'any_tx_x', 'time')]
  
  # Calculate the mean cuminc at each visit within each treated group
  results1 <- aggregate(ci ~ time + any_tx_x, FUN=mean, data=both1)
  results1$iter = j 
  return(results1)
}

### function to replace any failures with NAs
fail_row <- function(k){data.frame(time = NA_real_,
                                   high_dose_x  = NA_real_,
                                   ci= NA_real_,
                                   iter = k)}

## parallelize iterations
library("doSNOW")
cl <- parallel::makeCluster(12)
doParallel::registerDoParallel(cl)
registerDoSNOW(cl)

## Create progress bar
set.seed(1234)
n_iter <- 1000
iterations <- n_iter
pb <- txtProgressBar(min=0, max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## output arm-specific cumulative incidence by time for each replicate
bootres1 <- foreach(k = 1:iterations,
                    .combine = 'rbind',
                    .packages = c('tidyverse', 'splines'),
                    .options.snow = opts) %dopar% {
                      tryCatch({
                        out <- strap1(pdf_c1, k)   
                        out
                      }, error = function(e) fail_row(k)) ### add a fail row for any errors
                    }
close(pb)
stopCluster(cl)
