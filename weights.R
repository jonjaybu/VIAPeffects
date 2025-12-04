##################################################################################################
############################## WEIGHTS ####################################################
##################################################################################################
library(dplyr)

######################### TRIAL 1 #################################################

### add a squared time term for quadratic polynomial
pdf_c1$time2 <- pdf_c1$time * pdf_c1$time

############ a) probability of being protocol-censored by week
# rich denominator model (see documentation for variable definitions)
d_Cproto_mod2 <- glm(not_Cproto ~ time + time2 + #includes quadratic time term
                       black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                     data=pdf_c1, family=binomial())
# simple numerator model for stabilization
n_Cproto_mod2 <- glm(not_Cproto ~ time + time2,
                     data=pdf_c1, family=binomial())

## assign weights
# define risk window for protocol censoring
pdf_c1$risk_proto <- 0 #default
pdf_c1$risk_proto[pdf_c1$time<=4] <- 1 #only at risk of protocol censoring during weeks 0-4

#predict fitted models
pdf_c1$d_Cproto_mod2 <- predict(d_Cproto_mod2, type="response") #denominators
pdf_c1$d_Cproto_mod2[pdf_c1$risk_proto==0] <- 1 #set weights to exactly 1 outside of risk windows for protocol

pdf_c1$n_Cproto_mod2 <- predict(n_Cproto_mod2, type="response") #numerators
pdf_c1$n_Cproto_mod2[pdf_c1$risk_proto==0] <- 1 #set weights to exactly 1 outside of risk windows for protocol

pdf_c1$w_proto <- pdf_c1$n_Cproto_mod2/pdf_c1$d_Cproto_mod2  ##stabilized weights
pdf_c1 <- pdf_c1 %>% group_by(clientID) %>% mutate(sw_proto = cumprod(w_proto)) ##take cumulative product for IPW at each time point

#truncate at 1st and 99th percentiles
pdf_c1$sw_proto[pdf_c1$sw_proto < quantile(pdf_c1$sw_proto, 0.01)] <- quantile(pdf_c1$sw_proto, 0.01)
pdf_c1$sw_proto[pdf_c1$sw_proto > quantile(pdf_c1$sw_proto, 0.99)] <- quantile(pdf_c1$sw_proto, 0.99)

############ b) probability of being administratively censored by week 
### same general approach as above
d_Cadmin_mod2 <- glm(not_Cadmin ~ time + time2 + 
                       black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                     data=pdf_c1, family=binomial())
n_Cadmin_mod2 <- glm(not_Cadmin ~ time + time2,
                     data=pdf_c1, family=binomial())

## assign weights
pdf_c1$d_Cadmin_mod2 <- predict(d_Cadmin_mod2, type="response")
pdf_c1$d_Cadmin_mod2[pdf_c1$time<=52] <- 1 #outside risk window -- there's no risk of being admin censored < week 52
pdf_c1$n_Cadmin_mod2 <- predict(n_Cadmin_mod2, type="response")
pdf_c1$n_Cadmin_mod2[pdf_c1$time<=52] <- 1 #outside risk window

pdf_c1$w_admin <- pdf_c1$n_Cadmin_mod2/pdf_c1$d_Cadmin_mod2 #stabilized weights
pdf_c1 <- pdf_c1 %>% group_by(clientID) %>% mutate(sw_admin = cumprod(w_admin)) #take cumulative product

#truncate values
pdf_c1$sw_admin[pdf_c1$sw_admin < quantile(pdf_c1$sw_admin, 0.01)] <- quantile(pdf_c1$sw_admin, 0.01)
pdf_c1$sw_admin[pdf_c1$sw_admin > quantile(pdf_c1$sw_admin, 0.99)] <- quantile(pdf_c1$sw_admin, 0.99)

############ c) probability of being death-censored by week 
d_Cdeath_mod2 <- glm(not_Cdeath ~ time + time2 + 
                       black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                     data=pdf_c1, family=binomial())
n_Cdeath_mod2 <- glm(not_Cdeath ~ time + time2,
                     data=pdf_c1, family=binomial())

pdf_c1$d_Cdeath_mod2 <- predict(d_Cdeath_mod2, type="response")
pdf_c1$n_Cdeath_mod2 <- predict(n_Cdeath_mod2, type="response")
pdf_c1$w_death <- pdf_c1$n_Cdeath_mod2/pdf_c1$d_Cdeath_mod2 
pdf_c1 <- pdf_c1 %>% group_by(clientID) %>% mutate(sw_death = cumprod(w_death))

pdf_c1$sw_death[pdf_c1$sw_death < quantile(pdf_c1$sw_death, 0.01)] <- quantile(pdf_c1$sw_death, 0.01)
pdf_c1$sw_death[pdf_c1$sw_death > quantile(pdf_c1$sw_death, 0.99)] <- quantile(pdf_c1$sw_death, 0.99)

######################### combine all weights ######################### 
pdf_c1 <- pdf_c1 %>% mutate(sw_all = sw_proto*sw_admin*sw_death) ### outcome models will be weighted by sw_all

######################### TRIAL 2 #################################################
######################### SEE TRIAL 1 CODE FOR CLOSER ANNOTATION ##################

### add a squared time term for quadratic polynomial
pdf_c2$time2 <- pdf_c2$time * pdf_c2$time

####### PROTOCOL CENSORING
d_Cproto_mod2 <- glm(not_Cproto ~ time + time2 + 
                       black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                     data=pdf_c2, family=binomial())
n_Cproto_mod2 <- glm(not_Cproto ~ time + time2,
                     data=pdf_c2, family=binomial())

## at-risk periods differ between trials for treatment group; same for control group
pdf_c2$risk_proto <- 0 
pdf_c2$risk_proto[pdf_c2$high_dose_x==0 & pdf_c2$time<=4] <- 1 #control clones at risk until t=4
pdf_c2$risk_proto[pdf_c2$high_dose_x==1 & pdf_c2$time<=8] <- 1 #treatment clones at risk until t=8

pdf_c2$d_Cproto_mod2 <- predict(d_Cproto_mod2, type="response")
pdf_c2$d_Cproto_mod2[pdf_c2$risk_proto==0] <- 1

pdf_c2$n_Cproto_mod2 <- predict(n_Cproto_mod2, type="response")
pdf_c2$n_Cproto_mod2[pdf_c2$risk_proto==0] <- 1

pdf_c2$w_proto <- pdf_c2$n_Cproto_mod2/pdf_c2$d_Cproto_mod2 
pdf_c2 <- pdf_c2 %>% group_by(clientID) %>% mutate(sw_proto = cumprod(w_proto))

pdf_c2$sw_proto[pdf_c2$sw_proto < quantile(pdf_c2$sw_proto, 0.01)] <- quantile(pdf_c2$sw_proto, 0.01)
pdf_c2$sw_proto[pdf_c2$sw_proto > quantile(pdf_c2$sw_proto, 0.99)] <- quantile(pdf_c2$sw_proto, 0.99)

####### ADMIN CENSORING -- SAME AS TRIAL 1
d_Cadmin_mod2 <- glm(not_Cadmin ~ time + time2 + 
                       black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                     data=pdf_c2, family=binomial())
n_Cadmin_mod2 <- glm(not_Cadmin ~ time + time2,
                     data=pdf_c2, family=binomial())

pdf_c2$d_Cadmin_mod2 <- predict(d_Cadmin_mod2, type="response")
pdf_c2$d_Cadmin_mod2[pdf_c2$time<=52] <- 1 #outside risk window
pdf_c2$n_Cadmin_mod2 <- predict(n_Cadmin_mod2, type="response")
pdf_c2$n_Cadmin_mod2[pdf_c2$time<=52] <- 1 #outside risk window

pdf_c2$w_admin <- pdf_c2$n_Cadmin_mod2/pdf_c2$d_Cadmin_mod2 
pdf_c2 <- pdf_c2 %>% group_by(clientID) %>% mutate(sw_admin = cumprod(w_admin))

pdf_c2$sw_admin[pdf_c2$sw_admin < quantile(pdf_c2$sw_admin, 0.01)] <- quantile(pdf_c2$sw_admin, 0.01)
pdf_c2$sw_admin[pdf_c2$sw_admin > quantile(pdf_c2$sw_admin, 0.99)] <- quantile(pdf_c2$sw_admin, 0.99)

####### DEATH CENSORING -- SAME AS TRIAL 1
d_Cdeath_mod2 <- glm(not_Cdeath ~ time + time2 + 
                       black + latinx + age + cismale + factor(inj_type=="GSW") + admit + ISS_none + MH + factor(inj_year)+ ice,
                     data=pdf_c2, family=binomial())
n_Cdeath_mod2 <- glm(not_Cdeath ~ time + time2,
                     data=pdf_c2, family=binomial())

pdf_c2$d_Cdeath_mod2 <- predict(d_Cdeath_mod2, type="response")
pdf_c2$n_Cdeath_mod2 <- predict(n_Cdeath_mod2, type="response")
pdf_c2$w_death <- pdf_c2$n_Cdeath_mod2/pdf_c2$d_Cdeath_mod2 
pdf_c2 <- pdf_c2 %>% group_by(clientID) %>% mutate(sw_death = cumprod(w_death))

pdf_c2$sw_death[pdf_c2$sw_death < quantile(pdf_c2$sw_death, 0.01)] <- quantile(pdf_c2$sw_death, 0.01)
pdf_c2$sw_death[pdf_c2$sw_death > quantile(pdf_c2$sw_death, 0.99)] <- quantile(pdf_c2$sw_death, 0.99)

###################### COMBINE ALL WEIGHTS #########################
pdf_c2 <- pdf_c2 %>% mutate(sw_all = sw_proto*sw_admin*sw_death)
