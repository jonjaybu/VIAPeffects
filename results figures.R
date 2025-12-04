##################################################################################################
############################## MAIN RESULTS, FIGURES ####################################################
##################################################################################

################################### TRIAL 1

# call arm-specific cuminc at 1, 2, 3 years
results1 %>% group_by(any_tx_x) %>% summarise(ci[time %in% c(52, 52*2, 52*3)]) 

# call arm-specific 95% CIs at 1, 2, 3 years  
bootres1_sum <- bootres1 %>% filter(time %in% c(52, 52*2, 52*3)) %>% group_by(time, any_tx_x) %>% 
  summarise(lower=quantile(ci, c(0.025)), upper=quantile(ci, c(0.975)))

# call ratio CIs, convert to % reductions
bootres1 %>% filter(time %in% c(52, 52*2, 52*3)) %>%
  group_by(time) %>% summarize(ratio=ci[any_tx_x==1]/ci[any_tx_x==0]) %>% #calculate ratio of treatment to control arm cuminc
  group_by(time) %>% summarize(lb=quantile(ratio, 0.025), ub=quantile(ratio, 0.975))  %>% #find 95% CIs
  mutate(effect_lb=1-ub, effect_ub=1-lb) #convert to risk reduction

# figure 2 code
colors=c("#4D4D4D", "#1B9E77")
lines=c("longdash", "solid")

### Standardized marginal cumulative incidence curves
r1 <- ggplot(results1) +
  geom_errorbar(data=bootres1_sum, aes(x=time, group=factor(any_tx_x), color=factor(any_tx_x), 
                                       ymin=lower, ymax=upper, linetype=factor(any_tx_x)), 
                width=10, alpha=0.8, position=position_dodge(width=3)) +
  geom_line(aes(x=time, y=ci, group=factor(any_tx_x), color=factor(any_tx_x), linetype=factor(any_tx_x)), 
            size=1.2) +
  scale_color_manual(values=c(colors[1], colors[2]), labels=c("Nonadherers", "Adherers"), name="") +
  scale_linetype_manual(values = c(lines[1], lines[2]), labels=c("Nonadherers", "Adherers"), name="") +   
  labs(y="Cumulative incidence\n") +
  scale_x_continuous(breaks = seq(0, 156, by=26), labels = seq(0, 3, by=0.5), name="\nYears from index injury") +
  ylim(0, NA) +
  theme_minimal() 

################################### TRIAL 2
  
# call arm-specific cuminc at 1, 2, 3 years
results2 %>% group_by(high_dose_x)%>% summarise(ci[time %in% c(52, 52*2, 52*3)])  

# call arm-specific 95% CIs at 1, 2, 3 years  
bootres2_sum <- bootres2 %>% filter(time %in% c(52, 52*2, 52*3)) %>% group_by(time, high_dose_x) %>% 
  summarise(lower=quantile(ci, c(0.025)), upper=quantile(ci, c(0.975)))

# call ratio CIs, convert to % reductions
bootres2 %>% filter(time %in% c(52, 52*2, 52*3)) %>%
  group_by(time) %>% summarize(ratio=ci[high_dose_x==1]/ci[high_dose_x==0]) %>%
  group_by(time) %>% summarize(lb=quantile(ratio, 0.025), ub=quantile(ratio, 0.975))
  mutate(effect_lb=1-ub, effect_ub=1-lb) #convert to risk reduction
  
# figure 2 code
r2 <- ggplot(results2) +
  geom_errorbar(data=bootres2_sum, aes(x=time, group=factor(high_dose_x), color=factor(high_dose_x), 
                                       ymin=lower, ymax=upper, linetype=factor(high_dose_x)), 
                width=10, alpha=0.8, position=position_dodge(width=3)) +
  geom_line(aes(x=time, y=ci, group=factor(high_dose_x), color=factor(high_dose_x), linetype=factor(high_dose_x)), 
                size=1.2) +
  scale_color_manual(values=c(colors[1], colors[2]), labels=c("Nonadherers", "Adherers"), name="") +
  scale_linetype_manual(values = c(lines[1], lines[2]), labels=c("Nonadherers", "Adherers"), name="") +   
  labs(y="Cumulative incidence\n") +
  scale_x_continuous(breaks = seq(0, 156, by=26), labels = seq(0, 3, by=0.5), name="\nYears from index injury") +
  ylim(0, NA) +
  theme_minimal()
