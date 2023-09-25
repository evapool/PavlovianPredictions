################################################################################


#         Neural substrates of parallel devaluation-sensitive and 
#           devaluation-insensitive Pavlovian learning in humans 


#                        behavioral results

# last  modified on MAY 2023
################################################################################




#-------------------------------------------------------------------------------
#                   PRELIMINARY STUFF
#-------------------------------------------------------------------------------

# load libraries
if(!require(pacman)) {
  install.packages('pacman')
  library(pacman)
}


pacman::p_load(car, effects,afex, jtools, ez, emmeans,
               lme4, lmerTest,nlme, psych, rstudioapi, optimx,pbkrtest,
               tidyverse, plyr, dplyr, tidyr, reshape, reshape2, pastecs,Rmisc,
               devtools, 
               viridis, ggplot2, grid, gridExtra,cowplot,corrplot,DescTools,ggExtra,
               BayesFactor, bayestestR, DHARMa, brms)


require(lattice)


# set path
current_dir <- dirname(getActiveDocumentContext()$path)
setwd(current_dir)
home_path <- getwd()
setwd(home_path)

# run set up
source(file.path(home_path,'R','setup.R'), echo=F)# useful functions



#-------------------------------------------------------------------------------
#                          LEARNING PUPIL CS
#-------------------------------------------------------------------------------

# ------------------------------------------- statistical model PUPIL
mf_pupil = formula(CS_pupil ~ CS.value*run + (CS.value*run|ID))

pupil_mod = lmer(mf_pupil, data = LEARN, control = my_control2) 

summary(pupil_mod)

confint(pupil_mod, level = 0.95, method = "Wald") 


# ------------------------------------------ visualize assumptions check PUPIL
plot(fitted(pupil_mod),residuals(pupil_mod)) 
qqnorm(residuals(pupil_mod ))
hist(residuals(pupil_mod ))
simulateResiduals(fittedModel=pupil_mod, n=1000, plot=TRUE)
plot(density(LEARN$CS_pupil))


#-------------------------------------------- Bayes Factor 
LEARN.p.m$run_b = ifelse(LEARN.p.m$run  == "2", -1, 1) # we need to sum code
LEARN.p.m$run_b = factor(LEARN.p.m$run_b)


fit_pl_m <- brm( scale(CS_pupil) ~ CS.value * run_b + (CS.value * run_b |ID), 
                 prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                 data =LEARN.p.m, # on the aggregate data
                 iter = 40000, warmup=5000,
                 family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_pl_m, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_pl_m, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))




#-------------------------------------------------------------------------------
#                      LEARNING DWELL TIME ANTICIPATION
#-------------------------------------------------------------------------------


# ------------------------------------------- statistical model DW L
mf_dwl = formula(ANT_DW_L ~ CS.left*run + (CS.left*run|ID))

dwl_mod = lmer(mf_dwl, data = LEARN, control = my_control) 

summary(dwl_mod)

confint(dwl_mod, level = 0.95, method = "Wald") 


# ------------------------------------------ visualize assumptions check DW L
plot(fitted(dwl_mod),residuals(dwl_mod)) 
qqnorm(residuals(dwl_mod))
hist(residuals(dwl_mod))
simulateResiduals(fittedModel=dwl_mod, n=1000, plot=TRUE)


#-------------------------------------------- Bayes Factor 
LEARN.dwl.m$run_b = ifelse(LEARN.dwl.m$run  == "2", -1, 1)
LEARN.dwl.m$run_b = factor(LEARN.dwl.m$run_b)


fit_dwl_m <- brm( scale(ANT_DW_L) ~ CS.left * run_b + (CS.left * run_b |ID), 
                  prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                  data =LEARN.dwl.m, # on the aggregate data
                  iter = 40000, warmup=5000,
                  family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_dwl_m, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_dwl_m, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))




# ------------------------------------------- statistical model DW R
mf_dwr = formula(ANT_DW_R ~ CS.right*run + (CS.right*run|ID))

dwr_mod = lmer(mf_dwr, data = LEARN, control = my_control2) 

summary(dwr_mod)

confint(dwr_mod, level = 0.95, method = "Wald") 

# ------------------------------------------ visualize assumptions check DW R
plot(fitted(dwr_mod),residuals(dwr_mod)) 
qqnorm(residuals(dwr_mod))
hist(residuals(dwr_mod))
simulateResiduals(fittedModel=dwr_mod, n=1000, plot=TRUE)



#-------------------------------------------- Bayes Factor 
LEARN.dwr.m$run_b = ifelse(LEARN.dwr.m$run  == "2", -1, 1)
LEARN.dwr.m$run_b = factor(LEARN.dwr.m$run_b)

fit_dwr_m <- brm( scale(ANT_DW_R) ~ CS.right * run_b + (CS.right * run_b |ID), 
                  prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                  data =LEARN.dwr.m, # on the aggregate data
                  iter = 40000, warmup=5000,
                  family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_dwr_m, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_dwr_m, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


#-------------------------------------------------------------------------------
#                           LEARNING REACTION TIMES OUTCOME
#-------------------------------------------------------------------------------

# ------------------------------------------- statistical model SIDE
mf_side = formula(US_RT ~ congr_SIDE*run + (congr_SIDE*run|ID))

side_mod = lmer(mf_side, data = subset(LEARN, congr_ID == "congr" & CS_ID != "CSmi"), control = my_control) 

summary(side_mod)

confint(side_mod, level = 0.95, method = "Wald") 


# ------------------------------------------ visualize assumptions check SIDE
plot(fitted(side_mod),residuals(side_mod)) 
qqnorm(residuals(side_mod ))
hist(residuals(side_mod ))
simulateResiduals(fittedModel=side_mod, n=1000, plot=TRUE)


#-------------------------------------------- Bayes Factor 
SIDE.mean$run_b = ifelse(SIDE.mean$run  == "2", -1, 1)
SIDE.mean$run_b= factor(SIDE.mean$run_b)
SIDE.mean$congr_SIDE_b = ifelse(SIDE.mean$congr_SIDE == "congr", -1, 1)
SIDE.mean$congr_SIDE_b= factor(SIDE.mean$congr_SIDE_b)


fit_rt_side<- brm( scale(US_RT) ~ congr_SIDE * run_b + (congr_SIDE * run_b |ID), 
                   prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                   data =SIDE.mean, # on the aggregate data
                   iter = 40000, warmup=5000,
                   family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_rt_side, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_rt_side, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


# ------------------------------------------- statistical model ID
mf_id = formula(US_RT ~ congr_ID*run + (congr_ID*run|ID))

id_mod = lmer(mf_id, data = subset(LEARN, congr_SIDE == "congr" & CS_ID != "CSmi"), control = my_control) 

summary(id_mod) # this model is inclusive since the aggregated data do not provide the same output - do not drive conclusions form it

confint(id_mod, level = 0.95, method = "Wald") 


# ------------------------------------------ visualize assumptions check ID
plot(fitted(id_mod),residuals(id_mod )) 
qqnorm(residuals(id_mod ))
hist(residuals(id_mod ))
simulateResiduals(fittedModel=id_mod, n=1000, plot=TRUE)


#-------------------------------------------- Bayes Factor 
ID.mean$run_b = ifelse(ID.mean$run  == "2", -1, 1)
ID.mean$run_b= factor(ID.mean$run_b)
ID.mean$congr_ID_b = ifelse(ID.mean$congr_ID == "congr", -1, 1)
ID.mean$congr_ID_b= factor(ID.mean$congr_ID_b)


fit_rt_id<- brm( scale(US_RT) ~ congr_ID * run + (congr_ID * run |ID), 
                 prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                 data =ID.mean, # on the aggregate data
                 iter = 40000, warmup=5000,
                 family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_rt_id, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_rt_id, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


#-------------------------------------------------------------------------------
#         MANIPULATION CHECK FOR THE DEVALUATION PROCEUDRE
#-------------------------------------------------------------------------------


# ------------------------------------------- statistical model
mf = formula(outcome_liking ~ US_value * run + (US_value + run |ID)) # ! simplified error term

check_mod = lmer(mf, data = FOOD.c, control = my_control)

summary(check_mod)

confint(check_mod, level = 0.95, method = "Wald") 


# ------------------------------------------ visualize assumptions check
plot(fitted(check_mod ),residuals(check_mod )) 
qqnorm(residuals(check_mod ))
hist(residuals(check_mod ))
simulateResiduals(fittedModel=check_mod , n=1000, plot=TRUE)
plot(density(FOOD.c$outcome_liking))


#-------------------------------------------- Bayes Factor 
FOOD.c$run_b = ifelse(FOOD.c$run  == "Pre devaluation", -1, 1)
FOOD.c$US_value_b = ifelse(FOOD.c$US_value  == "val", -1, 1)

fit_mc <- brm( scale(outcome_liking) ~ US_value_b * run_b + (US_value_b * run_b |ID), 
               prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
               data =FOOD.c,
               iter = 40000, warmup=5000,
               family = gaussian(), save_pars = save_pars(all = TRUE))

describe_posterior(fit_mc, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_mc, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


# ------------------------------------------- plot raw data
p_deval_check <- ggplot(data = FOOD.c, aes(x = run, y = outcome_liking, fill = factor(US_value, levels = c("val","deval")), color = factor(US_value, levels = c("val","deval")))) +
  geom_abline(slope= 0, intercept=5, linetype = "dashed", color = "gray") +
  facet_grid(~run,scales = "free") +
  geom_point(aes(x = interaction(run, factor(US_value, levels = c("val","deval"))), y = outcome_liking, group = ID), alpha = 0.2, position = position_dodge(0.6))+
  geom_crossbar(data = FOOD.bg, aes(x = interaction(run, factor(US_value, levels = c("val","deval"))), y = outcome_liking, ymin=outcome_liking, ymax=outcome_liking), 
                width = 0.9 , alpha = 0) +
  geom_errorbar(data = FOOD.bg, aes(x = interaction(run, factor(US_value, levels = c("val","deval"))), ymin = outcome_liking-se, ymax = outcome_liking+se), width = 0.9)+
  geom_flat_violin(aes(x = 2.7), alpha = .5, position = position_nudge(x = -.15, y = 0), adjust = 1.8, trim = F, color = NA) +
  scale_x_discrete(labels = c("valued", "devalued")) +
  scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0, 12), expand = c(0, 0))+
  scale_fill_manual(values=c( "#66CC33","#3d3d3d")) +
  scale_color_manual(values=c( "#66CC33","#3d3d3d")) +
  theme_bw() +
  labs( x = "Outcome", y = "Pleasantness of the outcome")


p_deval_check = p_deval_check + timeline_theme + theme(legend.position="none") + theme( panel.grid.major.y = element_blank())

pdf(file.path(figures_path,'Figure_1C.pdf'))
print(p_deval_check)
dev.off() 




#-------------------------------------------------------------------------------
#                 OUTCOME DEVALUATION EFFECTS
#-------------------------------------------------------------------------------

# ------------------------------------------- statistical model PUPIL
mf = formula(CS_pupil_delta ~ CS_ID*run  + (CS_ID + run |ID)) 

deval_pupil_mod = lmer(mf, data = PUPIL.mean, control = my_control)

summary(deval_pupil_mod )

confint(deval_pupil_mod, level = 0.95, method = "Wald") 


# ------------------------------------------ visualize assumptions check 
plot(fitted(deval_pupil_mod,residuals(deval_pupil_mod)) 
     qqnorm(residuals(deval_pupil_mod))
     hist(residuals(deval_pupil_mod))
     simulateResiduals(fittedModel=deval_pupil_mod, n=1000, plot=TRUE)
     
     
#-------------------------------------------- Bayes Factor 
PUPIL.mean$run_b = ifelse(PUPIL.mean$run  == "2", -1, 1)
PUPIL.mean$run_b = factor(PUPIL.mean$run_b)
PUPIL.mean$CS_ID_b = ifelse(PUPIL.mean$CS_ID == "val", -1, 1)
PUPIL.mean$CS_ID_b = factor(PUPIL.mean$CS_ID_b)
     
    
fit_pupil_change<- brm(scale(CS_pupil_delta) ~ CS_ID_b * run_b + (CS_ID_b * run_b |ID), 
                       prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                       data =PUPIL.mean, # on the aggregate data
                       iter = 40000, warmup=5000,
                       family = gaussian(), save_pars = save_pars(all = TRUE))
     
     
describe_posterior(fit_pupil_change, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                    bf_prior = fit_pupil_change, diagnostic = "Rhat",  
                    test = c("p_direction", "bf"))
     
     
# ------------------------------------------- statistical model DW
mf = formula(DW_delta~ CS_ID*run  + (CS_ID + run |ID))

deval_dw_mod = lmer(mf, data = DW.mean, control = my_control)

summary(deval_dw_mod)

confint(deval_dw_mod, level = 0.95, method = "Wald") # anova

# ------------------------------------------ visualize assumptions check 
plot(fitted(deval_dw_mod,residuals(deval_dw_mod)) 
     qqnorm(residuals(deval_dw_mod))
     hist(residuals(deval_dw_mod))
     simulateResiduals(fittedModel=deval_dw_mod, n=1000, plot=TRUE)
     
     


#-------------------------------------------- Bayes Factor 
DW.mean$run_b = ifelse(DW.mean$run  == "2", -1, 1)
DW.mean$run_b = factor(DW.mean$run_b)
DW.mean$CS_ID_b = ifelse(DW.mean$CS_ID == "val", -1, 1)
DW.mean$CS_ID_b = factor(DW.mean$CS_ID_b)


fit_dw_change<- brm(scale(DW_delta) ~ CS_ID_b * run_b + (CS_ID_b * run_b |ID), 
                    prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                    data =DW.mean, # on the aggregate data
                    iter = 40000, warmup=5000,
                    family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_dw_change, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_pupil_change, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


     
     
# ------------------------------PLOT DEVALUATION IN EYES BEHAVIOR -----------
     
DEV_BEHAV.long <- gather(DEV_BEHAV, eyes_response, Mean_diff, dw_deval:pupil_deval, factor_key=TRUE)
     
     
DEV_BEHAV.long$eyes_response <- recode(DEV_BEHAV.long$eyes_response, pupil_deval = "Pupil",
                                            dw_deval = "Dwell time")
     
     
sumstat.eyes <- summarySE(DEV_BEHAV.long,
                               measurevar = c("Mean_diff"),
                               groupvars = "eyes_response")
     
     
     
eyes.deval.pp <-  ggplot(data = DEV_BEHAV.long, aes (x=eyes_response, y = Mean_diff, fill = eyes_response, color = eyes_response)) +
       geom_abline(slope=0, intercept=0, linetype = "dashed", color = "black") +
       geom_point(aes(color = eyes_response),position = position_jitterdodge(jitter.width = .1, jitter.height = 0),alpha = .3) +
       geom_flat_violin(scale = "count", trim = FALSE, alpha = .2, aes(fill = eyes_response, color = NA), color = NA) +
       geom_crossbar(data = sumstat.eyes, aes(y = Mean_diff, ymin=Mean_diff, ymax=Mean_diff), width = 0.5 , alpha = 0.0) +
       geom_errorbar(data = sumstat.eyes, aes(x = eyes_response, ymin = Mean_diff-ci, ymax = Mean_diff+ci), width = 0.3) +
       coord_flip() +
       theme_bw() +
       scale_fill_manual(values=c("#E69F00","#56B4E9")) + 
       scale_color_manual(values=c("#984C26","#000066")) +
       labs(
         title = '',
         x = '',
         y = "Devaluation effect (mean diff.)"
       ) 

eyes.deval.pp = eyes.deval.pp  + timeline_theme + theme(legend.position="none")  + 
  theme(axis.text.y = element_text(angle=90, vjust= 1, hjust=.5))


   
     
pdf(file.path(figures_path,'Figure_2C.pdf'),width=5,height=8)
print(eyes.deval.pp)
dev.off() 


#-------------------------------------------------------------------------------
# Additional analysis suggested by reviewers
#: Effect of outcome devaluation on RT during the test session                    
#-------------------------------------------------------------------------------

# Format database
TEST = subset(PAVMOD, run == "3")
TEST = subset(PAVMOD, CS_ID != "CSmi")
RT.test.m <- aggregate(US_RT ~ ID*CS_ID , data = TEST, function(x) mean(x, na.rm =T))

RT.test.m$CS_ID <- dplyr::recode(RT.test.m$CS_ID, deval = 'Devalued', val = 'Valued')

RT.test.m.bg = summarySEwithin(RT.test.m,
                               measurevar = c("US_RT"),
                               withinvars = c("CS_ID"),
                               idvar = "ID")


# statistical test frequentistic
mf_rt = formula(US_RT ~ CS_ID  + (CS_ID|ID))
rt_mod = lmer(mf_rt, data = TEST, control = my_control) 
summary(rt_mod)
confint(rt_mod, level = 0.95, method = "Wald") 

# Bayesian analysis
RT.test.m$CS_ID_b = ifelse(RT.test.m$CS_ID == "Valued", -1, 1)

fit_IDval <- brm( scale(US_RT) ~ CS_ID_b + (CS_ID_b|ID), 
                  prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                  data = RT.test.m,
                  iter = 40000, warmup=5000,
                  family = gaussian(), save_pars = save_pars(all = TRUE))

describe_posterior(fit_IDval, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_IDval, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))
