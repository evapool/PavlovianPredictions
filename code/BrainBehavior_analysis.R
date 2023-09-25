################################################################################


#         Neural substrates of parallel devaluation-sensitive and 
#           devaluation-insensitive Pavlovian learning in humans 


#                        brain behavior correlations

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
#                                   SPE ROI
#-------------------------------------------------------------------------------
# aggregate over left and right 
m.spe = aggregate(betas ~ ID + roi + SIDE_PE + ID_PE, data= db.BrainBehav.spe,  function(x) mean(x, na.rm =T))


#model RT_ID corr SPE
mf_spe_id = formula(betas ~ (ID_PE) + (1|ID))
spe_id_mod = lmer(mf_spe_id, data = m.spe, control = my_control) 
summary(spe_id_mod)

confint(spe_id_mod, level = 0.95, method = "Wald") 

# BF
spe_id.data = subset(DB.spe, RT_PE == "ID_PE")

fit_spe_id <- brm(scale(betas) ~ scale(ID_PE) + ( 1|ID), 
                  prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                  data = spe_id.data, # on the aggregate data
                  iter = 40000, warmup=5000, 
                  family = gaussian(), save_pars = save_pars(all = TRUE))

describe_posterior(fit_spe_id , estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_spe_id, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))




#model RT_SIDE corr SPE
mf_spe_side = formula(betas ~ (SIDE_PE) + (1|ID))
spe_side_mod = lmer(mf_spe_side, data =  db.BrainBehav.spe, control = my_control) 
summary(spe_side_mod)
confint(spe_side_mod, level = 0.95, method = "Wald") 


# BF
spe_side.data = subset(DB.spe, RT_PE == "SIDE_PE")

fit_spe_side <- brm(scale(betas) ~ scale(SIDE_PE) + ( 1|ID), 
                    prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                    data =subset(spe_side.data), # on the aggregate data
                    iter = 40000, warmup=5000,
                    family = gaussian(), save_pars = save_pars(all = TRUE))

describe_posterior(fit_spe_side, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_spe_side, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


#-------------------------------------------------------------------------------
#                                   RPE ROI
#-------------------------------------------------------------------------------
# model 1
mf_rpe_id = formula(betas ~ (ID_PE) + (1|ID))
rpe_id_mod = lmer(mf_rpe_id, db.BrainBehav.rpe, control = my_control) 
summary(rpe_id_mod)
confint(rpe_id_mod, level = 0.95, method = "Wald") 


# BF
rpe_id.data = subset(DB.rpe, RT_PE == "ID_PE")

fit_rpe_id <- brm(scale(betas) ~ scale(ID_PE) + ( 1|ID), 
                  prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                  data =subset(rpe_id.data), # on the aggregate data
                  iter = 40000, warmup=5000,
                  family = gaussian(), save_pars = save_pars(all = TRUE))

describe_posterior(fit_rpe_id, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_rpe_id , diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))

# model 2
mf_rpe_side = formula(betas ~ (SIDE_PE) + (1|ID))
rpe_side_mod = lmer(mf_rpe_side, data = db.BrainBehav.rpe, control = my_control) 
summary(rpe_side_mod)
confint(rpe_side_mod, level = 0.95, method = "Wald") 


#BF
rpe_side.data = subset(DB.rpe, RT_PE == "SIDE_PE")

fit_rpe_side <- brm(scale(betas) ~ scale(SIDE_PE) + ( 1|ID), 
                    prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                    data =subset(rpe_side.data), # on the aggregate data
                    iter = 40000, warmup=5000,
                    family = gaussian(), save_pars = save_pars(all = TRUE))

describe_posterior(fit_rpe_side, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_rpe_id , diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))



#-------------------------------------------------------------------------------
#                                   PLOT CORRELATIONS
#-------------------------------------------------------------------------------

# pannel A
DB.spe = subset(DB.all, feature == "SPE")


pp_lin_spe = ggplot(DB.spe, aes(x = RT, y = betas, fill = RT_PE, color = RT_PE)) +
  geom_point(alpha = .5, size = 3.5, position = position_jitterdodge(jitter.width = .0, jitter.height = 0)) +
  geom_smooth(method = lm, level = .95, alpha = .1, fullrange=TRUE) +
  ylab('Betas of SPE (a.u)')+
  xlab('RT [unexpected - expected] (a.u)')+
  scale_fill_manual(name="Outcome Feature",labels=c("Side", "Identity"), values=c("#FFCC33", "#990000")) +
  scale_color_manual(name="Outcome Feature",labels=c("Side", "Identity"), values=c("#FFCC33", "#990000")) +
  ylim(c(-1.5,2)) +
  xlim(c(-0.35, 0.35))+
  theme_bw()

pp_lin_spe = pp_lin_spe + timeline_theme 


pp_lin_spe <- ggMarginal(pp_lin_spe + theme(legend.position = "bottom"),type = "density", 
                         groupColour = TRUE, groupFill = T, alpha = 0.2)


pdf(file.path(figures_path,'Figure_4A.pdf'))
print(pp_lin_spe)
dev.off()


# pannel B
DB.rpe = subset(DB.all, feature == "RPE")

pp_lin_rpe = ggplot(DB.rpe, aes(x = RT, y = betas, fill = RT_PE, color = RT_PE)) +
  geom_point(alpha = .5, size = 3.5, position = position_jitterdodge(jitter.width = .0, jitter.height = 0)) +
  geom_smooth(method = lm, level = .95, alpha = .1, fullrange=TRUE) +
  ylab('Betas of RPE (a.u)')+
  xlab('RT [unexpected - expected] (a.u)')+
  scale_fill_manual(name="Outcome Feature",labels=c("Side", "Identity"), values=c("#FFCC33", "#990000")) +
  scale_color_manual(name="Outcome Feature",labels=c("Side", "Identity"), values=c("#FFCC33", "#990000")) +
  ylim(c(-1.5,2)) +
  xlim(c(-0.35, 0.35))+
  theme_bw()

pp_lin_rpe = pp_lin_rpe  + timeline_theme 

pp_lin_rpe  <- ggMarginal(pp_lin_rpe  + theme(legend.position = "bottom"),type = "density", 
                          groupColour = F, yparams = list(fill = "#56B4E9", color = "#000066",alpha=.2))


pdf(file.path(figures_path,'Figure_4B.pdf'))
print(pp_lin_rpe)
dev.off()




