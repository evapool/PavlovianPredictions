###############################################################################


#         Neural substrates of parallel devaluation-sensitive and 
#           devaluation-insensitive Pavlovian learning in humans 


#                        fMRI analysis

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





#------------------------------------------------------------------------------
#          Parallel Pavlovian predictions about affective value and
#                        perceptual attributes of the outcome
#-------------------------------------------------------------------------------


# ------------------- REWARD PREDICTION ERROR DISTRIBUTION --------------------

# --------------------- MIDBRAIN
t.test(LEARN_RPE$RPE_in_RPE_VTA_ROI); se(LEARN_RPE$RPE_in_RPE_VTA_ROI)
# BF
ttestBF(LEARN_RPE$RPE_in_RPE_VTA_ROI)


# -----------------------VS and Subcallosal
t.test(LEARN_RPE$RPE_in_RPE_VS_ROI); se(LEARN_RPE$RPE_in_RPE_VS_ROI)
# BF
ttestBF(LEARN_RPE$RPE_in_RPE_VS_ROI)


# ----------------------- vmPFC
t.test(LEARN_RPE$RPE_in_RPE_VMPF_ROI); se(LEARN_RPE$RPE_in_RPE_VMPF_ROI)
# BF
ttestBF(LEARN_RPE$RPE_in_RPE_VMPF_ROI)




# ------------------- REWARD PREDICTION ERROR DEVALUATION--------------------

#----------------------- Create variables for bayse factor computation
CHANGE_RPE_long$value_b = ifelse(CHANGE_RPE_long$value  == "devalued", -1, 1)
CHANGE_RPE_long$value_b = factor(CHANGE_RPE_long$value_b)
CHANGE_RPE_long$roi_b   = as.character(CHANGE_RPE_long$roi)


#-------------------------------------------------------------- VS subcallosal
vs.mdl = lme(betas ~ value,  data= subset(CHANGE_RPE_long, roi == "RPE_VS"), random= ~ value|ID)
summary(vs.mdl)
intervals(vs.mdl, which = "fixed")

# BF
fit_vs_dev<- brm( scale(betas) ~ value_b + (value_b |ID), 
                  prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                  data =subset(CHANGE_RPE_long, roi == "RPE_VS"), # on the aggregate data
                  iter = 40000, warmup=5000,
                  family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_vs_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_vs_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


#-------------------------------------------------------------- mid brain
vta.mdl = lme(betas ~ value,  data= subset(CHANGE_RPE_long, roi == "RPE_VTA"), random= ~ value|ID)
summary(vta.mdl)
intervals(vta.mdl, which = "fixed")

# BF
fit_vta_dev<- brm( scale(betas) ~ value_b + (value_b |ID), 
                   prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                   data =subset(CHANGE_RPE_long, roi == "RPE_VTA"), # on the aggregate data
                   iter = 40000, warmup=5000,
                   family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_vta_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_vta_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))



#------------------------- VMPFC
vpmfc.mdl = lme(betas ~ value,  data= subset(CHANGE_RPE_long, roi == "RPE_VMPF"), random= ~ value|ID)
summary(vpmfc.mdl)
intervals(vpmfc.mdl, which = "fixed")

# BF
fit_vmpf_dev <- brm( scale(betas) ~ value_b + (value_b |ID), 
                     prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                     data =subset(CHANGE_RPE_long, roi == "RPE_VMPF"), # on the aggregate data
                     iter = 40000, warmup=5000,
                     family = gaussian(), save_pars = save_pars(all = TRUE))

describe_posterior(fit_vmpf_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_vmpf_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))





# ------------------- STATE PREDICTION ERROR DISTRUBITION------------------------

# --------------------- OFC left
t.test(LEARN_SPE$SPE_in_SPE_left_OFC_ROI); se(LEARN_SPE$SPE_in_SPE_left_OFC_ROI)
# BF
ttestBF(LEARN_SPE$SPE_in_SPE_left_OFC_ROI)


# --------------------- OFC right
t.test(LEARN_SPE$SPE_in_SPE_right_OFC_ROI); se(LEARN_SPE$SPE_in_SPE_left_OFC_ROI)
# BF
ttestBF(LEARN_SPE$SPE_in_SPE_right_OFC_ROI)


# --------------------- MFG right
t.test(LEARN_SPE$SPE_in_SPE_right_MFG_ROI); se(LEARN_SPE$SPE_in_SPE_right_MFG_ROI)
# BF
ttestBF(LEARN_SPE$SPE_in_SPE_right_MFG_ROI)



# --------------------- SFG
t.test(LEARN_SPE$SPE_in_SPE_SFG_ROI); se(LEARN_SPE$SPE_in_SPE_SFG_ROI)
# BF
ttestBF(LEARN_SPE$SPE_in_SPE_SFG_ROI)


# --------------------- Midbrain

t.test(LEARN_SPE$SPE_in_SPE_VTA_ROI); se(LEARN_SPE$SPE_in_SPE_VTA_ROI)
# BF
ttestBF(LEARN_SPE$SPE_in_SPE_VTA_ROI)



# ------------------- STATE PREDICTION ERROR DEVALUATION------------------------


#-------------------------------------------- variables to compute bayse factors
CHANGE_SPE_long$value_b = ifelse(CHANGE_SPE_long$value  == "devalued", -1, 1)
CHANGE_SPE_long$value_b= factor(CHANGE_SPE_long$value_b)

CHANGE_SPE_long$side_b = ifelse(CHANGE_SPE_long$side == "right", -1, 1)
CHANGE_SPE_long$side_b= factor(CHANGE_SPE_long$side_b)


#--------------------------------------------------------------  MFG
mfg.mdl = lme(betas ~ value,  data= subset(CHANGE_SPE_long, roi == "SPE_MFG" & side == "right"), random= ~ value|ID)
summary(mfg.mdl)
intervals(mfg.mdl, which = "fixed")

#BF
fit_mfg_dev<- brm( scale(betas) ~ value_b + (value_b |ID), 
                   prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                   data =subset(CHANGE_SPE_long, roi == "SPE_MFG" & side == "right"), # on the aggregate data
                   iter = 40000, warmup=5000,                  
                   family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_mfg_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_mfg_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))



#-------------------------------------------------------------- SFG
sfg.mdl = lme(betas ~ value,  data= subset(CHANGE_SPE_long, roi == "SPE_SFG"), random= ~ value|ID)
summary(sfg.mdl)
intervals(sfg.mdl, which = "fixed")

#BF
fit_sfg_dev<- brm( scale(betas) ~ value_b + (value_b |ID), 
                   prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                   data =subset(CHANGE_SPE_long, roi == "SPE_SFG"), # on the aggregate data
                   iter = 40000, warmup=5000,                  
                   family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_sfg_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_sfg_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


#--------------------------------------------------------------  OFC
ofc.mdl = lme(betas ~ value*side,  data= subset(CHANGE_SPE_long, roi == "SPE_OFC"), random= ~ value/side|ID)
summary(ofc.mdl)
intervals(ofc.mdl, which = "fixed")


#BF
fit_ofc_dev<- brm( scale(betas) ~ value_b*side_b + (value_b*side_b |ID), 
                   prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                   data =subset(CHANGE_SPE_long, roi == "SPE_OFC"), # on the aggregate data
                   iter = 40000, warmup=5000, 
                   family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_ofc_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_ofc_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


#--------------------------------------------------------------  VTA
vta.mdl = lme(betas ~ value,  data= subset(CHANGE_SPE_long, roi == "SPE_VTA"), random= ~ value|ID)
summary(vta.mdl)
intervals(vta.mdl, which = "fixed")


#BF
fit_vta_dev<- brm( scale(betas) ~ value_b + (value_b |ID), 
                   prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                   data =subset(CHANGE_SPE_long, roi == "SPE_VTA"), # on the aggregate data
                   iter = 40000, warmup=5000,  
                   family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_vta_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_vta_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))







#-------------------------------------------------------------------------------
#                             PLOT FIGURE 3
#-------------------------------------------------------------------------------

#--------------------- PLOT RPE ------------------------------------------------
rpe.deval.pp <-  ggplot(data = rpe.bs, aes (x=roi, y = diff_betas, fill = roi, color = roi)) +
  geom_abline(slope=0, intercept=0, linetype = "dashed", color = "black") +
  geom_point(aes(color = roi),position = position_jitterdodge(jitter.width = .1, jitter.height = 0),alpha = .1) +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .2, aes(fill = roi, color = NA), color = NA) +
  geom_crossbar(data = sumstat.rpe, aes(y = diff_betas, ymin=diff_betas, ymax=diff_betas), width = 0.5 , alpha = 0.0) +
  geom_errorbar(data = sumstat.rpe, aes(x = roi, ymin = diff_betas-ci, ymax = diff_betas+ci), width = 0.3) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous( limits = c(-5,5)) +
  scale_fill_manual(values=c("#56B4E9","#56B4E9","#56B4E9")) + 
  scale_color_manual(values=c("#0e3c5a","#0e3c5a","#0e3c5a")) +
  labs(
    title = '',
    x = 'ROI based on Reward PE',
    y = "Betas [valued - devalued] (arb. units)"
  ) 


rpe.deval.pp = rpe.deval.pp + timeline_theme + theme(legend.position="none") 


pdf(file.path(figures_path,'Figure_3B_RPE.pdf'))
print(rpe.deval.pp)
dev.off() 

# write sourcefile
rpe.bs.source = rpe.bs[c('ID','roi', 'feature', 'diff_betas')]
write.csv(rpe.bs.source,file.path(figures_path, 'SourceData_3B_RPE_individual_estimates.csv'))

sumstat.rpe$mean = sumstat.rpe$diff_betas
rpe.bs.aggregate = sumstat.rpe[c('roi','N', 'mean', 'ci')]
write.csv(rpe.bs.aggregate,file.path(figures_path, 'SourceData_3B_RPE_aggregated_estimates.csv'))



#--------------------- PLOT SPE -------------------------------------------------
spe.deval.pp <- ggplot(data = spe.bs, aes (x=roi, y = diff_betas, fill = roi, color = roi)) +
  geom_abline(slope=0, intercept=0, linetype = "dashed", color = "black") +
  geom_point(aes(color = roi),position = position_jitterdodge(jitter.width = .1, jitter.height = 0),alpha = .1) +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .2, aes(fill = roi, color = NA), color = NA) +
  geom_crossbar(data = sumstat.spe, aes(y = diff_betas, ymin=diff_betas, ymax=diff_betas), width = 0.5 , alpha = 0.0) +
  geom_errorbar(data = sumstat.spe, aes(x = roi, ymin = diff_betas-ci, ymax = diff_betas+ci), width = 0.3) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous( limits = c(-5,5)) +
  scale_fill_manual(values=c("#E69F00","#E69F00","#E69F00","#E69F00")) + 
  scale_color_manual(values=c("#9D6C00","#9D6C00","#9D6C00","#9D6C00")) +
  theme_bw() +
  labs(
    title = '',
    x = 'ROI based on State PE',
    y = "Betas [valued - devalued] (arb. units)"
  ) 

spe.deval.pp = spe.deval.pp + timeline_theme + theme(legend.position="none") 


pdf(file.path(figures_path,'Figure_3B_SPE.pdf'))
print(spe.deval.pp)
dev.off() 


# write sourcefile
spe.bs.source = spe.bs[c('ID','roi', 'feature', 'diff_betas')]
write.csv(spe.bs.source,file.path(figures_path, 'SourceData_3B_SPE_individual_estimates.csv'))

sumstat.spe$mean = sumstat.spe$diff_betas
spe.bs.aggregate = sumstat.spe[c('roi','N', 'mean', 'ci')]
write.csv(spe.bs.aggregate,file.path(figures_path, 'SourceData_3B_SPE_aggregated_estimates.csv'))





#--------------------------------------------------------------------------------
#            Pavlovian predictions about spatial location and taste
#                      identity attributes of the outcome                   
#-------------------------------------------------------------------------------

#----------------------- DISTRIBUTION  ID  --------------------------------------


# --------------------- IFG

# MEAN ACC AND CI
t.test(LEARN_ID_ifg$IFG_R, mu = 0.5); se(LEARN_ID_ifg$IFG_R) 
# BF
ttestBF(LEARN_ID_ifg$IFG_R, mu = 0.5)

# --------------------- IPS

# MEAN ACC AND CI
t.test(LEARN_ID_ips$IPS_R, mu = 0.5); se(LEARN_ID_ips$IPS_R) 
# BF
ttestBF(LEARN_ID_ips$IPS_R, mu = 0.5)


# --------------------- PCL_L

# MEAN ACC AND CI
t.test(LEARN_ID_pcl$PCL_L, mu = 0.5); se(LEARN_ID_pcl$PCL_L) 
# BF
ttestBF(LEARN_ID_pcl$PCL_L, mu = 0.5)

# --------------------- PCG_L

# MEAN ACC AND CI
t.test(LEARN_ID_pcg$PCG_L, mu = 0.5); se(LEARN_ID_pcg$PCG_L) 
# BF
ttestBF(LEARN_ID_pcg$PCG_L, mu = 0.5)




#----------------------- DEVALUATION  ID  --------------------------------------


#-------------------------------------------- get variables to compute BF 
ID_long$CS_b = ifelse(ID_long$CS == "devalued", -1, 1)
ID_long$CS_b = factor(ID_long$CS_b)

ID_long = subset(ID_long, prepost == "post")

#------------------------------------------------ in each ROI: L_PCG 
pcg_dev_mod = lme(index ~ CS,  data= subset(ID_long, roi == "L_PCG"), random= ~ CS|ID)
summary(pcg_dev_mod )
intervals(pcg_dev_mod, which = "fixed" )


#BF
fit_lpcg_dev<- brm( scale(index) ~ CS_b + ( CS_b|ID), 
                    prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                    data = subset(ID_long, roi == "L_PCG"), # on the aggregate data
                    iter = 40000, warmup=5000, 
                    family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_lpcg_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_lpcg_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


#------------------------------------------------  in each ROI L_PCL 
pcl_dev_mod = lme(index ~ CS,  data= subset(ID_long, roi == "L_PCL"), random= ~ CS|ID)
summary(pcl_dev_mod )
intervals(pcl_dev_mod, which = "fixed")


#BF
fit_pcl_dev<- brm( scale(index) ~ CS_b + ( CS_b|ID), 
                   prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                   data = subset(ID_long, roi == "L_PCL"), # on the aggregate data
                   iter = 40000, warmup=5000, 
                   family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_pcl_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_pcl_dev,
                   diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))



#------------------------------------------------  in each ROI R_IPS
ips_dev_mod = lme(index ~ CS,  data= subset(ID_long, roi == "R_IPS"), random= ~ CS|ID)
summary(ips_dev_mod )
intervals(ips_dev_mod, which = "fixed")


#BF
fit_ips_dev<- brm(scale(index) ~ CS_b + ( CS_b|ID), 
                  prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                  data = subset(ID_long, roi == "R_IPS"), # on the aggregate data
                  iter = 40000, warmup=5000, 
                  family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_ips_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_ips_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))



#------------------------------------------------  in each ROI R_IFG
ifg_dev_mod = lme(index ~ CS,  data= subset(ID_long, roi == "R_IFG"), random= ~ CS|ID)
summary(ifg_dev_mod )
intervals(ifg_dev_mod, which = "fixed")


#BF
fit_ifg_dev<- brm(scale(index) ~ CS_b + ( CS_b|ID), 
                  prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                  data = subset(ID_long, roi == "R_IFG"), # on the aggregate data
                  iter = 40000, warmup=5000, 
                  family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_ifg_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_ifg_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))




#---------------------- DISTRUBITION SIDE --------------------------------------

# --------------------- cuneus
# MEAN ACC AND CI
t.test(LEARN_SIDE_cuneus$Cuneus, mu = 0.5); se(LEARN_SIDE_cuneus$Cuneus) 
# BF
ttestBF(LEARN_SIDE_cuneus$Cuneus, mu = 0.5)

# --------------------- ips
# MEAN ACC AND CI
t.test(LEARN_SIDE_ips$R_IPS, mu = 0.5); se(LEARN_SIDE_ips$R_IPS) 
# BF
ttestBF(LEARN_SIDE_ips$R_IPS, mu = 0.5)


# --------------------- stg  / smg l
# MEAN ACC AND CI
t.test(LEARN_SIDE_smg_l$STG_L, mu = 0.5); se(LEARN_SIDE_smg_l$STG_L) 
# BF
ttestBF(LEARN_SIDE_smg_l$STG_L, mu = 0.5)

# --------------------- stg / smg r
# MEAN ACC AND CI
t.test(LEARN_SIDE_smg_r$STG_R, mu = 0.5); se(LEARN_SIDE_smg_r$STG_R) 
# BF
ttestBF(LEARN_SIDE_smg_r$STG_R, mu = 0.5)



# --------------------- lat occ
# MEAN ACC AND CI
t.test(LEARN_SIDE_latocc$R_LATOCC, mu = 0.5); se(LEARN_SIDE_latocc$R_LATOCC) 
# BF
ttestBF(LEARN_SIDE_latocc$R_LATOCC, mu = 0.5)




#---------------------- DEVALUATION SIDE ---------------------------------------

#-------------------------------------------- compute variable for BF 
SIDE_long_m$CS_b = ifelse(SIDE_long_m$CS == "devalued", -1, 1)
SIDE_long_m$CS_b = factor(SIDE_long_m$CS_b)
SIDE_long_m = subset(SIDE_long_m, prepost == "post")


#------------------------------------------------ in each ROI: Cuneus 
Cuneus_dev_mod = lme(index ~ CS,  data= subset(SIDE_long_m, roi == "Cuneus"), random= ~ CS|ID)
summary(Cuneus_dev_mod  )
intervals(Cuneus_dev_mod, which = "fixed" )


#BF
fit_cuneus_dev<- brm(scale(index) ~ CS_b + ( CS_b|ID), 
                     prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                     data= subset(SIDE_long_m, roi == "Cuneus"), # on the aggregate data
                     iter = 40000, warmup=5000, 
                     family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_cuneus_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_cuneus_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


##------------------------------------------------ in each ROI: IPS 
ips_dev_mod = lme(index ~ CS,  data= subset(SIDE_long_m, roi == "IPS"), random= ~ CS|ID)
summary(ips_dev_mod)
intervals(ips_dev_mod, which = "fixed")


#BF
fit_ips_dev <- brm(scale(index) ~ CS_b + ( CS_b|ID), 
                   prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                   data= subset(SIDE_long_m, roi == "IPS"), # on the aggregate data
                   iter = 40000, warmup=5000, 
                   family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_ips_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_ips_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


##------------------------------------------------ in each ROI: LACTOCC 
lactocc_dev_mod  = lme(index ~ CS,  data= subset(SIDE_long_m, roi == "LACTOCC"), random= ~ CS|ID)
summary(lactocc_dev_mod) 
intervals(lactocc_dev_mod, which = "fixed")


#BF
fit_latocc_dev <- brm(scale(index) ~ CS_b + ( CS_b|ID), 
                      prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                      data= subset(SIDE_long_m, roi == "LACTOCC"), # on the aggregate data
                      iter = 40000, warmup=5000, 
                      family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(fit_latocc_dev, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_latocc_dev, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))



##------------------------------------------------ in each ROI: STG /SMG
smg_dev_mod = lme(index ~ CS,  data= subset(SIDE_long_m, roi == "SMG"), random= ~ CS|ID)
summary(smg_dev_mod) 
intervals(smg_dev_mod, which = "fixed")


#BF
smg_dev_mod  <- brm(scale(index) ~ CS_b + ( CS_b|ID), 
                    prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                    data= subset(SIDE_long_m, roi == "SMG"), # on the aggregate data
                    iter = 40000, warmup=5000, 
                    family = gaussian(), save_pars = save_pars(all = TRUE))


describe_posterior(smg_dev_mod , estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = smg_dev_mod, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))




#-------------------------------------------------------------------------------
#                             PLOT FIGURE 5
#-------------------------------------------------------------------------------

#--------------------------- PLOT ID -------------------------------------------
id.bs = ddply(ID_long,.(ID,roi),transform, index2=mean(index2))
id.bs = subset (id.bs, CS == "valued")
id.bs = subset (id.bs, prepost == "post")

id.bs$roi <- dplyr::recode(id.bs$roi, "L_PCG" = "PCG", "L_PCL" = "PCL", "R_IFG" = "IFG", "R_IPS" = "IPS")
id.bs$roi <- as.character((id.bs$roi))

id.bg = ddply(id.bs,.(roi),summarise, index2=mean(index2, na.rm = T))



sumstat.id <- summarySE(id.bs,
                        measurevar = c("index2"),
                        groupvars = "roi")


id.deval.pp <-  ggplot(data = id.bs, aes (x=roi, y = index2, fill = roi, color = roi)) +
  geom_abline(slope=0, intercept=0, linetype = "dashed", color = "black") +
  geom_point(aes(fill = NA, color = roi),position = position_jitterdodge(jitter.width = .1, jitter.height = 0),alpha = .2) +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .3, aes(fill = roi, color = NA), color = NA) +
  geom_crossbar(data = sumstat.id, aes(y = index2, ymin=index2, ymax=index2), width = 0.5 , alpha = 0.0) +
  geom_errorbar(data = sumstat.id, aes(x = roi, ymin = index2-ci, ymax = index2+ci), width = 0.3) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous( limits = c(-2,2)) +
  scale_fill_manual(values=c("#CC0000","#CC0000","#CC0000","#CC0000")) + 
  scale_color_manual(values=c("#660000","#660000","#660000","#660000")) +
  labs(
    title = '',
    x = 'Voxles decoding Identity',
    y = "Betas [mean diff.] "
  ) 


id.deval.pp = id.deval.pp + timeline_theme + theme(legend.position="none")  


pdf(file.path(figures_path,'Figure_5B_identity.pdf'))
print(id.deval.pp)
dev.off() 

# write sourcefile
id.bs$mean_diff = id.bs$index2
id.bs.source = id.bs[c('ID','roi', 'mean_diff')]
write.csv(id.bs.source,file.path(figures_path, 'SourceData_5B_identity_individual_estimates.csv'))


sumstat.id$mean = sumstat.id$index2 
id.aggregate = sumstat.id[c('roi','N', 'mean', 'ci')]
write.csv(id.aggregate,file.path(figures_path, 'SourceData_5B_identity_aggregated_estimates.csv'))




#--------------------------- PLOT SIDE -------------------------------------------------
side.bs = ddply(SIDE_long_m_index,.(ID,roi),transform, index2=mean(index2))
side.bs = subset (side.bs, CS == "valued")
side.bs = subset (side.bs, prepost == "post")

side.bs$roi <- dplyr::recode(side.bs$roi, "LACTOCC" = "LOC")
side.bs$roi <- as.character((side.bs$roi))

side.bg = ddply(side.bs,.(roi),summarise, index2=mean(index2, na.rm = T))


sumstat.side <- summarySE(side.bs,
                          measurevar = c("index2"),
                          groupvars = "roi")

side.deval.pp <-  ggplot(data = side.bs, aes (x=roi, y = index2, fill = roi, color = roi)) +
  geom_abline(slope=0, intercept=0, linetype = "dashed", color = "black") +
  geom_point(aes(fill = NA, color = roi),position = position_jitterdodge(jitter.width = .1, jitter.height = 0),alpha = .2) +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .3, aes(fill = roi, color = NA), color = NA) +
  geom_crossbar(data = sumstat.side, aes(y = index2, ymin=index2, ymax=index2), width = 0.5 , alpha = 0.0) +
  geom_errorbar(data = sumstat.side, aes(x = roi, ymin = index2-ci, ymax = index2+ci), width = 0.3) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous( limits = c(-2,2)) +
  scale_fill_manual(values=c("#FFCC33","#FFCC33","#FFCC33","#FFCC33")) + 
  scale_color_manual(values=c("#9D6C00","#9D6C00","#9D6C00","#9D6C00")) +
  labs(
    title = '',
    x = 'Voxels decoding Side',
    y = "Betas [mean diff.]"
  ) 


side.deval.pp = side.deval.pp + timeline_theme + theme(legend.position="none") 


pdf(file.path(figures_path,'Figure_5B_side.pdf'))
print(side.deval.pp)
dev.off()

# write sourcefile
side.bs$mean_diff =  side.bs$index2
side.bs.source =  side.bs[c('ID','roi', 'mean_diff')]
write.csv(side.bs.source,file.path(figures_path, 'SourceData_5B_side_individual_estimates.csv'))


sumstat.side$mean = sumstat.side$index2 
side.aggregate = sumstat.side[c('roi','N', 'mean', 'ci')]
write.csv(side.aggregate,file.path(figures_path, 'SourceData_5B_side_aggregated_estimates.csv'))



