###############################################################################


#         Neural substrates of parallel devaluation-sensitive and 
#           devaluation-insensitive Pavlovian learning in humans 


#                        analysis for supplementary materials

# last  modified on AUGUST 2023
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

# set up specific to revisions

pannel_theme <- theme(strip.text.x = element_text(size = 12, face = "bold"),
                      strip.background = element_rect(color="white", fill="white", linetype="solid"),
                      panel.grid.major.x = element_blank() ,
                      panel.grid.major.y = element_line(size=.2, color="lightgrey"),
                      panel.grid.minor = element_blank(),
                      legend.text  = element_text(size =  14),
                      legend.title = element_text(size =  14),
                      legend.position = "bottom",
                      axis.title.x =element_blank(),
                      axis.ticks.x=element_blank(), 
                      axis.text.x=element_blank(),
                      axis.title.y = element_text(size =  14),
                      axis.line = element_line(size = 0.5))



#-------------------------------------------------------------------------------
#                  REVIEWER 1
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# R1 Comment 2 :  Plot valued and devalued conditions separately                  
#-------------------------------------------------------------------------------


# ------------------------------------ Pupil -----------------------------------

#format database
db.pupil.change <- subset(PUPIL.mean, run == 2)

db.pupil.change$CS_ID <- dplyr::recode(db.pupil.change$CS_ID, deval = 'Devalued', 
                                       val = 'Valued')
db.pupil.change.bc = summarySEwithin(db.pupil.change,
                                     measurevar = c("pupil_change"),
                                     withinvars = c("CS_ID"),
                                     idvar = "ID")
# plot
pupil.deval.pp <-  ggplot(data = db.pupil.change, 
                          aes (x=CS_ID, y = pupil_change, fill = CS_ID, color = CS_ID)) +
  geom_abline(slope=0, intercept=0, linetype = "dashed", color = "black") +
  geom_point(aes(x = CS_ID, y = pupil_change, group = ID), 
             alpha = 0.3, position = position_dodge(0.6)) +
  geom_flat_violin(aes(x = 2.7), alpha = .5, position = position_nudge(x = -.15, y = 0), 
                   adjust = 1.8, trim = F, color = NA) +
  geom_crossbar(data = db.pupil.change.bc, 
                aes(x = CS_ID, y = pupil_change, ymin=pupil_change, ymax=pupil_change),
                width = 0.8) +
  geom_errorbar(data = db.pupil.change.bc, 
                aes(x = CS_ID, ymin = pupil_change-se, ymax = pupil_change+se), 
                width = 0.7) +
  theme_bw() +
  scale_fill_manual(values=c("black","#56B4E9")) + 
  scale_color_manual(values=c("black","#000066")) +
  annotate("segment", x = 1, xend = 2, y = 0.3, yend = 0.3) +
  annotate("text", x = 1.5, y = 0.33, label = "**", size = 6) +
  labs( title = '', x = 'CS', y = "Pupil (after - before)") 


pupil.deval.pp = pupil.deval.pp  + timeline_theme + theme(legend.position="none")  + 
  theme(axis.text.y = element_text(angle=90, vjust= 1, hjust=.5))

pdf(file.path(figures_path,'Supplementary_Figure_%_A.pdf'),width=5,height=8)
print(pupil.deval.pp )
dev.off() 


# ----------------------------------- Dwell Time -------------------------------

# format database
db.dw.change <- subset(DW.mean, run == 2)

db.dw.change$CS_ID <- dplyr::recode(db.dw.change$CS_ID, deval = 'Devalued', 
                                    val = 'Valued')

db.dw.change.bc = summarySEwithin(db.dw.change,
                                  measurevar = "dw_change",
                                  withinvar = "CS_ID",
                                  idvar = "ID")

# plot
dw.deval.pp <- ggplot(data =db.dw.change, 
                      aes(x = CS_ID, y = dw_change , fill = CS_ID, color = CS_ID)) +
  geom_abline(slope = 0, intercept = 0, linetype= "dashed", color="black") +
  geom_point(aes(x = CS_ID, y = dw_change, group = ID), 
             alpha = 0.3, position = position_dodge(0.6)) +
  geom_flat_violin(aes(x = 2.7), alpha = .5, position = position_nudge(x = -.15, y = 0), 
                   adjust = 1.8, trim = F, color = NA) +
  geom_crossbar(data = db.dw.change.bc, 
                aes(x = CS_ID, y = dw_change, ymin=dw_change, ymax=dw_change),
                width = 0.8, alpha = 0.0) +
  geom_errorbar(data = db.dw.change.bc, 
                aes(x = CS_ID, ymin = dw_change-se, ymax = dw_change+se), 
                width = 0.7) +
  theme_bw() +
  scale_fill_manual(values=c("black","#E69F00")) + 
  scale_color_manual(values=c("black","#984C26")) +
  labs( title = '',x = 'CS',y = "Dwell time (after - before)") 

#make it nice
dw.deval.pp =  dw.deval.pp + timeline_theme + theme(legend.position="none")  + 
  theme(axis.text.y = element_text(angle=90, vjust= 1, hjust=.5))

#print
pdf(file.path(figures_path,'Supplementary_Figure_%_B.pdf'),width=5,height=8)
print(dw.deval.pp )
dev.off()


# ----------------------------------- RPE ROI ----------------------------------

# Format database
supp.rpe.bs = ddply(CHANGE_RPE_long,.(ID,roi),transform, diff_betas=mean(betas[value == "valued"] - mean(betas[value == "devalued"])))
supp.rpe.bg = ddply(supp.rpe.bs,.(roi,value),summarise, betas=mean(betas, na.rm = T))

supp.rpe.bs$roi <- dplyr::recode(supp.rpe.bs$roi, RPE_VTA = 'midbrain', 
                                 RPE_VMPF = 'vmPFC',
                                 RPE_VS = 'VS')

supp.sumstat.rpe  <- summarySEwithin(supp.rpe.bs,
                                     measurevar = c("betas"),
                                     withinvar = c("roi","value"),
                                     idvar = "ID")

# plot
pp.rpe.roi.deval <- ggplot(data = supp.rpe.bs, 
                           aes(x = value, y = betas, fill = value, color = value)) +
  geom_abline(slope = 0, intercept = 0, linetype= "dashed", color="black") +
  geom_point(aes(x = value, y = betas, group = ID), 
             alpha = 0.3, position = position_dodge(0.6)) +
  geom_flat_violin(aes(x = 2.7), alpha = .5, position = position_nudge(x = -.15, y = 0), 
                   adjust = 1.8, trim = F, color = NA) +
  geom_crossbar(data = supp.sumstat.rpe, aes(x = value, y = betas, ymin=betas, ymax=betas),
                width = 0.8 ) +
  geom_errorbar(data = supp.sumstat.rpe, aes(x = value, y = betas, ymin = betas-se, ymax = betas+se),
                width = 0.7) +
  theme_bw() +
  facet_wrap(~roi, nrow = 1) +
  scale_y_continuous( limits = c(-6.5,3)) +
  scale_fill_manual(values=c("black","#56B4E9")) + 
  scale_color_manual(values=c("black","#000066")) +
  labs(title = '', x = 'CS', y = "Betas" ) 


pp.rpe.roi.deval <- pp.rpe.roi.deval + theme_bw( ) + pannel_theme

pdf(file.path(figures_path,'Supplementary_Figure_6.pdf'),width=8,height=5)
print(pp.rpe.roi.deval)
dev.off()



#------------------------------------ SPE ROI ----------------------------------

# format database
supp.spe.bs = ddply(CHANGE_SPE_long,.(ID,roi),transform, diff_betas=mean(betas[value == "valued"] - mean(betas[value == "devalued"])))

supp.spe.bg = ddply(supp.spe.bs,.(roi,value),summarise, betas=mean(betas, na.rm = T))

supp.spe.bs = subset (supp.spe.bs, side != "left") # let remove repetitions


supp.spe.bs$roi <- dplyr::recode(supp.spe.bs$roi, SPE_MFG = 'MFG', 
                                 SPE_OFC = 'OFC',
                                 SPE_SFG = 'SFG',
                                 SPE_VTA = 'midbrain')


supp.sumstat.spe  <- summarySEwithin(supp.spe.bs,
                                     measurevar = c("betas"),
                                     withinvar = c("roi","value"),
                                     idvar = "ID")

#plot
pp.spe.roi.deval <- ggplot(data = supp.spe.bs, aes(x = value, y = betas, fill = value, color = value)) +
  geom_abline(slope = 0, intercept = 0, linetype= "dashed", color="black") +
  geom_point(aes(x = value,
                 y = betas, group = ID), 
             alpha = 0.3, position = position_dodge(0.6)) +
  geom_flat_violin(aes(x = 2.7), alpha = .5, position = position_nudge(x = -.15, y = 0), 
                   adjust = 1.8, trim = F, color = NA) +
  geom_crossbar(data = supp.sumstat.spe, aes(x = value, y = betas, ymin=betas, ymax=betas), width = 0.8 , alpha = 0.0) +
  geom_errorbar(data = supp.sumstat.spe, aes(x = value, y = betas, ymin = betas-se, ymax = betas+se), width = 0.7) +
  theme_bw() +
  facet_wrap(~roi, nrow = 1) +
  scale_y_continuous( limits = c(-6.5,3)) +
  scale_fill_manual(values=c("black","#E69F00")) + 
  scale_color_manual(values=c("black","#984C26")) +
  labs(title = '', x = 'CS', y = "Betas") 

# make it nice
pp.spe.roi.deval <- pp.spe.roi.deval + theme_bw( ) + pannel_theme
# prints
pdf(file.path(figures_path,'Supplementary_Figure_7.pdf'),width=8,height=5)
print(pp.spe.roi.deval)
dev.off()



# ----------------------------------- ID ROI --------------------------------------

# format database
id.bs = subset (ID_long, prepost == "post")

id.bs$roi <- dplyr::recode(id.bs$roi, "L_PCG" = "PCG", "L_PCL" = "PCL", "R_IFG" = "IFG", "R_IPS" = "IPS")

id.bg = ddply(id.bs,.(roi,CS),summarise, index=mean(index, na.rm = T))

sumstat.id <- summarySEwithin(id.bs,
                              measurevar = c("index"),
                              withinvar = c("roi","CS"),
                              idvar = "ID")

#plot
pp.id.roi.deval <- ggplot(data = id.bs, aes(x = CS, y = index, fill = CS, color = CS)) +
  geom_abline(slope = 0, intercept = 0, linetype= "dashed", color="black") +
  geom_point(aes(x = CS,
                 y = index, group = ID), 
             alpha = 0.3, position = position_dodge(0.6)) +
  geom_flat_violin(aes(x = 2.7), alpha = .5, position = position_nudge(x = -.15, y = 0), 
                   adjust = 1.8, trim = F, color = NA) +
  geom_crossbar(data = sumstat.id , aes(x = CS, y = index, ymin=index, ymax=index), width = 0.8 , alpha = 0.0) +
  geom_errorbar(data = sumstat.id , aes(x = CS, y = index, ymin = index-se, ymax = index+se), width = 0.7) +
  theme_bw() +
  facet_wrap(~roi, nrow = 1) +
  scale_y_continuous( limits = c(-6.5,4.5)) +
  scale_fill_manual(values=c("black","#CC0000")) + 
  scale_color_manual(values=c("black","#660000")) +
  labs(title = '', x = 'CS', y = "Betas (after - before)") 

#make it nice
pp.id.roi.deval <- pp.id.roi.deval + theme_bw( ) + pannel_theme

#print
pdf(file.path(figures_path,'Supplementary_Figure_8.pdf'),width=8,height=5)
print(pp.id.roi.deval)
dev.off()


# ---------------------------------- SIDE ROI ---------------------------
side.bs = subset (SIDE_long, prepost == "post")
side.bs.m =  aggregate(index ~ID + roi + CS + prepost, data =side.bs  , function(x) mean(x, na.rm =T))


side.bs.m$roi <- dplyr::recode(side.bs.m$roi, "LACTOCC" = "LOC")

side.bg = ddply(side.bs,.(roi,CS),summarise, index=mean(index, na.rm = T))

sumstat.side <- summarySEwithin(side.bs.m,
                                measurevar = c("index"),
                                withinvar = c("roi","CS"),
                                idvar = "ID")


pp.side.roi.deval <- ggplot(data = side.bs, aes(x = CS, y = index, fill = CS, color = CS)) +
  geom_abline(slope = 0, intercept = 0, linetype= "dashed", color="black") +
  geom_point(aes(x = CS,
                 y = index, group = ID), 
             alpha = 0.3, position = position_dodge(0.6)) +
  geom_flat_violin(aes(x = 2.7), alpha = .5, position = position_nudge(x = -.15, y = 0), 
                   adjust = 1.8, trim = F, color = NA) +
  geom_crossbar(data = sumstat.side , aes(x = CS, y = index, ymin=index, ymax=index), width = 0.8 , alpha = 0.0) +
  geom_errorbar(data = sumstat.side , aes(x = CS, y = index, ymin = index-se, ymax = index+se), width = 0.7) +
  theme_bw() +
  facet_wrap(~roi, nrow = 1) +
  scale_y_continuous(limits = c(-6.5,4.5)) +
  scale_fill_manual(values=c("black","#FFCC33")) + 
  scale_color_manual(values=c("black","#9D6C00")) +
  labs(
    title = '',
    x = 'CS',
    y = "Betas (after - before)"
  ) 


pp.side.roi.deval <- pp.side.roi.deval + theme_bw( ) + pannel_theme

pdf(file.path(figures_path,'Supplementary_Figure_9.pdf'),width=8,height=5)
print(pp.side.roi.deval)
dev.off()







#-------------------------------------------------------------------------------
# R1  Comment 6: Distribution of the learning rates parameters                 
#-------------------------------------------------------------------------------

# ----------------------------format database ----------------------------------
free_parameters = aggregate(alpha ~ ID + eta, data = MDL,  function(x) mean(x, na.rm =T))

db.free_parameters <- gather(free_parameters,parameter_name, parameter_value, eta:alpha, factor_key=TRUE)

db.free_parameters$parameter_name <- recode(db.free_parameters$parameter_name, "alpha" = "\u03b1",
                                            "eta" = "\u03B7")

db.free_parameters.bg = summarySEwithin(db.free_parameters,
                                        measurevar = c("parameter_value"),
                                        withinvars = c("parameter_name"),
                                        idvar = "ID")

corr.test(free_parameters$eta,free_parameters$alpha)



# --------------------------- plot raw data ------------------------------------
parameter_distribution <- ggplot(data = db.free_parameters, 
                                 aes(x = 1, y = parameter_value,fill = parameter_name, 
                                                                color = parameter_name) ) +
  
  geom_point(aes(x = 1, y = parameter_value, group = ID), alpha = 0.9, position = position_dodge(0.6)) +
  
  geom_crossbar(data = db.free_parameters.bg , aes(x = 1,y = parameter_value, ymin=parameter_value, ymax=parameter_value), 
                width = 0.8 , alpha = 0) +
  geom_errorbar(data = db.free_parameters.bg , aes(x = 1, ymin = parameter_value-se, ymax = parameter_value+se), 
                width = 0.7) +

  facet_wrap(~parameter_name) +

  geom_flat_violin(aes(x = 1.8), alpha = .3, position = position_nudge(x = -.15, y = 0), 
                   adjust = 1.8, trim = F, color = NA) +
  scale_fill_viridis_d(aesthetics = c("color", "fill"),begin = 0, end = 0.5, option = "viridis") +
  theme_bw() +
  labs(x = '', y = "Estimated value")


parameter_distribution = parameter_distribution + timeline_theme +theme(legend.position="none") +
  theme( panel.grid.major.y = element_blank(),  axis.title.x =element_blank(),
         axis.ticks.x=element_blank(), 
         axis.text.x=element_blank())

cairo_pdf(file.path(figures_path,'Supplementary_Figure_3.pdf'))
print(parameter_distribution)
dev.off() 




#-------------------------------------------------------------------------------
# R1  Comment 6: model based predictors of the pupil                 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# R1 Minor Comment8 :   Control MVPA procedure with runs as foldes                  
#-------------------------------------------------------------------------------

#----------------------- DISTRIBUTION  ID  --------------------------------------


# --------------------- IFG
# MEAN ACC AND CI
t.test(R.LEARN_ID_ifg$IFG_R, mu = 0.5); se(R.LEARN_ID_ifg$IFG_R) 
# BF
ttestBF(R.LEARN_ID_ifg$IFG_R, mu = 0.5)

# --------------------- IPS

# MEAN ACC AND CI
t.test(R.LEARN_ID_ips$IPS_R, mu = 0.5); se(R.LEARN_ID_ips$IPS_R) 
# BF
ttestBF(R.LEARN_ID_ips$IPS_R, mu = 0.5)

# --------------------- PCL_L

# MEAN ACC AND CI
t.test(R.LEARN_ID_pcl$PCL_L, mu = 0.5); se(R.LEARN_ID_pcl$PCL_L) 
# BF
ttestBF(R.LEARN_ID_pcl$PCL_L, mu = 0.5)

# --------------------- PCG_L

# MEAN ACC AND CI
t.test(R.LEARN_ID_pcg$PCG_L, mu = 0.5); se(R.LEARN_ID_pcg$PCG_L) 
# BF
ttestBF(R.LEARN_ID_pcg$PCG_L, mu = 0.5)

#---------------------- DISTRUBITION SIDE --------------------------------------

# --------------------- cuneus
# MEAN ACC AND CI
t.test(R.LEARN_SIDE_cuneus$Cuneus, mu = 0.5); se(R.LEARN_SIDE_cuneus$Cuneus) 
# BF
ttestBF(R.LEARN_SIDE_cuneus$Cuneus, mu = 0.5)


# --------------------- ips
# MEAN ACC AND CI
t.test(R.LEARN_SIDE_ips$R_IPS, mu = 0.5); se(R.LEARN_SIDE_ips$R_IPS) 
# BF
ttestBF(R.LEARN_SIDE_ips$R_IPS, mu = 0.5)

# --------------------- stg l (supra marginal gyrus)
# MEAN ACC AND CI
t.test(R.LEARN_SIDE_stg_l$STG_L, mu = 0.5); se(R.LEARN_SIDE_stg_l$STG_L) 
# BF
ttestBF(R.LEARN_SIDE_stg_l$STG_L, mu = 0.5)


# --------------------- stg r (supra marginal gyrus)
# MEAN ACC AND CI
t.test(R.LEARN_SIDE_stg_r$STG_R, mu = 0.5); se(R.LEARN_SIDE_stg_r$STG_R) 
# BF
ttestBF(R.LEARN_SIDE_stg_r$STG_R, mu = 0.5)

# --------------------- stg average
# MEAN ACC AND CI
t.test(Plot.R.SIDE_acc$STG, mu = 0.5); se(Plot.R.SIDE_acc$STG) 
# BF
ttestBF(Plot.R.SIDE_acc$STG, mu = 0.5)


# --------------------- lat occ
# MEAN ACC AND CI
t.test(R.LEARN_SIDE_latocc$R_LATOCC, mu = 0.5); se(R.LEARN_SIDE_latocc$R_LATOCC) 
# BF
ttestBF(R.LEARN_SIDE_latocc$R_LATOCC, mu = 0.5)




#----------------------------------- PLOT --------------------------------------


# plot ID
sumstat.id.control.acc <- summarySE(Plot.R.ID_acc.long,
                                    measurevar = c("accuracy"),
                                    groupvars = "ROI")


ID.control.acc.pp <- ggplot(data = Plot.R.ID_acc.long, aes (x=ROI, y = accuracy, fill = ROI, color = ROI)) +
  geom_abline(slope=0, intercept=0.5, linetype = "dashed", color = "black") +
  geom_point(aes(color = ROI),position = position_jitterdodge(jitter.width = .5, jitter.height = 0),alpha = .2) +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .1, aes(fill = ROI, color = NA), color = NA) +
  geom_crossbar(data = sumstat.id.control.acc, aes(y = accuracy, ymin=accuracy, ymax=accuracy), width = 0.5 , alpha = 0.0) +
  geom_errorbar(data = sumstat.id.control.acc, aes(x = ROI, ymin = accuracy-se, ymax = accuracy+se), width = 0.3) +
  theme_bw() +
  scale_y_continuous( limits = c(0.3,0.7)) +
  scale_fill_manual(values=c("#CC0000","#CC0000","#CC0000","#CC0000")) + 
  scale_color_manual(values=c("#660000","#660000","#660000","#660000")) +
  theme_bw() +
  labs(
    title = '',
    x = 'ROI',
    y = "Decoding accuracy"
  ) 

ID.control.acc.pp  = ID.control.acc.pp + timeline_theme + theme(legend.position="none") 



# Plot SIDE

Plot.R.SIDE_acc.long$ROI <- dplyr::recode(Plot.R.SIDE_acc.long$ROI, "STG" = "SMG")

sumstat.side.control.acc <- summarySE(Plot.R.SIDE_acc.long,
                                      measurevar = c("accuracy"),
                                      groupvars = "ROI")

SIDE.control.acc.pp <- ggplot(data = Plot.R.SIDE_acc.long, aes (x=ROI, y = accuracy, fill = ROI, color = ROI)) +
  geom_abline(slope=0, intercept=0.5, linetype = "dashed", color = "black") +
  geom_point(aes(color = ROI),position = position_jitterdodge(jitter.width = .5, jitter.height = 0),alpha = .2) +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .1, aes(fill = ROI, color = NA), color = NA) +
  geom_crossbar(data = sumstat.side.control.acc, aes(y = accuracy, ymin=accuracy, ymax=accuracy), width = 0.5 , alpha = 0.0) +
  geom_errorbar(data = sumstat.side.control.acc, aes(x = ROI, ymin = accuracy-se, ymax = accuracy+se), width = 0.3) +
  theme_bw() +
  scale_fill_manual(values=c("#FFCC33","#FFCC33","#FFCC33","#FFCC33")) + 
  scale_color_manual(values=c("#9D6C00","#9D6C00","#9D6C00","#9D6C00")) +
  theme_bw() +
  scale_y_continuous( limits = c(0.3,0.9)) +
  labs(
    title = '',
    x = 'ROI',
    y = "Decoding accuracy"
  ) 


SIDE.control.acc.pp  = SIDE.control.acc.pp + timeline_theme + theme(legend.position="none") 


control.acc.pp = plot_grid(SIDE.control.acc.pp, ID.control.acc.pp, ncol = 1, labels=c('A', 'B'), scale = 0.9, label_size = 28)


pdf(file.path(figures_path,'Supplementary_Figure_2.pdf'),width=12,height=8)
print(control.acc.pp)
dev.off() 


# ----------------model based value  derived from the FWD model ---------------

# frequentistic approach
mdl_value = formula(pupil ~ MB_VV + ( MB_VV|ID))
value_mod = lmer(mdl_value, data = MDL, control = my_control) 
mb_value  = coef(value_mod)
summary(value_mod)
confint(value_mod, level = 0.95, method = "Wald")


# BF
fit_pupil_mb_val <- brm(scale(pupil) ~ MB_VV + (MB_VV |ID), 
                        prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                        data =MDL, # on the aggregate data
                        iter = 40000, warmup=5000,
                        family = gaussian(), save_pars = save_pars(all = TRUE))

describe_posterior(fit_pupil_mb_val , estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_pupil_mb_val, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


# ----------------model based value  derived from the RW model ---------------

# frequentistic approach
mf_mdl_value = formula(pupil ~ MF_VV + (MF_VV|ID))
mf_value_mod = lmer(mf_mdl_value, data = MDL, control = my_control) 
mf_value     = coef(mf_value_mod)

summary(mf_value_mod )
confint(mf_value_mod , level = 0.95, method = "Wald")

# BF
fit_pupil_mf_val <- brm(scale(pupil) ~ MF_VV + (MF_VV |ID), 
                        prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                        data =MDL, # on the aggregate data
                        iter = 40000, warmup=5000,
                        family = gaussian(), save_pars = save_pars(all = TRUE))

describe_posterior(fit_pupil_mf_val , estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_pupil_mf_val, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))





#-------------------------------------------------------------------------------
# R1  Minor Comment3 : Liking ratings of the sweet and the salty outcomes                 
#-------------------------------------------------------------------------------

# format dataset
FOOD.m <- aggregate(outcome_liking ~ US_ID*ID , data = FOOD, function(x) mean(x, na.rm =T))

FOOD.m$US_ID = factor(FOOD.m$US_ID)
FOOD.m.bg = summarySEwithin(FOOD.m,
                            measurevar = c("outcome_liking"),
                            withinvars = c("US_ID"),
                            idvar = "ID")


# statistical test frequentistic
mf = formula(outcome_liking ~ US_ID * run + (US_ID + run |ID)) # ! simplified error term
check_id_val = lmer(mf, data = FOOD, control = my_control)
summary(check_id_val)
confint(check_id_val, level = 0.95, method = "Wald") 

# Bayesian analysis
FOOD$US_ID_b = ifelse(FOOD$US_ID == "sweet", -1, 1)

fit_IDval <- brm( scale(outcome_liking) ~ US_ID_b * scale(run) + (US_ID_b * scale(run) |ID), 
                  prior =  c(prior(normal(0,0.5), class="b", coef=""),prior(cauchy(0,0.5), class="sd")),
                  data =FOOD,
                  family = gaussian(), save_pars = save_pars(all = TRUE))

describe_posterior(fit_IDval, estimate = "median", 
                   dispersion = T, ci = .9, ci_method = "hdi", 
                   bf_prior = fit_IDval, diagnostic = "Rhat",  
                   test = c("p_direction", "bf"))


# Plot
p_US_ID_check <- ggplot(data = FOOD.m, aes(x = US_ID, y = outcome_liking, 
                                           fill = US_ID, color = US_ID)) +
  geom_abline(slope= 0, intercept=5, linetype = "dashed", color = "gray") +
  geom_point(aes(x = US_ID, y = outcome_liking, group = ID), 
             alpha = 0.9, position = position_dodge(0.5)) +
  geom_crossbar(data = FOOD.m.bg, aes(x = US_ID, y = outcome_liking,ymin=outcome_liking, ymax=outcome_liking), 
                width = 0.7 , alpha = 0) +
  geom_errorbar(data = FOOD.m.bg, aes(x = US_ID, ymin = outcome_liking-se, ymax = outcome_liking+se), 
                width = 0.6) +
  geom_flat_violin(aes(x = 2.7), alpha = .5, position = position_nudge(x = -.15, y = 0), 
                   adjust = 1.8, trim = F, color = NA) +
  scale_x_discrete(labels = c("Sweet", "Salty")) +
  scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0, 12), expand = c(0, 0))+
  scale_fill_viridis_d(aesthetics = c("color", "fill"),begin = 0.3, end = 0.7, option = "viridis") +
  theme_bw() +
  labs( x = "Outcome", y = "Pleasantness of the outcome")


p_US_ID_check = p_US_ID_check+ timeline_theme + theme(legend.position="none") +
  theme( panel.grid.major.y = element_blank())

pdf(file.path(figures_path,'Supplementary_Figure_1.pdf'))
print(p_US_ID_check )
dev.off() 

#-------------------------------------------------------------------------------
#                  REVIEWER 2
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# R2 Comment4 : Role of hunger across sessions                
#-------------------------------------------------------------------------------

# format database
HUNGER <- aggregate (hunger ~ run*ID, data = FOOD, function(x) mean(x, na.rm =T))
HUNGER$ID= factor(HUNGER$ID)
HUNGER$run= factor(HUNGER$run)


HUNGER.bg = summarySEwithin(HUNGER,
                            measurevar = c("hunger"),
                            withinvars = c("run"),
                            idvar = "ID")

# create plot
hunger.pp <- ggplot(data = HUNGER, aes(x = run,  y = hunger, fill = run, color = run)) +
  geom_abline(slope= 0, intercept=5, linetype = "dashed", color = "black") +
  geom_point(aes(x = run, y = hunger, group = ID), alpha = 0.7, position = position_dodge(0.3)) +
  geom_crossbar(data = HUNGER.bg, aes(x = run,  y = hunger, ymin=hunger, ymax=hunger), 
                width = 0.6 , alpha = 0) +
  geom_errorbar(data = HUNGER.bg, aes(x = run, ymin = hunger-se, ymax = hunger+se), 
                width = 0.7) +
  geom_flat_violin(aes(x = 3.7), alpha = .5, position = position_nudge(x = -.15, y = 0), 
                   adjust = 1.8, trim = F, color = NA) +
  scale_y_continuous(breaks = seq(0, 10, 1), limits = c(-0.5, 12), expand = c(0, 0)) +
  scale_fill_viridis_d(aesthetics = c("color", "fill"),begin = 0, end = 0.75, option = "magma", direction = -1) +
  theme_bw() +
  annotate("rect", xmin=0.6, xmax=2.4, ymin=0, ymax=11, alpha=0.2, fill="gray") +
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=11, alpha=0.2, fill="gray") +
  annotate("text", x = 1.5, y = 0.5, label = "Acquisition",size = 5) +
  annotate("text", x = 3, y = 0.5, label = "Test",size = 5) +
  annotate("segment", x = 2, xend = 3, y = 11.1, yend = 11.1) +
  annotate("segment", x = 1, xend = 3, y = 11.6, yend = 11.6) +
  annotate("text", x = 2.5, y = 11.2, label = "**", size = 6) +
  annotate("text", x = 2, y = 11.7, label = "**", size = 6) +
  labs( x = "Run", y = "Hunger ratings")

hunger.pp = hunger.pp+ timeline_theme + theme(legend.position="none") +
  theme( panel.grid.major.y = element_blank())


pdf(file.path(figures_path,'Supplementary_Figure_4.pdf'))
print(hunger.pp)
dev.off() 

# quick test for plot
anova(lmer(hunger ~ run + (1| ID), data = subset(HUNGER, run != "1"), control = my_control))
anova(lmer(hunger ~ run + (1| ID), data = subset(HUNGER, run != "3"), control = my_control))
anova(lmer(hunger ~ run + (1| ID), data = subset(HUNGER, run != "2"), control = my_control))

