###############################################################################


#         Neural substrates of parallel devaluation-sensitive and 
#           devaluation-insensitive Pavlovian learning in humans 


#                     analysis for revisions only

# last  modified on JUNE 2023
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
               BayesFactor, bayestestR, DHARMa, brms,cowplot)


require(lattice)


# set path
current_dir <- dirname(getActiveDocumentContext()$path)
setwd(current_dir)
home_path <- getwd()
setwd(home_path)

# run set up
source(file.path(home_path,'R','setup.R'), echo=F)

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


pp_theme <- theme_bw(base_size = 16, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 16, face = "bold"),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(size=.2, color="lightgrey"),
        panel.grid.minor = element_blank(),
        legend.text  = element_text(size =  12),
        legend.title = element_text(size =  12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size =  16),
        axis.line = element_line(size = 0.5))



#-------------------------------------------------------------------------------
#                  REVIEWER 1
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# R1 Comment 3: create side and identity SPE and test the devaluation
#-------------------------------------------------------------------------------


#--------------------- brain behavior correlation SPE id voxels ----------------

#model idSPE voxels and ID PE in RT
vta.mdl.man = lm(betas ~ ID_PE,  data= db.BrainBehav.spe.id.man)
summary(vta.mdl.man)

# model idSPE voxels and SIDE PE in RT
vta.mdl.man = lm(betas ~ SIDE_PE,  data= db.BrainBehav.spe.side.man)
summary(vta.mdl.man)


#--------------------- brain behavior correlation SPE side voxels --------------  

#model sideSPE voxels and SIDE PE in RT
mf_spe_side = formula(betas ~ (SIDE_PE) + (1|ID))
side_spe_side_mod = lmer(mf_spe_side, data = db.BrainBehav.spe.side.man, control = my_control) 
summary(side_spe_side_mod)
confint(side_spe_side_mod, level = 0.95, method = "Wald") 

# model idSPE voxels and SIDE PE in RT
mf_spe_id = formula(betas ~ (ID_PE) + (1|ID))
side_spe_id_mod = lmer(mf_spe_id, data = db.BrainBehav.spe.side.man, control = my_control) 
summary(side_spe_id_mod)
confint(id_spe_side_mod, level = 0.95, method = "Wald") 



# -------------------- devaluation in state id voxels ------------------------------------------
id.vta.mdl.manual = lme(betas ~ value,  data= subset(CHANGE_idSPE.manual.long, roi == "idSPE_VTA"), random= ~ value|ID)
summary(id.vta.mdl.manual)
intervals(id.vta.mdl.manual, which = "fixed")

# -------------------- devaluation in state id voxels ------------------------------------------
side.mfg.mdl.manual = lme(betas ~ value,  data= subset(CHANGE_sideSPE.manual.long, roi == "sideSPE_MFG"), random= ~ value|ID)
summary(side.mfg.mdl.manual)
intervals(side.mfg.mdl.manual, which = "fixed")

side.ofc.mdl.manual = lme(betas ~ value*side,  data= subset(CHANGE_sideSPE.manual.long, roi == "sideSPE_OFC"), random= ~ value/side|ID)
summary(side.ofc.mdl.manual)
intervals(side.ofc.mdl.manual, which = "fixed")

side.vta.mdl.manual = lme(betas ~ value,  data= subset(CHANGE_sideSPE.manual.long, roi == "sideSPE_VTA"), random= ~ value|ID)
summary(side.vta.mdl.manual)
intervals(side.vta.mdl.manual, which = "fixed")



# -------------------- PLOT BB--------------------------------------------------

sup.plot.side = aggregate(betas ~ID + SIDE_PE, data = db.BrainBehav.spe.side.man, function(x) mean(x, na.rm =T))
sup.plot.side$SPE_type = "sideSPE"

sup.plot.id = aggregate(betas ~ID + ID_PE, data = db.BrainBehav.spe.id.man, function(x) mean(x, na.rm =T))
sup.plot.id$SPE_type = "idSPE"


pp_idRT_idSPE = ggplot(sup.plot.id , aes(x = ID_PE, y = betas, color = SPE_type , fill = SPE_type )) +
  geom_point(alpha = .5, size = 3.5) +
  geom_smooth(method = lm, level = .95, alpha = .1, fullrange=TRUE) +
  scale_fill_manual(values=c("#990000")) +
  scale_color_manual( values=c("#990000")) +
  ylab('Betas of identity SPE (a.u)')+
  xlab('RT [unexpected identity - expected identity] (a.u)') +
  theme_bw()

pp_idRT_idSPE  = pp_idRT_idSPE + pp_theme + theme(legend.position="none") 

pp_idRT_idSPE  <- ggMarginal( pp_idRT_idSPE   + theme(legend.position = "none"),type = "density", 
                                 groupColour = TRUE, groupFill = T, alpha = 0.2)




pp_sideRT_sideSPE = ggplot(sup.plot.side , aes(x = SIDE_PE, y = betas, color = SPE_type , fill = SPE_type )) +
  geom_point(alpha = .5, size = 3.5) +
  geom_smooth(method = lm, level = .95, alpha = .3, fullrange=TRUE) +
  scale_fill_manual(values=c("#FFCC33")) +
  scale_color_manual( values=c("#FFCC33")) +
  ylab('Betas of side SPE (a.u)')+
  xlab('RT [unexpected side - expected side] (a.u)') +
  theme_bw()



pp_sideRT_sideSPE  = pp_sideRT_sideSPE + pp_theme + theme(legend.position="none") 


pp_sideRT_sideSPE  <- ggMarginal(pp_sideRT_sideSPE  + theme(legend.position = "none"),type = "density", 
                         groupColour = TRUE, groupFill = T, alpha = 0.2)




# -------------------- PLOT devaluation effects-----------------------------------------------


# format the dataset
CHANGE_idSPE.manual.long = ddply(CHANGE_idSPE.manual.long,.(ID,roi),transform, diff_betas=mean(betas[value == "valued"] - mean(betas[value == "devalued"])))
CHANGE_idSPE.manual.long = subset (CHANGE_idSPE.manual.long, value == "valued") # let remove repetitions


CHANGE_sideSPE.manual.long = ddply(CHANGE_sideSPE.manual.long,.(ID,roi),transform, diff_betas=mean(betas[value == "valued"] - mean(betas[value == "devalued"])))
CHANGE_sideSPE.manual.long = subset (CHANGE_sideSPE.manual.long, value == "valued") # let remove repetitions

r.pp.deval.spe.id = aggregate (diff_betas ~ID +roi,  data = CHANGE_idSPE.manual.long, function(x) mean(x, na.rm =T))
r.pp.deval.spe.side = aggregate (diff_betas~ID +roi,  data = CHANGE_sideSPE.manual.long, function(x) mean(x, na.rm =T))
r.pp.deval.spes = join (r.pp.deval.spe.id, r.pp.deval.spe.side, type = "full")


r.pp.deval.spes$roi <- dplyr::recode(r.pp.deval.spes$roi, "sideSPE_VTA" = "VTA (side SPE)",
                                     "sideSPE_OFC" = "OFC (side SPE)",
                                     "sideSPE_SFG" = "SFG (side SPE)",
                                     "sideSPE_MFG" = "VTA (side SPE)",
                                     "idSPE_VTA" = "VTA (id SPE)")


sumstat.spes <- summarySE(r.pp.deval.spes,
                         measurevar = c("diff_betas"),
                         groupvars = "roi")



r.deval.pp <- ggplot(data = r.pp.deval.spes, aes (x=roi, y = diff_betas, fill = roi, color = roi)) +
  geom_abline(slope=0, intercept=0, linetype = "dashed", color = "black") +
  geom_point(aes(color = roi),position = position_jitterdodge(jitter.width = .1, jitter.height = 0),alpha = .1) +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .2, aes(fill = roi, color = NA), color = NA) +
  geom_crossbar(data = sumstat.spes, aes(y = diff_betas, ymin=diff_betas, ymax=diff_betas), width = 0.5 , alpha = 0.0) +
  geom_errorbar(data = sumstat.spes, aes(x = roi, ymin = diff_betas-ci, ymax = diff_betas+ci), width = 0.3) +
  coord_flip() +
  scale_y_continuous( limits = c(-5,5)) +
  scale_fill_manual(values=c("#CC0000","#E69F00","#E69F00","#E69F00","#E69F00")) + 
  scale_color_manual(values=c("#660000","#9D6C00","#9D6C00","#9D6C00","#9D6C00" )) +
  theme_bw() +
  labs(
    title = '',
    x = 'ROI based on identity and side PE',
    y = "Betas [valued - devalued]"
  ) 

r.deval.pp = r.deval.pp + pp_theme + theme(legend.position="none") 



first_colon = plot_grid(pp_sideRT_sideSPE, pp_idRT_idSPE, ncol = 1, labels=c('A', 'B'), scale = 0.8, label_size = 28)
side_id_spe = plot_grid(first_colon, r.deval.pp, nrow=1, labels=c('',  'C'), scale = 1,label_size = 28) #Or labels="AUTO"

pdf(file.path(figures_path,'Revisions_idSPE_sideSPE.pdf'),width=12,height=12)
print(side_id_spe)
dev.off() 





#-------------------------------------------------------------------------------
#                  REVIEWER 2
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# R2 Comment1 : Effect of outcome devaluation on RT during the test session                    
#-------------------------------------------------------------------------------

# plot
rt.test.deval.pp <- ggplot(data = RT.test.m, 
                           aes (x=CS_ID, y = US_RT, fill = CS_ID, color = CS_ID)) +
  geom_point(aes(x = CS_ID, y = US_RT, group = ID), position = position_dodge(0.5), 
             alpha = .7) +
  geom_flat_violin(aes(x = 2.7), alpha = .5, position = position_nudge(x = -.15, y = 0), 
                   adjust = 1.8, trim = F, color = NA) +
  geom_crossbar(data = RT.test.m.bg, aes(y = US_RT, ymin=US_RT, ymax=US_RT), 
                width = 0.6 , alpha = 0.0) +
  geom_errorbar(data = RT.test.m.bg,  aes(x = CS_ID, ymin = US_RT-se, ymax = US_RT+se),
                width = 0.7) +
  theme_bw() +
  scale_fill_manual(values=c("black","#56B4E9")) + 
  scale_color_manual(values=c("black","#000066")) +
  theme_bw() +
  labs(title = '', x = 'CS',y = "Normalized RT") 

rt.test.deval.pp  =  rt.test.deval.pp  + timeline_theme + theme(legend.position="none") 

pdf(file.path(figures_path,'Revisions_RT.pdf'))
print(rt.test.deval.pp)
dev.off() 
