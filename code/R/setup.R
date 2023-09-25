####################################################################################################

### Set Up for PAVMOD fmri 

# last  modified on June 2023
####################################################################################################




#------------------------ SETUP PATH ------------------------------------------------------

# set working directory for outputs
figures_path  <- file.path(home_path, 'figures')
dabases_path  <- file.path(home_path, "databases")




#------------------------ GET TOOLS --------------------------------------------------------------


devtools::source_gist("2a1bb0133ff568cbe28d", 
                      filename = "geom_flat_violin.R")

devtools::source_gist("383aa93ffa161665c0dca1103ef73d9d", 
                      filename = "effect_CI.R")





#------------------------ PLOT SETTINGS ---------------------------------------------------------

# set theme for plots
averaged_theme <- theme_bw(base_size = 20, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(size=.2, color="lightgrey"),
        panel.grid.minor = element_blank(),
        legend.text  = element_text(size =  14),
        legend.title = element_text(size =  14),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size =  20),
        axis.line = element_line(size = 0.5))
#panel.border = element_blank())



timeline_theme <- theme_bw(base_size = 20, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 20, face = "bold"),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(size=.2, color="lightgrey"),
        panel.grid.minor = element_blank(),
        legend.text  = element_text(size =  14),
        legend.title = element_text(size =  14),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size =  20),
        axis.line = element_line(size = 0.5))
#panel.border = element_blank())

timeline_theme_minor_grid <- theme_bw(base_size = 20, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 20, face = "bold"),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(size=.2, color="lightgrey"),
        panel.grid.minor = element_line(size=.1, color="lightgrey"),
        legend.text  = element_text(size =  14),
        legend.title = element_text(size =  14),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size =  20),
        axis.line = element_line(size = 0.5))

pal = viridis::inferno(n=5)




# ----------------------  ANALYSIS  SETTINGS ----------------------------------------


#my_control1 = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))


my_control  = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)) #!!


options(contrasts = c("contr.sum", "contr.poly")) # this is key


# ---------------------- OPEN AND CLEAN BEHAVIOR DATASETS --------------------------------------------------


# open behavioral datasets
PAV  <- read.delim(file.path(dabases_path,'behavior','Database-PAVMOD-behavior.txt'), header = T, sep ='') # read dataset
FOOD <- read.delim(file.path(dabases_path,'behavior','Database-PAVMOD-outcome_liking.txt'), header = T, sep ='') # read dataset
CS   <- read.delim(file.path(dabases_path,'behavior','Database-PAVMOD-liking.txt'), header = T, sep ='') # read dataset
MDL  <- read.delim(file.path(dabases_path,'behavior','Database-PAVMOD-modeling-revisions.txt'), header = T, sep ='') # read dataset

# define factors
fac <- c("ID") 
MDL[fac] <- lapply(MDL[fac], factor)

# join behavioral datasets
PAVMOD <- join(CS, PAV, type = 'full')
#PAVMOD <- join(PAVMOD1,FOOD, type = 'full')

# define factors
fac <- c("ID", "trial", "run", "bin", "congr_ID", "congr_SIDE", "CS_ID", "CS_SIDE", "CS_ALL") 
PAVMOD[fac] <- lapply(PAVMOD[fac], factor)


# FOOD evaluation for plots
FOOD.c = subset(FOOD, !run == '1')
FOOD.c$run = factor(FOOD.c$run)
FOOD.c$run <- dplyr::recode(FOOD.c$run, "2" = "Pre devaluation", "3" = "Post devaluation" )

FOOD.c$US_value = factor(FOOD.c$US_value)
FOOD.bg = summarySEwithin(FOOD.c,
                          measurevar = c("outcome_liking"),
                          withinvars = c("US_value","run"),
                          idvar = "ID")



# ---------------------- OPEN AND CLEAN MB ROI DATASETS --------------------------------------------------


# read roi database
LEARN_RPE    <- read.delim(file.path(dabases_path,'RL_ROI','Database-PAVMOD-betasRPE_ROI_glmRPE.txt'), header = T, sep ='') # read dataset 
LEARN_SPE    <- read.delim(file.path(dabases_path,'RL_ROI','Database-PAVMOD-betasSPE_ROI_glmSPE.txt'), header = T, sep ='') # read dataset 

# get it in long format RPE
LEARN_RPE_long <- gather(LEARN_RPE, ROI, betas, RPE_in_RPE_VMPF_ROI:RPE_in_RPE_VTA_ROI, factor_key=TRUE)

LEARN_RPE_long[grepl("VMPF", LEARN_RPE_long$ROI), "roi"] <- "RPE_VMPF"
LEARN_RPE_long[grepl("VS", LEARN_RPE_long$ROI), "roi"] <- "RPE_VS"
LEARN_RPE_long[grepl("VTA", LEARN_RPE_long$ROI), "roi"] <- "RPE_VTA"

LEARN_RPE_long$side = "none"
LEARN_RPE_long$feature = "RPE"


# get it in long format SPE
LEARN_SPE_long <- gather(LEARN_SPE, ROI, betas, SPE_in_SPE_SFG_ROI:SPE_in_SPE_right_OFC_ROI, factor_key=TRUE)

LEARN_SPE_long[grepl("SFG", LEARN_SPE_long$ROI), "roi"] <- "SPE_SFG"
LEARN_SPE_long[grepl("VTA", LEARN_SPE_long$ROI), "roi"] <- "SPE_VTA"
LEARN_SPE_long[grepl("MFG", LEARN_SPE_long$ROI), "roi"] <- "SPE_MFG"
LEARN_SPE_long[grepl("OFC", LEARN_SPE_long$ROI), "roi"] <- "SPE_OFC"

LEARN_SPE_long$side = "none"
LEARN_SPE_long[grepl("left", LEARN_SPE_long$ROI), "side"] <- "left"
LEARN_SPE_long[grepl("right", LEARN_SPE_long$ROI), "side"] <- "right"

LEARN_SPE_long$feature = "SPE"


CHANGE_RPE    <- read.delim(file.path(dabases_path,'RL_ROI','Database-PAVMOD_betasUSdeval_ROI_glmRPE.txt'), header = T, sep ='') # read dataset 
CHANGE_SPE    <- read.delim(file.path(dabases_path,'RL_ROI','Database-PAVMOD_betasUSdeval_ROI_glmSPE.txt'), header = T, sep ='') # read dataset 

# get it in long format RPE
CHANGE_RPE_long <- gather(CHANGE_RPE, ROI, betas,  value_corr_in_RPE_VMPF_ROI:devalue_corr_in_RPE_VTA_ROI, factor_key=TRUE)

CHANGE_RPE_long[grepl("VMPF", CHANGE_RPE_long$ROI), "roi"] <- "RPE_VMPF"
CHANGE_RPE_long[grepl("VS", CHANGE_RPE_long$ROI), "roi"] <- "RPE_VS"
CHANGE_RPE_long[grepl("VTA", CHANGE_RPE_long$ROI), "roi"] <- "RPE_VTA"

CHANGE_RPE_long$side    = "none"
CHANGE_RPE_long$feature = "RPE"
CHANGE_RPE_long$value   = "valued"
CHANGE_RPE_long[grepl("devalue", CHANGE_RPE_long$ROI), "value"] <- "devalued"

fac <- c("ID","side", "roi", "value", "feature") 
CHANGE_RPE_long[fac] <- lapply(CHANGE_RPE_long[fac], factor)

# get it in long format SPE
CHANGE_SPE_long <- gather(CHANGE_SPE, ROI, betas, value_corr_in_SPE_SFG_ROI:devalue_corr_in_SPE_right_OFC_ROI, factor_key=TRUE)

CHANGE_SPE_long[grepl("SFG", CHANGE_SPE_long$ROI), "roi"] <- "SPE_SFG"
CHANGE_SPE_long[grepl("VTA", CHANGE_SPE_long$ROI), "roi"] <- "SPE_VTA"
CHANGE_SPE_long[grepl("MFG", CHANGE_SPE_long$ROI), "roi"] <- "SPE_MFG"
CHANGE_SPE_long[grepl("OFC", CHANGE_SPE_long$ROI), "roi"] <- "SPE_OFC"

CHANGE_SPE_long$side = "none"
CHANGE_SPE_long[grepl("left", CHANGE_SPE_long$ROI), "side"] <- "left"
CHANGE_SPE_long[grepl("right", CHANGE_SPE_long$ROI), "side"] <- "right"

CHANGE_SPE_long$feature = "SPE"
CHANGE_SPE_long$value   = "valued"
CHANGE_SPE_long[grepl("devalue", CHANGE_SPE_long$ROI), "value"] <- "devalued"

fac <- c("ID","side", "roi", "value", "feature") 
CHANGE_SPE_long[fac] <- lapply(CHANGE_SPE_long[fac], factor)


CHANGE_PE_long = join (CHANGE_RPE_long, CHANGE_SPE_long, type = "full")


# REMOVE THE ROI THAT DID NOT SURVIVE CORRECTION
CHANGE_SPE_long = subset(CHANGE_SPE_long, ROI != "value_corr_in_SPE_left_MFG_ROI" & ROI != "devalue_corr_in_SPE_left_MFG_ROI")


# ---------------------- OPEN AND CLEAN MVPA ROI DATASETS --------------------------------------------------
ROI_ID    <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD_betasCSdeval_ROI_mvpaID.txt'), header = T, sep ='') # read dataset 
ROI_SIDE  <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD-betasCSdeval_ROI_mvpaSIDE.txt'), header = T, sep ='') # read dataset 



# ACC during learning in ROI
LEARN_ID_ifg <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD-mvpaID_IFG_R_mvpaID.txt'), header = T, sep ='') # read dataset 
LEARN_ID_ips <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD-mvpaID_IPS_R_mvpaID.txt'), header = T, sep ='') # read dataset 
LEARN_ID_pcg <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD-mvpaID_PCG_L_mvpaID.txt'), header = T, sep ='') # read dataset 
LEARN_ID_pcl <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD-mvpaID_PCL_L_mvpaID.txt'), header = T, sep ='') # read dataset 


LEARN_SIDE_cuneus <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD-mvpaSIDE_Cuneus_mvpaSIDE.txt'), header = T, sep ='') # read dataset 
LEARN_SIDE_ips    <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD-mvpaSIDE_R_IPS_mvpaSIDE.txt'), header = T, sep ='') # read dataset 
LEARN_SIDE_latocc <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD-mvpaSIDE_R_LATOCC_mvpaSIDE.txt'), header = T, sep ='') # read dataset 
LEARN_SIDE_stg_l  <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD-mvpaSIDE_STG_L_mvpaSIDE.txt'), header = T, sep ='') # read dataset 
LEARN_SIDE_stg_r  <- read.delim(file.path(dabases_path,'MVPA_ROI','Database-PAVMOD-mvpaSIDE_STG_R_mvpaSIDE.txt'), header = T, sep ='') # read dataset 


fac = c("ID", "roi","prepost", "CS")


# ............................... ID ROI ................................................................
ID_long <- gather(ROI_ID, ROI, betas, value_corr_pre_in_IFG_R:devalue_corr_post_in_PCL_L, factor_key=TRUE)

# ROI
ID_long[grepl("IFG_R", ID_long$ROI), "roi"] <- "R_IFG"
ID_long[grepl("IPS_R", ID_long$ROI), "roi"] <- "R_IPS"
ID_long[grepl("PCG_L", ID_long$ROI), "roi"] <- "L_PCG"
ID_long[grepl("PCL_L", ID_long$ROI), "roi"] <- "L_PCL"

ID_long[grepl("pre", ID_long$ROI), "prepost"] <- "pre"
ID_long[grepl("post", ID_long$ROI), "prepost"] <- "post"

ID_long[grepl("value", ID_long$ROI), "CS"] <- "valued"
ID_long[grepl("devalue", ID_long$ROI), "CS"] <- "devalued"

# merge and define factors
ID_long[fac] <- lapply(ID_long[fac], factor)
ID_long <- ddply(ID_long, .(ID,roi,CS), transform, index  = betas[prepost=="post"] - betas[prepost=="pre"]) 
ID_long <- ddply(ID_long, .(ID,roi), transform, index2  = index[CS=="valued"] - index[CS=="devalued"]) 

# label 
ID_long$decoding = "ID"

# ............................... SIDE ROI ................................................................

# get it in long format
SIDE_long <- gather(ROI_SIDE, ROI, betas, value_corr_pre_in_Cuneus:devalue_corr_post_in_STG_R, factor_key=TRUE)

# ROI
SIDE_long[grepl("Cuneus", SIDE_long$ROI), "roi"] <- "Cuneus"
SIDE_long[grepl("R_IPS", SIDE_long$ROI), "roi"] <- "IPS"
SIDE_long[grepl("LATOCC", SIDE_long$ROI), "roi"] <- "LACTOCC"
SIDE_long[grepl("STG_L", SIDE_long$ROI), "roi"] <- "STG"
SIDE_long[grepl("STG_R", SIDE_long$ROI), "roi"] <- "STG"

# side
SIDE_long[grepl("Cuneus", SIDE_long$ROI), "laterality"] <- "bilateral"
SIDE_long[grepl("R_IPS", SIDE_long$ROI), "laterality"] <- "right"
SIDE_long[grepl("LATOCC", SIDE_long$ROI), "laterality"] <- "right"
SIDE_long[grepl("STG_L", SIDE_long$ROI), "laterality"] <- "left"
SIDE_long[grepl("STG_R", SIDE_long$ROI), "laterality"] <- "right"

# prepost
SIDE_long[grepl("pre", SIDE_long$ROI), "prepost"] <- "pre"
SIDE_long[grepl("post", SIDE_long$ROI), "prepost"] <- "post"

# value
SIDE_long[grepl("value", SIDE_long$ROI), "CS"] <- "valued"
SIDE_long[grepl("devalue", SIDE_long$ROI), "CS"] <- "devalued"

# merge and define factors
SIDE_long[fac] <- lapply(SIDE_long[fac], factor)
SIDE_long <- ddply(SIDE_long, .(ID,roi,CS), transform, index  = betas[prepost=="post"] - betas[prepost=="pre"]) 
SIDE_long <- ddply(SIDE_long, .(ID,roi), transform, index2  = index[CS=="valued"] - index[CS=="devalued"]) 

#label kind of roi
SIDE_long$decoding = "SIDE"

# average across sides
SIDE_long_m = aggregate(index ~ID + roi + CS + prepost, data = SIDE_long, function(x) mean(x, na.rm =T))
SIDE_long_m_index = aggregate(index2 ~ID + roi + CS + prepost, data = SIDE_long, function(x) mean(x, na.rm =T))





# ---------------------------------- DEFINE CONTRASTS OF INTEREST --------------------


# Value contrast 1
PAVMOD$CS.value[PAVMOD$CS_ID== 'CSmi'] <- -2
PAVMOD$CS.value[PAVMOD$CS_ID== 'val'] <- 1
PAVMOD$CS.value[PAVMOD$CS_ID== 'deval'] <- 1

# Value contrast left
PAVMOD$CS.left[PAVMOD$CS_SIDE== 'CSmi'] <- -1
PAVMOD$CS.left[PAVMOD$CS_SIDE== 'CSpL'] <- 2
PAVMOD$CS.left[PAVMOD$CS_SIDE== 'CSpR'] <- -1

# Value contrast right
PAVMOD$CS.right[PAVMOD$CS_SIDE== 'CSmi'] <- -1
PAVMOD$CS.right[PAVMOD$CS_SIDE== 'CSpL'] <- -1
PAVMOD$CS.right[PAVMOD$CS_SIDE== 'CSpR'] <- 2


fac <- c("CS.value", "CS.left", "CS.right") 
PAVMOD[fac] <- lapply(PAVMOD[fac], factor)


#------------------------- SUBSET  DATABASES  ---------------------------------

# define change database
CHANGE = subset (PAVMOD, !run == '1')
# define learning database
LEARN  = subset (PAVMOD, !run == "3")

LEARN.p.m <- aggregate(CS_pupil ~ CS.value+run+ID , data = LEARN, function(x) mean(x, na.rm =T))
LEARN.dwl.m <- aggregate(ANT_DW_L ~ CS.left+run+ID , data = LEARN, function(x) mean(x, na.rm =T))
LEARN.dwr.m <- aggregate(ANT_DW_R ~ CS.right+run+ID , data = LEARN, function(x) mean(x, na.rm =T))


# --------------------------- COMPUTE LEARNING INDEXES -------------------------
RT.mean    <- aggregate(US_RT ~ID + congr_ID+ congr_SIDE  +run, data = subset(LEARN, CS_ID !="CSmi"), function(x) mean(x, na.rm =T))

# ID RT prediction error index
ID.mean = subset(RT.mean, congr_SIDE == "congr") # let's narrow down to ID incongruencies only
ID.mean   <- ddply(ID.mean, .(ID), transform, ID_PE = mean(US_RT[congr_ID == "incongr"]) - mean(US_RT[congr_ID == "congr"]))
PE.ID = aggregate(ID_PE ~ID, data = ID.mean, function(x) mean(x, na.rm =T))


# SIDE RT prediction error
SIDE.mean = subset(RT.mean, congr_ID == "congr")
SIDE.mean <- ddply(SIDE.mean, .(ID), transform, SIDE_PE = mean(US_RT[congr_SIDE == "incongr"]) - mean(US_RT[congr_SIDE == "congr"]))
PE.SIDE = aggregate(SIDE_PE ~ID, data = SIDE.mean, function(x) mean(x, na.rm =T))



# PUPIL INDEX
pupil.mean    <- aggregate(CS_pupil~ID + CS_ID, data = LEARN, function(x) mean(x, na.rm =T))
pupil.mean   <- ddply(pupil.mean, .(ID), transform, delta_pupil = CS_pupil - (CS_pupil[CS_ID == "CSmi"]))
pupil.mean = aggregate(delta_pupil ~ID, data = pupil.mean, function(x) mean(x, na.rm =T))


# DWELL INDEX
dwell.mean   <- aggregate(ANT_DW_congr~ID + CS_ID, data = LEARN, function(x) mean(x, na.rm =T))
dwell.mean   <- ddply(dwell.mean, .(ID), transform, delta_dwell = ANT_DW_congr - (ANT_DW_congr[CS_ID == "CSmi"]))
dwell.mean = aggregate(delta_dwell ~ID, data = dwell.mean, function(x) mean(x, na.rm =T))


# common database
LEARN_BEHAV = join(PE.ID, PE.SIDE, type = "full")
LEARN_BEHAV = join(LEARN_BEHAV, dwell.mean, type = "full")
LEARN_BEHAV = join(LEARN_BEHAV, pupil.mean, type = "full")

# average learning signal
LEARN_BEHAV  = ddply (LEARN_BEHAV, .(ID), transform, MEAN_PE = (SIDE_PE + ID_PE)/2)

# average


# --------------------- COMPUTE DEVALUATION EFFICACY INDEXES -------------------

# pupil
PUPIL.mean <- aggregate(CS_pupil ~ID + CS_ID +run, data = CHANGE, function(x) mean(x, na.rm =T))
PUPIL.mean <- ddply(PUPIL.mean, .(ID), transform, CS_pupil_delta = CS_pupil - CS_pupil[CS_ID=="CSmi"])
PUPIL.mean <- ddply(subset(PUPIL.mean, CS_ID !='CSmi'), .(ID), transform, pupil_change = CS_pupil_delta[run=="3"] - CS_pupil_delta[run=="2"]) 
PUPIL.mean <- ddply(PUPIL.mean, .(ID), transform, pupil_deval = pupil_change[CS_ID=="val"] - pupil_change[CS_ID=="deval"])
PU_devaluation <- aggregate(pupil_deval ~ ID, data = PUPIL.mean, FUN = mean)


# dwell
DW.mean   <- aggregate(ANT_DW_congr ~ ID + CS_ID + run, data = CHANGE, function(x) mean(x, na.rm = T))
DW.mean   <- ddply(DW.mean, .(ID), transform, DW_delta = ANT_DW_congr - ANT_DW_congr[CS_ID=="CSmi"])
DW.mean   <- ddply(subset(DW.mean, CS_ID !='CSmi'), .(ID), transform, dw_change = DW_delta[run=="3"] - DW_delta[run=="2"])
DW.mean   <- ddply(DW.mean, .(ID), transform, dw_deval = dw_change[CS_ID=="val"] - dw_change[CS_ID=="deval"])
DW_devaluation <- aggregate(dw_deval ~ ID, data = DW.mean, FUN = mean)


# create database
DEV_BEHAV = join (DW_devaluation, PU_devaluation, type = 'full')
DEV_BEHAV$ID <- factor(DEV_BEHAV$ID)


# combine data base
CHANGE_PE_long = join(CHANGE_PE_long, DEV_BEHAV, type = "full")





#----------------FORMAT DATASET FOR BRAIN BEHAVIOR ANALYSIS --------------------

LEARN = join (LEARN, LEARN_BEHAV, type = "full")

db.BrainBehav.rpe = join(LEARN_BEHAV, LEARN_RPE_long, type = "full")
db.BrainBehav.spe = join(LEARN_BEHAV, LEARN_SPE_long, type = "full")

db.BrainBehav = join(db.BrainBehav.rpe, db.BrainBehav.spe, type = "full")

fac <- c("roi", "side", "feature") 
db.BrainBehav[fac] <- lapply(db.BrainBehav[fac], factor)
db.BrainBehav.rpe[fac] <- lapply(db.BrainBehav.rpe[fac], factor)
db.BrainBehav.spe[fac] <-  lapply(db.BrainBehav.spe[fac], factor)


#---------------------for plots behavior ---------------------------------------
DB.pert = aggregate(betas ~ ID + SIDE_PE + ID_PE + feature, data = db.BrainBehav,  function(x) mean(x, na.rm =T))

DB.pert_long <- gather(DB.pert, RT_PE, RT, SIDE_PE:ID_PE, factor_key=TRUE)
DB.all = join(DB.pert, DB.pert_long, type = "full")
DB.rpe = subset(DB.all, feature == "RPE")
DB.spe = subset(DB.all, feature == "SPE")


# for plots fmri spe
spe.bs = ddply(CHANGE_SPE_long,.(ID,roi),transform, diff_betas=mean(betas[value == "valued"] - mean(betas[value == "devalued"])))
spe.bg = ddply(spe.bs,.(roi),summarise, diff_betas=mean(diff_betas, na.rm = T))

spe.bs = subset (spe.bs, value == "valued") # let remove repetitions
spe.bs = subset (spe.bs, side != "left") # let remove repetitions


spe.bs$roi <- dplyr::recode(spe.bs$roi, SPE_MFG = 'MFG', 
                            SPE_OFC = 'OFC',
                            SPE_SFG = 'SFG',
                            SPE_VTA = 'midbrain',)


sumstat.spe <- summarySE(spe.bs,
                         measurevar = c("diff_betas"),
                         groupvars = "roi")


# for plots fmri rpe
rpe.bs = ddply(CHANGE_RPE_long,.(ID,roi),transform, diff_betas=mean(betas[value == "valued"] - mean(betas[value == "devalued"])))

rpe.bg = ddply(rpe.bs,.(roi),summarise, diff_betas=mean(diff_betas, na.rm = T))

rpe.bs = subset (rpe.bs, value == "valued")

rpe.bs$roi <- dplyr::recode(rpe.bs$roi, RPE_VTA = 'midbrain', 
                            RPE_VMPF = 'vmPFC',
                            RPE_VS = 'VS')

sumstat.rpe <- summarySE(rpe.bs,
                         measurevar = c("diff_betas"),
                         groupvars = "roi")




# ------------------------------------------------------------------------------
#     open dataset for revisions identity and side prediction errors manual
# ------------------------------------------------------------------------------

# --------------------analysis for devaluation change --------------------------

CHANGE_idSPE.manual    <- read.delim(file.path(dabases_path,'revision_ROI','Database-PAVMOD-betasUSdeval_ROI_glm_idSPE-manual.txt'), header = T, sep ='') # read dataset 
CHANGE_sideSPE.manual  <- read.delim(file.path(dabases_path,'revision_ROI','Database-PAVMOD-betasUSdeval_ROI_glm_sideSPE-manual.txt'), header = T, sep ='') # read dataset 


# get idSPE in long format
CHANGE_idSPE.manual.long <- gather(CHANGE_idSPE.manual, ROI, betas, value_corr_in_idSPE_VTA_ROI:devalue_corr_in_idSPE_VTA_ROI, factor_key=TRUE)

CHANGE_idSPE.manual.long[grepl("VTA", CHANGE_idSPE.manual.long$ROI), "roi"] <- "idSPE_VTA"
CHANGE_idSPE.manual.long$side = "none"

CHANGE_idSPE.manual.long$feature = "idSPE"
CHANGE_idSPE.manual.long$value   = "valued"
CHANGE_idSPE.manual.long[grepl("devalue", CHANGE_idSPE.manual.long$ROI), "value"] <- "devalued"

fac <- c("ID","side", "roi", "value", "feature") 
CHANGE_idSPE.manual.long[fac] <- lapply(CHANGE_idSPE.manual.long[fac], factor)

# get sidSePE in long format
CHANGE_sideSPE.manual.long <- gather(CHANGE_sideSPE.manual, ROI, betas, value_corr_in_sideSPE_MFG_ROI:devalue_corr_in_sideSPE_right_OFC_ROI, factor_key=TRUE)

CHANGE_sideSPE.manual.long [grepl("SFG", CHANGE_sideSPE.manual.long$ROI), "roi"] <- "sideSPE_SFG"
CHANGE_sideSPE.manual.long [grepl("VTA", CHANGE_sideSPE.manual.long$ROI), "roi"] <- "sideSPE_VTA"
CHANGE_sideSPE.manual.long [grepl("MFG", CHANGE_sideSPE.manual.long$ROI), "roi"] <- "sideSPE_MFG"
CHANGE_sideSPE.manual.long [grepl("OFC", CHANGE_sideSPE.manual.long$ROI), "roi"] <- "sideSPE_OFC"

CHANGE_sideSPE.manual.long$side = "none"
CHANGE_sideSPE.manual.long[grepl("left", CHANGE_sideSPE.manual.long$ROI), "side"] <- "left"
CHANGE_sideSPE.manual.long[grepl("right", CHANGE_sideSPE.manual.long$ROI), "side"] <- "right"


CHANGE_sideSPE.manual.long$feature = "sideSPE"
CHANGE_sideSPE.manual.long$value   = "valued"
CHANGE_sideSPE.manual.long[grepl("devalue", CHANGE_sideSPE.manual.long$ROI), "value"] <- "devalued"

fac <- c("ID","side", "roi", "value", "feature") 
CHANGE_sideSPE.manual.long[fac] <- lapply(CHANGE_sideSPE.manual.long[fac], factor)


# --------------------analysis for brain behavior correlation --------------------------

LEARN_idSPE.man    <- read.delim(file.path(dabases_path,'revision_ROI','Database-PAVMOD-betasidSPE_ROI_glmidSPE-manual.txt'), header = T, sep ='') # read dataset 
LEARN_sideSPE.man  <- read.delim(file.path(dabases_path,'revision_ROI','Database-PAVMOD-betassideSPE_ROI_glmsideSPE-manual.txt'), header = T, sep ='') # read dataset 

# get idSPE in long format
LEARN_idSPE.man$roi = "VTA"
LEARN_idSPE.man$betas = LEARN_idSPE.man$idSPE_in_idSPE_VTA_ROI
LEARN_idSPE.man$side = "none"
LEARN_idSPE.man$feature = "idSPE"


# get sideSPE in long format
LEARN_sideSPE.man.long <- gather(LEARN_sideSPE.man, ROI, betas, sideSPE_in_sideSPE_MFG_ROI:sideSPE_in_sideSPE_right_OFC_ROI, factor_key=TRUE)

LEARN_sideSPE.man.long[grepl("VTA", LEARN_sideSPE.man.long$ROI), "roi"] <- "SPE_VTA"
LEARN_sideSPE.man.long[grepl("OFC", LEARN_sideSPE.man.long$ROI), "roi"] <- "SPE_OFC"
LEARN_sideSPE.man.long[grepl("SFG", LEARN_sideSPE.man.long$ROI), "roi"] <- "SPE_SFG"
LEARN_sideSPE.man.long[grepl("MFG", LEARN_sideSPE.man.long$ROI), "roi"] <- "SPE_MFG"


LEARN_sideSPE.man.long$side = "none"
LEARN_sideSPE.man.long[grepl("left", LEARN_sideSPE.man.long$ROI), "side"] <- "left"
LEARN_sideSPE.man.long[grepl("right", LEARN_sideSPE.man.long$ROI), "side"] <- "right"

LEARN_sideSPE.man.long$feature = "sideSPE"

# merge for brain behavior analysis
db.BrainBehav.spe.id.man   = join(LEARN_BEHAV, LEARN_idSPE.man, type = "full")
db.BrainBehav.spe.side.man = join(LEARN_BEHAV, LEARN_sideSPE.man.long, type = "full")


fac <- c("roi", "side", "feature") 
# merge for brain behavior analysis
db.BrainBehav.spe.id.man[fac]   <- lapply(db.BrainBehav.spe.id.man[fac], factor)
db.BrainBehav.spe.side.man[fac] <- lapply(db.BrainBehav.spe.side.man[fac], factor)

# For revisions open datasets of mpva on roi by folded by run

# ACC during learning in ROI
R.LEARN_ID_ifg <- read.delim(file.path(dabases_path,'MVPA_ROI','revisions','id','IFG_R_ACC.txt'), header = T, sep ='') # read dataset 
R.LEARN_ID_ips <- read.delim(file.path(dabases_path,'MVPA_ROI','revisions','id','IPS_R_ACC.txt'), header = T, sep ='') # read dataset 
R.LEARN_ID_pcg <- read.delim(file.path(dabases_path,'MVPA_ROI','revisions','id','PCG_L_ACC.txt'), header = T, sep ='') # read dataset 
R.LEARN_ID_pcl <- read.delim(file.path(dabases_path,'MVPA_ROI','revisions','id','PCL_L_ACC.txt'), header = T, sep ='') # read dataset 


R.LEARN_SIDE_cuneus <- read.delim(file.path(dabases_path,'MVPA_ROI','revisions','side','Cuneus_ACC.txt'), header = T, sep ='') # read dataset 
R.LEARN_SIDE_ips    <- read.delim(file.path(dabases_path,'MVPA_ROI','revisions','side','R_IPS_ACC.txt'), header = T, sep ='') # read dataset 
R.LEARN_SIDE_latocc <- read.delim(file.path(dabases_path,'MVPA_ROI','revisions','side','R_LATOCC_ACC.txt'), header = T, sep ='') # read dataset 
R.LEARN_SIDE_stg_l  <- read.delim(file.path(dabases_path,'MVPA_ROI','revisions','side','STG_L_ACC.txt'), header = T, sep ='') # read dataset 
R.LEARN_SIDE_stg_r  <- read.delim(file.path(dabases_path,'MVPA_ROI','revisions','side','STG_R_ACC.txt'), header = T, sep ='') # read dataset 


# create dabase for plots

Plot.R.ID_acc       <- R.LEARN_ID_ifg
Plot.R.ID_acc$IPS   <- R.LEARN_ID_ips$IPS_R
Plot.R.ID_acc$PCG   <- R.LEARN_ID_pcg$PCG_L
Plot.R.ID_acc$PCL   <- R.LEARN_ID_pcl$PCL_L

Plot.R.ID_acc.long <- gather(Plot.R.ID_acc, ROI, accuracy, IFG_R:PCL, factor_key=TRUE)

Plot.R.ID_acc.long$ROI <- dplyr::recode(Plot.R.ID_acc.long$ROI, 
                                        IFG_R = 'IFG')


Plot.R.SIDE_acc       <- R.LEARN_SIDE_cuneus
Plot.R.SIDE_acc$IPS   <- R.LEARN_SIDE_ips$R_IPS
Plot.R.SIDE_acc$LOC   <- R.LEARN_SIDE_latocc$R_LATOCC
Plot.R.SIDE_acc$STG   <- (R.LEARN_SIDE_stg_l$STG_L + R.LEARN_SIDE_stg_r$STG_R) /2

Plot.R.SIDE_acc.long <- gather(Plot.R.SIDE_acc, ROI, accuracy, Cuneus:STG, factor_key=TRUE)
