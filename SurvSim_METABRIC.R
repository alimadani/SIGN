
rm(list = ls())
library(GSVA)

DataDir <- "/mnt/work1/users/bhklab/Users/Ali/Projects/TSC_Clinical/Data/Metabrick/"
PathwayDir <- "/mnt/work1/users/bhklab/Users/Ali/Projects/TSC_Clinical/Data/Pathways_Dataset/"
###################
ExpMat_All <- readRDS(paste(DataDir, "METABRIC_ExpMat_tumor.rds", sep = "", collapse = ""))
colnames(ExpMat_All) <- gsub("_", "-", colnames(ExpMat_All))

MetaData_Discovery <- read.table(paste(DataDir, "table_Discovery.txt", sep = "", collapse = ""), fill=T, sep = "\t", header = T)
MetaData_Validation <- read.table(paste(DataDir, "table_Validation.txt", sep = "", collapse = ""), fill=T, sep = "\t", header = T)

MetaData_Discovery <- as.matrix(MetaData_Discovery)
MetaData_Discovery[which(MetaData_Discovery[,"Treatment"] == "CT"),"Treatment"] <- "CT/RT"
MetaData_Discovery[which(MetaData_Discovery[,"Treatment"] == "CT/HT"), "Treatment"] <- "CT/HT/RT"
MetaData_Validation <- as.matrix(MetaData_Validation)
MetaData_Validation[which(MetaData_Validation[,"Treatment"] == "CT"), "Treatment"] <- "CT/RT"
MetaData_Validation[which(MetaData_Validation[,"Treatment"] == "CT/HT"), "Treatment"] <- "CT/HT/RT"
#############
setwd("/mnt/work1/users/bhklab/Users/Ali/Projects/TSC_Clinical/Scripts/EnsembleApproach/")
source("ASSIGN_Calculation.R")
source("ASSIGN_Wrapper.R")
source("BubbleSort.R")
source("ExpPhen_Matching.R")
source("ExpPheno_Categorize.R")
source("ExpPhen_Subdividing.R")
source("GeneMatching.R")
source("Genes_SimCal.R")
source("GSVA_Calculation.R")
source("Pathway_Similarity.R")
source("SIGN_Ensemble_SimCal.R")
source("Similarities_Wrapper.R")
source("SimSummary_2Class.R")
source("TSC.R")
source("Pathway_Grouping.R")
source("EventRenaming.R")
source("SIGN_Aggregate.R")
source("Survival_Stats.R")
source("SurvivalStat_PostProcess.R")
#################
PathwaySets <- Pathway_Grouping(PathwayDir, Pattern = "_entrez_ProcessPathways.rds")
PathwaySets <- PathwaySets[-which(names(PathwaySets) %in% c("c5_all"))]
##############
MetaData_Discovery[,"last_follow_up_status"] <- EventRenaming(MetaData_Discovery[,"last_follow_up_status"], Censored_Annot = "a")
MetaData_Validation[,"last_follow_up_status"] <- EventRenaming(MetaData_Validation[,"last_follow_up_status"], Censored_Annot = "a")
################ Matching Sample IDs in Expression and Metadata 
ExpPhenoList_Discovery <- ExpPhen_Matching(ExpMat_All, MetaData_Discovery, "METABRIC_ID")
ExpPhenoList_Validation <- ExpPhen_Matching(ExpMat_All, MetaData_Validation, "METABRIC_ID")
###############
###############
ExpPhenoList_Discovery_Treatments <- ExpPhen_Subdividing(ExpPhenoList_Discovery, "Treatment")
ExpPhenoList_Validation_Treatments <- ExpPhen_Subdividing(ExpPhenoList_Validation, "Treatment")

############# Predicting survival based on similarity in Discovery cohort
ExpList_Discovery_Treatments <- ExpPhenoList_Discovery_Treatments$ExpList_Treatment
MetaList_Discovery_Treatments <- ExpPhenoList_Discovery_Treatments$MetaList_Treatment

ptm <- proc.time()

StatSum_List <- list()
for(SubDivIter in 4:4){ #length(ExpList_Discovery_Treatments)
 ExpMat_TMP <- ExpList_Discovery_Treatments[[SubDivIter]]
 MetaMat_TMP <- MetaList_Discovery_Treatments[[SubDivIter]]
 #####
 ExpMeta_List <- list(ExpMat_TMP, MetaMat_TMP)
 names(ExpMeta_List) <- c("Expression", "Meta")
 Time_ID <- "T"
 Event_ID <- "last_follow_up_status"
 ######
 ExpPhenoCats_List <- ExpPheno_Categorize(ExpMeta_List, Time_ID, Event_ID, Mad_Factor=1, MinNum_ExClass=10, Expression_Log2=FALSE)

 ExpList_TMP <- ExpPhenoCats_List$ExpList
 SIGN_Out <- SIGN_Ensemble_SimCal(ExpList_TMP, RefClassID = c("poor", "good"),
                                  TestClassID = c("poor", "intermediate", "good"), GeneID = "EntrezID", PathwaySets)
 ScoreMat <- SIGN_Aggregate(SIGN_Out, ExpPhenoCats_List$TimeList, ExpPhenoCats_List$EventList)

 Survival_Stat <- Survival_Stats(ScoreMat$Scores, ScoreMat$Time, ScoreMat$Event)
 SurStat_Summary <- SurvivalStat_PostProcess(Survival_Stat)

 StatSum_List[[SubDivIter]] <- SurStat_Summary
}
names(StatSum_List) <- names(ExpList_Discovery_Treatments)

proc.time() - ptm
###########################
##########################
########################## Predicting Survival in Validation Cohort Using Discovery cohort
ExpList_Validation_Treatments <- ExpPhenoList_Validation_Treatments$ExpList_Treatment
MetaList_Validation_Treatments <- ExpPhenoList_Validation_Treatments$MetaList_Treatment

ptm <- proc.time()

StatSum_List <- list()
for(SubDivIter in 6:6){ #length(ExpList_Discovery_Treatments)
 ExpMatDiscovery_TMP <- ExpList_Discovery_Treatments[[SubDivIter]]
 MetaMatDiscovery_TMP <- MetaList_Discovery_Treatments[[SubDivIter]]
 ###########
 ExpMatValidation_TMP <- ExpList_Validation_Treatments[[SubDivIter]]
 MetaMatValidation_TMP <- MetaList_Validation_Treatments[[SubDivIter]]
 ##########
 ExpMeta_List <- list(ExpMatDiscovery_TMP, MetaMatDiscovery_TMP)
 names(ExpMeta_List) <- c("Expression", "Meta")
 Time_ID <- "T"
 Event_ID <- "last_follow_up_status"
 #################
 Validation_Time <- as.numeric(MetaMatValidation_TMP[,Time_ID])
 Validation_Event <- as.character(MetaMatValidation_TMP[,Event_ID]) 
 ##############
 ExpPhenoCats_List <- ExpPheno_Categorize(ExpMeta_List, Time_ID, Event_ID, Mad_Factor=1, MinNum_ExClass=10, Expression_Log2=FALSE)

 ExpList_TMP <- ExpPhenoCats_List$ExpList
 ExpList_TMP[["Validation"]] <- ExpMatValidation_TMP

 SIGN_Out <- SIGN_Ensemble_SimCal(ExpList_TMP, RefClassID = c("poor", "good"), TestClassID = c("Validation"), GeneID = "EntrezID", PathwaySets)
 ScoreMat <- SIGN_Aggregate(SIGN_Out, list(Time = Validation_Time), list(Event = Validation_Event))

 Survival_Stat <- Survival_Stats(ScoreMat$Scores, ScoreMat$Time, ScoreMat$Event)
 SurStat_Summary <- SurvivalStat_PostProcess(Survival_Stat)

 StatSum_List[[SubDivIter]] <- SurStat_Summary
}
names(StatSum_List) <- names(ExpList_Discovery_Treatments)

proc.time() - ptm



