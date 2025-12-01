library(tidyverse)
library(cbioportalR)

#usethis::edit_r_environ()
#paste CBIOPORTAL_TOKEN = '8b7286ee-9372-483a-9c1b-64cfd4d3c2e3'
#restart

#set db (can config local copies, custom searches)
set_cbioportal_db("public")

#set specific databases
all_studies <- available_studies()
study_set <- all_studies %>% filter(grepl('_tcga_pan_can_atlas',studyId)) %>% pull(studyId)

#gather data if not downloaded
tdat <- data.frame(matrix(nrow=0, ncol=6))
colnames(tdat) <- c("hugoGeneSymbol", "alteration", "n","total_samp", "cancer_type","rel_incidence")
for (study in study_set){
  print(study)
  cna <- get_cna_by_study(study)
  gene_cna <- cna %>% 
    group_by(hugoGeneSymbol) %>% 
    count(alteration) %>% 
    ungroup() %>% 
    mutate(total_samp=(get_study_info(study_id = study) %>% pull(cnaSampleCount)), #get total number of samples in dataset
           cancer_type=get_study_info(study_id = study) %>% pull(cancerType.shortName), #labels
           rel_incidence=n/total_samp) #calculate Pr(SV|C)
  tdat <- rbind(tdat, gene_cna)
}

tdat <- tdat %>% mutate(alteration=case_when(alteration<0~'del', alteration>0~'dup'))
#save data if not downloaded
#saveRDS(file='data/tcga_pancan_cna.rds',tdat)

#load and add origimed pancan atlas
tdat <- readRDS(file='data/tcga_pancan_cna.rds')

study <- 'pan_origimed_2020'
cancer_types<-get_clinical_by_study('pan_origimed_2020')
cancer_types <- cancer_types %>% filter(clinicalAttributeId=='CANCER_TYPE') %>% select(uniquePatientKey, value)
cna <- get_cna_by_study(study)
temptdat <- data.frame(matrix(nrow=0, ncol=6))
colnames(temptdat) <- c("hugoGeneSymbol", "alteration", "n","total_samp", "cancer_type","rel_incidence")

#iter cancer types
for (i in unique(cancer_types$value)){
  patients <- cancer_types %>% filter(value==i) %>% pull(uniquePatientKey)
  tcna <- cna %>% filter(uniquePatientKey %in% patients)
  
  gene_cna <- tcna %>% 
    group_by(hugoGeneSymbol) %>% 
    count(alteration) %>% 
    ungroup() %>% 
    mutate(total_samp=length(patients), #get total number of samples in dataset
           cancer_type=i, #labels
           rel_incidence=n/length(patients)) #calculate Pr(SV|C)
  
  temptdat <- rbind(temptdat, gene_cna)
}

#save
saveRDS(file='data/origimed_pancan_cna.rds',temptdat)
