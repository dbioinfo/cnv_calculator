library(tidyverse)
library(cbioportalR)
library(Matrix)

#usethis::edit_r_environ()
#paste CBIOPORTAL_TOKEN = '8b7286ee-9372-483a-9c1b-64cfd4d3c2e3'
#restart

#set db (can config local copies, custom searches)
set_cbioportal_db("public")

#set specific databases
all_studies <- available_studies()
study_set <- all_studies %>% filter(grepl('_tcga_pan_can_atlas',studyId)) %>% pull(studyId)

#gather data in chunks (of roughly 10 studies)
tdat <- list()
for (study in study_set[30:32]){
  #progress
  print(study)
  
  #download data
  cna <- get_cna_by_study(study)

  #process into co-occurrence dups/dels
  dups <- cna %>% filter(alteration>0)
  dels <- cna %>% filter(alteration<0)
  dup_co_freq <- crossprod(table(dups[c('uniquePatientKey','hugoGeneSymbol')]))
  del_co_freq <- crossprod(table(dels[c('uniquePatientKey','hugoGeneSymbol')]))
  
  #save data in sparse matrices to conserve memory
  tdat[[study]] <- list(dups=Matrix(data=dup_co_freq,
                                    sparse=T), 
                        dels=Matrix(data=del_co_freq,
                                    sparse=T))
  
}

#savepoint for future merge
#saveRDS(tdat, file='data/cnvco_part1.rds')
#saveRDS(tdat, file='data/cnvco_part2.rds')
#saveRDS(tdat, file='data/cnvco_part3.rds')
saveRDS(tdat, file='data/cnvco_part4.rds')




#####
#add pan origimed
study <- 'pan_origimed_2020'
cancer_types<-get_clinical_by_study('pan_origimed_2020')
cancer_types <- cancer_types %>% filter(clinicalAttributeId=='CANCER_TYPE') %>% select(uniquePatientKey, value)
cna <- get_cna_by_study(study)
temptdat <- list()

#iter cancer types
for (i in unique(cancer_types$value)){
  #subset this cancer
  patients <- cancer_types %>% filter(value==i) %>% pull(uniquePatientKey)
  tcna <- cna %>% filter(uniquePatientKey %in% patients)
  
  #calculate co occurrence
  dups <- tcna %>% filter(alteration>0)
  dels <- tcna %>% filter(alteration<0)
  dup_co_freq <- crossprod(table(dups[c('uniquePatientKey','hugoGeneSymbol')]))
  del_co_freq <- crossprod(table(dels[c('uniquePatientKey','hugoGeneSymbol')]))
  diag(dup_co_freq) <- 0
  diag(del_co_freq) <- 0
  
  #store cancer type and matrices
  temptdat[[i]] <- list(dups=dup_co_freq, dels=del_co_freq)
}

#save
saveRDS(tdat, file='data/cnvco_part5.rds')
