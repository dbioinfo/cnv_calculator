#download data from cbioportal straight into postgresql

library(tidyverse)
library(cBioPortalData)
library(dbplyr)
library(RPostgreSQL)

#set up api
cbio <- cBioPortal()
studies <- getStudies(cbio, buildReport = TRUE)
studies <- studies %>% select(studyId) %>% filter(grepl('tcga_pan_can_atlas_2018', studyId))

#connect to server
#drv <- dbDriver('PostgreSQL')
conn <- dbConnect(drv = RPostgres::Postgres(),
                  dbname='postgres', 
                  host='127.0.0.1', 
                  user='postgres',
                  password='password', 
                  port=5432)

#check if tables exist, if so, delete them
#if (dbExistsTable(conn, 'cnv_main')){dbRemoveTable(conn, 'cnv_main')}


#translate gene symbol to entrez if needed
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
translate <- getBM(attributes = c('entrezgene_id','hgnc_symbol'), mart = ensembl)

#####
# cnv_main

#get max genes available
igenes <- getGenePanel(cbio, genePanelId='UCLA_1202')$entrezGeneId
jgenes <- getGenePanel(cbio, genePanelId='IMPACT505')$entrezGeneId
xlgenes <-  getGenePanel(cbio, genePanelId='Bait_UCSF500')$entrezGeneId
igenes <- unique(c(igenes, jgenes, xlgenes))

#create table with first round of data
study <- studies$studyId[1]
studies <- studies[2:nrow(studies),]


cprefix <- str_split_i(study,'_',1)
tmp <- getDataByGenes(cBioPortal(), 
               genes = igenes, 
               molecularProfileIds = paste0(cprefix,'_tcga_pan_can_atlas_2018_gistic'),  
               sampleListId = paste0(cprefix,'_tcga_pan_can_atlas_2018_all'))
tmp <- tmp[[1]] %>% select(uniquePatientKey, hugoGeneSymbol, value) %>% 
              rename(., value='alteration', hugoGeneSymbol='gene') %>% 
              mutate(cancer=str_to_upper(cprefix))
dbWriteTable(conn, 'cnv_main',tmp)

#append to existing for tcga
for (study in studies$studyId){
  cprefix <- str_split_i(study,'_',1)
  tmp <- getDataByGenes(cBioPortal(), 
                        genes = igenes, 
                        molecularProfileIds = paste0(cprefix,'_tcga_pan_can_atlas_2018_gistic'),  
                        sampleListId = paste0(cprefix,'_tcga_pan_can_atlas_2018_all'))
  tmp <- tmp[[1]] %>% select(uniquePatientKey, hugoGeneSymbol, value) %>% 
    rename(., value='alteration', hugoGeneSymbol='gene') %>% 
    mutate(cancer=str_to_upper(cprefix))
  dbAppendTable(conn, 'cnv_main',tmp)
  print(study)
}

#add origimed data
tmps <- getDataByGenes(cBioPortal(),
                       genes = igenes,
                       molecularProfileIds = 'pan_origimed_2020_cna',
                       sampleListId = 'pan_origimed_2020_all')
clindat <- clinicalData(cBioPortal(),'pan_origimed_2020')

#translate their annotation to TCGA format UGH
translate <-c(
  'Colorectal Carcinoma'='COADREAD',
  'Liver Hepatocellular Carcinoma'='HCC',
  'Soft Tissue Sarcoma'='SOFT_TISSUE',
  'Carcinoma of Uterine Cervix'='UCS',
  'Urothelial Carcinoma'='BLCA',
  'Thymic Tumor'='THYM',
  'Kidney Renal Cell Carcinoma'='CCRCC',
  'Ovarian Carcinoma'='OV',
  'Gastric Cancer'='STAD',
  'Pancreatic Cancer'='PAAD',
  'Breast Carcinoma'='BRCA',
  'Esophageal Carcinoma'='ESCA',
  'Small Cell Lung Cancer'='LUSC',
  'Head and Neck Carcinoma'='HNSC',
  'Uterine Corpus Endometrial Carcinoma'='UCES',
  'Non Small Cell Lung Cancer'='LUAD',
  'Thyroid Carcinoma'='THYM')

tmp <- tmps$pan_origimed_2020_cna %>% 
  merge(., clindat %>% select(patientId, CANCER_TYPE), by='patientId') %>% 
  select(uniquePatientKey, hugoGeneSymbol, CANCER_TYPE, value) %>% 
  mutate(cancer=recode(CANCER_TYPE,!!!translate)) %>% 
  select(!CANCER_TYPE) %>% 
  rename(value='alteration', hugoGeneSymbol='gene')

dbAppendTable(conn, 'cnv_main', tmp)




#####
#SNPS



#create table with first round of data
study <- studies$studyId[1]
studies <- studies[2:nrow(studies),]


cprefix <- str_split_i(study,'_',1)
tmp <- getDataByGenes(cBioPortal(), 
                      genes = igenes, 
                      molecularProfileIds = paste0(cprefix,'_tcga_pan_can_atlas_2018_mutations'),  
                      sampleListId = paste0(cprefix,'_tcga_pan_can_atlas_2018_cnaseq'))
tmp <- tmp[[1]] %>% 
  filter(variantType=='SNP') %>% 
  select(uniquePatientKey, hugoGeneSymbol, referenceAllele, variantAllele, chr, startPosition, endPosition) %>% 
  rename(., hugoGeneSymbol='gene') %>% 
  mutate(cancer=str_to_upper(cprefix))
dbWriteTable(conn, 'cnv_snp',tmp)

#append to existing for tcga
for (study in studies$studyId){
  cprefix <- str_split_i(study,'_',1)
  tmp <- getDataByGenes(cBioPortal(), 
                        genes= igenes, 
                        molecularProfileIds = paste0(cprefix,'_tcga_pan_can_atlas_2018_mutations'),  
                        sampleListId = paste0(cprefix,'_tcga_pan_can_atlas_2018_cnaseq'))
  tmp <- tmp[[1]] %>% 
    filter(variantType=='SNP') %>% 
    select(uniquePatientKey, hugoGeneSymbol, referenceAllele, variantAllele, chr, startPosition, endPosition) %>%  
    rename(., hugoGeneSymbol='gene') %>% 
    mutate(cancer=str_to_upper(cprefix))
  dbAppendTable(conn, 'cnv_snp',tmp)
  print(study)
}


#append origimed data
tmp <-  getDataByGenes(cBioPortal(),
                       genes = igenes,
                       molecularProfileIds = 'pan_origimed_2020_mutations',
                       sampleListId = 'pan_origimed_2020_all')
clindat <- clinicalData(cBioPortal(),'pan_origimed_2020')

#translate their annotation to TCGA format UGH
translate <-c(
  'Colorectal Carcinoma'='COADREAD',
  'Liver Hepatocellular Carcinoma'='HCC',
  'Soft Tissue Sarcoma'='SOFT_TISSUE',
  'Carcinoma of Uterine Cervix'='UCS',
  'Urothelial Carcinoma'='BLCA',
  'Thymic Tumor'='THYM',
  'Kidney Renal Cell Carcinoma'='CCRCC',
  'Ovarian Carcinoma'='OV',
  'Gastric Cancer'='STAD',
  'Pancreatic Cancer'='PAAD',
  'Breast Carcinoma'='BRCA',
  'Esophageal Carcinoma'='ESCA',
  'Small Cell Lung Cancer'='LUSC',
  'Head and Neck Carcinoma'='HNSC',
  'Uterine Corpus Endometrial Carcinoma'='UCES',
  'Non Small Cell Lung Cancer'='LUAD',
  'Thyroid Carcinoma'='THYM')

tmp <- tmp$pan_origimed_2020_mutations %>% 
  merge(., clindat %>% select(patientId, CANCER_TYPE), by='patientId') %>% 
  select(uniquePatientKey, hugoGeneSymbol, referenceAllele, variantAllele, chr, startPosition, endPosition, CANCER_TYPE) %>% 
  mutate(cancer=recode(CANCER_TYPE,!!!translate)) %>% 
  select(!CANCER_TYPE) %>% 
  rename(hugoGeneSymbol='gene')

dbAppendTable(conn, 'cnv_snp', tmp)



#####
# create meta table with sample numbers 
# except this time, we need an index of each cancer-gene pair


tdat <- tbl(conn, 'cnv_main') %>% 
  select(uniquePatientKey, cancer, gene) %>% 
  distinct() %>% 
  group_by(cancer, gene) %>% 
  summarize(nsamp=n()) %>%  
  as.data.frame()

dbWriteTable(conn, 'cnv_nsamps', tdat, overwrite = T)

###
#genie?

#cbioportal can't directly serve the data and boy did i try
#but I finally found somewhere that will.

#library(synapser)

#synapse honestly makes it easier to use the web interface after all 
# :| 
setwd('WorkForaging/technium/cnv_calculator/')

#data
cna_genie <- read_delim('data/data_CNA.txt', delim='\t')

#meta

#####
# cnv_main
meta_cna <- read_delim('data/case_lists/cases_cnaseq.txt', '\t', skip = 3, col_names = c('ID'), trim_ws = T)
meta_cna <- meta_cna %>% 
  mutate(ID=case_when(
    grepl('case_list_ids', ID) ~ gsub('case_list_ids: ','',ID),
    T ~ ID
  ))
meta_cna <- meta_cna[2:nrow(meta_cna),] #remove header

#cancer subtypes
cfiles <- list.files('data/case_lists/')[-c(5,26)]

translate <-c(
  'Adrenocortical_Carcinoma'='ACC',
  'Bladder_Cancer'='BLCA',
  'Breast_Cancer'='BRCA',
  'Kidney_Cancer'='CCRCC',
  'Cervical_Cancer'='CESC',
  'Colorectal_Cancer'='COADREAD',
  'Esophagogastric_Cancer'='ESCA',
  'Gastrointestinal_Neuroendocrine_Tumor'='Gastrointestinal Neuroendocrine Tumor',
  'Hepatobiliary_Cancer'='HCC',
  'Head_and_Neck_Cancer'='HNSC',
  'Kidney_Cancer'='KIRC',
  'Leukemia'='LAML',
  'Liver_Cancer'='LIHC',
  'Lung_Cancer'='LUAD',
  'Melanoma'='SKCM',
  'Mesothelioma'='MESO',
  'Ovarian_Cancer'='OV',
  'Pancreatic_Cancer'='PAAD',
  'Pheochromocytoma'='PCPG',
  'Prostate_Cancer'='PRAD',
  'Soft_Tissue_Cancer'='SOFT_TISSUE',
  'Thymic_Tumor'='THYM',
  'Thyroid_Cancer'='THCA',
  'Endometrial_Cancer'='UCEC',
  'Uterine_Sarcoma'='UCS'
  )

for (cf in cfiles){
  genie_cname <- str_split(cf, 'cases_|\\.txt')[[1]][2]
  
  samples <- read_delim(paste0('data/case_lists/', cf), '\t', skip = 3, col_names = c('ID'), trim_ws = T)
  samples <- samples %>% mutate(ID=case_when(
                                  grepl('case_list_ids', ID) ~ gsub('case_list_ids: ','',ID),
                                  T ~ ID
                                  ))
  samples <- samples[2:nrow(samples),]
  samples$cancer <- genie_cname
  
  meta_cna <- merge(meta_cna, samples, by='ID', all.x = T)
  if ('cancer.x' %in% colnames(meta_cna)){
    meta_cna <- meta_cna %>% mutate(cancer.x = case_when(!is.na(cancer.y) ~ cancer.y, 
                      T ~ cancer.x
                      ))
    meta_cna <-  meta_cna %>% select(-cancer.y) %>% rename('cancer'=cancer.x)
  }
}

meta_cna <- meta_cna %>% 
        mutate(cancer=recode(cancer,!!!translate)) %>% 
        select(ID, cancer) %>% 
        filter(!is.na(cancer))

#some genie IDs are not in the data??
gsamples <- intersect(meta_cna$ID, colnames(cna_genie))

#need impact 505 genes
genes <- read_delim('data/data_gene_panel_MSK-IMPACT505.txt', '\t', skip = 2, col_names = c('Hugo_Symbol'), trim_ws = T)
genes <- t(genes[1,])[2:length(genes)] #read in as a list


#before write to db, need to pivot to long format, merge with meta, filter for genes, and rename columns
cna_genie <- cna_genie %>% 
  select(gsamples, 'Hugo_Symbol') %>%
  filter(Hugo_Symbol %in% genes) %>%
  pivot_longer(cols = -c(Hugo_Symbol), values_to = 'alteration', names_to = 'ID') %>% 
  merge(., meta_cna, by='ID') %>%
  filter(!is.na(alteration)) %>% 
  rename(gene='Hugo_Symbol', uniquePatientKey='ID') %>%
  mutate(cancer=recode(cancer,!!!translate), alteration=as.integer(alteration))

#write to db
dbAppendTable(conn, 'cnv_main', cna_genie)

###
#SNPS

#data
snp_genie <- read_delim('data/data_mutations_extended.txt', delim='\t')
snp_genie <- snp_genie %>% 
  filter(Variant_Type=='SNP') %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode, Reference_Allele, Tumor_Seq_Allele1, Chromosome, Start_Position, End_Position) %>% 
  rename(ID='Tumor_Sample_Barcode') %>% 
  merge(., meta_cna, by='ID') %>%
  rename(uniquePatientKey='ID', gene='Hugo_Symbol', referenceAllele='Reference_Allele', variantAllele='Tumor_Seq_Allele1', chr='Chromosome', startPosition='Start_Position', endPosition='End_Position') %>% 
  mutate(cancer=recode(cancer,!!!translate))

#write to db  
dbAppendTable(conn, 'cnv_snp', snp_genie)

