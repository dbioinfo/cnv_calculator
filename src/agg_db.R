library(tidyverse)
library(data.table)
library(Matrix)
#####
#single cna

t1<-readRDS('data/tcga_pancan_cna.rds')
t2<-readRDS('data/origimed_pancan_cna.rds')

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

#merge two datasets  
t2 <- t2 %>% mutate(cancer_type = recode(cancer_type, !!!translate))
tdat <- rbind(t1, t2 %>% filter(cancer_type %in% t1$cancer_type))
tdat <- tdat %>% 
  group_by(cancer_type, hugoGeneSymbol, alteration) %>% 
  summarise(n=sum(n),total_samp=sum(unique(total_samp)), rel_incidence=n/total_samp) %>% 
  ungroup()

#save
saveRDS(tdat, file='data/cna_single_agg.rds')


#####
#co-occurrence cna

#holy smokes that's a lot of data
#we're gonna need some sparse matrix addition
add_matrices <- function(a) { #efficient matrix addition 
  #a <- list(...)
  cols <- sort(unique(unlist(lapply(a, colnames))))
  rows <- sort(unique(unlist(lapply(a, rownames))))
  nrows <- length(rows)
  ncols <- length(cols)
  newms <- lapply(a, function(m) {
    s <- summary(m)
    i <- match(rownames(m), rows)[s$i]
    j <- match(colnames(m), cols)[s$j]
    ilj <- i<j
    sparseMatrix(i=ifelse(ilj, i, j),
                 j=ifelse(ilj, j, i),
                 x=s$x,
                 dims=c(nrows, ncols),
                 dimnames=list(rows, cols), symmetric=TRUE)
  })
  Reduce(`+`, newms)
}


#####
#seq read in files, add cnv count to make freq matrix
#change filename for each part until tdupmax and tdelmax are full
#so we don't kill the session by overloading RAM
#32Gb RAM limit on test server
f1 <- readRDS('data/cnvco_part3.rds')
f1 <- unlist(f1, recursive = T)
dupnames <- names(f1)[grep('dups', names(f1))]
delnames <- names(f1)[grep('dels', names(f1))]


#chunk dups and dels into RAM, since it's a sum, no problem
tdupmax <- add_matrices(f1[dupnames[1:5]])
t1dupmax <- add_matrices(f1[dupnames[6:10]])
t0dupmax<-readRDS('data/temp_dup_co.rds')
tdupmax <- add_matrices(list(tdupmax, t1dupmax, t0dupmax))
#savepoint
saveRDS(tdupmax, file='data/temp_dup_co.rds')

rm('t1dupmax', 'tdupmax', 't0dupmax')#save some space

tdelmax <- add_matrices(f1[delnames[1:5]])
t1delmax <- add_matrices(f1[delnames[6:10]])
t0delmax <- readRDS('data/temp_del_co.rds')
tdelmax <- add_matrices(list(tdelmax, t1delmax, t0dupmax))
saveRDS(tdelmax, file='data/temp_del_co.rds')

rm('t1delmax', 'tdelmax', 't0delmax')#save some space



#####
#put it all together for the server

p1 <- readRDS('data/cna_single_agg.rds')
p2 <- readRDS('data/temp_dup_co.rds')
p3 <- readRDS('data/temp_del_co.rds')
saveRDS(list(single_cna=p1,dup_co=p2,del_co=p3),'data/cnv.v1.rds')
