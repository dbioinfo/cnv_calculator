library(RPostgreSQL)
library(tidyverse)

#read in data
data <- readRDS('~/WorkForaging/technium/cnv_calculator/data/cna_single_agg.rds')
codupdata <- readRDS('~/WorkForaging/technium/cnv_calculator/data/temp_dup_co.rds')
codeldata <- readRDS('~/WorkForaging/technium/cnv_calculator/data/temp_dup_co.rds')

#open connection to (active!) sql server
drv <- dbDriver('PostgreSQL')
conn <- dbConnect(drv = RPostgres::Postgres(),
                  dbname='postgres', 
                  host='127.0.0.1', 
                  user='postgres',
                  password='password', 
                  port=5432)

#write main db contents into a table
dbWriteTable(conn,'cnv_main',value = data)

# write co-cnv data into a table after reading the data from online into chunks
data <- readRDS('~/WorkForaging/technium/cnv_calculator/data/cnvco_part5.rds')
#create table if need be
#dbWriteTable(conn, 'co_cnv', value=out)

subtype_names <- names(data)
alterations <- c('dups','dels')

for (sname in subtype_names){
  for (alt in alterations){
    summ <- summary(data[[sname]][[alt]]) #pivot data for format
    icancer <- stringr::str_to_upper(strsplit(sname, split='_')[[1]][1])
    out <- data.frame(Gene1=rownames(data[[sname]][[alt]])[summ$i], 
                      Gene2=colnames(data[[sname]][[alt]])[summ$j], 
                      nsamp=summ$x,
                      cancer_subtype=icancer,
                      alteration=stringr::str_sub(alt, end = -2))
    dbAppendTable(conn, 'co_cnv', value=out)
  }
}

#dbWriteTable(conn, 'dup_co_cnv', value = tibble(codupdata))
