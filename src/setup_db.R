library(tidyverse)
library(dbplyr)
library(RPostgreSQL)

options(dplyr.summarise.inform = FALSE)


conn <- dbConnect(drv = RPostgres::Postgres(),
                  dbname='public_cnv', 
                  host='cnv_db', 
                  user='postgres',
                  password='password', 
                  port=5432)

#####
#create single_cnv table with rel_incidence and nsamp
#summarize each alteration within cancer types and genes
single_cnv <- tbl(conn, 'cnv_main') %>% 
  group_by(cancer, gene, alteration) %>% summarise(freq=n()) %>% 
  mutate(nsamp=sum(freq), rel_incidence=freq/nsamp)
dbWriteTable(conn, 'single_cnv', as.data.frame(single_cnv), overwrite=T)

#this table turns out a little smaller, so no partition necessary, only index
# on alteration, cancer, gene (in that order!)
dbGetQuery(conn, "CREATE INDEX single_cnv_idx ON single_cnv USING btree (alteration, cancer, gene)")

#####
#create co_cnv table with some extra rules
#like adding all combinations of alterations (no combining with 0)
#partition on cancer index on alt1, alt2, g1, g2

#this one's gonna be like tens of Gb probably, so iter over cancer types, genes, and optimize
cancers <- tbl(conn, 'single_cnv') %>% select(cancer) %>% distinct() %>% pull(cancer)
icancer <- cancers[1]
cancers <- cancers[2:length(cancers)]

genes <- tbl(conn, 'cnv_main') %>% filter(alteration!=0) %>% select(gene) %>% distinct() %>% pull(gene)
igene <- genes[1]
genes <- genes[2:length(genes)]

alts <- c(-2, -1, 1, 2) #no need to count 0s
ialt <- -1

db_init <- T

for (icancer in cancers){ #iter cancers
for (igene in genes){ #iter genes
  for (ialt in alts){ #each alteration needs its own relation 
    #first gather the patients who have this gene altered in this way
    alt_patients <- tbl(conn, 'cnv_main') %>% 
      filter(cancer==icancer, alteration==ialt, gene==igene) %>% 
      distinct() %>% 
      pull(uniquePatientKey)  
    
    #second gather those patients' co-altered genes
    tmp_co <- tbl(conn, 'cnv_main') %>% 
      filter(uniquePatientKey %in% alt_patients, gene!=igene, alteration!=0) %>% 
      group_by(alteration, gene) %>% 
      summarise(freq=n())
    
    #lastly, put it all together and write it to the table
    out <- data.frame(tmp_co) %>% 
      rename(gene2='gene', alteration2='alteration') %>% 
      mutate(gene1=igene, alteration1=ialt, cancer=icancer)
    
    #if the db needs to be built from start, initialize with db_init <- T
    if (db_init){
      dbWriteTable(conn, 'co_cnv', out, overwrite=T)
      db_init <- F
    }else {dbAppendTable(conn, 'co_cnv', out)}
    
    
  }  
  
  if (igene==genes[floor(length(genes)/2)]){
    print(paste0(icancer, ' 50% Uploaded'))
  }
  
  #print(igene)
  } 
print(paste0(icancer, ' 100% Uploaded'))
}

dbGetQuery(conn, "CREATE INDEX co_cnv_idx ON co_cnv USING btree(alteration1, alteration2, cancer, gene1, gene2)")
