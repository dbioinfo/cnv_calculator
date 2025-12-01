library(tidyverse)
library(DBI)

###
# this script is used to update the database with new statistics
# it is run after the data has been updated with src/download_db.R



# Connect to the database
conn <- dbConnect(drv = RPostgres::Postgres(),
                  dbname='postgres', 
                  host='127.0.0.1', 
                  user='postgres',
                  password='password', 
                  port=5432)
"
# Compute nsamp statistics
cnv_nsamp <- tbl(conn, 'cnv_main') %>% 
  group_by(cancer) %>%
  summarise(nsamp = n_distinct(uniquePatientKey)) %>% 
  mutate(alteration='CNV')

snp_nsamp <- tbl(conn, 'cnv_snp') %>% 
  group_by(cancer) %>%
  summarise(nsamp = n_distinct(uniquePatientKey)) %>% 
  mutate(alteration='SNP')

nsamp <- rbind(as.data.frame(cnv_nsamp), as.data.frame(snp_nsamp))

# overwrite table in db
dbWriteTable(conn, 'cnv_nsamps', nsamp, overwrite = TRUE)
"

nsamp <- tbl(conn, 'cnv_nsamps')

# Compute single cnv statistics
tdat <- tbl(conn, "cnv_main") %>% 
  group_by(cancer, gene, alteration) %>%
  summarise(freq = n()) %>%
  merge(., nsamp, by = c('cancer', 'gene') ) %>%
  mutate(rel_incidence = freq/nsamp) 


# overwrite table in db
dbWriteTable(conn, "single_cnv", tdat, overwrite = TRUE)

#this table turns out a little smaller, so no partition necessary, only index
# on alteration, cancer, gene (in that order!)
# also delete if exists
dbExecute(conn, "DROP INDEX IF EXISTS single_cnv_idx")
dbExecute(conn, "CREATE INDEX single_cnv_idx ON single_cnv USING btree (alteration, cancer, gene)")

# Compute cooccurence cnv statistics

#iter through each cancer type, and get the cooccurence of each pair of genes
cancers <- tbl(conn, 'single_cnv') %>% select(cancer) %>% distinct() %>% pull(cancer)

genes <- tbl(conn, 'cnv_main') %>% filter(alteration!=0) %>% select(gene) %>% distinct() %>% pull(gene)

alts <- c(-2, -1, 1, 2) #no need to count 0s

db_init <- T

for (icancer in cancers){ #iter cancers
  for (igene in genes){ #iter genes
    for (ialt in alts){ #each alteration needs its own relation 
      #first gather the patients who have this gene altered in this way
      alt_patients <- tbl(conn, 'cnv_main') %>% 
        filter(cancer==icancer, alteration==ialt, gene==igene) %>% 
        distinct() %>% 
        pull(uniquePatientKey)  
      
      #if alt_patients is empty, skip to next alt
      if (length(alt_patients)==0){next}
      
      #second gather those patients' co-altered genes
      tmp_co <- tbl(conn, 'cnv_main') %>% 
        filter(uniquePatientKey %in% alt_patients, gene!=igene, alteration!=0) %>% 
        group_by(alteration, gene) %>% 
        summarise(freq=n())
      
      #lastly, put it all together and write it to the table
      out <- data.frame(tmp_co) %>% 
        rename(gene ='gene2', alteration='alteration2') %>% 
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

#this table is large, so index by cancer, gene1, gene2, alteration1, alteration2
dbExecute(conn, "DROP INDEX IF EXISTS co_cnv_idx")
dbExecute(conn, "CREATE INDEX co_cnv_idx ON co_cnv USING btree ( cancer,  gene1,  gene2, alteration1, alteration2)")

#add index to snps
dbExecute(conn, "DROP INDEX IF EXISTS cnv_snp_idx")
dbExecute(conn, 'CREATE INDEX cnv_snp_idx ON cnv_snp USING btree ("uniquePatientKey", cancer, gene)')
