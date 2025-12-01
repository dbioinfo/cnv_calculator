
# Define server logic
server <- shinyServer(function(input, output, session) {
  conn <- reactiveVal()
  tdat <- reactiveVal()
  dels <-reactiveVal()
  dups <-reactiveVal()
  genes <- reactiveVal()
  cancers <- reactiveVal()
  seer <- reactiveVal()
  igene_hist <- reactiveVal()
  meta_graph <- reactiveVal()
  g1_alt <- reactiveVal()
  g2_alt <- reactiveVal()
  snpsub <- reactiveVal()
  tbl_out <- reactiveVal()
  co_tbl <- reactiveVal()
  single_tbl <- reactiveVal()
  
  
  
observeEvent(input$load_data,{
    #input data
    print('connecting to server')
    conn <<- dbConnect(drv = RPostgres::Postgres(),
                      dbname='public_cnv',  #public_cnv
                      host='cnv_db',  #cnv_db
                      user='postgres',
                      password='password', 
                      port=5432)
    tdat <<- tbl(conn, 'single_cnv') %>% mutate(CNV=case_when(alteration<0~'Deletion', alteration>0~'Duplication'))
    
   
    
    #match with SEER db
    seer <<- read_csv('/home/dylan/WorkForaging/technium/cnv_calculator/data/seerstats.csv') %>% #just /data/seertats.csv in prod
            mutate(cancer=cancer_type)
    
    tdat <<- tdat %>% merge(., seer, by='cancer')
    
    genes <<- tdat %>% select(gene) %>% unique() %>% pull(gene)
    cancers <<- seer %>% pull(cancer) %>% unique()
    
    #populate some selecters and lists
    output$gene_select <- renderUI({
      selectizeInput(inputId='igene',
                         'Choose Gene',
                         choices = genes,
                          selected='MTAP')
    })
    output$gene_co_list <- renderUI({
      selectizeInput(inputId='igene_co',
                     'Choose Gene',
                     choices = genes, 
                     selected='MTAP')
    })
    
    output$g1_co_list <- renderUI({
      selectizeInput(inputId='g1_igene_co',
                     'Choose Gene 1',
                     choices = genes, 
                     selected='MTAP')
    })
    
    output$g2_co_list <- renderUI({
      selectizeInput(inputId='g2_igene_co',
                     'Choose Gene 2',
                     choices = genes, 
                     selected='CDKN2A')
    })
    
    output$g1_tbl_list <- renderUI({
      selectizeInput(inputId='g1_tbl_igene',
                     'Choose Target Gene',
                     choices = genes, 
                     selected='MTAP')
    })
    
    output$g2_tbl_list <- renderUI({
      selectizeInput(inputId='g2_tbl_igene',
                     'Choose Gene 2',
                     choices = genes, 
                     selected='CDKN2A')
    })
    
    output$cancer_list_tbl <- renderUI({
      selectInput(inputId = 'tbl_icancer',
                  'Choose Cancer Type',
                  choices = cancers,
                  selected = 'ACC')
    })
    
    output$zygosity_single <- renderUI({
      radioButtons('zygosity_single', choices = c('True','False'), label = 'Show Zygosity?', selected = 'False')
    })
    
    output$alteration_single <- renderUI({
      radioButtons('alteration_single', choices = c('CNV', 'SNP'), label = 'Alteration', selected = 'CNV')
    })
    
   
    ### download handlers
    output$download_igene_hist <- downloadHandler(
      filename = function() { paste(input$igene, '_cnv_incidence.png', sep='') },
      content = function(file) {
        ggsave(file, plot = igene_hist, device = "png", width = 4000, height=3000, units='px')
      }
    )
    
    output$download_co_pair <- downloadHandler(
      filename = function() { paste(input$igene, '_cnv_co_incidence.png', sep='') },
      content = function(file) {
        ggsave(file, plot = geneco_graph_pair, device = "png", width = 4000, height=3000, units='px')
      }
    )
    
    output$download_co_single <- downloadHandler(
      filename = function() { paste(input$igene, '_cnv_by_cancer_incidence.png', sep='') },
      content = function(file) {
        ggsave(file, plot = geneco_graph_single, device = "png", width = 4000, height=3000, units='px')
      }
    )
    
    output$download_meta_graph <- downloadHandler(
      filename = function() { paste(input$sel_meta, '.png', sep='') },
      content = function(file) {
        ggsave(file, plot = meta_graph, device = "png", width = 4000, height=3000, units='px')
      }
    )
    
    output$download_cotable <- downloadHandler(
      filename = function(){
        if (input$co_gopt == 'gsingle'){
          paste(input$g1_igene_co, input$g1_alt, 'cotable.csv', sep='_') 
        } else if (input$co_gopt == 'gpair') {
          paste(input$g1_igene_co, input$g1_alt, input$g2_igene_co, input$g2_alt, 'cotable.csv', sep='_')
        }
        
        }, 
      content = function(fname){
          write.csv(co_tbl, fname)
      }
    )
    
    output$download_single_tbl <- downloadHandler(
      filename = function(){
        paste(input$igene, '_cnv_by_cancer_type.csv')
      }, 
      content = function(fname){
        write.csv(top_cnv, fname)
      }
    )
    
    output$download_tbl <- downloadHandler(
      filename = function(){"tbl_results.csv"}, 
      content = function(fname){
        write.csv(tbl_out, fname)
      }
    )
    
    
    #plot data
    print('rendering single cnv hist...')
    output$igene_hist <- renderPlot({
      req(input$igene)
      igene <- input$igene
      
      print(igene)
      
      
      if (input$alteration_single=='CNV'){
        if (input$zygosity_single=='False'){ #meaning we do NOT show zygosity
          gdat <- tdat %>% 
            filter(gene==igene, alteration!=0) %>% 
            group_by(cancer, CNV) %>% 
            summarise(rel_incidence=sum(rel_incidence), total_incidence=max(total_incidence)) %>% 
            mutate(text=round(rel_incidence*total_incidence, 0)) 
        
          gg <- ggplot(gdat)+
            geom_bar(mapping=aes(x=cancer, y=rel_incidence, fill=CNV), stat='identity', position = 'dodge')+
            geom_text(mapping=aes(x=cancer, y=rel_incidence, group=CNV, label=text), 
                  position = position_dodge(width=1),
                  vjust=-0.5)+
            scale_fill_manual(values = c('#f37a98','#7fb568'))+
            scale_y_continuous(expand=c(0,0), limits = c(0, 1))+
            theme_bw()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                text = element_text(family='Arial', size = 22),
                plot.title = element_text(hjust = 0.5))+
            xlab('')+
            ylab('Proportion of Samples')+
            labs(title = (paste0('Relative CNV Incidence - ', igene)))
        } else if (input$zygosity_single=='True'){
          gdat <- tdat %>% 
            filter(gene==igene, alteration!=0) %>% 
            group_by(cancer, alteration) %>% 
            summarise(rel_incidence=sum(rel_incidence), total_incidence=max(total_incidence)) %>% 
            mutate(text=round(rel_incidence*total_incidence, 0), alteration=factor(alteration)) 
          
          
          gg <- ggplot(gdat)+
            geom_bar(mapping=aes(x=cancer, y=rel_incidence, fill=alteration), stat='identity', position = 'dodge')+
            geom_text(mapping=aes(x=cancer, y=rel_incidence, group=alteration, label=text), 
                      position = position_dodge(width=1),
                      vjust=-0.5)+
            scale_fill_manual(values = c('#f26161','#f29696','#8abf88','#41b33b'))+
            scale_y_continuous(expand=c(0,0), limits = c(0, 1))+
            theme_bw()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                  text = element_text(family='Arial', size = 22),
                  plot.title = element_text(hjust = 0.5))+
            xlab('')+
            ylab('Proportion of Samples')+
            labs(title = (paste0('Relative CNV Incidence - ', igene)))+
            guides(fill=guide_legend(title="Alteration"))
          
        }
      } else if(input$alteration_single=='SNP') {
        
        
        snpdat <- tbl(conn, 'cnv_snp') %>% 
          group_by(cancer, gene, uniquePatientKey) %>% 
          distinct()
        
        snpdat <- snpdat %>% 
          group_by(cancer, gene) %>% 
          summarize(freq=n())
        
        print(snpdat)
        
        nsamps <- tbl(conn, 'cnv_nsamps') %>%
          filter(alteration=='SNP') %>% 
          select(cancer, nsamp) %>% 
          group_by(cancer) %>% 
          summarize(nsamp=max(nsamp))
        
        snpdat <- merge(snpdat, nsamps, by='cancer')
        
        print(snpdat)
        
        snpdat <- merge(snpdat, seer, by='cancer') %>% 
          filter(gene==igene) %>% 
          mutate(rel_incidence = freq / nsamp,
                 text=round(rel_incidence*total_incidence, 0))
        
        print(snpdat)
        
        gg <- ggplot(snpdat)+
          geom_bar(mapping=aes(x=cancer, y=rel_incidence, fill=cancer), stat='identity', position = 'dodge')+
          geom_text(mapping=aes(x=cancer, y=rel_incidence, label=text), 
                    position = position_dodge(width=1),
                    vjust=-0.5)+
          #scale_fill_manual(values = c('#f26161','#f29696','#8abf88','#41b33b'))+
          scale_y_continuous(expand=c(0,0), limits = c(0, 1))+
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                text = element_text(family='Arial', size = 22),
                plot.title = element_text(hjust = 0.5))+
          xlab('')+
          ylab('Proportion of Samples')+
          labs(title = (paste0('Relative SNP Incidence - ', igene)))+
          guides(fill='none')
        
      }
        
        igene_hist <<- gg
        gg
      
      
    }, width=1000, height=800)
    
    #output table of top 100 most cnvs
    output$top_cnv <- DT::renderDataTable({
      igene <- input$igene
      
      if (input$zygosity_single=='False'){
        subset <- tdat %>% 
          filter(gene==igene, alteration!=0) %>% 
          group_by(cancer, CNV) %>% 
          summarise(rel_incidence=sum(rel_incidence), total_incidence=max(total_incidence)) %>% 
          arrange(-rel_incidence)
      } else {
        subset <- tdat %>% 
          filter(gene==igene, alteration!=0) %>% 
          group_by(cancer, alteration) %>% 
          summarise(rel_incidence=sum(rel_incidence), total_incidence=max(total_incidence)) %>% 
          arrange(-rel_incidence)
      }
      
      
      top_cnv <<- subset
      DT::datatable(subset,  list(pageLength = 10, scrollX=T), rownames = F)
      
    })
      
     
    
    
    })
  
  #metadata handler
observeEvent(input$load_meta, {
  
  print('rendering meta plot...')
  output$meta_graph <- renderPlot({

  if (input$sel_meta=='Cancer Incidence'){
    gg <- ggplot(seer)+
      geom_bar(mapping=aes(x=cancer,y=total_incidence,fill=cancer), stat='identity')+
      scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(nrow(seer)))+
      geom_text(mapping=aes(x=cancer,y=total_incidence, label=male), color='#3b91ed',vjust=-0.5, size=4)+
      geom_text(mapping=aes(x=cancer,y=total_incidence, label=female), color='#f542d7',vjust=-2, size=4)+
      theme_bw()+
      scale_y_continuous(trans='log2', breaks=c(10,100,1000,10000, 100000), labels = label_comma())+
      xlab('')+
      ylab('Incidence')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=22),
            axis.title = element_text(family='Arial', face='bold', size=22))+
      guides(fill='none')
  } else if (input$sel_meta=='Cancer CNV Samples') {
    gdat <- tbl(conn, 'single_cnv') %>% 
      select(cancer, alteration, freq, nsamp) %>%
      filter(alteration!=0) %>% 
      mutate(alteration=case_when(
        alteration<0 ~'Deletion',
        alteration>0 ~'Duplication',
        .default='NA'
      )) %>% 
      group_by(cancer, alteration) %>% 
      summarise(freq=sum(freq), nsamp=max(nsamp))
    
    
    
    #gdat <- gdat %>% filter(nsamp!=253)
    gg<-ggplot(gdat)+
      geom_bar(mapping=aes(x=cancer,y=freq,fill=alteration), stat='identity', position='dodge')+
      scale_fill_manual(values = c('#b5265f','#449e26'))+
      theme_bw()+
      xlab('')+
      ylab('N Samples')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=22),
            axis.title = element_text(family='Arial', face='bold', size=22))+
      guides(fill='none')
  }
   meta_graph<<- gg
   gg
    }, width=1000, height=800)
  })
  
observeEvent(input$load_cograph, {
    print('rendering co-alteration graph...')
    #print(input$co_gopt)
    #known gene pair graph
    if (input$co_gopt=='gpair'){
      
      req(input$g1_igene_co)
      req(input$g2_igene_co)
      print('geneco_graph_pair')
      
      igene_g1 <- input$g1_igene_co
      igene_g2 <- input$g2_igene_co
      
      if (input$g1_alt=='SNP'){ #input for gene 1 alt is SNP
        
        g1_patients <- tbl(conn, 'cnv_snp') %>% 
          filter(gene==igene_g1) %>% 
          select(uniquePatientKey) %>%
          distinct() %>% 
          as.data.frame()
        
        nsamps <- tbl(conn, 'cnv_nsamps') %>%
          filter(alteration=='SNP') %>%
          select(cancer, nsamp) %>% 
          group_by(cancer) %>% 
          summarize(nsamp=max(nsamp)) %>% 
          as.data.frame()
        
        if (input$g2_tbl=='CNV'){ #input for g2 alt is CNV
          
          output$geneco_graph_pair<- renderPlot({

            subset <- tbl(conn, 'cnv_main') %>%
              filter(uniquePatientKey %in% !! g1_patients$uniquePatientKey, alteration!=0, cancer %in% cancers) %>% 
              select(cancer, gene, alteration) %>% 
              group_by(cancer, gene, alteration) %>% 
              summarize(freq=n()) %>% 
              as.data.frame() %>% 
              merge(., nsamps, by='cancer') %>% 
              mutate(rel_incidence = freq / nsamp)
          
            
            if (input$zygosity_co=='False'){ #graph without zygosity

              subset <- subset %>% mutate(alteration=case_when(
                                          alteration<0~'Deletion',
                                          alteration>0~'Duplication'))
            
              fill_colors<-c('#f26161','#41b33b')
            
            
            } else if (input$zygosity_co=='True'){ #graph with zygosity
            
              fill_colors<-c('#f26161','#f29696','#8abf88','#41b33b')
            
            }
          
            gg<-ggplot(subset)+
              geom_bar(mapping=aes(x=reorder(cancer, desc(rel_incidence)), y=rel_incidence, fill=factor(alteration) ), stat='identity', position = 'dodge')+
              scale_y_continuous(expand = c(0,0),limits=c(0,1))+
              scale_fill_manual(values=fill_colors)+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                    text = element_text(family='Arial', size = 18),
                    plot.title = element_text(hjust = 0.5))+
              xlab('')+
              ylab(paste0('Proportion of ',igene_g1, ' SNP ',
                          '& ',igene_g2,' Co-occurences') )+
              labs(title=paste0(igene_g1, ' Mutation Co-Occurrence'),
                   subtitle = paste0(
                     paste0(igene_g1,' ',input$g1_alt,' & '),
                     paste0(igene_g2)
                   ))+ 
              guides(fill=guide_legend(title="Gene 2 CNV"))  
            
          
            gg  
          }, width = 1000, height=800)
        } else if (input$g2_tbl=='SNP'){ #input for g2 alt is also SNP
          output$geneco_graph_pair<- renderPlot({
          
            subset <- tbl(conn, 'cnv_snp') %>% 
              filter(uniquePatientKey %in% !! g1_patients$uniquePatientKey, cancer %in% cancers, gene==igene_g2) %>% 
              select(cancer, gene) %>%
              group_by(cancer, gene) %>% 
              summarize(freq = n()) %>% 
              as.data.frame() %>% 
              merge(., nsamps, by='cancer') %>% 
              mutate(rel_incidence = freq / nsamp)
          
          
            topngenes <- subset %>% arrange(desc(rel_incidence)) %>% 
                        head(n=25) %>% select(gene) %>% distinct() %>% as.data.frame()
            
            subset <- subset %>% 
              filter(gene %in% topngenes$gene) 
            
            print(subset)
            
            gg <- ggplot(subset) + 
              geom_bar(mapping=aes(x=reorder(cancer, desc(rel_incidence)), y=rel_incidence, fill=cancer), stat='identity', position='dodge')+
              scale_y_continuous(expand = c(0,0),limits=c(0,1))+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                    text = element_text(family='Arial', size = 18),
                    plot.title = element_text(hjust = 0.5))+
              xlab('')+
              ylab(paste0('Proportion of ',igene_g1, ' SNP ',
                          '& ',igene_g2,' SNP Co-occurences') )+
              labs(title=paste0(igene_g1, ' Mutation Co-Occurrence'),
                   subtitle = paste0(
                     paste0(igene_g1,' & '),
                     paste0(igene_g2, ' SNP Co-Occurrence')
                   ))+
              guides(fill='none')
            
            gg  
          }, width = 1000, height=800)
        }
      } else { #input for gene 1 alt is CNV
     
        if (input$g2_tbl=='CNV'){
          output$geneco_graph_pair<- renderPlot({
      
          req(input$g1_igene_co)
          req(input$g2_igene_co)
          print('geneco_graph_pair')
          
          if (input$g1_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_alt=='Deletion') {g1_alt <<- c(-1,-2)}
          igene_g1 <- input$g1_igene_co
          igene_g2 <- input$g2_igene_co
      
          if (input$zygosity_co=='False'){ #create graph with or without zygosity info
            #get igene cancer freq
            igene_freq <- tbl(conn, 'single_cnv') %>% 
              filter(gene==igene_g1, alteration %in% g1_alt) %>% 
              mutate(alteration=case_when(alteration<0~'Deletion', alteration>0~'Duplication')) %>% 
              select(cancer, freq, alteration) %>% 
              group_by(cancer) %>% 
              summarise(tfreq=sum(freq)) %>% 
              as_data_frame()
      
            #filter for subset of samples
            subset <- tbl(conn, 'co_cnv') %>% 
              filter(gene1==igene_g1, gene2==igene_g2, alteration1 %in% g1_alt) %>% 
              mutate(alteration2=case_when(alteration2<0~'Deletion', alteration2>0~'Duplication')) %>% 
              group_by(cancer, gene1, gene2, alteration2) %>% 
              summarise(freq = sum(freq)) %>% 
              as_data_frame() %>% 
              merge(igene_freq, by='cancer') %>% 
              mutate(rel_incidence = freq/tfreq)
        
           
            gg<-ggplot(subset)+
              geom_bar(mapping=aes(x=reorder(cancer, desc(rel_incidence)), y=rel_incidence, fill=factor(alteration2) ), stat='identity', position = 'dodge')+
              scale_y_continuous(expand = c(0,0))+
              scale_fill_manual(values=c('#f26161','#41b33b'))+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                 text = element_text(family='Arial', size = 18),
                 plot.title = element_text(hjust = 0.5))+
              xlab('')+
              ylab(paste0('Proportion of ',igene_g1, ' ', substr(input$g1_alt, 1, 3),
                       '& ',igene_g2,' Co-occurences') )+
              labs(title=paste0(igene_g1, ' Co-Occurrence'),
                subtitle = paste0(
                  paste0(igene_g1,' ',input$g1_alt,' & '),
                  paste0(igene_g2)
                ))+ 
              guides(fill=guide_legend(title="Gene 2 CNV"))
          } else if (input$zygosity_co=='True'){
            #get igene cancer freq
            igene_freq <- tbl(conn, 'single_cnv') %>% 
              filter(gene==igene_g1, alteration %in% g1_alt) %>% 
              select(cancer, freq, alteration) %>% 
              group_by(cancer) %>% 
              summarise(tfreq=sum(freq)) %>% 
              as_data_frame()
            
            #filter for subset of samples
            subset <- tbl(conn, 'co_cnv') %>% 
              filter(gene1==igene_g1, gene2==igene_g2, alteration1 %in% g1_alt) %>% 
              group_by(cancer, gene1, gene2, alteration2) %>% 
              summarise(freq = sum(freq)) %>% 
              as_data_frame() %>% 
              merge(igene_freq, by='cancer') %>% 
              mutate(rel_incidence = freq/tfreq)
            
            
            gg<-ggplot(subset)+
              geom_bar(mapping=aes(x=reorder(cancer, desc(rel_incidence)), y=rel_incidence, fill=factor(alteration2) ), stat='identity', position = 'dodge')+
              scale_y_continuous(expand = c(0,0))+
              scale_fill_manual(values=c('#f26161','#f29696','#8abf88','#41b33b'))+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                    text = element_text(family='Arial', size = 18),
                    plot.title = element_text(hjust = 0.5))+
              xlab('')+
              ylab(paste0('Proportion of ',igene_g1, ' ', substr(input$g1_alt, 1, 3),
                          '& ',igene_g2,' Co-occurences') )+
              labs(title=paste0(igene_g1, ' Co-Occurrence'),
                   subtitle = paste0(
                     paste0(igene_g1,' ',input$g1_alt,' & '),
                     paste0(igene_g2)
                   ))+ 
              guides(fill=guide_legend(title="Gene 2 CNV"))
          }
      
        gg
        }, width = 1000, height=800)
        print('co-table')
      
        output$cotable <- DT::renderDataTable({
          req(input$g1_igene_co)
          req(input$g2_igene_co)
          igene_g1 <- input$g1_igene_co
          igene_g2 <- input$g2_igene_co
        
          subset <- tbl(conn, 'co_cnv') %>% 
            filter(alteration1 %in% g1_alt,gene1==igene_g1,gene2==igene_g2) %>% 
            arrange(desc(freq)) %>% head(n=50)
          
          if (input$zygosity_co=='False'){
              subset <- subset %>% 
              mutate(alteration1=case_when(alteration1<0~'Deletion', alteration1>0~'Duplication')) %>% 
              mutate(alteration2=case_when(alteration2<0~'Deletion', alteration2>0~'Duplication')) 
          }
          
          co_tbl <<- data.frame(subset)
          DT::datatable(co_tbl %>% relocate(gene1, alteration1, gene2, alteration2), 
                      list(pageLength = 10, scrollX=T),rownames = F)
          
          
        
      })
      
      } else if (input$g2_tbl=='SNP'){
        output$geneco_graph_pair<- renderPlot({
          
          req(input$g1_igene_co)
          req(input$g2_igene_co)
          print('geneco_graph_pair_snp')
          if (input$g1_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_alt=='Deletion') {g1_alt <<- c(-1,-2)}
          igene_g1 <- input$g1_igene_co
          igene_g2 <- input$g2_igene_co
          
          #get igene cancer freq
          igene_freq <- tbl(conn, 'single_cnv') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            select(cancer, freq, alteration) %>% 
            group_by(cancer) %>% 
            summarise(tfreq=sum(freq)) %>% 
            as_data_frame()
          
          igene_samps <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            pull(uniquePatientKey)
          
          
          subset <- tbl(conn, 'cnv_snp') %>% 
            filter(uniquePatientKey %in% igene_samps, gene == igene_g2) %>% 
            select(uniquePatientKey, cancer, gene) %>% 
            distinct() %>% #gets rid of double counted snp (ie gene has 2+ mutations)
            group_by(cancer, gene) %>% 
            summarise(snpfreq=n()) %>% 
            merge(igene_freq, by='cancer') %>% 
            mutate(rel_incidence=snpfreq/tfreq) %>% 
            arrange(-rel_incidence) 
          
          snpsub <<- subset
          
          gg<-ggplot(subset)+
            geom_bar(mapping=aes(x=reorder(cancer, desc(rel_incidence)), y=rel_incidence, fill=cancer), stat='identity', position = 'dodge')+
            scale_y_continuous(expand = c(0,0))+
            theme_bw()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                  text = element_text(family='Arial', size = 18),
                  plot.title = element_text(hjust = 0.5))+
            xlab('')+
            ylab(paste0('Proportion of ',igene_g1, ' ', substr(input$g1_alt, 1, 3),
                        '& ',igene_g2,' SNP Co-occurences') )+
            labs(title=paste0(igene_g1, ' SNP Co-Occurrence'),
                 subtitle = paste0(
                   paste0(igene_g1,' ',input$g1_alt,' & '),
                   paste0(igene_g2)
                 ))
          
          gg
        }, width = 1000, height=800)
        
        output$cotable <<- DT::renderDataTable({
          req(input$g1_igene_co)
          req(input$g2_igene_co)
          igene_g1 <- input$g1_igene_co
          igene_g2 <- input$g2_igene_co
          co_tbl <<- data.frame(snpsub)
          DT::datatable(co_tbl %>% arrange(-rel_incidence), 
                        list(pageLength = 10, scrollX=T),rownames = F)
          
          
          
        })  
        
      }
     
    }
    ##single gene explorer
    } else if (input$co_gopt=='gsingle'){
      print('gsingle')
      
      req(input$igene_co)
      igene<- input$igene_co
      
      if (input$g1_alt=='SNP'){ #input for gene 1 alt is SNP
        print('snp')
        if (input$g2_alt!='SNP'){ #input for g2 alt is CNV
          print('snp_cnv')
          output$geneco_graph_single <- renderPlot({
          
            g_patients <- tbl(conn, 'cnv_snp') %>% 
              filter(gene==igene) %>% 
              select(uniquePatientKey) %>%
              distinct() %>% 
              as.data.frame()
            
            nsamps <- tbl(conn, 'cnv_nsamps') %>%
              select(cancer, nsamp) %>% 
              group_by(cancer) %>% 
              summarize(nsamp=max(nsamp)) %>% 
              as.data.frame()
            
            if (input$g1_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_alt=='Deletion') {g1_alt <<- c(-1,-2)}
            
            subset <- tbl(conn, 'cnv_main') %>%
              filter(uniquePatientKey %in% !! g_patients$uniquePatientKey, alteration %in% g1_alt, cancer %in% cancers) %>% 
              select(cancer, gene, alteration) %>% 
              group_by(cancer, gene, alteration) %>% 
              summarize(freq=n()) %>% 
              as.data.frame() %>% 
              merge(., nsamps, by='cancer') %>% 
              mutate(rel_incidence = freq / nsamp) 

            topngenes <- subset %>% arrange(desc(rel_incidence)) %>% 
              head(n=25) %>% select(gene) %>% distinct() %>% as.data.frame()
            
            subset <- subset %>% 
              filter(gene %in% topngenes$gene) 
            
            gg<-ggplot(data.frame(subset))+
              geom_bar(mapping=aes(x=reorder(gene, desc(rel_incidence)), y=rel_incidence, fill=cancer), stat='identity', position = 'dodge')+
              scale_y_continuous(expand = c(0,0))+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                    text = element_text(family='Arial', size = 18),
                    plot.title = element_text(hjust = 0.5))+
              xlab('')+
              ylab(paste0('Proportion of ',input$igene_co, ' Mutation ',' Co-occurences') )+
              labs(title=paste0(input$igene_co, ' Co-Occurrence'),
                   subtitle = paste0(
                     paste0(igene,' SNP & '),
                     paste0('Gene 2 ',input$g2_alt)
                   ))
            
            gg
          }, width = 1000, height=800)
          
        } else if (input$g2_alt=='SNP'){ #input for g2 alt is SNP
          print('SNP v SNP')
          output$geneco_graph_single <- renderPlot({
            g_patients <- tbl(conn, 'cnv_snp') %>% 
              filter(gene==igene) %>% 
              select(uniquePatientKey) %>%
              distinct() %>% 
              as.data.frame()
            
            nsamps <- tbl(conn, 'cnv_nsamps') %>%
              select(cancer, nsamp) %>% 
              group_by(cancer) %>% 
              summarize(nsamp=max(nsamp)) %>% 
              as.data.frame()
            
            subset <- tbl(conn, 'cnv_snp') %>%
              filter(uniquePatientKey %in% !! g_patients$uniquePatientKey, cancer %in% cancers) %>% 
              select(cancer, gene) %>% 
              group_by(cancer, gene) %>% 
              summarize(freq=n()) %>% 
              as.data.frame() %>% 
              merge(., nsamps, by='cancer') %>% 
              mutate(rel_incidence = freq / nsamp)
            print(cancers)
            topngenes <- subset %>% arrange(desc(rel_incidence)) %>% 
              head(n=25) %>% select(gene) %>% distinct() %>% as.data.frame()
            
            subset <- subset %>% 
              filter(gene %in% topngenes$gene) 
            
            gg<-ggplot(data.frame(subset))+
              geom_bar(mapping=aes(x=reorder(gene, desc(rel_incidence)), y=rel_incidence, fill=cancer), stat='identity', position = 'dodge')+
              scale_y_continuous(expand = c(0,0))+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                    text = element_text(family='Arial', size = 18),
                    plot.title = element_text(hjust = 0.5))+
              xlab('')+
              ylab(paste0('Proportion of ',input$igene_co, ' Mutation ',' Co-occurences') )+
              labs(title=paste0(input$igene_co, ' Co-Occurrence'),
                   subtitle =  paste0(igene,' SNP & Gene 2 SNP')
                   )
            
            gg
          }, width = 1000, height=800)
        }
      } else { #input for gene 1 alt is CNV
        print('cnv')
        
          if (input$g2_alt!='SNP'){
            output$geneco_graph_single <- renderPlot({
              
              if (input$g1_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_alt=='Deletion') {g1_alt <<- c(-1,-2)}
              if (input$g2_alt=='Duplication'){g2_alt <<- c(1,2)}else if (input$g2_alt=='Deletion') {g2_alt <<- c(-1,-2)}
              
              
              #need to get 30 most frequently co-duped genes
              co_alt_genes <- tbl(conn, 'co_cnv') %>% 
                filter(alteration1 %in% g1_alt,alteration2 %in% g2_alt,gene1==igene,gene2!=igene) %>% 
                group_by(gene2) %>% #group by gene 2, summing all cancer types into 1 measure
                summarise(tfreq = sum(freq)) %>% 
                as.data.frame() %>% 
                arrange(desc(tfreq)) %>% #arrange by sum
                head(n=15) %>%  #top 15
                pull(gene2)
              
              igene_freq <- tbl(conn, 'single_cnv') %>% 
                filter(gene==igene, alteration %in% g1_alt) %>% 
                select(cancer, freq, alteration) %>% 
                rename('tfreq'=freq, 'alteration1'=alteration) %>% 
                as_data_frame()
              
              #then get their most common cancer types
              subset <- tbl(conn, 'co_cnv') %>% 
                filter(alteration1 %in% g1_alt,alteration2 %in% g2_alt,gene1==igene,gene2 %in% co_alt_genes) %>% 
                arrange(desc(freq)) %>% head(n=150) %>% 
                merge(igene_freq, by=c('cancer', 'alteration1'))
              #print(subset)
              gg<-ggplot(data.frame(subset))+
                geom_bar(mapping=aes(x=reorder(gene2, desc(freq/tfreq)), y=freq/tfreq, fill=cancer), stat='identity', position = 'dodge')+
                scale_y_continuous(expand = c(0,0))+
                theme_bw()+
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                      text = element_text(family='Arial', size = 18),
                      plot.title = element_text(hjust = 0.5))+
                xlab('')+
                ylab(paste0('Proportion of ',input$igene_co, ' ', substr(input$g1_alt, 1, 3),' Co-occurences') )+
                labs(title=paste0(input$igene_co, ' Co-Occurrence'),
                     subtitle = paste0(
                       paste0(igene,' ',input$g1_alt,' & '),
                       paste0('Gene 2 ',input$g2_alt)
                     ))
              
              gg
            }, width = 1000, height=800)
            print('co-table')
            
            output$cotable <<- DT::renderDataTable({
              igene <- input$igene
              subset <- tbl(conn, 'co_cnv') %>% 
                filter(alteration1 %in% g1_alt,alteration2 %in% g2_alt,gene1==igene,gene2!=igene) %>% 
                arrange(desc(freq)) %>% head(n=50)
              
              co_tbl <<- data.frame(subset)
              DT::datatable(co_tbl %>% relocate(gene1, alteration1, gene2, alteration2), 
                            list(pageLength = 10, scrollX=T),rownames = F)
              
              
            })
          } else if (input$g2_alt=='SNP'){
            req(input$igene_co)
            igene<- input$igene_co
            print('snp')
            output$geneco_graph_single <- renderPlot({
              
              if (input$g1_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_alt=='Deletion') {g1_alt <<- c(-1,-2)}
            
              igene_samps <- tbl(conn, 'cnv_main') %>% 
                filter(gene==igene, alteration %in% g1_alt) %>% 
                pull(uniquePatientKey)
              
              igene_freq <- tbl(conn, 'single_cnv') %>% 
                filter(gene==igene, alteration %in% g1_alt) %>% 
                select(cancer, freq, alteration) %>%
                group_by(cancer) %>% 
                summarise(tfreq=sum(freq)) %>% 
                as_data_frame()
              
              selgenes <- tbl(conn, 'cnv_snp') %>% 
                filter(uniquePatientKey %in% igene_samps) %>% 
                select(uniquePatientKey, cancer, gene) %>% 
                distinct() %>% 
                group_by(cancer, gene) %>% 
                summarise(snpfreq=n()) %>% 
                merge(igene_freq, by='cancer') %>% 
                mutate(rel_incidence=snpfreq/tfreq) %>% 
                arrange(-rel_incidence) %>% 
                head(n=50) %>% 
                pull(gene) %>% 
                unique()
              
              subset <- tbl(conn, 'cnv_snp') %>% 
                filter(uniquePatientKey %in% igene_samps, gene %in% selgenes, cancer %in% cancers) %>% 
                select(uniquePatientKey, cancer, gene) %>% 
                distinct() %>% 
                group_by(cancer, gene) %>% 
                summarise(snpfreq=n()) %>% 
                merge(igene_freq, by='cancer') %>% 
                mutate(rel_incidence=snpfreq/tfreq) %>% 
                arrange(-rel_incidence) 
              
              snpsub <<- subset
              #print(head(subset))
              gg<-ggplot(data.frame(subset))+
                geom_bar(mapping=aes(x=reorder(gene, desc(rel_incidence)), y=rel_incidence, fill=cancer), stat='identity', position = 'dodge')+
                scale_y_continuous(expand = c(0,0))+
                theme_bw()+
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                      text = element_text(family='Arial', size = 18),
                      plot.title = element_text(hjust = 0.5))+
                xlab('')+
                ylab(paste0('Proportion of ',input$igene_co, ' ', substr(input$g1_alt, 1, 3),' SNP Co-occurences') )+
                labs(title=paste0(input$igene_co, ' Co-Occurrence'),
                     subtitle = paste0(
                       paste0(igene,' ',input$g1_alt,' & '),
                       paste0('Gene 2 SNP')
                     ))
              
              gg
            }, width = 1000, height=800)
            print('co-table')
            
            output$cotable <<- DT::renderDataTable({
              co_tbl <<- data.frame(snpsub)
              DT::datatable(co_tbl %>% arrange(-rel_incidence), 
                            list(pageLength = 10, scrollX=T),rownames = F)
              
              
              
            })
            
            
          }
          
          
        }#no indent
      

   }
  })
  
  observeEvent(input$load_tbl,{
    
    if (input$tbl_mode == 'Gene-Specific'){
        if (input$g2_tbl_alt!='SNP'){
          
          if (input$g1_tbl_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_tbl_alt=='Deletion') {g1_alt <<- c(-1,-2)}
          if (input$g2_tbl_alt=='Duplication'){g2_alt <<- c(1,2)}else if (input$g2_tbl_alt=='Deletion') {g2_alt <<- c(-1,-2)}
          
          igene_g1 <- input$g1_tbl_igene
          igene_g2 <- input$g2_tbl_igene
          alt1 <- input$g1_tbl_alt
          
          #calculate the incidence of the target gene alteration
          g1_single <- tbl(conn, 'single_cnv') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            group_by(cancer, gene, nsamp) %>% 
            summarise(g1_freq=sum(freq)) %>% 
            mutate(g1_rel_incidence=g1_freq/nsamp, alteration1=alt1) %>% 
            rename('gene1'=gene)
          
          
          #calculate the co-occuring incidence of gene 2
          out <- tbl(conn, 'co_cnv') %>% 
            filter(alteration1 %in% g1_alt,alteration2 %in% g2_alt,gene1==igene_g1,gene2==igene_g2) %>% 
            group_by(gene1, gene2, cancer, alteration2) %>% 
            summarise(co_freq=sum(freq)) %>% 
            merge(g1_single, by=c('cancer','gene1')) %>% 
            mutate(co_rel_incidence=co_freq/nsamp) %>% 
            merge(seer, by='cancer') %>% 
            mutate(est_co_incidence=round(co_rel_incidence*total_incidence, 0)) %>% 
            select(-cancer_type) %>% 
            relocate(cancer, 
                     gene1, alteration1, 
                     gene2, alteration2, 
                     g1_freq, nsamp, g1_rel_incidence,
                     co_freq, co_rel_incidence)
            
            
          
        } else if (input$g2_tbl_alt=='SNP'){
          
          if (input$g1_tbl_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_tbl_alt=='Deletion') {g1_alt <<- c(-1,-2)}
          
          igene_g1 <- input$g1_tbl_igene
          igene_g2 <- input$g2_tbl_igene
          alt1 <- input$g1_tbl_alt
          
          #calculate the incidence of the target gene alteration
          g1_single <- tbl(conn, 'single_cnv') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            group_by(cancer, gene, nsamp) %>% 
            summarise(g1_freq=sum(freq)) %>% 
            mutate(g1_rel_incidence=g1_freq/nsamp, alteration1=alt1) %>% 
            rename('gene1'=gene)
          
          g1_samps <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            pull(uniquePatientKey) %>% 
            unique()
          
          out <- tbl(conn, 'cnv_snp') %>% 
            filter(uniquePatientKey %in% g1_samps, gene==igene_g2) %>% 
            select(cancer, gene, uniquePatientKey) %>% 
            distinct() %>%
            group_by(cancer, gene) %>% 
            summarise(co_freq = n()) %>% 
            merge(g1_single, by='cancer') %>%
            merge(seer, by='cancer') %>% 
            mutate(co_rel_incidence=co_freq/nsamp, alteration2='SNP') %>%
            mutate(est_co_incidence=round(co_rel_incidence*total_incidence, 0)) %>% 
            rename('gene2'=gene) %>% 
            select(-cancer_type) %>% 
            relocate(cancer, 
                     gene1, alteration1, 
                     gene2, alteration2, 
                     g1_freq, nsamp, g1_rel_incidence,
                     co_freq, co_rel_incidence)
          
          
          
        }
       
    } else if (input$tbl_mode == 'Cancer-Specific'){
        if (input$g2_tbl_mode=='CNV'){
          
          if (input$g1_tbl_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_tbl_alt=='Deletion') {g1_alt <<- c(-1,-2)}
          
          igene_g1 <- input$g1_tbl_igene
          icancer <- input$tbl_icancer
          alt1 <- input$g1_tbl_alt
          
          #calculate the incidence of the target gene alteration
          g1_single <- tbl(conn, 'single_cnv') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            group_by(cancer, gene, nsamp) %>% 
            summarise(g1_freq=sum(freq)) %>% 
            mutate(g1_rel_incidence=g1_freq/nsamp, alteration1=alt1) %>% 
            rename('gene1'=gene)
          
          
          #calculate the co-occuring incidence of gene 2
          out <- tbl(conn, 'co_cnv') %>% 
            filter(alteration1 %in% g1_alt,gene1==igene_g1,cancer==icancer) %>% 
            group_by(gene1, gene2, cancer, alteration2) %>% 
            summarise(co_freq=sum(freq)) %>% 
            merge(g1_single, by=c('cancer','gene1')) %>% 
            mutate(co_rel_incidence=co_freq/nsamp) %>% 
            merge(seer, by='cancer') %>% 
            mutate(est_co_incidence=round(co_rel_incidence*total_incidence, 0)) %>% 
            select(-cancer_type) %>% 
            relocate(cancer, 
                     gene1, alteration1, 
                     gene2, alteration2, 
                     g1_freq, nsamp, g1_rel_incidence,
                     co_freq, co_rel_incidence)
          
          
        
          
        } else if (input$g2_tbl_mode=='SNP'){
          
          if (input$g1_tbl_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_tbl_alt=='Deletion') {g1_alt <<- c(-1,-2)}
          
          igene_g1 <- input$g1_tbl_igene
          icancer <- input$tbl_icancer
          alt1 <- input$g1_tbl_alt
          
          #calculate the incidence of the target gene alteration
          g1_single <- tbl(conn, 'single_cnv') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            group_by(cancer, gene, nsamp) %>% 
            summarise(g1_freq=sum(freq)) %>% 
            mutate(g1_rel_incidence=g1_freq/nsamp, alteration1=alt1) %>% 
            rename('gene1'=gene)
          
          g1_samps <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            pull(uniquePatientKey) %>% 
            unique()
          
          out <- tbl(conn, 'cnv_snp') %>% 
            filter(uniquePatientKey %in% g1_samps, cancer==icancer) %>% 
            select(cancer, gene, uniquePatientKey) %>% 
            distinct() %>%
            group_by(cancer, gene) %>% 
            summarise(co_freq = n()) %>% 
            merge(g1_single, by='cancer') %>%
            merge(seer, by='cancer') %>% 
            mutate(co_rel_incidence=co_freq/nsamp, alteration2='SNP') %>%
            mutate(est_co_incidence=round(co_rel_incidence*total_incidence, 0)) %>% 
            rename('gene2'=gene) %>% 
            select(-cancer_type) %>% 
            relocate(cancer, 
                     gene1, alteration1, 
                     gene2, alteration2, 
                     g1_freq, nsamp, g1_rel_incidence,
                     co_freq, co_rel_incidence)
          
        }
      
    }
    
    out <- out %>% 
      rename('cancer_incidence'=total_incidence)
    
    output$usr_tbl <- DT::renderDataTable({
      
      DT::datatable(data.frame(out), 
                    list(pageLength = 15, scrollX=T),rownames = F)
      
    })
    
    tbl_out <<- out
    
  })
})
