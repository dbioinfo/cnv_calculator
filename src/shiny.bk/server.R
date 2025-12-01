
# Define server logic
server <- shinyServer(function(input, output, session) {
  tdat <- reactiveVal()
  dels <-reactiveVal()
  dups <-reactiveVal()
  seer <- reactiveVal()
  igene_hist <- reactiveVal()
  meta_graph <- reactiveVal()
  cnv_samps <- reactiveVal()
  observeEvent(input$load_data,{
    #input data
    req(input$datfile)
    print(input$datfile$datapath)
    print('input tdat')
    tdat <<- readRDS(input$datfile$datapath)
    print('input dels')
    dels <<- tdat$del_co
    print('input dups')
    dups <<- tdat$dup_co
    tdat <<- tdat$single_cna %>% mutate(CNV=case_when(alteration=='del'~'Deletion', alteration=='dup'~'Duplication'))
    cnv_samps <<- tdat %>% 
      select(cancer_type, total_samp, hugoGeneSymbol, n, alteration) %>% 
      group_by(cancer_type, hugoGeneSymbol, alteration) %>% 
      summarize(n=sum(n), total_samp=unique(total_samp)) %>% 
      ungroup() %>% 
      group_by(cancer_type, alteration) %>% 
      summarize(n=length(unique(hugoGeneSymbol)), total_samp=sum(unique(total_samp))) %>%
      pivot_wider(names_from = 'alteration', values_from = 'n')
    
    #match with SEER db
    seer <<- read_csv('/data/seerstats.csv')
    tdat <<- tdat %>% merge(., seer, by='cancer_type')
    
    genes <- tdat %>% pull(hugoGeneSymbol) %>% unique()
    #populate gene choices
    output$gene_select <- renderUI({
      selectizeInput(inputId='igene',
                         'Choose Gene',
                         choices = genes,
                          selected='MTAP')
    })
    
    output$gene_co_list <- renderUI({
      selectizeInput(inputId='igene_co',
                     'Choose Gene',
                     choices = genes, #specifically rownames of dels/dups
                     selected='MTAP')
    })
    
    #plot data
    output$igene_hist <- renderPlot({
      req(input$igene)
      igene<- input$igene
      gg <- ggplot(tdat %>% filter(hugoGeneSymbol==igene))+
        geom_bar(mapping=aes(x=cancer_type, y=rel_incidence, fill=CNV), stat='identity', position = 'dodge')+
        geom_text(mapping=aes(x=cancer_type, y=rel_incidence, group=CNV, label=round(rel_incidence*total_incidence, 0) ), 
                  position = position_dodge(width=1),
                  vjust=-0.5)+
        scale_fill_manual(values = c('#f37a98','#7fb568'))+
        scale_y_continuous(expand=c(0,0), limits = c(0, 1.05*max(tdat %>% filter(hugoGeneSymbol==igene) %>% pull(rel_incidence))))+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
              text = element_text(family='Arial', size = 22),
              plot.title = element_text(hjust = 0.5))+
        xlab('')+
        ylab('Proportion of Samples')+
        labs(title = (paste0('Relative CNV Incidence - ', igene)))
      igene_hist <<- gg
      gg
        
    
    }, width=1000, height=800)
      
    output$download_igene_hist <- downloadHandler(
      filename = function() { paste(input$igene, '_cnv_incidence.png', sep='') },
      content = function(file) {
        ggsave(file, plot = igene_hist, device = "png", width = 4000, height=3000, units='px')
      }
    )
    
    output$download_meta_graph <- downloadHandler(
      filename = function() { paste(input$sel_meta, '.png', sep='') },
      content = function(file) {
        ggsave(file, plot = meta_graph, device = "png", width = 4000, height=3000, units='px')
      }
    )
    #output table of top 100 most cnvs
    output$top_cnv <- DT::renderDataTable({
      DT::datatable((tdat %>% group_by(hugoGeneSymbol, CNV) %>% summarize(Total_samples=sum(n))), 
                    list(pageLength = 10, scrollX=T),rownames = F)
    })
      
      
    })
  
  #metadata handler
  observeEvent(input$load_meta, {
  req(input$datfile)
  output$meta_graph <- renderPlot({

  if (input$sel_meta=='Cancer Incidence'){
    gg <- ggplot(seer)+
      geom_bar(mapping=aes(x=cancer_type,y=total_incidence,fill=cancer_type), stat='identity')+
      scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(nrow(seer)))+
      geom_text(mapping=aes(x=cancer_type,y=total_incidence, label=male), color='#3b91ed',vjust=-0.5, size=4)+
      geom_text(mapping=aes(x=cancer_type,y=total_incidence, label=female), color='#f542d7',vjust=-2, size=4)+
      theme_bw()+
      scale_y_continuous(trans='log2', breaks=c(10,100,1000,10000, 100000), labels = label_comma())+
      xlab('')+
      ylab('Incidence')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=22),
            axis.title = element_text(family='Arial', face='bold', size=22))+
      guides(fill='none')
  } else if (input$sel_meta=='Cancer CNV Samples') {
    gdat <- cnv_samps
    
    colourCount <- nrow(gdat)
    getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
    gdat <- gdat %>% filter(total_samp!=253)
    gg<-ggplot(gdat)+
      geom_bar(mapping=aes(x=cancer_type,y=total_samp,fill=cancer_type), stat='identity')+
      scale_fill_manual(values = getPalette(colourCount))+
      geom_text(mapping=aes(x=cancer_type,y=total_samp, label=del), color='#b5265f', size=4, vjust=-0.5)+
      geom_text(mapping=aes(x=cancer_type,y=total_samp, label=dup), color='#449e26',size=4, vjust =-2)+
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
    req(input$igene_co)
    
    output$geneco_graph <- renderPlot({
     if (input$del_or_dup =='Duplication'){
       subset <- data.frame(dups[,input$igene_co]) %>% 
         rownames_to_column('gene')
       colnames(subset) <- c('gene','freq')
      subset <- subset %>% arrange(-freq) %>% 
         head(n=30) %>% 
         filter(gene!=input$igene_co)
       gg<-ggplot(subset)+
         geom_bar(mapping=aes(x=reorder(gene, -freq), y=freq), stat='identity', fill='salmon')+
         theme_bw()+
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
               text = element_text(family='Arial', size = 18),
               plot.title = element_text(hjust = 0.5))+
         xlab('')+
         ylab('Number of Co=occurences')+
         ggtitle(paste0(input$igene_co, ' Co-Occurrence'))
     } else if (input$del_or_dup =='Deletion'){
       subset <- data.frame(dels[,input$igene_co]) %>% 
         rownames_to_column('gene')
       colnames(subset) <- c('gene','freq')
       subset <- subset %>% arrange(-freq) %>% 
         head(n=30) %>% 
         filter(gene!=input$igene_co)
       gg<-ggplot(subset)+
         geom_bar(mapping=aes(x=reorder(gene, -freq), y=freq), stat='identity', fill='salmon')+
         theme_bw()+
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
               text = element_text(family='Arial', size = 18),
               plot.title = element_text(hjust = 0.5))+
         xlab('')+
         ylab('Number of Co-occurences')+
         ggtitle(paste0(input$igene_co, ' Co-Occurrence'))
       
     }
    gg
    }, width = 1000, height=800)
    
    output$cotable <- DT::renderDataTable({
      if (input$del_or_dup =='Duplication'){
        subset <- data.frame(dups[,input$igene_co]) %>% 
          rownames_to_column('gene')
      } else if (input$del_or_dup =='Deletion'){
        subset <- data.frame(dels[,input$igene_co]) %>% 
          rownames_to_column('gene')
      }
      colnames(subset) <- c('gene','freq')
      subset <- subset %>% arrange(-freq) %>% 
        filter(gene!=input$igene_co)
      DT::datatable(subset, 
                    list(pageLength = 10, scrollX=T),rownames = F)
    })
  })
})
