library(shiny)
library(shinythemes)
library(RColorBrewer)
library(Matrix)
library(plotly)
library(scales)
library(DT)
library(readr)
library(stringr)
library(tidyr)
library(shinyjs)
library(sortable)
library(tibble)
library(rclipboard)
library(dplyr)
library(DBI)
library(cBioPortalData)
library(biomaRt)

options(shiny.maxRequestSize = 10000*1024^2)
tags <- htmltools::tags
select <- dplyr::select
rename <- dplyr::rename
ui <- shinyUI(navbarPage(title = "CNV Incidence Calculator v2.3", 
                         theme=shinytheme('flatly'),
                         tabsetPanel(
                           tabPanel("Single Gene",
                            sidebarLayout(
                              sidebarPanel(
                                actionButton('load_data','Load Database'), 
                                hr(),
                                uiOutput('gene_select'), 
                                hr(),
                                fluidRow(
                                  column(6,uiOutput('zygosity_single') ),
                                  column(6, uiOutput('alteration_single'))
                                  ),
                                hr(),
                                DT::dataTableOutput("top_cnv"),
                                hr(),
                                fluidRow(
                                  column(6,
                                         actionButton('load_single',label='Refresh Plot'),
                                         align='left'),
                                  column(6,
                                         downloadButton('download_single_tbl',label='Download Table'),
                                         align='right')
                                )
                              ),
                              mainPanel(hr(), fluidRow(align='center', 
                                                       plotOutput("igene_hist")
                                                       ),
                                              fluidRow(align='center', 
                                                      br(), br(), br(), br(),
                                                      br(), br(), br(), br(),
                                                      br(), br(), br(), br(),
                                                      br(), br()
                                              ), 
                                              fluidRow(align='right',
                                                       downloadButton('download_igene_hist','Download',align='bottom')
                                                       )
                                        )
                              )
                           ),
                           tabPanel("Co-occurrence",
                                    sidebarLayout(
                                      sidebarPanel(
                                      
                                        
                                        #choice of which graph to generate dynamically changes input options
                                        h2("Select Graphing Options"),
                                        fluidRow(
                                          column(6,
                                                radioButtons('co_gopt', label='',
                                                     choices= c('Known Gene Pair'='gpair',
                                                                'Single Gene Explorer'='gsingle'),
                                                     selected = 'Known Gene Pair')
                                                ), 
                                          conditionalPanel(condition ="input.co_gopt == 'gpair'",
                                              column(6,
                                                 radioButtons('zygosity_co', choices = c('True','False'), label = 'Show Zygosity?', selected = 'True'),
                                                 hr(),
                                                 )
                                          )
                                        ),
                                        
                                        
                                        #dynamic ui
                                        conditionalPanel(condition = "input.co_gopt == 'gpair'",
                                                         fluidRow(
                                                           column(3, radioButtons('g1_alt','Gene 1 Alteration',choices=c('Deletion','Duplication','SNP'),selected = 'Deletion')),
                                                           column(3, radioButtons('g2_tbl','Gene 2 Alteration',choices=c('Deletion','Duplication','SNP'),selected = 'Deletion')),
                                                          column(3,uiOutput('g1_co_list')),
                                                          column(3,uiOutput('g2_co_list'))
                                                         )
                                        ),
                                        conditionalPanel(condition = "input.co_gopt == 'gsingle'",
                                                         fluidRow(
                                                         column(3,radioButtons('g1_alt','Gene 1 Alteration',choices=c('Deletion','Duplication','SNP'),selected = 'Deletion')),
                                                         column(3,radioButtons('g2_alt','Gene 2 Alteration',choices=c('Deletion','Duplication','SNP'),selected = 'Deletion')),
                                                         column(6, uiOutput('gene_co_list'))
                                                         )
                                                         
                                        ),
                                        
                                        
                                        
                                        br(),
                                        fluidRow(
                                          column(12, actionButton('load_cograph','Load Co-occurrence Plot'), align='right')
                                        ),
                                        br(),
                                        DT::dataTableOutput('cotable'), #data table,
                                        fluidRow(
                                          column(12,
                                                 downloadButton('download_cotable',label='Download Table'),
                                                 align='right')
                                        )
                                      ),
                                      mainPanel(
                                        fluidRow(align='center',
                                                 conditionalPanel(condition = "input.co_gopt == 'gpair'",
                                                                  plotOutput('geneco_graph_pair')
                                                                  
                                                 ),
                                                 conditionalPanel(condition = "input.co_gopt == 'gsingle'",
                                                                   plotOutput('geneco_graph_single'),
                                                                  
                                                 )
                                                 
                                        ), 
                                        fluidRow(align='right', 
                                                 br(), br(), br(), br(),br(), br(), br(), br(),br(), br(), br(), br(),br(), br(),
                                                 conditionalPanel(condition = "input.co_gopt == 'gpair'",
                                                                  column(12,
                                                                         downloadButton('download_co_pair',label='Download'),
                                                                         align='right')
                                                 ),
                                                 conditionalPanel(condition = "input.co_gopt == 'gsingle'",
                                                                  column(12,
                                                                         downloadButton('download_co_single',label='Download'),
                                                                         align='right')
                                                 )
                                        )
                                      )
                                    )   
                           ),
                           tabPanel("Table Editor",
                                    sidebarLayout(
                                      sidebarPanel(
                                        h2("Select a Gene and Filtering Option"),
                                        br(),br(),
                                        
                                        fluidRow(
                                          column(6,uiOutput('g1_tbl_list')),
                                          column(6,selectInput('g1_tbl_alt','Target Gene Alteration', choices = c('Deletion','Duplication'),selected = 'Deletion'))
                                          ),
                                        br(),
                                        fluidRow( 
                                                 column(12,
                                                        radioButtons('tbl_mode',label = 'Query Mode', choices = c('Gene-Specific','Cancer-Specific'), selected='Gene-Specific'),
                                                        align='center')
                                                 ),
                                        br(),
                                        conditionalPanel(condition = "input.tbl_mode == 'Gene-Specific'",
                                                         fluidRow(
                                                           column(6,uiOutput('g2_tbl_list')),
                                                           column(6,selectInput('g2_tbl_alt','Gene 2 Alteration', choices = c('Deletion','Duplication','SNP'),selected = 'Deletion'))
                                                         ),         
                                        ), 
                                        
                                        conditionalPanel(condition = "input.tbl_mode == 'Cancer-Specific'",
                                                        fluidRow(
                                                          column(6, radioButtons('g2_tbl_mode', label = 'CNV or SNP', choices = c('CNV','SNP'), selected = 'CNV')),
                                                          column(6, uiOutput('cancer_list_tbl'))
                                                        )                
                                        ),
                                        br(),
                                        fluidRow(
                                          column(12, actionButton('load_tbl','Load Table'), align='right')
                                          )
                                      ),
                                      mainPanel(
                                        br(), br(),
                                        DT::dataTableOutput('usr_tbl'),
                                        br(), 
                                        fluidRow(
                                                column(12,
                                                       downloadButton('download_tbl',"Download Table"), 
                                                       align='right'
                                                )
                                        )
                                      )
                                    )
                           ),
                           tabPanel("Metadata",
                            sidebarLayout(
                              sidebarPanel(
                                selectInput('sel_meta', 'Select Data', choices = c('Cancer Incidence','Cancer CNV Samples')),
                                actionButton('load_meta','Load Metadata Graph')),
                              mainPanel(
                                fluidRow(align='center',
                                         plotOutput('meta_graph')
                                ), 
                                fluidRow(align='center', 
                                         br(), br(), br(), br(),
                                         br(), br(), br(), br(),
                                         br(), br(), br(), br(),
                                         br(), br()
                                ), 
                                fluidRow(align='right',
                                         downloadButton('download_meta_graph','Download Figure',align='bottom')
                                )
                              )
                            )   
                           ),
			  tabPanel("Update Database",
			    fluidPage(
			      #this custom CSS bit allows us to make the area of the list a bit longer (or, it should)
			      tags$style(HTML("
				    .selectize-dropdown-content {
				      max-height: 500px !important; /* Adjust dropdown height */
				      overflow-y: auto !important; /* Add scroll for long lists */
				    }
				    .selectize-input {
				      min-height: 30px; /* Adjust the height of the selection box */
				    }
				  ")),
			      fluidRow(align='center',br(),br()),
			      fluidRow(align='center',
				       column(6, uiOutput('extant_gene_select')),
				       column(6, uiOutput('new_gene_select'))
			      ),
			      fluidRow(align='center',
				       column(6, textOutput("ngenes") ),
				       column(6, textOutput("noffgenes"))
			      ),
			      fluidRow(align='center',
				       column(12, actionButton('update_db',label='Update Database'))
			      ),
			      fluidRow(align='center',
				       column(12, textOutput("update_progress") )
			      )
			   )
                         )
		)
))




# Define server logic
server <- shinyServer(function(input, output, session) {
  #log_connection <- file("/logs/temp.log", open = "a") 
  #sink(log_connection, type='message')
  #sink(log_connection, type-'output')
  #on.exit({sink(type='message');sink(type='output')})
  conn <- reactiveVal()
  tdat <- reactiveVal()
  tdat_single <- reactiveVal()
  dels <-reactiveVal()
  dups <-reactiveVal()
  genes <- reactiveVal()
  cancers <- reactiveVal()
  subset <- reactiveVal()
  seer <- reactiveVal()
  igene_hist <- reactiveVal()
  meta_graph <- reactiveVal()
  g1_alt <- reactiveVal()
  g2_alt <- reactiveVal()
  snpsub <- reactiveVal()
  tbl_out <- reactiveVal()
  co_tbl <- reactiveVal()
  single_tbl <- reactiveVal()
  nsamps <- reactiveVal()
  ngenes <- reactiveVal()
  noffgenes <- reactiveVal()
  fill_colors <- reactiveVal()
  gdat <- reactiveVal()
  #select function
  print('update selecter')
  select <- dplyr::select

  #params
  cnv_colors <<- c('#f37a98','#7fb568')
  alteration_colors <<- c('#f26161','#f29696','#8abf88','#41b33b')
  
observeEvent(input$load_data,{
    #input data
    print('connecting to server')
    conn <<- dbConnect(drv = RPostgres::Postgres(),
                      dbname='public_cnv',  #public_cnv
                      host='cnv_db',  #cnv_db
                      user='postgres',
                      password='password', 
                      port=5432)
    dbExecute(conn, "SET work_mem = '1.5GB';") #utilize more ram
    tdat <<- tbl(conn, 'cnv_main') %>% mutate(CNV=case_when(alteration<0~'Deletion', alteration>0~'Duplication'))
    
   
    
    #seer is now a table within our db
    seer <<- tbl(conn, 'seer_main')
    nsamps <<- tbl(conn, 'cnv_nsamps') %>% select(cancer, nsamp) %>% distinct()
    genes <<- tbl(conn, 'genes') %>% pull(genes)
    offsite_genes <<- tbl(conn, 'offsite_genes') %>% pull(newgenes)
    ngenes <<- length(genes)
    noffgenes<<-length(offsite_genes)
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
     
    output$new_gene_select <- renderUI({
      selectizeInput(inputId='offsite_genes',
                         'Choose Offsite Genes to Download',
                         choices = offsite_genes,
                         multiple=T)
    })

    output$extant_gene_select <- renderUI({
      selectizeInput(inputId='extant_genes',
                         'Browse Extant Genes',
                         choices = genes,
                         multiple=T)
    })
    output$ngenes <- renderText({ paste0("Number of extant genes: ", ngenes )})
    output$noffgenes <- renderText({ paste0("Number of offsite genes: ",  noffgenes )})
   
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
    
    

    
    })

#####
# Single Graphs
#####
observeEvent(input$load_single, {
  print('rendering single cnv hist...')
    req(input$igene)
    igene <- input$igene
    
    print(igene)
    
    
    
    if (input$alteration_single=='CNV'){
      
      #once igene is chosen, it makes the database much easier to nav
      gdat <<- tdat %>%   group_by(cancer, gene, alteration) %>% summarise(freq=n()) %>% 
        mutate(nsamp=sum(freq), rel_incidence=freq/nsamp)
      gdat <<- gdat %>%   filter(gene==igene) %>% group_by(cancer, alteration) %>% 
        left_join(., seer, by='cancer') %>% 
        summarise(rel_incidence=sum(rel_incidence), total_incidence=max(total_incidence)) %>% 
        mutate(text=round(rel_incidence*total_incidence, 0)) 
      gdat <<- gdat %>% filter(alteration!=0) %>% 
        group_by(cancer, alteration) %>% 
        mutate(text=round(rel_incidence*total_incidence, 0)) 
      
      
      if (input$zygosity_single=='False'){ #meaning we do NOT show zygosity
        gdat <<- gdat %>% 
          mutate(CNV=case_when(
            alteration<0~'Deletion',
            alteration>0~'Duplication'
          )) %>%
          ungroup() %>% 
          select(-alteration) %>%
          group_by(cancer, CNV) %>% 
          summarize(rel_incidence=sum(rel_incidence), total_incidence=max(total_incidence)) %>%
          mutate(text=round(rel_incidence*total_incidence, 0))
        
        ### ~graph~ ###
        gg <- ggplot(gdat)+
          geom_bar(mapping=aes(x=cancer, y=rel_incidence, fill=CNV), stat='identity', position = 'dodge')+
          geom_text(mapping=aes(x=cancer, y=rel_incidence, group=CNV, label=text), 
                    position = position_dodge(width=1),
                    vjust=-0.5)+
          scale_fill_manual(values = cnv_colors)+
          scale_y_continuous(expand=c(0,0), limits = c(0, 1))+
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                text = element_text(family='Arial', size = 22),
                plot.title = element_text(hjust = 0.5))+
          xlab('')+
          ylab('Proportion of Samples')+
          labs(title = (paste0('Relative CNV Incidence - ', igene)))
      } else if (input$zygosity_single=='True'){ #this is the default, meaning we DO show zygosity
        
        ### ~graph~ ###
        gg <- ggplot(gdat)+
          geom_bar(mapping=aes(x=cancer, y=rel_incidence, fill=factor(alteration)), stat='identity', position = 'dodge')+
          geom_text(mapping=aes(x=cancer, y=rel_incidence, group=factor(alteration), label=text), 
                    position = position_dodge(width=1),
                    vjust=-0.5)+
          scale_fill_manual(values = alteration_colors)+
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
    } else if(input$alteration_single=='SNP') { #access different table if snp data is needed
      
      snpdat <- tbl(conn, 'cnv_main') %>% #flatten multiple snps in a single patient
        group_by(cancer, gene, uniquePatientKey) %>% 
        filter(gene==igene) %>% 
        distinct()
      
      snpdat <- snpdat %>% #get freq of snps
        group_by(cancer, gene) %>% 
        summarize(freq=n())
      
      nsamps <- tbl(conn, 'cnv_nsamps') %>% #get total samples per cancer
        select(cancer, nsamp) %>%
        group_by(cancer) %>% 
        summarize(nsamp=max(nsamp))
      
      gdat <<- left_join(snpdat, nsamps, by='cancer') %>% #paste data all together 
        mutate(rel_incidence = as.numeric(freq )/ nsamp) %>%
        left_join(., seer, by='cancer') %>% 
        mutate(text=round(rel_incidence*total_incidence, 0))
      
      
      ### ~graph~ ###
      gg <- ggplot(gdat)+
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
    output$igene_hist <- renderPlot({
      igene_hist
       
     }, width=1000, height=800)
  
  #output table of top 100 most cnvs
  output$top_cnv <- DT::renderDataTable({
    igene <- input$igene
    
    if (input$zygosity_single=='False'){
      subset <- tdat %>% 
        filter(gene==igene) %>%
	as.data.frame() %>%
        mutate(CNV=case_when(alteration<0~'Deletion', alteration>0~'Duplication', alteration==0~'No Change')) %>% 
        select(-alteration) %>%
        group_by(cancer, CNV) %>% 
        summarise(freq=n()) %>% 
        mutate(nsamp=sum(freq), rel_incidence=freq/nsamp) %>% 
        filter(CNV!='No Change') %>% 
        ungroup() %>% 
        relocate(cancer, CNV) 
        
    } else {
      subset <- tdat %>% 
        filter(gene==igene) %>%
        group_by(cancer, alteration) %>% 
        summarise(freq=n()) %>% 
        mutate(nsamp=sum(freq), rel_incidence=freq/nsamp) %>% 
        filter(alteration!=0) #%>%  
        #
    }
    
    
    
    DT::datatable(as.data.frame(subset) %>% arrange(-rel_incidence),  list(pageLength = 10, scrollX=T), rownames = F)
    
  })
  
})


#####
# Metadata graphs
#####
observeEvent(input$load_meta, {
  
  print('rendering meta plot...')
  output$meta_graph <- renderPlot({

  if (input$sel_meta=='Cancer Incidence'){
    gg <- ggplot(seer)+
      geom_bar(mapping=aes(x=cancer,y=total_incidence,fill=cancer), stat='identity')+
      scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(22))+
      geom_text(mapping=aes(x=cancer, y=total_incidence, label=male), color='#3b91ed',vjust=-0.5, size=4)+
      geom_text(mapping=aes(x=cancer, y=total_incidence, label=female), color='#f542d7',vjust=-2, size=4)+
      theme_bw()+
      scale_y_continuous(trans='log2', breaks=c(10,100,1000,10000, 100000), labels = label_comma())+
      xlab('')+
      ylab('Incidence')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=22),
            axis.title = element_text(family='Arial', face='bold', size=22))+
      guides(fill='none')
  
      
  } else if (input$sel_meta=='Cancer CNV Samples') {
    gdat <- tdat %>% select(gene, alteration, cancer) %>% group_by(cancer, gene, alteration) %>% summarise(freq=n())
    gdat <- gdat %>% 
      ungroup() %>%
      select(cancer, alteration, freq) %>%
      filter(alteration!=0) %>% 
      mutate(alteration=case_when(
        alteration<0 ~'Deletion',
        alteration>0 ~'Duplication',
        .default='NA'
      )) %>% 
      group_by(cancer, alteration) %>% 
      summarise(freq=sum(freq))
    
    
    
    gg<-ggplot(gdat)+
      geom_bar(mapping=aes(x=cancer,y=freq,fill=alteration), stat='identity', position='dodge')+
      scale_fill_manual(values = cnv_colors)+
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
  
#####
# Co-occurence graphs and tables
#####
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
          filter(gene==igene_g1) %>% #filter(alteration=='SNP') %>%
          select(cancer, nsamp) %>% 
          group_by(cancer) %>% 
          summarize(nsamp=max(nsamp)) %>% 
          as.data.frame()
        
        if (input$g2_tbl!='SNP'){ #input for g2 alt is either del or dup
          
	    if (input$g2_alt=='Duplication'){g2_alt <<- c(1,2)}else if (input$g2_alt=='Deletion') {g2_alt <<- c(-1,-2)}
            subset <<- tbl(conn, 'cnv_main') %>%
              filter(uniquePatientKey %in% !! g1_patients$uniquePatientKey, gene==igene_g2, alteration %in% g2_alt, cancer %in% cancers) %>% 
              select(cancer, gene, alteration) %>% 
              group_by(cancer, gene, alteration) %>% 
              summarize(freq=n()) %>% 
              as.data.frame() %>% 
              merge(., nsamps, by='cancer') %>% 
              mutate(rel_incidence = freq / nsamp) %>% 
	      arrange(-rel_incidence)
          
            
            if (input$zygosity_co=='False'){ #graph without zygosity

              subset <<- subset %>% mutate(alteration=case_when(
                                          alteration<0~'Deletion',
                                          alteration>0~'Duplication'))

              fill_colors<-cnv_colors
            
            
            } else if (input$zygosity_co=='True'){ #graph with zygosity
            
              fill_colors<-alteration_colors
            
            }
	    
          #plot render 
          output$geneco_graph_pair<- renderPlot({
            gg<-ggplot(subset)+
              geom_bar(mapping=aes(x=cancer, y=rel_incidence, fill=factor(alteration) ), stat='identity', position = 'dodge')+
              scale_y_continuous(expand = c(0,0),limits=c(0,1))+
              scale_fill_manual(values=fill_colors)+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                    text = element_text(family='Arial', size = 18),
                    plot.title = element_text(hjust = 0.5))+
              xlab('')+
              ylab(paste0('Proportion of ',igene_g1, ' SNP ',
                          '& ',igene_g2, ' ',input$g2_alt,' Co-occurences') )+
              labs(title=paste0(igene_g1 , ' Mutation Co-Occurrence'),
                   subtitle = paste0(
                     paste0(igene_g1,' SNP  & '),
                     paste0(igene_g2, ' ', input$g2_alt)
                   ))+ 
              guides(fill=guide_legend(title=paste0(igene_g2," CNV") ))  
            
          
            gg  
          }, width = 1000, height=800)

	  #table render
	  
	  output$cotable <- DT::renderDataTable({

	    co_tbl <<- data.frame(subset)
	    DT::datatable(co_tbl, 
		      list(pageLength = 10, scrollX=T),rownames = F)

	  })

        } else if (input$g2_tbl=='SNP'){ #input for g2 alt is also SNP
          
          subset <<- tbl(conn, 'cnv_snp') %>% 
              filter(uniquePatientKey %in% !! g1_patients$uniquePatientKey, cancer %in% cancers, gene==igene_g2) %>% 
              select(cancer, gene) %>%
              group_by(cancer, gene) %>% 
              summarize(freq = n()) %>% 
              as.data.frame() %>% 
              merge(., nsamps, by='cancer') %>% 
              mutate(rel_incidence = freq / nsamp) %>%
	      arrange(-rel_incidence)
          
          
            
          
	  #plot render
          output$geneco_graph_pair<- renderPlot({
            gg <- ggplot(subset) + 
              geom_bar(mapping=aes(x=cancer, y=rel_incidence, fill=cancer), stat='identity', position='dodge')+
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

	  #table render

	  output$cotable <- DT::renderDataTable({

	    co_tbl <<- data.frame(subset)
	    DT::datatable(co_tbl, 
		      list(pageLength = 10, scrollX=T),rownames = F)

	  })

        }
      } else { #input for gene 1 alt is CNV
     
        if (input$g2_tbl!='SNP'){
      
          req(input$g1_igene_co)
          req(input$g2_igene_co)
          print('geneco_graph_pair')
          
          if (input$g1_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_alt=='Deletion') {g1_alt <<- c(-1,-2)}
	  if (input$g2_alt=='Duplication'){g2_alt <<- c(1,2)}else if (input$g2_alt=='Deletion') {g2_alt <<- c(-1,-2)}
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
            subset <<- tbl(conn, 'co_cnv') %>% 
              filter(gene1==igene_g1, gene2==igene_g2, alteration1 %in% g1_alt, alteration2 %in% g2_alt) %>% 
              mutate(alteration2=case_when(alteration2<0~'Deletion', alteration2>0~'Duplication')) %>% 
              group_by(cancer, gene1, gene2, alteration2) %>% 
              summarise(freq = sum(freq)) %>% 
              as_data_frame() %>% 
              merge(igene_freq, by='cancer') %>% 
              mutate(rel_incidence = freq/tfreq)
        
            
            fill_colors<<-cnv_colors
           

          } else if (input$zygosity_co=='True'){
            #get igene cancer freq
            igene_freq <- tbl(conn, 'single_cnv') %>% 
              filter(gene==igene_g1, alteration %in% g1_alt) %>% 
              select(cancer, freq, alteration) %>% 
              group_by(cancer) %>% 
              summarise(tfreq=sum(freq)) %>% 
              as_data_frame()
            
            #filter for subset of samples
            subset <<- tbl(conn, 'co_cnv') %>% 
              filter(gene1==igene_g1, gene2==igene_g2, alteration1 %in% g1_alt, alteration2 %in% g2_alt) %>% 
              group_by(cancer, gene1, gene2, alteration2) %>% 
              summarise(freq = sum(freq)) %>% 
              as_data_frame() %>% 
              merge(igene_freq, by='cancer') %>% 
              mutate(rel_incidence = freq/tfreq, alteration2=factor(alteration2)) %>%
	      arrange(-rel_incidence)
            
            
            fill_colors<<-alteration_colors



          }
          #render plot 
          output$geneco_graph_pair<- renderPlot({
	  if (input$zygosity_co=='True'){
          gg<-ggplot(subset)+
              geom_bar(mapping=aes(x=cancer, y=rel_incidence, fill=alteration2 ), stat='identity', position = 'dodge')+
              scale_y_continuous(expand = c(0,0), limits=c(0,1))+
              scale_fill_manual(values=fill_colors)+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                    text = element_text(family='Arial', size = 18),
                    plot.title = element_text(hjust = 0.5))+
              xlab('')+
              ylab(paste0('Proportion of ',igene_g1, ' ', substr(input$g1_alt, 1, 3),
                          ' & ',igene_g2, ' ', substr(input$g2_alt, 1,3)  ,' Co-occurences') )+
              labs(title=paste0(igene_g1, ' Co-Occurrence'),
                   subtitle = paste0(
                     paste0(igene_g1,' ',input$g1_alt,' & '),
                     paste0(igene_g2, ' ', input$g2_alt)
                   ))+ 
              guides(fill=guide_legend(title=paste0(igene_g2, ' ', input$g2_alt)))
          gg
	  } else {

          gg<-ggplot(subset)+
              geom_bar(mapping=aes(x=cancer, y=rel_incidence, fill=cancer), stat='identity', position = 'dodge')+
              scale_y_continuous(expand = c(0,0), limits=c(0,1))+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                    text = element_text(family='Arial', size = 18),
                    plot.title = element_text(hjust = 0.5))+
              xlab('')+
              ylab(paste0('Proportion of ',igene_g1, ' ', substr(input$g1_alt, 1, 3),
                          ' & ',igene_g2, ' ', substr(input$g2_alt, 1,3)  ,' Co-occurences') )+
              labs(title=paste0(igene_g1, ' Co-Occurrence'),
                   subtitle = paste0(
                     paste0(igene_g1,' ',input$g1_alt,' & '),
                     paste0(igene_g2, ' ', input$g2_alt)
                   ))+ 
              guides(fill='none')
          gg

	  }
        }, width = 1000, height=800)

        ##render table 
	output$cotable <- DT::renderDataTable({

	  co_tbl <<- data.frame(subset)
	  DT::datatable(co_tbl ,#%>% relocate(gene1, alteration1, gene2, alteration2, rel_incidence), 
		      list(pageLength = 10, scrollX=T),rownames = F)

	})
      
      
      } else if (input$g2_tbl=='SNP'){
          
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
          
          
          subset <<- tbl(conn, 'cnv_snp') %>% 
            filter(uniquePatientKey %in% igene_samps, gene == igene_g2) %>% 
            select(uniquePatientKey, cancer, gene) %>% 
            distinct() %>% #gets rid of double counted snp (ie gene has 2+ mutations)
            group_by(cancer, gene) %>% 
            summarise(snpfreq=n()) %>% 
            merge(igene_freq, by='cancer') %>% 
            mutate(rel_incidence=snpfreq/tfreq) %>% 
            arrange(-rel_incidence) 
          
          
      #render plot
        output$geneco_graph_pair<- renderPlot({
          gg<-ggplot(subset)+
            geom_bar(mapping=aes(x=cancer, y=rel_incidence, fill=cancer), stat='identity', position = 'dodge')+
            scale_y_continuous(expand = c(0,0))+
            theme_bw()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'bold'),
                  text = element_text(family='Arial', size = 18),
                  plot.title = element_text(hjust = 0.5))+
            xlab('')+
            ylab(paste0('Proportion of ',igene_g1, ' ', substr(input$g1_alt, 1, 3),
                        ' & ',igene_g2,' SNP Co-occurences') )+
            labs(title=paste0(igene_g1, ' SNP Co-Occurrence'),
                 subtitle = paste0(
                   paste0(igene_g1,' ',input$g1_alt,' & '),
                   paste0(igene_g2)
                 ))
          
          gg
        }, width = 1000, height=800)
        #render table 
        output$cotable <<- DT::renderDataTable({
          req(input$g1_igene_co)
          req(input$g2_igene_co)
          igene_g1 <- input$g1_igene_co
          igene_g2 <- input$g2_igene_co
          co_tbl <<- data.frame(subset)
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
            
            subset <<- tbl(conn, 'cnv_main') %>%
              filter(uniquePatientKey %in% !! g_patients$uniquePatientKey, alteration %in% g1_alt, cancer %in% cancers) %>% 
              select(cancer, gene, alteration) %>% 
              group_by(cancer, gene, alteration) %>% 
              summarize(freq=n()) %>% 
              as.data.frame() %>% 
              merge(., nsamps, by='cancer') %>% 
              mutate(rel_incidence = freq / nsamp) 

            topngenes <- subset %>% arrange(desc(rel_incidence)) %>% 
              head(n=100) %>% select(gene) %>% distinct() %>% head(n=25) %>% as.data.frame()
            
            subset <<- subset %>% 
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
            
            subset <<- tbl(conn, 'cnv_snp') %>%
              filter(uniquePatientKey %in% !! g_patients$uniquePatientKey, cancer %in% cancers) %>% 
              select(cancer, gene) %>% 
              group_by(cancer, gene) %>% 
              summarize(freq=n()) %>% 
              as.data.frame() %>% 
              merge(., nsamps, by='cancer') %>% 
              mutate(rel_incidence = freq / nsamp)
            print(cancers)
            topngenes <- subset %>% arrange(desc(rel_incidence)) %>% 
              head(n=100) %>% select(gene) %>% distinct() %>% head(n=25) %>% as.data.frame()
            
            subset <<- subset %>% 
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
              subset <<- tbl(conn, 'co_cnv') %>% 
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
            
            #output$cotable <<- DT::renderDataTable({
            #  igene <- input$igene
            #  subset <- tbl(conn, 'co_cnv') %>% 
            #    filter(alteration1 %in% g1_alt,alteration2 %in% g2_alt,gene1==igene,gene2!=igene) %>% 
            #    arrange(desc(freq)) %>% head(n=50)
              
            #  co_tbl <<- data.frame(subset)
            #  DT::datatable(co_tbl %>% relocate(gene1, alteration1, gene2, alteration2), 
            #                list(pageLength = 10, scrollX=T),rownames = F)
              
              
            #})
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
              
              subset <<- tbl(conn, 'cnv_snp') %>% 
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
  


#####
# Table editor
#####

observeEvent(input$load_tbl,{
    
    if (input$tbl_mode == 'Gene-Specific'){
        if (input$g2_tbl_alt!='SNP'){
          
          if (input$g1_tbl_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_tbl_alt=='Deletion') {g1_alt <<- c(-1,-2)}
          if (input$g2_tbl_alt=='Duplication'){g2_alt <<- c(1,2)}else if (input$g2_tbl_alt=='Deletion') {g2_alt <<- c(-1,-2)}
          
          igene_g1 <- input$g1_tbl_igene
          igene_g2 <- input$g2_tbl_igene
          alt1 <- input$g1_tbl_alt
          alt2 <- input$g2_tbl_alt
          
	  nsamps <- tbl(conn, 'cnv_nsamps') %>% filter(gene==igene_g1) %>% select(cancer, nsamp) %>% summarize(nsamp=sum(nsamp))
          #calculate the incidence of the target gene alteration
          g1_single <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            group_by(cancer, gene) %>% 
            summarise(g1_freq=n()) %>% 
            left_join(., nsamps, by='cancer') %>% 
            mutate(g1_rel_incidence=as.double(g1_freq)/nsamp, alteration1=alt1) %>% 
            rename('gene1'=gene)
          
          
          #calculate the co-occuring incidence of gene 2 (needs rewriting)
          g1_patients <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            select(uniquePatientKey) %>% 
            distinct() %>% 
            pull(uniquePatientKey)
          
          g2_single <- tbl(conn, 'cnv_main') %>% 
            filter(uniquePatientKey %in% g1_patients, gene==igene_g2, alteration %in% g2_alt) %>%
            group_by(cancer, gene) %>%
            summarise(co_freq=n()) %>%
            left_join(., nsamps, by='cancer') %>%
            mutate(co_rel_incidence=as.double(co_freq)/nsamp, alteration2=alt2) %>%
            rename('gene2'=gene)
            
          out <- left_join(g1_single, g2_single, by=c('cancer','nsamp')) %>% 
            relocate(cancer, 
                     gene1, alteration1, 
                     gene2, alteration2, 
                     nsamp,
                     g1_freq,  g1_rel_incidence,
                     co_freq, co_rel_incidence)
            
            
          
        } else if (input$g2_tbl_alt=='SNP'){
          
          if (input$g1_tbl_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_tbl_alt=='Deletion') {g1_alt <<- c(-1,-2)}
          
          igene_g1 <- input$g1_tbl_igene
          igene_g2 <- input$g2_tbl_igene
          alt1 <- input$g1_tbl_alt
          
	  nsamps <- tbl(conn, 'cnv_nsamps') %>% filter(gene==igene_g1) %>% group_by(cancer, nsamp) %>% summarize(nsamp=sum(nsamp))
          #calculate the incidence of the target gene alteration
          g1_single <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            group_by(cancer, gene) %>% 
            summarise(g1_freq=n()) %>% 
            left_join(., nsamps, by='cancer') %>% 
            mutate(g1_rel_incidence=as.double(g1_freq)/nsamp, alteration1=alt1) %>% 
            rename('gene1'=gene)
          
          g1_patients <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt) %>% 
            select(uniquePatientKey) %>% 
            distinct() %>% 
            pull(uniquePatientKey)
          
          g2_single <- tbl(conn, 'cnv_snp') %>% 
            filter(uniquePatientKey %in% g1_patients, gene==igene_g2) %>% 
	    mutate(location = paste0('chr',chr, ':', startPosition)) %>%
            select(cancer, gene, uniquePatientKey, referenceAllele, variantAllele, location) %>% 
            distinct() %>%
            group_by(cancer, gene, referenceAllele, variantAllele, location) %>% 
            summarise(co_freq = n()) %>% 
            left_join(., nsamps, by='cancer') %>% 
            mutate(co_rel_incidence=as.double(co_freq)/nsamp, alteration2='SNP') %>% 
            rename('gene2'=gene)
          
          out <- left_join(g1_single, g2_single, by=c('cancer', 'nsamp')) %>% 
            relocate(cancer, 
                     gene1, alteration1, 
                     gene2, referenceAllele, variantAllele, location, 
                     nsamp,
                     g1_freq,  g1_rel_incidence,
                     co_freq, co_rel_incidence) %>%
	  filter(gene2!="")
          

        }
       
    } else if (input$tbl_mode == 'Cancer-Specific'){
        if (input$g2_tbl_mode=='CNV'){
          
          if (input$g1_tbl_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_tbl_alt=='Deletion') {g1_alt <<- c(-1,-2)}
          
          igene_g1 <- input$g1_tbl_igene
          icancer <- input$tbl_icancer
          alt1 <- input$g1_tbl_alt
          
          #calculate the incidence of the target gene alteration
          g1_single <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt, cancer==icancer) %>% 
            group_by(cancer, gene) %>% 
            summarise(g1_freq=n()) %>% 
            left_join(., nsamps, by='cancer') %>% 
            mutate(g1_rel_incidence=as.double(g1_freq)/nsamp, alteration1=alt1) %>% 
            rename('gene1'=gene)
          
          
          #calculate the co-occuring incidence within cancer
          g1_patients <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt, cancer==icancer) %>% 
            select(uniquePatientKey) %>% 
            distinct() %>% 
            pull(uniquePatientKey)
          
          cancer_cnv <- tbl(conn, 'cnv_main') %>% 
            filter(cancer==icancer, uniquePatientKey %in% g1_patients, alteration!=0) %>%
            mutate(alteration2=
                     case_when(alteration<0 ~ 'Deletion',
                               alteration>0 ~ 'Duplication')) %>%
            rename(gene2='gene') %>% 
            group_by(cancer, gene2, alteration2) %>%
            summarise(co_freq=n()) 
          
          out <- left_join(g1_single, cancer_cnv, by='cancer') %>% 
            mutate(co_rel_incidence = as.double(co_freq)/nsamp) %>% 
            relocate(cancer, 
                     gene1, alteration1, 
                     gene2, alteration2, 
                     nsamp,
                     g1_freq,  g1_rel_incidence,
                     co_freq, co_rel_incidence)
          

        } else if (input$g2_tbl_mode=='SNP'){
          
          if (input$g1_tbl_alt=='Duplication'){g1_alt <<- c(1,2)}else if (input$g1_tbl_alt=='Deletion') {g1_alt <<- c(-1,-2)}
          
          igene_g1 <- input$g1_tbl_igene
          icancer <- input$tbl_icancer
          alt1 <- input$g1_tbl_alt
          
          #calculate the incidence of the target gene alteration
          g1_single <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt, cancer==icancer) %>% 
            group_by(cancer, gene) %>% 
            summarise(g1_freq=n()) %>% 
            left_join(., nsamps, by='cancer') %>% 
            mutate(g1_rel_incidence=as.double(g1_freq)/nsamp, alteration1=alt1) %>% 
            rename('gene1'=gene)
          
          
          #calculate the co-occuring incidence within cancer
          g1_patients <- tbl(conn, 'cnv_main') %>% 
            filter(gene==igene_g1, alteration %in% g1_alt, cancer==icancer) %>% 
            select(uniquePatientKey) %>% 
            distinct() %>% 
            pull(uniquePatientKey)
          
          cancer_snp <- tbl(conn, 'cnv_snp') %>% 
            filter(cancer==icancer, uniquePatientKey %in% g1_patients) %>%
            rename(gene2='gene') %>% 
            select(cancer, gene2, uniquePatientKey) %>% 
            distinct() %>%
            group_by(cancer, gene2) %>% 
            summarise(co_freq = n()) %>% 
            left_join(., nsamps, by='cancer') %>% 
            mutate(co_rel_incidence=as.double(co_freq)/nsamp, alteration2='SNP') 
          
          out <- left_join(g1_single, cancer_snp, by=c('cancer','nsamp') ) %>% 
            relocate(cancer, 
                     gene1, alteration1, 
                     gene2, alteration2, 
                     nsamp,
                     g1_freq,  g1_rel_incidence,
                     co_freq, co_rel_incidence)
          
        }
      
    }
    
    
    out <- left_join(out, seer, by='cancer') %>% 
      mutate(co_incidence = round(total_incidence * co_rel_incidence, 0) )
    
    output$usr_tbl <- DT::renderDataTable({
      
      DT::datatable(data.frame(out), 
                    list(pageLength = 15, scrollX=T),rownames = F)
      
    })
    
    tbl_out <<- out
    
  })

observeEvent(input$update_db , {
    
    dbBegin(conn) #start transaction, to keep db safe from disconnects and overflow
    
    ingenename <- input$offsite_genes
    output$update_progress <- renderText({ paste(c('Updating Database with following genes:',ingenename), sep=' ') })


    tryCatch({

	cbio <- cBioPortal()
	studies <- getStudies(cbio, buildReport = TRUE)
	studies <- studies %>% select(studyId) %>% filter(grepl('tcga_pan_can_atlas_2018', studyId))
	ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	translate <- getBM(attributes = c('entrezgene_id','hgnc_symbol'), mart = ensembl)
	incancers <- str_extract(studies$studyId, "^[^_]+") #every cancer we need to download from tcga

	###cnv_main
	output$update_progress <- renderText({'Step 1/8: Downloading data from cbioportal'})
	ingene <- translate %>% filter(hgnc_symbol %in% ingenename) %>% pull(entrezgene_id)
	indat <- c()
	#tcga
	for (incancer in incancers){
	  tryCatch({
	    tmp <- getDataByGenes(cBioPortal(),
				  genes = ingene,
				  molecularProfileIds = paste0(incancer,'_tcga_pan_can_atlas_2018_gistic'),
				  sampleListId = paste0(incancer,'_tcga_pan_can_atlas_2018_all'))
	    tmp <- tmp[[1]] %>% select(uniquePatientKey, hugoGeneSymbol, value) %>%
	      rename('alteration'=value, 'gene'=hugoGeneSymbol) %>%
	      mutate(cancer=str_to_upper(incancer))
	    indat <- rbind(indat,tmp)
	    print(paste0(incancer,': cnv download successful'))
	  }, error=function(e){
		print(e)
	  })
	}
	output$update_progress <- renderText({'Step 2/8: Updating cnv_main'})
	if (nrow(indat)>0){
	    dbAppendTable(conn, 'cnv_main',indat)
	} else {
            stop('no data found??? we should just exit then')
	}
	#origimed
	#TODO: add origimed data

	###cnv_snp

	output$update_progress <-  renderText({'Step 3/8: Updating cnv_snp'})
	indat <- c()
	for (incancer in incancers){
	  try({
	    tmp <- getDataByGenes(cBioPortal(),
				genes = ingene,
				molecularProfileIds = paste0(incancer,'_tcga_pan_can_atlas_2018_mutations'),
				sampleListId = paste0(incancer,'_tcga_pan_can_atlas_2018_cnaseq'))
	    tmp <- tmp[[1]] %>%
	      filter(variantType=='SNP') %>%
	      select(uniquePatientKey, hugoGeneSymbol, referenceAllele, variantAllele, chr, startPosition, endPosition) %>%
	      rename('gene'=hugoGeneSymbol) %>%
	      mutate(cancer=str_to_upper(incancer))

	    if (nrow(tmp)>0){
	      indat <- rbind(indat, tmp)
	    }
	    print(paste0(incancer,': snp download successful'))
	  })
	}
	dbAppendTable(conn, 'cnv_snp',indat)

	###cnv_nsamps
	output$update_progress <- renderText({'Step 4/8: Updating cnv_nsamps'})
	tdat <- tbl(conn, 'cnv_main') %>%
	  select(uniquePatientKey, cancer, gene) %>%
	  filter(gene %in% ingenename) %>%
	  distinct() %>%
	  group_by(cancer, gene) %>%
	  summarize(nsamp=n()) %>%
	  as.data.frame()
	dbAppendTable(conn, 'cnv_nsamps', tdat)

	###single_cnv
	output$update_progress <- renderText({'Step 5/8: Updating single_cnv'})
	tdat <- tbl(conn, 'cnv_main') %>%
	  filter(gene %in% ingenename) %>%
	  group_by(cancer, gene, alteration) %>% summarise(freq=n()) %>%
	  mutate(nsamp=sum(freq), rel_incidence=freq/nsamp) %>%
	  as.data.frame()
	dbAppendTable(conn, 'single_cnv', tdat)

	###co_cnv
	output$update_progress <- renderText({'Step 6/8: Updating co_cnv'})
	out <- c()
	alts <- c(-2, -1, 1, 2) #no need to count 0s
	for (incancer in incancers){ #iter cancers
	  for (igene in ingenename){ #iter genes
	    for (ialt in alts){ #each alteration needs its own relation
	      #first gather the patients who have this gene altered in this way
	      alt_patients <- tbl(conn, 'cnv_main') %>%
		filter(cancer==toupper(incancer), alteration==ialt, gene==igene) %>%
		distinct() %>%
		pull(uniquePatientKey)

	      #second gather those patients' co-altered genes
	      tmp_co <- tbl(conn, 'cnv_main') %>%
		filter(uniquePatientKey %in% alt_patients, gene!=igene, alteration!=0) %>%
		group_by(alteration, gene) %>%
		summarise(freq=n())

	      #lastly, put it all together and write it to the table
	      #forwards
	      tmp <- data.frame(tmp_co) %>%
		rename('gene2'=gene, 'alteration2'=alteration) %>%
		mutate(gene1=igene, alteration1=ialt, cancer=toupper(incancer))
	      #backwards
	      tmp2 <- tmp %>%
		rename('tgene'=gene1, 'talt'=alteration1) %>%
		rename('gene1'=gene2, 'alteration1'=alteration2) %>%
		rename('gene2'=tgene, 'alteration2'=talt)

	      tmp <- rbind(tmp, tmp2)
	      if (nrow(tmp)>0){
		out <- rbind(out, tmp)
	      }
	    }

	    try({if (igene==ingenename[floor(length(ingenename)/2)]){
	      print(paste0(incancer, ' 50% Uploaded'))
	    } })

	  }
	  print(paste0(incancer, ' 100% Uploaded'))
	}
	dbAppendTable(conn, 'co_cnv', out)

	###genes table
	output$update_progress <- renderText({'Step 7/8: Updating gene lists'})
	dbAppendTable(conn, 'genes', tibble(genes=ingenename))

	###offsite_genes
	print('Updating offsite_genes')
	gene_list <- paste0("'", ingenename, "'", collapse = ", ")
	query <- sprintf(
	  "DELETE FROM offsite_genes WHERE newgenes IN (%s);",
	  gene_list
	)
	dbExecute(conn, query)

      output$update_progress <- renderText({'Step 8/8: Committing changes to database'})
      dbCommit(conn) #if everything ran, save the changes
      output$update_progress <- renderText({'Update Complete!'})
    }, error = function(e){
      stop(paste0("Update Database FAILED: ", e))
      dbRollback(conn) #if some error occurred, return to before changes
      output$update_progress <- renderText({'Update Failed :c'})
    })
    #update counters
    ngenes <<- tbl(conn, 'genes') %>% pull(genes) %>% length()    
    noffgenes <<- tbl(conn, 'offsite_genes') %>% pull(newgenes) %>% length()    
    output$ngenes <- renderText({ paste0("Number of extant genes: ", ngenes )})
    output$noffgenes <- renderText({ paste0("Number of offsite genes: ",  noffgenes )})

  })

})


shinyApp(ui = ui, server = server)
