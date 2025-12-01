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


options(shiny.maxRequestSize = 10000*1024^2)


ui <- shinyUI(navbarPage(title = "CNV Incidence Calculator 1.3", 
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
                                  column(12,
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
                                                           column(3, radioButtons('g2_tbl','Gene 2 Variant',choices=c('CNV','SNP'),selected = 'CNV')),
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
                           )
                         )
))


