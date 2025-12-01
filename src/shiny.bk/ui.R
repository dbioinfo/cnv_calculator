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


options(shiny.maxRequestSize = 10000*1024^2)


ui <- shinyUI(navbarPage(title = "CNV Incidence Calculator",
                         theme=shinytheme('flatly'),
                         tabsetPanel(
                           tabPanel("Single Gene",
                            sidebarLayout(
                              sidebarPanel(
                                fileInput("datfile", label = "Locate RDS object"),
                                actionButton('load_data','Load Database'), 
                                hr(),
                                uiOutput('gene_select'), 
                                hr(),
                                DT::dataTableOutput("top_cnv")
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
                                                       downloadButton('download_igene_hist','Download Figure',align='bottom')
                                                       )
                                        )
                              )
                           ),
                           tabPanel("Co-occurrence",
                                    sidebarLayout(
                                      sidebarPanel(
                                        radioButtons('del_or_dup','Choose Alteration',choices=c('Deletion','Duplication'),selected = 'Deletion'),
                                        uiOutput('gene_co_list'),
                                        actionButton('load_cograph','Load Co-occurrence Plot'),
                                        DT::dataTableOutput('cotable')#data table
                                        ),
                                      mainPanel(
                                        fluidRow(align='center',
                                                 plotOutput('geneco_graph')
                                        ), 
                                        fluidRow(align='center', 
                                                 br(), br(), br(), br(),
                                                 br(), br(), br(), br(),
                                                 br(), br(), br(), br(),
                                                 br(), br()
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