#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

fluidPage(
  
  titlePanel("Chr6 Project"),
  
  fluidRow(
    column(2,
           fileInput(
             "nodefile",
             "Choose node file",
             accept = c("text/csv", ".csv")
           ),
           fileInput(
             "edgefile",
             "Choose edge file",
             accept = c("text/csv", ".csv")
           ),
           hr(),
           uiOutput("patient_focus"),
           uiOutput("gene_focus"),
           uiOutput("range_focus"),
           uiOutput("range_submit"),
           # uiOutput("edge_stuff"),
           # selectInput("Focus", "Focus on node:",
           #             c(1:3)),
           uiOutput("gene_sim_filter"),
           uiOutput("hi_gene_sim_filter"),
           hr(),
           uiOutput("hi_filter"),
           hr(),
           uiOutput("zoom_toggle"),
           uiOutput("dominant_toggle"),
           uiOutput("update_button"),
           uiOutput("generate_button"),
           style="overflow-x: scroll; overflow-y: scroll"),
    column(7,
           tabsetPanel(
             tabPanel("Network", visNetworkOutput("network_proxy", height = "800px")),
             tabPanel("Patient Table", dataTableOutput("patient_node_data")),
             tabPanel("Gene Table", dataTableOutput("gene_node_data")),
             type = "tabs"
           )),
    column(3,
           dataTableOutput("node_subnet_data"),
           # textInput("test", "Test"),
           style="overflow-x: scroll; overflow-y: scroll")
    )
)



# shinyUI(fluidUI(
#   
#   
#   
#   
#   
#   pageWithSidebar(
#   headerPanel("Chr6 Project"),
# 
#   sidebarPanel(
#     fileInput(
#       "nodefile",
#       "Choose node file",
#       accept = c("text/csv", ".csv")
#     ),
#     fileInput(
#       "edgefile",
#       "Choose edge file",
#       accept = c("text/csv", ".csv")
#     ),
#     hr(),
#     uiOutput("patient_focus"),
#     uiOutput("gene_focus"),
#     uiOutput("range_focus"),
#     uiOutput("range_submit"),
#     # uiOutput("edge_stuff"),
#     # selectInput("Focus", "Focus on node:",
#     #             c(1:3)),
#     uiOutput("sim_filter"),
#     hr(),
#     uiOutput("hi_filter"),
#     hr(),
#     uiOutput("zoom_toggle"),
#     uiOutput("dominant_toggle"),
#     uiOutput("update_button"),
#     uiOutput("generate_button")
#   ),
#   
#   mainPanel(
#     tabsetPanel(
#       tabPanel("Network", visNetworkOutput("network_proxy", height = "800px")),
#       tabPanel("Patient Table", dataTableOutput("patient_node_data")),
#       tabPanel("Gene Table", dataTableOutput("gene_node_data")),
#       type = "tabs"
#     ),
#     # visNetworkOutput("network_proxy_nodes", height = "600px")
#   )
# ))
#
# # Define UI for application that draws a histogram
# shinyUI(fluidPage(
#
#   # Application title
#   titlePanel("The Chromosome 6 Project"),
#
#   # # Sidebar with a slider input for number of bins
#   # sidebarLayout(
#   #     sidebarPanel(
#   #         sliderInput("bins",
#   #                     "Number of bins:",
#   #                     min = 1,
#   #                     max = 50,
#   #                     value = 30)
#   #     ),
#   #
#   #     # Show a plot of the generated distribution
#   #     mainPanel(
#   #         plotOutput("distPlot")
#   #     )
#   # ),
#
#   fluidRow(
#     column(
#       width = 4,
#       # selectInput("color", "Color: ", c("blue", "red", "green")),
#       # selectInput("Focus", "Focus on node: ", sort(nodes$id))
#     ),
#     column(
#       width = 8,
#       visNetworkOutput("network_proxy_nodes", height = "400px")
#     )
#   )
# ))
# # visNetworkOutput("network_proxy_nodes", height = "400px")
