#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(visNetwork)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  # DATA UPLOAD ==================================================================
  nodedata <- reactive({
    infile <- input$nodefile
    if (is.null(infile)) {
      return(NULL)
    }
    read.csv(infile$datapath)
  })

  edgedata <- reactive({
    infile <- input$edgefile
    if (is.null(infile)) {
      return(NULL)
    }
    read.csv(infile$datapath)
  })

  is_data_loaded <- function() {
    net_nodes <- nodedata()
    net_edges <- edgedata()
    if (is.null(net_nodes) | is.null(net_edges)) {
      return(FALSE)
    }
    else {
      return(TRUE)
    }
  }


  # UI ===========================================================================
  # __Focus select ----
  output$patient_focus <- renderUI({
    if (!is_data_loaded()) {
      return(NULL)
    }
    net_nodes <- nodedata()
    node_labels <- c("", net_nodes[net_nodes$group != "HI Gene", "label"])
    names(node_labels) <- node_labels
    selectInput("patient_pick", paste("Focus on patient: (", length(node_labels), ")"), node_labels)
  })

  output$gene_focus <- renderUI({
    if (!is_data_loaded()) {
      return(NULL)
    }
    net_nodes <- nodedata()
    node_labels <- c("", net_nodes[net_nodes$group == "HI Gene", "label"])
    names(node_labels) <- node_labels
    selectInput("gene_pick", paste("Focus on gene: (", length(node_labels), ")"), node_labels)
  })

  output$zoom_toggle <- renderUI({
    if (!is_data_loaded()) {
      return(NULL)
    }
    checkboxInput("zoom_button", "Zoom on focus", FALSE)
  })
  
  output$dominant_toggle <- renderUI({
    if (!is_data_loaded()) {
      return(NULL)
    }
    checkboxInput("dominant_button", "Toggle patients with dominant genes", FALSE)
  })

  # __Jaccard Gene Slider ----
  output$sim_filter <- renderUI({
    if (!is_data_loaded()) {
      return(NULL)
    }
    df <- edgedata()
    max_val <- ceiling(max(df$gene_sim) * 100)
    sliderInput(
      inputId = "sim_slider",
      label = "Gene Similarity Filter:",
      min = 0,
      max = 100,
      value = c(0, 100),
      round = TRUE,
      ticks = TRUE,
    )
  })

  # __HI Score Slider ----
  output$hi_filter <- renderUI({
    if (!is_data_loaded()) {
      return(NULL)
    }
    df <- edgedata()
    sliderInput(
      inputId = "hi_slider",
      label = "HI Score Gene Filter:",
      min = 0,
      max = 100,
      value = c(0, 100),
      round = TRUE,
      ticks = TRUE,
    )
  })

  # __Update button ----
  output$update_button <- renderUI({
    if (!is_data_loaded()) {
      return(NULL)
    }
    actionButton("update", "Update")
  })

  # __Generate button ----
  output$generate_button <- renderUI({
    if (!is_data_loaded()) {
      return(NULL)
    }
    actionButton("generate", "Generate")
  })
  
  # __Range look-up ----
  output$range_focus <- renderUI({
    if (!is_data_loaded()) {
      return(NULL)
    }
    textInput("range_field", "Filter by locus:", placeholder = "6:1000000-3000000")
  })
  
  # __Range submit ----
  output$range_submit <- renderUI({
    if (!is_data_loaded()) {
      return(NULL)
    }
    actionButton("range_button", "Go")
  })

  # OBSERVERS ====================================================================
  # __Initialize Observers ----
  observeEvent(input$generate, {
      focus_pat_obs$resume()
      focus_gene_obs$resume()
    }, once = TRUE)

  # __Patient Focus Observer ----
  focus_pat_obs <- observe(
    {
      # check_ready_obs$destroy()
      if (input$zoom_button) {
        visNetworkProxy("network_proxy") %>%
          visFocus(id = input$patient_pick, scale = 1) %>%
          visSelectNodes(id = input$patient_pick)
      }
      else {
        visNetworkProxy("network_proxy") %>%
          visSelectNodes(id = input$patient_pick)
      }
    },
    suspended = TRUE
  )

  # __Gene Focus Observer ----
  focus_gene_obs <- observe(
    {
      # check_ready_obs$destroy()
      if (input$zoom_button) {
        visNetworkProxy("network_proxy") %>%
          visFocus(id = input$gene_pick, scale = 1) %>%
          visSelectNodes(id = input$gene_pick)
      }
      else {
        visNetworkProxy("network_proxy") %>%
          visSelectNodes(id = input$gene_pick)
      }
    },
    suspended = TRUE
  )
  
  split_locus <- function(locus) {
    locus <- strsplit(locus, ":")[[1]]
    chromosome <- locus[1]
    start <- as.numeric(strsplit(locus[2], "-")[[1]][1])
    stop <- as.numeric(strsplit(locus[2], "-")[[1]][2])
    return(list(chromosome, start, stop))
  }
  
  # __Range update ----
  # observeEvent(input$range_button, {
  #   if (!is.null(input$range_field)) {
  #     locus_search <- split_locus(input$range_field)
  #     nodes <- nodedata()
  #     node_ids <- 
  #     for (node in nodes) {
  #       loci <- strsplit(node$ranges, 
  #     }
  #   }
  # })
  
  # __Dominant update ----
  observeEvent(input$dominant_button, {
    nodes <- nodedata()
    edges <- edgedata()
    
    dom_nodes <- nodes[nodes$group != "HI Gene" & nodes$hi > 0,]
    
    if (input$dominant_button) {
      visNetworkProxy("network_proxy") %>%
        visRemoveNodes(dom_nodes$id)
      new_choices <- nodes$id[nodes$group != "HI Gene" & !nodes$id %in% dom_nodes$id]
      updateSelectInput(session, "patient_pick",
                        choices = c("", new_choices),
                        label = paste("Focus on patient: (", length(new_choices), ")")
      )
    }
    else {
      visNetworkProxy("network_proxy") %>%
        visUpdateNodes(dom_nodes)
      new_choices <- nodes$id[nodes$group != "HI Gene"]
      updateSelectInput(session, "patient_pick",
                        choices = c("", new_choices),
                        label = paste("Focus on patient: (", length(new_choices), ")")
      )
    }
  }, ignoreInit = TRUE)
      
    
  # __Update function ----
  observeEvent(input$update, {

    # Get slider values.
    gene_sim <- input$sim_slider
    gene_hi <- input$hi_slider
    
    # Get other values.
    # remove_dominant <- input$dominant_button

    # Get original data.
    nodes <- nodedata()
    edges <- edgedata()

    # List components to remove.
    remove_edges <- edges$id[(edges$gene_sim * 100) < gene_sim[1] | gene_sim[2] < (edges$gene_sim * 100)]
    gene_nodes <- nodes[nodes$group == "HI Gene",]
    remove_nodes <- gene_nodes$id[(gene_nodes$hi) < gene_hi[1] | gene_hi[2] < (gene_nodes$hi)]
    # if (remove_dominant) {
      # remove_patients <- nodes$id[nodes$group != "HI Gene" & nodes$hi > 0]
      # remove_nodes <- c(remove_nodes, remove_patients)
    # }

    # Make new dropdown list.
    new_choices <- nodes[nodes$group == "HI Gene" & !nodes$id %in% remove_nodes, "label"]

    # Update dropdowns.
    updateSelectInput(session, "gene_pick",
      choices = c("", new_choices),
      label = paste("Focus on gene: (", length(new_choices), ")")
    )

    # Get current data.
    # current_nodes <- visGetNodes(visNetworkProxy("network_proxy"))
    # current_edges <- visGetEdges(visNetworkProxy("network_proxy"))
    # print(class(current_nodes))
    # View(current_nodes)
    # View(current_edges)

    # Update.
    visNetworkProxy("network_proxy") %>%
      visUpdateNodes(nodes = nodes) %>%
      visUpdateEdges(edges = edges) %>%
      visRemoveNodes(remove_nodes) %>%
      visRemoveEdges(remove_edges)
  })

  # check_ready_obs
  focus_pat_obs
  focus_gene_obs

  # MAIN NETWORK =================================================================

  main_graph <- eventReactive(input$generate, {
    net_nodes <- nodedata()
    net_edges <- edgedata()
    visNetwork(net_nodes, net_edges) %>%
      visOptions(
        highlightNearest = list(enabled = T, hover = F),
        nodesIdSelection = T, selectedBy = "group"
      ) %>%
      visEdges(smooth = F) %>%
      visPhysics(solver = "hierarchicalRepulsion", stabilization = F) %>%
      # visIgraphLayout() %>%
      visInteraction(navigationButtons = TRUE)
  })

  output$network_proxy <- renderVisNetwork({
    main_graph()
  })


  observeEvent(input$update | input$generate, {
    visNetworkProxy("network_proxy") %>%
      visGetNodes(addCoordinates = FALSE)
  })

  output$patient_node_data <- renderDataTable({
    if (!is.null(input$network_proxy_nodes)) {
      info <- data.frame(matrix(unlist(input$network_proxy_nodes),
        ncol = ncol(nodedata()),
        byrow = T
      ), stringsAsFactors = FALSE)
      colnames(info) <- colnames(nodedata())
      info <- info[info$group != "HI Gene", -8 ]
      # View(info)
      info
    }
  })

  output$gene_node_data <- renderDataTable({
    if (!is.null(input$network_proxy_nodes)) {
      info <- data.frame(matrix(unlist(input$network_proxy_nodes),
        ncol = ncol(nodedata()),
        byrow = T
      ), stringsAsFactors = FALSE)
      colnames(info) <- colnames(nodedata())
      info <- info[info$group == "HI Gene", -8 ]
      info
    }
  })
})
