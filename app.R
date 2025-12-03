options(warn = -1)
set.seed(001203)

library(shiny)
library(shinyjs)
library(Seurat)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(purrr)
library(tools)
library(patchwork)
library(cowplot)
library(harmony)

ui <- fluidPage(
  tags$head(
      tags$style(HTML("
    .harmony-checkbox label {
      color: red !important;
      font-weight: bold !important;
      margin-left: 4px;
    }
  ")),
    tags$style(HTML("
      .tooltip-circle {
        position: relative;
        display: inline-flex;
        justify-content: center;
        align-items: center;
        text-align: justify;
        width: 18px;
        height: 18px;
        background-color: grey;
        border-radius: 50%;
        color: white;
        font-weight: normal;
        font-size: 13px;
        cursor: pointer;
        margin: 10px;
      }

      .tooltip-circle .tooltip-text {
        visibility: hidden;
        background-color: white;
        color: black;
        border: 1px solid #d3d3d3;
        padding: 6px;
        border-radius: 6px;
        position: absolute;
        z-index: 1000;
        top: -5px;
        left: 30px;
        white-space: normal;
        width: 300px;
      }

      /* title */
      .tooltip-title {
        text-transform: uppercase;
        font-weight: bold;
        margin-bottom: 4px;
      }

      /* line */
      .tooltip-divider {
        height: 1px;
        background-color: #B0D4FF; 
        margin: 6px 0;
      }

      /* show */
      .tooltip-circle:hover .tooltip-text {
        visibility: visible;
      }
    "))
  ),
  titlePanel("Clustering Analysis for Single-Cell RNA-seq Data"),
  h4("Produced by Bingtian Ye (btye@bu.edu) in Waxman's Lab"),
  hr(),
  h4("Enter path to RDS (Seurat object) file"),
  fluidRow(
    column(8, uiOutput("data_path_input")),
    column(4, uiOutput("read_reset_ui"))
  ),
  hr(),
  uiOutput("control_ui"),
  tags$div(style = "margin-bottom: 5em;"),
  uiOutput("save_btn_panel")
)

server <- function(input, output, session) {
  seuratObj <- reactiveVal(NULL)
  variable_method <- reactiveVal(NULL)
  resolution_search_results <- reactiveVal(NULL)
  umap_plot <- reactiveVal(NULL)
  heatmap_plot <- reactiveVal(NULL)
  cluster_table <- reactiveVal(NULL)
  dotplot1_plot <- reactiveVal(NULL)
  dotplot2_plot <- reactiveVal(NULL)
  current_ident <- reactiveVal(NULL)
  subset_data_path <- reactiveVal(NULL)
  plot_title <- reactiveVal(NULL)
  harmony_state <- reactiveVal(NULL)
  harmony_meta_state <- reactiveVal(NULL)
  
  # CONTROL STEP STATUS
  status <- reactiveValues(
    normalized = FALSE,
    scaled = FALSE,
    pca = FALSE,
    neighbor = FALSE,
    clustered = FALSE
  )
  
  output$data_path_input <- renderUI({
    if (is.null(seuratObj())) {
      textInput(
        "data_path", NULL,
        value = "",
        width = "100%"
      )
    } else {
      tags$div(
        style = "padding: 8px; background-color: #e8f5e8; border-radius: 4px; border: 1px solid #c8e6c9;",
        tags$p(
          "📁 ",
          tags$strong("Loaded: "),
          if (is.null(subset_data_path())) {
            # normal: use input$data_path
            tags$code(input$data_path)
          } else {
            # after subset: show modified message
            tagList(
              "Subset data from ",
              tags$code(subset_data_path())
            )
          },
          style = "margin: 0; color: #2e7d32;"
        )
      )
    }
  })
  
  
  # READ or RESET BUTTON UI
  output$read_reset_ui <- renderUI({
    if (is.null(seuratObj())) {
      actionButton("read_btn_ui", "READ DATA", width = "100%", class = "btn-primary")
    } else {
      actionButton("reset_btn_ui", "RESET", width = "100%", class = "btn-warning")
    }
  })
  
  # READ DATA
  observeEvent(input$read_btn_ui, {
    req(input$data_path)
    path <- input$data_path
    if (!nzchar(path) || !file.exists(path)) {
      showModal(modalDialog(title = "Error", "Please enter a valid existing path!", easyClose = TRUE))
      return()
    }
    showModal(modalDialog(title = "Please wait", 
                          tagList(
                            p("Reading..."),
                            p("The waiting time is related to the size of the RDS file and usually takes one minute.")
                          ), footer = NULL, easyClose = FALSE))
    tryCatch({
      obj <- readRDS(path)
      DefaultAssay(obj) <- "RNA"  
      seuratObj(obj) 
      removeModal()
      showNotification("RDS file loaded.", type = "message")
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
      seuratObj(NULL)
    })
  })
  
  # RESET DATA
  observeEvent(input$reset_btn_ui, {
    seuratObj(NULL)
    status$normalized <- FALSE
    status$scaled <- FALSE
    status$pca <- FALSE
    status$neighbor <- FALSE
    status$clustered <- FALSE
    showNotification("Reset completed. You can now load a new dataset.", type = "message")
  })
  
  # save data
  output$save_btn_panel <- renderUI({
    req(seuratObj())
    absolutePanel(
      fixed = TRUE,
      draggable = FALSE,
      bottom = 20,
      left = 30,
      width = 120,
      height = "auto",
      style = "z-index: 1000;",
      actionButton(
        "save_data_btn",
        "SAVE DATA",
        width = "100%",
        style = "
        background-color: green;
        color: white;
        border: none;
        font-weight: bold;
      "
      )
    )
  })
  
  observeEvent(input$save_data_btn, {
    showModal(modalDialog(
      size = "m",
      title = "Save Data",
      fluidRow(column(12,
          textInput("save_path", "Input save path on SCC or local file name (enter the save path using forward slashes (/). The filename should end with '.rds')",
                    value = paste0(file_path_sans_ext(input$data_path), "_", Sys.Date(), ".rds"), width = "100%"))),
      footer = tagList(
        div(style = "display: flex; width: 100%; gap: 8px;",
            actionButton("btn_save", "Save", style = "
            flex: 1;
            background-color: #007bff; 
            color: white;
            border: none;
          "),
           downloadButton("downloadData", "Download Locally", style = "
            flex: 1;
            background-color: green;  
            color: white;
            border: none;
          "),
           actionButton("close_modal", "Close", style = "
          flex: 1;
          background-color: #dc3545;  
          color: white;
          border: none;
        ")
        )
      )
    ))
  })
  
  observeEvent(input$btn_save, {
    req(seuratObj(), input$save_path)
    save_path <- input$save_path
    showModal(modalDialog(title = "Please wait", 
                          tagList(
                            p("Saving..."),
                            p("The waiting time is related to the size of the dataset. Be patience...")
                          ), footer = NULL, easyClose = FALSE))
    tryCatch({
      saveRDS(seuratObj(), file = save_path)
      removeModal()
      showNotification(paste("Data saved to", save_path), type = "message")
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      basename(input$save_path) %||% paste0("data_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(seuratObj())
      saveRDS(seuratObj(), file)
    },
    contentType = "application/octet-stream"
  )
  observeEvent(input$close_modal, {
    removeModal()
  })
  
  
  # CONTROL ALL STEPS UI
  output$control_ui <- renderUI({
    req(seuratObj())
    
    layers <- Layers(seuratObj(), assay = "RNA")
    
    stop_flag <- FALSE 
    
    ui_list <- list()
    
    # --- Normalization ---
    if (is.null(seuratObj())) {
      ui_list <- append(ui_list, list(
        tags$p("🔴 DATA NORMALIZED", style = "color:red; font-weight:bold; font-size:14px;")
      ))
      stop_flag <- TRUE
    } else {
      ui_list <- append(ui_list, list(
        div(
          style = "display:flex; align-items:center; gap:6px;",
          if (isTRUE(status$normalized)) {
            tags$p("✅ DATA NORMALIZED", 
                   style = "color:black; font-weight:bold; font-size:14px; margin:0;")
          } else {
            tags$p("🔘 DATA NORMALIZED", 
                   style = "color:black; font-weight:bold; font-size:14px; margin:0;")
          },
          div(class = "tooltip-circle",
              "?",  
              span(class = "tooltip-text", 
                   div(class = "tooltip-title", "Data Normalized"),
                   div(class = "tooltip-divider"),
                   div("This step normalizes the gene expression measurements for each cell by
                       the total expression, multiplies this by a scale factor (default is 10,000),
                       and log-transforms the result. Normalization is essential to make gene expression
                       levels comparable across cells.")
              )
          )
        ),
        
        uiOutput("normalized_data_ui")
      ))
    }
    
    ui_list <- append(ui_list, list(hr()))
    
    # --- Scaling ---
    if (stop_flag || !"data" %in% layers) {
      ui_list <- append(ui_list, list(
        tags$p("🔴 DATA SCALED", style = "color:red; font-weight:bold; font-size:14px; margin:0;")
      ))
      stop_flag <- TRUE
    } else {
      ui_list <- append(ui_list, list(
        div(
          style = "display:flex; align-items:center; gap:6px;",
        if (isTRUE(status$scaled)) {
          tags$p("✅ DATA SCALED", style = "color:black; font-weight:bold; font-size:14px; margin:0;")
        } else {
          tags$p("🔘 DATA SCALED", style = "color:black; font-weight:bold; font-size:14px; margin:0;")
        },
        div(class = "tooltip-circle",
            "?",  
            span(class = "tooltip-text", 
                 div(class = "tooltip-title", "Data Scaled"),
                 div(class = "tooltip-divider"),
                 div("This step scales and centers the gene expression data. Scaling ensures that
                     each gene contributes equally to downstream analyses, such as PCA, by transforming
                     the data to have a mean of zero and a standard deviation of one.")
            )
        )
        ),
        uiOutput("scaled_data_ui")
      ))
    }
    
    ui_list <- append(ui_list, list(hr()))
    
    # --- PCA ---
    if (stop_flag || !"scale.data" %in% layers) {
      ui_list <- append(ui_list, list(
        tags$p("🔴 RUN PCA", style = "color:red; font-weight:bold; font-size:14px; margin:0;")
      ))
      stop_flag <- TRUE
    } else {
      ui_list <- append(ui_list, list(
         div(
          style = "display:flex; align-items:center; gap:6px;",
        if (isTRUE(status$pca)) {
          tags$p("✅ RUN PCA", style = "color:black; font-weight:bold; font-size:14px; margin:0;")
        } else {
          tags$p("🔘 RUN PCA", style = "color:black; font-weight:bold; font-size:14px; margin:0;")
        },
        div(class = "tooltip-circle",
            "?",  
            span(class = "tooltip-text", 
                 div(class = "tooltip-title", "Principal Component Analysis (PCA)"),
                 div(class = "tooltip-divider"),
                 div("PCA (Principal Component Analysis) is a dimensionality reduction method that transforms
                     high-dimensional gene expression data into a lower-dimensional space. This process identifies
                     the principal components (PCs) that capture the greatest variance in the data, which helps
                     with visualizing and clustering cells. The “Number of Principal Components” setting determines
                     how many PCs will be retained. It’s generally best to keep the default value, but if later
                     recommended dimensions suggest same as the PC number you input here, you can consider increasing
                     this number to 50 or even higher.")
            )
        )
        ),
        uiOutput("pca_ui")
      ))
    }
    
    ui_list <- append(ui_list, list(hr()))
    
    # --- Find Neighbors ---
    if (stop_flag || is.null(seuratObj()@reductions$pca) || seuratObj()@reductions$pca@assay.used != "RNA") {
      ui_list <- append(ui_list, list(
        tags$p("🔴 FIND NEIGHBORS", style = "color:red; font-weight:bold; font-size:14px;")
      ))
      stop_flag <- TRUE
    } else {
      ui_list <- append(ui_list, list(
        div(
          style = "display:flex; align-items:center; gap:6px;",
        if (isTRUE(status$neighbor)) {
          tags$p("✅ FIND NEIGHBORS", style = "color:black; font-weight:bold; font-size:14px; margin:0;")
        } else {
          tags$p("🔘 FIND NEIGHBORS", style = "color:black; font-weight:bold; font-size:14px; margin:0;")
        },
        div(class = "tooltip-circle",
            "?",  
            span(class = "tooltip-text", 
                 div(class = "tooltip-title", "Find Neighbors"),
                 div(class = "tooltip-divider"),
                 div("This step constructs a K-nearest neighbor graph based on the PCA-reduced data.
                     It identifies neighboring cells for each cell, which is essential for subsequent
                     clustering analyses. The program will provide a recommended number of PCs for follow-up
                     analysis; it is recommended to use the recommended value.")
            )
        )
        ),
        uiOutput("neighbor_ui")
      ))
    }
    
    ui_list <- append(ui_list, list(hr()))
    
    # --- Find Clusters ---
    if (stop_flag || is.null(seuratObj()@graphs$RNA_nn)) {
      ui_list <- append(ui_list, list(
        tags$p("🔴 FIND CLUSTERS", style = "color:red; font-weight:bold; font-size:14px;")
      ))
      stop_flag <- TRUE
    } else {
      ui_list <- append(ui_list, list(
        div(
          style = "display:flex; align-items:center; gap:6px;",
        if (isTRUE(status$clustered)) {
          tags$p("✅ FIND CLUSTERS", style = "color:black; font-weight:bold; font-size:14px; margin:0;")
        } else {
          tags$p("🔘 FIND CLUSTERS", style = "color:black; font-weight:bold; font-size:14px; margin:0;")
        },
        div(class = "tooltip-circle",
            "?",  
            span(class = "tooltip-text", 
                 div(class = "tooltip-title", "Find Clusters"),
                 div(class = "tooltip-divider"),
                 div("This step performs clustering on the K-nearest neighbor graph constructed in the previous step.
                     It groups cells into clusters based on their similarities in gene expression profiles.
                     After that, it will provide visualization for the clustering result.")
            )
        )
        ),
        uiOutput("clustering_ui")
      ))
    }
    
    ui_list <- append(ui_list, list(hr()))
    
    # --- Rename and Subset ---
    if (stop_flag|| is.null(seuratObj()@meta.data$shiny_clusters)) {
      ui_list <- append(ui_list, list(
        tags$p("🔴 RENAME AND SUBSET", style = "color:red; font-weight:bold; font-size:14px;")
      ))
    } else {
      ui_list <- append(ui_list, list(
        div(
          style = "display:flex; align-items:center; gap:6px;",
          tags$p(
            "🔘 RENAME AND SUBSET",
            style = "color:black; font-weight:bold; font-size:14px; margin:0;"
          ),
          div(
            class = "tooltip-circle",
            "?",
            span(
              class = "tooltip-text",
              div(class = "tooltip-title", "RENAME AND SUBSET"),
              div(class = "tooltip-divider"),
              div(HTML("\"Metadata column name\" is the name in seurat, which is the identification of cluster.<br>
             Rename the corresponding input box, then select only the desired samples in the checkboxes and click \"subset\".<br>
             After subsetting, click \"Analysis on Subset Data\". The system will reset ShinyApp and use the subset data as 
             input for analysis.<br>
             It is recommended to save the data after both renaming and subsetting."))
            )
          )
        ),
        
        uiOutput("final_ui")
      ))
    }
    
    
    tagList(ui_list)
  })

  
  # NORMALIZATION STEP UI
  output$normalized_data_ui <- renderUI({
    req(seuratObj())
    fluidRow(
      column(4,actionButton("normalize_btn_ui", "Normalize Data", width = "200px", class = "btn-primary")))
  })
  
  output$scaled_data_ui <- renderUI({
    tagList(
      div(
        style = "display:flex; align-items:center; gap:6px;",
        tags$h5("Gene Selection Method:", style = "color:black; font-weight:bold; font-size:14px; margin:0;"),
        div(class = "tooltip-circle",
            "?",  
            span(class = "tooltip-text", 
                 div(class = "tooltip-title", "Gene Selection Method"),
                 div(class = "tooltip-divider"),
                 div("This option lets you choose whether to scale only the highly variable genes or all genes. 
                      Using variable features allows the program to focus on the genes that show the strongest and 
                      most informative expression patterns, which generally speeds up computation and highlights 
                      biologically meaningful signals. Scaling all genes offers a more comprehensive approach but 
                      may require significantly more processing time. Note that the scaling step only affects 
                      downstream PCA results.")
            )
        )
      ),
      fluidRow(
        tags$div(
          style = "margin-left:20px;",
          radioButtons(
            "findvariable",
            NULL,
            choices = c(
              "Use Variable Features (Recommended)" = "variable",
              "Use All Genes (Longer Processing Time)" = "all"
            ),
            selected = input$findvariable %||% "variable",
            inline = TRUE
          )
        )
      ),
      fluidRow(
        column(4,
               actionButton("scaling_btn_ui", "Scale Data",
                            width = "200px", class = "btn-primary"))
      )
    )
  })
  
  
  # PCA STEP UI
  output$pca_ui <- renderUI({
    req(seuratObj())
    gene_method <- if (!is.null(input$findvariable)) {
      if (input$findvariable == "variable") "Variable Features" else "All Genes"
    } else {
      "Variable Features" 
    }
    tagList(
      fluidRow(
        tags$div(
          style = "margin-left:20px;",
          h5("Gene Selection Method:", 
             gene_method,
             style = "font-weight: bold;")
        )
      ),
      fluidRow(
        div(
          style = "display: inline-block; vertical-align: top; margin-left:30px;",
          h5("Number of Principal Components (PCs):")
        ),
        div(
          style = "display: inline-block; vertical-align: middle; margin-left:10px;",
          numericInput("npc", NULL, value = 30, min = 1, max = 50, step = 1, width = "100px")
        )
      ),
      fluidRow(
        column(4,
               actionButton("pca_btn_ui", "Run PCA",
                            width = "200px", class = "btn-primary"))
      )
    )
  })
  

  
  # NEIGHBOR STEP UI
  output$neighbor_ui <- renderUI({ 
    req(seuratObj())
    
    recommand_value <- if(!is.null(seuratObj()@reductions$pca)) {
      recommended_dimensions()$recommended
    } else {
      NULL
    }
    
    tagList(
      if (!is.null(seuratObj()@reductions$pca) && isTRUE(status$pca)) {
        tagList(
          div(
            style = "display: flex; align-items: flex-start; gap: 30px; width: 850px;",  
            div(
              style = "width: 400px; height: 400px;",
              h5("Elbow Plot with Recommended Dimensions:"),
              plotOutput("elbow_plot", width = "400px", height = "400px")
            ),
            div(
              style = "width: 400px;",  
              h5("Dimension Selection Criteria"),
              wellPanel(
                style = "background-color: #f7f7f7; width: 400px;",
                tags$p(tags$strong("Recommended dimensions based on:")),
                tags$ul(
                  tags$li(tags$strong("Priority 1:"), " PCs after the elbow point, where variance drop < 0.5%."),
                  tags$li(tags$strong("Priority 2:"), " Cumulative variance ≥ 90%.")
                ),
                tags$p("Elbow Point: ",
                       tags$strong(textOutput("elbow_point", inline = TRUE),
                                   style = "color: grey; font-size: 16px;")),
                tags$p("Last variance drop < 0.5% point: ",
                       tags$strong(textOutput("drop_point", inline = TRUE),
                                   style = "color: grey; font-size: 16px;")),
                tags$p("Cumulative variance ≥ 90% point: ",
                       tags$strong(textOutput("cum_point", inline = TRUE),
                                   style = "color: grey; font-size: 16px;")),
                tags$p("Recommended dimensions: ",
                       tags$strong(textOutput("recommended_dims", inline = TRUE),
                                   style = "color: red; font-size: 16px;"))
              )
            )
          ),br(),br(),br())
      },
      fluidRow(
        div(
          style = "display: inline-block; vertical-align: top; margin-left:20px;",
          h5("Select the number of dimensions (PCs) to use for finding neighbors:",
             style = "font-weight: bold;")
        ),
        div(
          style = "display: inline-block; vertical-align: middle; margin-left:10px;",
          numericInput("nn_dims", NULL,
                       value = recommand_value,
                       min = 1, max = input$npc %||% 50, step = 1, width = "100%"))
      ),
      fluidRow(
        div(
          style = "display: flex; margin-left:30px;",
          div(
            class = "harmony-checkbox",
            checkboxInput(
              inputId = "integration_harmony",
              label = "Integration (Harmony)",
              value = ifelse(!is.null(harmony_state()), harmony_state(), FALSE),
              width = "200px"
            ) %>% 
              tagAppendAttributes(style = "white-space: nowrap; display: inline-flex; align-items:center;")
          ),
          # Tooltip
          div(
            class = "tooltip-circle",
            style = "margin-left:10px;",
            "?",
            span(
              class = "tooltip-text",
              div(class = "tooltip-title", "INTEGRATION SAMPLE (HARMONY)"),
              div(class = "tooltip-divider"),
              div("If your Seurat object contains multiple samples or batches, enabling this option will
             show the Harmony integration parameter. This corrects batch effects and ensures proper
             integration based on biology rather than technical variation.")
            )
          )
        )
      ),
      
      conditionalPanel(
        condition = "input.integration_harmony == true",
          div(
            style = "display: inline-block; vertical-align: middle; margin-left:20px;",
            uiOutput("harmony_metadata_ui"))),
      fluidRow(
        column(4,
               actionButton("neighbor_btn_ui", "Find Neighbors",
                            width = "200px", class = "btn-primary"))
      )
    )
    })
  
  # CLUSTERING STEP UI
  output$clustering_ui <- renderUI({ 
    req(seuratObj())
    tagList(
      tags$div(style = "margin-bottom: 5px;",
               fluidRow(
                 div(style = "display: inline-block; vertical-align: top; margin-left:20px;margin-right:0;",
                     h5("Upload Gene Marker Excel file for cluster annotation (optional):", 
                        style = "font-weight: bold;")),
                 div(class = "tooltip-circle",
                     "?",  
                     span(class = "tooltip-text", 
                          div(class = "tooltip-title", "Customize Your Gene Marker"),
                          div(class = "tooltip-divider"),
                          div("You can upload an Excel file containing your own gene markers for cluster annotation.
                              First, you need download the template file and fill in your gene markers for each 
                              cluster in the provided format. The excel file included three sheets: \"HEATMAP\",
                              \"DOTPLOT1\", \"DOTPLOT2\". \"HEATMAP\" marker will generate heatmap and also 
                              define the cluster name in umap. Please ensure that all markers in the \"HEATMAP\"
                              table exist in the data, and that no duplicate markers appear in different columns
                              within the two \"DOTPLOT\". The first row of each column lists the marker's class,
                              followed by the marker's name. Note the capitalization.")
                     )
                 )
               ),
               fluidRow(
                 div(
                   style = "display: inline-block; vertical-align: top; margin-left:30px;",
                   downloadButton("download_template", "Download Gene Marker Template", 
                                  style = "background-color: #4CAF50; color: white;")
                 ),
                 div(
                   style = "display: inline-block; vertical-align: top; margin-left:40px;",
                   tags$h5("Upload Excel File:", style = "font-weight: bold;")
                 ),
                 div(
                   style = "display: inline-block; vertical-align: top; margin-left:10px; margin-bottom: -40px;",
                   fileInput("upload_excel", NULL,
                             accept = c(".xlsx", ".xls"),
                             buttonLabel = "Select File...",
                             placeholder = "No file selected, Using default",
                             width = "330px")
                 )
               ),
               hr(),
               conditionalPanel(
                 condition = "input.search_resolution == true",
                 tagList(
                   fluidRow(
                     div(
                       style = "display: flex; align-items: center; gap: 6px;",
                       h5("Search resolutions in range for cluster numbers (Longer Processing Time!!!!!):",
                          style = "font-weight: bold; margin-left: 20px; margin-right: 0px;"),
                       div(class = "tooltip-circle",
                           "?",  
                           span(class = "tooltip-text", 
                                div(class = "tooltip-title", "Search resolution range"),
                                div(class = "tooltip-divider"),
                                div("You can specify a range of resolution values to search for optimal clustering.
                   The 'Start' value sets the beginning of the range, 'End' sets the maximum resolution,
                   and 'Step' determines the increment between each resolution value tested.
                   For example, with a Start of 0.01, End of 1, and Step of 0.01, the program will
                   evaluate resolutions from 0.01 to 1 in increments of 0.01. Click 'Search in Shiny' 
                   you will get a list of cluster number under different resolution. Click 'Search on SCC background' 
                   will launch a qsub mission on the background, which will inclueded a list of cluster number and 
                   under visualization result different resolution.")
                   )))),
                   fluidRow(
                     div(
                       style = "display: inline-block; vertical-align: middle; margin-left:20px;",
                       h5("Start :", style = "font-weight: bold; display: inline-block; margin-right:8px;"),
                       numericInput("start_resolution", NULL, 
                                    value = 0.01, min = 0, max = 5.0, step = 0.01, width = "100px")
                     ),
                     div(
                       style = "display: inline-block; vertical-align: middle; margin-left:20px;",
                       h5("Step :", style = "font-weight: bold; display: inline-block; margin-right:8px;"),
                       numericInput("step_resolution", NULL, 
                                    value = 0.01, min = 0.001, max = 1, step = 0.001, width = "100px")
                     ),
                     div(
                       style = "display: inline-block; vertical-align: middle; margin-left:20px;",
                       h5("End :", style = "font-weight: bold; display: inline-block; margin-right:8px;"),
                       numericInput("end_resolution", NULL, 
                                    value = 1, min = 0, max = 5.0, step = 0.01, width = "100px")
                     ),
                     div(
                       style = "display: inline-block; vertical-align: middle; margin-left:30px;",
                       actionButton("search_resolution_btn_ui", 
                                    "Search in Shiny",
                                    width = "150px", class = "btn-primary")
                     ),
                     div(
                       style = "display: inline-block; vertical-align: middle; margin-left:30px;",
                       actionButton("search_on_scc", 
                                    "Search on SCC background",
                                    width = "200px", class = "btn-danger")
                     )
                   ),
                   hr()
                 )
               ),
               fluidRow(
                 div(style = "display: inline-block; vertical-align: top; margin-left:20px;",
                     h5("Select resolution for clustering:",
                        style = "font-weight: bold;")),
                 div(style = "display: inline-block; vertical-align: top; margin-left:10px;",
                     uiOutput("resolution_ui")),
                 div(class = "tooltip-circle",
                     "?",  
                     span(class = "tooltip-text", 
                          div(class = "tooltip-title", "Clustering Resolution"),
                          div(class = "tooltip-divider"),
                          div("The resolution parameter influences the granularity of the clustering results.
                              A higher resolution leads to more clusters (finer granularity), while a lower
                              resolution results in fewer clusters (coarser granularity). It is recommended
                              to start with a resolution of 0.3 to 0.5 and adjust based on the specific dataset
                              and biological context. If you want to see different resolution's cluster number, you 
                              can check \"Search Resolutions\" option to get more detail.")
                     )
                 ),
                 div(style = "display: inline-block; vertical-align: top; margin-left:30px;",
                     checkboxInput("search_resolution", "Search Resolutions", width = "100%")),
                 uiOutput("show_search_result_ui", inline = TRUE)
                 ),
               fluidRow(
                 column(4,
                        actionButton("clustering_btn_ui", "Find Clusters",
                                     width = "200px", class = "btn-primary"))
               ),
               fluidRow(
                 column(12,
                        conditionalPanel(
                          condition = "output.clustered_status == true",
                          uiOutput("umap_ui")
                        ))
               )
      )
    )
  })
  
  output$resolution_ui <- renderUI({
    numericInput("resolution", NULL, 
                 value = isolate(input$resolution) %||% 0.03, 
                 min = 0.01, max = 5.0, step = 0.001, width = "100%")
  })
  
  
  # FINAL STEP UI
  output$final_ui <- renderUI({
    req(seuratObj(),input$resolution)
    if(is.null(current_ident())){
      val <- "shiny_clusters"
      current_ident(val)
    }
    tagList(
      fluidRow(
        div(style = "display: inline-block; vertical-align: top; margin-left:20px;",
            h5("Metadata column name:",
               style = "font-weight: bold;")),
        div(style = "display: inline-block; vertical-align: top; margin-left:20px;",
            h5(current_ident(),
               style = "font-weight: bold; color: red;")),
        div(style = "display: inline-block; vertical-align: top; margin-left:10px;",
            textInput("rename_idents_input", NULL, value = current_ident(),
              width = "200px"
            )),
        div(style = "display: inline-block; vertical-align: top; margin-left:30px;",
            actionButton("rename_idents_btn","Rename Column",
              width = "150px",class = "btn-primary"))
      ),
      uiOutput("sample_subset_rename_ui"),
      h5("Attention: Please using \"SAVE DATA\" to save your current result before subset!!!", style = "font-weight: bold; color:red;"),
      fluidRow(
        div(style = "display: inline-block; vertical-align: top; margin-left:30px;",
          actionButton("rename_btn_ui", "Rename",
                       width = "200px", class = "btn-primary")),
      div(style = "display: inline-block; vertical-align: top; margin-left:30px;",
          actionButton("subset_btn_ui", "Subset",
                       width = "200px", class = "btn-danger")),
      div(style = "display: inline-block; vertical-align: top; margin-left:30px;",
          conditionalPanel(condition = "input.confirm_subset > 0",
            actionButton("keep_analysis", "Analysis on Subset Data",
                       width = "200px", class = "btn-success")))
    ))
  })
  
  observeEvent(input$search_on_scc, {
    showModal(
      modalDialog(
        title = tags$strong("Run SCC Job in Background"),
        easyClose = TRUE,
        
        # ---- Main text ----
        tags$p(
          "This will run an SCC background job. After submission, SCC will generate analysis results and figures for multiple resolutions.",
          style = "font-size:15px; font-weight:500;"
        ),
        tags$p(
          "Search range use the parameter in Shiny, please edit shiny start, step and end step first.",
          style = "font-size:15px; font-weight:500;"
        ),
        # Save path
        textInput(
          "scc_save_path", 
          "Save path on scc:", 
          placeholder = "/projectnb/your_path/..."
        ),
        
        # Project name
        textInput(
          "scc_project_name", 
          "Project name:", 
          placeholder = "e.g., My_SingleCell_Project"
        ),
        
        # Runtime (hours) + number of cores
        fluidRow(
          column(
            6,
            numericInput(
              "scc_runtime",
              "Runtime (hours):",
              value = 12,
              min = 1,
              step = 1
            )
          ),
          column(
            6,
            numericInput(
              "scc_cores",
              "Number of cores:",
              value = 8,
              min = 1,
              step = 1
            )
          )
        ),
        
        # Optional email
        textInput(
          "scc_email", 
          "Email (optional):", 
          placeholder = "your_email@bu.edu"
        ),
        
        br(),
        
        # ---- Footer buttons ----
        footer = tagList(
          modalButton("Cancel"),
          actionButton("run_scc_job", "Run", class = "btn-primary")
        )
      )
    )
  })
  
  observeEvent(input$run_scc_job, {
    req(seuratObj(), input$scc_save_path, input$scc_project_name, input$scc_runtime, input$scc_cores, input$start_resolution, input$step_resolution, input$end_resolution)
    
    showModal(modalDialog(
      title = "Please wait",
      tagList(
        p("Submitting SCC job..."),
        p("You will receive an email when the job is completed.")
      ),
      footer = NULL,
      easyClose = FALSE
    ))
    
  
    
    tryCatch({
      # Call function to submit SCC job
      source("/projectnb/wax-es/00_shinyapp/Clustering/Clustering_shinyapp/scc_job_submission.R")
      submit_scc_job(
        seuratObj(),
        save_path = input$scc_save_path,
        project_name = input$scc_project_name,
        start_resolution = input$start_resolution,
        step_resolution = input$step_resolution,
        end_resolution = input$end_resolution,
        runtime = input$scc_runtime,
        cores = input$scc_cores,
        email = input$scc_email
      )
      
      removeModal()
      showNotification("SCC job submitted successfully.", type = "message")
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  
                 
  observeEvent(input$rename_idents_btn,{
    req(seuratObj(), input$rename_idents_input)
    current_ident <- current_ident()
    new_ident <- input$rename_idents_input
    tryCatch({
      if (new_ident == "") {
        showNotification("Please enter a new ident name", type = "error")
        return()
      }
      if (new_ident %in% names(seuratObj()@meta.data)) {
        showNotification(paste("Ident name", new_ident, "already exists"), type = "error")
        return()
      }
      seurat_obj <- seuratObj()
      seurat_obj@meta.data[[new_ident]] <- seurat_obj@meta.data[[current_ident]]
      seurat_obj@meta.data[[current_ident]] <- NULL
      Idents(seurat_obj) <- seurat_obj@meta.data[[new_ident]]
      seuratObj(seurat_obj)
      showNotification(paste("Ident renamed to", new_ident), type = "message")
      current_ident(new_ident) 

    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
    
  })
  
  # NORMALIZATION STEP SERVER
  observeEvent(input$normalize_btn_ui,{
    req(seuratObj())
    showModal(modalDialog(title = "Please wait", 
                          tagList(
                            p("Normalizing..."),
                            p("The waiting time is related to the size of the dataset and usually takes one minute.")
                          ), footer = NULL, easyClose = FALSE))
    tryCatch({
      seuratObj(NormalizeData(seuratObj(), verbose = FALSE))
      status$normalized <- TRUE
      status$scaled <- FALSE
      status$pca <- FALSE
      status$neighbor <- FALSE
      status$clustered <- FALSE
      removeModal()
      showNotification("Data normalized.", type = "message")
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  
  observeEvent(input$scaling_btn_ui,{
    req(seuratObj())
    showModal(modalDialog(title = "Please wait", 
                          tagList(
                            p("Scaling..."),
                            p("The waiting time is related to the size of the dataset and usually takes one minute.")
                          ), footer = NULL, easyClose = FALSE))
    tryCatch({
      if (input$findvariable == "variable") {
        seuratObj(FindVariableFeatures(seuratObj(), selection.method = "vst", nfeatures = 2000, verbose = FALSE))
        seuratObj(ScaleData(seuratObj(), features = VariableFeatures(seuratObj()), verbose = FALSE))
      } else {
        seuratObj(ScaleData(seuratObj(), features = rownames(seuratObj()), verbose = FALSE))
      }
      status$normalized <- TRUE
      status$scaled <- TRUE
      status$pca <- FALSE
      status$neighbor <- FALSE
      status$clustered <- FALSE
      variable_method(if(input$findvariable == "variable") "Variable Features" else "All Genes")
      removeModal()
      showNotification("Data scaled.", type = "message")
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  
  # PCA STEP SERVER
  observeEvent(input$pca_btn_ui,{
    req(seuratObj())
    showModal(modalDialog(title = "Please wait", 
                          tagList(
                            p("Running PCA..."),
                            p("The waiting time is related to the size of the dataset, Gene Selection Method and the number of PCs.")
                          ), footer = NULL, easyClose = FALSE))
    tryCatch({
      if (input$findvariable == "variable") {
        seuratObj(suppressWarnings(RunPCA(seuratObj(), features = VariableFeatures(seuratObj()), npcs = input$npc, verbose = FALSE)))
      } else {
        seuratObj(suppressWarnings(RunPCA(seuratObj(), features = rownames(seuratObj()), npcs = input$npc, verbose = FALSE)))
      }
      status$normalized <- TRUE
      status$scaled <- TRUE
      status$pca <- TRUE
      status$neighbor <- FALSE
      status$clustered <- FALSE
      removeModal()
      showNotification("PCA completed.", type = "message")
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  
  # ======== PCA recommendation calculation ========
  
  recommended_dimensions <- reactive({
    req(seuratObj(), seuratObj()@reductions$pca)
    
    
    pca_stdev <- seuratObj()@reductions$pca@stdev
    variance <- pca_stdev^2
    total_variance <- sum(variance)
    cumulative_var <- cumsum(variance) / total_variance
    percent_per_pc <- variance / total_variance * 100
    
    # 1 Variance drop between PCs ≥0.5%
    variance_diff <- diff(percent_per_pc)
    drop_point <- tail(which(variance_diff < -0.5), 1) + 1
    if (is.na(drop_point)) drop_point <- length(pca_stdev)
    
    # 2 Elbow point (using 2nd derivative)
    find_elbow_point <- function(variance) {
      first_diff <- diff(variance)
      second_diff <- diff(first_diff)
      elbow_point <- which.min(second_diff) + 1
      return(min(elbow_point, length(variance)))
    }
    elbow_point <- find_elbow_point(variance) + 1
    
    # 3 Cumulative variance ≥90%
    cum_point <- which(cumulative_var >= 0.9)[1]
    if (is.na(cum_point)) cum_point <- length(pca_stdev)
    
    # 4 Final recommended PCs
    recommended_point <- min(max(drop_point, elbow_point), cum_point)
    recommended_point <- min(recommended_point, length(pca_stdev))
    
    list(
      recommended = recommended_point,
      elbow_point = elbow_point,
      drop_point = drop_point,
      cum_point = cum_point,
      cumulative_var = cumulative_var
    )
  })
  
  output$elbow_point <- renderText({
    rec <- recommended_dimensions()
    paste0(rec$elbow_point, " PCs")
  })
  
  output$drop_point <- renderText({
    rec <- recommended_dimensions()
    paste0(rec$drop_point, " PCs")
  })
  
  output$cum_point <- renderText({
    rec <- recommended_dimensions()
    paste0(rec$cum_point, " PCs")
  })
  
  output$recommended_dims <- renderText({
    rec <- recommended_dimensions()
    paste0(rec$recommended, " PCs")
  })
  
  
  # PCA ELBOW PLOT
  output$elbow_plot <- renderPlot({
    req(seuratObj(), seuratObj()@reductions$pca)
    
    rec <- recommended_dimensions()
    pca_stdev <- seuratObj()@reductions$pca@stdev
    variance <- pca_stdev^2
    total_variance <- sum(variance)
    
    plot_data <- data.frame(
      PC = 1:length(variance),
      Variance = variance,
      Percent = variance / total_variance * 100,
      CumulativePercent = rec$cumulative_var * 100
    )
    
    ggplot(plot_data, aes(x = PC, y = Percent)) +
      geom_line(color = "blue", size = 1) +
      geom_point(color = "blue", size = 1.5) +
      geom_vline(xintercept = rec$recommended, 
                 color = "red", linetype = "dashed", size = 1.5) +
      geom_point(data = plot_data[rec$recommended, ], 
                 color = "red", size = 5, shape = 18) +
      annotate("text", x = rec$recommended, y = max(plot_data$Percent) * 0.9,
               label = paste("Recommended:", rec$recommended, "PCs"),
               color = "red", hjust = -0.1, size = 5, fontface = "bold") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "red", size = 11)
      )
  })
  
  observeEvent(input$integration_harmony, {
    harmony_state(input$integration_harmony)
  })
  
  observeEvent(input$harmony_metadata, {
    harmony_meta_state(input$harmony_metadata)
  })
  
  # Harmony metadata UI
  output$harmony_metadata_ui <- renderUI({
    req(seuratObj())
    metadata_cols <- colnames(seuratObj()@meta.data)
    tagList( 
    div(
        style = "display: inline-block; vertical-align: middle; margin-left:20px;",
        h5("Select metadata column for Harmony integration (sample identification):", style = "font-weight: bold; display: inline-block; margin-right:8px;"),
        div(
          style = "margin-right:10px; display:inline-block;",
          selectInput(
            "harmony_metadata", NULL,
            choices = metadata_cols,
            selected = harmony_meta_state() %||% metadata_cols[1],
            width = "200px"
          )
        )), 
      div(
        style = "margin-left:40px;",
        textOutput("harmony_metadata_preview")
      ),
    br()
    )
  })
  
  
  
  output$harmony_metadata_preview <- renderText({
    req(input$harmony_metadata, seuratObj())
    vals <- seuratObj()@meta.data[[input$harmony_metadata]]
    unique_vals <- unique(vals)
    preview_vals <- unique_vals[1:min(10, length(unique_vals))]
    result <- paste(preview_vals, collapse = ", ")
    if(length(unique_vals) > 10){
      result <- paste0(result, ", ...")
    }
    return(paste("Preview values:", result))
  })
  
  # NEIGHBOR STEP SERVER
  observeEvent(input$neighbor_btn_ui,{
    req(seuratObj())
    showModal(modalDialog(title = "Please wait", 
                          tagList(
                            p("Finding Neighbors..."),
                            p("The waiting time is related to the size of the dataset and number of PCs.")
                          ), footer = NULL, easyClose = FALSE))
    tryCatch({
      if(input$integration_harmony){
        req(input$harmony_metadata)
        seuratObj <- RunHarmony(seuratObj(), group.by.vars = input$harmony_metadata, reduction.use = "pca",
                                dims.use = 1:input$npc, assay.use = "RNA", reduction.save = "harmony", 
                                project.dim = TRUE, verbose = FALSE)
        seuratObj(FindNeighbors(seuratObj, reduction = "harmony", dims = 1:input$nn_dims, verbose = FALSE))
      } else {
        seuratObj(FindNeighbors(seuratObj(), reduction = "pca", dims = 1:input$nn_dims, verbose = FALSE))
      }
      updateCheckboxInput(session, "integration_harmony",
                          value = harmony_state())
      updateSelectInput(session, "harmony_metadata",
                        selected = harmony_meta_state())
      status$normalized <- TRUE
      status$scaled <- TRUE
      status$pca <- TRUE
      status$neighbor <- TRUE
      status$clustered <- FALSE
      removeModal()
      showNotification("Neighbors found.", type = "message")
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  
  # CLUSTERING STEP SERVER
  observeEvent(input$clustering_btn_ui,{
    req(seuratObj())
    showModal(modalDialog(title = "Please wait", 
                          tagList(
                            p("Finding Clusters..."),
                            p("The waiting time is related to the size of the dataset and resolution.")
                          ), footer = NULL, easyClose = FALSE))
    tryCatch({
      sobj <- seuratObj()
      sobj <- FindClusters(sobj, resolution = input$resolution, verbose = FALSE)
      
      # Annotate clusters if gene markers provided
      FindCellTypesByMarkers <- function(sobj, biomarkers = NULL) {
        
        if(is.null(biomarkers)) {
          stop("parameter biomarkers should be specied (now NULL)")
        }
        
        DefaultAssay(sobj) <- "RNA"
        
        tmp <- map_dfr(names(biomarkers), function(name) {
          mm <- FindAllMarkers(sobj,features = biomarkers[[name]] ,logfc.threshold = -Inf, min.pct = 0.0) 
          
          if (is_empty(mm)){
            mm <- tibble(pval = 1, 
                         avg_log2FC = 0, 
                         pct.1 = 0, 
                         pct.2 = 0, 
                         p_val_adj = 1,
                         cluster = as.factor(0), 
                         gene = biomarkers[[name]][1])
          }
          mm %>% 
            group_by(cluster) %>% 
            summarise(score = mean(avg_log2FC)) %>%
            mutate(celltype = name)
        }) 
        
        tmp <- left_join(tmp %>% tidyr::expand(cluster, celltype), tmp) %>%
          mutate(score = coalesce(score, 0))
        
        labels <- tmp %>% 
          group_by(cluster) %>%
          summarise(prediction = ifelse(score[which.max(score)]> 0, celltype[which.max(score)], "Unknown" )) %>% 
          ungroup() %>% 
          mutate(prediction = paste0(prediction,"(",cluster,")"))
        
        new_labels <- labels$prediction
        names(new_labels) <- labels$cluster
        
        ## swap names and values
        old_labels <- setNames(names(new_labels), new_labels)
        
        list(new_labels = new_labels, old_labels = old_labels, heatmap_table = tmp)   
      }
      
      # Convert Excel to Marker List
      
      convert_marker_excel_to_list <- function(path, sheet_name, allowed_genes = NULL) {
        tryCatch({
          df <- read.xlsx(path, sheet = sheet_name, colNames = TRUE)
          
          if (is.null(df) || ncol(df) == 0) {
            stop(paste0("The '", sheet_name, "' sheet appears empty or has no valid columns."))
          }
          
          if (any(is.na(names(df))) || any(names(df) == "")) {
            stop(paste0("Some columns in '", sheet_name, "' have no names. Please check header row."))
          }
          
          biomarkers <- list()
          
          for (col in names(df)) {
            genes <- df[[col]]
            genes <- genes[!is.na(genes) & genes != ""]
            genes <- as.character(genes)
            
            if (!is.null(allowed_genes)) {
              genes <- genes[genes %in% allowed_genes]
            }
            
            if (length(genes) > 0) {
              biomarkers[[col]] <- genes
            }
          }
          
          if (length(biomarkers) == 0) {
            stop(paste0("No valid markers found in '", sheet_name, "' sheet."))
          }
          
          return(biomarkers)
          
        }, error = function(e) {
          stop(paste("❌ Failed to convert Excel to marker list. Please check Excel format.\nDetails:", e$message))
        })
      }
      
      # Create Heatmap Function
      CreateCellTypesHeatmap <- function(df, labels_text_size = 6, xaxis_text_size = 12, yaxis_text_size = 12, rotate_x = FALSE){
        
        df$cluster <- factor(df$cluster, levels = sort(unique(df$cluster),TRUE), ordered = TRUE)
        
        tmp <- ggplot(df, aes(x = celltype, y = as.factor(cluster))) +
          geom_tile(aes(fill = score),color= "gray50",size = 0.1)+
          scale_fill_gradient2(low = "blue", mid="white", high = "tomato")+
          geom_text(aes(label=round(score,1)), size = labels_text_size)+
          scale_x_discrete(position = "top") +
          ylab("Clusters")+
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position = "bottom",
                axis.text.x = element_text(color = "black", size = xaxis_text_size),
                axis.text.y = element_text(color = "black", size = yaxis_text_size)) +
          coord_fixed()
        
        if (rotate_x){
          # tmp <- tmp+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color = "black", size = xaxis_text_size))
          tmp <- tmp+theme(axis.text.x = element_text(angle = 90, hjust=0, color = "black", size = xaxis_text_size),
                           axis.text.x.top = element_text(vjust = 0.5))
        }
        tmp
      }
      
      custom_dotplot <- function(seurat_obj, gene_list, assay = "RNA", title = "", cols = c("lightgrey", "blue")) {
        
        all_genes <- unlist(gene_list, use.names = FALSE)
        DefaultAssay(seurat_obj) <- assay
        clusters <- unique(Idents(seurat_obj))
        results <- data.frame()
        for(cluster in clusters) {
          cluster_cells <- WhichCells(seurat_obj, idents = cluster)
          cluster_data <- GetAssayData(seurat_obj, assay = assay, slot = "data")[all_genes, cluster_cells, drop = FALSE]
          
          for(gene in all_genes) {
            if(gene %in% rownames(cluster_data)) {
              avg_exp <- mean(expm1(cluster_data[gene, ]))
              pct_exp <- sum(cluster_data[gene, ] > 0) / length(cluster_cells) * 100
              results <- rbind(results, data.frame(
                gene = gene,
                cluster = cluster,
                avg_exp = avg_exp,
                pct_exp = pct_exp
              ))
            }
          }
        }
        
        results <- results %>%
          group_by(gene) %>%
          mutate(avg_exp_scaled = scale(avg_exp)[,1]) %>%
          ungroup()
        
        results$cluster <- factor(results$cluster, levels = sort(unique(results$cluster),TRUE), ordered = TRUE)
        
        gene_groups <- data.frame()
        for(group_name in names(gene_list)) {
          group_genes <- gene_list[[group_name]]
          gene_groups <- rbind(gene_groups, 
                               data.frame(gene = group_genes, group = group_name))
        }
        
        gene_groups$group <- factor(gene_groups$group, levels = names(gene_list))
        
        results <- results %>%
          left_join(gene_groups, by = "gene")
        
        results$gene <- factor(results$gene, levels = all_genes)
        
        p <- ggplot(results, aes(x = gene, y = cluster)) +
          geom_point(aes(size = pct_exp, color = avg_exp_scaled)) +
          scale_color_gradientn(
            colors = c("lightgrey", "blue"),
            name = "Avg Exp",
            guide = guide_colorbar(barwidth = 5, barheight = 0.5)
          ) +
          scale_size(
            range = c(0, 6),  
            name = "Per Exp",
            limits = c(0, 100),
            breaks = c(0, 25, 50, 75, 100)
          ) +
          facet_grid(. ~ group, scales = "free_x", space = "free_x") +
          theme_cowplot() +
          theme(
            axis.text.x = element_text(
              angle = 90, 
              hjust = 1, 
              vjust = 0.5,
              size = 8
            ),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.direction = "horizontal", 
            legend.position = "bottom",
            legend.box = "horizontal", 
            legend.justification = "center",
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9),
            panel.grid.major = element_blank(),  
            panel.grid.minor = element_blank(),  
            panel.spacing.y = unit(0.1, "lines"),  
            panel.spacing.x = unit(1, "lines"),   
            strip.text = element_text(face = "bold", size = 9),
            strip.background = element_rect(fill = "grey95", color = NA)
          ) +
          labs(x = "Genes", y = "Clusters")
        
        return(p)
      }
      
      # Read Gene Markers
        path <- if (!is.null(input$upload_excel)) {
          input$upload_excel$datapath
        } else {
          "/projectnb/wax-es/00_shinyapp/Clustering/Clustering_shinyapp/Gene_Markers.xlsx"
        }
        
        allowed_genes <- rownames(sobj)
        tryCatch({
          biomarkers_heatmap <- convert_marker_excel_to_list(path, sheet_name = "HEATMAP", allowed_genes)
        }, error = function(e) {
          showModal(modalDialog(
            title = "Error",
            paste("Failed to read gene marker file:", e$message),
            easyClose = TRUE
          ))
        })
        
        tryCatch({
          biomarkers_dotplot1 <- convert_marker_excel_to_list(path, sheet_name = "DOTPLOT1", allowed_genes)
        }, error = function(e) {
          biomarkers_dotplot1 <- NULL
        })
        
        tryCatch({
          biomarkers_dotplot2 <- convert_marker_excel_to_list(path, sheet_name = "DOTPLOT2", allowed_genes)
        }, error = function(e) {
          biomarkers_dotplot1 <- NULL
        })
      
        if (!is.null(biomarkers_dotplot1)) {
          dotplot1_plot(custom_dotplot(sobj, biomarkers_dotplot1))
        }
        
        if (!is.null(biomarkers_dotplot2)) {
          dotplot2_plot(custom_dotplot(sobj, biomarkers_dotplot2))
        }
        
      # Annotate clusters
        findcelltypes_result <- FindCellTypesByMarkers(sobj, biomarkers_heatmap)
      heatmap <- CreateCellTypesHeatmap(findcelltypes_result$heatmap_table,
                                 labels_text_size = 4,
                                 xaxis_text_size = 10,
                                 yaxis_text_size = 10,
                                 rotate_x = TRUE)
      heatmap_plot(heatmap)
      
      sobj <- RenameIdents(sobj, findcelltypes_result$new_labels)
      sobj$shiny_clusters <- Idents(sobj)
      Idents(sobj) <- sobj$shiny_clusters
      if(input$integration_harmony){sobj <- RunUMAP(
        sobj,
        dims = 1:input$nn_dims,
        min.dist = 0.3,
        reduction = "harmony",
        verbose = FALSE
      )}else{
        sobj <- RunUMAP(
          sobj,
          dims = 1:input$nn_dims,
          min.dist = 0.3,
          reduction = "pca",
          verbose = FALSE
        )
      }
      if(harmony_state()){
        plot_title(paste0("resolution: ", input$resolution, ", min_dist: ", 0.3, ", PCA:", input$nn_dims, ", Harmony"))
      } else {
        plot_title(paste0("resolution: ", input$resolution, ", min_dist: ", 0.3, ", PCA:", input$nn_dims))
      }
      umap_plot(DimPlot(sobj, reduction = "umap", label = FALSE) + 
                  ggtitle("") +
                  theme(plot.title = element_text(face = "bold", size = 10)))
      seuratObj(sobj)
      status$normalized <- TRUE
      status$scaled <- TRUE
      status$pca <- TRUE
      status$neighbor <- TRUE
      status$clustered <- TRUE
      removeModal()
      showNotification("Clustering completed.", type = "message")
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  
  output$clustered_status <- reactive({
    isTRUE(status$clustered)
  })
  outputOptions(output, "clustered_status", suspendWhenHidden = FALSE)
  
  
  # UMAP PLOT UI
  output$umap_ui <- renderUI(
    tagList(
      hr(),
      fluidRow(
        div(style = "display: inline-block; vertical-align: top; margin-left:20px;",
            h5("UMAP Minimum Distance:",style = "font-weight: bold;")),
        div(style = "display: inline-block; vertical-align: top; margin-left:10px;",
            uiOutput("umap_dist_ui")),
        div(class = "tooltip-circle",
            "?",  
            span(class = "tooltip-text", 
                 div(class = "tooltip-title", "Minimum Distance"),
                 div(class = "tooltip-divider"),
                 div("Minimum distance parameter controls how tightly UMAP packs points together.
                     Lower values result in a more clustered embedding, while higher values
                     produce a more spread-out representation. A typical starting value is 0.3.
                     Adjusting this parameter can help reveal different structures in the data.")
            )
        ),
        div(style = "display: inline-block; vertical-align: top; margin-left:30px;",
            actionButton("update_umap_btn_ui", "Update UMAP", width = "120px", class = "btn-primary"))
      ),
      fluidRow(
        column(12,
               uiOutput("cluster_plot_ui")))
    )
  )
  
  output$umap_dist_ui <- renderUI({
    numericInput("umap_dist", NULL, 
                 value = isolate(input$umap_dist) %||% 0.3, 
                 min = 0.01, max = 0.99, step = 0.01, width = "100%")
  })
  
  output$download_template <- downloadHandler(
    filename = function() {
      "Gene_Markers.xlsx"
    },
    content = function(file) {
      template_path <- "/projectnb/wax-es/00_shinyapp/Clustering/Clustering_shinyapp/Gene_Markers.xlsx" 
      
      if (file.exists(template_path)) {
        file.copy(template_path, file)
      }
    }
  )
  
  # upload excel
  observeEvent(input$upload_excel, {
    req(input$upload_excel)
    
    ext <- tools::file_ext(input$upload_excel$name)
    validate(
      need(ext %in% c("xlsx", "xls"), "Please upload .xlsx or .xls file")
    )
    
    tryCatch({
      
      sheets <- openxlsx::getSheetNames(input$upload_excel$datapath)
      seurat_genes <- rownames(seuratObj())
      
      missing_list <- list()
      
      for (sheet in sheets) {
        df <- openxlsx::read.xlsx(input$upload_excel$datapath, sheet = sheet)
        
        for (colname in colnames(df)) {
          col_values <- df[[colname]] |>
            as.character() |>
            unique() |>
            setdiff(c(NA, "", " "))
          
          if (length(col_values) == 0) next
          
          missing <- setdiff(col_values, seurat_genes)
          
          if (length(missing) > 0) {
            missing_list[[paste0(sheet, "_", colname)]] <- list(
              sheet = sheet,
              column = colname,
              missing_markers = missing
            )
          }
        }
      }
      
      sheets_to_check <- intersect(c("DOTPLOT1", "DOTPLOT2"), sheets)
      
      for (sheet in sheets_to_check) {
        df <- openxlsx::read.xlsx(input$upload_excel$datapath, sheet = sheet)
        
        markers_by_group <- list()
        for (colname in colnames(df)) {
          genes <- df[[colname]] |>
            as.character() |>
            unique() |>
            setdiff(c(NA, "", " "))
          if (length(genes) > 0) {
            markers_by_group[[colname]] <- genes
          }
        }
        
        gene_groups_map <- list()
        for (group in names(markers_by_group)) {
          for (gene in markers_by_group[[group]]) {
            if (!is.null(gene_groups_map[[gene]])) {
              gene_groups_map[[gene]] <- c(gene_groups_map[[gene]], group)
            } else {
              gene_groups_map[[gene]] <- group
            }
          }
        }
        
        duplicates <- Filter(function(g) length(g) > 1, gene_groups_map)
        
        if (length(duplicates) > 0) {
          dup_msg <- paste(
            sapply(names(duplicates), function(gene) {
              paste0(
                "Marker '", gene, "' appears in columns:<br>",
                paste(duplicates[[gene]], collapse = ", ")
              )
            }), collapse = "<br><br>"
          )
          
          showModal(modalDialog(
            title = paste0("🚫 Error: Duplicate markers detected in ", sheet, " sheet"),
            HTML(
              paste0(
                "The following markers appear in multiple columns within ", sheet, 
                ". Please remove duplicates before proceeding:<br>", dup_msg
              )
            ),
            easyClose = FALSE,
            footer = modalButton("Close")
          ))
          
          file.remove(input$upload_excel$datapath)
          return(NULL)
        }
      }
      
      if (length(missing_list) > 0) {
        msg <- ""
        for (item in missing_list) {
          msg <- paste0(
            msg,
            "<b>Sheet:</b> ", item$sheet, "<br>",
            "<b>Column:</b> ", item$column, "<br>",
            "<b>Missing markers:</b> ", paste(item$missing_markers, collapse = ", "), 
            "<br><br>"
          )
        }
        
        showModal(modalDialog(
          title = "⚠️ Warning: Some markers are not found in Seurat object",
          HTML(paste0(
            msg,
            "<i>If you proceed, these markers will be automatically ignored in downstream analysis.</i>"
          )),
          easyClose = TRUE
        ))
      }
      
      showNotification("File uploaded successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      file.remove(input$upload_excel$datapath)
    })
    
  })
  
  
  # SEARCH RESOLUTION BUTTON
  observeEvent(input$search_resolution_btn_ui, {
    req(seuratObj(), status$neighbor, input$start_resolution, input$end_resolution, input$step_resolution)
    
    seq_res <- seq(input$start_resolution, input$end_resolution, by = input$step_resolution)
    cluster_numbers <- numeric(length(seq_res))
    
    showModal(modalDialog(
      title = "Please wait",
      tagList(
        p("Searching resolutions..."),
        p("The waiting time depends on dataset size and resolution range.")
      ),
      footer = NULL, easyClose = FALSE
    ))
    
    tryCatch({
      for (i in seq_along(seq_res)) {
        res <- seq_res[i]
        tmp_obj <- FindClusters(seuratObj(), resolution = res, verbose = FALSE)
        num_clusters <- length(unique(Idents(tmp_obj)))
        cluster_numbers[i] <- num_clusters
      }
      
      result_df <- data.frame(
        Resolution = seq_res,
        Number_of_Clusters = cluster_numbers
      )
      resolution_search_results(result_df)
      removeModal()
      
      showModal(modalDialog(
        title = "Resolution Search Results",
        tagList(
          p("Search completed! You can click 'Show search result' later to view again."),
          renderTable(result_df, align = 'c', bordered = TRUE, striped = TRUE, hover = TRUE,
                      width = 'auto', rownames = FALSE)
        ),
        easyClose = TRUE
      ))
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(
        title = "Resolution search failed",
        p(e$message),
        easyClose = TRUE
      ))
    })
  })
  
  output$show_search_result_ui <- renderUI({
    req(resolution_search_results())
    div(style = "display: inline-block; vertical-align: top; margin-left:10px;",
        actionButton("show_search_result_btn", "Show search result", 
                     class = "btn-info", width = "150px"))
  })
  
  observeEvent(input$show_search_result_btn, {
    req(resolution_search_results())
    showModal(modalDialog(
      title = "Resolution Search Results",
      renderTable(resolution_search_results(), align = 'c', bordered = TRUE, striped = TRUE, hover = TRUE,
                  width = 'auto', rownames = FALSE),
      easyClose = TRUE
    ))
  })
  
  
  observeEvent(input$update_umap_btn_ui, {
    req(seuratObj(), status$clustered, input$umap_dist)
    
    if (input$umap_dist >= 0.01 && input$umap_dist <= 1.0) {
      showModal(modalDialog("Updating UMAP with new min.dist...", footer = NULL, easyClose = FALSE))
      
      tryCatch({
        if(harmony_state()){
          plot_title(paste0("resolution: ", input$resolution, ", min_dist: ", input$umap_dist %||% 0.3, ", PCA:", input$nn_dims, ", Harmony"))
        } else {
          plot_title(paste0("resolution: ", input$resolution, ", min_dist: ", input$umap_dist %||% 0.3, ", PCA:", input$nn_dims))
        }
        seuratObj(RunUMAP(seuratObj(), 
                          dims = 1:input$nn_dims, 
                          min.dist = input$umap_dist, 
                          verbose = FALSE))
        umap_plot(DimPlot(seuratObj(), reduction = "umap", label = FALSE) + 
                    ggtitle("") +
                    theme(plot.title = element_text(face = "bold", size = 10)))
        removeModal()
        showNotification("UMAP updated with new min.dist!", type = "message")
      }, error = function(e) {
        removeModal()
        showModal(modalDialog("UMAP update failed:", e$message, easyClose = TRUE))
      })
    }
  })
  
  output$cluster_umap_ui <- renderPlot({
    req(!is.null(umap_plot()))
    umap_plot()
  })
  
  output$cluster_table_ui <- renderTable({
    req(seuratObj())
    get_cluster_summary_table <- function(seurat_obj) {
      all_cells <- seurat_obj@meta.data %>%
        select(seurat_clusters, counts = nCount_RNA, genes = nFeature_RNA) %>%
        group_by(seurat_clusters) %>%
        summarise(ncells = n(),
                  avg.counts = as.integer(round(mean(counts))), 
                  avg.genes = as.integer(round(mean(genes)))) %>%  
        ungroup()
      
      result <- all_cells %>%
        mutate(all = sum(ncells),
               cluster = as.character(seurat_clusters),
               pct = round(100 * ncells / all, 2)) %>%
        select(cluster, ncells, pct, avg.counts, avg.genes)
      
      result <- bind_rows(
        result,
        result %>%
          summarise(cluster = "total",
                    ncells = sum(ncells),
                    pct = sum(pct),
                    avg.counts = as.integer(round(mean(seuratObj()@meta.data$nCount_RNA))),
                    avg.genes = as.integer(round(mean(seuratObj()@meta.data$nFeature_RNA))))
      )
      
      return(result)
    }
    cluster_table(get_cluster_summary_table(seuratObj()))
    cluster_table()
  }, 
  align = 'c',
  bordered = FALSE,
  striped = FALSE,
  hover = FALSE,
  width = 'auto',
  rownames = FALSE,
  colnames = TRUE
  )
  
  output$cluster_heatmap_ui <- renderPlot({
    req(!is.null(heatmap_plot()))
    heatmap_plot()
  })
  
  output$cluster_dotplot1_ui <- renderPlot({
    req(!is.null(dotplot1_plot()))
    dotplot1_plot()
  })
  
  output$cluster_dotplot2_ui <- renderPlot({
    req(!is.null(dotplot2_plot()))
    dotplot2_plot()
  })
  
  output$plot_title <- renderText({
    req(plot_title())
    plot_title()
  })
  
  output$cluster_plot_ui <- renderUI({
    tagList(
      div(style = "display: inline-block; vertical-align: top; margin-left:30px; font-weight: bold;",fluidRow(textOutput("plot_title"))),
      fluidRow(
        column(4,
               plotOutput("cluster_umap_ui")),
        column(3,
               tags$style(HTML("
                            table {
                              border: none !important;
                              font-size: 12px !important;
                            }
                            table th, table td {
                              border: none !important;
                              font-size: 12px !important;
                            }
                            table tr {
                              border: none !important;
                            }
                          ")),
               tableOutput("cluster_table_ui")
        ),
        column(5,
               plotOutput("cluster_heatmap_ui"))
      ),
      hr(),
      fluidRow( 
        column(12, plotOutput("cluster_dotplot1_ui"))
      ),
      fluidRow( 
        column(12, plotOutput("cluster_dotplot2_ui"))
      ),
      fluidRow( 
        column(12, textInput("title_for_download", "Title for Downloaded Plots:", 
                             value = plot_title() %||% "Cluster Plots", width = "800px"))
      ),
      fluidRow(
        div(style = "display: inline-block; vertical-align: top; margin-left:10px;", 
            downloadButton("download_cluster_plots", "Download Cluster Plots", 
            style = "background-color: #4CAF50; color: white;")),
        div(class = "tooltip-circle",
            "?",  
            span(class = "tooltip-text", 
                 div(class = "tooltip-title", "DOWNLOAD VISUALIZATION RESULT"),
                 div(class = "tooltip-divider"),
                 div("Edit \"Title for Downloaded Plots:\" to change the title of visualization file. 
                     Click the 'Download Cluster Plots' button to download a ZIP file containing the UMAP plot, 
                     heatmap, dot plots, and a combined plot with the specified title.")
            )
        )
    ))
  })  
  
  # DOWNLOAD CLUSTER PLOTS
  output$download_cluster_plots <- downloadHandler(
    filename = function() {
      paste("Cluster_Plots_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      temp_dir <- tempdir()
      umap_file <- file.path(temp_dir, "Cluster_UMAP.pdf")
      heatmap_file <- file.path(temp_dir, "Cluster_Heatmap.pdf")
      dotplot1_file <- file.path(temp_dir, "Cluster_Dotplot1.pdf")
      dotplot2_file <- file.path(temp_dir, "Cluster_Dotplot2.pdf")
      combined_file <- file.path(temp_dir, "Cluster_Combined.pdf")
      
      create_combined_cluster_plot_patchwork <- function(umap_plot, table_data, heatmap_plot, dotplot1 = NULL, dotplot2 = NULL) {
        
        table_plot <- function(table_data) {
          if (is.null(table_data)) return(ggplot() + theme_void())
          
          table_grob <- gridExtra::tableGrob(
            table_data,
            rows = NULL,
            theme = gridExtra::ttheme_minimal(
              base_size = 8,
              padding = unit(c(2, 2), "mm")
            )
          )
          
          ggplot() + 
            annotation_custom(table_grob) + 
            theme_void() +
            theme(plot.margin = margin(5, 5, 5, 5))
        }
        
        title_text <- input$title_for_download %||% "Cluster Plots"
        top_row <- (umap_plot | table_plot(table_data) | heatmap_plot) + 
          plot_layout(widths = c(7, 4, 7))
        title_plot <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, 
                   label = title_text, 
                   size = 10,
                   hjust = 0.5, vjust = 0.5, 
                   fontface = "bold") +
          theme_void() +
          theme(
            plot.margin = margin(5, 5, 5, 5),
            plot.background = element_rect(fill = "white", color = NA)
          )
        combined_plot <- title_plot / top_row
        
        bottom_plots <- list()
        
        if (!is.null(dotplot1)) {
          bottom_plots[[length(bottom_plots) + 1]] <- dotplot1
        }
        if (!is.null(dotplot2)) {
          bottom_plots[[length(bottom_plots) + 1]] <- dotplot2
        }
        
        if (length(bottom_plots) > 0) {
          for (plot in bottom_plots) {
            combined_plot <- combined_plot / plot
          }
        }
        
        n_parts <- 1 + 1 + length(bottom_plots)  
        heights <- c(0.2, 2, rep(1.5, length(bottom_plots)))
        
        combined_plot <- combined_plot + 
          plot_layout(heights = heights[1:n_parts], ncol = 1)
        
        return(combined_plot)
      }
      
      combined_plot <- create_combined_cluster_plot_patchwork(
        umap_plot(),
        cluster_table(),
        heatmap_plot(),
        dotplot1_plot(),
        dotplot2_plot()
      )
      
      ggsave(heatmap_file, plot = heatmap_plot(), width = 7, height = 8)
      ggsave(umap_file, plot = umap_plot(), width = 7, height = 8)
      ggsave(combined_file, plot = combined_plot, width = 18, height = 20)
      
      if (!is.null(dotplot1_plot())) {
        ggsave(dotplot1_file, plot = dotplot1_plot(), width = 18, height = 6)
      }
      
      if (!is.null(dotplot2_plot())) {
        ggsave(dotplot2_file, plot = dotplot2_plot(), width = 18, height = 6)
      }
      
      zip_files <- c(heatmap_file)
      if (file.exists(dotplot1_file)) {
        zip_files <- c(zip_files, dotplot1_file)
      }
      if (file.exists(dotplot2_file)) {
        zip_files <- c(zip_files, dotplot2_file)
      }
      zip_files <- c(zip_files, umap_file, combined_file)
      zip::zip(zipfile = file, files = zip_files, mode = "cherry-pick")
    }
  )
  
  output$sample_subset_rename_ui <- renderUI({
    req(seuratObj())
    vals <- unique(Idents(seuratObj()))
    vals <- vals[order(as.numeric(gsub(".*\\((\\d+)\\).*", "\\1", vals)))]
    
    tagList(
      lapply(seq_along(vals), function(i) {
        val <- vals[i]
        fluidRow(
          column(
            4,
            checkboxInput(
              paste0("sample_cb_", val),
              label = paste0(i - 1, ": ", val),
              value = TRUE
            )
          ),
          column(
            8,
            textInput(
              paste0("sample_rename_", val),
              label = NULL,
              value = val,
              placeholder = "New name"
            )
          )
        )
      })
    )
  })
  
  
  observeEvent(input$rename_btn_ui,{
    req(seuratObj())
    seurat_object <- seuratObj()
    get_selected_values <- function(obj, col, prefix) {
      vals <- as.character(unique(obj@meta.data[[col]]))
      selected <- sapply(vals, function(val) {
        if(isTRUE(input[[paste0(prefix, "_", val)]])) val else NA
      })
      na.omit(selected)
    }
    
    apply_renaming <- function(obj, vals, prefix) {
      ident_vec <- as.character(Idents(obj))
      for (val in vals) {
        new_name <- input[[paste0(prefix, "_", val)]]
        if (!is.null(new_name) && nzchar(new_name)) {
          ident_vec[ident_vec == val] <- new_name
        }
      }
      Idents(obj) <- ident_vec
      return(obj)
    }
    tryCatch({
      # Get selected values with validation
      sample_vals <- get_selected_values(seuratObj(), current_ident(), "sample_cb")
      sample_table <- data.frame(
        `Selected Sample` = sample_vals,
        `Renamed Value` = sapply(sample_vals, function(val) {
          new_name <- input[[paste0("sample_rename_", val)]]
          ifelse(!is.null(new_name) && nzchar(new_name), new_name, val)
        }),
        check.names = FALSE
      )
      
      seurat_object <- apply_renaming(seurat_object, sample_vals, "sample_rename")
      
      # Update reactive values
      seuratObj(seurat_object)
      showNotification("Renaming completed successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  observeEvent(input$subset_btn_ui, {
    showModal(
      modalDialog(
        title = "Confirm Subset",
        tagList(
          h4("Are you sure you want to subset the data?"),
          br(),
          p("⚠️ Please make sure you have saved the data before subsetting."),
          p("After subsetting, some data may no longer be recoverable.")
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_subset", "Confirm", class = "btn-danger")
        ),
        easyClose = FALSE
      )
    )
  })
  
  observeEvent(input$keep_analysis, {
    showModal(
      modalDialog(
        title = "Confirm analysising on subset data",
        tagList(
          h4("Are you sure you want to keep analysis on subset data?"),
          br(),
          p("⚠️ Please ensure that you save the data after subsetting before analyzing the subset data.."),
          p("Click this buttom, some data may no longer be recoverable.")
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_keep_analysis", "Confirm", class = "btn-danger")
        ),
        easyClose = FALSE
      )
    )
  })
  
  observeEvent(input$confirm_subset, {
    req(seuratObj())
    showModal(
      modalDialog(
        title = "Please wait",
        tagList(
          p("Subsetting data..."),
          p("The waiting time depends on the size of the dataset and number of selected samples.")
        ),
        footer = NULL, easyClose = FALSE
      )
    )
    seuratobj <- seuratObj()
    tryCatch({
      get_selected_values <- function(obj, col, prefix) {
        vals <- as.character(unique(obj@meta.data[[col]]))
        selected <- sapply(vals, function(val) {
          if(isTRUE(input[[paste0(prefix, "_", val)]])) val else NA
        })
        na.omit(selected)
      }
      sample_vals <- get_selected_values(seuratObj(), current_ident(), "sample_cb")
      
      validate(
        need(length(sample_vals) > 0, "Please select at least one sample")
      )
      keep_cells <- names(Idents(seuratobj))[(Idents(seuratobj) %in% sample_vals) ]
      validate(
        need(length(keep_cells) > 0, "No cells match selected samples")
      )
      seuratobj <- subset(seuratobj, cells = keep_cells)
      seuratObj(seuratobj)
      
      showNotification("Subsetting completed successfully!", type = "message")
      removeModal()
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  observeEvent(input$confirm_keep_analysis, {
    req(seuratObj())
    tryCatch({
      
      subset_data_path(input$data_path)
      
      status$normalized <- FALSE
      status$scaled <- FALSE
      status$pca <- FALSE
      status$neighbor <- FALSE
      status$clustered <- FALSE
      variable_method(NULL)
      resolution_search_results(NULL)
      umap_plot(NULL)
      heatmap_plot(NULL)
      cluster_table(NULL)
      dotplot1_plot(NULL)
      dotplot2_plot(NULL)
      current_ident(NULL)
      subset_data_path(NULL)
      plot_title(NULL)
      harmony_state(NULL)
      harmony_meta_state(NULL)
      shinyjs::runjs("location.reload();")
      removeModal()
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
}  

shinyApp(ui = ui, server = server)
