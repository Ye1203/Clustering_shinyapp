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

ui <- fluidPage(
  titlePanel("Clustering"),
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
        value = if (is.null(input$data_path)) {
          "/projectnb/wax-es/00_shinyapp/Clustering/forBingtian/example_data.rds"
        } else {
          input$data_path
        },
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
  
  
  # READ  or RESET BUTTON UI
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
          textInput("save_path", "Input save path on SCC or local file name", value = paste0(file_path_sans_ext(input$data_path), "_", Sys.Date(), ".rds"), width = "100%"))),
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
        if (isTRUE(status$normalized)) {
          tags$p("✅ DATA NORMALIZED", style = "color:black; font-weight:bold; font-size:14px;")
        } else {
          tags$p("🔘 DATA NORMALIZED", style = "color:black; font-weight:bold; font-size:14px;")
        },
        uiOutput("normalized_data_ui")
      ))
    }
    
    ui_list <- append(ui_list, list(hr()))
    
    # --- Scaling ---
    if (stop_flag || !"data" %in% layers) {
      ui_list <- append(ui_list, list(
        tags$p("🔴 DATA SCALED", style = "color:red; font-weight:bold; font-size:14px;")
      ))
      stop_flag <- TRUE
    } else {
      ui_list <- append(ui_list, list(
        if (isTRUE(status$scaled)) {
          tags$p("✅ DATA SCALED", style = "color:black; font-weight:bold; font-size:14px;")
        } else {
          tags$p("🔘 DATA SCALED", style = "color:black; font-weight:bold; font-size:14px;")
        },
        uiOutput("scaled_data_ui")
      ))
    }
    
    ui_list <- append(ui_list, list(hr()))
    
    # --- PCA ---
    if (stop_flag || !"scale.data" %in% layers) {
      ui_list <- append(ui_list, list(
        tags$p("🔴 RUN PCA", style = "color:red; font-weight:bold; font-size:14px;")
      ))
      stop_flag <- TRUE
    } else {
      ui_list <- append(ui_list, list(
        if (isTRUE(status$pca)) {
          tags$p("✅ RUN PCA", style = "color:black; font-weight:bold; font-size:14px;")
        } else {
          tags$p("🔘 RUN PCA", style = "color:black; font-weight:bold; font-size:14px;")
        },
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
        if (isTRUE(status$neighbor)) {
          tags$p("✅ FIND NEIGHBORS", style = "color:black; font-weight:bold; font-size:14px;")
        } else {
          tags$p("🔘 FIND NEIGHBORS", style = "color:black; font-weight:bold; font-size:14px;")
        },
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
        if (isTRUE(status$clustered)) {
          tags$p("✅ FIND CLUSTERS", style = "color:black; font-weight:bold; font-size:14px;")
        } else {
          tags$p("🔘 FIND CLUSTERS", style = "color:black; font-weight:bold; font-size:14px;")
        },
        uiOutput("clustering_ui")
      ))
    }
    
    ui_list <- append(ui_list, list(hr()))
    
    # --- Rename and Subset ---
    if (stop_flag) {
      ui_list <- append(ui_list, list(
        tags$p("🔴 RENAME AND SUBSET", style = "color:red; font-weight:bold; font-size:14px;")
      ))
    } else {
      ui_list <- append(ui_list, list(
        tags$p("🔘 RENAME AND SUBSET", style = "color:black; font-weight:bold; font-size:14px;"),
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
      fluidRow(
        tags$div(
          style = "margin-left:20px;",
          radioButtons(
            "findvariable",
            "Gene Selection Method:",
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
      10
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
          )
        )
      },
      br(),br(),br(),
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
                 div(style = "display: inline-block; vertical-align: top; margin-left:20px;",
                     h5("Upload Gene Marker Excel file for cluster annotation (optional):", 
                        style = "font-weight: bold;"))
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
                   h5("Search resolutions in range for cluster numbers (Longer Processing Time!!!!!):",
                      style = "font-weight: bold; margin-left:0px;"),
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
                                    "Search",
                                    width = "100px", class = "btn-primary")
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
                     numericInput("resolution", NULL, 
                                  value = 0.03, min = 0.01, max = 5.0, step = 0.001, width = "100%")),
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
  
  
  # FINAL STEP UI
  output$final_ui <- renderUI({
    req(seuratObj(),input$resolution, status$clustered)
    if(is.null(current_ident())){
      val <- "shiny_clusters"
      current_ident(val)
    }
    tagList(
      fluidRow(
        div(style = "display: inline-block; vertical-align: top; margin-left:20px;",
            h5("Current ident name:",
               style = "font-weight: bold;")),
        div(style = "display: inline-block; vertical-align: top; margin-left:20px;",
            h5(current_ident(),
               style = "font-weight: bold; color: red;")),
        div(style = "display: inline-block; vertical-align: top; margin-left:10px;",
            textInput("rename_idents_input", NULL, value = current_ident(),
              width = "200px"
            )),
        div(style = "display: inline-block; vertical-align: top; margin-left:30px;",
            actionButton("rename_idents_btn","Rename Ident",
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
  
  # NEIGHBOR STEP SERVER
  observeEvent(input$neighbor_btn_ui,{
    req(seuratObj())
    showModal(modalDialog(title = "Please wait", 
                          tagList(
                            p("Finding Neighbors..."),
                            p("The waiting time is related to the size of the dataset and number of PCs.")
                          ), footer = NULL, easyClose = FALSE))
    tryCatch({
      seuratObj(FindNeighbors(seuratObj(), dims = 1:input$nn_dims, verbose = FALSE))
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
      seuratObj(FindClusters(seuratObj(), resolution = input$resolution, verbose = FALSE))
      
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
      convert_marker_excel_to_list <- function(path, sheet_name) {
        tryCatch({
          df <- read.xlsx(path, sheet = sheet_name, colNames = TRUE)
          
          if (is.null(df) || ncol(df) == 0) {
            stop(paste0("The \'", sheet_name,"\' sheet appears empty or has no valid columns."))
          }
          
          if (any(is.na(names(df))) || any(names(df) == "")) {
            stop(paste0("Some columns in \'", sheet_name,"\' have no names. Please check header row."))
          }
          biomarkers <- list()
          for (col in names(df)) {
            genes <- df[[col]]
            genes <- genes[!is.na(genes) & genes != ""]
            if (length(genes) > 0) {
              biomarkers[[col]] <- as.character(genes)
            }
          }
          
          if (length(biomarkers) == 0) {
            stop(paste0("No valid markers found in \'", sheet_name,"\' sheet."))
          }
          
          return(biomarkers)
        }, error = function(e) {
          stop(paste("❌ Failed to convert Excel to marker list. Please check Excel format.\nDetails:", e$message))
        })
      }
      
      # Create Heatmap Function
      CreateCellTypesHeatmap <- function(df, labels_text_size = 6, xaxis_text_size = 12, yaxis_text_size = 12, rotate_x = FALSE, celltype_order = NULL){
        
        if(!is.null(celltype_order)){
          df$celltype = factor(df$celltype, levels = celltype_order, ordered = TRUE)
        }
        
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
                axis.text.y = element_text(color = "black", size = yaxis_text_size))
        
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
            name = "Average Expression",
            guide = guide_colorbar(barwidth = 5, barheight = 0.5)
          ) +
          scale_size(
            range = c(0, 6),  
            name = "Percent Expressed",
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
        
        
        tryCatch({
          biomarkers_heatmap <- convert_marker_excel_to_list(path, sheet_name = "HEATMAP")
        }, error = function(e) {
          showModal(modalDialog(
            title = "Error",
            paste("Failed to read gene marker file:", e$message),
            easyClose = TRUE
          ))
        })
        
        tryCatch({
          biomarkers_dotplot1 <- convert_marker_excel_to_list(path, sheet_name = "DOTPLOT1")
        }, error = function(e) {
          biomarkers_dotplot1 <- NULL
        })
        
        tryCatch({
          biomarkers_dotplot2 <- convert_marker_excel_to_list(path, sheet_name = "DOTPLOT2")
        }, error = function(e) {
          biomarkers_dotplot1 <- NULL
        })
      
      # Annotate clusters
      findcelltypes_result <- FindCellTypesByMarkers(seuratObj(), biomarkers_heatmap)
      heatmap <- CreateCellTypesHeatmap(findcelltypes_result$heatmap_table,
                                 labels_text_size = 5,
                                 xaxis_text_size = 12,
                                 yaxis_text_size = 12,
                                 rotate_x = TRUE,
                                 celltype_order = names(biomarkers_heatmap))
      heatmap_plot(heatmap)
      
      if (!is.null(biomarkers_dotplot1)) {
        dotplot1 <- custom_dotplot(seuratObj(), biomarkers_dotplot1, title = NULL)
        dotplot1_plot(dotplot1)
      }
      
      if (!is.null(biomarkers_dotplot2)) {
        dotplot2 <- custom_dotplot(seuratObj(), biomarkers_dotplot2, title = NULL)
        dotplot2_plot(dotplot2)
      }
      seurat_object <- seuratObj()
      seurat_object <- RenameIdents(seurat_object, findcelltypes_result$new_labels)
      seurat_object <- AddMetaData(seurat_object, metadata = Idents(seurat_object), col.name = "shiny_clusters")
      Idents(seurat_object) <- seurat_object@meta.data$shiny_clusters
      seurat_object <- RunUMAP(seurat_object, 
                        dims = 1:input$nn_dims, 
                        min.dist = 0.3, 
                        verbose = FALSE)
      umap_plot(DimPlot(seurat_object, reduction = "umap", label = FALSE) + 
                  ggtitle(paste0("resolution: ", input$resolution, ", min_dist: ", 0.3, ", PCA:", input$nn_dims)) +
                  theme(plot.title = element_text(face = "bold", size = 10)))
      seuratObj(seurat_object)
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
          numericInput("umap_dist", NULL, value = input$umap_dist %||% 0.3, min = 0.01, max = 0.99, step = 0.01, width = "100%")),
        div(style = "display: inline-block; vertical-align: top; margin-left:30px;",
            actionButton("update_umap_btn_ui", "Update UMAP", width = "120px", class = "btn-primary"))
        ),
      fluidRow(
        column(12,
               uiOutput("cluster_plot_ui")))
    )
  )
  
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
    
    # Read Upload file
    tryCatch({
      library(readxl)
      uploaded_data <- read.xlsx(input$upload_excel$datapath)
      
      showNotification("File upload successfully！", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
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
        seuratObj(RunUMAP(seuratObj(), 
                          dims = 1:input$nn_dims, 
                          min.dist = input$umap_dist, 
                          verbose = FALSE))
        umap_plot(DimPlot(seuratObj(), reduction = "umap", label = FALSE) + 
                    ggtitle(paste0("resolution = ", input$resolution, ", min_dist = ", input$umap_dist %||% 0.3, ", PCA:", input$nn_dims)) +
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
  
  output$cluster_plot_ui <- renderUI({
    tagList(
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
        column(12, downloadButton("download_cluster_plots", "Download Cluster Plots", 
                                   style = "background-color: #4CAF50; color: white;"))
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
        
        top_row <- (umap_plot | table_plot(table_data) | heatmap_plot) + 
          plot_layout(widths = c(7, 4, 7))
        
        bottom_plots <- list()
        
        if (!is.null(dotplot1)) {
          bottom_plots[[length(bottom_plots) + 1]] <- dotplot1
        }
        if (!is.null(dotplot2)) {
          bottom_plots[[length(bottom_plots) + 1]] <- dotplot2
        }
        if (length(bottom_plots) == 0) {
          return(top_row)
        }
        combined_plot <- top_row
        for (plot in bottom_plots) {
          combined_plot <- combined_plot / plot
        }
        heights <- c(1, rep(1, length(bottom_plots)))
        combined_plot <- combined_plot + 
          plot_layout(heights = heights, ncol = 1)
        
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
      shinyjs::runjs("location.reload();")
      removeModal()
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
}  

shinyApp(ui = ui, server = server)
