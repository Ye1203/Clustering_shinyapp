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
library(tidyr)
library(qpdf)
library(loupeR)
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
    .dimplot-modal {
      padding: 0;
    }

    .dimplot-modal {
      width: 100%;
    }

    .dimplot-modal {
      max-width: 100%;
    }

    .modal-dialog.modal-xl:has(.dimplot-modal) {
      width: 80% !important;
      max-width: 80% !important;
    }

    .modal-header .close {
      display: none;
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
  # REACTIVE VALUES
  seuratObj <- reactiveVal(NULL)
  variable_method <- reactiveVal(NULL)
  resolution_search_results <- reactiveVal(NULL)
  umap_plot <- reactiveVal(NULL)
  heatmap_plot <- reactiveVal(NULL)
  cluster_table <- reactiveVal(NULL)
  table_data <- reactiveVal(NULL)
  dotplot1_plot <- reactiveVal(NULL)
  dotplot2_plot <- reactiveVal(NULL)
  current_ident <- reactiveVal(NULL)
  subset_data_path <- reactiveVal(NULL)
  plot_title <- reactiveVal(NULL)
  harmony_state <- reactiveVal(NULL)
  harmony_meta_state <- reactiveVal(NULL)
  user_nn_dims <- reactiveVal(NULL)
  uploaded_filename <- reactiveVal(NULL)
  user_resolution <- reactiveVal(0.10)
  temp_seuratObj <- reactiveVal(NULL)
  min_dist_default <- reactiveVal(NULL)
  umap_dist_state <- reactiveVal(NULL)
  combined_plot <- reactiveVal(NULL)
  uploaded_file_info <- reactiveVal(NULL)
  sample_ident_analysis_path <- reactiveValues(paths = list())
  show_download_print <- reactiveVal(FALSE)
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
      auto_min_dist <- function(n_cells) {
        md <- if (n_cells < 10000) 0.3
        else if (n_cells < 15000) 0.2
        else if (n_cells < 30000) 0.1
        else if (n_cells < 50000) 0.05
        else if (n_cells < 100000) 0.03
        else 0.01
        md
      }
      min_dist_default(auto_min_dist(ncol(obj)))
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
    variable_method(NULL)
    resolution_search_results(NULL)
    umap_plot(NULL)
    heatmap_plot(NULL)
    cluster_table(NULL)
    table_data(NULL)
    dotplot1_plot(NULL)
    dotplot2_plot(NULL)
    current_ident(NULL)
    subset_data_path(NULL)
    plot_title(NULL)
    harmony_state(NULL)
    harmony_meta_state(NULL)
    user_nn_dims(NULL)
    umap_dist_state <- reactiveVal(NULL)
    uploaded_filename(NULL)
    status$normalized <- FALSE
    status$scaled <- FALSE
    status$pca <- FALSE
    status$neighbor <- FALSE
    status$clustered <- FALSE
    temp_seuratObj(NULL)
    combined_plot(NULL)
    show_download_print(FALSE)
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
  
  # SAVE DATA MODAL
  observeEvent(input$save_data_btn, {
    if(is.null(input$title_for_download)){
      save_value <- paste0(tools::file_path_sans_ext(input$data_path), "_", Sys.Date(), ".rds")
      }else{
      save_value <- file.path(dirname(input$data_path), paste0(input$title_for_download, "_", Sys.Date(), ".rds"))}
    showModal(modalDialog(
      size = "l",
      title = "Save Data",
      fluidRow(
        column(12,
               checkboxInput("rds_file_chk", "RDS File", value = TRUE),
               checkboxInput("loupe_file_chk", "Loupe Browser File", value = TRUE),
               checkboxInput("cmd_history_chk", "Command History File", value = TRUE)
        ),
        uiOutput("optional_save_options_ui"),
        column(12,
               textInput(
                 "save_path",
                 "Input save path on SCC or local file name (use forward slashes (/). The filename should end with '.rds')",
                 value = save_value,
                 width = "100%"
               )
        ),
        column(12,
               tags$p("The \"Download Locally\" button will not display any pop-up window or notification. When downloading 
                      files (especially large files such as RDS or Loupe Browser files), please be patient and wait until the 
                      download is fully completed. Make sure the file has been successfully saved to your local system before 
                      performing any further actions.", style = "color: rgb(215, 47, 40); font-weight:bold; font-size:14px;"))
      ),
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
    
    output$optional_save_options_ui <- renderUI({
      opts <- c()
      if (!is.null(input$integration_harmony) && input$integration_harmony && !is.null(status$neighbor) && status$neighbor) {
        opts <- c(opts, "Harmony Plot (Integration)" = "harmony_plot")
      }
      if (!is.null(status$clustered) && status$clustered) {
        opts <- c(opts, "Combined Plot (umap, heatmap, dotplot)" = "combined_plot", "Gene Marker (for clustering)" = "gene_marker")
      }
      if (length(opts) == 0) return(NULL) 
      
      tags$div(
        style = "margin-left: 15px;", 
        checkboxGroupInput(
          "save_options_optional",
          "Other files to save:",
          choices = opts,
          selected = NULL
        )
      )
    })
  })
  
  
  observeEvent(input$btn_save, {
    req(seuratObj(), input$save_path)
    save_path <- input$save_path
    selected_files <- c()
    if (isTRUE(input$rds_file_chk)) selected_files <- c(selected_files, "rds_file")
    if (isTRUE(input$loupe_file_chk)) {selected_files <- c(selected_files, "loupe_file")}
    if (isTRUE(input$cmd_history_chk)) selected_files <- c(selected_files, "command_history")
    if (!is.null(input$save_options_optional)) selected_files <- c(selected_files, input$save_options_optional)
    
    showModal(modalDialog(title = "Please wait", 
                          tagList(
                            p("Saving..."),
                            p("The waiting time depends on dataset size. Please be patient...")
                          ), footer = NULL, easyClose = FALSE))
    
    tryCatch({
      if ("rds_file" %in% selected_files) {
        if (!is.null(temp_seuratObj())) saveRDS(temp_seuratObj(), file = save_path)
        else saveRDS(seuratObj(), file = save_path)
      }
      
      if ("loupe_file" %in% selected_files) {
        
        loupe_base <- tools::file_path_sans_ext(save_path)
        setup()
        create_loupe_from_seurat(
          if (!is.null(temp_seuratObj())) temp_seuratObj() else seuratObj(),
          output_name = loupe_base
        )
      }
      
      if ("command_history" %in% selected_files) {
        cmd_lines <- c(
          paste0("findvariable: ", input$findvariable %||% "variable"),
          paste0("npc: ", input$npc),
          paste0("nn_dims: ", input$nn_dims),
          paste0("integration_harmony: ", input$integration_harmony %||% FALSE),
          paste0("harmony_metadata: ", if (input$integration_harmony) input$harmony_metadata else "NULL"),
          paste0("resolution: ", input$resolution),
          paste0("umap_dist: ", input$umap_dist)
        )
        if (!is.null(input$sample_select)) {
          vals <- sort(unique(seuratObj()@meta.data[[input$sample_select]]))
          for (val in vals) {
            sel <- input[[paste0("sample_cb_", val)]] %||% FALSE
            rename_val <- input[[paste0("sample_rename_", val)]] %||% val
            cmd_lines <- c(cmd_lines, paste0("sample_", val, ": ", if (sel) paste0(val,"->",rename_val) else "unselected"))
          }
        }
        if (!is.null(current_ident())) {
          vals <- sort(unique(seuratObj()@meta.data[[current_ident()]]))
          for (val in vals) {
            sel <- input[[paste0("cluster_cb_", val)]] %||% FALSE
            rename_val <- input[[paste0("cluster_rename_", val)]] %||% val
            cmd_lines <- c(cmd_lines, paste0("cluster_", val, ": ", if (sel) paste0(val,"->",rename_val) else "unselected"))
          }
        }
        cmd_file <- paste0(tools::file_path_sans_ext(save_path), "_CommandHistory_", Sys.Date(), ".txt")
        writeLines(cmd_lines, cmd_file)
      }
      if ("harmony_plot" %in% selected_files) {
        harmony_plot <- DimPlot(
          seuratObj(),
          reduction = "umap",
          group.by = input$harmony_metadata,
          split.by = input$harmony_metadata
        ) + ggtitle(NULL)
        harmony_file <- paste0(file_path_sans_ext(save_path),paste0("HarmonyDimPlotSplit_", tools::file_path_sans_ext(basename(input$data_path)), "_", Sys.Date(), ".png"))
        ggsave(harmony_file, plot = harmony_plot, width = 20, height = 8, dpi = 300)
      }
      
      if ("combined_plot" %in% selected_files) {

        comb_plot <- create_combined_cluster_plot_patchwork(
          umap_plot(),
          table_data(),
          heatmap_plot(),
          dotplot1_plot(),
          dotplot2_plot()
        )
        comb_file <- paste0(file_path_sans_ext(save_path),paste0(input$title_for_download, ".pdf"))
        height_combined <- if(length(unique(heatmap_plot()$data$cluster)) <= 10){20}else{20+0.72*length(unique(heatmap_plot()$data$cluster))}
        ggsave(comb_file, plot = comb_plot, width = 18, height = height_combined)
      }
      
      if ("gene_marker" %in% selected_files) {
        gene_file <- if (!is.null(uploaded_file_info()) && file.exists(uploaded_file_info()$datapath)) {
          uploaded_file_info()$datapath
        } else {
          file.path(getwd(), "Gene_Markers.xlsx")
        }
        if (!file.exists(gene_file)) {
          showNotification("Gene Marker file not found!", type = "error")
        } else {
          gene_dest <- paste0(tools::file_path_sans_ext(save_path), "_GeneMarkers.xlsx")
          file.copy(gene_file, gene_dest, overwrite = TRUE)
        }
      }
      
      showNotification("Selected files saved successfully.", type = "message")
      removeModal()
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("Clustering_",file_path_sans_ext(basename(input$data_path)), "_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(seuratObj())
      tmp_dir <- tempdir()
      files_to_zip <- c()
      selected_files <- c()
      if (isTRUE(input$rds_file_chk)) selected_files <- c(selected_files, "rds_file")
      if (isTRUE(input$loupe_file_chk)) {selected_files <- c(selected_files, "loupe_file")}
      if (isTRUE(input$cmd_history_chk)) selected_files <- c(selected_files, "command_history")
      if (!is.null(input$save_options_optional)) selected_files <- c(selected_files, input$save_options_optional)
      
      if ("rds_file" %in% selected_files) {
        rds_path <- file.path(tmp_dir, basename(input$save_path))
        if (!is.null(temp_seuratObj())) saveRDS(temp_seuratObj(), file = rds_path)
        else saveRDS(seuratObj(), file = rds_path)
        files_to_zip <- c(files_to_zip, rds_path)
      }
      
      if ("loupe_file" %in% selected_files) {
        
        loupe_output <- file.path(
          tmp_dir,
          paste0(tools::file_path_sans_ext(basename(input$data_path)))
        )
        setup()
        create_loupe_from_seurat(
          if (!is.null(temp_seuratObj())) temp_seuratObj() else seuratObj(),
          output_name = loupe_output
        )
        
        loupe_file <- paste0(loupe_output, ".cloupe")
        
        if (file.exists(loupe_file)) {
          files_to_zip <- c(files_to_zip, loupe_file)
        } else {
          showNotification("Loupe file not found after creation!", type = "error")
        }
      }
      
      if ("command_history" %in% selected_files) {
        cmd_lines <- c(
          paste0("findvariable: ", input$findvariable %||% "variable"),
          paste0("npc: ", input$npc),
          paste0("nn_dims: ", input$nn_dims),
          paste0("integration_harmony: ", input$integration_harmony %||% FALSE),
          paste0("harmony_metadata: ", if (input$integration_harmony) input$harmony_metadata else "NULL"),
          paste0("resolution: ", input$resolution),
          paste0("umap_dist: ", input$umap_dist)
        )
        if (!is.null(input$sample_select)) {
          vals <- sort(unique(seuratObj()@meta.data[[input$sample_select]]))
          for (val in vals) {
            sel <- input[[paste0("sample_cb_", val)]] %||% FALSE
            rename_val <- input[[paste0("sample_rename_", val)]] %||% val
            cmd_lines <- c(cmd_lines, paste0("sample_", val, ": ", if (sel) paste0(val,"->",rename_val) else "unselected"))
          }
        }
        if (!is.null(current_ident())) {
          vals <- sort(unique(seuratObj()@meta.data[[current_ident()]]))
          for (val in vals) {
            sel <- input[[paste0("cluster_cb_", val)]] %||% FALSE
            rename_val <- input[[paste0("cluster_rename_", val)]] %||% val
            cmd_lines <- c(cmd_lines, paste0("cluster_", val, ": ", if (sel) paste0(val,"->",rename_val) else "unselected"))
          }
        }
        cmd_file <- file.path(tmp_dir, paste0(tools::file_path_sans_ext(basename(input$save_path)), "_CommandHistory_", Sys.Date(), ".txt"))
        writeLines(cmd_lines, cmd_file)
        files_to_zip <- c(files_to_zip, cmd_file)
      }
      
      if ("harmony_plot" %in% selected_files) {
        harmony_plot <- DimPlot(
          seuratObj(),
          reduction = "umap",
          group.by = input$harmony_metadata,
          split.by = input$harmony_metadata
        ) + ggtitle(NULL)
        harmony_file <- file.path(tmp_dir, paste0("HarmonyDimPlotSplit_", tools::file_path_sans_ext(basename(input$data_path)), "_", Sys.Date(), ".png"))
        ggsave(harmony_file, plot = harmony_plot, width = 20, height = 8, dpi = 300)
        files_to_zip <- c(files_to_zip, harmony_file)
      }
      
      if ("combined_plot" %in% selected_files) {
        comb_plot <- create_combined_cluster_plot_patchwork(
          umap_plot(),
          table_data(),
          heatmap_plot(),
          dotplot1_plot(),
          dotplot2_plot()
        )
        comb_file <- file.path(tmp_dir, paste0(input$title_for_download, ".pdf"))
        height_combined <- if(length(unique(heatmap_plot()$data$cluster)) <= 10){20}else{20+0.72*length(unique(heatmap_plot()$data$cluster))}
        ggsave(comb_file, plot = comb_plot, width = 18, height = height_combined)
        files_to_zip <- c(files_to_zip, comb_file)
      }
      if ("gene_marker" %in% selected_files) {
        gene_file <- if (!is.null(uploaded_file_info()) && file.exists(uploaded_file_info()$datapath)) {
          uploaded_file_info()$datapath
        } else {
          file.path(getwd(), "Gene_Markers.xlsx")
        }
        if (!file.exists(gene_file)) {
          showNotification("Gene Marker file not found for download!", type = "error")
        } else {
          gene_dest <- file.path(tmp_dir, paste0(tools::file_path_sans_ext(basename(input$data_path)), "_Gene_Marker.xlsx"))
          file.copy(gene_file, gene_dest, overwrite = TRUE)
          files_to_zip <- c(files_to_zip, gene_dest)
        }
      }
      
      zip::zip(zipfile = file, files = files_to_zip, mode = "cherry-pick")
    },
    contentType = "application/zip"
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
        tags$p("🔴 SELECT FEATURE AND DATA SCALED", style = "color:red; font-weight:bold; font-size:14px; margin:0;")
      ))
      stop_flag <- TRUE
    } else {
      ui_list <- append(ui_list, list(
        div(
          style = "display:flex; align-items:center; gap:6px;",
        if (isTRUE(status$scaled)) {
          tags$p("✅ SELECT FEATURE AND DATA SCALED", style = "color:black; font-weight:bold; font-size:14px; margin:0;")
        } else {
          tags$p("🔘 SELECT FEATURE AND DATA SCALED", style = "color:black; font-weight:bold; font-size:14px; margin:0;")
        },
        div(class = "tooltip-circle",
            "?",  
            span(class = "tooltip-text", 
                 div(class = "tooltip-title", "SELECT FEATURE AND DATA SCALED"),
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
                     After that, it will provide visualization for the clustering result. Note that the cluster data 
                     is based on the neighbor data from the previous step. In other words, if you selected Harmony 
                     for your neighbor data, then the cluster data will also be generated using Harmony's data.")
            )
        )
        ),
        uiOutput("clustering_ui")
      ))
    }
    
    ui_list <- append(ui_list, list(hr()))
    
    # --- Rename and Subset ---
    if (stop_flag|| status$clustered == FALSE) {
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
             After subsetting, click \"Subset and Recluster Data\". The system will reset ShinyApp and use the subset data as 
             input for analysis.<br>
             It is recommended to save the data after both renaming and subsetting."))
            )
          )
        ),
        
        uiOutput("final_ui")
      ))
    }
    ui_list <- append(ui_list, list( div(style = "height: 50px;")))
    
    tagList(ui_list)
  })

  # NORMALIZATION STEP UI
  output$normalized_data_ui <- renderUI({
    req(seuratObj())
    div(
      style = "display: inline-block; vertical-align: middle;",
      actionButton(
        "normalize_btn_ui",
        "Normalize Data",
        width = "200px",
        class = "btn-primary"
      ))
  })
  
  # SCALING STEP UI
  output$scaled_data_ui <- renderUI({
    tagList(
      div(
        style = "display:flex; align-items:center; gap:6px;",
        tags$h5("Feature Gene Selection Method:", style = "color:black; font-weight:bold; font-size:14px; margin:0;"),
        div(class = "tooltip-circle",
            "?",  
            span(class = "tooltip-text", 
                 div(class = "tooltip-title", "Feature Gene Selection Method"),
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
            style = "margin-left:20px;width = 100%;",
            uiOutput("findvariable_ui")
          )
      ),
      div(
        style = "display: inline-block; vertical-align: middle;",
        actionButton(
          "scaling_btn_ui",
          "Scale Data",
          width = "200px",
          class = "btn-primary"
        ))
    )
  })
  
  output$findvariable_ui <- renderUI({
    
    has_variable <- !is.null(seuratObj()) && length(VariableFeatures(seuratObj())) > 0
    
    choices <- c(
      "Use Top 2,000 Most Variable Features (Recommended)" = "variable",
      "Use All Genes (Long Processing Time!!!)" = "all",
      "Use Top 2,000 Integration Features (Multi Samples)" = "integration"
    )
    
    if (has_variable) {
      choices <- c(
        choices,
        "Use Original Set of Variable Features (Useful when Subsetting, to maintain Harmony across samples)" = "original"
      )
    }
    
    tagList(
      radioButtons(
        "findvariable",
        NULL,
        choices = choices,
        selected = input$findvariable %||% "variable",
        inline = FALSE,
        width = "100%"
      ),
      
      conditionalPanel(
        condition = "input.findvariable == 'integration'",
        
        div(
          style = "margin-left: 20px; margin-top: 5px;",
          
          selectInput(
            "integration_meta",
            label = "Select Sample column name",
            choices = colnames(seuratObj()@meta.data),
            selected = NULL
          )
        )
      )
    )
  })
  
  # PCA STEP UI
  output$pca_ui <- renderUI({
    req(seuratObj())
    gene_method <- if (!is.null(input$findvariable)) {
      if (input$findvariable == "variable") {
        "Variable Features"
      } else if (input$findvariable == "original") {
        "Original Variable Features"
      } else if (input$findvariable == "integration") {
        "Integration Features"}else {
        "All Genes"
      }
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
      div(
        style = "display: inline-block; vertical-align: middle;",
        actionButton(
          "pca_btn_ui",
          "Run PCA",
          width = "200px",
          class = "btn-primary"
        ))
    )
  })
  
  debounced_nn_dims <- reactive(input$nn_dims) %>% debounce(3000)
  
  observeEvent(debounced_nn_dims(), {
    user_nn_dims(debounced_nn_dims())
  })
  
  # NEIGHBOR STEP UI
  output$neighbor_ui <- renderUI({ 
    req(seuratObj())
    
    recommend_value <- if(!is.null(seuratObj()@reductions$pca)) {
      recommended_dimensions()$recommended
    } else {
      NULL
    }
    current_value <- user_nn_dims() %||% recommend_value %||% 1
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
                       value = current_value,
                       min = 1, max = input$npc %||% 50, step = 1, width = "100%"))
      ),
      fluidRow(
        div(
          style = "display: flex; margin-left:30px;",
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
          ),
          div(
            class = "harmony-checkbox",
            checkboxInput(
              inputId = "integration_harmony",
              label = "Harmony Integration (If you have multi samples in the data)",
              value = ifelse(!is.null(harmony_state()), harmony_state(), FALSE),
              width = "350px"
            ) %>% 
              tagAppendAttributes(style = "white-space: nowrap; display: inline-flex; align-items:center;")
          )
          
        )
      ),
      
      conditionalPanel(
        condition = "input.integration_harmony == true",
          div(
            style = "display: inline-block; vertical-align: middle; margin-left:20px;",
            uiOutput("harmony_metadata_ui"))),
      
      fluidRow(
        div(
          style = "display: inline-block; vertical-align: middle; margin-left: 20px;",
          actionButton(
            "neighbor_btn_ui",
            "Find Neighbors",
            width = "200px",
            class = "btn-primary"
          )
        ),
        div(
          style = "display: inline-block; vertical-align: middle; margin-left: 20px;",
          uiOutput("harmony_check_ui")
        )
      )
      
      )
    })
  
  output$harmony_check_ui <- renderUI({
    if (isTRUE(status$neighbor)) {
        actionButton(
          "open_dimplot_modal",
          "Check Harmony Result",
          width = "200px",
          class = "btn-danger"
      )
    }
  })
  
  observeEvent(input$open_dimplot_modal, {
    
    showModal(
      modalDialog(
        title = "Harmony DimPlot",
        size = "xl",
        easyClose = FALSE,
        footer = NULL,
        
        tags$div(
          class = "dimplot-modal",
          ## ===== plot =====
          plotOutput("dimplot_modal", height = "600px"),
          hr(),
          ## ===== split.by + close =====
          tags$div(
            style = "
            display: flex;
            align-items: center;
            justify-content: space-between;
          ",
            
            tags$div(
              style = "display: flex; align-items: center;",
              checkboxInput(
                "split_by",
                label = paste0("Split by ", input$harmony_metadata),
                value = TRUE
              )
            ),
            
            downloadButton(
              "download_dimplot_modal",
              "Download Plot",
              style = "background-color: #4CAF50; color: white;"
            ),
            
            modalButton("Close")
          )
        )))})
  
  output$dimplot_modal <- renderPlot({
    
    if (input$split_by) {
      DimPlot(
        seuratObj(),
        reduction = "umap",
        group.by = input$harmony_metadata,
        split.by = input$harmony_metadata
      ) + ggtitle(NULL)
    } else {
      DimPlot(
        seuratObj(),
        reduction = "umap",
        group.by = input$harmony_metadata,
      ) + ggtitle(NULL)
    }
  })

  output$download_dimplot_modal <- downloadHandler(
    filename = function() {
      paste0("HarmonyDimPlotSplit_", tools::file_path_sans_ext(basename(input$data_path)), "_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(seuratObj(), input$harmony_metadata)
      
      if (isTRUE(input$split_by)) {
        plot <- DimPlot(
          seuratObj(),
          reduction = "umap",
          group.by = input$harmony_metadata,
          split.by = input$harmony_metadata
        ) + ggtitle(NULL)
        ggsave(
          filename = file,
          plot = plot,
          width = 20,
          height = 8,
          dpi = 300
        )
      } else {
        plot <- DimPlot(
          seuratObj(),
          reduction = "umap",
          group.by = input$harmony_metadata
        ) + ggtitle(NULL)
        ggsave(
          filename = file,
          plot = plot,
          width = 12,
          height = 8,
          dpi = 300
        )
      }
      
      showNotification("DimPlot saved to local", type = "message")
    }
  )
  
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
                             placeholder = if (is.null(uploaded_file_info())) {
                               "No file selected, Using default"
                             } else {
                               paste("Selected:", uploaded_file_info()$name)
                             },
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
                                    value = 0.01, min = 0.0001, max = 5.0, step = 0.0001, width = "100px")
                     ),
                     div(
                       style = "display: inline-block; vertical-align: middle; margin-left:20px;",
                       h5("End :", style = "font-weight: bold; display: inline-block; margin-right:8px;"),
                       numericInput("end_resolution", NULL, 
                                    value = 0.5, min = 0, max = 5.0, step = 0.01, width = "100px")
                     ),
                     div(
                       style = "display: inline-block; vertical-align: middle; margin-left:20px;",
                       h5("Step :", style = "font-weight: bold; display: inline-block; margin-right:8px;"),
                       numericInput("step_resolution", NULL, 
                                    value = 0.02, min = 1e-10, max = 1, step = 0.001, width = "100px")
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
               div(
                 style = "display: inline-block; vertical-align: middle;",
                 actionButton(
                   "clustering_btn_ui",
                   "Find Clusters",
                   width = "200px",
                   class = "btn-primary"
                 )
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
  
  debounced_resolution <- reactive(input$resolution) %>% debounce(3000)
  
  observeEvent(debounced_resolution(), {
    user_resolution(debounced_resolution())
  })
  
  # RESOLUTION INPUT UI
  output$resolution_ui <- renderUI({
    numericInput("resolution", NULL, 
                 value = user_resolution() %||% 0.10, 
                 min = 0.001, max = 5.000, step = 0.001, width = "100%")
  })
  
  
  # FINAL STEP UI
  output$final_ui <- renderUI({
    req(seuratObj(),input$resolution)
    if(is.null(current_ident())){
      val <- "shiny_clusters"
      current_ident(val)
    }
    mdcols <- names(sapply(seuratObj()@meta.data, function(x) !is.numeric(x))[sapply(seuratObj()@meta.data, function(x) !is.numeric(x))])
    mdcols <- setdiff(mdcols, "CB")
    tagList(
      fluidRow(
        div(style = "display: inline-block; vertical-align: top; margin-left:20px;",
            h5("Cluster metadata column name:",
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
      fluidRow(
        column(6,
               h4("Subsetting Samples"),
               fluidRow(
               div(style = "display: inline-block; vertical-align: top; margin-left:20px;",
                   h5("Sample:",style = "font-weight: bold;")),
               div(style = "display: inline-block; vertical-align: top; margin-left:20px;",
                   selectInput("sample_select", NULL, choices = mdcols, selected = if(!is.null(input$harmony_metadata)){input$harmony_metadata}else{NULL}))),
               plotOutput("sample_umap_plot", width = "80%"),
               uiOutput("sample_subset_rename_ui"),  # Combined UI
               fluidRow(
                 div(
                   style = "display: flex; flex-direction: column; gap: 20px; width: 500px;",
                   div(
                     style = "display: flex; align-items: center; gap: 10px;",
                     div(
                       class = "tooltip-circle",
                       "?",
                       span(
                         class = "tooltip-text",
                         div(class = "tooltip-title", 
                             "Analyze Individual Samples (No Reclustering)"),
                         div(class = "tooltip-divider"),
                         div(
                           "12. Clicking \"Analyze Individual Samples (No Reclustering)\" will analyze the selected sample(s)
                           without reclustering and give the opportunity to \"Download Result\" and print to screen UMAP,
                           data table, heatmap, and dotplot for each sample. You can download the results to your computer.
                           This allows the user to observe the data for each sample within the original dataset. If you want
                           to recluster and then reanalyze one or more samples individually, select the sample you want and
                           click on \"Subset & Rename\" button, and then, if desired, \"Save Data\"."
                         )
                       )
                     ),
                     actionButton(
                       "sample_ident_analysis",
                       "Analyze Individual Samples (No Reclustering)",
                       width = "70%",
                       class = "btn-success"
                     )
                   ),
                   div(
                     style = "margin-left: 50px;",
                     uiOutput("download_ident_sample_ui")
                   )
                 )
               )),
               column(6,
                      h4("Subsetting Celltype Clusters"),
                      h4(paste0("Cluster: ", current_ident())),
                      h4(""),
                      plotOutput("cluster_umap_plot", width = "80%"),
                      uiOutput("cluster_subset_rename_ui")  # Combined UI
               )),
            fluidRow(
              column(4),
              
              column(
                4,
                h5(
                  "Attention: If you want to retain the original data, please use \"SAVE DATA\" to save your current result before subset.!!!",
                  style = "font-weight: bold; color:red;"
                ),
                div(
                  style = "display: flex; flex-direction: column; gap: 10px;",
                  actionButton(
                    "subset_btn_ui",
                    "Rename & Subset",
                    width = "100%",
                    class = "btn-danger"
                  ),
                  conditionalPanel(
                    condition = "input.confirm_subset > 0",
                    actionButton(
                      "keep_analysis",
                      "Subset and Recluster Data",
                      width = "100%",
                      class = "btn-success"
                    )
                  )
                )
              ),
              
              column(4)
            ))
  })
  output$download_ident_sample_ui <- renderUI({
    if (show_download_print()) {
      tagList(
        div(
          style = "display: flex; flex-direction: column; gap: 5px; width: 350px;",
          downloadButton(
            "download_ident_sample_analysis_combined",
            "Download (combined file)",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          ),
          downloadButton(
            "download_ident_sample_analysis_separate",
            "Download (separate file for each sample)",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          ),
          actionButton(
            "print_ident_sample_analysis",
            "Print on Screen",
            class = "btn-warning",
            style = "width: 100%;"
          )
        )
      )
    } else {
      NULL
    }
  })
  # SCC JOB MODAL
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
          "Project name (project name on scc, like wax-es or wax-dk):", 
          placeholder = "e.g., wax-es"
        ),
        
        # Project name
        textInput(
          "scc_pca_num", 
          "PCA number (Multiple PCAs are separated by commas(,).:", 
          value = input$nn_dims,
          placeholder = "e.g., 7,8,9"
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
              value = if(length(seq(input$step_resolution, input$end_resolution, by = input$step_resolution)) <= 10) {
                1
              }else{8},
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
  
  # SUBMIT SCC JOB ACTION
  observeEvent(input$run_scc_job, {
    req(seuratObj(), input$scc_save_path, input$scc_project_name, input$scc_runtime, input$scc_cores, input$start_resolution, input$step_resolution, input$end_resolution)
    
    showModal(modalDialog(
      title = "Please wait",
      tagList(
        p("Submitting SCC job..."),
        p("Don't close the window." )
      ),
      footer = NULL,
      easyClose = FALSE
    ))
    marker_path <- if (!is.null(input$upload_excel)) {
      input$upload_excel$datapath
    } else {
      file.path(getwd(), "Gene_Markers.xlsx")
    }
    tryCatch({
      # Call function to submit SCC job

      source(file.path(getwd(), "scc_job_submission.R"))
      min.dist <- if(is.null(input$umap_dist)){input$umap_dist}else{min_dist_default()}
      qsub_job_number <- submit_scc_job(
        seuratObj(),
        save_path = input$scc_save_path,
        project_name = input$scc_project_name,
        start_resolution = input$start_resolution,
        step_resolution = input$step_resolution,
        end_resolution = input$end_resolution,
        integration_harmony = input$integration_harmony,
        pca_num = input$scc_pca_num,
        harmony_metadata = if (input$integration_harmony) {
          input$harmony_metadata
        } else {
          NULL
        },
        runtime = input$scc_runtime,
        cores = input$scc_cores,
        email = input$scc_email,
        marker_path = marker_path
      )
      removeModal()
      showModal(
        modalDialog(
          title = "Clustering Visualization Job Submitted",
          tagList(
            p("The Clustering analysis is now running in the background."),
            p("Your job ID is:"),
            tags$pre(style = "color: red; user-select: text;", qsub_job_number),
            p("You can copy the following command in the SCC terminal to check the job status:"),
            tags$pre(style = "color: red; user-select: text;", paste0("qstat -j ", qsub_job_number)),
            p("You can now close the shinyapp and wait for the results."),
            p("For subsequent analysis, it is recommended to use the following path of seurat rds file as input:"),
            tags$pre(style = "color: blue; user-select: text;", 
                     file.path(input$scc_save_path, "input_seurat_file.rds")),
            if (nchar(input$scc_email) > 0) {
              p("The results will be sent to your email address: ", input$scc_email)
            } else {
              NULL
            }
          ),
          footer = modalButton("Close"),
          easyClose = TRUE
        )
      )
      
      showNotification("SCC job submitted successfully.", type = "message")
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  
  # RENAME IDENTS SERVER               
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
      temp_seuratObj(NULL)
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
  
  # SCALING STEP SERVER
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
      } else if(input$findvariable == "all")  {
        seuratObj(ScaleData(seuratObj(), features = rownames(seuratObj()), verbose = FALSE))
      } else if(input$findvariable == "integration"){
        obj <- seuratObj()  
        seurat_list <- SplitObject(obj, split.by = input$integration_meta)
        seurat_list <- lapply(seurat_list, FindVariableFeatures)
        features <- SelectIntegrationFeatures(seurat_list)
        VariableFeatures(obj) <- features 
        obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
        seuratObj(obj)
        }else if(input$findvariable == "original"){
        seuratObj(ScaleData(seuratObj(), features = VariableFeatures(seuratObj()), verbose = FALSE))
      }
      status$normalized <- TRUE
      status$scaled <- TRUE
      status$pca <- FALSE
      status$neighbor <- FALSE
      status$clustered <- FALSE
      variable_method(switch(input$findvariable,
                             "variable" = "Variable Features",
                             "all" = "All Genes",
                             "integration" = "Integration Features",
                             "original" = "Original Features"))
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
      if (input$findvariable == "all") {
        seuratObj(suppressWarnings(RunPCA(seuratObj(), features = rownames(seuratObj()), npcs = input$npc, verbose = FALSE)))
      } else {
        seuratObj(suppressWarnings(RunPCA(seuratObj(), features = VariableFeatures(seuratObj()), npcs = input$npc, verbose = FALSE)))
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
  
  # DISPLAY RECOMMENDATION RESULTS
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
  
  # HARMONY INTEGRATION VALUE UPDATE
  observeEvent(input$integration_harmony, {
    harmony_state(input$integration_harmony)
  })
  
  observeEvent(seuratObj(), {
    obj <- seuratObj()
    has_harmony <- "harmony" %in% names(obj@reductions)
    if (is.null(harmony_state())) {
      harmony_state(has_harmony)
    }
  })
  
  observeEvent(input$harmony_metadata, {
    harmony_meta_state(input$harmony_metadata)
  })
  
  # Harmony metadata UI
  output$harmony_metadata_ui <- renderUI({
    req(seuratObj())
    
    md <- seuratObj()@meta.data
    metadata_cols <- colnames(md)
    
    metadata_cols <- metadata_cols[!grepl("CB", metadata_cols, ignore.case = TRUE)]
    numeric_filter <- sapply(md[ , metadata_cols, drop = FALSE], function(x) {
      is.numeric(x) && length(unique(x)) > 100
    })
    metadata_cols <- metadata_cols[!numeric_filter]
    default_col <- if (any(grepl("^SampleID", metadata_cols, ignore.case = TRUE))) {
      metadata_cols[grepl("^SampleID", metadata_cols, ignore.case = TRUE)][1]
    } else {
      metadata_cols[1]
    }
    selected_col <- isolate(harmony_meta_state() %||% default_col)
    tagList(
      div(
        style = "display: inline-block; vertical-align: middle; margin-left:20px;",
        h5("Select metadata column for Harmony integration (sample identification):",
           style = "font-weight: bold; display: inline-block; margin-right:8px;"),
        div(
          style = "margin-right:10px; display:inline-block;",
          selectInput(
            "harmony_metadata", NULL,
            choices = metadata_cols,
            selected = selected_col,
            width = "200px"
          )
        )
      ),
      div(
        style = "margin-left:40px;",
        textOutput("harmony_metadata_preview")
      ),
      br()
    )
  })
  
  # Harmony metadata preview
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
      if(input$integration_harmony == TRUE){
        req(input$harmony_metadata)
        seuratObj <- RunHarmony(seuratObj(), group.by.vars = input$harmony_metadata, reduction.use = "pca",
                                dims.use = 1:input$nn_dims, assay.use = "RNA", reduction.save = "harmony", 
                                project.dim = TRUE, verbose = FALSE)
        seuratObj <- FindNeighbors(seuratObj, reduction = "harmony", dims = 1:input$nn_dims, verbose = FALSE)
        seuratObj <- RunUMAP(seuratObj, min.dist = min_dist_default(), dims = 1:input$nn_dims, reduction = "harmony", verbose = FALSE)
      } else {
        seuratObj <- FindNeighbors(seuratObj(), reduction = "pca", dims = 1:input$nn_dims, verbose = FALSE)
        seuratObj <- RunUMAP(seuratObj, min.dist = min_dist_default(), dims = 1:input$nn_dims, reduction = "pca", verbose = FALSE)
      }
      seuratObj(seuratObj)
      
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
  
  create_combined_cluster_plot_patchwork <- function(umap_plot, table_data, heatmap_plot, dotplot1 = NULL, dotplot2 = NULL, title = NA) {
    
    table_plot <- function(table_data) {
      if (is.null(table_data)) return(ggplot() + theme_void())
      
      table_grob <- gridExtra::tableGrob(
        table_data,
        rows = NULL,
        theme = gridExtra::ttheme_minimal(
          base_size = 11,
          padding = unit(c(2, 2), "mm")
        )
      )
      
      ggplot() + 
        annotation_custom(table_grob) + 
        theme_void() +
        theme(plot.margin = margin(5, 5, 5, 5))
    }
    
    cluster_n <- sum(table_data$cluster != "total", na.rm = TRUE)
    top_row <- (umap_plot | table_plot(table_data) | heatmap_plot) + 
      plot_layout(widths = c(4, 3, 5))
    title_pdf <- if(is.na(title)){input$title_for_download %||% plot_title()}else{title}
    title_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = title_pdf,
               size = 8,
               hjust = 0.5, vjust = 0.5) +
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
    height_heatmap <- if(length(unique(heatmap_plot()$data$cluster)) <= 10){2}else{length(unique(heatmap_plot()$data$cluster))*0.2}
    heights <- c(0.15, height_heatmap, rep(1.5, length(bottom_plots)))
    
    combined_plot <- combined_plot + 
      plot_layout(heights = heights[1:n_parts], ncol = 1)
    
    return(combined_plot)
  }
  
  # Annotate clusters if gene markers provided
  FindCellTypesByMarkers <- function(sobj, biomarkers = NULL) {
    if (is.null(biomarkers)) {
      stop("parameter biomarkers should be specified (now NULL)")
    }
    
    DefaultAssay(sobj) <- "RNA"
    
    clusters <- levels(Idents(sobj))
    all_genes <- unique(unlist(biomarkers))
    
    fc_all <- purrr::map_dfr(clusters, function(clust) {
      fc <- tryCatch({
        FoldChange(sobj, ident.1 = clust, features = all_genes, assay = "RNA", slot = "data")
      }, error = function(e) {
        data.frame(avg_log2FC = rep(NA_real_, length(all_genes)), row.names = all_genes)
      })
      
      fc_df <- as.data.frame(fc)
      fc_df$gene <- rownames(fc_df)
      fc_df$cluster <- clust
      fc_df[, c("cluster", "gene", "avg_log2FC")]
    })
    
    celltypes_ordered <- names(biomarkers)
    
    scores <- purrr::map_dfr(seq_along(biomarkers), function(i) {
      celltype <- celltypes_ordered[i]
      genes <- biomarkers[[i]]
      tmp <- fc_all %>%
        filter(gene %in% genes) %>%
        group_by(cluster) %>%
        summarise(score = mean(avg_log2FC, na.rm = TRUE), .groups = "drop") %>%
        mutate(celltype = celltype)
      tmp
    })
    
    scores <- scores %>%
      mutate(
        celltype = factor(celltype, levels = celltypes_ordered),
        cluster = factor(cluster, levels = clusters)
      )
    
    scores <- tidyr::expand(scores, cluster, celltype) %>%
      left_join(scores, by = c("cluster", "celltype")) %>%
      mutate(score = coalesce(score, 0)) %>%
      arrange(cluster, celltype)
    
    labels <- scores %>%
      group_by(cluster) %>%
      summarise(prediction = ifelse(max(score) > 0, as.character(celltype[which.max(score)]), "Unknown"), .groups = "drop") %>%
      mutate(
        prediction = paste0(prediction, "(", cluster, ")"),
        cluster = factor(cluster, levels = clusters)
      ) %>%
      arrange(cluster)
    
    new_labels <- labels$prediction
    names(new_labels) <- labels$cluster
    
    old_labels <- setNames(names(new_labels), new_labels)
    
    list(new_labels = new_labels, old_labels = old_labels, heatmap_table = scores)
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
      stop(paste("Failed to convert Excel to marker list. Please check Excel format.
Details:", e$message))
    })
  }
  
  # Create Heatmap Function
  CreateCellTypesHeatmap <- function(df, labels_text_size = 6, xaxis_text_size = 12, yaxis_text_size = 12, rotate_x = FALSE){
    
    clusters_numeric <- sort(as.numeric(levels(df$cluster)), decreasing = TRUE)
    df$cluster <- factor(df$cluster,
                         levels = as.character(clusters_numeric),
                         ordered = TRUE)
    
    tmp <- ggplot(df, aes(x = celltype, y = as.factor(cluster))) +
      geom_raster(aes(fill = score))+
      scale_fill_gradient2(low = "blue", mid="white", high = "tomato") +
      geom_text(aes(label=round(score,1)), size = labels_text_size) +
      scale_x_discrete(position = "top") +
      ylab("Clusters") +
      theme_minimal() +
      theme(legend.title = element_blank(),
            legend.position = "bottom",
            axis.text.x = element_text(color = "black", size = xaxis_text_size),
            axis.text.y = element_text(color = "black", size = yaxis_text_size)) +
      coord_fixed()
    
    if (rotate_x){
      tmp <- tmp + theme(axis.text.x = element_text(angle = 90, hjust=0, color = "black", size = xaxis_text_size),
                         axis.text.x.top = element_text(vjust = 0.5))
    }
    tmp
  }
  
  custom_dotplot <- function(
    seurat_obj,
    gene_list,
    all_gene_list,
    assay = "RNA",
    title = NULL,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5
  ) {
    
    plot_genes <- unlist(gene_list, use.names = FALSE)
    all_genes <- unique(unlist(all_gene_list, use.names = FALSE))
    
    DefaultAssay(seurat_obj) <- assay
    clusters <- unique(Idents(seurat_obj))
    
    results <- data.frame()
    
    for(cluster in clusters) {
      
      cluster_cells <- WhichCells(seurat_obj, idents = cluster)
      
      cluster_data <- GetAssayData(
        seurat_obj,
        assay = assay,
        slot = "data"
      )[all_genes, cluster_cells, drop = FALSE]
      
      for(gene in all_genes) {
        
        if(gene %in% rownames(cluster_data)) {
          
          avg_exp <- mean(expm1(cluster_data[gene, ]))
          pct_exp <- sum(cluster_data[gene, ] > 0) / length(cluster_cells) * 100
          
          results <- rbind(
            results,
            data.frame(
              gene = gene,
              cluster = cluster,
              avg_exp = avg_exp,
              pct_exp = pct_exp
            )
          )
        }
      }
    }
    results <- results %>%
      group_by(gene) %>%
      mutate(
        avg_exp_scaled = scale(log1p(avg_exp))[,1]
      ) %>%
      ungroup()
    
    results$avg_exp_scaled <- pmin(
      pmax(results$avg_exp_scaled, col.min),
      col.max
    )
    
    results_plot <- results %>%
      filter(gene %in% plot_genes)
    
    clusters_numeric <- sort(as.numeric(unique(results_plot$cluster)), decreasing = TRUE)
    
    results_plot$cluster <- factor(
      results_plot$cluster,
      levels = as.character(clusters_numeric),
      ordered = TRUE
    )
    
    gene_groups <- data.frame()
    
    for(group_name in names(gene_list)) {
      
      group_genes <- gene_list[[group_name]]
      
      gene_groups <- rbind(
        gene_groups,
        data.frame(gene = group_genes, group = group_name)
      )
    }
    
    gene_groups$group <- factor(gene_groups$group, levels = names(gene_list))
    
    results_plot <- results_plot %>%
      left_join(gene_groups, by = "gene")
    
    results_plot$gene <- factor(results_plot$gene, levels = plot_genes)
    
    p <- ggplot(results_plot, aes(x = gene, y = cluster)) +
      geom_point(aes(size = pct_exp, color = avg_exp_scaled)) +
      scale_color_gradientn(
        colors = cols,
        name = "Avg Exp",
        limits = c(col.min, col.max),
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
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
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
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "grey95", color = NA)
      ) +
      labs(x = "Genes", y = "Clusters", title = title)
    
    return(p)
  }
  get_cluster_summary_table <- function(seurat_obj, meta_col = NA) {

    if(is.na(meta_col)){    
      meta <- data.frame(
      cluster = as.character(Idents(seurat_obj)),
      counts = seurat_obj$nCount_RNA,
      genes = seurat_obj$nFeature_RNA)
      }else{
      meta <- data.frame(
      cluster = as.character(seurat_obj@meta.data[[meta_col]]),
      counts = seurat_obj$nCount_RNA,
      genes = seurat_obj$nFeature_RNA)
    }
    all_cells <- meta %>%
      group_by(cluster) %>%
      summarise(
        ncells = n(),
        avg.counts = as.integer(round(mean(counts))),
        avg.genes = as.integer(round(mean(genes)))
      ) %>%
      ungroup()
    
    result <- all_cells %>%
      mutate(all = sum(ncells),
             pct = round(100 * ncells / all, 2)) %>%
      select(cluster, ncells, pct, avg.counts, avg.genes)
    
    result <- bind_rows(
      result,
      result %>%
        summarise(
          cluster = "total",
          ncells = sum(ncells),
          pct = sum(pct),
          avg.counts = as.integer(round(mean(seurat_obj@meta.data$nCount_RNA))),
          avg.genes = as.integer(round(mean(seurat_obj@meta.data$nFeature_RNA)))
        )
    )
    
    return(result)
  }
  
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
      Idents(sobj) <- factor(Idents(sobj), 
                             levels = sort(as.numeric(levels(Idents(sobj)))))
      # Read Gene Markers
        path <- if (!is.null(input$upload_excel)) {
          input$upload_excel$datapath
        } else {
          file.path(getwd(), "Gene_Markers.xlsx")
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
          biomarkers_dotplot2 <- NULL
        })
        biomarkers_dotplot <- c(biomarkers_dotplot1, biomarkers_dotplot2)
        
        if (!is.null(biomarkers_dotplot1)) {
          dotplot1_plot(custom_dotplot(sobj, biomarkers_dotplot1, biomarkers_dotplot))
        }
        
        if (!is.null(biomarkers_dotplot2)) {
          dotplot2_plot(custom_dotplot(sobj, biomarkers_dotplot2, biomarkers_dotplot))
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
      labs <- levels(Idents(sobj))
      cluster_num <- as.numeric(sub(".*\\((\\d+)\\).*", "\\1", labs))
      new_levels <- labs[order(cluster_num)]
      Idents(sobj) <- factor(Idents(sobj), levels = new_levels)
      sobj$shiny_clusters <- Idents(sobj)
      base_name <- if (!is.null(subset_data_path())) {
        paste0("Subset: ", file_path_sans_ext(basename(subset_data_path())))
      } else {
        file_path_sans_ext(basename(input$data_path))
      }
      
      title <- paste0(
        base_name,
        ", resolution_", input$resolution,
        ", min_dist_", min_dist_default(),
        if ("harmony" %in% names(seuratObj()@reductions)) ", Harmony" else ""
      )
      
      plot_title(title)
      
      umap_plot(DimPlot(sobj, reduction = "umap", label = FALSE) + 
                  ggtitle("") +
                  theme(plot.title = element_text(face = "bold", size = 10),
                        legend.text = element_text(size = 10))+
                  guides(color = guide_legend(ncol = 1, 
                                              override.aes = list(size = 3))))
      df <- get_cluster_summary_table(sobj, "seurat_clusters")
      
      colnames(df) <- c("cluster", "ncell", "pct", "avg.counts", "avg.genes")
      table_data(df)
      combined_plot <- create_combined_cluster_plot_patchwork(
        umap_plot(),
        table_data(),
        heatmap_plot(),
        dotplot1_plot(),
        dotplot2_plot()
      )
      combined_plot(combined_plot)
      
      seuratObj(sobj)
      status$normalized <- TRUE
      status$scaled <- TRUE
      status$pca <- TRUE
      status$neighbor <- TRUE
      status$clustered <- TRUE
      temp_seuratObj(NULL)
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
                     produce a more spread-out representation. 
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
  
  debounced_dist <- reactive(input$umap_dist) %>% debounce(3000)
  
  observeEvent(debounced_dist(), {
    umap_dist_state(input$umap_dist)
  })
  
  #   UMAP DISTANCE INPUT UI
  output$umap_dist_ui <- renderUI({
    numericInput("umap_dist", NULL, 
                 value = umap_dist_state() %||% min_dist_default(), 
                 min = 0.001, max = 1, step = 0.001, width = "100%")
  })
  
  # DOWNLOAD EXCEL TEMPLATE
  output$download_template <- downloadHandler(
    filename = function() {
      "Gene_Markers.xlsx"
    },
    content = function(file) {
      template_path <- file.path(getwd(), "Gene_Markers.xlsx")
      if (file.exists(template_path)) {
        file.copy(template_path, file)
      }
    }
  )
  
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
      duplicates_within <- list()
      duplicates_across <- list()
      
      sheets_to_check <- sheets[grepl("^DOTPLOT", sheets, ignore.case = TRUE)]
      
      for (sheet in sheets_to_check) {
        df <- openxlsx::read.xlsx(input$upload_excel$datapath, sheet = sheet)
        gene_groups_map_temp <- list()
        
        for (colname in colnames(df)) {
          col_genes <- as.character(df[[colname]])
          col_genes <- trimws(col_genes)
          col_genes <- col_genes[!is.na(col_genes) & col_genes != ""]
          
          dup_within_col <- col_genes[duplicated(col_genes)]
          if (length(dup_within_col) > 0) {
            duplicates_within[[paste0(sheet, "_", colname)]] <- list(
              sheet = sheet,
              column = colname,
              duplicate_markers = unique(dup_within_col)
            )
          }
          
          for (gene in unique(col_genes)) {
            if (!is.null(gene_groups_map_temp[[gene]])) {
              gene_groups_map_temp[[gene]] <- c(gene_groups_map_temp[[gene]], colname)
            } else {
              gene_groups_map_temp[[gene]] <- colname
            }
          }
          
          missing <- setdiff(col_genes, seurat_genes)
          if (length(missing) > 0) {
            missing_list[[paste0(sheet, "_", colname)]] <- list(
              sheet = sheet,
              column = colname,
              missing_markers = missing
            )
          }
        }
        
        dup_across_cols <- Filter(function(g) length(g) > 1, gene_groups_map_temp)
        if (length(dup_across_cols) > 0) {
          for (gene in names(dup_across_cols)) {
            duplicates_across[[paste0(sheet, "_", gene)]] <- list(
              sheet = sheet,
              gene = gene,
              columns = dup_across_cols[[gene]]
            )
          }
        }
      }
      
      error_msg <- ""
      if (length(duplicates_within) > 0) {
        for (d in duplicates_within) {
          error_msg <- paste0(
            error_msg,
            "<b>Sheet:</b> ", d$sheet, "<br>",
            "<b>Column:</b> ", d$column, "<br>",
            "<b>Duplicate markers in column:</b> ", paste(d$duplicate_markers, collapse = ", "), 
            "<br><br>"
          )
        }
      }
      
      if (length(duplicates_across) > 0) {
        for (d in duplicates_across) {
          error_msg <- paste0(
            error_msg,
            "<b>Sheet:</b> ", d$sheet, "<br>",
            "<b>Marker:</b> ", d$gene, "<br>",
            "<b>Appears in multiple columns:</b> ", paste(d$columns, collapse = ", "),
            "<br><br>"
          )
        }
      }
      
      if (error_msg != "") {
        file.remove(input$upload_excel$datapath)
        uploaded_filename(NULL)
        showModal(modalDialog(
          title = "🚫 Error: Duplicate markers detected",
          HTML(error_msg),
          easyClose = FALSE,
          footer = tagList(
            modalButton("Close"),
            downloadButton("download_error_msg", "Download Error", class = "btn btn-success")
          )
        ))
        
        output$download_error_msg <- downloadHandler(
          filename = function() {
            paste0("duplicate_markers_error_", Sys.Date(), ".txt")
          },
          content = function(file) {
            txt <- gsub("<br>", "\n", error_msg)
            txt <- gsub("<.*?>", "", txt) 
            writeLines(txt, file)
          }
        )
        
        return(NULL)
      }
      
      if (length(missing_list) > 0) {
        warning_msg <- ""
        for (item in missing_list) {
          warning_msg <- paste0(
            warning_msg,
            "<b>Sheet:</b> ", item$sheet, "<br>",
            "<b>Column:</b> ", item$column, "<br>",
            "<b>Missing markers:</b> ", paste(item$missing_markers, collapse = ", "), 
            "<br><br>"
          )
        }
        showModal(modalDialog(
          title = "⚠️ Warning: Some markers are not found in Seurat object",
          HTML(warning_msg),
          easyClose = TRUE
        ))
      }
      file_info <- list(
        name = input$upload_excel$name,
        datapath = input$upload_excel$datapath,
        size = input$upload_excel$size,
        type = input$upload_excel$type
      )
      uploaded_file_info(file_info)
      showNotification("File uploaded successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      uploaded_filename(NULL)
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
  
  # SHOW SEARCH RESULT BUTTON UI AND MODAL
  output$show_search_result_ui <- renderUI({
    req(resolution_search_results())
    div(style = "display: inline-block; vertical-align: top; margin-left:10px;",
        actionButton("show_search_result_btn", "Show search result", 
                     class = "btn-info", width = "150px"))
  })
  
  # SHOW SEARCH RESULT MODAL
  observeEvent(input$show_search_result_btn, {
    req(resolution_search_results())
    showModal(modalDialog(
      title = "Resolution Search Results",
      renderTable(resolution_search_results(), align = 'c', bordered = TRUE, striped = TRUE, hover = TRUE,
                  width = 'auto', rownames = FALSE),
      easyClose = TRUE
    ))
  })
  
  # UPDATE UMAP BUTTON
  observeEvent(input$update_umap_btn_ui, {
    req(seuratObj(), status$clustered, input$umap_dist)
    
    if (input$umap_dist >= 0.001 && input$umap_dist <= 1.0) {
      showModal(modalDialog("Updating UMAP with new min.dist...", footer = NULL, easyClose = FALSE))
      
      tryCatch({
        base_name <- if (!is.null(subset_data_path())) {
          paste0("Subset: ", file_path_sans_ext(basename(subset_data_path())))
        } else {
          file_path_sans_ext(basename(input$data_path))
        }
        
        title <- paste0(
          base_name,
          ", resolution: ", input$resolution,
          ", min_dist: ", input$umap_dist %||% min_dist_default(),
          if ("harmony" %in% names(seuratObj()@reductions)) ", Harmony" else ""
        )
        
        plot_title(title)
        if("harmony" %in% names(seuratObj()@reductions)){
          seuratObj(RunUMAP(seuratObj(), 
                            dims = 1:input$nn_dims, 
                            min.dist = input$umap_dist, 
                            reduction = "harmony",
                            verbose = FALSE))
        } else {
          seuratObj(RunUMAP(seuratObj(), 
                            dims = 1:input$nn_dims, 
                            min.dist = input$umap_dist, 
                            reduction = "pca",
                            verbose = FALSE))
        }
        umap_plot(DimPlot(seuratObj(), reduction = "umap", label = FALSE) + 
                    ggtitle("") +
                    theme(plot.title = element_text(face = "bold", size = 10),
                          legend.text = element_text(size = 10)))
        removeModal()
        showNotification("UMAP updated with new min.dist!", type = "message")
      }, error = function(e) {
        removeModal()
        showModal(modalDialog("UMAP update failed:", e$message, easyClose = TRUE))
      })
    }
  })
  
  # CLUSTER PLOTS UI
  output$cluster_umap_ui <- renderPlot({
    req(!is.null(umap_plot()))
    umap_plot()
  })
  
  # CLUSTER SUMMARY TABLE UI
  output$cluster_table_ui <- renderPlot({
    req(seuratObj(), table_data())
    
    df <- table_data()
    tg <- gridExtra::tableGrob(
      df,
      rows = NULL,
      theme = gridExtra::ttheme_minimal(
        base_size = 10,
        padding = grid::unit(c(3, 3), "mm"),
        core = list(
          fg_params = list(hjust = 0.5, x = 0.5)
        ),
        colhead = list(
          fg_params = list(fontface = "bold", hjust = 0.5, x = 0.5)
        )
      )
    )
    grid::grid.newpage()
    grid::grid.draw(tg)
  })
  
  # CLUSTER HEATMAP UI
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
    height_heatmap <- if(length(unique(heatmap_plot()$data$cluster)) <= 10){400}else{length(unique(heatmap_plot()$data$cluster))*40}

    tagList(
      div(style = "display: inline-block; vertical-align: top; margin-left:30px; font-weight: bold;",fluidRow(textOutput("plot_title"))),
      fluidRow(
        column(4,
               plotOutput("cluster_umap_ui", height = height_heatmap)),
        column(3,
               plotOutput("cluster_table_ui", height = height_heatmap)),
        column(5,
               plotOutput("cluster_heatmap_ui", height = height_heatmap))
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
  get_combined_data <- function(df_table, df_heat, df_dot1 = NULL, df_dot2 = NULL, file_path = "cluster_analysis_results.xlsx"){
    
    wb <- createWorkbook()
    
    addWorksheet(wb, "table_data")
    writeData(wb, "table_data", df_table)
    
    addWorksheet(wb, "heatmap_data")
    
    heat_wide <- df_heat %>%
      mutate(cluster = as.character(cluster)) %>%
      pivot_wider(
        id_cols = cluster,
        names_from = celltype,
        values_from = score
      ) %>%
      arrange(cluster)
    
    heat_mat <- as.data.frame(heat_wide)
    rownames(heat_mat) <- heat_mat$cluster
    heat_mat$cluster <- NULL
    
    writeData(wb, "heatmap_data", heat_mat, rowNames = TRUE)
    
    num_style <- createStyle(numFmt = "0.00")
    addStyle(
      wb, "heatmap_data", num_style,
      rows = 2:(nrow(heat_mat)+1),
      cols = 2:(ncol(heat_mat)+1),
      gridExpand = TRUE
    )
    
    write_dotplot_sheet <- function(wb, df, sheet_name, value_col){
      
      # -----------------------------
      # 确保 cluster 是数字排序
      # -----------------------------
      df <- df %>%
        mutate(cluster = as.numeric(as.character(cluster))) %>%
        arrange(cluster, group, gene)
      
      df_wide <- df %>%
        select(cluster, gene, group, !!sym(value_col)) %>%
        pivot_wider(
          id_cols = cluster,
          names_from = gene,
          values_from = !!sym(value_col)
        )
      
      mat <- as.data.frame(df_wide)
      rownames(mat) <- mat$cluster
      mat$cluster <- NULL
      
      gene_order <- colnames(mat)
      group_map <- df %>%
        distinct(gene, group) %>%
        arrange(match(gene, gene_order))
      
      group_row <- as.character(group_map$group)
      
      addWorksheet(wb, sheet_name)
      
      writeData(wb, sheet_name,
                t(group_row),
                startRow = 1,
                startCol = 2,
                colNames = FALSE)
      
      writeData(wb, sheet_name,
                t(gene_order),
                startRow = 2,
                startCol = 2,
                colNames = FALSE)
      
      writeData(wb, sheet_name,
                mat,
                startRow = 3,
                startCol = 1,
                rowNames = TRUE,
                colNames = FALSE)
      
      writeData(wb, sheet_name,
                "cluster",
                startRow = 2,
                startCol = 1,
                colNames = FALSE)
      
      rle_group <- rle(group_row)
      start_col <- 2
      center_style <- createStyle(
        halign = "center",
        valign = "center",
        textDecoration = "bold"
      )
      
      for(i in seq_along(rle_group$lengths)){
        len <- rle_group$lengths[i]
        end_col <- start_col + len - 1
        
        if(len > 1){
          mergeCells(wb, sheet = sheet_name,
                     cols = start_col:end_col,
                     rows = 1)
        }
        addStyle(
          wb,
          sheet = sheet_name,
          style = center_style,
          rows = 1,
          cols = start_col:end_col,
          gridExpand = TRUE,
          stack = TRUE
        )
        
        start_col <- start_col + len
      }
      
      num_style <- createStyle(numFmt = "0.00")
      addStyle(
        wb, sheet_name, num_style,
        rows = 3:(nrow(mat)+2),
        cols = 2:(ncol(mat)+1),
        gridExpand = TRUE
      )
    }
    
    if(!is.null(df_dot1)){
      write_dotplot_sheet(wb, df_dot1, "dotplot1_avg_exp", "avg_exp")
      write_dotplot_sheet(wb, df_dot1, "dotplot1_pct_exp", "pct_exp")
    }
    if(!is.null(df_dot2)){
      write_dotplot_sheet(wb, df_dot2, "dotplot2_avg_exp", "avg_exp")
      write_dotplot_sheet(wb, df_dot2, "dotplot2_pct_exp", "pct_exp")
    }
    
    saveWorkbook(wb, file_path, overwrite = TRUE)
    
  }
  # DOWNLOAD CLUSTER PLOTS
  output$download_cluster_plots <- downloadHandler(
    filename = function() {
      paste("Cluster_Plots_", file_path_sans_ext(basename(input$data_path)), Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      temp_dir <- tempdir()
      umap_file <- file.path(temp_dir, paste0("Cluster_UMAP_",input$title_for_download,".pdf"))
      heatmap_file <- file.path(temp_dir, paste0("Cluster_Heatmap_",input$title_for_download,".pdf"))
      dotplot1_file <- file.path(temp_dir, paste0("Cluster_Dotplot1_",input$title_for_download,".pdf"))
      dotplot2_file <- file.path(temp_dir, paste0("Cluster_Dotplot2_",input$title_for_download,".pdf"))
      combined_file <- file.path(temp_dir, paste0("Cluster_Combined_",input$title_for_download,".pdf"))
      combined_data <- file.path(temp_dir, paste0("Cluster_Data_",input$title_for_download,".xlsx"))
      combined_plot <- create_combined_cluster_plot_patchwork(
        umap_plot(),
        table_data(),
        heatmap_plot(),
        dotplot1_plot(),
        dotplot2_plot()
      )
      height_heatmap <- length(unique(heatmap_plot()$data$cluster))*0.8
      ggsave(heatmap_file, plot = heatmap_plot(), width = 7, height = height_heatmap)
      ggsave(umap_file, plot = umap_plot(), width = 7, height = 8)
      height_combined <- if(length(unique(heatmap_plot()$data$cluster)) <= 10){20}else{20+0.72*length(unique(heatmap_plot()$data$cluster))}
      ggsave(combined_file, plot = combined_plot, width = 18, height = height_combined)
      if (!is.null(dotplot1_plot())) {
        ggsave(dotplot1_file, plot = dotplot1_plot(), width = 18, height = 6)
      }
      
      if (!is.null(dotplot2_plot())) {
        ggsave(dotplot2_file, plot = dotplot2_plot(), width = 18, height = 6)
      }
      
      get_combined_data(table_data(),
                        heatmap_plot()$data,
                        if (!is.null(dotplot1_plot()$data)) dotplot1_plot()$data else NULL,
                        if (!is.null(dotplot2_plot()$data)) dotplot2_plot()$data else NULL,
                        file_path = combined_data)
      
      zip_files <- c(heatmap_file)
      if (file.exists(dotplot1_file)) {
        zip_files <- c(zip_files, dotplot1_file)
      }
      if (file.exists(dotplot2_file)) {
        zip_files <- c(zip_files, dotplot2_file)
      }
      zip_files <- c(zip_files, umap_file, combined_file, combined_data)
      zip::zip(zipfile = file, files = zip_files, mode = "cherry-pick")
    }
  )
  
  output$sample_subset_rename_ui <- renderUI({
    req(seuratObj(),input$sample_select)
    vals <- sort(unique(seuratObj()@meta.data[[input$sample_select]]))

    tagList(
      lapply(vals, function(val) {
        fluidRow(
          column(4,
                 checkboxInput(paste0("sample_cb_", val), label = val, value = TRUE)
          ),
          column(8,
                 textInput(paste0("sample_rename_", val), label = NULL, value = val,
                           placeholder = val)
          )
        )
      }))
  })
  output$cluster_subset_rename_ui <- renderUI({
    req(seuratObj())
    vals <- sort(unique(seuratObj()@meta.data[[current_ident()]]))
    
    tagList(
      lapply(vals, function(val) {
        fluidRow(
          column(4,
                 checkboxInput(paste0("cluster_cb_", val), label = val, value = TRUE)
          ),
          column(8,
                 textInput(paste0("cluster_rename_", val), label = NULL, value = val,
                           placeholder = val)
          )
        )
      }))
  })
  
  
  # SUBSET BUTTON
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
  
  # KEEP ANALYSIS BUTTON
  observeEvent(input$keep_analysis, {
    showModal(
      modalDialog(
        title = "Confirm Analyzing on subset data",
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
    showModal(modalDialog(
      title = "Processing Subset",
      "Applying subset and rename operations...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      get_selected_values <- function(obj, col, prefix) {
        vals <- as.character(unique(obj@meta.data[[col]]))
        selected <- sapply(vals, function(val) {
          if(isTRUE(input[[paste0(prefix, "_", val)]])) val else NA
        })
        na.omit(selected)
      }
      apply_renaming <- function(obj, col, vals, prefix) {
        
        original_levels <- levels(factor(obj@meta.data[[col]]))
        obj@meta.data[[col]] <- as.character(obj@meta.data[[col]])
        for (val in vals) {
          new_name <- input[[paste0(prefix, "_", val)]]
          if (!is.null(new_name) && nzchar(new_name)) {
            obj@meta.data[[col]][obj@meta.data[[col]] == val] <- new_name
          }
        }
        new_levels <- sapply(original_levels, function(val) {
          new_name <- input[[paste0(prefix, "_", val)]]
          if (!is.null(new_name) && nzchar(new_name)) {
            new_name
          } else {
            val
          }
        })
        new_levels <- new_levels[new_levels %in% unique(obj@meta.data[[col]])]
        
        obj@meta.data[[col]] <- factor(
          obj@meta.data[[col]],
          levels = new_levels
        )
        
        return(obj)
      }
      
      # Get selected values with validation
      sample_vals <- get_selected_values(seuratObj(), input$sample_select, "sample_cb")
      cluster_vals <- get_selected_values(seuratObj(), current_ident(), "cluster_cb")
      
      # Validate at least one sample and cluster is selected
      validate(
        need(length(sample_vals) > 0, "Please select at least one sample"),
        need(length(cluster_vals) > 0, "Please select at least one cluster")
      )
      
      # Get intersecting cells with validation
      sample_cells <- colnames(seuratObj())[seuratObj()@meta.data[[input$sample_select]] %in% sample_vals]
      cluster_cells <- colnames(seuratObj())[seuratObj()@meta.data[[current_ident()]] %in% cluster_vals]
      keep_cells <- intersect(sample_cells, cluster_cells)

      validate(
        need(length(keep_cells) > 0, 
             "No cells match the selected criteria. Please adjust your subsetting options.")
      )
      # Create subset
      obj_sub <- subset(seuratObj(), cells = keep_cells)
      
      # Apply renaming
      obj_sub <- apply_renaming(obj_sub, input$sample_select, sample_vals, "sample_rename")
      obj_sub <- apply_renaming(obj_sub, current_ident(), cluster_vals, "cluster_rename")
      temp_seuratObj(obj_sub)
      removeModal()
      showNotification("Subset and renaming completed successfully!", type = "message")
      show_download_print(FALSE)
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$cluster_umap_plot <- renderPlot({
    req(seuratObj(), current_ident())
    sobj <- if(!is.null(temp_seuratObj())){temp_seuratObj()}else{seuratObj()}
    DimPlot(sobj, reduction = "umap", label = FALSE, group.by = current_ident()) + 
      ggtitle("") +
      theme(
        plot.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 8)   
      )
  })
  
  output$sample_umap_plot <- renderPlot({
    req(seuratObj(), input$sample_select)
    sobj <- if(!is.null(temp_seuratObj())){temp_seuratObj()}else{seuratObj()}
    DimPlot(sobj, reduction = "umap", label = FALSE, group.by = input$sample_select) + 
      ggtitle("") +
      theme(plot.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 10),
            axis.text = element_text(size = 8))
  })
  
  observeEvent(input$sample_ident_analysis, {
    req(seuratObj())
    req(input$sample_select)
    sobj <- if(!is.null(temp_seuratObj())){temp_seuratObj()}else{seuratObj()}
    showModal(modalDialog(
      title = "Processing Analyzing",
      "Please wait...",
      footer = NULL,
      easyClose = FALSE
    ))
    #MARKER
    get_selected_values <- function(obj, col, checkbox_prefix = "sample_cb", rename_prefix = "sample_rename") {
      vals <- as.character(unique(obj@meta.data[[col]]))
      mapping <- sapply(vals, function(val) {
        checked <- isTRUE(input[[paste0(checkbox_prefix, "_", val)]])
        if (!checked) return(NULL)
        new_name <- input[[paste0(rename_prefix, "_", val)]]
        if (is.null(new_name) || new_name == "") val else new_name
      }, simplify = FALSE)
      unlist(mapping)
    }
    
    sample_vals <- get_selected_values(seuratObj(), input$sample_select)
    if (length(sample_vals) == 0) {
      removeModal()
      showModal(modalDialog(
        title = "Warning",
        "Please select at least one sample.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return(NULL)
    }
    pdf_dir <- file.path(tempdir(), "sample_analysis_pdfs")
    dir.create(pdf_dir, showWarnings = FALSE)
    
    pdf_paths <- list()
    path <- if (!is.null(input$upload_excel)) {
      input$upload_excel$datapath
    } else {
      file.path(getwd(), "Gene_Markers.xlsx")
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
      biomarkers_dotplot2 <- NULL
    })
    
    generate_heatmap_data <- function(sobj, biomarkers = NULL) {
      if (is.null(biomarkers)) {
        stop("parameter biomarkers should be specified (now NULL)")
      }
      
      DefaultAssay(sobj) <- "RNA"
      
      clusters <- levels(Idents(sobj))
      all_genes <- unique(unlist(biomarkers))
      
      fc_all <- purrr::map_dfr(clusters, function(clust) {
        fc <- tryCatch({
          FoldChange(sobj, ident.1 = clust, features = all_genes, assay = "RNA", slot = "data")
        }, error = function(e) {
          data.frame(avg_log2FC = rep(NA_real_, length(all_genes)), row.names = all_genes)
        })
        
        fc_df <- as.data.frame(fc)
        fc_df$gene <- rownames(fc_df)
        fc_df$cluster <- clust
        fc_df[, c("cluster", "gene", "avg_log2FC")]
      })
      
      celltypes_ordered <- names(biomarkers)
      
      scores <- purrr::map_dfr(seq_along(biomarkers), function(i) {
        celltype <- celltypes_ordered[i]
        genes <- biomarkers[[i]]
        tmp <- fc_all %>%
          filter(gene %in% genes) %>%
          group_by(cluster) %>%
          summarise(score = mean(avg_log2FC, na.rm = TRUE), .groups = "drop") %>%
          mutate(celltype = celltype)
        tmp
      })
      
      scores <- scores %>%
        mutate(
          celltype = factor(celltype, levels = celltypes_ordered),
          cluster = factor(cluster, levels = clusters)
        )
      
      scores <- tidyr::expand(scores, cluster, celltype) %>%
        left_join(scores, by = c("cluster", "celltype")) %>%
        mutate(score = coalesce(score, 0)) %>%
        arrange(cluster, celltype)
      return(scores)
    }
    sample_df <- FetchData(sobj, vars = input$sample_select)
    
    all_sample_plots <- list()
    all_sample_tables <- list()
    
    for (sample in sample_vals) {
      removeModal()
      showModal(modalDialog(
        title = "Analyzing",
        paste0("Analysis on sample: ", sample, ". Please wait..."),
        footer = NULL,
        easyClose = FALSE
      ))
      cells_in_sample <- rownames(sample_df)[as.character(sample_df[[input$sample_select]]) == sample]
      sobj_subset <- subset(sobj, cells = cells_in_sample)
      Idents(sobj_subset) <- "seurat_clusters"
      
      heatmap_df <- generate_heatmap_data(sobj_subset, biomarkers_heatmap)
      heatmap <- CreateCellTypesHeatmap(heatmap_df,
                                        labels_text_size = 4,
                                        xaxis_text_size = 10,
                                        yaxis_text_size = 10,
                                        rotate_x = TRUE)
      biomarkers_dotplot <- c(biomarkers_dotplot1, biomarkers_dotplot2)
      dp1 <- if (!is.null(biomarkers_dotplot1)) custom_dotplot(sobj_subset, biomarkers_dotplot1, biomarkers_dotplot) else NULL
      dp2 <- if (!is.null(biomarkers_dotplot2)) custom_dotplot(sobj_subset, biomarkers_dotplot2, biomarkers_dotplot) else NULL
      
      umap <- DimPlot(sobj_subset, reduction = "umap", label = FALSE, group.by = current_ident()) +
        ggtitle("") + 
        theme(plot.title = element_text(face = "bold", size = 10),
              legend.text = element_text(size = 10)) +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
      
      df <- get_cluster_summary_table(sobj_subset)
      colnames(df) <- c("cluster", "ncell", "pct", "avg.counts", "avg.genes")
      combined_plot <- create_combined_cluster_plot_patchwork(
        umap, df, heatmap, dp1, dp2, title = paste0(sample, " from ", input$title_for_download)
      )
      
      all_sample_plots[[sample]] <- combined_plot
      all_sample_tables[[sample]] <- list(
        summary = df,
        heatmap = heatmap$data,
        dotplot1 = dp1,
        dotplot2 = dp2
      )
      
      excel_file <- file.path(pdf_dir, paste0(sample, "_cluster_data.xlsx"))
      
      get_combined_data(
        df_table = df,
        df_heat  = heatmap$data,
        df_dot1  = if (!is.null(dp1)) dp1$data else NULL,
        df_dot2  = if (!is.null(dp2)) dp2$data else NULL,
        file_path = excel_file)
      
      height_combined <- if(length(unique(heatmap$data$cluster)) <= 10){20}else{20+0.72*length(unique(heatmap$data$cluster))}
      
      pdf_file <- file.path(pdf_dir, paste0(sample, "_Clustering_combined.pdf"))
      ggsave(pdf_file, plot = combined_plot, width = 18, height = height_combined)
      pdf_paths[[sample]] <- list(
        pdf   = pdf_file,
        excel = excel_file
      )
    }
    
    if (length(all_sample_plots) > 0) {
      combined_pdf_file <- file.path(pdf_dir, "All_Samples_Combined.pdf")
      pdf(combined_pdf_file, width = 18, height = 20) 
      for (sample in names(all_sample_plots)) {
        print(all_sample_plots[[sample]])
      }
      dev.off()
      
      combined_excel_file <- file.path(pdf_dir, "All_Samples_Data.xlsx")
      showModal(modalDialog(
        title = "Analyzing",
        paste0("Generate combiend file..."),
        footer = NULL,
        easyClose = FALSE
      ))
      # Create a combined Excel workbook for multiple samples
      wb <- openxlsx::createWorkbook()
      
      # Center style for headers and merged cells
      center_style <- createStyle(
        halign = "center",
        valign = "center",
        textDecoration = "bold"
      )
      
      # Number style with 2 decimal places
      num_style <- createStyle(numFmt = "0.00")
      
      # -------------------------
      # 1. Write table_data sheet (summary per sample)
      # -------------------------
      openxlsx::addWorksheet(wb, "table_data")
      current_row <- 1
      for (sample in names(all_sample_tables)) {
        # Write sample header
        openxlsx::writeData(
          wb, "table_data", 
          x = data.frame(Sample = paste0("=== ", sample, " ==="), stringsAsFactors = FALSE),
          startCol = 1, startRow = current_row
        )
        current_row <- current_row + 2
        
        # Write summary table and apply number formatting
        summary_df <- all_sample_tables[[sample]]$summary
        openxlsx::writeData(
          wb, "table_data", 
          x = summary_df,
          startCol = 1, startRow = current_row, colNames = TRUE
        )
        
        # Apply number formatting to numeric columns in summary
        if (!is.null(summary_df) && nrow(summary_df) > 0) {
          # Find numeric columns (skip first column if it's cluster/celltype)
          num_cols <- which(sapply(summary_df, is.numeric))
          if (length(num_cols) > 0) {
            addStyle(
              wb, "table_data", num_style,
              rows = (current_row + 1):(current_row + nrow(summary_df)),
              cols = num_cols,
              gridExpand = TRUE
            )
          }
        }
        
        current_row <- current_row + nrow(summary_df) + 2
      }
      
      # -------------------------
      # 2. Write heatmap_data sheet
      # -------------------------
      openxlsx::addWorksheet(wb, "heatmap_data")
      current_row <- 1
      
      for (sample in names(all_sample_tables)) {
        # Sample header
        openxlsx::writeData(
          wb, "heatmap_data", 
          x = data.frame(Sample = paste0("=== ", sample, " ==="), stringsAsFactors = FALSE),
          startCol = 1, startRow = current_row
        )
        current_row <- current_row + 2
        
        # Convert heatmap data to wide format with proper cluster ordering
        heat_wide <- all_sample_tables[[sample]]$heatmap %>%
          mutate(
            cluster = as.numeric(as.character(cluster)),  # ensure numeric
            score = round(score, 2)                       # round to 2 decimals
          ) %>%
          arrange(cluster) %>%  # sort by numeric cluster
          pivot_wider(
            id_cols = cluster,
            names_from = celltype,
            values_from = score
          ) %>%
          arrange(cluster)  # ensure final sorting
        
        heat_mat <- as.data.frame(heat_wide)
        rownames(heat_mat) <- heat_mat$cluster
        heat_mat$cluster <- NULL
        
        # Write heatmap matrix with row names (cluster numbers)
        openxlsx::writeData(
          wb, "heatmap_data", x = heat_mat,
          startCol = 1, startRow = current_row, rowNames = TRUE
        )
        
        # Apply number formatting - fix row range calculation
        if (nrow(heat_mat) > 0 && ncol(heat_mat) > 0) {
          data_start_row <- current_row + 1  # +1 because rowNames=TRUE puts rownames in first column
          data_end_row <- current_row + nrow(heat_mat)
          data_start_col <- 2  # first column is rownames (cluster)
          data_end_col <- 1 + ncol(heat_mat)
          
          addStyle(
            wb, "heatmap_data", num_style,
            rows = data_start_row:data_end_row,
            cols = data_start_col:data_end_col,
            gridExpand = TRUE
          )
        }
        
        current_row <- current_row + nrow(heat_mat) + 2
      }
      
      # -------------------------
      # 3. Function to write combined dotplot sheets
      # -------------------------
      write_combined_dotplot_sheet <- function(wb, sheet_name, all_tables, value_col, dotplot_type = "dotplot1") {
        openxlsx::addWorksheet(wb, sheet_name)
        current_row <- 1
        
        for (sample in names(all_tables)) {
          df_dot <- all_tables[[sample]][[dotplot_type]]
          if (!is.null(df_dot) && !is.null(df_dot$data)) {
            # Process data with proper cluster ordering
            df <- df_dot$data %>%
              mutate(
                cluster = as.numeric(as.character(cluster)),  # convert to numeric for proper sorting
                !!sym(value_col) := round(!!sym(value_col), 2)
              ) %>%
              arrange(cluster, group, gene)  # sort by numeric cluster
            
            # Convert to wide format
            df_wide <- df %>%
              select(cluster, gene, group, !!sym(value_col)) %>%
              pivot_wider(
                id_cols = cluster,
                names_from = gene,
                values_from = !!sym(value_col)
              ) %>%
              arrange(cluster)  # ensure clusters are in numeric order
            
            mat <- as.data.frame(df_wide)
            rownames(mat) <- mat$cluster
            mat$cluster <- NULL
            
            # Get gene order and group mapping
            gene_order <- colnames(mat)
            group_map <- df %>%
              distinct(gene, group) %>%
              arrange(match(gene, gene_order))
            group_row <- as.character(group_map$group)
            
            # Write sample header
            openxlsx::writeData(
              wb, sheet_name, 
              x = data.frame(Sample = paste0("=== ", sample, " ==="), stringsAsFactors = FALSE),
              startCol = 1, startRow = current_row
            )
            current_row <- current_row + 1
            
            # Write group row
            openxlsx::writeData(
              wb, sheet_name, t(group_row),
              startRow = current_row, startCol = 2, colNames = FALSE
            )
            
            # Write gene row
            openxlsx::writeData(
              wb, sheet_name, t(gene_order),
              startRow = current_row + 1, startCol = 2, colNames = FALSE
            )
            
            # Write cluster header
            openxlsx::writeData(
              wb, sheet_name, "cluster",
              startRow = current_row + 1, startCol = 1, colNames = FALSE
            )
            
            # Write matrix
            openxlsx::writeData(
              wb, sheet_name, mat,
              startRow = current_row + 2, startCol = 1,
              rowNames = TRUE, colNames = FALSE
            )
            
            # Merge group header cells
            rle_group <- rle(group_row)
            start_col <- 2
            for(i in seq_along(rle_group$lengths)) {
              len <- rle_group$lengths[i]
              end_col <- start_col + len - 1
              
              if(len > 1) {
                mergeCells(wb, sheet = sheet_name,
                           cols = start_col:end_col,
                           rows = current_row)
              }
              addStyle(wb, sheet = sheet_name, style = center_style,
                       rows = current_row, cols = start_col:end_col,
                       gridExpand = TRUE, stack = TRUE)
              
              start_col <- start_col + len
            }
            
            # Apply number formatting - FIXED row range calculation
            if (nrow(mat) > 0 && ncol(mat) > 0) {
              data_start_row <- current_row + 2
              data_end_row <- current_row + 1 + nrow(mat)
              data_start_col <- 2
              data_end_col <- 1 + ncol(mat)
              
              addStyle(
                wb, sheet_name, num_style,
                rows = data_start_row:data_end_row,
                cols = data_start_col:data_end_col,
                gridExpand = TRUE
              )
            }
            
            # Update row for next sample
            current_row <- current_row + 3 + nrow(mat) + 1
          }
        }
      }
      
      # -------------------------
      # 4. Write dotplot sheets if available
      # -------------------------
      has_dotplot1 <- any(sapply(all_sample_tables, function(x) !is.null(x$dotplot1)))
      if (has_dotplot1) {
        write_combined_dotplot_sheet(wb, "dotplot1_avg_exp", all_sample_tables, "avg_exp", "dotplot1")
        write_combined_dotplot_sheet(wb, "dotplot1_pct_exp", all_sample_tables, "pct_exp", "dotplot1")
      }
      
      has_dotplot2 <- any(sapply(all_sample_tables, function(x) !is.null(x$dotplot2)))
      if (has_dotplot2) {
        write_combined_dotplot_sheet(wb, "dotplot2_avg_exp", all_sample_tables, "avg_exp", "dotplot2")
        write_combined_dotplot_sheet(wb, "dotplot2_pct_exp", all_sample_tables, "pct_exp", "dotplot2")
      }
      
      # -------------------------
      # 5. Save workbook
      # -------------------------
      openxlsx::saveWorkbook(wb, combined_excel_file, overwrite = TRUE)
      pdf_paths[["All"]] <- list(
        pdf = combined_pdf_file,
        excel = combined_excel_file
      )
    }
    
    show_download_print(TRUE)
    sample_ident_analysis_path$paths <- pdf_paths
    removeModal()
    showNotification(paste("Generated PDFs for", length(sample_vals), "samples (paths stored)."), type = "message")
    
  })
  
  output$download_ident_sample_analysis_combined <- downloadHandler(
    filename = function() {
      paste0("combined_sample_analysis_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(sample_ident_analysis_path$paths)
      req(sample_ident_analysis_path$paths[["All"]]) 
      
      tmp_dir <- tempdir()
      combined_pdf <- sample_ident_analysis_path$paths[["All"]]$pdf
      combined_excel <- sample_ident_analysis_path$paths[["All"]]$excel
      
      pdf_dest <- file.path(tmp_dir, basename(combined_pdf))
      excel_dest <- file.path(tmp_dir, basename(combined_excel))
      
      file.copy(combined_pdf, pdf_dest, overwrite = TRUE)
      file.copy(combined_excel, excel_dest, overwrite = TRUE)
      
      zip::zipr(zipfile = file, files = c(pdf_dest, excel_dest), recurse = FALSE)
    },
    contentType = "application/zip"
  )
  
  output$download_ident_sample_analysis_separate <- downloadHandler(
    filename = function() {
      paste0("separate_sample_analysis_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(sample_ident_analysis_path$paths)
      
      tmp_dir <- tempdir()
      all_files <- c()
      
      for (sample in names(sample_ident_analysis_path$paths)) {
        if (sample == "All") next  
        pdf_src   <- sample_ident_analysis_path$paths[[sample]]$pdf
        excel_src <- sample_ident_analysis_path$paths[[sample]]$excel
        
        if (file.exists(pdf_src) && file.exists(excel_src)) {
          pdf_dest <- file.path(tmp_dir, basename(pdf_src))
          excel_dest <- file.path(tmp_dir, basename(excel_src))
          
          file.copy(pdf_src, pdf_dest, overwrite = TRUE)
          file.copy(excel_src, excel_dest, overwrite = TRUE)
          
          all_files <- c(all_files, pdf_dest, excel_dest)
        }
      }
      
      if (length(all_files) == 0) {
        showNotification("No separate sample files found to download.", type = "warning")
        return(NULL)
      }
      
      zip::zipr(zipfile = file, files = all_files, recurse = FALSE)
    },
    contentType = "application/zip"
  )
  observe({
    req(sample_ident_analysis_path$paths)
    
    sample_names <- names(sample_ident_analysis_path$paths)
    
    for (sample in sample_names) {
      safe_sample_id <- URLencode(sample, reserved = TRUE)
      
      output_id <- paste0("pdf_preview_", safe_sample_id)
      
      local({
        current_sample <- sample
        current_output_id <- paste0("pdf_preview_", URLencode(current_sample, reserved = TRUE))
        
        output[[current_output_id]] <- renderUI({
          pdf_path <- sample_ident_analysis_path$paths[[current_sample]]$pdf
          
          if (!is.null(pdf_path) && file.exists(pdf_path)) {

            tryCatch({
              base64_pdf <- base64enc::base64encode(pdf_path)
              tags$iframe(
                src = paste0("data:application/pdf;base64,", base64_pdf),
                style = "width: 100%; height: 1000px; border: 1px solid #ddd;"
              )
            }, error = function(e) {
              div(
                class = "alert alert-danger",
                icon("exclamation-triangle"),
                paste("Error loading PDF:", e$message)
              )
            })
          } else {
            div(
              class = "alert alert-warning",
              icon("exclamation-triangle"),
              paste("PDF file not found or inaccessible:", current_sample)
            )
          }
        })
      })
    }
  })
  
  observeEvent(input$print_ident_sample_analysis, {
    
    req(sample_ident_analysis_path$paths)
    
    sample_names <- names(sample_ident_analysis_path$paths)
    if ("All" %in% sample_names) {
      sample_names <- c("All", sample_names[sample_names != "All"])
    }
    
    showModal(
      modalDialog(
        title = "Sample Clustering Result",
        
        do.call(tabsetPanel, 
                c(list(id = "pdf_tabs", type = "pills"),
                  lapply(sample_names, function(sample) {
                    tabPanel(
                      title = sample,
                      value = sample,
                      div(style = "height: 1000px; margin-top: 20px;",
                          uiOutput(paste0("pdf_preview_", gsub("[^[:alnum:]]", "", sample))))
                    )
                  })
                )
        ),
        
        easyClose = TRUE,
        size = "l",
        footer = modalButton("Close")
      )
    )
  })
  # CONFIRM KEEP ANALYSIS BUTTON
  observeEvent(input$confirm_keep_analysis, {
    req(seuratObj())
    tryCatch({
      if(!is.null(temp_seuratObj())){
        seuratObj(temp_seuratObj())
      }
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
      table_data(NULL)
      dotplot1_plot(NULL)
      dotplot2_plot(NULL)
      current_ident(NULL)
      plot_title(NULL)
      harmony_state(NULL)
      umap_dist_state <- reactiveVal(NULL)
      harmony_meta_state(NULL)
      user_nn_dims(NULL)
      uploaded_filename(NULL)
      temp_seuratObj(NULL)
      combined_plot(NULL)
      show_download_print(FALSE)
      shinyjs::runjs("location.reload();")
      removeModal()
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
}  

shinyApp(ui = ui, server = server)
