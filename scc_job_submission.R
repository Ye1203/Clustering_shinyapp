submit_scc_job <- function(seuratObj,
                           save_path,
                           project_name,
                           start_resolution,
                           step_resolution,
                           end_resolution,
                           integration_harmony = FALSE,
                           nn_dims = 30,
                           runtime = 12,
                           cores = 8,
                           email = NA,
                           marker_path = NULL) {
  
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  saveRDS(seuratObj, file = file.path(save_path, "input_seurat_file.rds"))
  
  file.copy(from = marker_path, to = file.path(save_path, "Gene_Markers.xlsx"), overwrite = TRUE)
  
  if (!file.exists(file.path(save_path, "input_seurat_file.rds"))) {
    stop(paste("Failed to create file:", file.path(save_path, "input_seurat_file.rds")))
  }
  
  vis_dir <- file.path(save_path, "Visualization_results")
  if (!dir.exists(vis_dir)) {
    dir.create(vis_dir, recursive = TRUE)
  }
  
  r_script_path <- file.path(save_path, "scc_job_script.R")
  
  r_code <- paste0('
dyn.load("/projectnb/wax-es/00_shinyapp/Clustering/conda_env/lib/libicui18n.so.75")
# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(cowplot)
library(future.apply)
library(gridExtra)
library(tidyr)
library(patchwork)

# Configuration parameters
output_file_path <- "', save_path, '"
email_address <- "', ifelse(is.na(email), "", email), '"
start_resolution <- ', start_resolution, '
step_resolution <- ', step_resolution, '
end_resolution <- ', end_resolution, '
cores <- ', cores, '
integration_harmony <- ', ifelse(integration_harmony, "TRUE", "FALSE"), '
nn_dims <- ', nn_dims, '

# Email function
send_email <- function(to, subject, body) {
  if (!is.na(to) && nzchar(to)) {
    f <- tempfile(fileext = ".txt")
    mail_content <- c(
      paste0("To: ", to),
      paste0("Subject: ", subject),
      "",
      body
    )
    writeLines(mail_content, f)
    cmd <- paste0("/usr/sbin/sendmail -t < ", shQuote(f))
    ret <- system(cmd)
    unlink(f)
    if (ret != 0) {
      warning("sendmail command failed with exit code ", ret)
    }
  }
}

tryCatch({
  # Load data
  seurat_obj <- readRDS(file.path(output_file_path, "input_seurat_file.rds"))
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Set up parallel processing
  plan(multithreaded, workers = cores)
  options(future.globals.maxSize = 10 * 1024^3)
  set.seed(001203)
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
    
    scores <- purrr::map_dfr(names(biomarkers), function(celltype) {
      genes <- biomarkers[[celltype]]
      tmp <- fc_all %>%
        filter(gene %in% genes) %>%
        group_by(cluster) %>%
        summarise(score = mean(avg_log2FC, na.rm = TRUE)) %>%
        mutate(celltype = celltype)
      tmp
    })
    
    scores <- left_join(tidyr::expand(scores, cluster, celltype), scores, by = join_by(cluster, celltype)) %>%
      mutate(score = coalesce(score, 0))
    
    labels <- scores %>%
      group_by(cluster) %>%
      summarise(prediction = ifelse(max(score) > 0, celltype[which.max(score)], "Unknown")) %>%
      ungroup() %>%
      mutate(prediction = paste0(prediction, "(", cluster, ")"))
    
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
        stop(paste0("The \'", sheet_name, "\' sheet appears empty or has no valid columns."))
      }
      
      if (any(is.na(names(df))) || any(names(df) == "")) {
        stop(paste0("Some columns in \'", sheet_name, "\' have no names. Please check header row."))
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
        stop(paste0("No valid markers found in \'", sheet_name, "\' sheet."))
      }
      
      return(biomarkers)
      
    }, error = function(e) {
      stop(paste("Failed to convert Excel to marker list. Please check Excel format.\\nDetails:", e$message))
    })
  }
  
  # Create Heatmap Function
  CreateCellTypesHeatmap <- function(df, labels_text_size = 6, xaxis_text_size = 12, yaxis_text_size = 12, rotate_x = FALSE){
    
    clusters_numeric <- sort(as.numeric(unique(df$cluster)))
    df$cluster <- factor(df$cluster, 
                        levels = as.character(clusters_numeric),
                        ordered = TRUE)
    
    tmp <- ggplot(df, aes(x = celltype, y = as.factor(cluster))) +
      geom_tile(aes(fill = score), color= "gray50", size = 0.1) +
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
  
  # Custom Dotplot Function
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
    
    clusters_numeric <- sort(as.numeric(unique(results$cluster)))
    results$cluster <- factor(results$cluster, 
                             levels = as.character(clusters_numeric),
                             ordered = TRUE)
    
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
        colors = cols,
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
  
  # Cluster Summary Table Function
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
                  avg.counts = as.integer(round(mean(seurat_obj@meta.data$nCount_RNA))),
                  avg.genes = as.integer(round(mean(seurat_obj@meta.data$nFeature_RNA))))
    )
    
    return(result)
  }
  
  # Create Combined Cluster Plot Function
  create_combined_cluster_plot_patchwork <- function(umap_plot, table_data, heatmap_plot, dotplot1 = NULL, dotplot2 = NULL, resolution) {
    
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
    
    title_text <- paste("Resolution: ", resolution)
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
  
  # Main processing
  markers_path <- file.path(output_file_path, "Gene_Markers.xlsx")
  
  resolutions <- seq(start_resolution, end_resolution, by = step_resolution)
  
  combined_folder <- file.path(output_file_path, "Visualization_results", "Combined")
  dir.create(combined_folder, recursive = TRUE, showWarnings = FALSE)
  umap_folder <- file.path(output_file_path, "Visualization_results", "UMAPs")
  dir.create(umap_folder, recursive = TRUE, showWarnings = FALSE)
  heatmap_folder <- file.path(output_file_path, "Visualization_results", "Heatmaps")
  dir.create(heatmap_folder, recursive = TRUE, showWarnings = FALSE)
  if("DOTPLOT1" %in% openxlsx::getSheetNames(markers_path)) {
    dotplot1_folder <- file.path(output_file_path, "Visualization_results", "Dotplots1")
    dir.create(dotplot1_folder, recursive = TRUE, showWarnings = FALSE)
  }
  if("DOTPLOT2" %in% openxlsx::getSheetNames(markers_path)) {
    dotplot2_folder <- file.path(output_file_path, "Visualization_results", "Dotplots2")
    dir.create(dotplot2_folder, recursive = TRUE, showWarnings = FALSE)
  }
  
  message("Starting cluster analysis...")
  message("Resolutions to process: ", paste(resolutions, collapse = ", "))
  

  future.apply::future_lapply(resolutions, function(resolution) {
    
    umap_file <- file.path(umap_folder, paste0("Resolution: ", resolution, ".pdf"))
    heatmap_file <- file.path(heatmap_folder, paste0("Resolution: ", resolution, ".pdf"))
    dotplot1_file <- file.path(dotplot1_folder, paste0("Resolution: ", resolution, ".pdf"))
    dotplot2_file <- file.path(dotplot2_folder, paste0("Resolution: ", resolution, ".pdf"))
    combined_file <- file.path(combined_folder, paste0("Resolution: ", resolution, ".pdf"))
    
    sobj <- seurat_obj
    
    # Find clusters
    sobj <- FindClusters(sobj, resolution = resolution, verbose = FALSE)
    Idents(sobj) <- factor(Idents(sobj), 
                             levels = sort(as.numeric(levels(Idents(sobj)))))
    allowed_genes <- rownames(sobj)
    
    # Load biomarkers
    biomarkers_heatmap <- NULL
    if (file.exists(markers_path)) {
      tryCatch({
        biomarkers_heatmap <- convert_marker_excel_to_list(markers_path, sheet_name = "HEATMAP", allowed_genes)
      }, error = function(e) {
        message("Failed to read heatmap markers: ", e$message)
      })
    }
    
    biomarkers_dotplot1 <- NULL
    if (file.exists(markers_path)) {
      tryCatch({
        biomarkers_dotplot1 <- convert_marker_excel_to_list(markers_path, sheet_name = "DOTPLOT1", allowed_genes)
      }, error = function(e) {
        biomarkers_dotplot1 <- NULL
      })
    }
    
    biomarkers_dotplot2 <- NULL
    if (file.exists(markers_path)) {
      tryCatch({
        biomarkers_dotplot2 <- convert_marker_excel_to_list(markers_path, sheet_name = "DOTPLOT2", allowed_genes)
      }, error = function(e) {
        biomarkers_dotplot2 <- NULL
      })
    }
    
    dotplot1_plot <- NULL
    dotplot2_plot <- NULL
    
    if (!is.null(biomarkers_dotplot1)) {
      dotplot1_plot <- custom_dotplot(sobj, biomarkers_dotplot1)
    }
    
    if (!is.null(biomarkers_dotplot2)) {
      dotplot2_plot <- custom_dotplot(sobj, biomarkers_dotplot2)
    }
    
    findcelltypes_result <- NULL
    if (!is.null(biomarkers_heatmap)) {
      findcelltypes_result <- FindCellTypesByMarkers(sobj, biomarkers_heatmap)
    } else {
      stop("biomarkers_heatmap is NULL, cannot annotate clusters")
    }
    
    heatmap_plot <- CreateCellTypesHeatmap(findcelltypes_result$heatmap_table,
                                           labels_text_size = 4,
                                           xaxis_text_size = 10,
                                           yaxis_text_size = 10,
                                           rotate_x = TRUE)
    
    sobj <- RenameIdents(sobj, findcelltypes_result$new_labels)
    sobj$shiny_clusters <- Idents(sobj)
    Idents(sobj) <- sobj$shiny_clusters
    
    # Run UMAP
    if (integration_harmony) {
      sobj <- RunUMAP(sobj,
                      dims = 1:nn_dims,
                      min.dist = 0.3,
                      reduction = "harmony",
                      verbose = FALSE)
    } else {
      sobj <- RunUMAP(sobj,
                      dims = 1:nn_dims,
                      min.dist = 0.3,
                      reduction = "pca",
                      verbose = FALSE)
    }
    
    umap_plot <- DimPlot(sobj, reduction = "umap", label = FALSE) + 
      ggtitle("") +
      theme(plot.title = element_text(face = "bold", size = 10))+
      guides(color = guide_legend(ncol = 1,
                                  override.aes = list(size = 3)))
    
    table_data <- get_cluster_summary_table(sobj)
    
    combined_plot <- create_combined_cluster_plot_patchwork(
      umap_plot = umap_plot,
      table_data = table_data,
      heatmap_plot = heatmap_plot,
      dotplot1 = dotplot1_plot,
      dotplot2 = dotplot2_plot,
      resolution = resolution
    )
    
    # Save plots
    ggsave(heatmap_file, plot = heatmap_plot, width = 7, height = 8)
    ggsave(umap_file, plot = umap_plot, width = 7, height = 8)
    ggsave(combined_file, plot = combined_plot, width = 18, height = 20)
    
    if (!is.null(dotplot1_plot)) {
      ggsave(dotplot1_file, plot = dotplot1_plot, width = 18, height = 6)
    }
    
    if (!is.null(dotplot2_plot)) {
      ggsave(dotplot2_file, plot = dotplot2_plot, width = 18, height = 6)
    }
    
    message("Finished resolution: ", resolution)
    
    return(NULL)  
  })
  
  message("All resolutions completed successfully!")
  
  # Send success email
  if (!is.na(email_address) && nzchar(email_address)) {
    body <- c(
      "Hi,",
      "",
      "The cluster analysis procedure has been completed successfully.",
      paste0("Results are available in folder: ", output_file_path),
      "",
      "Best,",
      "Bingtian"
    )
    send_email(to = email_address, subject = "Cluster Analysis - COMPLETED", body = body)
  }
  
  return(list(success = TRUE, message = "Analysis completed successfully"))
  
}, error = function(e) {
  err_msg <- conditionMessage(e)
  message("ERROR: ", err_msg)
  
  # Send error email
  if (!is.na(email_address) && nzchar(email_address)) {
    body <- c(
      "Hi,",
      "",
      "There was an error in the cluster analysis procedure.",
      paste0("Error message: ", err_msg),
      "",
      paste0("Output folder: ", output_file_path),
      "",
      "Please check the log files for more details. If you need assistance, feel free to contact me (btye@bu.edu).",
      "",
      "Best,",
      "Bingtian"
    )
    send_email(to = email_address, subject = "Cluster Analysis - ERROR", body = body)
  }
  
})
')
  
  writeLines(r_code, con = r_script_path)
  
  # Generate qsub file
  qsub_file_path <- file.path(save_path, "launch_job.qsub")
  log_path <- file.path(save_path, "cluster_job.log")
  qsub_content <- paste0(
    "#!/bin/bash\n",
    "#$ -N clustering\n",
    "#$ -cwd\n",
    "#$ -j y\n",
    "#$ -o ", log_path, "\n",
    "#$ -pe omp ", cores, "\n",
    "#$ -l h_rt=", runtime, ":00:00\n",
    "#$ -V\n",
    "#$ -P ", project_name, "\n",
    "\n",
    "echo \"==========================================================\"\n",
    "echo \"Starting on       : $(date)\"\n",
    "echo \"Running on node   : $(hostname)\"\n",
    "echo \"Current job ID    : $JOB_ID\"\n",
    "echo \"Current job name  : $JOB_NAME\"\n",
    "echo \"==========================================================\"\n",
    "\n",
    "module load R/4.4.3\n",
    "\n",
    "Rscript ", r_script_path, "\n",
    "\n",
    "echo \"==========================================================\"\n",
    "echo \"Ending on         : $(date)\"\n",
    "echo \"==========================================================\"\n"
  )
  writeLines(qsub_content, con = qsub_file_path)
  
  qsub_output <- system(paste("qsub", qsub_file_path), intern = TRUE)
  job_id <- sub(".*?([0-9]+).*", "\\1", qsub_output)
  return(job_id)
}