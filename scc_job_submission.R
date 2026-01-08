submit_scc_job <- function(seuratObj,
                           save_path,
                           project_name,
                           start_resolution,
                           step_resolution,
                           end_resolution,
                           integration_harmony = FALSE,
                           pca_num,
                           harmony_metadata,
                           runtime = 12,
                           cores = 8,
                           email = NA,
                           marker_path = NULL) {
  if(start_resolution <= 0){start_resolution = 0.0001}
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  writeLines(capture.output(print(seuratObj@commands)), file.path(save_path, "command_info.txt"))
  
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
library(gridExtra)
library(tidyr)
library(patchwork)
library(future.apply)
library(grid)
library(harmony)

# Configuration parameters
output_file_path <- "', save_path, '"
email_address <- "', ifelse(is.na(email), NULL, email), '"
start_resolution <- ', start_resolution, '
step_resolution <- ', step_resolution, '
end_resolution <- ', end_resolution, '
pca_num <- "', pca_num, '"
harmony_metadata <- "', harmony_metadata, '"
cores <- ', cores, '
integration_harmony <- ', ifelse(integration_harmony, "TRUE", "FALSE"), '
pca_num <- as.numeric(strsplit(pca_num, ",")[[1]])
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
  set.seed(001203)
  
  check_last_scale_and_pca <- function(commands_file) {
  if (!file.exists(commands_file)) {
    return(list(
      scale_used_variable_features = NA,
      pca_npcs = NA
    ))
  }
  
  lines <- tryCatch(readLines(commands_file), error = function(e) NULL)
  if (is.null(lines)) {
    return(list(
      scale_used_variable_features = NA,
      pca_npcs = NA
    ))
  }
  
  scale_idx <- grep("^\\\\$Scale", lines)
  if (length(scale_idx) == 0) {
    scale_used_variable_features <- NA
  } else {
    last_scale_start <- scale_idx[length(scale_idx)]
    next_block <- grep("^\\\\$", lines)
    next_block <- next_block[next_block > last_scale_start]
    last_scale_end <- if (length(next_block) > 0) min(next_block) - 1 else length(lines)
    scale_block <- lines[last_scale_start:last_scale_end]
    scale_command <- scale_block[grep("^Command:", scale_block)]
    scale_used_variable_features <- if (length(scale_command) == 0) {
      NA
    } else {
      grepl("VariableFeatures", scale_command)
    }
  }
  
  pca_idx <- grep("^\\\\$RunPCA", lines)
  if (length(pca_idx) == 0) {
    pca_npcs <- NA
  } else {
    last_pca_start <- pca_idx[length(pca_idx)]
    next_block <- grep("^\\\\$", lines)
    next_block <- next_block[next_block > last_pca_start]
    last_pca_end <- if (length(next_block) > 0) min(next_block) - 1 else length(lines)
    pca_block <- lines[last_pca_start:last_pca_end]
    npcs_line <- pca_block[grep("npcs\\\\s*:", pca_block)]
    pca_npcs <- if (length(npcs_line) == 0) {
      NA
    } else {
      suppressWarnings(as.numeric(sub(".*npcs\\\\s*:\\\\s*", "", npcs_line)))
    }
  }
  
  list(
    scale_used_variable_features = scale_used_variable_features,
    pca_npcs = pca_npcs
  )
}

  ana_info  <- check_last_scale_and_pca(file.path(output_file_path, "command_info.txt"))

  message("Generate different PCA UMAP embeddings...")
  for(pca in pca_num){
  if(integration_harmony){
    seurat_obj <- RunHarmony(seurat_obj, group.by.vars = harmony_metadata, reduction.use = "pca",
                            dims.use = 1:pca, assay.use = "RNA", reduction.save = "harmony", 
                            project.dim = TRUE, verbose = FALSE)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:pca, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, min.dist = 0.3, dims = 1:pca, reduction = "harmony", verbose = FALSE)
  } else {
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:pca, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, min.dist = 0.3, dims = 1:pca, reduction = "pca", verbose = FALSE)
  }
  saveRDS(seurat_obj, file.path(output_file_path, paste0("input_seurat_file_pca_", pca, ".rds")))}
  remove(seurat_obj)
  
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
        FoldChange(sobj, ident.1 = clust, features = all_genes, assay = "RNA", layer = "data")
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
      stop(paste("Failed to convert Excel to marker list. Please check Excel format.\nDetails:", e$message))
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
  
  # Custom Dotplot Function
  custom_dotplot <- function(seurat_obj, gene_list, assay = "RNA", title = NULL, cols = c("lightgrey", "blue")) {
    
    all_genes <- unlist(gene_list, use.names = FALSE)
    DefaultAssay(seurat_obj) <- assay
    clusters <- unique(Idents(seurat_obj))
    results <- data.frame()
    
    for(cluster in clusters) {
      cluster_cells <- WhichCells(seurat_obj, idents = cluster)
      cluster_data <- GetAssayData(seurat_obj, assay = assay, layer = "data")[all_genes, cluster_cells, drop = FALSE]
      
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
    
    clusters_numeric <- sort(as.numeric(unique(results$cluster)), decreasing = TRUE)
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
          size = 12
        ),
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
  create_combined_cluster_plot_patchwork <- function(umap_plot, table_data, heatmap_plot, dotplot1 = NULL, dotplot2 = NULL, title_text) {
    
    table_plot <- function(table_data) {
      if (is.null(table_data)) return(ggplot() + theme_void())
      
      table_grob <- gridExtra::tableGrob(
        table_data,
        rows = NULL,
        theme = gridExtra::ttheme_minimal(
          base_size = 12,
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
    
    title_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = title_text, 
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
    heights <- c(0.15, 2, rep(1.5, length(bottom_plots)))
    
    combined_plot <- combined_plot + 
      plot_layout(heights = heights[1:n_parts], ncol = 1)
    
    return(combined_plot)
  }
  
  # Combined Data Frame
  process_single_result <- function(x) {
    label_df <- data.frame(
      cluster    = names(x$label_name),
      label_name = as.vector(x$label_name),
      stringsAsFactors = FALSE
    )
    
    base_df <- data.frame(
      index          = x$index,
      cluster_number = x$cluster_number,
      stringsAsFactors = FALSE
    )
    base_df <- base_df[rep(1, nrow(label_df)), , drop = FALSE]
    
    df <- cbind(base_df, label_df)
    df <- dplyr::left_join(df, x$table_data, by = "cluster")
    heatmap_wide <- x$heatmap_data %>%
      tidyr::pivot_wider(
        names_from  = celltype,
        values_from = score,
        values_fill = 0
      ) %>%
      dplyr::arrange(match(cluster, df$cluster)) %>%
      dplyr::select(-cluster)
    
    df <- cbind(df, heatmap_wide)
    
    df
  }
  
  # Main processing
  markers_path <- file.path(output_file_path, "Gene_Markers.xlsx")
  resolutions <- seq(start_resolution, end_resolution, by = step_resolution)
  all_png_files <- character()
  output_folder <- file.path(output_file_path, "Visualization_results")
  temp_png_folder <- file.path(output_folder, "temp_pngs")
  dir.create(temp_png_folder, recursive = TRUE, showWarnings = FALSE)
  
  message("Starting cluster analysis...")
  message("Resolutions to process: ", paste(resolutions, collapse = ", "))
  message("PCA to process: ", paste(pca_num, collapse = ", "))
  cores_n <- if(cores > 1) cores - 1 else 1
  plan(multisession, workers = cores_n)
  options(future.globals.maxSize = 30 * 1024^3)
  
  resolution_index_table <- data.frame(
    resolution = sort(resolutions),
    pca = rep(pca_num, each = length(resolutions)),
    stringsAsFactors = FALSE) %>%
    dplyr::arrange(resolution, pca)%>%
    dplyr::mutate(
      index_int = dplyr::row_number(),
      index = ifelse(index_int < 10,
                     sprintf("%02d", index_int),
                     as.character(index_int))
    ) %>%
    dplyr::select(index, resolution, pca)


results <- future.apply::future_lapply(seq_len(nrow(resolution_index_table)), function(i) {
  resolution <- resolution_index_table$resolution[i]
  pca        <- resolution_index_table$pca[i]
  index      <- resolution_index_table$index[i]
  dyn.load("/projectnb/wax-es/00_shinyapp/Clustering/conda_env/lib/libicui18n.so.75")
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(cowplot)
  library(gridExtra)
  library(tidyr)
  library(patchwork)
  message("Processing resolution: ", resolution)

  sobj <- readRDS(file.path(output_file_path, paste0("input_seurat_file_pca_", pca, ".rds")))
  sobj <- FindClusters(sobj, resolution = resolution, verbose = FALSE)

  cluster_number <- length(levels(Idents(sobj)))
    
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
    labs <- levels(Idents(sobj))
    cluster_num <- as.numeric(sub(".*\\\\(([0-9]+)\\\\).*", "\\\\1", labs))
    new_levels <- labs[order(cluster_num)]
    Idents(sobj) <- factor(Idents(sobj), levels = new_levels)
    sobj$shiny_clusters <- Idents(sobj)
    Idents(sobj) <- sobj$shiny_clusters
    
    umap_plot <- DimPlot(sobj, reduction = "umap", label = FALSE) + 
      ggtitle(NULL) +
      theme(plot.title = element_text(face = "bold", size = 10)) +
      guides(color = guide_legend(ncol = 1,
                                  override.aes = list(size = 3)))
    
    table_data <- get_cluster_summary_table(sobj)

    combined_plot <- create_combined_cluster_plot_patchwork(
      umap_plot = umap_plot,
      table_data = table_data,
      heatmap_plot = heatmap_plot,
      dotplot1 = dotplot1_plot,
      dotplot2 = dotplot2_plot,
      title_text = paste0("[", index, "] ","Resolution: ", resolution, ", PCA: ", pca, ", Cluster Number: ", cluster_number)
    )
    
    png_file <- file.path(
    temp_png_folder,
    sprintf("Index:%s_Resolution:%.2f_PCA:%d_ClusterNumber:%d.png", index, resolution, pca, cluster_number)
  )

  png(png_file, width = 18 * 200, height = 18 * 200, res = 200, type = "cairo")
  print(combined_plot)
  dev.off()

  rm(sobj, combined_plot)
  gc()

  list(
    index = index,
    cluster_number = cluster_number,
    png_file = png_file,
    label_name = findcelltypes_result$new_labels,
    table_data = table_data,
    heatmap_data = findcelltypes_result$heatmap_table
  )
}, future.seed = 001203)

  
  message("All resolutions processed! Creating final PDF...")
    resolution_summary <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      index = x$index,
      cluster_number = x$cluster_number
    )
  }))%>%
    left_join(resolution_index_table, by = "index") %>%
    relocate(cluster_number, .before = resolution)
  
  final_df <- do.call(rbind, lapply(results, process_single_result)) %>%
    left_join(resolution_index_table, by = "index") %>%
    relocate(index, resolution, pca) %>%   
    mutate(index = as.numeric(index)) %>%
    select(-cluster)
    
  all_png_files <- setNames(
    vapply(results, function(x) x$png_file, ""),
    vapply(results, function(x) x$index, "")
  )
  final_pdf <- file.path(output_folder, "All_Resolutions_Combined.pdf")
  
  pdf(final_pdf, width = 18, height = 18)
  n_rows <- nrow(resolution_summary)
  gc()
  # First Page
  grid::grid.newpage()
  
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow = 2, 
        ncol = 1,
        heights = grid::unit(c(0.1, 0.9), "npc")
      )
    )
  )
  
  # =========================
  # Title block (row 1)
  # =========================
  grid::pushViewport(grid::viewport(layout.pos.row = 1))
  
  grid::grid.text(
    "Cluster Analysis Summary", 
    x = 0.5, y = 0.75,
    gp = grid::gpar(fontsize = 36, fontface = "bold")
  )
  
  grid::grid.text(
    paste("Generated:", Sys.time()),
    x = 0.5, y = 0.5,
    gp = grid::gpar(fontsize = 20)
  )
  
  grid::grid.text(
    paste("Folder:", output_file_path),
    x = 0.5, y = 0.3,
    gp = grid::gpar(fontsize = 20)
  )
  
  grid::popViewport()  # title block
  
  # =========================
  # Tables block (row 2)
  # =========================
  grid::pushViewport(grid::viewport(layout.pos.row = 2))
  
  part_size <- ceiling(n_rows / 3)
  
  tables <- list()
  for (i in 1:3) {
    start_row <- (i - 1) * part_size + 1
    end_row <- min(i * part_size, n_rows)
  
    if (start_row <= end_row) {
      sub_data <- resolution_summary[start_row:end_row, c("index", "resolution", "pca", "cluster_number")]
      colnames(sub_data) <- c("Index" ,"Resolution", "PCA", "Clusters")
  
      tables[[i]] <- gridExtra::tableGrob(
        sub_data,
        rows = NULL,
        theme = gridExtra::ttheme_minimal(
          base_size = 16,
          padding = grid::unit(c(3, 3), "mm"),
          core = list(
            fg_params = list(hjust = 0.5, x = 0.5)
          ),
          colhead = list(
            fg_params = list(fontface = "bold", hjust = 0.5, x = 0.5)
          )
        )
      )
    } else {
      tables[[i]] <- grid::nullGrob()
    }
  }
  
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow = 1,
        ncol = 3,
        widths = grid::unit(c(0.33, 0.34, 0.33), "npc"),
        respect = TRUE
      )
    )
  )
  
  for (i in 1:3) {
    grid::pushViewport(grid::viewport(layout.pos.col = i))
    if (!inherits(tables[[i]], "null")) {
      grid::grid.draw(tables[[i]])
    }
    grid::popViewport()
  }
  
  grid::popViewport()  # table layout
  grid::popViewport()  # row 2
  grid::popViewport()  # page
  plots_per_page <- 4
  total_pages <- ceiling(length(all_png_files) / plots_per_page)
  
  for (page in 1:total_pages) {
    
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, length(all_png_files))
    
    grid::grid.newpage()
    
    pushViewport(viewport(
      width  = unit(1, "npc"),
      height = unit(1, "npc"),
      layout = grid.layout(
        nrow = 2,
        ncol = 2,
        heights = unit(rep(1, 2), "null"),
        widths  = unit(rep(1, 2), "null"),
        respect = TRUE
      )
    ))
    
    
    cell_positions <- list(
      c(1, 1),  
      c(1, 2), 
      c(2, 1), 
      c(2, 2)  
    )
    
    for (i in seq_along(start_idx:end_idx)) {
      idx <- start_idx + i - 1
      if (idx <= length(all_png_files)) {
        png_file <- all_png_files[idx]
        
        if (file.exists(png_file)) {
          pushViewport(viewport(layout.pos.row = cell_positions[[i]][1],
                                layout.pos.col = cell_positions[[i]][2]))
          
          img <- png::readPNG(png_file)
          grid::grid.raster(
            img,
            width  = unit(1, "npc"),
            height = unit(1, "npc"),
            interpolate = FALSE
          )
          
          popViewport()
        }
      }
    }
    
    popViewport()
    
    if (end_idx - start_idx + 1 < plots_per_page) {
      empty_cells <- plots_per_page - (end_idx - start_idx + 1)
    }
  }
  
  dev.off()
  
  summary_file <- file.path(output_folder, "Resolution_ClusterNumber_Summary.xlsx")
  
  # SAVE EXCEL
  wb <- createWorkbook()
  addWorksheet(wb, "Sheet1")
  
  info_row        <- 1
  max_row         <- 2
  min_row         <- 3
  header_row      <- 4
  data_start_row  <- 5
  
  start_col <- 10
  end_col   <- ncol(final_df)
  info <- paste(
    paste0("Generated: ", Sys.time()),
    paste0("Folder: ", output_file_path),
    paste0("Scale used variable features: ", ana_info$scale_used_variable_features),
    paste0("PCA npcs: ", ana_info$pca_npcs),
    paste0("harmony integration: ", integration_harmony),
    if(integration_harmony){ paste0("harmony metadata: ", harmony_metadata) },
    sep = ", "
  )
  writeData(wb, "Sheet1", x = info, startRow = info_row, startCol = 1)
  writeData(wb, "Sheet1", x = "MIN:", startRow = min_row, startCol = 9)
  writeData(wb, "Sheet1", x = "MAX:", startRow = max_row, startCol = 9)
  
  num_mat <- final_df[, start_col:end_col, drop = FALSE]
  
  num_mat[] <- lapply(num_mat, function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[!is.finite(x)] <- NA  
    x
  })
  
  max_vals <- vapply(num_mat, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  }, numeric(1))
  
  min_vals <- vapply(num_mat, function(x) {
    if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  }, numeric(1))
  
  writeData(
    wb, "Sheet1",
    x = t(as.matrix(max_vals)),
    startRow = max_row,
    startCol = start_col,
    colNames = FALSE
  )
  
  writeData(
    wb, "Sheet1",
    x = t(as.matrix(min_vals)),
    startRow = min_row,
    startCol = start_col,
    colNames = FALSE
  )
  
  writeData(
    wb, "Sheet1",
    x = final_df,
    startRow = header_row,
    startCol = 1,
    colNames = TRUE
  )
  
  data_end_row <- 4 + nrow(final_df)
  
  conditionalFormatting(
    wb, "Sheet1",
    cols = start_col:end_col,
    rows = 2:data_end_row,
    style = c("#003CB8", "#FFFFFF", "#FF6347"),
    type = "colorScale"
  )
  
  header_style <- createStyle(
    border = "TopBottomLeftRight",
    borderStyle = "thin",
    textDecoration = "bold"
  )
  
  addStyle(
    wb, "Sheet1",
    style = header_style,
    rows = header_row,
    cols = 1:ncol(final_df),
    gridExpand = TRUE
  )
  
  saveWorkbook(wb, summary_file, overwrite = TRUE)
  
  message("Summary saved to: ", summary_file)
  
  unlink(temp_png_folder, recursive = TRUE)
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
    "#$ -l mem_per_core=16G\n",
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