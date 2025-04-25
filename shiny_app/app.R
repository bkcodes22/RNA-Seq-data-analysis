library(shiny)
suppressWarnings(suppressPackageStartupMessages(library(shinyjs)))
suppressWarnings(suppressPackageStartupMessages(library(DT)))
suppressWarnings(suppressPackageStartupMessages(library(limma)))
library(edgeR)
library(ggplot2)
library(FactoMineR)
library(RColorBrewer)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Differential Gene Expression Analysis"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("input_mode", "Select input source:",
                   choices = c("From raw FASTQ files" = "fastq",
                               "From existing count file" = "counts")),
      textInput("dirpath", "Directory containing *.tab files:",
                placeholder = "/path/to/bam_files/"),
      selectInput("col_choice", "Choose count column to keep:",
                  choices = c("unstranded", "same_stranded", "reverse_stranded")),
      actionButton("go_fastq", "Load, Merge & Save"),
      hr(),
      textInput("csvpath", "Path to existing counts CSV:",
                placeholder = "/path/to/raw_counts.csv"),
      actionButton("go_counts", "Load Counts File"),
      actionButton("normalize_tmm", "Normalize Counts (TMM)", icon = icon("magic")),
      actionButton("do_diffexp", "Run DE Analysis", icon = icon("chart-line")),
      actionButton("do_pca", "Run PCA Analysis", icon = icon("project-diagram")),
      selectInput("pca_format", "Download PCA as:", choices = c("PDF", "PNG", "SVG")),
      downloadButton("downloadPCA", "Download PCA Plot"),
      verbatimTextOutput("status")
    ),
    mainPanel(
      DTOutput("mergedTable"),
      plotOutput("pcaPlot")
    )
  )
)

server <- function(input, output, session) {
  observe({
    toggleState("dirpath",   input$input_mode == "fastq")
    toggleState("col_choice", input$input_mode == "fastq")
    toggleState("go_fastq",  input$input_mode == "fastq")
    toggleState("csvpath",   input$input_mode == "counts")
    toggleState("go_counts", input$input_mode == "counts")
  })
  
  fastq_result <- eventReactive(input$go_fastq, {
    req(input$input_mode == "fastq")
    validate(need(dir.exists(input$dirpath), "❌ Directory not found."))
    fns <- list.files(input$dirpath, pattern = "ReadsPerGene\\.out\\.tab$", full.names = TRUE)
    validate(need(length(fns) > 0, "❌ No .tab files found."))
    proc <- function(p) {
      d <- read.csv(p, sep = "\t", header = FALSE)[-1:-4, ]
      colnames(d) <- c("gene", "unstranded", "same_stranded", "reverse_stranded")
      s <- sub("ReadsPerGene\\.out\\.tab$", "", basename(p))
      o <- d[, c("gene", input$col_choice)]; colnames(o)[2] <- s; o
    }
    L <- lapply(fns, proc)
    M <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), L)
    M <- M[order(M$gene), ]
    od <- file.path(input$dirpath, "diff_exp")
    if (!dir.exists(od)) dir.create(od, recursive = TRUE)
    wf <- file.path(od, "raw_counts.csv")
    write.csv(M, wf, row.names = FALSE, quote = TRUE)
    list(data = M, msg = paste0("✅ raw_counts.csv saved to:\n  ", wf))
  })
  
  loaded <- reactiveVal(NULL)
  kept_cols <- reactiveVal(NULL)
  metadata <- reactiveVal(NULL)
  design_matrix <- reactiveVal(NULL)
  dge_obj <- reactiveVal(NULL)
  
  observeEvent(input$go_counts, {
    req(input$input_mode == "counts")
    validate(need(file.exists(input$csvpath), "❌ CSV not found."))
    df <- read.csv(input$csvpath, header = TRUE, check.names = FALSE)
    if (colnames(df)[1] != "gene") colnames(df)[1] <- "gene"
    loaded(df)
    samples <- setdiff(colnames(df), "gene")
    showModal(modalDialog(
      title = "Select Columns to Retain",
      checkboxGroupInput("keep_cols", "Choose sample columns:",
                         choices = samples, selected = samples),
      footer = tagList(modalButton("Cancel"), actionButton("submit_keep", "OK")),
      easyClose = FALSE
    ))
  })
  
  observeEvent(input$submit_keep, {
    req(loaded(), input$keep_cols)
    removeModal()
    df <- loaded()[, c("gene", input$keep_cols), drop = FALSE]
    kept_cols(input$keep_cols)
    loaded(df)
    showModal(modalDialog(
      title = "Rename Kept Columns?",
      "Would you like to rename the selected sample columns?",
      footer = tagList(modalButton("No, thanks"), actionButton("start_rename", "Yes, rename"))
    ))
  })
  
  observeEvent(input$start_rename, {
    req(loaded(), kept_cols())
    removeModal()
    cols <- kept_cols()
    showModal(modalDialog(
      title = "Enter New Names",
      lapply(cols, function(cn) {
        textInput(paste0("nm_", cn), label = cn, value = cn)
      }),
      footer = tagList(modalButton("Cancel"), actionButton("apply_rename", "Apply")),
      easyClose = FALSE
    ))
  })
  
  observeEvent(input$apply_rename, {
    req(loaded(), kept_cols())
    df <- loaded()
    for (cn in kept_cols()) {
      newn <- input[[paste0("nm_", cn)]]
      if (nzchar(newn) && !newn %in% colnames(df)) {
        colnames(df)[colnames(df) == cn] <- newn
      }
    }
    loaded(df)
    removeModal()
    showModal(modalDialog(
      title = "Save Processed File?",
      "Save results to diff_exp/raw_counts_for_diff_exp_final.csv?",
      footer = tagList(modalButton("No"), actionButton("save_final", "Yes"))
    ))
  })
  
  observeEvent(input$save_final, {
    req(loaded())
    removeModal()
    base_dir <- if (input$input_mode == "fastq") input$dirpath else dirname(input$csvpath)
    out_dir <- file.path(base_dir, "diff_exp")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_file <- file.path(out_dir, "raw_counts_for_diff_exp_final.csv")
    write.csv(loaded(), out_file, row.names = FALSE)
    output$status <- renderText(paste0("✅ Final counts saved to:\n  ", out_file))
  })
  
  observeEvent(input$normalize_tmm, {
    req(loaded())
    df <- loaded()
    samples <- colnames(df)[-1]
    showModal(modalDialog(
      title = "Assign Metadata Groups to Samples",
      lapply(samples, function(smp) {
        textInput(paste0("grp_", smp), label = paste("Group for", smp), placeholder = "e.g., Control / Treated")
      }),
      footer = tagList(modalButton("Cancel"), actionButton("submit_groups", "Normalize Now")),
      easyClose = FALSE
    ))
  })
  
  observeEvent(input$submit_groups, {
    req(loaded())
    removeModal()
    df <- loaded()
    samples <- colnames(df)[-1]
    genes <- df$gene
    counts <- df[, -1, drop = FALSE]
    grp_vector <- sapply(samples, function(smp) input[[paste0("grp_", smp)]])
    metadata(data.frame(sample = samples, group = grp_vector, stringsAsFactors = FALSE))
    groups <- factor(grp_vector)
    design <- model.matrix(~0 + groups)
    colnames(design) <- levels(groups)
    design_matrix(design)
    
    count_matrix <- as.matrix(counts)
    rownames(count_matrix) <- genes
    dge <- DGEList(counts = count_matrix)
    dge <- calcNormFactors(dge, method = "TMM")
    dge_obj(dge)
    tmm_counts <- cpm(dge, log = FALSE)
    norm_df <- data.frame(gene = rownames(tmm_counts), round(tmm_counts, 3), check.names = FALSE)
    
    base_dir <- if (input$input_mode == "fastq") input$dirpath else dirname(input$csvpath)
    out_dir <- file.path(base_dir, "diff_exp")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    tmm_file <- file.path(out_dir, "tmm_normalized_counts.csv")
    write.csv(norm_df, tmm_file, row.names = FALSE)
    
    loaded(norm_df)
    output$status <- renderText(paste0("✅ TMM-normalized counts saved to:\n  ", tmm_file))
  })
  
  observeEvent(input$do_diffexp, {
    req(dge_obj(), design_matrix(), metadata())
    showModal(modalDialog(
      title = "Specify Contrast",
      textInput("contrast_text", "Enter contrast (e.g., Control - Treated):", ""),
      footer = tagList(modalButton("Cancel"), actionButton("run_contrast", "Run")),
      easyClose = FALSE
    ))
  })
  
  observeEvent(input$run_contrast, {
    req(input$contrast_text, dge_obj(), design_matrix())
    removeModal()
    logCPM <- cpm(dge_obj(), log = TRUE, prior.count = 3)
    fit <- lmFit(logCPM, design_matrix())
    cont.mat <- makeContrasts(contrasts = input$contrast_text, levels = design_matrix())
    fit2 <- contrasts.fit(fit, cont.mat)
    fit2 <- eBayes(fit2, trend = TRUE)
    res <- topTable(fit2, n = nrow(logCPM))
    
    # Show modal to ask for logFC cutoff
    showModal(modalDialog(
      title = "Specify logFC Cutoff",
      numericInput("logfc_cutoff", "Enter logFC cutoff:", value = 1, min = 0.1, step = 0.1),
      footer = tagList(modalButton("Cancel"), actionButton("apply_cutoff", "Apply")),
      easyClose = FALSE
    ))
    
    # Save the results and apply padj < 0.05 and user-defined logFC cutoff
    observeEvent(input$apply_cutoff, {
      removeModal()
      
      logfc_cutoff <- input$logfc_cutoff
      
      # Filter results based on padj < 0.05 and logFC cutoff
      filtered_res <- res[res$adj.P.Val < 0.05 & abs(res$logFC) >= logfc_cutoff, ]
      
      # Count upregulated and downregulated genes
      upregulated <- sum(filtered_res$logFC > 0)
      downregulated <- sum(filtered_res$logFC < 0)
      
      output$status <- renderText({
        paste0("✅ Differential expression results saved to:\n  ", result_file,
               "\n\nUpregulated genes (logFC > ", logfc_cutoff, "): ", upregulated,
               "\nDownregulated genes (logFC < -", logfc_cutoff, "): ", downregulated)
      })
      
      # Save filtered results
      base_dir <- if (input$input_mode == "fastq") input$dirpath else dirname(input$csvpath)
      out_dir <- file.path(base_dir, "diff_exp")
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      contrast_name <- gsub("[[:space:]]+", "_", input$contrast_text)
      result_file <- file.path(out_dir, paste0("diff_exp_results_", contrast_name, ".csv"))
      write.csv(filtered_res, result_file, row.names = TRUE)
    })
  })
  
  
  observeEvent(input$do_pca, {
  req(dge_obj(), metadata())
  logCPM <- cpm(dge_obj(), log = TRUE, prior.count = 3)
  pca_res <- PCA(t(logCPM), scale.unit = TRUE, ncp = 5, graph = FALSE)
  var_explained <- round((pca_res$eig[1:2, "percentage of variance"]), 2)
  pca_df <- data.frame(pca_res$ind$coord)
  pca_df$Group <- metadata()$group
  
  # Simplified color palette
  colors <- c("red", "yellow", "green", "black", "blue", "orange", "purple", "pink", "brown", "gray")
  
  # Ensure the number of colors matches the number of groups
  if (length(unique(pca_df$Group)) > length(colors)) {
    colors <- rep(colors, length.out = length(unique(pca_df$Group)))
  }
  
  p <- ggplot(pca_df, aes(x = Dim.1, y = Dim.2, color = Group)) +
    geom_point(size = 4, alpha = 0.85) +  # Simplified with just dots
    scale_color_manual(values = colors) +  # Apply the color scheme
    labs(
      title = "PCA on logCPM Counts",
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_blank()
    )
  
  output$pcaPlot <- renderPlot(p)
  assign("current_pca_plot", p, envir = .GlobalEnv)
})

  
  output$downloadPCA <- downloadHandler(
    filename = function() {
      paste0("PCA_plot.", tolower(input$pca_format))
    },
    content = function(file) {
      ggsave(file, plot = get("current_pca_plot", envir = .GlobalEnv), device = tolower(input$pca_format), width = 8, height = 6)
    }
  )
  
  output$mergedTable <- renderDT({
    if (input$input_mode == "fastq") {
      req(fastq_result()); fastq_result()$data
    } else {
      req(loaded()); loaded()
    }
  }, options = list(pageLength = 10))
  
  output$status <- renderText({
    if (input$input_mode == "fastq" && !is.null(fastq_result())) fastq_result()$msg
    else ""
  })
}

shinyApp(ui, server)
