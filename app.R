#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(renv)
library(shiny)
library(shinythemes)
library(readxl)
# library(IOBR) ## an error has occur
library(immunedeconv)
# library(EPIC)
# library(xCell)
library(ggplot2)
library(reshape2)
library(DESeq2)
### fins aqui finsionava.
library(pheatmap)
##
library(PoiClaClu)
library(grid)
library(stats)
library(grDevices)

library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(ggpubr)
library(maxstat)
library(org.Hs.eg.db)
library(stats)
library(ggupset)



ui <- fluidPage(
  theme = shinytheme("united"),
  tags$div(
    tags$img(src = "HTGAnalyzer_logo.png", height = "150px", style = "display: block; margin-left: auto; margin-right: auto;"),
    style = "text-align: center;"
  ),

  titlePanel("HTGAnalyzer Shiny App"),

  sidebarLayout(
    sidebarPanel(
      helpText("Select your data files to begin the analysis. Upload counts and annotation files."),
      fileInput("counts_file", "Upload counts file (.xlsx)", accept = c(".xlsx")),
      helpText("Upload your annotation file with sample data. The 'id' column should be specified."),
      fileInput("annot_file", "Upload annotation file (.xlsx)", accept = c(".xlsx")),
      helpText("Select the type of data you're working with. Choose 'HTG' for HTG EdgeSeq data or 'RNAseq' for RNA-seq data."),
      selectInput("file_type", "Select file type", choices = c("HTG", "RNAseq"), selected = "HTG"),
      helpText("Choose if you want to perform quality control on the data. This step is available for HTG data."),
      checkboxInput("QC", "Quality Control (QC)", FALSE),

      conditionalPanel(
        condition = "input.file_type == 'HTG' && input.QC == true",
        helpText("Set thresholds for various quality control parameters. Adjust as needed based on the data."),
        numericInput("threshold_superior_pos", "Threshold superior pos", 5),
        numericInput("threshold_line_pos", "Threshold line pos", 4),
        numericInput("threshold_inferior_lib", "Threshold inferior lib", 5e+06),
        numericInput("threshold_lib", "Threshold lib", 7e+06),
        numericInput("threshold_superior_nc", "Threshold superior nc", 0.05),
        numericInput("threshold_line_nc", "Threshold line nc", 0.045),
        numericInput("threshold_superior_gdna", "Threshold superior gdna", 0.025),
        numericInput("threshold_line_gdna", "Threshold line gdna", 0.02),
        numericInput("threshold_superior_ercc", "Threshold superior ercc", 0.03),
        numericInput("threshold_line_ercc", "Threshold line ercc", 0.025),
        numericInput("threshold_inferior_median", "Threshold inferior median", 3),
        numericInput("threshold_line_median", "Threshold line median", 5),
        helpText("Use a regular expression pattern to identify control probes in the count data."),
        textInput("pattern", "Pattern", "^NC-|^POS-|^GDNA-|^ERCC-"),
        helpText("Choose whether to save QC results as a CSV file. Specify the file name."),
        checkboxInput("save_csv", "Save results in CSV", TRUE),
        conditionalPanel(
          condition = "input.save_csv == true",
          textInput("csv_file", "CSV file name", "QC_results.csv")
        ),
        helpText("Choose whether to remove outliers during analysis."),
        checkboxInput("remove_outliers", "Remove Outliers", TRUE)
      ),

      helpText("Select if you want to perform differential expression analysis (DEA), GSEA and heatmap."),
      checkboxInput("DEA", "Differential Expression Analysis (DEA), GSEA and Heatmap", FALSE),

      conditionalPanel(
        condition = "input.DEA == true",
        helpText("Specify the design formula and contrast for DEA. Example: Design formula could be 'HPV_status'."),
        textInput("design_formula", "Design Formula", "HPV_status"),
        textInput("contrast_input", "Contrast", 'c("HPV_status", "Positive", "Negative")'),
        helpText("Set thresholds for gene expression filtering."),
        numericInput("percentage_gene", "Percentage gene", 0.2),
        numericInput("threshold_gene", "Threshold gene", 200),
        numericInput("threshold_subject", "Threshold subject", 10),
        helpText("Choose if you want to perform Gene Set Enrichment Analysis (GSEA) as part of DEA."),
        checkboxInput("GSEA", "Gene Set Enrichment Analysis (GSEA)", FALSE)
      ),

      conditionalPanel(
        condition = "input.DEA == true",
        helpText("Choose if you want to generate a heatmap for differential expression analysis."),
        checkboxInput("generate_heatmap", "Generate Heatmap", FALSE),
        conditionalPanel(
          condition = "input.generate_heatmap == true",
          textInput("heatmap_columns", "Heatmap columns", "HPV_status,Ciclina_D1")
        )
      ),

      helpText("Select if you want to perform tumor microenvironment profiling (TME)."),
      checkboxInput("TME", "Tumor Microenvironment Profile (TME)", FALSE),

      helpText("Select if you want to perform survival analysis on the data."),
      checkboxInput("survival_analysis", "Survival Analysis", FALSE),

      conditionalPanel(
        condition = "input.survival_analysis == true",
        helpText("Specify the variable and time for survival analysis. Example: 'Recurrence_01' and 'Time_to_death_surv'."),
        textInput("variable_01", "Variable 01", "Recurrence_01"),
        textInput("time", "Time", "Time_to_death_surv"),
        textInput("genes_to_use", "Genes to use", "CCND1,MMP10,CTTN")
      ),

      actionButton("run_analysis", "Run Analysis")
    ),

    mainPanel(
      h2("DESCRIPTION HTG ANALYZER SHINY APP"),
      p("This app enables easy analysis of HTG Edge and RNAseq data through the ", strong("HTG_auto"), " function."),
      tags$ul(
        tags$li(strong("Upload Counts & Annotation Files:"), " Add input files to begin the analysis."),
        tags$li(strong("File Type:"), " Choose HTG or RNAseq as your data source."),
        tags$li(strong("Quality Control (QC):"), " Available for HTG data to filter and clean data, with customizable thresholds specifically designed for the HTG Transcriptomic Panel. QC not available for RNAseq data"),
        tags$li(strong("Differential Expression Analysis (DEA):"), " Identify differentially expressed genes with optional GSEA and heatmap generation."),
        tags$li(strong("Tumor Microenvironment (TME):"), " Profile the tumor microenvironment in your samples."),
        tags$li(strong("Survival Analysis:"), " Analyze survival data with selected clinical variables.")
      ),
      p("This app automates quality control, DEA, and more to facilitate HTG and RNAseq data analysis. It uses the same variable names as the ", strong("HTG_auto"), " function from the HTGAnalyzer package. For more information, visit the ", a("HTGAnalyzer GitHub Site", href = "https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer"), "."),
      h3("Estimated Analysis Times"),
      tags$ul(
        tags$li("Quality Control: 5 min (not required for RNAseq)"),
        tags$li("Differential Expression Analysis: 5–10 min"),
        tags$li("GSEA: 20–30 min (optional in DEA)"),
        tags$li("Tumor Microenvironment Analysis: 5–10 min"),
        tags$li("Survival Analysis: 5–10 min")
      ),
      textOutput("dimension_output"),
      textOutput("dimension_output1"),
      # RESULTS QC
      tabsetPanel(
        tabPanel("Quality Control",
      tableOutput("summary_stats"),
      tableOutput("ratiosb"),
      downloadButton("download_summary", "Download Summary Stats"),
      downloadButton("download_ratiosb", "Download QC Results"),
      plotOutput("pos_genes_plot", width = "100%", height = "400px"),
      downloadButton("download_pos_genes_plot", "Download Positive Control Plot (PDF)"),
      plotOutput("size_lib", width = "100%", height = "400px"),
      downloadButton("download_size_lib_plot", "Download Library Size Plot (PDF)"),
      plotOutput("neg", width = "100%", height = "400px"),
      downloadButton("download_neg_plot", "Download Negative Control Plot (PDF)"),
      plotOutput("gdna", width = "100%", height = "400px"),
      downloadButton("download_gdna_plot", "Download Genomic DNA Plot (PDF)"),
      plotOutput("ERCC", width = "100%", height = "400px"),
      downloadButton("download_ercc_plot", "Download ERCC Plot (PDF)"),
      plotOutput("med", width = "100%", height = "400px"),
      downloadButton("download_med_plot", "Download Median Plot (PDF)"),
      plotOutput("qc_plot_heatmap", width = "100%", height = "400px"),
      downloadButton("download_qc_plot_heatmap", "Download QC Heatmap PDF"),
      textOutput("outliers")
      ),
      # RESULTS DEA
      tabPanel("Differential Expression Analisis",
      textOutput("error_design_formula"),
      textOutput("filtering_data"),
      textOutput("removing_outliers"),
      textOutput("missatge2"),
      textOutput("dds_missatge"),
      verbatimTextOutput("dds_output"),
      textOutput("results_VST"),
      verbatimTextOutput("vsd_output"),
      textOutput("missatge3"),
      plotOutput("DEA_heatmap", width = "100%", height = "400px"),
      downloadButton("DEA_heatmap_download", "Download HEATMAP (PDF)"),
      textOutput("missatge_contrast"),
      textOutput("summary_results"),
      verbatimTextOutput("printsummary_results"),
      textOutput("top10_contrast"),
      verbatimTextOutput("printtop10_contrast"),
      downloadButton("download_res", "Download DEA results"),
      textOutput("missatge_size_factors_plot"),
      plotOutput("size_factors_plot", width = "100%", height = "400px"),
      downloadButton("download_DEA_plot", "Download Size Factors Plot"),
      plotOutput("correlation_heatmap", width = "100%", height = "400px"),
      textOutput("correlation_heatmap_mes"),
      downloadButton("download_correlation_heatmap", "Download correlation heatmap"),
      plotOutput("poisson_heatmap", width = "100%", height = "400px"),
      textOutput("missatge4"),
      plotOutput("poison", width = "100%", height = "400px"),
      downloadButton("download_poison", "Download poisson distance"),
      textOutput("missatge5"),
      plotOutput("distance_heatmap", width = "100%", height = "400px"),
      downloadButton("download_distance_heatmap", "DEA Sample Distance Heatmap"),
      textOutput("missatge6"),
      plotOutput("vsd_boxplot", width = "100%", height = "400px"),
      downloadButton("download_vsd_boxplot", "Download Boxplot VSD"),
      textOutput("topgenes"),
      textOutput("missatge7"),
      plotOutput("pca_plot", width = "100%", height = "400px"),
      downloadButton("download_pca_plot", "Download PCA plot"),
      textOutput("missatge8"),
      textOutput("missatge9"),
      plotOutput("maPlot", width = "100%", height = "400px"),
      downloadButton("downloadMAPlot", "Download MA plot"),
      plotOutput("volvano", width = "100%", height = "400px"),
      downloadButton("downloadvolvano", "Download Volcano plot"),
      textOutput("missatge10"),
      plotOutput("dispersion", width = "100%", height = "400px"),
      downloadButton("downloaddispersion", "Download Dispersion plot")
      ),
      # RESULTS GSEA
      tabPanel("GSEA",
      plotOutput("dotplot1", width = "100%", height = "400px"),
      downloadButton("download_dotplot1", "Download Dotplot GSEA 1"),
      plotOutput("dotplot2", width = "100%", height = "400px"),
      downloadButton("download_dotplot2", "Download Dotplot GSEA 2"),
      plotOutput("emapplot1"),
      downloadButton("download_emapplot", "Download Emapplot GSEA"),
      plotOutput("ridgeplot1"),
      downloadButton("download_ridgeplot", "Download Ridgeplot GSEA"),
      plotOutput("heatplot1"),
      downloadButton("download_heatplot", "Download Heatplot GSEA"),
      plotOutput("treeplot1"),
      downloadButton("download_treeplot1", "Download Treeplot GSEA"),
      plotOutput("a"),
      downloadButton("download_a", "Download GSEA Plot"),
      plotOutput("b"),
      downloadButton("download_b", "Download Top 5 Gene Sets"),
      plotOutput("gsea_dotplot3"),
      downloadButton("download_dotplot3", "Download KEGG Dotplot"),
      plotOutput("gsea_emapplot2"),
      downloadButton("download_emapplot2", "Download KEGG Emapplot"),
      plotOutput("gsea_ridgeplot2"),
      downloadButton("download_ridgeplot2", "Download KEGG Ridgeplot"),
      plotOutput("gsea_heatplot2"),
      downloadButton("download_heatplot2", "Download KEGG Heatplot"),
      plotOutput("gsea_treeplot2"),
      downloadButton("download_treeplot2", "Download KEGG Treeplot"),
      plotOutput("go_bar_plot"),
      downloadButton("download_go_bar_plot", "Download GO Bar Plot"),
      plotOutput("go_bar_plot2"),
      downloadButton("download_go_bar_plot2", "Download GO Bar Plot")
      ),
      # RESULTS TPM
      tabPanel("Tumor Microenvironment",
               textOutput("genes_omitted"),
      textOutput("initial_message"),
      textOutput("normalization_message"),
      tableOutput("tpm_table"),
      downloadButton("download_tpm", "Download TPM CSV"),
      textOutput("cat_message"),
      tableOutput("imm_epic_table"),
      downloadButton("download_epic", "Download imm_epic.csv"),
      tableOutput("imm_qti_table"),
      downloadButton("download_qti", "Download imm_qti.csv"),
      tableOutput("imm_xcell_table"),
      downloadButton("download_xcell", "Download imm_xcell.csv"),
      verbatimTextOutput("TME"),
      plotOutput("plot_cell_fraction_EPIC"),
      downloadButton("download_plot_EPIC", "Descarregar EPIC"),
      plotOutput("plot_cell_fraction_quanTIseq"),
      downloadButton("download_plot_quanTIseq", "Descarregar quanTIseq"),
      plotOutput("plot_cell_fraction_xCell"),
      downloadButton("download_plot_xCell", "Descarregar xCell"),
      verbatimTextOutput("resultados_norm_imm_qtiA"),
      verbatimTextOutput("parametric_results_imm_qtiA"),
      verbatimTextOutput("resultados_norm_epicA"),
      verbatimTextOutput("parametric_results_epicA"),
      verbatimTextOutput("resultados_norm_xcellA"),
      verbatimTextOutput("parametric_results_xcellA"),
      plotOutput("heatmap_qti"),
      downloadButton("download_heatmap_qti", "Download Heatmap qti"),
      plotOutput("heatmap_EPIC"),
      downloadButton("download_heatmap_EPIC", "Download Heatmap EPIC"),
      plotOutput("heatmap_xcell"),
      downloadButton("download_heatmap_xcell", "Download Heatmap xcell")
      ),
      # SURVIVAL
      tabPanel("Survival Analysis",
               textOutput("normalizing"),
               textOutput("logrank"),
               textOutput("checkduplicates"),
               textOutput("duplicatesfound"),
               textOutput("noduplicates"),
               textOutput("top10"),
               verbatimTextOutput("genes_used"),
               verbatimTextOutput("genes_used2"),
               textOutput("startingsurvival"),
               textOutput("cutoffmessage"),
               textOutput("summaryfit"),
               tableOutput("fit1_df"),
               downloadButton("download_fit1_df", "Download Fit1 Table"),
               verbatimTextOutput("surv_diff_output"),
               textOutput("survival")

      )
    )
    )
  )
)





# Servidor para la aplicación
server <- function(input, output, session) {
  observeEvent(input$file_type, {
    if (input$file_type == "RNAseq") {
      updateCheckboxInput(session, "QC", value = FALSE)
    }
  })

  output$file_message <- renderText({
    req(input$run_analysis)
    paste("Note: The analysis will generate PDF reports, which will be saved in the working directory.")
  })

  observeEvent(input$run_analysis, {
    req(input$counts_file)
    req(input$annot_file)

    # Verificación de variables en design_formula
    if (input$DEA) {
      design_formula_vars <- input$design_formula
    }

    # Convertir el input de contrast en un vector
    contrast_vector <- if (input$DEA) {
      eval(parse(text = input$contrast_input))
    } else {
      NULL
    }

    # Convertir heatmap_columns a vector
    heatmap_columns_vector <- if (input$generate_heatmap) {
      strsplit(input$heatmap_columns, ",")[[1]]
    } else {
      NULL
    }

    # Convertir genes_to_use a vector
    genes_to_use_vector <- strsplit(input$genes_to_use, ",")[[1]]

    # Llamar a HTG_auto
    counts_file_path = input$counts_file$datapath
    AnnotData_file_path = input$annot_file$datapath
    file_type = input$file_type
    QC = input$QC
    threshold_superior_pos = input$threshold_superior_pos
    threshold_line_pos = input$threshold_line_pos
    threshold_inferior_lib = input$threshold_inferior_lib
    threshold_lib = input$threshold_lib
    threshold_superior_nc = input$threshold_superior_nc
    threshold_line_nc = input$threshold_line_nc
    threshold_superior_gdna = input$threshold_superior_gdna
    threshold_line_gdna = input$threshold_line_gdna
    threshold_superior_ercc = input$threshold_superior_ercc
    threshold_line_ercc = input$threshold_line_ercc
    threshold_inferior_median = input$threshold_inferior_median
    threshold_line_median = input$threshold_line_median
    pattern = input$pattern
    save_csv = input$save_csv
    csv_file = if (input$save_csv) input$csv_file else "QC_results.csv"
    DEA = input$DEA
    GSEA = input$DEA && input$GSEA
    design_formula = input$design_formula
    contrast = contrast_vector
    percentage_gene = if (input$DEA) input$percentage_gene else NULL
    threshold_gene = if (input$DEA) input$threshold_gene else NULL
    threshold_subject = if (input$DEA) input$threshold_subject else NULL
    generate_heatmap = input$generate_heatmap
    heatmap_columns = heatmap_columns_vector
    TME = input$TME
    survival_analysis = input$survival_analysis
    variable_01 = if (input$survival_analysis) input$variable_01 else NULL
    time = if (input$survival_analysis) input$time else NULL
    genes_to_use = genes_to_use_vector
    remove_outliers = input$remove_outliers

    # Importar los datos
    if (file_type == "HTG") {
      # Leer y procesar archivo HTG
      htg_db <- readxl::read_excel(counts_file_path)
      htg_db <- as.data.frame(htg_db)
      rownames_db <- htg_db[[1]]
      htg_db <- htg_db[, -1]
      rownames(htg_db) <- rownames_db

      cts <- htg_db
      if (grepl("^Sample ID", rownames(cts)[1])) {
        cts <- cts[-c(1:4), ]
        cts <- apply(cts, 2, as.numeric)
        rownames(cts) <- rownames(htg_db)[-c(1:4)]
        cts <- as.data.frame(cts)
      }
      counts_data <- cts
    } else if (file_type == "RNAseq") {
      # Leer y procesar archivo RNAseq
      rna_db <- readxl::read_excel(counts_file_path)
      rna_db <- as.data.frame(rna_db)
      counts_data <- rna_db
    } else {
      stop("file_type has to be 'HTG' or 'RNAseq'")
    }
    ###########################
    ###########################

    counts_data <- as.data.frame(counts_data)
    col_data <- readxl::read_excel(AnnotData_file_path)
    col_data <- as.data.frame(col_data)

    # Normalizar los nombres de las columnas en col_data para que sean válidos en R
    colnames(col_data) <- gsub(" ", "_", colnames(col_data))

    # Obtener los nombres de las columnas en counts_data
    counts_colnames <- colnames(counts_data)
    rownames(col_data) <- col_data$id
    col_data_rownames <- rownames(col_data)

    common_ids <- base::intersect(counts_colnames, col_data_rownames)

    counts_data <- counts_data[, counts_colnames %in% common_ids]
    col_data <- col_data[common_ids, , drop = FALSE]

    dims_message <- paste("Dimensions of filtered counts file:", paste(dim(counts_data), collapse = " x "))
    output$dimension_output <- renderText({dims_message})

    dims_message1 <- paste("Dimensions of filtered annotation files:", paste(dim(col_data), collapse = " x "))
    output$dimension_output1 <- renderText({dims_message1})

    if (ncol(counts_data) == nrow(col_data)) {
      cat("The filtered data sets are now aligned.")
    } else {
      cat("There is a mismatch in the number of columns in counts_data_filtered and rows in col_data_filtered.")
    }

    # Perform quality control if QC is TRUE
    if (QC) {
      # Filter counts_data data
      counts_filtered <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
      min_values <- apply(counts_filtered, 2, min)
      max_values <- apply(counts_filtered, 2, max)
      mean_values <- apply(counts_filtered, 2, mean)
      median_values <- apply(counts_filtered, 2, median)
      mode_values <- apply(counts_filtered, 2, function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      })
      sd_values <- apply(counts_filtered, 2, sd)
      var_values <- apply(counts_filtered, 2, var)
      range_values <- apply(counts_filtered, 2, function(x) max(x) - min(x))
      quartile_1 <- apply(counts_filtered, 2, function(x) quantile(x, 0.25))
      quartile_3 <- apply(counts_filtered, 2, function(x) quantile(x, 0.75))
      iqr_values <- quartile_3 - quartile_1
      skewness_values <- apply(counts_filtered, 2, function(x) {
        n <- length(x)
        mean_x <- mean(x)
        sd_x <- sd(x)
        sum((x - mean_x)^3) / ((n - 1) * (sd_x^3))
      })
      kurtosis_values <- apply(counts_filtered, 2, function(x) {
        n <- length(x)
        mean_x <- mean(x)
        sd_x <- sd(x)
        sum((x - mean_x)^4) / ((n - 1) * (sd_x^4)) - 3
      })
      missing_values <- apply(counts_filtered, 2, function(x) sum(is.na(x)))
      cv_values <- sd_values / mean_values

      summary_stats <- data.frame(
        Min = min_values,
        Max = max_values,
        Mean = mean_values,
        Median = median_values,
        Mode = mode_values,
        SD = sd_values,
        Variance = var_values,
        Range = range_values,
        Q1 = quartile_1,
        Q3 = quartile_3,
        IQR = iqr_values,
        Skewness = skewness_values,
        Kurtosis = kurtosis_values,
        Missing = missing_values,
        CV = cv_values
      )
      summary_stats$ID <- rownames(summary_stats)
      output$summary_stats <- renderTable({summary_stats})

      # Configurar el archivo de descarga
      output$download_summary <- downloadHandler(
        filename = function() { "summary_stats.csv" },
        content = function(file) {
          write.csv(summary_stats, file, row.names = FALSE)
        })
      # Subsets
      cts_ERCC <- as.data.frame(subset(counts_data, grepl("^ERCC-", rownames(counts_data))))
      cts_NC <- as.data.frame(subset(counts_data, grepl("^NC-", rownames(counts_data))))
      cts_POS <- as.data.frame(subset(counts_data, grepl("^POS-", rownames(counts_data))))
      cts_GDNA <- as.data.frame(subset(counts_data, grepl("^GDNA-", rownames(counts_data))))

      # Ratios
      total_genes <- colSums(counts_filtered)
      total_POS <- colSums(cts_POS)
      total_NC <- colSums(cts_NC)
      total_GDNA <- colSums(cts_GDNA)
      total_ERCC <- colSums(cts_ERCC)

      # Calculate ratios
      ratios <- data.frame(
        total_POS = total_POS,
        total_GDNA = total_GDNA,
        total_genes = total_genes,
        total_NC = total_NC,
        total_ERCC = total_ERCC
      )
      ratios$`pos/genes` <- (ratios$total_POS / ratios$total_genes) * 100
      ratios$`gdna/genes` <- (ratios$total_GDNA / ratios$total_genes) * 100
      ratios$`nc/genes` <- (ratios$total_NC / ratios$total_genes) * 100
      ratios$`ERCC/genes` <- (ratios$total_ERCC / ratios$total_genes) * 100

      # Add median column
      ratios$median <- summary_stats$Median
      ratiosb<-ratios
      ratiosb$min<- summary_stats$Min
      ratiosb$max<- summary_stats$Max
      ratiosb$mean<- summary_stats$Mean

      # Add sample names as a factor column
      ratios$samples <- factor(rownames(ratios))

      # Calcula los valores de las columnas Q0 a Q5 usando los umbrales definidos
      ratiosb$Q0 <- ifelse(ratiosb$`pos/genes` < threshold_line_pos, "PASS", "FAIL")
      ratiosb$Q1 <- ifelse(ratiosb$total_genes > threshold_lib, "PASS", "FAIL")
      ratiosb$Q2 <- ifelse(ratiosb$total_NC > threshold_line_nc, "PASS", "FAIL")
      ratiosb$Q3 <- ifelse(ratiosb$`gdna/genes` < threshold_line_gdna, "PASS", "FAIL")
      ratiosb$Q4 <- ifelse(ratiosb$`ERCC/genes` < threshold_line_ercc, "PASS", "FAIL")
      ratiosb$Q5 <- ifelse(ratiosb$median > threshold_line_median, "PASS", "FAIL")

      # Agrega los rownames como una columna 'ID'
      ratiosb$ID <- rownames(ratiosb)

      # Extrae el número de paciente y ordena por este número
      ratiosb$ID_num <- as.numeric(gsub("patient_", "", ratiosb$ID))  # Extrae el número de ID
      ratiosb <- ratiosb[order(ratiosb$ID_num), ]  # Ordena por la columna ID_num
      ratiosb$ID_num <- NULL  # Elimina la columna auxiliar

      # Elimina los rownames originales y reorganiza las columnas según el orden deseado
      rownames(ratiosb) <- NULL
      ratiosb <- ratiosb[, c("ID", "total_POS", "total_NC", "total_GDNA", "total_ERCC",
                             "pos/genes", "Q0", "total_genes", "Q1",
                             "nc/genes", "Q2", "gdna/genes", "Q3",
                             "ERCC/genes", "Q4", "median", "Q5",
                             "min", "max", "mean")]
      ratiosb$QC_status <- ifelse(
        apply(ratiosb[, c("Q0", "Q1", "Q2", "Q3", "Q4", "Q5")], 1, function(x) any(x == "FAIL")),
        "FAIL",
        "PASS"
      )

      output$ratiosb <- renderTable({ratiosb})

      output$download_ratiosb <- downloadHandler(
        filename = function() { "QC_results_QC.csv" },
        content = function(file) {
          write.csv(ratiosb, file, row.names = FALSE)
        }
      )

      # Guardar ratiosb como CSV opcionalmente
      if (input$save_csv) {
        write.csv(ratiosb, input$csv_file, row.names = FALSE)
        cat("QC DATA SAVED AS '", input$csv_file, "'\n")
      }

      # Optionally save as CSV
      if (save_csv) {
        write.csv(ratiosb, csv_file, row.names = FALSE)
        cat("QC DATA SAVED AS '", csv_file, "'")
      }
      ###
      library_size <- colSums(counts_filtered)

      # Create dataframe of library size
      lib_s2 <- data.frame(Sample = colnames(counts_filtered), Size = library_size)
      ratios_heat <- as.matrix(ratios)

      # Add a fourth column to ratios_heat with library sizes from lib_s2
      ratios_heat <- cbind(ratios_heat, Size = "")
      ratios_heat[, "Size"] <- lib_s2[match(rownames(ratios_heat), rownames(lib_s2)), "Size"]
      ratios_heat <- as.data.frame(ratios_heat)

      # Convert values of ratios_heat to numeric
      cols_to_convert <- c("total_POS", "total_GDNA", "total_genes", "total_NC", "total_ERCC",
                           "pos/genes", "gdna/genes", "nc/genes", "ERCC/genes", "median", "Size")
      ratios_heat[cols_to_convert] <- lapply(ratios_heat[cols_to_convert], as.numeric)
      str(ratios_heat)

      assign_01_QC <- function(valor, threshold) {
        ifelse(valor < threshold, 0, 1)
      }

      # Function to assign 0 or 1 according to library size value
      assign_01_size <- function(valor, threshold) {
        ifelse(valor > threshold, 0, 1)
      }

      # Create binary matrix for the heatmap
      bin_matrix <- matrix(0, nrow = nrow(ratios_heat), ncol = 6)
      for (i in 1:nrow(ratios_heat)) {
        bin_matrix[i, 1] <- assign_01_QC(ratios_heat[i, "pos/genes"], threshold_line_pos)
        bin_matrix[i, 2] <- assign_01_size(ratios_heat[i, "Size"], threshold_lib)
        bin_matrix[i, 3] <- assign_01_QC(ratios_heat[i, "nc/genes"], threshold_line_nc)
        bin_matrix[i, 4] <- assign_01_QC(ratios_heat[i, "gdna/genes"], threshold_line_gdna)
        bin_matrix[i, 5] <- assign_01_QC(ratios_heat[i, "ERCC/genes"], threshold_line_ercc)
        bin_matrix[i, 6] <- assign_01_size(ratios_heat[i, "median"], threshold_line_median)
      }

      # Row and column names
      rownames(bin_matrix) <- rownames(ratios_heat)
      colnames(bin_matrix) <- c("QC0", "QC1", "QC2", "QC3", "QC4", "QC5")

      # Convert the matrix to a data frame for ggplot2
      bin_df <- as.data.frame(bin_matrix)
      bin_df$Sample <- rownames(bin_df)
      # Melt the data frame
      bin_df_melted <- reshape2::melt(bin_df, id.vars = "Sample")

      max_value_pos <- max(ratios$`pos/genes`, threshold_line_pos)
      colores_pos <- ifelse(ratios$`pos/genes` <= threshold_line_pos, "#4793AF",
                            ifelse(ratios$`pos/genes` <= threshold_superior_pos, "#FFC470", "red"))
      min_size <- min(ratios$`pos/genes`, threshold_line_pos)
      colores <- ifelse(lib_s2$Size < threshold_inferior_lib, "red",
                        ifelse(lib_s2$Size <= threshold_lib, "#FFC470", "#4793AF"))
      max_size_lib <- max(lib_s2$Size, threshold_lib)
      min_size_lib <- min(lib_s2$Size, threshold_lib)
      max_value_neg <- max(ratios$`nc/genes`, threshold_line_nc)
      min_value_neg <- min(ratios$`nc/genes`, threshold_line_nc)
      colores_nc <- ifelse(ratios$`nc/genes` <= threshold_line_nc, "#4793AF",
                           ifelse(ratios$`nc/genes` <= threshold_superior_nc, "#FFC470", "red"))
      max_value_gdna <- max(ratios$`gdna/genes`, threshold_line_gdna)
      colores_gdna <- ifelse(ratios$`gdna/genes` <= threshold_line_gdna, "#4793AF",
                             ifelse(ratios$`gdna/genes` <= threshold_superior_gdna, "#FFC470", "red"))
      max_value_ERCC <- max(ratios$`ERCC/genes`, threshold_line_ercc)
      colores_ercc <- ifelse(ratios$`ERCC/genes` <= threshold_line_ercc, "#4793AF",
                             ifelse(ratios$`ERCC/genes` <= threshold_superior_ercc, "#FFC470", "red"))
      max_value_med <- max(ratios$median, threshold_line_median)
      colores_med <- ifelse(ratios$median < threshold_inferior_median, "red",
                            ifelse(ratios$median <= threshold_line_median, "#FFC470", "#4793AF"))

      ######
      # Plotting
      output$pos_genes_plot <- renderPlot({
        plot(ratios$`pos/genes`, xlab = "", ylab = "pos/genes", col = colores_pos,
             xaxt = "n", pch = 19, main = "Positive control 4% (QC0)", ylim = c(0, max_value_pos))
        axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
        abline(h = threshold_line_pos, col = "red")
      })
      output$size_lib <- renderPlot({
        plot(lib_s2$Size, xlab = "", ylab = "Library Size", col = colores,
             xaxt = "n", pch = 19, main = "Library Size per Sample (QC1)", cex.axis = 0.8,
             ylim = c(min_size_lib, max_size_lib))
        axis(1, at = 1:length(lib_s2$Sample), labels = lib_s2$Sample, las = 2, cex.axis = 0.8)
        abline(h = threshold_lib, col = "red")
      })
      output$neg <- renderPlot({
        plot(ratios$`nc/genes`, xlab = "", ylab = "nc/genes", col = colores_nc,
             xaxt = "n", pch = 19, main = "Negative Control (QC2)", ylim = c(min_value_neg, max_value_neg))
        axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
        abline(h = threshold_line_nc, col = "red")
      })
      output$gdna <- renderPlot({
        plot(ratios$`gdna/genes`, xlab = "", ylab = "gdna/genes", col = colores_gdna,
             xaxt = "n", pch = 19, main = "Genomic DNA (QC3)", ylim = c(0, max_value_gdna))
        axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
        abline(h = threshold_line_gdna, col = "red")
      })
      output$ERCC <- renderPlot({
        plot(ratios$`ERCC/genes`, xlab = "", ylab = "ERCC/genes", col = colores_ercc,
             xaxt = "n", pch = 19, main = "ERCC (QC4)", ylim = c(0, max_value_ERCC))
        axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
        abline(h = threshold_line_ercc, col = "red")
      })
      output$med <- renderPlot({
        plot(ratios$median, xlab = "", ylab = "Median", col = colores_med,
             xaxt = "n", pch = 19, main = "Median (QC5)", ylim = c(0, max_value_med))
        axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
        abline(h = threshold_line_median, col = "red")
      })

      #####
      output$download_pos_genes_plot <- downloadHandler(
        filename = function() { "positive_control_plot.pdf" },
        content = function(file) {
          pdf(file)
          plot(ratios$`pos/genes`, xlab = "", ylab = "pos/genes", col = colores_pos,
               xaxt = "n", pch = 19, main = "Positive control 4% (QC0)", ylim = c(0, max_value_pos))
          axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
          abline(h = threshold_line_pos, col = "red")
          dev.off()
        }
      )

      # Descargar el gráfico de Library Size
      output$download_size_lib_plot <- downloadHandler(
        filename = function() { "library_size_plot.pdf" },
        content = function(file) {
          pdf(file)
          plot(lib_s2$Size, xlab = "", ylab = "Library Size", col = colores,
               xaxt = "n", pch = 19, main = "Library Size per Sample (QC1)", cex.axis = 0.8,
               ylim = c(min_size_lib, max_size_lib))
          axis(1, at = 1:length(lib_s2$Sample), labels = lib_s2$Sample, las = 2, cex.axis = 0.8)
          abline(h = threshold_lib, col = "red")
          dev.off()
        }
      )

      # Descargar el gráfico de Negative Control
      output$download_neg_plot <- downloadHandler(
        filename = function() { "negative_control_plot.pdf" },
        content = function(file) {
          pdf(file)
          plot(ratios$`nc/genes`, xlab = "", ylab = "nc/genes", col = colores_nc,
               xaxt = "n", pch = 19, main = "Negative Control (QC2)", ylim = c(0, max_value_neg))
          axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
          abline(h = threshold_line_nc, col = "red")
          dev.off()
        }
      )

      # Descargar el gráfico de Genomic DNA
      output$download_gdna_plot <- downloadHandler(
        filename = function() { "gdna_plot.pdf" },
        content = function(file) {
          pdf(file)
          plot(ratios$`gdna/genes`, xlab = "", ylab = "gdna/genes", col = colores_gdna,
               xaxt = "n", pch = 19, main = "Genomic DNA (QC3)", ylim = c(0, max_value_gdna))
          axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
          abline(h = threshold_line_gdna, col = "red")
          dev.off()
        }
      )

      # Descargar el gráfico de ERCC
      output$download_ercc_plot <- downloadHandler(
        filename = function() { "ercc_plot.pdf" },
        content = function(file) {
          pdf(file)
          plot(ratios$`ERCC/genes`, xlab = "", ylab = "ERCC/genes", col = colores_ercc,
               xaxt = "n", pch = 19, main = "ERCC (QC4)", ylim = c(0, max_value_ERCC))
          axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
          abline(h = threshold_line_ercc, col = "red")
          dev.off()
        }
      )

      # Descargar el gráfico de Median
      output$download_med_plot <- downloadHandler(
        filename = function() { "median_plot.pdf" },
        content = function(file) {
          pdf(file)
          plot(ratios$median, xlab = "", ylab = "Median", col = colores_med,
               xaxt = "n", pch = 19, main = "Median Values (QC6)", ylim = c(0, max_value_med))
          axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
          abline(h = threshold_line_median, col = "red")
          dev.off()
        }
      )



      # Create the heatmap with ggplot2
      a<- ggplot2::ggplot(bin_df_melted, ggplot2::aes(x = variable, y = Sample, fill = factor(value))) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_fill_manual(values = c("0" = "#FFF9D0", "1" = "red"), labels = c("OK", "Outlier")) +
        ggplot2::labs(x = "QC Metrics", y = " ", fill = "QC Status") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       axis.text.y = ggplot2::element_text(size = 7),
                       legend.position = "bottom")

      output$qc_plot_heatmap <- renderPlot({print(a) })
      output$download_qc_plot_heatmap <- downloadHandler(
        filename = function() {
          paste("QC_plots_heatmap.pdf")
        },
        content = function(file) {
          pdf(file)
          print(a)
          dev.off()}
      )

      rows_with_1 <- suppressWarnings(rownames(bin_matrix)[apply(bin_matrix, 1, any)])
      outliers<- rows_with_1
      output$outliers <- renderText({
        paste("Los outliers son los siguientes:", paste(outliers, collapse = "/"))
      })

      #############################fi de QC

    } else {
      outliers <- NULL
    }

    # if (!is.null(design_formula)) {
    #   if (!all(all.vars(design_formula) %in% colnames(col_data))) {
    #     stop("Variables in design_formula must be columns in col_data.")
    #   }
    # }
    if (design_formula == "") {
      output$error_design_formula <- renderText("Desing formula is necesary")
      return()
    }


    if (pattern == "") {

      output$filtering_data <- renderText({"FILTERING THE COUNT DATA. DELETING THE PROVES..."})
      counts_data <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
    }

    # Removing outliers if specified
    if (remove_outliers && !is.null(outliers) && length(outliers) > 0) {
      output$removing_outliers <- renderText({"REMOVING OUTLIERS..."})
      cat("REMOVING OUTLIERS...\n")
      counts_filtered <- counts_data[, !colnames(counts_data) %in% outliers]
      AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
    } else {
      counts_filtered <- counts_data
      AnnotData <- col_data
    }
    ############################
    #####
    ####
    ####
    if (DEA) {
       output$missatge2 <- renderText({"STARTING THE DIFERENTIAN EXPRESSION ANALYSIS."})
       cat("STARTING THE DIFERENTIAN EXPRESSION ANALYSIS.\n")


       design_formul <- as.formula(paste("~", " ", design_formula))


       ### Variables should not have spaces
       colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))

       ### CHECK IF COLNAMES OF SAMPLE IDs IN col_data ARE THE NAMES OF COL IN COUNT DATA
       col_data <- AnnotData[order(AnnotData$id), ]
       counts_filtered <- counts_filtered[, order(colnames(counts_filtered))]

       # Create DESeqDataSet object
       rownames(col_data)<- col_data$id  #important that the columns is called id
       col_data[[contrast[1]]] <- as.factor(col_data[[contrast[1]]])
       counts_filtered<- as.data.frame(counts_filtered)

       dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = design_formul)
       n_genes <- nrow(DESeq2::counts(dds))
       n_subj <- ncol(DESeq2::counts(dds))
       zero_threshold <- ceiling(n_subj * 0.2)
       keep_genes <- rowSums(DESeq2::counts(dds) == 0) <= zero_threshold
       smallest_group_size <- ceiling(n_subj * percentage_gene)
       keep_genes <- keep_genes & (rowSums(DESeq2::counts(dds) >= threshold_gene) >= smallest_group_size)
       # Apply gene filters to the DESeq2 object
       dds <- dds[keep_genes, ]
       output$dds_missatge <- renderText({"Preprocessed Data for Differential Analysis"})
       output$dds_output <- renderPrint({dds})

       # Perform DESeq2 analysis
       vsd <- DESeq2::vst(dds, blind = FALSE)
       output$results_VST <- renderText({"RESULTS OF VST NORMALIZATION"})
        output$vsd_output <- renderPrint({vsd})
       dds <- DESeq2::DESeq(dds)

       # Generate heatmap if generate_heatmap is TRUE
       if (generate_heatmap) {
         output$missatge3 <- renderText({"HEATMAP"})

         # Seleccionar las filas con los 500 genes más expresados
         selecto <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)), decreasing = TRUE)[1:500]

         # Crear un data frame para las anotaciones de columnas
         df <- as.data.frame(SummarizedExperiment::colData(dds)[, heatmap_columns])
             # Generar el heatmap con nombres de columnas (muestras) y nombres de filas si se desea


         lola<-pheatmap::pheatmap(SummarizedExperiment::assay(vsd)[selecto,],
                               cluster_rows = FALSE,
                               show_rownames = FALSE,
                               show_colnames = TRUE,
                               cluster_cols = TRUE,
                               annotation_col = df)

         output$DEA_heatmap <- renderPlot({print(lola)})
         output$DEA_heatmap_download <- downloadHandler(
           filename = function() {
             paste("DEA_heatmap.pdf")
           },
           content = function(file) {
             pdf(file)
             print(lola)
             dev.off()}
         )

       }

       # Results contrast
       output$missatge_contrast <- renderText({"GENERATING CONTRAST RESULTS"})
       res <- DESeq2::results(dds, contrast = contrast, cooksCutoff = TRUE)
       output$summary_results<- renderText({"SUMMARY OF RESULT OF THE CONTRAST"})
       output$printsummary_results <- renderPrint({ print(summary(res))})

       output$top10_contrast<- renderText({"RESULTS OF THE CONTRAST (TOP 10 p-adj)"})
       output$printtop10_contrast <- renderPrint({print(head(res[order(res$padj), ], 10))})

    output$download_res <- downloadHandler(
      filename = function() { "Results_contrast_DEA.csv" },
      content = function(file) {
        write.csv(res, file, row.names = TRUE)
      }
    )


      dds <- DESeq2::estimateSizeFactors(dds)


      output$missatge_size_factors_plot <- renderText({"SIZE FACTORS PLOT"})
      output$size_factors_plot <- renderPlot({plot(DESeq2::sizeFactors(dds), colSums(DESeq2::counts(dds)),
           xlab = "Size Factors", ylab = "Column Sums of Counts",
           main = "Size Factors vs. Column Sums")
      text(DESeq2::sizeFactors(dds), colSums(DESeq2::counts(dds)), labels = colnames(DESeq2::counts(dds)), pos = 3, cex = 0.7)
      abline(lm(colSums(DESeq2::counts(dds)) ~ DESeq2::sizeFactors(dds) + 0))})

      output$download_DEA_plot <- downloadHandler(
      filename = function() { "DEA_plot.pdf" },
      content = function(file) {
        pdf(file)
        plot(DESeq2::sizeFactors(dds), colSums(DESeq2::counts(dds)),
               xlab = "Size Factors", ylab = "Column Sums of Counts",
               main = "Size Factors vs. Column Sums")
        text(DESeq2::sizeFactors(dds), colSums(DESeq2::counts(dds)), labels = colnames(DESeq2::counts(dds)), pos = 3, cex = 0.7)
        abline(lm(colSums(DESeq2::counts(dds)) ~ DESeq2::sizeFactors(dds) + 0))
        dev.off()
      }
    )
      vsd_cor <- cor(SummarizedExperiment::assay(vsd))
      sample_ids <- vsd$id
      annotation_col <- data.frame(SampleID = sample_ids)
      rownames(annotation_col) <- colnames(vsd_cor)

      aa<- pheatmap::pheatmap(vsd_cor,
                             main = "Sample-to-Sample Correlation Heatmap",
                             display_numbers = FALSE,
                             annotation_col = annotation_col,
                             annotation_legend = FALSE)
      output$correlation_heatmap_mes <- renderText({"CORRELATION HEATMAP"})
      output$correlation_heatmap <- renderPlot({print(aa)})

      output$download_correlation_heatmap <- downloadHandler(
        filename = function() { "Sample_Correlation_Heatmap.pdf" },
        content = function(file) {
          pdf(file)
          print(aa)
          dev.off()
        }
      )

      # Calcular la distancia de Poisson
      pois_distance <- PoiClaClu::PoissonDistance(t(DESeq2::counts(dds, normalized = TRUE)))
      samplePoisDistMatrix <- as.matrix(pois_distance$dd)
      sample_ids <- dds$id
      colnames(samplePoisDistMatrix) <- sample_ids
      rownames(samplePoisDistMatrix) <- sample_ids
      colors <- grDevices::colorRampPalette(c("white", "#4793AF", "#013649"))(255)

      output$missatge4 <- renderText({"GENERATING POISSON DISTANCES PLOT"})

      output$poison <- renderPlot({print(b)})

      output$download_poison <- downloadHandler(
        filename = function() { "DEA_poisson_distance.pdf" },
        content = function(file) {
          pdf(file, width = 10, height = 8)
          print(b)
          dev.off()
        }
      )

#### proba de posar en comptes de stats:: sese stats::
      sampleDists <- dist(t(SummarizedExperiment::assay(vsd)))
      sampleDistMatrix <- as.matrix(sampleDists)
      rownames(sampleDistMatrix) <- paste(vsd$id, sep = " - ")
      colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
      # colors <- grDevices::colorRampPalette(c("white", "#4793AF"))(255)
      output$missatge5 <- renderText({"GENERATING VSD DISTANCE HEATMAP"})

      c<- pheatmap::pheatmap(
        sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        #col = colors,
        fontsize = 8,
        main = "Sample Distance Heatmap",
        display_numbers = FALSE)

      output$distance_heatmap <- renderPlot({print(c)})

      output$download_distance_heatmap <- downloadHandler(
        filename = function() { "DEA_Sample_Distance_Heatmap.pdf" },
        content = function(file) {
          pdf(file, width = 10, height = 8)
          print(c)
          dev.off()
        }
      )

    output$missatge6 <- renderText({"GENERATING VSD BOXPLOT"})

      d<- graphics::boxplot(SummarizedExperiment::assay(vsd), las = 2, main = "vsd", cex.axis = 0.6)
      ## Renderizar el gráfico de boxplot
      output$vsd_boxplot <- renderPlot({
        graphics::boxplot(SummarizedExperiment::assay(vsd), las = 2, main = "vsd", cex.axis = 0.6)
      })

      res_sorted <- res[order(res$padj), ]
      top_genes_indices <- head(row.names(res_sorted), 10)

      output$topgenes <- renderText({
        paste("Plots for individual plot will be store in", paste(getwd(), "DEA_ gene_expression_plots.pdf", sep = "/"))
      })

      pdf("DEA_ gene_expression_plots.pdf")

      for (gene_index in top_genes_indices) {
        gen_a2m <- as.data.frame(SummarizedExperiment::assay(vsd)[gene_index, ])
        rownames(gen_a2m)
        gen_a2m$status <- col_data[[design_formula]]
        gen_a2m_ordered <- gen_a2m[order(gen_a2m$status), ]


        # Define colors based on status
        levels_design_formula <- unique(col_data[[design_formula]])

        num_colors <- length(levels_design_formula)
        palette <- colorRampPalette(c("red", "orange", "yellow", "green", "purple", "#4793AF"))(num_colors)
        colors <- palette

        group_colors <- colors[as.numeric(factor(col_data[[design_formula]], levels = levels_design_formula))]

        # Graficar
        plot(gen_a2m_ordered$`SummarizedExperiment::assay(vsd)[gene_index, ]`,
             xlab = "",
             ylab = gene_index,
             col = group_colors,
             xaxt = "n",
             pch = 19)
        axis(1, at = 1:nrow(gen_a2m_ordered), labels = rownames(gen_a2m_ordered), las = 2, cex.axis = 0.6)
        legend("topright", legend = levels_design_formula, fill = colors)
      }
      dev.off()




    # Función para descargar el boxplot
    output$download_vsd_boxplot <- downloadHandler(
      filename = function() { "VSD_Boxplot.pdf" },
      content = function(file) {
        pdf(file, width = 10, height = 8)
        graphics::boxplot(SummarizedExperiment::assay(vsd), las = 2, main = "vsd", cex.axis = 0.6)
        dev.off()
      }
    )
       output$missatge7 <- renderText({"GENERATING PCA PLOT"})

      normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
      pca_result <- prcomp(t(normalized_counts), scale. = TRUE)
      pca_data <- as.data.frame(pca_result$x)
      pca_data$Sample <- SummarizedExperiment::colData(dds)$id
      pca_data$Condition_Group <- SummarizedExperiment::colData(dds)[[design_formula]]

      color_palette <- c("#4793AF", "#E57373")
      names(color_palette) <- unique(pca_data$Condition_Group)

      pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = Condition_Group)) +
        ggplot2::geom_point(size = 3) +
        ggrepel::geom_text_repel(ggplot2::aes(label = Sample), size = 8, max.overlaps = 15) +  # Aumenta el tamaño de las etiquetas de muestra
        ggplot2::labs(title = "PCA of Normalized Counts",
                      x = "Principal Component 1",
                      y = "Principal Component 2") +
        ggplot2::theme_minimal() +
        ggplot2::scale_color_manual(values = color_palette) +
        ggplot2::theme(
          legend.position = "right",
          plot.title = ggplot2::element_text(size = 25),          # Título del gráfico
          axis.title.x = ggplot2::element_text(size = 25),        # Título del eje X
          axis.title.y = ggplot2::element_text(size = 25),        # Título del eje Y
          axis.text = ggplot2::element_text(size = 20),           # Texto de los ejes
          legend.title = ggplot2::element_text(size = 25),        # Título de la leyenda
          legend.text = ggplot2::element_text(size = 20)          # Texto de los elementos de la leyenda
        )

    output$pca_plot <- renderPlot({print(pca_plot)})
    output$download_pca_plot <- downloadHandler(
      filename = function() { "DEA_PCA.pdf" },
      content = function(file) {
        pdf(file, width = 10, height = 8)
        print(pca_plot)
        dev.off()
      }
    )
      output$missatge8 <- renderText({"GENERATING VOLCANO PLOT"})
      volcan<- EnhancedVolcano::EnhancedVolcano(res,
                                              lab = rownames(res),
                                              x = 'log2FoldChange',
                                              y = 'padj',
                                              pCutoff = 5e-2)
      output$volvano <- renderPlot({volcan})
      output$downloadvolvano <- downloadHandler(
        filename = function() {
          paste("DEA_volcano.pdf")
        },
        content = function(file) {
          pdf
          volcan
          dev.off()
        }
      )

      output$downloadvolvano <- downloadHandler(
        filename = function() {
          paste("DEA_volcano.pdf")
        },
        content = function(file) {
          pdf(file)
          print(volcan)

          dev.off()
        })


      output$missatge9 <- renderText({"GENERATING MA PLOT"})
      output$maPlot <- renderPlot({DESeq2::plotMA(res, main = "MA Plot of Results")})
    output$downloadMAPlot <- downloadHandler(
      filename = function() {
        paste("MA_plot.pdf")
      },
      content = function(file) {
        # Generar y guardar el gráfico en un archivo PDF
        pdf(file)
        DESeq2::plotMA(res, main = "MA Plot of Results")
        dev.off()
      }
    )


      output$missatge10 <- renderText({"GENERATING DISPERSION PLOT"})


    output$dispersion <- renderPlot({DESeq2::plotDispEsts(dds, main = "DIPERSION PLOT")})
    output$downloaddispersion<- downloadHandler(
      filename = function() {
        paste("DEA_dispersion_plot.pdf")
      },
      content = function(file) {
        pdf(file)
        DESeq2::plotDispEsts(dds, main = "DIPERSION PLOT")
        dev.off()
      })


     } else {cat("Skipping Diferential expresion analysis\n")}
    ############################
    #####
    ####
    ####
    if (GSEA) {
      if (!exists("res")) {
        cat("Warning: The variable 'res' does not exist. GSEA analysis cannot proceed without it.")
      }
      else if (!inherits(res, "DESeqResults")) {
        cat("Warning: The variable 'res' exists but is not of class 'DESeqResults'. GSEA analysis cannot proceed.")
      }
      else {
        cat("\033[33mSTARTING GSEA")

        # gseGO  Analysis
        cat("Performing gseGO Analysis")
        cat("Preparing gene list for GSEA")
        original_gene_list <- res$log2FoldChange
        names(original_gene_list) <- rownames(res)
        gene_list <- na.omit(original_gene_list)
        gene_list <- sort(gene_list, decreasing = TRUE)
        gse2 <- clusterProfiler::gseGO(geneList = gene_list, ont = "BP", keyType = "SYMBOL", nPermSimple = 500000,
                                       minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, eps = 0,
                                       OrgDb = org.Hs.eg.db::org.Hs.eg.db, pAdjustMethod = "bonferroni")

        # Dotplot for gseGO
        cat("Creating Dotplot for gseGO ")
        dotplot1 <- clusterProfiler::dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                                             title = "gseGO Enrichment Results: Pathways", color = "p.adjust", size = "Count")

        dotplot2 <- clusterProfiler::dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                                             title = "gseGO Enrichment Results: Pathways", color = "p.adjust", size = "Count") + ggplot2::facet_grid(.~.sign)

        # Emaplot for gseGO
        cat("Creating Emaplot for gseGO ")
        x2 <- enrichplot::pairwise_termsim(gse2)
        #emapplot1 <- enrichplot::emapplot(x2, max.overlaps = 70, min.segment.length = 0.3, point_size = 0.3, font.size = 5) +   ggplot2::ggtitle("Enrichment Map gseGO ")
        emapplot1 <- enrichplot::emapplot(x2) + ggplot2::ggtitle("Enrichment Map gseGO ")


        # Ridgeplot for gseGO
        cat("Creating Ridgeplot for gseGO ")
        ridgeplot1 <- enrichplot::ridgeplot(gse2)  +  ggplot2::labs(x = "gseGO enrichment distribution", font.size = 7) +  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))




        # Ridgeplot for gseGO
        ridgeplot1 <- enrichplot::ridgeplot(gse2)  +  ggplot2::labs(x = "gseGO enrichment distribution", font.size = 7) +  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))

        # Heatplot for gseGO
        heatplot1 <- enrichplot::heatplot(gse2, showCategory = 10) + ggplot2::ggtitle("gseGO Heatplot")

        # Treeplot for gseGO
        treeplot1 <- suppressWarnings(enrichplot::treeplot(x2) + ggplot2::ggtitle("gseGO Treeplot"))

        # Create gseaplot2 plots with titles
        a <- enrichplot::gseaplot2(gse2, geneSetID = 1, title = paste("GSEA Plot:", gse2$Description[1]))
        b <- enrichplot::gseaplot2(gse2, geneSetID = 1:5, pvalue_table = TRUE, title = "GSEA: Top 5 Gene Sets")
      output$gsea_dotplot <- renderPlot({print(dotplot1)})
      output$gsea_dotplot2 <- renderPlot({print(dotplot2)})
      output$gsea_emapplot <- renderPlot({print(emapplot1)})
      output$gsea_ridgeplot <- renderPlot({print(ridgeplot1)})
      output$gsea_heatplot <- renderPlot({print(heatplot1)})
      output$gsea_treeplot1 <- renderPlot({print(treeplot1)})
      output$gsea_a <- renderPlot({print(a)})
      output$gsea_b <- renderPlot({print(b)})
      output$download_dotplot1 <- downloadHandler(
        filename = function() { paste("gsea_dotplot.pdf") },
        content = function(file) {
          pdf(file)
          print(dotplot1)
          dev.off()
        }
      )
      output$download_dotplot2 <- downloadHandler(
        filename = function() { paste("gsea_dotplot2.pdf") },
        content = function(file) {
          pdf(file)
          print(dotplot2)
          dev.off()
        }
      )
      output$download_emapplot <- downloadHandler(
        filename = function() { paste("gsea_emapplot.pdf") },
        content = function(file) {
          pdf(file)
          print(emapplot1)
          dev.off()
        }
      )
      output$download_ridgeplot <- downloadHandler(
        filename = function() { paste("gsea_ridgeplot.pdf") },
        content = function(file) {
          pdf(file)
          print(ridgeplot1)
          dev.off()
        }
      )
      output$download_heatplot <- downloadHandler(
        filename = function() { paste("gsea_heatplot.pdf") },
        content = function(file) {
          pdf(file)
          print(heatplot1)
          dev.off()
        }
      )
      output$download_treeplot1 <- downloadHandler(
        filename = function() { paste("gsea_treeplot1.pdf") },
        content = function(file) {
          pdf(file)
          print(treeplot1)
          dev.off()
        }
      )

      output$download_a <- downloadHandler(
        filename = function() { paste("gsea_GSEA_Plot.pdf") },
        content = function(file) {
          pdf(file)
          print(a)
          dev.off()
        }
      )
      output$download_b <- downloadHandler(
        filename = function() { paste("gsea_GSEA_Top_5_Gene_Sets.pdf") },
        content = function(file) {
          pdf(file)
          print(b)
          dev.off()
        }
      )

        # KEGG Analysis
        ids <- clusterProfiler::bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb =  org.Hs.eg.db::org.Hs.eg.db)
        dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]
        df2 <- res[rownames(res) %in% dedup_ids$SYMBOL, ]
        df2$Y <- dedup_ids$ENTREZID
        kegg_gene_list <- df2$log2FoldChange
        names(kegg_gene_list) <- df2$Y
        kegg_gene_list <- na.omit(kegg_gene_list)
        kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

        kk2 <- clusterProfiler::gseKEGG(geneList = kegg_gene_list, organism = "hsa", minGSSize = 3, maxGSSize = 800,
                                        pvalueCutoff = 0.05, pAdjustMethod = "none", keyType = "ncbi-geneid", nPermSimple = 100000)

        # Dotplot for KEGG
        dotplot3 <- clusterProfiler::dotplot(kk2, showCategory = 10, title = "Enriched Pathways for KEGG", split = ".sign", font.size = 9) + ggplot2::facet_grid(.~.sign)

        # Emaplot for KEGG
        cat("a")
        x3 <- enrichplot::pairwise_termsim(kk2)
        emapplot2 <- enrichplot::emapplot(x3) + ggplot2::ggtitle("KEGG Enrichment Map")

        # Ridgeplot for KEGG
        ridgeplot2 <- enrichplot::ridgeplot(kk2) +  ggplot2::labs(x = "KEGG enrichment distribution", font.size = 6) +  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
        cat("b")

        # Heatplot for KEGG
        heatplot2 <- enrichplot::heatplot(kk2, showCategory = 10) + ggplot2::ggtitle("KEGG Heatplot")
        cat("c")

        # Treeplot for KEGG
        cat("Creating Treeplot for KEGG ")
        treeplot2 <- suppressWarnings(enrichplot::treeplot(x3) + ggplot2::ggtitle("KEGG Treeplot"))

        upset_plot <- enrichplot::upsetplot(kk2) + ggplot2::labs(title = "Up set plot for KEGG")
      output$gsea_dotplot3 <- renderPlot({print(dotplot3)})
      output$gsea_emapplot2 <- renderPlot({print(emapplot2)})
      output$gsea_ridgeplot2 <- renderPlot({print(ridgeplot2)})
      output$gsea_heatplot2 <- renderPlot({print(heatplot2)})
      output$gsea_treeplot2 <- renderPlot({print(treeplot2)})

      output$download_dotplot3 <- downloadHandler(
        filename = function() { paste("KEGG_dotplot.pdf") },
        content = function(file) {
          pdf(file)
          print(dotplot3)
          dev.off()
        }
      )
      output$download_emapplot2 <- downloadHandler(
        filename = function() { paste("KEGG_emapplot.pdf") },
        content = function(file) {
          pdf(file)
          print(emapplot2)
          dev.off()
        }
      )
      output$download_ridgeplot2 <- downloadHandler(
    filename = function() { paste("KEGG_ridgeplot.pdf") },
    content = function(file) {
    pdf(file)
    print(ridgeplot2)
    dev.off()
  }
)
output$download_heatplot2 <- downloadHandler(
  filename = function() { paste("KEGG_heatplot.pdf") },
  content = function(file) {
    pdf(file)
    print(heatplot2)
    dev.off()
  }
)
output$download_treeplot2 <- downloadHandler(
  filename = function() { paste("KEGG_treeplot1.pdf") },
  content = function(file) {
    pdf(file)
    print(treeplot2)
    dev.off()
  }
)


  # enrichGO Analysis
  cat("Performing GO Enrichment Analysis")
  sig_genes_df <- subset(res, padj < 0.05)

  if (nrow(sig_genes_df) > 0) {
    # Prepare the gene list for enrichment analysis
    genes <- sig_genes_df$log2FoldChange
    names(genes) <- rownames(sig_genes_df)

    # Perform GO enrichment analysis
    cat("Performing GO Enrichment Analysis")
    go_enrich <- clusterProfiler::enrichGO(
      gene = names(genes),
      universe = names(gene_list),
      OrgDb =  org.Hs.eg.db::org.Hs.eg.db,
      keyType = 'SYMBOL',
      readable = TRUE,
      ont = "BP",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.10
    )

          # Create a bar plot for significant GO terms
          cat("Creating Bar Plot for GO Enrichment")
          go_results <- go_enrich@result
          significant_terms <- go_results[go_results$qvalue < 0.05, ]
          significant_terms <- significant_terms[order(significant_terms$qvalue), ]

          # Check if there are significant terms to plot
          if (nrow(significant_terms) > 0) {
            bar_plot1 <- ggplot2::ggplot(significant_terms, ggplot2::aes(x = reorder(Description, -Count), y = Count)) +
              ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
              ggplot2::labs(title = "Significant GO Terms", x = "GO Term", y = "Count") +
              ggplot2::theme_minimal() +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

            output$go_bar_plot <- renderPlot({
              print(bar_plot1)
            })

            output$download_go_bar_plot <- downloadHandler(
              filename = function() { paste("GO_bar_plot.pdf") },
              content = function(file) {
                pdf(file)
                print(bar_plot1)
                dev.off()
              }
            )

          } else {
            cat("No significant GO terms found to plot.")
          }

          # # Save the GO enrichment results to a CSV file
          # write.csv(go_results, "go_results.csv")
          # cat("GO enrichment results saved to 'go_results.csv'")
        } else {
          cat("No significant genes found for GO enrichment analysis.")
        }


        if (nrow(sig_genes_df) > 0) {
          genes <- sig_genes_df$log2FoldChange
          names(genes) <- rownames(sig_genes_df)
          go_enrich <- clusterProfiler::enrichGO(gene = names(genes), universe = names(gene_list), OrgDb =  org.Hs.eg.db::org.Hs.eg.db,
                                                 keyType = 'SYMBOL', readable = TRUE, ont = "BP",
                                                 pvalueCutoff = 0.05, qvalueCutoff = 0.10)
          go_results <- go_enrich@result
          significant_terms <- go_results[go_results$qvalue < 0.05, ]
          significant_terms <- significant_terms[order(significant_terms$qvalue), ]

          bar_plot2 <- ggplot2::ggplot(significant_terms, ggplot2::aes(x = reorder(Description, -Count), y = Count)) +
            ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
            ggplot2::labs(title = "Significant GO Terms", x = "GO Term", y = "Count") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

          output$go_bar_plot2 <- renderPlot({
            print(bar_plot2)
          })

          output$download_go_bar_plot2 <- downloadHandler(
            filename = function() { paste("GO_bar_plot2.pdf") },
            content = function(file) {
              pdf(file)
              print(bar_plot2)
              dev.off()
            }
          )

        }
    }
    }
    ############################
    #####
    ####
    ####
    if (TME) {
      cat("STARTING TME")

      if (!is.null(pattern)) {
        # Remove outliers based on pattern
        filtered <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
      } else {
        # If no pattern is provided, do not filter based on pattern
        filtered <- counts_data
      }

      if (remove_outliers) {
        # Remove columns corresponding to outliers
        counts_data <- filtered[, !colnames(filtered) %in% outliers]
        AnnotData <- AnnotData[!AnnotData[["id"]] %in% outliers, ]
      } else {
        counts_data <- counts_data
        AnnotData <- AnnotData
      }

      # print(paste0("Inicial gene number: ", dim(counts_data)[1]))
      # ## Normalización TPM
      # cat("TPM normalization performed and stored on tpm_counts.csv")

       tpm_counts <- IOBR::count2tpm(counts_data,
                               idType = "Symbol",
                               org = "hsa",
                               source = "biomart")



      # Se almacenan los genes omitidos en la normalización TPM
      genes_omitidos <- base::setdiff(rownames(counts_data), rownames(tpm_counts))
      print(paste0("Number of genes omitted during TPM normalization due to their length not being available in Biomart:   ", dim(counts_data)[1]-dim(tpm_counts)[1]))
      tpm_counts <- as.data.frame(tpm_counts)
      output$genes_omitted <- renderText({
        paste0("Number of genes omitted during TPM normalization due to their length not being available in Biomart: ",
               dim(counts_data)[1] - dim(tpm_counts)[1])
      })


      # Mostrar el número inicial de genes
      output$initial_message <- renderText({
        paste0("Inicial gene number: ", dim(counts_data)[1])
      })

      # Mostrar el mensaje de normalización
      output$normalization_message <- renderText({
        "TPM normalization performed and stored on tpm_counts.csv"
      })

      # Mostrar el head de tpm_counts como tabla
      output$tpm_table <- renderTable({
        head(tpm_counts)
      })

      # Descargar el archivo tpm_counts.csv
      output$download_tpm <- downloadHandler(
        filename = function() {
          paste("tpm_counts.csv")
        },
        content = function(file) {
          write.csv(tpm_counts, file)
        }
      )

      if (!is.null(DEA) && DEA) {
        cat("We are going to use information from DEA\n")
        dds <- DESeq2::estimateSizeFactors(dds)
        normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
      } else {
        cat("Performing normalization.\n")
        design_formul <- as.formula("~ 1")
        colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))
        AnnotData <- AnnotData[order(AnnotData$id), ]
        counts_data <- counts_data[, order(colnames(counts_data))]
        if (!identical(colnames(counts_data), AnnotData$id)) {
          stop("Column names of counts_data and IDs in AnnotData do not match.\n")
        }
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_data, colData = AnnotData, design = design_formul)
        dds <- DESeq2::estimateSizeFactors(dds)
        normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
      }

      ## Deconvolución
      imm_epic <- immunedeconv::deconvolute(tpm_counts, method = "epic")
      imm_qti <- immunedeconv::deconvolute(tpm_counts, method = "quantiseq")
      imm_xcell <- immunedeconv::deconvolute(tpm_counts, method = "xcell")
      cat("results of the TME will be stored in imm_epic.csv, imm_qti.csv and imm_xcell.csv \n")
      # write.csv(imm_epic, file = "imm_epic.csv")
      # write.csv(imm_qti, file = "imm_qti.csv")
      # write.csv(imm_xcell, file = "imm_xcell.csv")
      output$cat_message <- renderText({
        "Results of the TME will be stored in imm_epic.csv, imm_qti.csv and imm_xcell.csv"
      })

      # Mostrar las primeras filas de los resultados como tabla
      output$imm_epic_table <- renderTable({head(imm_epic)})

      output$imm_qti_table <- renderTable({head(imm_qti)})

      output$imm_xcell_table <- renderTable({head(imm_xcell)})

      # Permitir la descarga de los archivos CSV
      output$download_epic <- downloadHandler(
        filename = function() { "TMEimm_epic.csv" },
        content = function(file) {
          write.csv(imm_epic, file)
        }
      )

      output$download_qti <- downloadHandler(
        filename = function() { "TMEimm_qti.csv" },
        content = function(file) {
          write.csv(imm_qti, file)
        }
      )

      output$download_xcell <- downloadHandler(
        filename = function() { "TMEimm_xcell.csv" },
        content = function(file) {
          write.csv(imm_xcell, file)
        }
      )

      # Transponer y cambiar colnames
      std.im.df <- function(imm_df){
        imm_df <- as.data.frame(t(imm_df))
        celltype_imm <- imm_df[1,]
        imm_df <- imm_df[-1,]
        colnames(imm_df) <- celltype_imm
        imm_df[-1,]
        rn_imm <- rownames(imm_df)
        imm_df <- as.data.frame(sapply(imm_df, as.numeric))
        rownames(imm_df) <- rn_imm
        return(imm_df)
      }

      imm_epic <- std.im.df(imm_epic)
      imm_qti <- std.im.df(imm_qti)
      imm_xcell <- std.im.df(imm_xcell)

      # Se incluye la variable AnnotData$design_formula en los dataframes
      cat("Are they in the same order?\n")

      rownames(AnnotData)<- AnnotData$id
      AnnotData <- AnnotData[order(rownames(AnnotData)), ]
      imm_epic <- imm_epic[order(rownames(imm_epic)), ]
      imm_qti <- imm_qti[order(rownames(imm_qti)), ]
      imm_xcell <- imm_xcell[order(rownames(imm_xcell)), ]
      AnnotData$id
      rownames(imm_epic)
      rownames(imm_qti)
      rownames(imm_xcell)

      cat("Have to be true.\n")
      cat("EPIC\n")
      print(all(rownames(AnnotData)==rownames(imm_epic)))
      cat("qti\n")
      print(all(rownames(AnnotData)==rownames(imm_qti)))
      cat("xcell\n")
      print(all(rownames(AnnotData)==rownames(imm_xcell)))


      imm_epic[[design_formula]] <- factor(AnnotData[[design_formula]])
      imm_qti[[design_formula]] <- factor(AnnotData[[design_formula]])
      imm_xcell[[design_formula]] <- factor(AnnotData[[design_formula]])

      replace_space_with_underscore <- function(df) {
        colnames(df) <- gsub(" ", "_", colnames(df))
        return(df)
      }
      imm_epic <- replace_space_with_underscore(imm_epic)
      imm_qti <- replace_space_with_underscore(imm_qti)
      imm_xcell <- replace_space_with_underscore(imm_xcell)

      ############# imm_epic
      design_formula_sym <- rlang::sym(design_formula)
      # Vector with column names except the last one
      column_names <- names(imm_epic)[-ncol(imm_epic)]
      # Loop to generate plots for each column

      pdf("TME_plots_imm_EPIC.pdf")
      for (col_name in column_names) {
        # Calculate means by design formula
        grouped_data <- dplyr::group_by(imm_epic, !!design_formula_sym)
        means_df <- dplyr::summarize(grouped_data, mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

        # Perform t-test or ANOVA
        formula <- as.formula(paste0("`", col_name, "` ~ ", rlang::as_label(design_formula_sym)))
        p_value <- tryCatch({
          if (dplyr::n_distinct(imm_epic[[rlang::as_label(design_formula_sym)]]) == 2) {
            t.test(formula, data = imm_epic)$p.value
          } else {
            summary(stats::aov(formula, data = imm_epic))[[1]]$`Pr(>F)`[1]
          }
        }, error = function(e) {
          NA  # In case of error, return NA for p-value
        })
        # Generate plot
        plot <- ggplot2::ggplot(imm_epic, ggplot2::aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                                       fill = !!design_formula_sym, color = !!design_formula_sym)) +
          ggplot2::geom_jitter(alpha = 1, width = 0.3, height = 0) +
          ggplot2::geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
          ggplot2::geom_point(data = means_df, ggplot2::aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                              shape = 22, color = "black", size = 3, stroke = 1.5,
                              show.legend = F) +
          ggplot2::labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_epic") +
          ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1),
                                      limits = c(0, NA)) +
          ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
          ggplot2::annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                            hjust = 1.1, vjust = 1.1, size = 5, color = "red")
        print(plot)
      }
      dev.off()



      ############# imm_qti
      # Vector with column names except the last one
      column_names <- names(imm_qti)[-ncol(imm_qti)]

      pdf("TME_plots_imm_qti.pdf")
      for (col_name in column_names) {
        # Calculate means by design formula
        grouped_data <- dplyr::group_by(imm_qti, !!design_formula_sym)
        means_df <- dplyr::summarize(grouped_data, mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

        # Perform t-test or ANOVA
        formula <- as.formula(paste0("`", col_name, "` ~ ", rlang::as_label(design_formula_sym)))
        p_value <- tryCatch({
          if (dplyr::n_distinct(imm_qti[[rlang::as_label(design_formula_sym)]]) == 2) {
            t.test(formula, data = imm_qti)$p.value
          } else {
            summary(stats::aov(formula, data = imm_qti))[[1]]$`Pr(>F)`[1]
          }
        }, error = function(e) {
          NA  # In case of error, return NA for p-value
        })

        # Generate plot
        plot <- ggplot2::ggplot(imm_qti, ggplot2::aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                                      fill = !!design_formula_sym, color = !!design_formula_sym)) +
          ggplot2::geom_jitter(alpha = 1, width = 0.3, height = 0) +
          ggplot2::geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
          ggplot2::geom_point(data = means_df, ggplot2::aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                              shape = 22, color = "black", size = 3, stroke = 1.5,
                              show.legend = F) +
          ggplot2::labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_qti") +
          ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1),
                                      limits = c(0, NA)) +
          ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
          ggplot2::annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                            hjust = 1.1, vjust = 1.1, size = 5, color = "red")
        print(plot)
      }
      dev.off()

      ############# imm_xcell
      # Vector with column names except the last one
      column_names <- names(imm_xcell)[-ncol(imm_xcell)]
      pdf("plots_imm_xcell.pdf")
      for (col_name in column_names) {
        # Calculate means by design formula
        grouped_data <- dplyr::group_by(imm_xcell, !!design_formula_sym)
        means_df <- dplyr::summarize(grouped_data, mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

        # Perform t-test or ANOVA
        formula <- as.formula(paste0("`", col_name, "` ~ ", rlang::as_label(design_formula_sym)))
        p_value <- tryCatch({
          if (dplyr::n_distinct(imm_xcell[[rlang::as_label(design_formula_sym)]]) == 2) {
            t.test(formula, data = imm_xcell)$p.value
          } else {
            summary(stats::aov(formula, data = imm_xcell))[[1]]$`Pr(>F)`[1]
          }
        }, error = function(e) {
          NA  # In case of error, return NA for p-value
        })
        # Generate plot
        plot <- ggplot2::ggplot(imm_xcell, ggplot2::aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                                        fill = !!design_formula_sym, color = !!design_formula_sym)) +
          ggplot2::geom_jitter(alpha = 1, width = 0.3, height = 0) +
          ggplot2::geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
          ggplot2::geom_point(data = means_df, ggplot2::aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                              shape = 22, color = "black", size = 3, stroke = 1.5,
                              show.legend = F) +
          ggplot2::labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_xcell") +
          ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1),
                                      limits = c(0, NA)) +
          ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
          ggplot2::annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                            hjust = 1.1, vjust = 1.1, size = 5, color = "red")
        print(plot)
      }
      dev.off()

      output$TME <- renderText({
        paste("Other Plots saved in ", paste(getwd(), sep = "/"))
      })



      # Composición celular del TME i por grupo

      # Function to generate the bar plot for each sample
      plot_bar <- function(df, paleta, titulo, legend.position) {
        # Convertir las filas en una columna llamada "Sample"
        df <- tibble::rownames_to_column(df, var = "Sample")
        df <- tidyr::pivot_longer(
          df, cols = colnames(df)[2:(ncol(df) - 1)],
          names_to = "Cell_Type", values_to = "Value")
        df <- dplyr::mutate(df,
                            Sample = factor(Sample, levels = rev(unique(Sample))),
                            Cell_Type = factor(Cell_Type, levels = rev(unique(Cell_Type))))

        p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Value, fill = Cell_Type)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::labs(title = titulo,
                        x = "Samples",
                        y = "Cell Fraction (%)") +
          ggplot2::coord_flip() +
          ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
          ggplot2::scale_fill_manual(values = paleta) +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = legend.position,
                         axis.text.y = ggplot2::element_text(size = 5)) +
          ggplot2::scale_y_continuous(labels = scales::percent)

        return(p)
      }

      plot_bar_group <- function(df, paleta, titulo, design_formula, legend_position = "right") {
        suppressWarnings({
          design_formula_sym <- rlang::sym(design_formula)
          niveles_tipo_cel <- colnames(df)[1:(ncol(df) - 1)]
          df_rownames <- tibble::rownames_to_column(df, var = "Sample")
          df_long <- tidyr::pivot_longer(
            df_rownames,
            cols = niveles_tipo_cel,
            names_to = "Cell_Type",
            values_to = "Value")
          df_grouped <- dplyr::group_by(df_long, !!design_formula_sym, Cell_Type)
          df_summarised <- dplyr::summarise(df_grouped,
                                            Average = mean(Value, na.rm = TRUE),
                                            .groups = "drop")
          df_ungrouped <- dplyr::ungroup(df_summarised)
          promedios <- dplyr::mutate(df_ungrouped,
                                     !!design_formula_sym := factor(!!design_formula_sym, levels = rev(unique(!!design_formula_sym))),
                                     Cell_Type = factor(Cell_Type, levels = rev(niveles_tipo_cel)))
          p <- ggplot2::ggplot(promedios, ggplot2::aes(x = !!design_formula_sym, y = Average, fill = Cell_Type)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::labs(title = titulo,
                          x = "Samples",
                          y = "Cell Fraction (%)") +
            ggplot2::coord_flip() +
            ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
            ggplot2::scale_fill_manual(values = paleta) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = legend_position,
                           axis.text.y = ggplot2::element_text(size = 5)) +
            ggplot2::scale_y_continuous(labels = scales::percent)
          return(p)
        })
      }

      # Function to combine both plots into one
      plot_combined <- function(df, paleta, titulo_individual, titulo_grupo, design_formula, legend_position = "right") {
        p1 <- plot_bar(df, paleta, titulo_individual, legend_position)
        p2 <- plot_bar_group(df, paleta, titulo_grupo, design_formula, legend_position)

        combined_plot <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1))
        return(combined_plot)
      }

      # Example usage with extended palette for larger datasets
      paleta_imm <- c("grey95","#FB8072","#FFED6F","#6F6C87","#94DFD1","#FDB462", "#B3DE69", "#FFB1D9")
      paleta_qti <- c("grey95","#B3DE69","#6F6C87","#94DFD1","#FDB462", "#FB8072","#FFFFB3","#8BB07A","#FFED6F","#80B1D3","#FFB1D9")
      paleta_extendida <- c("grey95","#2CA02C","#E6F8E0","#7F7F7F","#FF8000","#DF7401","#6F6C80","#F5A9E1",
                            "#F7D358","#6F6C99","#FB8072","#FFFFB3","#F5A9F2",
                            "#8BB07A","#E377C2","#FE2E2E","#FFED6F","#1F77B4","#80B1D3","#0489B1","#FF0000","#8BB07A",
                            "#0000FF","#D6616B","#58FA82","#98DF8A",
                            "#8C564B", "#C49C94","#CEE3F6","#E0F8EC","#58FAF4", "#ADD8E6","#9EDAE5","#94DFD1","#FFA500",
                            "#FF7F0E","#FDB462","#FFB1D9","#B3DE69")
      # Assuming 'imm_epic', 'imm_qti', and 'imm_xcell' dataframes are already loaded
      combined_plot_EPIC <- plot_combined(imm_epic, paleta_imm, "EPIC Individual", "EPIC Average", design_formula, "right")
      combined_plot_quanTIseq <- plot_combined(imm_qti, paleta_qti, "quanTIseq Individual", "quanTIseq Average", design_formula, "right")
      plot_bar <- function(df, paleta, titulo, legend.position) {
        # Convertir las filas en una columna llamada "Sample"
        df <- tibble::rownames_to_column(df, var = "Sample")
        df <- tidyr::pivot_longer(
          df, cols = colnames(df)[2:(ncol(df) - 1)],
          names_to = "Cell_Type", values_to = "Value")
        df <- dplyr::mutate(df,
                            Sample = factor(Sample, levels = rev(unique(Sample))),
                            Cell_Type = factor(Cell_Type, levels = rev(unique(Cell_Type))))

        p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Value, fill = Cell_Type)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::labs(title = titulo,
                        x = "Samples",
                        y = "Enrichment Scores") +
          ggplot2::coord_flip() +
          ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
          ggplot2::scale_fill_manual(values = paleta) +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = legend.position,
                         axis.text.y = ggplot2::element_text(size = 5)) +
          ggplot2::scale_y_continuous(labels = scales::percent)

        return(p)
      }

      plot_bar_group <- function(df, paleta, titulo, design_formula, legend_position = "right") {
        suppressWarnings({
          design_formula_sym <- rlang::sym(design_formula)
          niveles_tipo_cel <- colnames(df)[1:(ncol(df) - 1)]
          df_rownames <- tibble::rownames_to_column(df, var = "Sample")
          df_long <- tidyr::pivot_longer(
            df_rownames,
            cols = niveles_tipo_cel,
            names_to = "Cell_Type",
            values_to = "Value")
          df_grouped <- dplyr::group_by(df_long, !!design_formula_sym, Cell_Type)
          df_summarised <- dplyr::summarise(df_grouped,
                                            Average = mean(Value, na.rm = TRUE),
                                            .groups = "drop")
          df_ungrouped <- dplyr::ungroup(df_summarised)
          promedios <- dplyr::mutate(df_ungrouped,
                                     !!design_formula_sym := factor(!!design_formula_sym, levels = rev(unique(!!design_formula_sym))),
                                     Cell_Type = factor(Cell_Type, levels = rev(niveles_tipo_cel)))
          p <- ggplot2::ggplot(promedios, ggplot2::aes(x = !!design_formula_sym, y = Average, fill = Cell_Type)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::labs(title = titulo,
                          x = "Samples",
                          y = "Enrichment Scores") +
            ggplot2::coord_flip() +
            ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
            ggplot2::scale_fill_manual(values = paleta) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = legend_position,
                           axis.text.y = ggplot2::element_text(size = 5)) +
            ggplot2::scale_y_continuous(labels = scales::percent)
          return(p)
        })
      }

      # Function to combine both plots into one
      plot_combined <- function(df, paleta, titulo_individual, titulo_grupo, design_formula, legend_position = "right") {
        p1 <- plot_bar(df, paleta, titulo_individual, legend_position)
        p2 <- plot_bar_group(df, paleta, titulo_grupo, design_formula, legend_position)

        combined_plot <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1))
        return(combined_plot)
      }
      combined_plot_xCell <- plot_combined(imm_xcell, paleta_extendida, "xCell Individual", "xCell Average", design_formula, "right")
      combined_plot_xCell <- plot_combined(imm_xcell, paleta_extendida, "xCell Individual", "xCell Average", design_formula, "right")

      #pdf("plot_cell_fraction_Average_cell_fraction_EPIC.pdf", width = 11, height = 14)
      #print(combined_plot_EPIC)
      #dev.off()
      #pdf("plot_cell_fraction_Average_cell_fraction_quanTIseq.pdf", width = 11, height = 14)
      #print(combined_plot_quanTIseq)
      #dev.off()
      #pdf("plot_cell_fraction_Average_cell_fraction_xCell.pdf", width = 11, height = 14)
      #print(combined_plot_xCell)
      #dev.off()

      # Per a la visualització dels gràfics
      output$plot_cell_fraction_EPIC <- renderPlot({
        print(combined_plot_EPIC)
      })

      output$plot_cell_fraction_quanTIseq <- renderPlot({
        print(combined_plot_quanTIseq)
      })

      output$plot_cell_fraction_xCell <- renderPlot({
        print(combined_plot_xCell)
      })

      # Per a la descàrrega dels gràfics en format PDF
      output$download_plot_EPIC <- downloadHandler(
        filename = function() { "plot_cell_fraction_Average_cell_fraction_EPIC.pdf" },
        content = function(file) {
          pdf(file, width = 11, height = 14)
          print(combined_plot_EPIC)
          dev.off()
        }
      )

      output$download_plot_quanTIseq <- downloadHandler(
        filename = function() { "plot_cell_fraction_Average_cell_fraction_quanTIseq.pdf" },
        content = function(file) {
          pdf(file, width = 11, height = 14)
          print(combined_plot_quanTIseq)
          dev.off()
        }
      )

      output$download_plot_xCell <- downloadHandler(
        filename = function() { "plot_cell_fraction_Average_cell_fraction_xCell.pdf" },
        content = function(file) {
          pdf(file, width = 11, height = 14)
          print(combined_plot_xCell)
          dev.off()
        }
      )

      ########################
      trans_formato_largo <- function(df, design_formula_sym) {
        df_largo <- tidyr::pivot_longer(df, cols = -dplyr::all_of(design_formula_sym), names_to = "Cell_Type",
                                        values_to = "Fraction")
        return(df_largo)
      }

      # Function to perform Shapiro-Wilk normality test
      prueba_norm <- function(df, design_formula_sym) {
        df_largo <- trans_formato_largo(df, design_formula_sym)

        resultados_normalidad <- dplyr::group_by(df_largo, Cell_Type)
        resultados_normalidad <- dplyr::summarise(resultados_normalidad,
                                                  shapiro_test = list(stats::shapiro.test(Fraction)))
        resultados_normalidad <- dplyr::mutate(resultados_normalidad,
                                               p.value = purrr::map_dbl(shapiro_test, "p.value"))


        return(resultados_normalidad)
      }

      # Function to check if all values in a column are the same
      check_column_equal <- function(column) {
        all_equal <- length(unique(column)) == 1
        return(all_equal)
      }

      # Function to filter valid columns
      filter_valid_columns <- function(df, design_formula_sym, dataset_name) {
        cat("Checking if All Values in Each Column Are Equal in dataset:", dataset_name, "")

        equality_results <- base::sapply(df, check_column_equal)
        unequal_columns_count <- base::sum(!equality_results)

        if (unequal_columns_count > 0) {
          cat("In dataset", dataset_name, "the following columns have the same value in all rows and will be excluded from the analysis")
          print(names(df)[equality_results])
        } else {
          cat("All columns in dataset", dataset_name, "have variability. No columns to exclude.")
        }

        valid_columns <- !equality_results
        df_filtered <- df[, valid_columns]

        return(df_filtered)
      }

      # Function to perform parametric tests (t-test and ANOVA)
      perform_parametric_tests <- function(df, design_formula_sym) {
        df_largo <- trans_formato_largo(df, design_formula_sym)
        cell_types <- unique(df_largo$Cell_Type)

        results <- list()

        for (cell_type in cell_types) {

          df_cell_type <- dplyr::filter(df_largo, Cell_Type == cell_type)
          groups <- unique(df_cell_type[[design_formula_sym]])


          if (length(groups) == 2) {
            # Perform t-test
            t_test_result <- t.test(Fraction ~ df_cell_type[[design_formula_sym]], data = df_cell_type)
            results[[cell_type]] <- list(Test = "t-test", p.value = t_test_result$p.value)
          } else {
            # Perform ANOVA
            anova_result <- stats::aov(Fraction ~ df_cell_type[[design_formula_sym]], data = df_cell_type)
            p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]
            results[[cell_type]] <- list(Test = "ANOVA", p.value = p_value)
          }
        }

        results_df <- dplyr::bind_rows(lapply(names(results), function(cell_type) {
          result <- results[[cell_type]]
          tibble::tibble(Cell_Type = cell_type, Test = result$Test, p.value = result$p.value)
        }))

        return(results_df)
      }

      # Apply the filtering and normality test for each dataset
      cat("Apply the filtering and normality test for each dataset\n")
      imm_qti_filtered <- filter_valid_columns(imm_qti, design_formula_sym, "quanTIseq")
      imm_epic_filtered <- filter_valid_columns(imm_epic, design_formula_sym, "EPIC")
      imm_xcell_filtered <- filter_valid_columns(imm_xcell, design_formula_sym, "xcell")

      # Perform the normality test
      if (ncol(imm_qti_filtered) > 1) {
        cat("Performing Shapiro-Wilk test for quanTIseq\n")
        resultados_norm_imm_qti <- prueba_norm(imm_qti_filtered, design_formula_sym)
        print(resultados_norm_imm_qti)
        output$resultados_norm_imm_qtiA <- renderPrint({
          cat("Performing Shapiro-Wilk test for quanTIseq\n")
          print(resultados_norm_imm_qti)
        })
        # Perform parametric tests if normality is satisfied
        cat("Performing parametric tests for quanTIseq\n")
        parametric_results_imm_qti <- perform_parametric_tests(imm_qti_filtered, design_formula_sym)
        print(parametric_results_imm_qti)
        output$parametric_results_imm_qtiA <- renderPrint({
          cat("Performing parametric tests for quanTIseq\n")
          print(parametric_results_imm_qti)
        })
      } else {
        cat("Can't perform Shapiro-Wilk test for quanTIseq. No valid columns available.\n")
      }

      if (ncol(imm_epic_filtered) > 1) {
        cat("Performing Shapiro-Wilk test for EPIC\n")
        resultados_norm_epic <- prueba_norm(imm_epic_filtered, design_formula_sym)
        print(resultados_norm_epic)
        output$resultados_norm_epicA <- renderPrint({
          cat("Performing Shapiro-Wilk test for EPIC\n")
          print(resultados_norm_epic)
        })

        # Perform parametric tests if normality is satisfied
        parametric_results_epic <- perform_parametric_tests(imm_epic_filtered, design_formula_sym)
        print(parametric_results_epic)
        output$parametric_results_epicA <- renderPrint({
          cat("Performing Shapiro-Wilk test for EPIC\n")
          print(parametric_results_epic)
        })
      } else {
        cat("Can't perform Shapiro-Wilk test for EPIC. No valid columns available.\n")
      }

      if (ncol(imm_xcell_filtered) > 1) {
        cat("Performing Shapiro-Wilk test for xcell\n")
        resultados_norm_xcell <- prueba_norm(imm_xcell_filtered, design_formula_sym)
        print(resultados_norm_xcell)
        output$resultados_norm_xcellA <- renderPrint({
          cat("Performing Shapiro-Wilk test for xcell\n")
          print(resultados_norm_xcell)
        })
        # Perform parametric tests if normality is satisfied
        cat("Performing parametric tests for xcell\n")
        parametric_results_xcell <- perform_parametric_tests(imm_xcell_filtered, design_formula_sym)
        output$parametric_results_xcellA <- renderPrint({
          cat("Performing Shapiro-Wilk test for xcell\n")
          print(parametric_results_xcell)
        })
      } else {
        cat("Can't perform Shapiro-Wilk test for xcell. No valid columns available.\n")
      }

      ############# Heatmaps
      # Transponer y estandarizar por filas
      h_imm_epic <- as.data.frame(t(imm_epic))
      ## Delete the row uncharacterized cell
      h_imm_epic <- head(h_imm_epic, -1)
      rn_himmepic <- rownames(h_imm_epic)
      # "Cancer associated fibroblast" will become CAFs
      rn_himmepic[2] <- "CAFs"
      cl_himmepic <- colnames(h_imm_epic)
      h_imm_epic <- apply(h_imm_epic, 2, as.numeric)
      h_imm_epic <- t(apply(h_imm_epic,1,scale))
      rownames(h_imm_epic) <- rn_himmepic
      colnames(h_imm_epic) <- cl_himmepic
      h_imm_epic <- as.data.frame(h_imm_epic)


      # Transponer y estandarizar por filas
      h_imm_qti <- as.data.frame(t(imm_qti))
      h_imm_qti <- head(h_imm_qti, -1)
      rn_himmqti <- rownames(h_imm_qti)
      rn_himmqti[7] <- "T cell CD4+"
      rn_himmqti[9] <- "T cell regulatory"
      cl_himmqti <- colnames(h_imm_qti)
      h_imm_qti <- apply(h_imm_qti, 2, as.numeric)
      h_imm_qti <- t(apply(h_imm_qti,1,scale))
      rownames(h_imm_qti) <- rn_himmqti
      colnames(h_imm_qti) <- cl_himmqti
      h_imm_qti <- as.data.frame(h_imm_qti)

      # Crear dataframe para heatmap
      h_imm_xcell <- as.data.frame(imm_xcell)
      # Transponer y estandarizar por filas
      h_imm_xcell <- as.data.frame(t(h_imm_xcell))
      h_imm_xcell <- head(h_imm_xcell, -4)

      # Poblaciones celulares no interesantes
      xcell_row_delete <- c('Common lymphoid progenitor', 'Common myeloid progenitor',
                            'Granulocyte-monocyte progenitor', 'Hematopoietic stem cell')

      h_imm_xcell <- dplyr::filter(h_imm_xcell, !rownames(h_imm_xcell) %in% xcell_row_delete)

      rn_himmxcell <- rownames(h_imm_xcell)
      cl_himmxcell <- colnames(h_imm_xcell)
      h_imm_xcell <- apply(h_imm_xcell, 2, as.numeric)
      h_imm_xcell <- t(apply(h_imm_xcell,1,scale))
      rownames(h_imm_xcell) <- rn_himmxcell
      colnames(h_imm_xcell) <- cl_himmxcell
      h_imm_xcell <- as.data.frame(h_imm_xcell)

      # Factorización de variables y renombrado en AnnotData
      AnnotData[[design_formula]] <- as.factor(AnnotData[[design_formula]])

      ################################### HEATMAP
      cat("GENERATING HEATMAPS")
      remove_nan_rows <- function(df) {
        # Remove rows with any NaN values
        df_cleaned <- df[!apply(df, 1, function(row) any(is.nan(row))), ]
        return(df_cleaned)
      }

      # Generar el heatmap para h_imm_epic
      combined_data <- remove_nan_rows(h_imm_epic)

      if (ncol(combined_data) < 3) {
        cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
      } else {
        if (any(is.na(combined_data))) {
          cat("Cannot generate the heatmap because there are NaN values in the data.\n")
        } else {
          col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

          if (ncol(combined_data) == nrow(col_ann_data)) {
            col_ann_data <- as.data.frame(col_ann_data)
          } else {
            stop("Dimensions of col_ann_data and combined_data do not match")
          }

          annotation_col <- as.data.frame(col_ann_data[, design_formula, drop = FALSE])
          rownames(annotation_col) <- colnames(combined_data)

          HeatmapEPIC <- pheatmap::pheatmap(
            as.matrix(combined_data),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            annotation_col = annotation_col,
            fontsize = 9,
            color = colorRampPalette(c("#4793AF", "white", "#013649"))(50),
            border_color = "grey60",
            main = "Heatmap EPIC",
            legend = TRUE,
            angle_col = 45,
            silent = TRUE
          )
        }
      }

      # Generar el heatmap para h_imm_qti
      combined_data <- remove_nan_rows(h_imm_qti)

      if (ncol(combined_data) < 3) {
        cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
      } else {
        if (any(is.na(combined_data))) {
          cat("Cannot generate the heatmap because there are NaN values in the data.\n")
        } else {
          col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

          if (ncol(combined_data) == nrow(col_ann_data)) {
            col_ann_data <- as.data.frame(col_ann_data)
          } else {
            stop("Dimensions of col_ann_data and combined_data do not match\n")
          }

          annotation_col <- as.data.frame(col_ann_data[, design_formula, drop = FALSE])
          rownames(annotation_col) <- colnames(combined_data)

          Heatmap_qti <- pheatmap::pheatmap(
            as.matrix(combined_data),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            annotation_col = annotation_col,
            fontsize = 9,
            color = colorRampPalette(c("#4793AF", "white", "#013649"))(50),
            border_color = "grey60",
            main = "Heatmap QTI",
            legend = TRUE,
            angle_col = 45,
            silent = FALSE
          )
        }
      }

      # Generar el heatmap para h_imm_xcell
      combined_data <- remove_nan_rows(h_imm_xcell)

      if (ncol(combined_data) < 3) {
        cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
      } else {
        if (any(is.na(combined_data))) {
          cat("Cannot generate the heatmap because there are NaN values in the data.\n")
        } else {
          col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

          if (ncol(combined_data) == nrow(col_ann_data)) {
            col_ann_data <- as.data.frame(col_ann_data)
          } else {
            stop("Dimensions of col_ann_data and combined_data do not match")
          }

          annotation_col <- as.data.frame(col_ann_data[, design_formula, drop = FALSE])
          rownames(annotation_col) <- colnames(combined_data)

          Heatmap_xcell <- pheatmap::pheatmap(
            as.matrix(combined_data),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            annotation_col = annotation_col,
            fontsize = 9,
            color = colorRampPalette(c("#4793AF", "white", "#013649"))(50),
            border_color = "grey60",
            main = "Heatmap XCELL",
            legend = TRUE,
            angle_col = 45,
            silent = TRUE
          )
        }
      }

      cat("Heatmap of qti, EPIC and xcell will be stored on plots_TME_heatmap.pdf\n")

      if (exists("Heatmap_qti")) {
        #pdf("plots_TME_Heatmap_qti.pdf", width = 15, height = 14)
        #print(Heatmap_qti)
        #dev.off()
        output$heatmap_qti <- renderPlot({print(Heatmap_qti)})
        output$download_heatmap_qti <- downloadHandler(
          filename = function() { paste("plots_TME_Heatmap_qti.pdf") },
          content = function(file) {
            pdf(file, width = 15, height = 14)
            if (exists("Heatmap_qti")) {
              print(Heatmap_qti)
            } else {
              cat("Heatmap_qti object does not exist.\n")
            }
            dev.off()
          }
        )
      } else {
        cat("Heatmap_qti object does not exist.\n")
      }

      if (exists("HeatmapEPIC")) {
        output$heatmap_EPIC <- renderPlot({print(HeatmapEPIC)})
        output$download_heatmap_EPIC <- downloadHandler(
          filename = function() { paste("plots_TME_HeatmapEPIC.pdf") },
          content = function(file) {
            pdf(file, width = 15, height = 14)
            if (exists("HeatmapEPIC")) {
              print(HeatmapEPIC)
            } else {
              cat("HeatmapEPIC object does not exist.\n")
            }
            dev.off()
          }
        )

        #pdf("plots_TME_HeatmapEPIC.pdf", width = 15, height = 14)
        #print(HeatmapEPIC)
        #dev.off()
      } else {
        cat("HeatmapEPIC object does not exist.\n")
      }
      if (exists("Heatmap_xcell")) {
        #pdf("plots_TME_Heatmap_xcell.pdf", width = 15, height = 14)
        #print(Heatmap_xcell)
        #dev.off()
        output$heatmap_xcell <- renderPlot({print(Heatmap_xcell)})
        output$download_heatmap_xcell <- downloadHandler(
          filename = function() { paste("plots_TME_Heatmap_xcell.pdf") },
          content = function(file) {
            pdf(file, width = 15, height = 14)
            if (exists("Heatmap_xcell")) {
              print(Heatmap_xcell)
            } else {
              cat("Heatmap_xcell object does not exist.\n")
            }
            dev.off()
          }
        )
      } else {
        cat("Heatmap_xcell object does not exist.\n")
      }


    }else {cat("TME analysis skipped.")}

    #####
    #####
    #####
    if (survival_analysis) {
      cat("STARTING SURVIVAL ANALYSIS")
      if (is.null(col_data[[time]]) || is.null(col_data[[variable_01]])) {
        stop("Variables for survival analysis are required.")
      }

      if (!is.null(pattern)) {
        # Remove outliers based on pattern
        filtered <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
      } else {
        filtered <- counts_data
      }

      if (remove_outliers) {
        counts_filtered <- filtered[, !colnames(filtered) %in% outliers]
        AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
      } else {
        counts_filtered <- counts_data
        AnnotData <- col_data
      }

      output$normalizing <- renderText({"Normalizing data..."})
      colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))
      col_data <- AnnotData[order(AnnotData$id), ]
      col_data<- as.data.frame(col_data)
      counts_filtered <- counts_filtered[, order(colnames(counts_filtered))]
      if (!identical(colnames(counts_filtered), col_data$id)) {
        stop("Column names of counts_filtered and IDs in col_data do not match.\n")
      }

      # Clean column names to avoid issues with special characters
      clean_column_names <- function(names) {
        gsub("[^[:alnum:]_]", "_", names)
      }

      # Create DESeqDataSet object
      if (!exists("dds")) {
        rownames(col_data) <- col_data$id
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = ~ 1)
        dds <- DESeq2::estimateSizeFactors(dds)
      }
      normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
      df_t <- t(normalized_counts)

      # Checking for duplicate columns
      output$checkduplicates <- renderText({"Checking for duplicate columns"})
      duplicated_columns <- colnames(df_t)[duplicated(colnames(df_t))]
      if (length(duplicated_columns) > 0) {
        output$duplicatesfound <- renderText({paste("Duplicate columns found:", paste(duplicated_columns, collapse = ", "))})
        } else {
         output$noduplicates <- renderText({"No duplicate columns."})
      }

      rownames(col_data) <- col_data$id
      ids_data <- rownames(df_t)

      # Filter data
      cat("Filtering data.")
      subset_data <- dplyr::filter(col_data, id %in% ids_data)
      df_ta <- as.data.frame(df_t)
      df_ta$id <- rownames(df_ta)

      # Selecting genes
      if (!is.null(genes_to_use)) {
        cat("Using provided genes\n")
        output$genes_used <- renderText({"Using provided genes"})
        top_genes <- genes_to_use
      } else if (exists("res")) {
        cat("Selecting TOP 10 genes with the lowest padj\n")
        output$genes_used2 <- renderText({"Selecting TOP 10 genes with the lowest padj"})

        output$top10 <- renderText({"Selecting TOP 10 genes with the lowest padj"})
        top_genes <- rownames(head(res[order(res$padj), ], 10))
        top_genes <- rownames(head(res[order(res$padj), ], 10))
      } else {
        stop("No genes provided or found in res object.")
      }
      selected_df_t <- df_t[, top_genes, drop = FALSE]
      selected_df_t <- as.data.frame(selected_df_t)
      selected_df_t$id <- rownames(selected_df_t)
      col_data$id <- rownames(col_data)
      merged_data <- merge(col_data, selected_df_t, by = "id")
      merged_data[[time]] <- as.numeric(merged_data[[time]])

      cat("Starting survival analysis.\n")
      output$startingsurvival <- renderText({"Starting survival analysis."})

      pdf("survival_analysis_plots_HTG_analyzer.pdf")
      top_genes_clean <- clean_column_names(top_genes)

      # Perform survival analysis for each gene
      for (i in top_genes_clean) {
        if (!is.numeric(merged_data[[i]])) {
          next
        }
        cat("Performing analysis for column:\033[0m ", i, "\n")
        # Perform MAXSTAT test
        merged_data$time <- merged_data[[time]]
        merged_data$variable_01 <- merged_data[[variable_01]]
        gene_column <- merged_data[[i]]    #get(i, merged_data)
        MAXSTAT <- maxstat::maxstat.test(survival::Surv(time, variable_01) ~ gene_column, data = merged_data,
                                         smethod = "LogRank", pmethod = "Lau92", iscores = TRUE, minprop = 0.45, maxprop = 0.55)
        cut.off <- MAXSTAT$estimate
        output$cutoffmessage <- renderText({
          paste("CUT OFF using MAXSTAT:", cut.off)
        })

        # Create a new variable based on the cutoff
        new_column_name <- paste0(i, "_mRNA_expression")
        merged_data[[new_column_name]] <- ifelse(gene_column > cut.off, "High", "Low")
        merged_data[[new_column_name]] <- factor(merged_data[[new_column_name]])

        # Fit survival model
        cat("Fitting survival model")
        surv_object <- survival::Surv(merged_data$time, merged_data$variable_01)
        surv_formula <- as.formula(paste("surv_object ~", new_column_name))

        fit1 <- survival::survfit(surv_formula, data = merged_data)

        # Summary of the fit
        cat("Summary of the fit")
        output$summaryfit <- renderText({"Summary of the fit"})
        print(summary(fit1))
        fit1_df <- data.frame(
          time = fit1$time,
          n_risk = fit1$n.risk,
          n_event = fit1$n.event,
          n_censor = fit1$n.censor,
          survival = fit1$surv,
          std_err = fit1$std.err,
          cumhaz = fit1$cumhaz,
          std_chaz = fit1$std.chaz,
          lower = fit1$lower,
          upper = fit1$upper
        )
        fit1_df$strata <- rep(names(fit1$strata), times = fit1$strata)
        merged_data_subset <- data.frame(
          time = merged_data$time,
          n_event = merged_data$variable_01,
          id = merged_data$id
        )
        fit1_df <- merge(fit1_df, merged_data_subset, by = c("time", "n_event"), all.x = TRUE)
        csv_filename <- paste0("Summary_of_the_fit_", new_column_name, ".csv")
        # write.csv(fit1_df, file = csv_filename, row.names = FALSE)
        output$fit1_df <- renderTable({head(fit1_df)})
        output$download_fit1_df <- downloadHandler(
          filename = function() {
            paste0("Summary_of_the_fit_", Sys.Date(), ".csv")
          },
          content = function(file) {
            write.csv(fit1_df, file, row.names = FALSE)
          }
        )


        # Log-rank test and p-value
        output$logrank <- renderText({"Performing log-rank test and obtaining p-value"})
        surv_diff <- survival::survdiff(surv_formula, data = merged_data)
        p_value <- 1 - stats::pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

        print(surv_diff)
        output$surv_diff_output <- renderPrint({
          print(surv_diff)
        })
        surv_diff_df <- data.frame(
          group = attr(surv_diff$n, "dimnames")[[1]],
          N = as.vector(surv_diff$n),
          Observed = surv_diff$obs,
          Expected = surv_diff$exp,
          `O-E^2/E` = (surv_diff$obs - surv_diff$exp)^2 / surv_diff$exp,
          `Chisq` = rep(surv_diff$chisq, 2),  # chi-cuadrado repetido dos veces
          p_value = rep(surv_diff$pvalue, 2)    # p-valor repetido dos veces
        )
        csv_filename <- paste0("surv_diff_summary_", new_column_name, ".csv")
        write.csv(surv_diff_df, file = csv_filename, row.names = FALSE)
        cat("P-value")
        print(p_value)

        # Generate Kaplan-Meier plot
        cat("Generating Kaplan-Meier plot")
        palette <- c("#9A3449", "#D4A8B1")
        plot(fit1, lty = 1, col = palette, lwd = 4, main = paste("Survival analysis for", i, "\n", "p-value =", format(p_value, digits = 3)))

        # Add a legend
        legend("topright",
               legend = c("High", "Low"),
               lty = 1,
               col = palette,
               lwd = 4)
        cat("Plots and .csv saved in ", paste(getwd(), "survival_analysis_plots_top10.pdf", sep = "/"))
        output$survival <- renderText({
          paste("Plots and .csv saved in ", paste(getwd(), "survival_analysis_plots_top10.pdf", sep = "/"))
        })

      }
      dev.off()

    } else {
      cat("Skipping survival analysis.")
    }
    ### HTG_analysis
    ###################################



    # Muestra el resultado en la consola laia
    output$console_output <- renderPrint({
      result
    })
  })
}

shinyApp(ui = ui, server = server)
