#' SpliceImpactR Shiny application
#'
#' A self-contained Shiny interface that exposes the most common parts of the
#' SpliceImpactR workflow. The app bundles a fast demo mode that uses the
#' package's test annotation and toy rMATS/HIT Index outputs, while still
#' allowing users to point the interface at cached or locally downloaded
#' resources for larger analyses. To launch the app from an installed package
#' run:
#'
#' ```r
#' shiny::runApp(system.file("shiny", package = "SpliceImpactR"))
#' ```
#'
library(shiny)
library(data.table)
library(ggplot2)
library(SpliceImpactR)

# Prepares the toy manifest distributed with the package for instant runs.
demo_sample_frame <- function() {
  data.frame(
    path = c(
      check_extdata_dir("rawData/control_S5/"),
      check_extdata_dir("rawData/control_S6/"),
      check_extdata_dir("rawData/control_S7/"),
      check_extdata_dir("rawData/control_S8/"),
      check_extdata_dir("rawData/case_S1/"),
      check_extdata_dir("rawData/case_S2/"),
      check_extdata_dir("rawData/case_S3/"),
      check_extdata_dir("rawData/case_S4/")
    ),
    sample_name = c("S5", "S6", "S7", "S8", "S1", "S2", "S3", "S4"),
    condition = c(rep("control", 4), rep("case", 4)),
    stringsAsFactors = FALSE
  )
}

supported_events <- c("ALE", "AFE", "MXE", "SE", "A3SS", "A5SS", "RI")

ui <- fluidPage(
  titlePanel("SpliceImpactR Shiny"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      actionButton("load_demo_all", "Load bundled demo (annotations + events + proteins + PPI)", width = "100%"),
      hr(),
      h4("Annotations"),
      radioButtons(
        "annotation_mode",
        "Annotation source",
        choices = c(
          "Bundled test data (fast)" = "test",
          "Cached RDS directory" = "cached",
          "Download from GENCODE" = "link",
          "Local files" = "path"
        ),
        selected = "test"
      ),
      selectInput("annotation_species", "Species", c("human", "mouse"), selected = "human"),
      numericInput("annotation_release", "GENCODE release", value = 45, min = 1),
      textInput("annotation_base", "Cache/output directory", value = "annotation_cache"),
      conditionalPanel(
        condition = "input.annotation_mode == 'path'",
        textInput("gtf_path", "GTF path"),
        textInput("tx_path", "Transcript FASTA path"),
        textInput("aa_path", "Protein FASTA path")
      ),
      actionButton("load_annotations", "Load annotations", width = "100%"),
      hr(),
      h4("Splicing data"),
      checkboxInput("use_demo_data", "Use bundled toy rMATS + HIT Index output", value = TRUE),
      fileInput("sample_manifest", "Upload manifest (path, sample_name, condition)", accept = c(".csv", ".tsv")),
      selectInput("event_types", "Event types", choices = supported_events, selected = supported_events, multiple = TRUE),
      selectInput("rmats_use", "rMATS counts", choices = c("JCEC", "JC"), selected = "JCEC"),
      checkboxInput("keep_first_last", "Filter for annotated first/last exons", value = TRUE),
      actionButton("load_splicing", "Load splicing data", width = "100%"),
      hr(),
      h4("Protein features"),
      checkboxInput("use_demo_proteins", "Use bundled InterPro/SignalP demo domains", value = TRUE),
      fileInput("protein_file", "Upload protein features (RDS/CSV)", accept = c(".rds", ".csv", ".tsv")),
      actionButton("load_proteins", "Load protein features", width = "100%"),
      hr(),
      h4("PPI interactions"),
      checkboxInput("use_demo_ppi", "Use bundled PPIDM demo set", value = TRUE),
      actionButton("load_ppi", "Load PPI maps", width = "100%"),
      hr(),
      h4("Result filters"),
      textInput("gene_filter", "Global gene filter", placeholder = "ENSG... or gene symbol"),
      textInput("transcript_filter", "Filter by transcript ID", placeholder = "ENST..."),
      textInput("protein_filter", "Filter by protein ID", placeholder = "ENSP..."),
      actionButton("clear_filters", "Clear filters"),
      hr(),
      h4("Differential inclusion"),
      numericInput("min_reads", "Minimum total reads per site/sample", value = 10, min = 0),
      numericInput("padj_thr", "FDR cutoff", value = 0.05, min = 0, step = 0.01),
      numericInput("dpsi_thr", "|Delta PSI| cutoff", value = 0.1, min = 0, step = 0.01),
      actionButton("run_di", "Run differential inclusion", width = "100%"),
      hr(),
      h4("Annotation mapping"),
      helpText("Match events to the loaded annotation to highlight likely transcripts."),
      actionButton("map_events", "Map events to transcripts", width = "100%")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel(
          "Status",
          h4("Annotation status"),
          tableOutput("annotation_summary"),
          h4("Splicing inputs"),
          tableOutput("sample_preview"),
          h4("Event counts"),
          tableOutput("event_summary")
        ),
        tabPanel(
          "Differential inclusion",
          textInput("di_gene_filter", "Gene filter (DI tab)", value = "", placeholder = "Overrides global gene filter"),
          plotOutput("volcano", height = 400),
          downloadButton("download_di_plot", "Download volcano plot"),
          verbatimTextOutput("di_colnames"),
          h4("Significant events"),
          tableOutput("di_table")
        ),
        tabPanel(
          "Event mapping",
          textInput("map_gene_filter", "Gene filter (mapping tab)", value = "", placeholder = "Overrides global gene filter"),
          h4("Transcript matches for loaded events"),
          tableOutput("mapping_summary"),
          tableOutput("mapping_gene_table"),
          actionButton("show_mapping_table", "Show full mapping table")
        ),
        tabPanel(
          "Protein consequences",
          textInput("prot_gene_filter", "Gene filter (protein tab)", value = "", placeholder = "Overrides global gene filter"),
          plotOutput("protein_plot", height = 350),
          h4("Domain overlaps with inclusion exons"),
          tableOutput("protein_summary"),
          tableOutput("protein_gene_table"),
          actionButton("show_protein_table", "Show full protein consequence table")
        ),
        tabPanel(
          "Transcript pairs",
          helpText("Lookup domain/consequence summaries for two transcripts or proteins."),
          selectizeInput("pair_gene", "Gene", choices = NULL, options = list(create = TRUE)),
          selectizeInput("tx_a", "Transcript A", choices = NULL, options = list(create = TRUE)),
          selectizeInput("tx_b", "Transcript B", choices = NULL, options = list(create = TRUE)),
          textInput("prot_a", "Protein A (optional)", placeholder = "ENSP..."),
          textInput("prot_b", "Protein B (optional)", placeholder = "ENSP..."),
          actionButton("run_pair_query", "Query pair"),
          h4("Pair summary"),
          tableOutput("pair_summary"),
          h4("Feature differences (by database)"),
          tableOutput("pair_diff_summary"),
          h4("Transcript-level feature summary"),
          tableOutput("pair_features")
        ),
        tabPanel(
          "Seq/Domain plots",
          actionButton("run_downstream", "Run sequence/domain plots", width = "100%"),
          h4("Alignment summary"),
          plotOutput("alignment_plot", height = 350),
          downloadButton("download_alignment_plot", "Download alignment plot"),
          h4("Length comparison"),
          plotOutput("length_plot", height = 350),
          downloadButton("download_length_plot", "Download length plot"),
          h4("Domain enrichment"),
          plotOutput("domain_plot", height = 400),
          downloadButton("download_domain_plot", "Download domain plot"),
          tableOutput("domain_table")
        ),
        tabPanel(
          "Integrative summary",
          actionButton("run_integrative", "Run integrative summary", width = "100%"),
          plotOutput("integrative_plot", height = 600),
          downloadButton("download_integrative_plot", "Download integrative plot"),
          h4("Event counts by type"),
          tableOutput("integrative_by_type"),
          h4("PPI/domain prevalence"),
          tableOutput("integrative_domains"),
          h4("Relative use (pre vs post filter)"),
          tableOutput("integrative_relative_use")
        ),
        tabPanel(
          "Downloads",
          helpText("Export full tables for follow-up in your own scripts."),
          downloadButton("download_di", "Download differential inclusion results"),
          downloadButton("download_mapping", "Download event-to-transcript map"),
          downloadButton("download_proteins", "Download protein consequence table")
        ),
        tabPanel(
          "PPI",
          textInput("ppi_gene_filter", "Gene filter (PPI tab)", value = "", placeholder = "Overrides global gene filter"),
          helpText("n_ppi = total unique partner transcripts per event (from PPIDM). partners = sample of partner transcript IDs."),
          h4("PPI impact summary"),
          plotOutput("ppi_plot", height = 350),
          downloadButton("download_ppi_plot", "Download PPI plot"),
          tableOutput("ppi_summary"),
          tableOutput("ppi_gene_table"),
          actionButton("show_ppi_table", "Show full PPI table")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    annotations = NULL,
    sample_frame = NULL,
    splicing = NULL,
    di = NULL,
    mapped = NULL,
    protein_features = NULL,
    exon_features = NULL,
    ppi_pairs = NULL
  )
  
  output$di_colnames <- renderPrint({
    if (is.null(rv$di)) return("DI results not loaded yet.")
    colnames(as.data.table(rv$di))
  })
  apply_filters <- function(dt, gene_col = NULL, tx_col = NULL, prot_col = NULL, gene_override = NULL) {
    out <- as.data.table(dt)
    gene_col <- if (!is.null(gene_col) && gene_col %in% names(out)) gene_col else NULL
    if (is.null(gene_col)) {
      if ("gene_id" %in% names(out)) gene_col <- "gene_id"
      else if ("gene_id.x" %in% names(out)) gene_col <- "gene_id.x"
    }
    active_gene_filter <- if (is.null(gene_override) || !nzchar(gene_override)) input$gene_filter else gene_override
    if (!is.null(gene_col) && nzchar(active_gene_filter)) {
      filt <- grepl(active_gene_filter, out[[gene_col]], ignore.case = TRUE) |
        ("gene_name" %in% names(out) && grepl(active_gene_filter, out$gene_name, ignore.case = TRUE))
      out <- out[filt]
    }
    if (!is.null(tx_col) && nzchar(input$transcript_filter) && tx_col %in% names(out)) {
      out <- out[grepl(input$transcript_filter, out[[tx_col]], ignore.case = TRUE)]
    }
    if (!is.null(prot_col) && nzchar(input$protein_filter) && prot_col %in% names(out)) {
      out <- out[grepl(input$protein_filter, out[[prot_col]], ignore.case = TRUE)]
    }
    out
  }
  
  normalize_di_cols <- function(dt) {
    D <- as.data.table(dt)
    dpsi_cols <- intersect(c("delta_psi", "delta.psi", "delta_psi.x", "delta.psi.x"), names(D))
    if (length(dpsi_cols)) data.table::setnames(D, dpsi_cols[1], "delta_psi")
    padj_cols <- intersect(c("padj", "padj.x", "FDR", "fdr", "FDR.x"), names(D))
    if (length(padj_cols)) data.table::setnames(D, padj_cols[1], "padj")
    missing_cols <- setdiff(c("delta_psi", "padj"), names(D))
    if (length(missing_cols)) {
      warn_msg <- paste("DI table missing columns:", paste(missing_cols, collapse = ", "))
      message(warn_msg)
      showNotification(warn_msg, type = "error")
      return(NULL)
    }
    D
  }
  
  observeEvent(input$load_annotations, {
    withProgress(message = "Loading annotations", value = 0, {
      mode <- input$annotation_mode
      species <- input$annotation_species
      release <- input$annotation_release
      base_dir <- input$annotation_base
      incProgress(0.2, detail = paste("Mode:", mode))
      
      ann <- switch(
        mode,
        test = get_annotation(load = "test"),
        cached = get_annotation(load = "cached", base_dir = base_dir, species = species, release = release),
        link = get_annotation(load = "link", base_dir = base_dir, species = species, release = release),
        path = {
          req(input$gtf_path, input$tx_path, input$aa_path)
          get_annotation(
            load = "path",
            base_dir = base_dir,
            species = species,
            release = release,
            gtf_path = input$gtf_path,
            transcript_path = input$tx_path,
            translation_path = input$aa_path
          )
        }
      )
      
      rv$annotations <- ann
      incProgress(1, detail = "Annotations ready")
    })
  })
  
  observeEvent(input$load_demo_all, {
    withProgress(message = "Loading bundled demo pipeline", value = 0, {
      incProgress(0.1, detail = "Annotations")
      ann <- get_annotation(load = "test")
      rv$annotations <- ann
      
      incProgress(0.2, detail = "Splicing (demo)")
      manifest <- demo_sample_frame()
      data <- get_rmats_hit(
        manifest,
        event_types = input$event_types,
        use = input$rmats_use,
        keep_annotated_first_last = isTRUE(input$keep_first_last)
      )
      rv$splicing <- data
      rv$sample_frame <- manifest
      
      incProgress(0.35, detail = "Proteins (demo InterPro/SignalP)")
      pf <- get_protein_features(
        biomaRt_databases = c("interpro", "signalp"),
        gtf_df = ann$annotations,
        test = TRUE
      )
      pf_dt <- as.data.table(pf)
      if (!"transcript_id" %in% names(pf_dt) && "ensembl_transcript_id" %in% names(pf_dt)) pf_dt[, transcript_id := ensembl_transcript_id]
      if (!"gene_id" %in% names(pf_dt) && "gene_id.x" %in% names(pf_dt)) pf_dt[, gene_id := gene_id.x]
      rv$protein_features <- pf_dt
      ef <- get_exon_features(ann$annotations, pf_dt)
      if (!"transcript_id" %in% names(ef) && "ensembl_transcript_id" %in% names(ef)) ef[, transcript_id := ensembl_transcript_id]
      if ("gene_id.x" %in% names(ef) && !"gene_id" %in% names(ef)) setnames(ef, "gene_id.x", "gene_id")
      rv$exon_features <- ef
      
      incProgress(0.5, detail = "Event mapping")
      rv$mapped <- match_events_to_annotations_vec(rv$splicing, rv$annotations$annotations)
      
      incProgress(0.65, detail = "Differential inclusion (demo)")
      rv$di <- get_differential_inclusion(
        rv$splicing,
        min_total_reads = input$min_reads,
        cooks_cutoff = "4/n",
        parallel_glm = FALSE,
        verbose = FALSE
      )
      
      incProgress(0.8, detail = "PPI (demo PPIDM)")
      ppidm <- get_ppidm(test = TRUE)
      ppi_raw <- get_isoform_interactions(rv$protein_features, ppidm, init = TRUE, save = FALSE)
      ppi_dt <- as.data.table(ppi_raw)
      tx_a <- intersect(c("txA", "transcript_id_A", "transcriptA"), names(ppi_dt))
      tx_b <- intersect(c("txB", "transcript_id_B", "transcriptB"), names(ppi_dt))
      if (length(tx_a) == 0 || length(tx_b) == 0) {
        tx_a <- intersect(c("transcript_id", "ensembl_transcript_id"), names(ppi_dt))[1]
        tx_b <- intersect(c("i.transcript_id", "partner_transcript_id"), names(ppi_dt))[1]
      }
      part1 <- ppi_dt[, .(transcript_id = .SD[[tx_a]], partner_id = .SD[[tx_b]])]
      part2 <- ppi_dt[, .(transcript_id = .SD[[tx_b]], partner_id = .SD[[tx_a]])]
      partners <- data.table::rbindlist(list(part1, part2), use.names = TRUE, fill = TRUE)
      rv$ppi_pairs <- unique(partners[nzchar(transcript_id) & nzchar(partner_id)])
      
      incProgress(1, detail = "Demo ready")
    })
  })
  
  observeEvent(input$load_splicing, {
    manifest <- NULL
    if (isTruthy(input$use_demo_data)) {
      manifest <- demo_sample_frame()
    }
    
    if (!is.null(input$sample_manifest)) {
      ext <- tools::file_ext(input$sample_manifest$name)
      delim <- if (tolower(ext) == "tsv") "\t" else ","
      manifest <- fread(input$sample_manifest$datapath, sep = delim, data.table = FALSE)
    }
    
    validate(need(!is.null(manifest), "Provide a manifest or enable the demo dataset."))
    validate(need(all(c("path", "sample_name", "condition") %in% names(manifest)),
                  "Manifest must include path, sample_name, and condition columns."))
    
    rv$sample_frame <- manifest
    
    withProgress(message = "Loading splicing data", value = 0, {
      data <- get_rmats_hit(
        manifest,
        event_types = input$event_types,
        use = input$rmats_use,
        keep_annotated_first_last = isTRUE(input$keep_first_last)
      )
      rv$splicing <- data
      rv$di <- NULL
      rv$mapped <- NULL
      rv$protein_features <- NULL
      rv$exon_features <- NULL
      incProgress(1, detail = "Events loaded")
    })
  })
  
  observeEvent(input$load_proteins, {
    req(rv$annotations)
    withProgress(message = "Loading protein features", value = 0, {
      pf <- NULL
      if (isTruthy(input$use_demo_proteins)) {
        pf <- get_protein_features(
          biomaRt_databases = c("interpro", "signalp"),
          gtf_df = rv$annotations$annotations,
          test = TRUE
        )
      }
      
      if (!is.null(input$protein_file)) {
        ext <- tolower(tools::file_ext(input$protein_file$name))
        pf <- if (ext == "rds") {
          readRDS(input$protein_file$datapath)
        } else {
          delim <- if (ext == "tsv") "\t" else ","
          fread(input$protein_file$datapath, sep = delim, data.table = TRUE)
        }
      }
      
      validate(need(!is.null(pf), "Provide a protein feature file or enable demo domains."))
      
      pf_dt <- as.data.table(pf)
      if (!"transcript_id" %in% names(pf_dt) && "ensembl_transcript_id" %in% names(pf_dt)) {
        pf_dt[, transcript_id := ensembl_transcript_id]
      }
      if (!"gene_id" %in% names(pf_dt) && "gene_id.x" %in% names(pf_dt)) {
        pf_dt[, gene_id := gene_id.x]
      }
      
      rv$protein_features <- pf_dt
      ef <- get_exon_features(rv$annotations$annotations, pf_dt)
      if (!"transcript_id" %in% names(ef) && "ensembl_transcript_id" %in% names(ef)) {
        ef[, transcript_id := ensembl_transcript_id]
      }
      if ("gene_id.x" %in% names(ef) && !"gene_id" %in% names(ef)) {
        setnames(ef, "gene_id.x", "gene_id")
      }
      rv$exon_features <- ef
      incProgress(1, detail = "Protein features ready")
    })
  })
  
  observeEvent(input$load_ppi, {
    req(rv$protein_features)
    withProgress(message = "Loading PPI interactions", value = 0, {
      ppidm <- if (isTRUE(input$use_demo_ppi)) {
        get_ppidm(test = TRUE)
      } else {
        get_ppidm(download = TRUE)
      }
      ppi_raw <- get_isoform_interactions(rv$protein_features, ppidm, init = TRUE, save = FALSE)
      ppi_dt <- as.data.table(ppi_raw)
      
      tx_a <- intersect(c("txA", "transcript_id_A", "transcriptA"), names(ppi_dt))
      tx_b <- intersect(c("txB", "transcript_id_B", "transcriptB"), names(ppi_dt))
      if (length(tx_a) == 0 || length(tx_b) == 0) {
        tx_a <- intersect(c("transcript_id", "ensembl_transcript_id"), names(ppi_dt))[1]
        tx_b <- intersect(c("i.transcript_id", "partner_transcript_id"), names(ppi_dt))[1]
      }
      validate(need(!is.null(tx_a) && !is.null(tx_b), "Could not detect transcript columns in PPI map"))
      
      part1 <- ppi_dt[, .(transcript_id = .SD[[tx_a]], partner_id = .SD[[tx_b]])]
      part2 <- ppi_dt[, .(transcript_id = .SD[[tx_b]], partner_id = .SD[[tx_a]])]
      partners <- data.table::rbindlist(list(part1, part2), use.names = TRUE, fill = TRUE)
      partners <- partners[nzchar(transcript_id) & nzchar(partner_id)]
      rv$ppi_pairs <- unique(partners)
      incProgress(1, detail = "PPI ready")
    })
  })
  
  observeEvent(input$clear_filters, {
    updateTextInput(session, "gene_filter", value = "")
    updateTextInput(session, "transcript_filter", value = "")
    updateTextInput(session, "protein_filter", value = "")
  })
  
  observeEvent(input$run_di, {
    req(rv$splicing)
    withProgress(message = "Running quasi-binomial models", value = 0, {
      di <- get_differential_inclusion(
        rv$splicing,
        min_total_reads = input$min_reads,
        cooks_cutoff = "4/n",
        parallel_glm = FALSE,
        verbose = FALSE
      )
      rv$di <- di
      incProgress(1, detail = "Differential inclusion complete")
    })
  })
  
  observeEvent(input$map_events, {
    req(rv$splicing, rv$annotations)
    withProgress(message = "Matching events to transcripts", value = 0, {
      mapped <- match_events_to_annotations_vec(rv$splicing, rv$annotations$annotations)
      rv$mapped <- mapped
      incProgress(1, detail = "Mapping complete")
    })
  })
  
  output$annotation_summary <- renderTable({
    req(rv$annotations)
    ann <- as.data.table(rv$annotations$annotations)
    data.frame(
      Genes = uniqueN(na.omit(ann$gene_id)),
      Transcripts = ann[type == "transcript", uniqueN(transcript_id)],
      Exons = ann[type == "exon", .N],
      Protein.coding.transcripts = ann[!is.na(protein_id), uniqueN(transcript_id)]
    )
  })
  
  output$sample_preview <- renderTable({
    if (is.null(rv$sample_frame)) return(NULL)
    head(rv$sample_frame, 10)
  })
  
  output$event_summary <- renderTable({
    if (is.null(rv$splicing)) return(NULL)
    rv$splicing[, .(events = .N, samples = uniqueN(sample), genes = uniqueN(gene_id)), by = event_type]
  })
  
  di_plot_obj <- reactive({
    req(rv$di)
    di_dt <- normalize_di_cols(rv$di)
    if (is.null(di_dt)) {
      showNotification("Differential inclusion results missing padj or delta_psi.", type = "error")
      return(NULL)
    }
    plot_di_volcano_dt(di_dt, padj_thr = input$padj_thr, dpsi_thr = input$dpsi_thr)
  })
  
  output$volcano <- renderPlot({
    di_plot_obj()
  })
  
  output$di_table <- renderTable({
    req(rv$di)
    sig <- normalize_di_cols(rv$di)
    if (is.null(sig)) {
      showNotification("Differential inclusion results missing padj or delta_psi.", type = "error")
      return(NULL)
    }
    sig <- sig[is.finite(padj) & padj <= input$padj_thr &
                 is.finite(delta_psi) & abs(delta_psi) >= input$dpsi_thr]
    sig <- apply_filters(sig, gene_col = "gene_id", gene_override = input$di_gene_filter)
    sig[order(padj)][1:min(.N, 25)]
  })
  
  output$download_di_plot <- downloadHandler(
    filename = function() paste0("spliceimpactr_di_volcano_", Sys.Date(), ".png"),
    content = function(file) {
      plt <- di_plot_obj()
      req(!is.null(plt))
      ggsave(file, plot = plt, width = 7, height = 5, dpi = 300)
    }
  )
  
  output$mapping_summary <- renderTable({
    req(rv$mapped)
    keep <- c("event_type", "gene_id", "transcript_id")
    tmp <- rv$mapped[, ..keep]
    tmp <- apply_filters(tmp, gene_col = "gene_id", tx_col = "transcript_id", gene_override = input$map_gene_filter)
    tmp[, .(n_events = .N, n_transcripts = uniqueN(transcript_id)), by = event_type][order(-n_events)]
  })
  
  output$mapping_gene_table <- renderTable({
    req(rv$mapped)
    gene_filter_active <- nzchar(input$map_gene_filter) || nzchar(input$gene_filter)
    if (!gene_filter_active) return(NULL)
    keep <- c("event_id", "event_type", "gene_id", "transcript_id", "exons", "inc_exons_by_idx")
    dat <- apply_filters(rv$mapped[, ..keep], gene_col = "gene_id", tx_col = "transcript_id", gene_override = input$map_gene_filter)
    if (!nrow(dat)) return(NULL)
    dat[order(event_id)][1:min(.N, 50)]
  })
  
  observeEvent(input$show_mapping_table, {
    req(rv$mapped)
    keep <- c("event_id", "event_type", "gene_id", "transcript_id", "exons", "inc_exons_by_idx")
    dat <- apply_filters(rv$mapped[, ..keep], gene_col = "gene_id", tx_col = "transcript_id", gene_override = input$map_gene_filter)
    showModal(modalDialog(
      title = "Mapping table",
      size = "l",
      easyClose = TRUE,
      div(style = "max-height: 500px; overflow-y: auto;",
          dataTableOutput("mapping_full"))
    ))
    output$mapping_full <- renderDataTable(dat, options = list(pageLength = 25))
  })
  
  protein_consequences <- reactive({
    req(rv$mapped, rv$exon_features)
    m <- as.data.table(rv$mapped)
    # expand inclusion exon IDs
    inc_long <- m[, .(exon_id = unlist(strsplit(inc_exons_by_idx, ";", fixed = TRUE))),
                  by = .(event_id, event_type, gene_id, transcript_id)]
    inc_long <- inc_long[nzchar(exon_id)]
    if (!nrow(inc_long)) return(data.table())
    cons <- merge(
      inc_long,
      rv$exon_features,
      by = "exon_id",
      allow.cartesian = TRUE
    )
    if ("gene_id.x" %in% names(cons) && !"gene_id" %in% names(cons)) {
      cons[, gene_id := gene_id.x]
    }
    if (!"transcript_id" %in% names(cons) && "ensembl_transcript_id" %in% names(cons)) {
      cons[, transcript_id := ensembl_transcript_id]
    }
    cons <- unique(cons, by = c("event_id", "exon_id", "feature_id", "database", "transcript_id"))
    cons[]
  })
  
  output$protein_plot <- renderPlot({
    cons <- protein_consequences()
    if (!nrow(cons)) return(NULL)
    cons <- apply_filters(cons, gene_col = "gene_id", tx_col = "transcript_id", prot_col = "ensembl_peptide_id", gene_override = input$prot_gene_filter)
    agg <- cons[, .N, by = .(event_type, database)][order(event_type)]
    ggplot(agg, aes(x = event_type, y = N, fill = database)) +
      geom_col(position = "stack") +
      labs(x = "Event type", y = "Overlapping domain count", fill = "Database") +
      theme_classic()
  })
  
  output$protein_summary <- renderTable({
    cons <- protein_consequences()
    if (!nrow(cons)) return(NULL)
    cons <- apply_filters(cons, gene_col = "gene_id", tx_col = "transcript_id", prot_col = "ensembl_peptide_id", gene_override = input$prot_gene_filter)
    cons[, .(
      events = uniqueN(event_id),
      transcripts = uniqueN(transcript_id),
      domains = uniqueN(feature_id)
    ), by = database][order(-domains)]
  })
  
  output$protein_gene_table <- renderTable({
    cons <- protein_consequences()
    gene_filter_active <- nzchar(input$prot_gene_filter) || nzchar(input$gene_filter)
    if (!gene_filter_active || !nrow(cons)) return(NULL)
    cons <- apply_filters(cons, gene_col = "gene_id", tx_col = "transcript_id", prot_col = "ensembl_peptide_id", gene_override = input$prot_gene_filter)
    if (!nrow(cons)) return(NULL)
    cols <- c("event_id", "event_type", "gene_id", "transcript_id", "database", "name", "alt_name")
    cons[, ..cols][order(event_id)][1:min(.N, 50)]
  })
  
  observeEvent(input$show_protein_table, {
    cons <- protein_consequences()
    req(nrow(cons))
    cons <- apply_filters(cons, gene_col = "gene_id", tx_col = "transcript_id", prot_col = "ensembl_peptide_id", gene_override = input$prot_gene_filter)
    cols <- c("event_id", "event_type", "gene_id", "transcript_id",
              "exon_id", "database", "name", "alt_name")
    dat <- cons[, ..cols][order(event_id)]
    showModal(modalDialog(
      title = "Protein consequences",
      size = "l",
      easyClose = TRUE,
      div(style = "max-height: 500px; overflow-y: auto;",
          dataTableOutput("protein_full"))
    ))
    output$protein_full <- renderDataTable(dat, options = list(pageLength = 25))
  })
  
  output$download_di <- downloadHandler(
    filename = function() paste0("spliceimpactr_di_results_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$di)
      fwrite(rv$di, file)
    }
  )
  
  output$download_mapping <- downloadHandler(
    filename = function() paste0("spliceimpactr_event_mapping_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$mapped)
      fwrite(rv$mapped, file)
    }
  )
  
  output$download_proteins <- downloadHandler(
    filename = function() paste0("spliceimpactr_protein_consequences_", Sys.Date(), ".csv"),
    content = function(file) {
      cons <- protein_consequences()
      req(cons)
      fwrite(cons, file)
    }
  )
  
  ppi_for_events <- reactive({
    req(rv$mapped, rv$ppi_pairs)
    m <- as.data.table(rv$mapped)
    m <- apply_filters(m, gene_col = "gene_id", tx_col = "transcript_id", gene_override = input$ppi_gene_filter)
    if (!nrow(m)) return(data.table())
    partners <- unique(as.data.table(rv$ppi_pairs))
    partner_list <- partners[, .(
      n_ppi = uniqueN(partner_id),
      partner_sample = paste(head(unique(partner_id), 10), collapse = ";")
    ), by = transcript_id]
    
    per_tx <- merge(m, partner_list, by = "transcript_id", all.x = TRUE)
    per_tx[is.na(n_ppi), `:=`(n_ppi = 0L, partner_sample = "")]
    
    per_ev <- per_tx[, .(
      transcripts = uniqueN(transcript_id),
      n_ppi_total = sum(n_ppi, na.rm = TRUE),
      impacted = any(n_ppi > 0),
      partners = {
        pvec <- unique(unlist(strsplit(partner_sample, ";", fixed = TRUE)))
        pvec <- pvec[nzchar(pvec)]
        if (length(pvec)) paste(head(pvec, 20), collapse = ";") else "none"
      }
    ), by = .(event_id, event_type, gene_id)]
    
    per_ev[, `:=`(
      n_inc_ppi = as.integer(n_ppi_total),
      n_exc_ppi = 0L,
      n_ppi = as.integer(n_ppi_total)
    )]
    
    list(per_tx = per_tx, per_event = per_ev)
  })
  
  ppi_plot_obj <- reactive({
    res <- ppi_for_events()
    if (is.null(res) || !nrow(res$per_event)) return(NULL)
    plot_ppi_summary(res$per_event)
  })
  
  output$ppi_plot <- renderPlot({
    ppi_plot_obj()
  })
  
  output$ppi_summary <- renderTable({
    res <- ppi_for_events()
    if (is.null(res) || !nrow(res$per_event)) return(NULL)
    dt <- res$per_event
    dt[, .(
      events = .N,
      impacted = sum(impacted),
      mean_ppi = round(mean(n_ppi, na.rm = TRUE), 2)
    ), by = event_type][order(-events)]
  })
  
  output$ppi_gene_table <- renderTable({
    res <- ppi_for_events()
    if (is.null(res) || !nrow(res$per_event)) return(NULL)
    dt <- res$per_event
    gene_filter_active <- nzchar(input$ppi_gene_filter) || nzchar(input$gene_filter)
    if (!gene_filter_active) return(NULL)
    dt[order(event_id), .(event_id, event_type, gene_id, n_ppi, partners)][1:min(.N, 50)]
  })
  
  observeEvent(input$show_ppi_table, {
    res <- ppi_for_events()
    req(!is.null(res), nrow(res$per_event))
    dat <- res$per_event[, .(event_id, event_type, gene_id, n_ppi, partners)]
    showModal(modalDialog(
      title = "PPI per event",
      size = "l",
      easyClose = TRUE,
      div(style = "max-height: 500px; overflow-y: auto;",
          dataTableOutput("ppi_full"))
    ))
    output$ppi_full <- renderDataTable(dat, options = list(pageLength = 25))
  })
  
  output$download_ppi_plot <- downloadHandler(
    filename = function() paste0("spliceimpactr_ppi_plot_", Sys.Date(), ".png"),
    content = function(file) {
      plt <- ppi_plot_obj()
      req(!is.null(plt))
      ggsave(file, plot = plt, width = 9, height = 5, dpi = 300)
    }
  )
  
  downstream_results <- eventReactive(input$run_downstream, {
    req(rv$mapped, rv$annotations, rv$exon_features)
    if (is.null(rv$annotations$sequences)) {
      showNotification("Annotation sequences are required for downstream plots.", type = "error")
      return(NULL)
    }
    withProgress(message = "Running sequence/domain plots", value = 0, {
      map_dt <- unique(as.data.table(rv$mapped), by = c("event_id", "transcript_id"))
      # attach delta_psi if present elsewhere
      if (!"delta_psi" %in% names(map_dt)) {
        di_dt <- normalize_di_cols(rv$di)
        if (!is.null(di_dt) && "event_id" %in% names(di_dt)) {
          map_dt <- merge(map_dt, di_dt[, .(event_id, delta_psi)], by = "event_id", all.x = TRUE)
        }
      }
      if (!"delta_psi" %in% names(map_dt) || all(is.na(map_dt$delta_psi))) {
        showNotification("Downstream plots require delta_psi on mapped events; run DI first or supply a table with delta_psi.", type = "error")
        return(NULL)
      }
      
      incProgress(0.2, detail = "Attaching sequences")
      x_seq <- tryCatch(
        attach_sequences(map_dt, rv$annotations$sequences),
        error = function(e) {
          showNotification(paste("attach_sequences failed:", e$message), type = "error")
          return(NULL)
        }
      )
      if (is.null(x_seq)) return(NULL)
      pairs <- tryCatch(get_pairs(x_seq, source = "multi"), error = function(e) {
        showNotification(paste("get_pairs failed:", e$message), type = "error"); NULL
      })
      if (is.null(pairs)) return(NULL)
      
      incProgress(0.4, detail = "Comparing sequences")
      seq_compare <- tryCatch(compare_sequence_frame(pairs, rv$annotations$annotations),
                              error = function(e) {showNotification(paste("compare_sequence_frame failed:", e$message), type="error"); NULL})
      if (is.null(seq_compare)) return(NULL)
      alignment_plot <- tryCatch(plot_alignment_summary(seq_compare), error = function(e) NULL)
      length_plot <- tryCatch(plot_length_comparison(seq_compare), error = function(e) NULL)
      
      domain_plot <- NULL
      enriched_domains <- NULL
      incProgress(0.6, detail = "Domain overlaps")
      hits_domain <- tryCatch(get_domains(seq_compare, rv$exon_features),
                              error = function(e) {showNotification(paste("get_domains failed:", e$message), type="error"); NULL})
      if (!is.null(hits_domain) && !is.null(rv$sample_frame) && !is.null(rv$protein_features)) {
        bg <- tryCatch(get_background(
          source = "hit_index",
          input = rv$sample_frame,
          annotations = rv$annotations$annotations,
          protein_features = rv$protein_features
        ), error = function(e) {showNotification(paste("get_background failed:", e$message), type="error"); NULL})
        if (!is.null(bg)) {
          enriched_domains <- tryCatch(enrich_domains_hypergeo(hits_domain, bg, db_filter = "interpro"),
                                       error = function(e) NULL)
          if (!is.null(enriched_domains) && nrow(enriched_domains)) {
            domain_plot <- tryCatch(plot_enriched_domains_counts(enriched_domains, top_n = 20),
                                    error = function(e) NULL)
          }
        }
      }
      incProgress(1, detail = "Done")
      
      list(
        seq_compare = seq_compare,
        alignment_plot = alignment_plot,
        length_plot = length_plot,
        enriched_domains = enriched_domains,
        domain_plot = domain_plot,
        hits_domain = hits_domain
      )
    })
  })
  
  output$alignment_plot <- renderPlot({
    res <- downstream_results()
    if (is.null(res) || is.null(res$alignment_plot)) return(NULL)
    res$alignment_plot
  })
  
  output$length_plot <- renderPlot({
    res <- downstream_results()
    if (is.null(res) || is.null(res$length_plot)) return(NULL)
    res$length_plot
  })
  
  output$domain_plot <- renderPlot({
    res <- downstream_results()
    if (is.null(res) || is.null(res$domain_plot)) return(NULL)
    res$domain_plot
  })
  
  output$domain_table <- renderTable({
    res <- downstream_results()
    if (is.null(res) || is.null(res$enriched_domains) || inherits(res$enriched_domains, "try-error")) return(NULL)
    head(res$enriched_domains, 20)
  })
  
  output$download_alignment_plot <- downloadHandler(
    filename = function() paste0("spliceimpactr_alignment_plot_", Sys.Date(), ".png"),
    content = function(file) {
      res <- downstream_results()
      if (is.null(res) || is.null(res$alignment_plot)) stop("No alignment plot available")
      ggsave(file, plot = res$alignment_plot, width = 7, height = 5, dpi = 300)
    }
  )
  
  output$download_length_plot <- downloadHandler(
    filename = function() paste0("spliceimpactr_length_plot_", Sys.Date(), ".png"),
    content = function(file) {
      res <- downstream_results()
      if (is.null(res) || is.null(res$length_plot)) stop("No length plot available")
      ggsave(file, plot = res$length_plot, width = 7, height = 5, dpi = 300)
    }
  )
  
  output$download_domain_plot <- downloadHandler(
    filename = function() paste0("spliceimpactr_domain_plot_", Sys.Date(), ".png"),
    content = function(file) {
      res <- downstream_results()
      if (is.null(res) || is.null(res$domain_plot)) stop("No domain plot available")
      ggsave(file, plot = res$domain_plot, width = 7, height = 5, dpi = 300)
    }
  )
  
  integrative_results <- eventReactive(input$run_integrative, {
    dr <- downstream_results()
    if (is.null(dr) || is.null(dr$hits_domain)) {
      showNotification("Run sequence/domain plots first to generate hits_domain.", type = "error")
      return(NULL)
    }
    if (is.null(rv$ppi_pairs)) {
      showNotification("Load PPI maps before running integrative summary.", type = "error")
      return(NULL)
    }
    di_dt <- normalize_di_cols(rv$di)
    if (is.null(di_dt)) return(NULL)
    
    ppi_dt <- as.data.table(rv$ppi_pairs)
    if (!("ensembl_transcript_id" %in% names(ppi_dt))) {
      if ("transcript_id" %in% names(ppi_dt)) setnames(ppi_dt, "transcript_id", "ensembl_transcript_id")
    }
    partner_col <- NULL
    for (cand in c("i.ensembl_transcript_id", "partner_id", "partner_transcript_id")) {
      if (cand %in% names(ppi_dt)) { partner_col <- cand; break }
    }
    if (is.null(partner_col)) {
      showNotification("Could not detect partner transcript column in PPI map.", type = "error")
      return(NULL)
    }
    setnames(ppi_dt, partner_col, "i.ensembl_transcript_id")
    
    hits_final <- tryCatch(get_ppi_switches(dr$hits_domain, ppi_dt),
                           error = function(e) {
                             showNotification(paste("get_ppi_switches failed:", e$message), type = "error")
                             NULL
                           })
    if (is.null(hits_final)) return(NULL)
    
    int <- tryCatch(integrated_event_summary(hits_final, di_dt),
                    error = function(e) {
                      showNotification(paste("integrated_event_summary failed:", e$message), type = "error")
                      NULL
                    })
    int
  })
  
  output$integrative_plot <- renderPlot({
    res <- integrative_results()
    if (is.null(res) || is.null(res$plot)) return(NULL)
    res$plot
  })
  
  output$download_integrative_plot <- downloadHandler(
    filename = function() paste0("spliceimpactr_integrative_plot_", Sys.Date(), ".png"),
    content = function(file) {
      res <- integrative_results()
      if (is.null(res) || is.null(res$plot)) stop("No integrative plot available")
      ggsave(file, plot = res$plot, width = 10, height = 10, dpi = 300)
    }
  )
  
  output$integrative_by_type <- renderTable({
    res <- integrative_results()
    if (is.null(res) || is.null(res$summaries$by_type)) return(NULL)
    res$summaries$by_type
  })
  
  output$integrative_domains <- renderTable({
    res <- integrative_results()
    if (is.null(res) || is.null(res$summaries$domain_prevalence)) return(NULL)
    res$summaries$domain_prevalence
  })
  
  output$integrative_relative_use <- renderTable({
    res <- integrative_results()
    if (is.null(res) || is.null(res$summaries$relative_use)) return(NULL)
    res$summaries$relative_use
  })
  
  observeEvent(rv$annotations, {
    ann <- as.data.table(rv$annotations$annotations)
    gene_choices <- unique(ann[type == "gene", gene_id])
    updateSelectizeInput(session, "pair_gene", choices = gene_choices, server = TRUE)
    updateSelectizeInput(session, "tx_a", choices = character(0), server = TRUE)
    updateSelectizeInput(session, "tx_b", choices = character(0), server = TRUE)
  }, ignoreNULL = TRUE)
  
  observeEvent(input$pair_gene, {
    req(rv$annotations)
    ann <- as.data.table(rv$annotations$annotations)
    tx_choices <- ann[type == "transcript" & gene_id == input$pair_gene, unique(transcript_id)]
    updateSelectizeInput(session, "tx_a", choices = tx_choices, server = TRUE, selected = NULL)
    updateSelectizeInput(session, "tx_b", choices = tx_choices, server = TRUE, selected = NULL)
  }, ignoreInit = TRUE)
  
  pair_query <- eventReactive(input$run_pair_query, {
    if (is.null(rv$exon_features) || !nrow(rv$exon_features)) {
      showNotification("Load protein features before querying transcript pairs.", type = "warning")
      return(data.table())
    }
    if (!nzchar(input$pair_gene)) {
      showNotification("Select a gene before querying transcript pairs.", type = "warning")
      return(data.table())
    }
    tx_ids <- c(input$tx_a, input$tx_b)
    tx_ids <- tx_ids[nzchar(tx_ids)]
    if (length(tx_ids) != 2) {
      showNotification("Pick exactly two transcripts to compare.", type = "warning")
      return(data.table())
    }
    ef <- as.data.table(rv$exon_features)
    ef <- ef[gene_id == input$pair_gene]
    if (!nrow(ef)) {
      showNotification("No features found for this gene in loaded annotations.", type = "message")
      return(data.table())
    }
    prot_ids <- c(input$prot_a, input$prot_b)
    prot_ids <- prot_ids[nzchar(prot_ids)]
    if (length(prot_ids)) {
      ef <- ef[ensembl_peptide_id %in% prot_ids | transcript_id %in% tx_ids]
    } else {
      ef <- ef[transcript_id %in% tx_ids]
    }
    if (!nrow(ef)) {
      showNotification("No overlapping features for the selected transcripts/proteins.", type = "message")
      return(data.table())
    }
    unique(ef)
  })
  
  output$pair_summary <- renderTable({
    ef <- pair_query()
    if (!nrow(ef)) return(NULL)
    ef[, .(
      transcripts = uniqueN(transcript_id),
      proteins = uniqueN(ensembl_peptide_id),
      domains = uniqueN(feature_id)
    ), by = .(database)][order(-domains)]
  })
  
  output$pair_features <- renderTable({
    ef <- pair_query()
    if (!nrow(ef)) return(NULL)
    cols <- c("transcript_id", "ensembl_peptide_id", "database", "name", "alt_name")
    unique(ef[, ..cols])[order(transcript_id, database, name)][1:min(.N, 50)]
  })
  
  pair_diff_tables <- reactive({
    ef <- pair_query()
    if (!nrow(ef) || !nzchar(input$tx_a) || !nzchar(input$tx_b)) return(NULL)
    cols <- c("transcript_id", "database", "feature_id", "name", "alt_name")
    ft <- unique(ef[, ..cols])
    key_cols <- c("database", "feature_id", "name", "alt_name")
    a_feats <- ft[transcript_id == input$tx_a]
    b_feats <- ft[transcript_id == input$tx_b]
    
    shared <- merge(a_feats[, ..key_cols], b_feats[, ..key_cols], by = key_cols)
    if (nrow(shared)) shared[, status := "shared"]
    
    only_a <- fsetdiff(a_feats[, ..key_cols], b_feats[, ..key_cols])
    if (nrow(only_a)) only_a[, status := "only A"]
    
    only_b <- fsetdiff(b_feats[, ..key_cols], a_feats[, ..key_cols])
    if (nrow(only_b)) only_b[, status := "only B"]
    
    out <- rbindlist(list(shared, only_a, only_b), use.names = TRUE, fill = TRUE)
    if (!nrow(out)) return(NULL)
    out[, status := factor(status, levels = c("shared", "only A", "only B"))]
    summary_tbl <- out[, .(
      shared = sum(status == "shared"),
      only_A = sum(status == "only A"),
      only_B = sum(status == "only B"),
      total_A = sum(status %in% c("shared", "only A")),
      total_B = sum(status %in% c("shared", "only B"))
    ), by = database][order(-shared - only_A - only_B)]
    
    detail_diff <- unique(out[status != "shared"])
    if (nrow(detail_diff)) {
      detail_diff[, transcript_id := ifelse(status == "only A", input$tx_a, input$tx_b)]
    }
    
    list(
      summary = summary_tbl,
      differences = detail_diff
    )
  })
  
  output$pair_diff_summary <- renderTable({
    tbls <- pair_diff_tables()
    if (is.null(tbls)) return(NULL)
    tbls$summary
  })
  
  output$pair_diff_full <- renderTable({
    tbls <- pair_diff_tables()
    if (is.null(tbls)) return(NULL)
    diffs <- tbls$differences
    if (!nrow(diffs)) return(NULL)
    diffs[, .(status, transcript_id, database, name, alt_name)][
      order(status, database, name)
    ][1:min(.N, 100)]
  })
}

shinyApp(ui, server)