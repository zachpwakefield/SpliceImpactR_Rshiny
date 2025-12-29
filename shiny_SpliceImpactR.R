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
      actionButton("map_events", "Map events to transcripts", width = "100%"),
      hr(),
      h4("Sequence/domain + PPI"),
      helpText("Run after mapping and loading proteins/PPIs to populate Protein consequences, PPI, and Integrative tabs."),
      actionButton("run_downstream", "Run sequence/domain + PPI summary", width = "100%"),
      hr(),
      h4("HIT/PSI comparisons"),
      helpText("Compute HIT index and PSI comparisons across conditions."),
      actionButton("run_hit_overview_sidebar", "Run HIT/PSI summaries", width = "100%")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel(
          "How to use",
          h3("What SpliceImpactR does"),
          p("SpliceImpactR maps alternative splicing events to coding consequences, domain changes, and isoform-specific protein-protein interaction (PPI) switches. The Shiny app mirrors the published workflow described in the bioRxiv preprint and the package README, so you can reproduce figures directly from curated rMATS/HIT index outputs."),
          tags$ul(
            tags$li("Quantifies differential inclusion (|\u0394PSI| and FDR)."),
            tags$li("Aligns inclusion vs. exclusion isoforms and flags frame changes or rescues."),
            tags$li("Calls domain gains/losses with InterPro/SignalP overlap via ", code("get_domains(seq_compare, exon_features)"), "."),
            tags$li("Calls PPI rewiring with ", code("hits_final <- get_ppi_switches(hits_domain, ppi)"), " using PPIDM-derived isoform pairs.")
          ),
          h4("Data you need"),
          tags$ul(
            tags$li("rMATS junction counts (JCEC/JC) plus HIT index output for your samples (or load the bundled toy set)."),
            tags$li("GENCODE annotations (bundled test, cached RDS, GENCODE download, or local GTF/FASTA)."),
            tags$li("Protein features (demo InterPro/SignalP or your own overlaps)."),
            tags$li("PPI map (demo PPIDM predictions or downloaded PPIDM).")
          ),
          h4("Step-by-step workflow in this app"),
          tags$ol(
            tags$li("Load annotations for your species/release."),
            tags$li("Load splicing inputs (toggle demo or upload manifest)."),
            tags$li("Load protein features, then load PPIs (needed for domain + PPI tabs)."),
            tags$li("Run differential inclusion, then map events to transcripts."),
            tags$li("Click the sidebar button ", strong("Run sequence/domain + PPI summary"), " to compute sequence alignments, domain hits, and PPI switches; all downstream tabs will populate after this."),
            tags$li("Explore Protein consequences, PPI, Seq/Domain plots, and Integrative summary tabs."),
            tags$li("Use Downloads to export DI, mapping, and consequence tables.")
          ),
          h4("Per-tab guidance"),
          tags$ul(
            tags$li(strong("Differential inclusion:"), " tune FDR and |\u0394PSI| thresholds; the table reflects active filters."),
            tags$li(strong("Event mapping:"), " links each event to transcripts before any domain/PPI calls."),
            tags$li(strong("Protein consequences:"), " populated after sequence/domain run; shows inclusion vs. exclusion domain gains/losses by database."),
            tags$li(strong("PPI:"), " populated after sequence/domain run with PPIs loaded; summarizes isoform-specific partners from ", code("get_ppi_switches"), "."),
            tags$li(strong("Seq/Domain plots:"), " alignment, length, and enrichment summaries reused by other tabs."),
            tags$li(strong("Integrative summary:"), " aggregates event counts, domain prevalence, and PPI rewiring across types."),
            tags$li(strong("Transcript pairs:"), " quick lookup of domain differences for two isoforms of the same gene."),
            tags$li(strong("Filters:"), " global gene/transcript/protein filters apply across tabs; per-tab gene filters override the global setting.")
          ),
          h4("Learn more"),
          tags$ul(
            tags$li(tags$a(href = "https://www.biorxiv.org/content/10.1101/2025.06.20.660706v1", target = "_blank", "bioRxiv preprint")),
            tags$li(tags$a(href = "https://github.com/fiszbein-lab/SpliceImpactR/tree/main", target = "_blank", "SpliceImpactR package repository")),
            tags$li(tags$a(href = "https://github.com/fiszbein-lab/SpliceImpactR/blob/main/README.md", target = "_blank", "README workflow and examples")),
            tags$li(tags$a(href = "https://www.fiszbeinlab.com/", target = "_blank", "Fiszbein Lab")),
            tags$li(tags$a(href = "https://github.com/Xinglab/rmats-turbo/tree/v4.2.0", target = "_blank", "rMATS"), " and ",
                    tags$a(href = "https://github.com/thepailab/HITindex", target = "_blank", "HITindex"), " references for input generation")
          )
        ),
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
          helpText("Run sequence/domain plots after loading proteins to summarize inclusion vs. exclusion domain changes."),
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
          helpText("Use the sidebar button to run sequence/domain + PPI summary, then view plots below."),
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
          "PPI",
          textInput("ppi_gene_filter", "Gene filter (PPI tab)", value = "", placeholder = "Overrides global gene filter"),
          helpText("Run sequence/domain plots with PPIs loaded to populate this tab."),
          helpText("n_ppi = total unique partner transcripts per event (from PPIDM). partners = sample of partner transcript IDs."),
          h4("PPI impact summary"),
          plotOutput("ppi_plot", height = 350),
          downloadButton("download_ppi_plot", "Download PPI plot"),
          tableOutput("ppi_summary"),
          tableOutput("ppi_gene_table"),
          actionButton("show_ppi_table", "Show full PPI table")
        ),
        tabPanel(
          "Event probe",
          helpText("Select an event ID from DI results to probe PSI by sample/condition."),
          selectizeInput("probe_event", "Event ID", choices = NULL, options = list(placeholder = "Select event_id")),
          actionButton("run_probe", "Probe event", width = "100%"),
          h4("PSI by sample/condition"),
          plotOutput("probe_plot", height = 350)
        ),
        tabPanel(
          "HIT/PSI summaries",
          helpText("Run HIT index and PSI overview comparisons across conditions."),
          actionButton("run_hit_overview", "Run HIT/PSI summaries", width = "100%"),
          h4("HIT index comparison"),
          plotOutput("hit_compare_plot", height = 400),
          h4("PSI overview (AFE example)"),
          plotOutput("psi_overview_plot", height = 400),
          h4("Proximal vs distal (AFE/ALE)"),
          plotOutput("proximal_plot", height = 300)
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
    di_norm = NULL,
    di_sig = NULL,
    matched = NULL,
    protein_features = NULL,
    exon_features = NULL,
    ppidm = NULL
  )
  downstream_results <- reactiveVal(NULL)
  
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
      f1 <- grepl(active_gene_filter, out[[gene_col]], ignore.case = TRUE)
      f1[is.na(f1)] <- FALSE
      if ("gene_name" %in% names(out)) {
        f2 <- grepl(active_gene_filter, out$gene_name, ignore.case = TRUE)
        f2[is.na(f2)] <- FALSE
        f1 <- f1 | f2
      }
      out <- out[f1]
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
  
  first_available <- function(dt, cols) {
    available <- cols[cols %in% names(dt)]
    if (!length(available)) return(rep(NA_character_, nrow(dt)))
    for (c in available) {
      candidate <- dt[[c]]
      if (!all(is.na(candidate))) return(candidate)
    }
    rep(NA_character_, nrow(dt))
  }
  
  collapse_list_vals <- function(x, n = 5) {
    vals <- as.character(unlist(x))
    vals <- vals[nzchar(vals) & !is.na(vals)]
    if (!length(vals)) return("")
    paste(head(vals, n), collapse = ";")
  }
  has_df_rows <- function(x) {
    is.data.frame(x) && nrow(x) > 0
  }
  as_dt_or_null <- function(x) {
    if (is.null(x) || !is.data.frame(x)) return(NULL)
    as.data.table(x)
  }
  
  build_domain_long <- function(hits_domain_dt, exon_features_dt) {
    hd <- as.data.table(hits_domain_dt)
    
    hd[, event_type   := first_available(hd, c("event_type_inc", "event_type_exc", "event_type"))]
    hd[, gene_id_plot := first_available(hd, c("gene_id", "gene_id_inc", "gene_id_exc"))]
    
    for (c in c("event_id", "transcript_id_inc", "transcript_id_exc")) {
      if (!c %in% names(hd)) hd[, (c) := NA_character_]
    }
    for (lc in c("inc_only_domains_list", "exc_only_domains_list", "either_domains_list")) {
      if (!lc %in% names(hd)) hd[, (lc) := vector("list", nrow(hd))]
    }
    
    expand_dir <- function(list_col, tx_col, direction_label) {
      out <- hd[, .(
        event_id,
        event_type,
        gene_id = gene_id_plot,
        transcript_id = get(tx_col),
        domain_vec = get(list_col)
      )]
      
      out <- out[!vapply(domain_vec, is.null, logical(1))]
      out <- out[, .(
        domain_id = as.character(unlist(domain_vec, use.names = FALSE))
      ), by = .(event_id, event_type, gene_id, transcript_id)]
      
      out <- out[!is.na(domain_id) & nzchar(domain_id)]
      out[, direction := direction_label]
      out[]
    }
    
    inc_long    <- expand_dir("inc_only_domains_list", "transcript_id_inc", "inclusion")
    exc_long    <- expand_dir("exc_only_domains_list", "transcript_id_exc", "exclusion")
    shared_long <- expand_dir("either_domains_list",   "transcript_id_inc", "shared")
    
    dom <- rbindlist(list(inc_long, exc_long, shared_long), use.names = TRUE, fill = TRUE)
    if (!nrow(dom)) return(dom)
    
    ef <- as.data.table(exon_features_dt)
    feat_lookup <- unique(ef[, .(feature_id, database, name, alt_name)])
    dom <- merge(dom, feat_lookup, by.x = "domain_id", by.y = "feature_id", all.x = TRUE)
    
    miss <- is.na(dom$database) | !nzchar(dom$database) | is.na(dom$name) | !nzchar(dom$name)
    if (any(miss)) {
      parsed <- tstrsplit(dom$domain_id[miss], ";", fixed = TRUE)
      if (length(parsed) >= 2) {
        dom[miss & (is.na(database) | !nzchar(database)), database := parsed[[1]]]
        dom[miss & (is.na(name) | !nzchar(name)), name := parsed[[2]]]
      }
    }
    
    dom[, database := fifelse(is.na(database) | !nzchar(database), "unlabeled", database)]
    dom[, name     := fifelse(is.na(name)     | !nzchar(name),     domain_id,   name)]
    dom[, alt_name := fifelse(is.na(alt_name), "", alt_name)]
    dom[]
  }
  
  standardize_ppi_cols <- function(ppi_dt) {
    dt <- as.data.table(ppi_dt)
    col_a <- intersect(c("ensembl_transcript_id", "txA", "transcript_id_A", "transcriptA"), names(dt))
    col_b <- intersect(c("i.ensembl_transcript_id", "partner_id", "partner_transcript_id", "txB", "transcript_id_B", "transcriptB"), names(dt))
    if (!length(col_a) || !length(col_b)) return(NULL)
    setnames(dt, col_a[1], "ensembl_transcript_id")
    setnames(dt, col_b[1], "i.ensembl_transcript_id")
    dt[, `:=`(
      ensembl_transcript_id = as.character(ensembl_transcript_id),
      i.ensembl_transcript_id = as.character(i.ensembl_transcript_id)
    )]
    dt <- dt[nzchar(ensembl_transcript_id) & nzchar(i.ensembl_transcript_id)]
    unique(dt[, .(ensembl_transcript_id, i.ensembl_transcript_id)])
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
      downstream_results(NULL)
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
      interpro_features <- get_protein_features(c("interpro"), ann$annotations, test = TRUE)
      signalp_features <- get_protein_features(c("signalp"), ann$annotations, test = TRUE)
      pf <- get_comprehensive_annotations(list(signalp_features, interpro_features))
      pf_dt <- as.data.table(pf)
      if (!"ensembl_transcript_id" %in% names(pf_dt) && "transcript_id" %in% names(pf_dt)) {
        pf_dt[, ensembl_transcript_id := transcript_id]
      }
      if (!"transcript_id" %in% names(pf_dt) && "ensembl_transcript_id" %in% names(pf_dt)) {
        pf_dt[, transcript_id := ensembl_transcript_id]
      }
      if (!"gene_id" %in% names(pf_dt) && "gene_id.x" %in% names(pf_dt)) pf_dt[, gene_id := gene_id.x]
      rv$protein_features <- pf_dt
      ef <- get_exon_features(ann$annotations, pf_dt)
      if (!"transcript_id" %in% names(ef) && "ensembl_transcript_id" %in% names(ef)) ef[, transcript_id := ensembl_transcript_id]
      if ("gene_id.x" %in% names(ef) && !"gene_id" %in% names(ef)) setnames(ef, "gene_id.x", "gene_id")
      rv$exon_features <- ef
      
      incProgress(0.5, detail = "Event mapping")
      rv$di_sig <- NULL
      rv$matched <- NULL
      
      incProgress(0.65, detail = "Differential inclusion (demo)")
      rv$di <- get_differential_inclusion(
        rv$splicing,
        min_total_reads = input$min_reads,
        cooks_cutoff = "4/n",
        parallel_glm = FALSE,
        verbose = FALSE
      )
      rv$di_norm <- normalize_di_cols(rv$di)
      rv$di_sig <- if (!is.null(rv$di_norm)) keep_sig_pairs(rv$di_norm, padj_thr = input$padj_thr, dpsi_thr = input$dpsi_thr) else NULL
      if (!is.null(rv$di_sig)) {
        matched <- get_matched_events_chunked(rv$di_sig, rv$annotations$annotations, chunk_size = 2000)
        rv$matched <- as.data.table(matched)
      }
      
      incProgress(0.8, detail = "PPI (demo PPIDM)")
      rv$ppidm <- get_ppidm(test = TRUE)
      
      res <- tryCatch(
        compute_downstream(show_progress = FALSE),
        error = function(e) {
          showNotification(paste("Downstream results error:", e$message), type = "error")
          NULL
        }
      )
      downstream_results(res)
      
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
      rv$di_norm <- NULL
      rv$di_sig <- NULL
      rv$matched <- NULL
      rv$protein_features <- NULL
      rv$exon_features <- NULL
      rv$ppidm <- NULL
      downstream_results(NULL)
      incProgress(1, detail = "Events loaded")
    })
  })
  
  observeEvent(input$load_proteins, {
    req(rv$annotations)
    withProgress(message = "Loading protein features", value = 0, {
      pf <- NULL
      if (isTruthy(input$use_demo_proteins)) {
        interpro_features <- get_protein_features(
          biomaRt_databases = c("interpro"),
          gtf_df = rv$annotations$annotations,
          test = TRUE
        )
        signalp_features <- get_protein_features(
          biomaRt_databases = c("signalp"),
          gtf_df = rv$annotations$annotations,
          test = TRUE
        )
        pf <- get_comprehensive_annotations(list(signalp_features, interpro_features))
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
      if (!"ensembl_transcript_id" %in% names(pf_dt) && "transcript_id" %in% names(pf_dt)) {
        pf_dt[, ensembl_transcript_id := transcript_id]
      }
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
      rv$ppidm <- NULL
      downstream_results(NULL)
      incProgress(1, detail = "Protein features ready")
    })
  })
  
  observeEvent(input$load_ppi, {
    withProgress(message = "Loading PPI interactions", value = 0, {
      ppidm <- if (isTRUE(input$use_demo_ppi)) {
        get_ppidm(test = TRUE)
      } else {
        get_ppidm(download = TRUE)
      }
      rv$ppidm <- ppidm
      incProgress(1, detail = "PPI metadata ready")
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
      rv$di_norm <- normalize_di_cols(di)
      rv$di_sig <- if (!is.null(rv$di_norm)) keep_sig_pairs(rv$di_norm, padj_thr = input$padj_thr, dpsi_thr = input$dpsi_thr) else NULL
      rv$matched <- NULL
      incProgress(1, detail = "Differential inclusion complete")
    })
  })
  
  observeEvent(input$map_events, {
    req(rv$di_norm, rv$annotations)
    withProgress(message = "Matching events to transcripts", value = 0, {
      rv$di_sig <- keep_sig_pairs(rv$di_norm, padj_thr = input$padj_thr, dpsi_thr = input$dpsi_thr)
      req(rv$di_sig)
      matched <- get_matched_events_chunked(rv$di_sig, rv$annotations$annotations, chunk_size = 2000)
      rv$matched <- as.data.table(matched)
      updateSelectizeInput(session, "probe_event", choices = unique(rv$di_norm$event_id), server = TRUE)
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
    req(rv$di_norm)
    sig <- rv$di_norm
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
    req(rv$matched)
    tmp <- unique(as.data.table(rv$matched)[, .(event_id, event_type, gene_id, transcript_id)])
    tmp <- apply_filters(tmp, gene_col = "gene_id", tx_col = "transcript_id", gene_override = input$map_gene_filter)
    tmp[, .(
      n_events = uniqueN(event_id),
      n_transcripts = uniqueN(transcript_id),
      n_event_tx_pairs = .N
    ), by = event_type][order(-n_events)]
  })
  
  output$mapping_gene_table <- renderTable({
    req(rv$matched)
    gene_filter_active <- nzchar(input$map_gene_filter) || nzchar(input$gene_filter)
    if (!gene_filter_active) return(NULL)
    keep <- c("event_id", "event_type", "gene_id", "transcript_id", "exons", "inc_exons_by_idx")
    dat <- apply_filters(rv$matched[, ..keep], gene_col = "gene_id", tx_col = "transcript_id", gene_override = input$map_gene_filter)
    if (!nrow(dat)) return(NULL)
    dat[order(event_id)][1:min(.N, 50)]
  })
  
  observeEvent(input$show_mapping_table, {
    req(rv$matched)
    keep <- c("event_id", "event_type", "gene_id", "transcript_id", "exons", "inc_exons_by_idx")
    dat <- apply_filters(rv$matched[, ..keep], gene_col = "gene_id", tx_col = "transcript_id", gene_override = input$map_gene_filter)
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
    res <- tryCatch(
      downstream_results(),
      error = function(e) {
        showNotification(paste("Downstream results error:", e$message), type = "error")
        NULL
      }
    )
    if (is.null(res) || is.null(res$domain_long) || !nrow(res$domain_long)) {
      return(data.table())
    }
    dt <- as.data.table(res$domain_long)
    needed <- c("event_id", "event_type", "gene_id", "transcript_id", "direction", "domain_id", "database", "name", "alt_name")
    missing_cols <- setdiff(needed, names(dt))
    if (length(missing_cols)) dt[, (missing_cols) := NA_character_]
    out <- dt[, ..needed]
    for (col in needed) set(out, j = col, value = as.character(out[[col]]))
    as.data.table(out)
  })
  
  output$protein_plot <- renderPlot({
    cons <- protein_consequences()
    if (!nrow(cons)) {
      return(NULL)
    }
    cons <- apply_filters(cons, gene_col = "gene_id", tx_col = "transcript_id", prot_col = "ensembl_peptide_id", gene_override = input$prot_gene_filter)
    if (!nrow(cons)) return(NULL)
    agg <- cons[, .N, by = .(event_type, direction, database)][order(event_type)]
    ggplot(agg, aes(x = event_type, y = N, fill = database)) +
      geom_col(position = "stack") +
      facet_wrap(~direction, ncol = 1) +
      labs(x = "Event type", y = "Domain changes", fill = "Database") +
      theme_classic()
  })
  
  output$protein_summary <- renderTable({
    cons <- protein_consequences()
    if (!nrow(cons)) return(NULL)
    cons <- apply_filters(cons, gene_col = "gene_id", tx_col = "transcript_id", prot_col = "ensembl_peptide_id", gene_override = input$prot_gene_filter)
    if (!nrow(cons)) return(NULL)
    cons[, .(
      events = uniqueN(event_id),
      transcripts = uniqueN(transcript_id),
      domains = uniqueN(domain_id)
    ), by = .(database, direction)][order(-domains)]
  })
  
  output$protein_gene_table <- renderTable({
    cons <- protein_consequences()
    gene_filter_active <- nzchar(input$prot_gene_filter) || nzchar(input$gene_filter)
    if (!gene_filter_active || !nrow(cons)) return(NULL)
    cons <- apply_filters(cons, gene_col = "gene_id", tx_col = "transcript_id", prot_col = "ensembl_peptide_id", gene_override = input$prot_gene_filter)
    if (!nrow(cons)) return(NULL)
    cols <- c("event_id", "event_type", "gene_id", "transcript_id", "direction", "database", "name", "alt_name")
    cons[, ..cols][order(event_id)][1:min(.N, 50)]
  })
  
  observeEvent(input$show_protein_table, {
    cons <- protein_consequences()
    req(nrow(cons))
    cons <- apply_filters(cons, gene_col = "gene_id", tx_col = "transcript_id", prot_col = "ensembl_peptide_id", gene_override = input$prot_gene_filter)
    cols <- c("event_id", "event_type", "gene_id", "transcript_id",
              "direction", "database", "name", "alt_name")
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
      req(rv$di_norm)
      fwrite(rv$di_norm, file)
    }
  )
  
  output$download_mapping <- downloadHandler(
    filename = function() paste0("spliceimpactr_event_mapping_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$matched)
      fwrite(rv$matched, file)
    }
  )
  
  output$download_proteins <- downloadHandler(
    filename = function() paste0("spliceimpactr_protein_consequences_", Sys.Date(), ".csv"),
    content = function(file) {
      cons <- protein_consequences()
      req(nrow(cons))
      fwrite(cons, file)
    }
  )
  
  ppi_results <- reactive({
    res <- downstream_results()
    if (is.null(res) || is.null(res$hits_final) || !nrow(res$hits_final)) {
      return(data.table(
        event_id = character(), event_type = character(), gene_for_plot = character(),
        transcript_id_inc = character(), transcript_id_exc = character(),
        n_inc_ppi = integer(), n_exc_ppi = integer(), n_ppi = integer(),
        inc_ppi = I(list()), exc_ppi = I(list())
      ))
    }
    hf <- as.data.table(res$hits_final)
    hf[, event_type := first_available(hf, c("event_type_inc", "event_type_exc", "event_type"))]
    hf[, gene_for_plot := first_available(hf, c("gene_id", "gene_id_inc", "gene_id_exc"))]
    if (!"gene_id" %in% names(hf)) hf[, gene_id := gene_for_plot]
    if (nzchar(input$transcript_filter)) {
      hf <- hf[grepl(input$transcript_filter, transcript_id_inc, ignore.case = TRUE) |
                 grepl(input$transcript_filter, transcript_id_exc, ignore.case = TRUE)]
    }
    hf <- apply_filters(hf, gene_col = "gene_for_plot", tx_col = NULL, gene_override = input$ppi_gene_filter)
    hf
  })
  
  ppi_plot_obj <- reactive({
    hf <- ppi_results()
    if (!nrow(hf)) return(NULL)
    plot_ppi_summary(hf)
  })
  
  output$ppi_plot <- renderPlot({
    plt <- ppi_plot_obj()
    if (is.null(plt)) {
      showNotification("Run sequence/domain plots after loading PPIs to summarize PPI switches.", type = "message")
      return(NULL)
    }
    plt
  })
  
  output$ppi_summary <- renderTable({
    hf <- ppi_results()
    if (!nrow(hf)) {
      showNotification("Run sequence/domain plots with PPIs loaded to populate PPI results.", type = "message")
      return(NULL)
    }
    hf[, .(
      events = .N,
      impacted = sum(n_ppi > 0, na.rm = TRUE),
      mean_ppi = round(mean(n_ppi, na.rm = TRUE), 2)
    ), by = event_type][order(-events)]
  })
  
  output$ppi_gene_table <- renderTable({
    hf <- ppi_results()
    gene_filter_active <- nzchar(input$ppi_gene_filter) || nzchar(input$gene_filter)
    if (!gene_filter_active || !nrow(hf)) return(NULL)
    hf[, .(
      event_id,
      event_type,
      gene_id = gene_for_plot,
      transcript_id_inc,
      transcript_id_exc,
      n_inc_ppi,
      n_exc_ppi,
      n_ppi,
      inc_ppi = vapply(inc_ppi, collapse_list_vals, character(1)),
      exc_ppi = vapply(exc_ppi, collapse_list_vals, character(1))
    )][order(event_id)][1:min(.N, 50)]
  })
  
  observeEvent(input$show_ppi_table, {
    hf <- ppi_results()
    req(nrow(hf))
    dat <- hf[, .(
      event_id,
      event_type,
      gene_id = gene_for_plot,
      transcript_id_inc,
      transcript_id_exc,
      n_inc_ppi,
      n_exc_ppi,
      n_ppi,
      inc_ppi = vapply(inc_ppi, function(x) collapse_list_vals(x, n = 10), character(1)),
      exc_ppi = vapply(exc_ppi, function(x) collapse_list_vals(x, n = 10), character(1))
    )][order(event_id)]
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
  
  compute_downstream <- function(show_progress = TRUE) {
    if (is.null(rv$matched) || is.null(rv$annotations) || is.null(rv$exon_features) || is.null(rv$protein_features)) {
      showNotification("Load annotations, proteins, and mapped events before running downstream analyses.", type = "warning")
      return(NULL)
    }
    if (is.null(rv$annotations$sequences)) {
      showNotification("Annotation sequences are required for downstream plots.", type = "error")
      return(NULL)
    }
    
    run_block <- function(inc_fn) {
      map_dt <- unique(as.data.table(rv$matched), by = c("event_id", "transcript_id"))
      if (!is.null(rv$protein_features) && "ensembl_transcript_id" %in% names(rv$protein_features)) {
        keep_tx <- unique(rv$protein_features$ensembl_transcript_id)
        map_dt <- map_dt[transcript_id %in% keep_tx]
      }
      if (!nrow(map_dt)) {
        showNotification("No mapped events available for downstream analysis after filtering to protein-annotated transcripts.", type = "warning")
        return(NULL)
      }
      
      inc_fn(0.2, detail = "Attaching sequences")
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
      
      inc_fn(0.4, detail = "Comparing sequences")
      seq_compare <- tryCatch(
        compare_sequence_frame(pairs, rv$annotations$annotations),
        error = function(e) {showNotification(paste("compare_sequence_frame failed:", e$message), type = "error"); NULL})
      seq_compare <- as_dt_or_null(seq_compare)
      if (!has_df_rows(seq_compare)) return(NULL)
      proximal_output <- tryCatch(get_proximal_shift_from_hits(pairs),
                                  error = function(e) {showNotification(paste("get_proximal_shift_from_hits failed:", e$message), type = "error"); NULL})
      proximal_plot <- if (!is.null(proximal_output) && nrow(proximal_output) > 0) {
        tryCatch(plot_prox_dist(proximal_output), error = function(e) NULL)
      } else NULL
      alignment_plot <- tryCatch(plot_alignment_summary(seq_compare), error = function(e) NULL)
      length_plot <- tryCatch(plot_length_comparison(seq_compare), error = function(e) NULL)
      
      domain_plot <- NULL
      enriched_domains <- NULL
      domain_summary <- NULL
      domain_long <- data.table()
      hits_domain <- NULL
      inc_fn(0.6, detail = "Domain overlaps")
      hits_domain <- tryCatch(
        get_domains(seq_compare, rv$exon_features),
        error = function(e) {showNotification(paste("get_domains failed:", e$message), type = "error"); NULL})
      hits_domain <- as_dt_or_null(hits_domain)
      if (has_df_rows(hits_domain) && !is.null(rv$sample_frame) && has_df_rows(rv$protein_features)) {
        bg <- tryCatch(get_background(
          source = "annotated",
          annotations = rv$annotations$annotations,
          protein_features = rv$protein_features
        ), error = function(e) {showNotification(paste("get_background failed:", e$message), type = "error"); NULL})
        bg <- as_dt_or_null(bg)
        if (has_df_rows(bg)) {
          enriched_domains <- tryCatch(enrich_domains_hypergeo(hits_domain, bg, db_filter = "interpro"),
                                       error = function(e) NULL)
          enriched_domains <- as_dt_or_null(enriched_domains)
          if (has_df_rows(enriched_domains)) {
            domain_plot <- tryCatch(plot_enriched_domains_counts(enriched_domains, top_n = 20),
                                    error = function(e) NULL)
          }
        }
      }
      
      if (has_df_rows(hits_domain)) {
        hd <- as.data.table(hits_domain)
        hd[, event_type := first_available(hd, c("event_type_inc", "event_type_exc", "event_type"))]
        hd[, gene_for_plot := first_available(hd, c("gene_id", "gene_id_inc", "gene_id_exc"))]
        
        domain_long <- tryCatch(
          build_domain_long(hd, rv$exon_features),
          error = function(e) {
            showNotification(paste("Domain-long construction failed:", e$message), type = "error")
            data.table()
          }
        )
        domain_summary <- hd[, .(
          events = .N,
          with_domain_change = sum(as.numeric(diff_n) > 0, na.rm = TRUE),
          inc_only = sum(as.numeric(inc_only_n), na.rm = TRUE),
          exc_only = sum(as.numeric(exc_only_n), na.rm = TRUE)
        ), by = event_type][order(-with_domain_change)]
      }
      
      hits_final <- NULL
      ppi_plot <- NULL
      if (!is.null(rv$ppidm) && has_df_rows(hits_domain)) {
        restrict_ids <- unique(c(hits_domain$transcript_id_inc, hits_domain$transcript_id_exc))
        restrict_pf <- as.data.table(rv$protein_features)[ensembl_transcript_id %in% restrict_ids]
        if (!nrow(restrict_pf)) {
          showNotification("No protein features matched domain hits; cannot compute PPIs.", type = "warning")
        } else {
          ppi_raw <- tryCatch(
            get_isoform_interactions(restrict_pf, rv$ppidm, init = TRUE, save = FALSE),
            error = function(e) {showNotification(paste("get_isoform_interactions failed:", e$message), type = "error"); NULL}
          )
          ppi_std <- if (!is.null(ppi_raw)) standardize_ppi_cols(ppi_raw) else NULL
          if (!is.null(ppi_std) && nrow(ppi_std) > 0) {
            hits_final <- tryCatch(get_ppi_switches(hits_domain, ppi_std),
                                   error = function(e) {
                                     showNotification(paste("get_ppi_switches failed:", e$message), type = "error")
                                     NULL
                                   })
            hits_final <- as_dt_or_null(hits_final)
            if (has_df_rows(hits_final)) {
              ppi_plot <- tryCatch(plot_ppi_summary(hits_final), error = function(e) NULL)
            }
          }
        }
      } else if (is.null(rv$ppidm)) {
        showNotification("Load PPI interactions to compute PPI switches alongside domain hits.", type = "message")
      }
      
      inc_fn(1, detail = "Done")
      
      list(
        seq_compare = seq_compare,
        alignment_plot = alignment_plot,
        length_plot = length_plot,
        enriched_domains = enriched_domains,
        domain_plot = domain_plot,
        hits_domain = hits_domain,
        domain_long = domain_long,
        domain_summary = domain_summary,
        hits_final = hits_final,
        ppi_plot = ppi_plot,
        pairs = pairs,
        proximal_output = proximal_output,
        proximal_plot = proximal_plot
      )
    }
    
    if (isTRUE(show_progress)) {
      withProgress(message = "Running sequence/domain plots", value = 0, {
        run_block(incProgress)
      })
    } else {
      inc_stub <- function(amount = 0, detail = NULL) {}
      run_block(inc_stub)
    }
  }
  
  observeEvent(input$run_downstream, {
    res <- tryCatch(
      compute_downstream(show_progress = TRUE),
      error = function(e) {
        showNotification(paste("Downstream results error:", e$message), type = "error")
        NULL
      }
    )
    downstream_results(res)
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
  
  output$proximal_plot <- renderPlot({
    res <- downstream_results()
    if (is.null(res) || is.null(res$proximal_plot)) return(NULL)
    res$proximal_plot
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
    if (is.null(dr) || is.null(dr$hits_final) || !nrow(dr$hits_final)) {
      showNotification("Run sequence/domain plots with PPIs loaded to compute hits_final before the integrative summary.", type = "error")
      return(NULL)
  }
  di_dt <- normalize_di_cols(rv$di)
  if (is.null(di_dt)) return(NULL)
  
  int <- tryCatch(integrated_event_summary(dr$hits_final, di_dt),
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
  
  hit_overview <- eventReactive(input$run_hit_overview + input$run_hit_overview_sidebar, {
    req(rv$sample_frame, rv$splicing)
    hc <- tryCatch(
      compare_hit_index(rv$sample_frame, condition_map = c(control = "control", test = "case")),
      error = function(e) {showNotification(paste("compare_hit_index failed:", e$message), type = "error"); NULL}
    )
    ov <- tryCatch(
      overview_splicing_comparison_fixed(rv$splicing, rv$sample_frame, depth_norm = "exon_files", event_type = "AFE"),
      error = function(e) {showNotification(paste("overview_splicing_comparison_fixed failed:", e$message), type = "error"); NULL}
    )
    list(hit = hc, overview = ov)
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
    updateSelectizeInput(session, "probe_event", choices = character(0), server = TRUE)
  }, ignoreNULL = TRUE)
  
  probe_choices <- reactive({
    di_events <- if (!is.null(rv$di_norm)) unique(rv$di_norm$event_id) else character(0)
    hit_events <- if (!is.null(rv$splicing)) unique(rv$splicing$event_id) else character(0)
    choices <- intersect(di_events, hit_events)
    if (!length(choices) && length(di_events)) choices <- di_events
    sort(unique(choices))
  })
  
  observeEvent(probe_choices(), {
    choices <- probe_choices()
    current <- isolate(input$probe_event)
    selected <- if (!is.null(current) && nzchar(current) && current %chin% choices) current else NULL
    updateSelectizeInput(session, "probe_event", choices = choices, server = TRUE, selected = selected)
  }, ignoreInit = TRUE)
  
  observeEvent(input$pair_gene, {
    req(rv$annotations)
    ann <- as.data.table(rv$annotations$annotations)
    tx_choices <- ann[type == "transcript" & gene_id == input$pair_gene, unique(transcript_id)]
    updateSelectizeInput(session, "tx_a", choices = tx_choices, server = TRUE, selected = NULL)
    updateSelectizeInput(session, "tx_b", choices = tx_choices, server = TRUE, selected = NULL)
  }, ignoreInit = TRUE)
  
  probe_event_plot <- eventReactive(input$run_probe, {
    req(rv$di_norm, rv$splicing)
    evt <- input$probe_event
    if (!nzchar(evt)) {
      showNotification("Select an event ID to probe.", type = "warning")
      return(NULL)
    }
    allowed <- intersect(unique(rv$di_norm$event_id), unique(rv$splicing$event_id))
    if (!evt %chin% allowed) {
      showNotification("Pick an event present in both DI results and the splicing table.", type = "warning")
      return(NULL)
    }
    tryCatch(
      probe_individual_event(rv$splicing, event = evt),
      error = function(e) {
        showNotification(paste("probe_individual_event failed:", e$message), type = "error")
        NULL
      }
    )
  })
  
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
  
  output$hit_compare_plot <- renderPlot({
    res <- hit_overview()
    if (is.null(res) || is.null(res$hit$plot)) return(NULL)
    res$hit$plot
  })
  
  output$psi_overview_plot <- renderPlot({
    res <- hit_overview()
    if (is.null(res) || is.null(res$overview)) return(NULL)
    res$overview
  })
  
  output$probe_plot <- renderPlot({
    probe_event_plot()
  })
}

shinyApp(ui, server)
