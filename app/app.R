library(shiny)
library(shinydashboard)
library(shinythemes)
library(plotly)
library(reactable)
library(shiny)
library(tidyverse)
library(reactable)
library(plotly)
library(here)
library(data.table)
library(scales)
library(speedyseq)
library(crosstalk)
library(gggenomes)

# --------  Read data ----------
source(here("script/import_data.R"))
source(here("script/utils.R"))
source(here("script/annotation_prot.R"))
load(here("data/vpf.RData"))

ctgmeta <- ctganno_merged
# filter(!is.na(virsorter_category))
smeta <- data.frame(sample_data(pseq)) %>% 
  mutate(samplename=rownames(.))
abundance <- data.frame(otu_table(pseq))

minlen <- min(ctganno_merged$contig_length)
maxlen <- max(ctganno_merged$contig_length)


dpath <- here("data/hms")
fp_amg <- here("data/raw/amg_database.tsv")
fp_anno_prot <- here(dpath, "annotations.tsv")
ganno <- read_anno_vprot(fp_anno_prot, fp_amg)


# ---------- UI -----------
# header ----
header <- dashboardHeader(title = "P0062: VirOil", titleWidth = 300)
sidebar <- dashboardSidebar(
  fileInput("upload", "Upload vpf.RData file"),
  uiOutput("contig_len"),
  radioButtons("logic", "Filter logic", 
               choices = c("AND", "OR"),
               selected = "AND"),
  uiOutput("vs2cat"),
  uiOutput("checkv_quality"),
  uiOutput("ctg_feature_col"),
  uiOutput("prot_feature_col"),
  uiOutput("sample_metadata"),
  radioButtons("abundance_scale", "Absolute/relative abundance",
               choices = c("Relative", "Absolute"),
               selected = "Absolute"),
  numericInput("abundance_min_threshold", "Minimum abundance to show",
               value = 0.001, min=0),
  width=300
)


# box ---------
box2 <- box(
  plotlyOutput("histplot"), 
  title = "Viral quality",
  status = "primary", 
  solidHeader = TRUE, 
  collapsible = TRUE,
  width = 3)
box3 <- box(
  plotlyOutput("taxa_abundance"), 
  title = "Reads abundance", 
  status = "primary", 
  solidHeader = TRUE, 
  collapsible = TRUE,
  width = 3)
box4 <- box(
  plotOutput("alpha_diversity"), 
  title = "Diversity",
  status = "primary", 
  solidHeader = TRUE, 
  collapsible = TRUE)
box_genome1 <- box(
  uiOutput("genomeSel"),
  uiOutput("geneFeature"),
  uiOutput("geneLabel"),
  title = "Select Genome",
  collapsible = TRUE,
  width = 2)
box_genome2 <- box(
  plotOutput("genomePlot"),
  reactableOutput("geneAnno"),
  title = "Genome Annotation",
  status = "primary", 
  solidHeader = TRUE, 
  collapsible = TRUE,
  width = 10)
tabBox_abundance <- tabBox(
  title = "Abundance",
  id = "tset_abundance",
  tabPanel("Count", reactableOutput("count")),
  tabPanel("TPM", reactableOutput("tpm"))
)
tabBox_ctganno <- tabBox(
  title = "Contig annotation",
  id = "tset_ctganno",
  tabPanel("ctganno", reactableOutput("results"))
)
tabBox_anno_prot <- tabBox(
  title = "Protein annotation",
  id = "tset_anno_prot",
  tabPanel("annovprot", reactableOutput("anno_vprot")),
  tabPanel("annoprot", reactableOutput("anno_prot"))
)


# --------- body -----------
body <- dashboardBody(
  fluidRow(tabBox_anno_prot, tabBox_abundance, box2, box3, box4, box_genome1, box_genome2)
)




# ------- Server ---------
server <- function(input, output) {
  # render UI ---------
  output$contig_len <- renderUI({
    sliderInput("contig_len", "Contig length",
                min=minlen, 
                max=maxlen, 
                step = 500,
                ticks = FALSE,
                post = " bp",
                value=c(minlen,maxlen))
  })
  
  output$checkv_quality <- renderUI({
    selectInput("checkv_quality", "CheckV quality",
                sort(unique(ctgmeta$checkv_quality)),
                selected = c("Complete", "High-quality"),
                multiple = TRUE)
  })
  
  output$vs2cat <- renderUI({
    selectInput("vs2cat", "VirSorter2 category",
                sort(unique(ctgmeta$virsorter_category)),
                selected = c("1", "2"),
                multiple = TRUE)
  })
  
  output$sample_metadata <- renderUI({
    selectInput("sample_metadata", "Sample metadata",
                colnames(smeta),
                selected = "samplename")
  })
  
  output$ctg_feature_col <- renderUI({
    selectInput("ctg_feature_col", "Feature column",
                colnames(ctgmeta),
                selected = "Family")
  })
  
  output$prot_feature_col <- renderUI({
    selectInput("prot_feature_col", "Protein annotation column",
                colnames(anno_vprot),
                selected = "ko_id")
  })
  
  output$genomeSel <- renderUI({
    selectInput("genomeSel", "Select viral contig",
                fgbks,
                selected = "DB_S2__NODE_13_length_111139_cov_10.702044-cat_1")
  })
  
  output$geneFeature <- renderUI({
    selectInput("geneFeature", "Select gene annotation",
                colnames(ganno),
                selected = "auxiliary_score")
  })
  
  output$geneLabel <- renderUI({
    selectInput("geneLabel", "Show gene label",
                colnames(ganno),
                selected = "kegg_hit")
  })
  
  # --------- reactive --------------
  
  # -------- filter contigs ---------
  ctgmeta_filtered <- reactive({
    if (is.null(input$checkv_quality)) {
      return(NULL)
    }
    if (input$logic=="AND") {
      ctgmeta %>% 
        filter(contig_length >= input$contig_len[1],
               contig_length <= input$contig_len[2]) %>% 
        filter(virsorter_category %in% input$vs2cat) %>%
        filter(checkv_quality %in% input$checkv_quality)
    } else if (input$logic=="OR") {
      ctgmeta %>% 
        filter(contig_length >= input$contig_len[1],
               contig_length <= input$contig_len[2]) %>% 
        filter(virsorter_category %in% input$vs2cat | 
                 checkv_quality %in% input$checkv_quality)
    }
  })
  
  # --------- load genbank --------------
  gbk <- reactive({
    fgbk <- here(fpath_gbk, paste0(input$genomeSel, ".gbk"))
    gggenomes::read_gbk(fgbk) %>% 
      mutate(gene=name) %>% 
      left_join(ganno, by="gene")
  })
  
  # ------ checkv quality --------
  output$histplot <- renderPlotly({
    if (is.null(ctgmeta_filtered())) {
      return()
    }
    p <- ctgmeta_filtered() %>% 
      ggplot(aes(x=contig_length, color=checkv_quality, fill=checkv_quality)) +
      geom_histogram(alpha=0.5, position = "identity", bins=30) +
      scale_x_continuous(trans=log10_trans(),
                         breaks = trans_breaks("log10",
                                               function(x) as.integer(10^x))) +
      theme_bw()
    ggplotly(p)
  })
  
  # ------ subset Phyloseq -------
  pseq_subset <- reactive({
    filter_tax_table(pseq, Contig %in% ctgmeta_filtered()$Contig)
    # filter_tax_table(pseq, Contig %in% ctgmeta_shared$key())
  })
  
  # --------- taxonomy bar plot ---------
  output$taxa_abundance <- renderPlotly({
    # if (is.null(ctgmeta_filtered())) {
    #   return()
    # }
    taxa_collapse <-  tax_glom(pseq_subset(), taxrank = input$ctg_feature_col)
    if (input$abundance_scale=="Relative") {
      taxa_collapse <- taxa_collapse %>% 
        transform_sample_counts(function(x) {x/sum(x)})
    }
    p <- taxa_collapse %>% 
      psmelt() %>% 
      filter(Abundance > input$abundance_min_threshold) %>%
      ggplot(aes(x = Sample, y = Abundance, fill=.data[[input$ctg_feature_col]])) +    # Color by Phylum
      geom_bar(stat = "identity", position="stack") +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab(paste0(input$abundance_scale, " abundance (> ", input$abundance_min_threshold, ")")) +
      theme_bw()
    
    # theme_bw()
    ggplotly(p)
  })
  
  # ------- feature annotation table --------
  output$results <- renderReactable({
    reactable(ctgmeta_filtered(),
              selection = "multiple", 
              filterable = TRUE,
              onClick = "select",
              searchable = TRUE)
  })
  
  # -------- alpha diversity -----------
  output$alpha_diversity <- renderPlot({
    plot_richness(pseq_subset())
    # ggplotly(p)
  })
  
  # ----------- genome plot -----------
  output$genomePlot <- renderPlot({
    gbk() %>% 
      gggenomes() +
      geom_seq() +
      geom_seq_label(size=6) +
      geom_gene(aes(fill=.data[[input$geneFeature]], label=input$geneLabel), position="strand") +
      geom_gene_tag(aes(label=.data[[input$geneLabel]]), size=4) +
      geom_text(aes(label=input$geneLabel, x=-25, y=1.1), size=6) +
      theme(legend.position = "bottom",
            legend.text=element_text(size=16),
            text=element_text(size=16),
            axis.text.x = element_text(size=16))
  })
  
  # ----------- gene annotation ------------
  output$geneAnno <- renderReactable({
    reactable(ganno,
              selection = "multiple", 
              filterable = TRUE,
              onClick = "select")
  })
  
  output$anno_vprot <- renderReactable({
    reactable(ganno,
              selection = "multiple", 
              filterable = TRUE,
              onClick = "select",
              columns = list(
              ko_id = colDef(html = TRUE, cell = function(x) {sprintf('<a href="https://www.genome.jp/dbget-bin/www_bget?ko:%s" target="_blank">%s</a>', x, x)}),
              viral_id = colDef(html = TRUE, cell = function(x) {sprintf('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=%s" target="_blank">%s</a>', x, x)}),
              pfam_hits = colDef(html = TRUE, cell = function(x) {sprintf('<a href="http://pfam.xfam.org/family/%s" target="_blank">%s</a>', str_extract(x, "PF[0-9.]+"), x)}),
              vogdb_id = colDef(html = TRUE, cell = function(x) {sprintf('<a href="https://vogdb.org/reports/vog_report?id=%s" target="_blank">%s</a>', x, x)}),
              vogdb_hits = colDef(html = TRUE, cell = function(x) {sprintf('<a href="https://www.uniprot.org/uniprot/%s" target="_blank">%s</a>', str_replace_all(str_extract(x, "\\|.*\\|"), "\\|", ""), x)})
              ))
  })
  
  # --------- abundance ---------
  output$count <- renderReactable({
    reactable(df_count, filterable = TRUE, onClick = "select", resizable = TRUE)
  })
  
  output$tpm <- renderReactable({
    reactable(df_tpm, filterable = TRUE, onClick = "select", resizable = TRUE)
  })
}

# Run the application 
shinyApp(
  ui = dashboardPage(header, sidebar, body), 
  server = server)
