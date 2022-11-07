library(shiny)
library(here)
library(tidyverse)
library(reactable)
library(data.table)
library(shinythemes)
library(plotly)
library(heatmaply)
library(pheatmap)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
source(here("R/utils.R"))

dpath <- here("pipeline")

fp_clst_ctg <- here(dpath, "wtp/probe/clstmeta_p95_s95_c80_cluster_ani_raw.tsv")
fp_taxa_krakenstd <- here(dpath, "metaall/taxonomy/viroil_taxa_krakenstd.tsv")
fp_taxa_krakenvir <- here(dpath, "metaall/taxonomy/viroil_taxa_krakenvir.tsv")
fp_taxa_mvref <- here(dpath, "metaall/taxonomy/viroil_taxa_mmseqsvref.tsv")
fp_taxa_vc2 <- here(dpath, "metaall/taxonomy/taxa_vc2.tsv")
fp_taxa_vc2ref <- here(dpath, "metaall/taxonomy/taxa_vc2_ref.tsv")
fp_taxa_vc2all <- here(dpath, "metaall/taxonomy/genome_by_genome_overview.csv")
fp_virus_virify <- here(dpath, "metaall/virify/08-final/virus_virify.tsv")
fp_virus_wtp <- here(dpath, "metaall/wtp/final_report.utf8.csv")
fp_clsf_pprmeta <- here(dpath, "metaall/virify/01-viruses/pprmeta/repctg_pprmeta.csv")
fp_map_virify <- here(dpath, "metaall/virify/virify_mapping.tsv")
fp_map_gene <- here(dpath, "metaall/ganno/qry_cluster.tsv")
fp_ganno_canthyd <- here(dpath, "metaall/ganno/canthyd.tsv")
fp_ganno_nr <- here(dpath, "metaall/ganno/rstAnnoNR.tsv")
fp_host_genome <- here(dpath, "metaall/gbk/Host_prediction_to_genome_m90.csv")
fp_host_genus <- here(dpath, "metaall/gbk/Host_prediction_to_genus_m90.csv")

df_clst_ctg <- fread(fp_clst_ctg, header = F, col.names = c("repctg", "members"))
df_taxa_krakenstd <- fread(fp_taxa_krakenstd)
df_taxa_krakenvir <- fread(fp_taxa_krakenvir)
df_taxa_mvref <- fread(fp_taxa_mvref)
df_taxa_vc2 <- fread(fp_taxa_vc2)
df_taxa_vc2ref <- fread(fp_taxa_vc2ref)
df_taxa_vc2all <- fread(fp_taxa_vc2all)
df_virus_virify <- fread(fp_virus_virify)
df_virus_wtp <- fread(fp_virus_wtp) %>%
  dplyr::select(-V1)
df_clsf_pprmeta <- parse_pprmeta(fp_clsf_pprmeta, fp_map_virify)
df_ganno_canthyd <- parse_ganno(fp_ganno_canthyd, fp_map_gene)
df_ganno_nr <- parse_ganno(fp_ganno_nr, fp_map_gene)
df_host_genome <- fread(fp_host_genome)
df_host_genus <- fread(fp_host_genus)


df_clst_ctg_seprows <- df_clst_ctg %>%
  mutate(tplg = ifelse(str_detect(members, "circular"), "circular", "linear")) %>%
  mutate(repctg = ifelse(str_detect(repctg, "circular"), repctg, paste0(repctg, "_", tplg))) %>%
  mutate(repctg = ifelse(str_detect(repctg, "_linear"), str_replace(repctg, "_linear", ""), repctg)) %>%
  dplyr::select(-tplg) %>%
  separate_rows(members, sep = ",") %>%
  mutate(ctgid = str_replace(members, "_length.*", "")) %>%
  dplyr::select(-members)

a <- df_ganno_canthyd %>%
  dplyr::select(ctgid) %>%
  distinct() %>%
  left_join(df_clst_ctg_seprows, by = "ctgid") %>%
  mutate(ctgid = str_replace(repctg, "_length.*", "")) %>%
  distinct() %>%
  inner_join(df_taxa_krakenstd, by = "ctgid") %>%
  mutate(ctglen = str_replace(repctg, ".*length", "")) %>%
  mutate(ctglen = str_replace(ctglen, "^[_=]", "")) %>%
  mutate(ctglen = str_replace(ctglen, "_.*", "")) %>%
  mutate(ctglen = str_replace(ctglen, "NZ", "3000000")) %>%
  mutate(ctglen = as.integer(ctglen)) %>%
  group_by(Species) %>%
  summarise(across(ctglen, sum))


df_clst_ctg_heatmap <- df_clst_ctg %>%
  mutate(tplg = ifelse(str_detect(members, "circular"), "circular", "linear")) %>%
  mutate(repctg = ifelse(str_detect(repctg, "circular"), repctg, paste0(repctg, "_", tplg))) %>%
  mutate(repctg = ifelse(str_detect(repctg, "_linear"), str_replace(repctg, "_linear", ""), repctg)) %>%
  dplyr::select(-tplg) %>%
  separate_rows(members, sep = ",") %>%
  mutate(sampleid = str_replace(members, "_.*", "")) %>%
  dplyr::select(-members) %>%
  distinct() %>%
  mutate(exist = 1) %>%
  pivot_wider(names_from = repctg, values_from = exist) %>%
  column_to_rownames("sampleid") %>%
  mutate_all(~replace_na(., 0))




# Define UI for application that draws a histogram
side_p <- sidebarPanel(
    sliderInput("height", "Height:", min = 600, max = 5000, value = 600),
    sliderInput("fontsize_row", "Row Font Size:", min = 1, max = 20, value = 10),
    selectInput(
       inputId = "topology",
       label = "Topology",
       choices = c("circular", "linear"),
       selected = "circular",
       multiple = FALSE,
       selectize = TRUE,
       width = NULL,
       size = NULL
    ),
    uiOutput("ui_sel_repctg"),
    width = 2
)

main_p <- mainPanel(
    tabsetPanel(
        tabPanel("PAM", fluidRow(uiOutput("ui_hm_clst_ctg"))),
        tabPanel("Contig Members", fluidRow(reactableOutput("ctg_members"),
                                            h4("pprmeta"),
                                            reactableOutput('clsf_pprmeta'),
                                            h4("Gene annotation: CANT-HYD"),
                                            reactableOutput('ganno_canthyd'),
                                            h4("Gene annotation: NR"),
                                            reactableOutput('ganno_nr'),
                                            h4("Kraken standard DB"),
                                            reactableOutput('taxa_kstd'),
                                            h4("Kraken virus DB"),
                                            reactableOutput('taxa_kvir'),
                                            h4("MMseqs virus RefSeq DB"),
                                            reactableOutput('taxa_mvref'),
                                            h4("vConTACT2 reference"),
                                            reactableOutput('taxa_vc2ref'),
                                            h4("vConTACT2 all"),
                                            reactableOutput('taxa_vc2all'),
                                            )),
        tabPanel("CANT-HYD", fluidRow(reactableOutput("tbl_ganno_canthyd"),
                                      h4("All contig clusters"),
                                      reactableOutput("tbl_ctg_members")
                                      )),
        tabPanel("Host", fluidRow(reactableOutput("host_genome"),
                                  reactableOutput("host_genus")))
    ),
    width = 10
)


navp_data <- tabPanel("Data", sidebarLayout(side_p, main_p))
# navp_genes <- tabPanel("Analyses", sidebarLayout(side_p, main_p))
ui <- navbarPage("P0062: VirOil", navp_data,
                theme = shinytheme("united"))

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$ui_hm_clst_ctg <- renderUI ({
    plotlyOutput("hm_clst_ctg", height = input$height, width = 1400)
  })

  output$ui_sel_repctg <- renderUI({
    selectInput("sel_repctg", "Select RepContig", choices = unique(df_clst_ctg$repctg), selected = NULL, multiple = T, selectize = T)
  })

  pam_mx <- reactive({
    if (input$topology == "circular") {
      df <- df_clst_ctg_heatmap[, str_detect(colnames(df_clst_ctg_heatmap), input$topology)]
    } else {
      idx <- colSums(df_clst_ctg_heatmap)
      df <- df_clst_ctg_heatmap[, idx > 1]
    }

    df <- t(as.matrix(df))
  })

  clst_all <- reactive({
    df_clst_ctg %>%
      separate_rows(members, sep = ",") %>%
      mutate(ctgid = str_replace(members, "_length.*", "")) %>%
      mutate(repid = str_replace(repctg, "_length.*", "")) %>%
      left_join(df_taxa_vc2, by = c("repid" = "contig_id")) %>%
      left_join(df_virus_virify, by = c("repid" = "ctgid")) %>%
      left_join(df_virus_wtp, by = c("repid" = "contig_name"))
  })

  clst_sel <- reactive({
    clst_all() %>%
      dplyr::filter(repctg %in% input$sel_repctg)
    # clst_all() %>%
    #   dplyr::filter(ctgid %in% unique(df_ganno_canthyd$ctgid))
  })

  output$hm_clst_ctg <- renderPlotly({
    heatmaply(pam_mx(),
              colors = c("#FFFFFF", "red"),
              # width = 1400,
              height = input$height,
              cluster_rows = T,
              cluster_cols = T,
              grid_size = 1,
              grid_color = "grey",
              fontsize_row = input$fontsize_row,
              fontsize_col = 8,
              )})

  output$ctg_members <- renderReactable({
    clst_sel() %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$clsf_pprmeta <- renderReactable({
    df_clsf_pprmeta %>%
      dplyr::filter(repid %in% clst_sel()$ctgid) %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$taxa_kstd <- renderReactable({
    df_taxa_krakenstd %>%
      dplyr::filter(ctgid %in% clst_sel()$ctgid) %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$taxa_kvir <- renderReactable({
    df_taxa_krakenvir %>%
      dplyr::filter(ctgid %in% clst_sel()$ctgid) %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$taxa_mvref <- renderReactable({
    df_taxa_mvref %>%
      dplyr::filter(ctgid %in% clst_sel()$ctgid) %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$taxa_vc2ref <- renderReactable({
    df_taxa_vc2ref %>%
      dplyr::filter(VC %in% clst_sel()$VC) %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$taxa_vc2all <- renderReactable({
    df_taxa_vc2all %>%
      dplyr::filter(VC %in% clst_sel()$VC) %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$ganno_canthyd <- renderReactable({
    df_ganno_canthyd %>%
      dplyr::filter(ctgid %in% clst_sel()$ctgid) %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$ganno_nr <- renderReactable({
    df_ganno_nr %>%
      dplyr::filter(ctgid %in% clst_sel()$ctgid) %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  # CANT-HYD
  output$tbl_ganno_canthyd <- renderReactable({
    df_ganno_canthyd %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$tbl_ctg_members <- renderReactable({
    df_ganno_canthyd %>%
      left_join(clst_all(), by = "ctgid") %>%
      dplyr::select(c(virify_quality, repctg)) %>%
      distinct() %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$host_genome <- renderReactable({
    df_host_genome %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })

  output$host_genus <- renderReactable({
    df_host_genus %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })


}

# Run the application
options(shiny.port = 8062)
shinyApp(ui = ui, server = server)
library(shinydashboard)
library(shinythemes)
library(here)
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
prj_path <- normalizePath(".")
here(prj_path, "R/import_data.R")
source(here(prj_path, "R/import_data.R"))
source(here(prj_path, "R/utils.R"))
source(here(prj_path, "R/annotation_prot.R"))
load(here(prj_path, "data/vpf.RData"))

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
  uiOutput("sample_metadata"),
  uiOutput("contig_len"),
  radioButtons("logic", "Filter logic",
               choices = c("AND", "OR"),
               selected = "AND"),
  uiOutput("vs2cat"),
  uiOutput("checkv_quality"),
  uiOutput("ctg_feature_col"),
  uiOutput("prot_feature_col"),
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
    selectInput("ctg_feature_col", "Contig annotations",
                colnames(ctgmeta),
                selected = "Family")
  })

  output$prot_feature_col <- renderUI({
    selectInput("prot_feature_col", "Protein annotations",
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
