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
fp_virus_virify <- here(dpath, "metaall/virify/08-final/virus_virify.tsv")
fp_virus_wtp <- here(dpath, "metaall/wtp/final_report.utf8.csv")
fp_clsf_pprmeta <- here(dpath, "metaall/virify/01-viruses/pprmeta/repctg_pprmeta.csv")
fp_map_virify <- here(dpath, "metaall/virify/virify_mapping.tsv")
fp_map_gene <- here(dpath, "metaall/ganno/qry_cluster.tsv")
fp_ganno_canthyd <- here(dpath, "metaall/ganno/canthyd.tsv")

df_clst_ctg <- fread(fp_clst_ctg, header = F, col.names = c("repctg", "members"))
df_taxa_krakenstd <- fread(fp_taxa_krakenstd)
df_taxa_krakenvir <- fread(fp_taxa_krakenvir)
df_taxa_mvref <- fread(fp_taxa_mvref)
df_taxa_vc2 <- fread(fp_taxa_vc2)
df_taxa_vc2ref <- fread(fp_taxa_vc2ref)
df_virus_virify <- fread(fp_virus_virify)
df_virus_wtp <- fread(fp_virus_wtp) %>%
  dplyr::select(-V1)
df_clsf_pprmeta <- parse_pprmeta(fp_clsf_pprmeta, fp_map_virify)
df_ganno_canthyd <- parse_ganno(fp_ganno_canthyd, fp_map_gene)


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
    sliderInput("width",
                "Width:",
                min = 500,
                max = 5000,
                value = 1000),
    sliderInput("height", "Height:", min = 500,
                max = 5000,
                value = 1200),
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
                                            h4("Kraken standard DB"),
                                            reactableOutput('taxa_kstd'),
                                            h4("Kraken virus DB"),
                                            reactableOutput('taxa_kvir'),
                                            h4("MMseqs virus RefSeq DB"),
                                            reactableOutput('taxa_mvref'),
                                            h4("vConTACT2 reference"),
                                            reactableOutput('taxa_vc2ref')))
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
    plotlyOutput("hm_clst_ctg", height = input$height, width = input$width)
  })

  output$ui_sel_repctg <- renderUI({
    selectInput("sel_repctg", "Select RepContig", choices = unique(df_clst_ctg$repctg), selected = NULL, multiple = T, selectize = T)
  })

  pam_mx <- reactive({
    if (input$topology == "circular") {
      df <- df_clst_ctg_heatmap[, str_detect(colnames(df_clst_ctg_heatmap), input$topology)]
    } else {
      df <- df_clst_ctg_heatmap
    }

    df <- t(as.matrix(df))
  })

  clst_sel <- reactive({
    df_clst_ctg %>%
      dplyr::filter(repctg %in% input$sel_repctg) %>%
      separate_rows(members, sep = ",") %>%
      mutate(ctgid = str_replace(members, "_length.*", "")) %>%
      mutate(repid = str_replace(repctg, "_length.*", "")) %>%
      left_join(df_taxa_vc2, by = c("repid" = "contig_id")) %>%
      left_join(df_virus_virify, by = c("repid" = "ctgid")) %>%
      left_join(df_virus_wtp, by = c("repid" = "contig_name"))
  })

  output$hm_clst_ctg <- renderPlotly({
    heatmaply(pam_mx(),
              colors = c("#FFFFFF", "red"),
              width = input$width,
              height = input$height,
              cluster_rows = T,
              cluster_cols = T,
              fontsize_col = 8,
              fontsize_row = 8,
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

  output$ganno_canthyd <- renderReactable({
    df_ganno_canthyd %>%
      dplyr::filter(ctgid %in% clst_sel()$repctg) %>%
      reactable(resizable = T, highlight = T, filterable = T, sortable = T, wrap = F)
  })


}

# Run the application
options(shiny.port = 8062)
shinyApp(ui = ui, server = server)
