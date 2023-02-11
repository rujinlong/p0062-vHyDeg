#' vhydeg2 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_vhydeg2_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarLayout(
      # Side panel
      sidebarPanel(
        sliderInput(ns("bitscore"), "Minimum bitscore", min = 150, max = 650, value = 300, step=10),
        uiOutput(ns("hydg")),
        width = 3
      ),

      # Main panel
      mainPanel(
        tabsetPanel(
          tabPanel("IMG/VR",
                   shinycssloaders::withSpinner(reactableOutput(ns("tbl_imgvr")), type=8)
          ),
          tabPanel("PHROGs",
                   shinycssloaders::withSpinner(reactableOutput(ns("tbl_phrog")), type=8)
          ),
        ),
        width = 9
      )
    )
  )
}

#' vhydeg2 Server Functions
#'
#' @noRd
mod_vhydeg2_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    data("vhydeg_db")
    dataobj <- vhydeg_db
    # ========== UI =============
    output$hydg <- renderUI({
      selectizeInput(ns('hydg'), label= "Filter by gene name", dataobj[["gene_names"]], selected = NULL, multiple=TRUE, options = list(placeholder = 'select multiple genes'))
    })

    # ========== load data =============
    df_imgvr <- reactive({
      df_show <- dataobj[["df_imgvr"]]
      if (length(input$hydg)==0) {
        df_show <- df_show %>%
          dplyr::filter(bitscore >= input$bitscore)
      } else {
        df_show <- df_show %>%
          dplyr::filter(gene_name %in% input$hydg) %>%
          dplyr::filter(bitscore >= input$bitscore)
      }
    }) %>% bindCache(input$hydg, input$bitscore)

    # ----------- show table ----------
    output$tbl_imgvr <- renderReactable({
      df_imgvr_parent <- df_imgvr() %>%
        dplyr::select(-c(ipr_db_name, ipr_db_acc, ipr_db_desc, ipr_acc, ipr_desc, ipr_evalue)) %>%
        dplyr::distinct()

      df_imgvr_children <- df_imgvr() %>%
        dplyr::select(c(gene_id, ipr_db_name, ipr_db_acc, ipr_db_desc, ipr_acc, ipr_desc, ipr_evalue)) %>%
        dplyr::distinct()

      reactable(df_imgvr_parent,
                searchable = TRUE,
                defaultPageSize = 25,
                sortable = TRUE,
                filterable = TRUE,
                highlight = TRUE,
                wrap = FALSE,
                resizable = TRUE,
                defaultColDef = colDef(minWidth = 120),
                columns = list(
                  gene_id = colDef(minWidth = 300, cell = function(value, index) {
                    url <- sprintf('https://img.jgi.doe.gov/cgi-bin/vr/main.cgi?section=MetaGeneDetail&page=geneDetail&data_type=assembled&gene_oid=%s&taxon_oid=%s', str_replace_all(value, "^.*\\|", ""), str_replace(str_replace(value, "^[^|]+[|]", ""), "[|].*", ""))
                    htmltools::tags$a(href = url, target = "_blank", as.character(value))
                  }),

                  genome_id = colDef(minWidth = 200, cell = function(value, index) {
                    url <- sprintf('https://img.jgi.doe.gov/cgi-bin/vr/main.cgi?section=ViralBrowse&page=uviginfo&uvig_id=%s', value)
                    htmltools::tags$a(href = url, target = "_blank", as.character(value))
                  }),

                  gene_name = colDef(minWidth = 130),
                  gene_description = colDef(minWidth = 250),
                  Ecosystem = colDef(minWidth = 200),
                  vOTU = colDef(minWidth = 150),
                  Topology = colDef(minWidth = 180),
                  geNomad_score = colDef(minWidth = 150),
                  MIUVIG_quality = colDef(minWidth = 150),
                  Taxa_classification = colDef(minWidth = 250),
                  Host_taxonomy_prediction = colDef(minWidth = 250),
                  Sequence_origin = colDef(minWidth = 150),

                  Enzyme = colDef(minWidth = 150),
                  Broad_Enzymatic_Group = colDef(minWidth = 150),
                  Closest_false_positive = colDef(minWidth = 150),

                  trusted_idx = colDef(
                    style = function(value) {
                      bar_style(width = value, fill = "#09A472", color = "#cfcfcf")
                    },
                    align = "left",
                    format = colFormat(digits = 2)
                  ),

                  noise_idx = colDef(
                    style = function(value) {
                      bar_style(width = value, fill = "#829C08", color = "#cfcfcf")
                    },
                    align = "left",
                    format = colFormat(digits = 2)
                  ),

                  phrog_id = colDef(cell = function(value, index) {
                    url <- sprintf('https://phrogs.lmge.uca.fr/cgi-bin/script_mega_2018.py?mega=%s', str_replace_all(value, "^phrog_", ""))
                    htmltools::tags$a(href = url, target = "_blank", as.character(value))
                  })
                ),

                details = function(index) {
                  children = df_imgvr_children %>%
                    dplyr::filter(gene_id == df_imgvr_parent$gene_id[index]) %>%
                    dplyr::select(-gene_id)
                  tbl <- reactable(children, resizable = T, wrap = F, highlight = T, sortable = T, filterable = T)
                  htmltools::div(style = "padding: 1rem", tbl)
                }
      )
    })

    output$tbl_phrog <- renderReactable({
      dataobj[["df_phrog"]] %>%
        reactable(searchable = TRUE,
                  defaultPageSize = 50,
                  sortable = TRUE,
                  filterable = TRUE,
                  highlight = TRUE,
                  wrap = TRUE,
                  resizable = TRUE,
                  defaultColDef = colDef(minWidth = 50),
                  columns = list(
                    phrog_id = colDef(cell = function(value, index) {
                      url <- sprintf('https://phrogs.lmge.uca.fr/cgi-bin/script_mega_2018.py?mega=%s', str_replace_all(value, "^phrog_", ""))
                      htmltools::tags$a(href = url, target = "_blank", as.character(value))
                    }),
                    gene_name = colDef(maxWidth = 150),
                    bitscore_vmax = colDef(maxWidth = 150)
                  ))
    })
  })
}

## To be copied in the UI
# mod_vhydeg2_ui("vhydeg2_1")

## To be copied in the server
# mod_vhydeg2_server("vhydeg2_1")
