#' tutorial UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_tutorial_ui <- function(id){
  ns <- NS(id)
  tagList(
    fixedPage(
      includeHTML("workflow/91-tutorial.html")
  ))
}

#' tutorial Server Functions
#'
#' @noRd
mod_tutorial_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_tutorial_ui("tutorial_1")

## To be copied in the server
# mod_tutorial_server("tutorial_1")
