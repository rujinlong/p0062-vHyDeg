#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # mod_vhydeg_server("vhydeg_1")
  mod_vhydeg2_server("vhydeg2_1")
}
