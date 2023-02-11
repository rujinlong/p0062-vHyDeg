attachment::att_amend_desc()
golem::add_module(name = "vhydeg", with_test = TRUE)
golem::add_module(name = "vhydeg2", with_test = TRUE)
golem::add_fct("helpers", with_test = TRUE)
golem::add_utils("helpers", with_test = TRUE)
golem::add_js_file("script")
golem::add_js_handler("handlers")
golem::add_css_file("custom")
golem::add_sass_file("custom")
usethis::use_data_raw(name = "my_dataset", open = F)
usethis::use_test("app")

# change vignette name
usethis::use_vignette("vhydeg")
devtools::build_vignettes()
usethis::use_github(private = T)
usethis::use_coverage("codecov")

# finishing setup github
rstudioapi::navigateToFile("dev/03_deploy.R")
