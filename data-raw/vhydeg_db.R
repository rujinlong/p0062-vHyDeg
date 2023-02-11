## code to prepare `vhydeg_db` dataset goes here

vhydeg_db <- readRDS("workflow/data/31-vhydeg_db/vhydeg_db.rds")
usethis::use_data(vhydeg_db, overwrite = TRUE)
