library(here)
source(here("script/utils.R"))

dpath <- here("data/hms")
fp_amg <- here("data/raw/amg_database.tsv")
fp_anno_prot <- here(dpath, "annotations.tsv")

# --------- AMGs ------------
read_amgdb <- function(fp_amg_db) {
  amgs <- fread(fp_amg_db) %>% 
    filter(KO!="") %>% 
    dplyr::rename(ko_id=KO, gene_func=gene) %>% 
    dplyr::select(c(ko_id, gene_func)) %>% 
    mutate(amg=1)
  return(amgs)
}

# ------- prot annotations ---------
read_anno_vprot <- function(fp_anno_prot, fp_amg) {
  amgs <- read_amgdb(fp_amg)
  anno_prot <- fread(fp_anno_prot) %>% 
    left_join(amgs, by = "ko_id")
}

anno_vprot <- read_anno_vprot(fp_anno_prot, fp_amg)


x <- "Glycosyl transferases group 1 [PF00534.23]; Glycosyl transferases group 1 [PF13692.9]; Glycosyl transferases group 1 [PF13524.9]"

