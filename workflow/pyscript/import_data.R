library(here)
library(tidyverse)
library(data.table)
library(speedyseq)
library(readxl)
library(gggenomes)


source(here("script/utils.R"))
dpath <- here("data/hms")

# ------- bins (optional) --------
bins <- read_bins(here(dpath, "clusters.tsv"), 
                  here(dpath, "vambbins_RF_predictions.txt"))

# -------  contig annotations --------
checkv <- read_CheckV(here(dpath, "quality_summary_viruses.tsv"))
raw_taxa <- read_taxonomy(here(dpath, "taxonomy.tsv"))
catbat <- read_CATBAT(here(dpath, "out.CAT.contig2classification_official.txt"))
vcontact <- read_vContact2(here(dpath, "genome_by_genome_overview.csv"), "__")
ctganno <- list("checkv"=checkv,
                "catbat"=catbat,
                "vcontact"=vcontact,
                "taxonomy"=raw_taxa,
                "bins"=bins)

ctganno_merged <- ctganno$checkv %>% 
  left_join(ctganno$taxonomy, by = "Contig") %>% 
  left_join(ctganno$catbat, by = "contig_id") %>% 
  left_join(ctganno$vcontact, by = "contig_id") %>% 
  left_join(ctganno$bins, by = "Contig") %>% 
  mutate(rowid = Contig) %>% 
  column_to_rownames("rowid")

# -------- gene annotation -----------
ganno <- fread(here(dpath, "annotations.tsv")) %>% 
  mutate(auxiliary_score=as.factor(auxiliary_score))

# ------- sample metadata ---------
sample_metadata <- read_excel(here("data/raw/metadata.xlsx"))

# ---------- abundance ------------
fp_count <- here(dpath, "abundance_contigs_count.tsv")
fp_tpm <- here(dpath, "abundance_contigs_tpm.tsv")
df_count <- fread(fp_count)
df_tpm <- fread(fp_tpm)


# Phyloseq
phy_count <- df_count %>% 
  column_to_rownames("Contig") %>% 
  as.matrix() %>% 
  otu_table(taxa_are_rows = T)
phy_meta <- sample_metadata %>% 
  column_to_rownames("sampleid") %>% 
  sample_data()
phy_taxa <- ctganno_merged %>% 
  as.matrix() %>%
  tax_table()
pseq <- phyloseq(phy_count, phy_taxa, phy_meta)
save(ctganno, ctganno_merged, ganno, pseq, file = here("data/vpf.RData"))






# Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
# # read gbk
# fgbk <- here("data/hms/test.gbk")
# fgbk <- here("data/hms/genbank.gbk")
# a <- gggenomes::read_gbk(fgbk)
#   mutate(gene=name) %>% 
#   left_join(ganno, by="gene")
# 
# class(a)
