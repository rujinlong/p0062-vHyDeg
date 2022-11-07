-   [Introduction](#introduction)
    -   [Tasks](#tasks)
        -   [Protein annotation](#protein-annotation)
        -   [Tree visualization](#tree-visualization)
    -   [Files written](#files-written)
        -   [References](#references)

Updated at 2022-11-07 15:58:25 UTC.

# Introduction

We have three phylogenetic trees for three protein families,

1.  fam0023_tree.nwk
2.  fam0026_tree.nwk
3.  fam0253_tree.nwk

There is a table contains annotation of all the proteins in these trees.
Annotations include taxonomy, environment etc. Expected output of this
task is three phylogenetic trees, with taxonomy and environment
annotations.

Trees can be plotted using `ggtree`. Annotations can be added using
`ggtreeExtra`. Tutorial can be found
[here](https://yulab-smu.top/treedata-book/chapter10.html). One example
is [this](https://www.nature.com/articles/s41586-020-2007-4/figures/2).

## Tasks

TODO:

-   [ ] Change layout of the tree if necessary

<div class="cell">

``` r
library(here)
library(conflicted)
library(tidyverse)
```

<div class="cell-output-stderr">

    ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ✔ ggplot2 3.3.6      ✔ purrr   0.3.5 
    ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ✔ readr   2.1.3      ✔ forcats 0.5.2 

</div>

``` r
library(treedataverse)
```

<div class="cell-output-stderr">

    ── Attaching packages ─────────────────────────────────── treedataverse 0.0.1 ──
    ✔ ape         5.6.2      ✔ ggtree      3.4.4 
    ✔ tidytree    0.4.1      ✔ ggtreeExtra 1.6.1 
    ✔ treeio      1.20.2     

</div>

``` r
library(data.table)
library(tidytree)
library(ape)
library(ggtree)

dir_prj <- normalizePath("..")
wd <- here(dir_prj, "data/fig1_HYDG_tree")

# function to read the ".nwk" tree file, plot using ggtree, and save to pdf file
create_tree <- function(ftree, df_anno_imgvr, df_anno_nr) {
  tree_tbl <- read.tree(ftree) %>%
    as_tibble() %>%
    mutate(UVIG = ifelse(str_detect(label, "^IMGVR"), label, "")) %>%
    mutate(UVIG = str_replace_all(UVIG, "\\|.*", ""))

  tree_IMGVR <- tree_tbl %>% 
    dplyr::filter(UVIG != "") %>%
    left_join(df_anno_imgvr[, c("UVIG", "Kingdom", "Phylum", "Class", "eco3")], by = "UVIG") %>%
    mutate(Kingdom = "Viruses")

  tree_nr <- tree_tbl %>%
    dplyr::filter(UVIG == "") %>%
    left_join(df_anno_nr[, c("protid", "Kingdom", "Phylum", "Class")], by = c("label" = "protid")) %>%
    mutate(eco3 = "Unknown")

  tree <- rbind(tree_IMGVR, tree_nr) %>%
    as.treedata()

  p <- ggtree(tree, layout="equal_angle", color="grey") +
    geom_tiplab2(aes(label=Kingdom, color = Kingdom), size=0, color="black", offset=0) +
    geom_tippoint(aes(color=Kingdom), size=0.2, alpha=0.7) 
  # output file name is basename of ftree without extension
  ggsave(filename=path_target(paste0(basename(ftree), ".pdf")), plot=p, width = 20, height = 13)
}
```

</div>

### Protein annotation

Proteins from three sources,

1.  IMGVR, which include query viral proteins from IMGVR database.
2.  Custom data, which include query viral proteins form other sources.
3.  NR, which include proteins of bacteria, archaea, eukaryota and
    viruses in NCBI NR database that is homologous of viral proteins in
    1 and 2.

Annotations of NR are stored in `df_anno_nr`, and saved to file
`annotation_nr.tsv`.

<div class="cell">

``` r
# taxid2lineage <- fread(here(wd, "nf_tmp_taxaid2lineage.tsv.gz"), header=F) %>%
#   setnames(colnames(.), c("taxid", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
# 
# name2taxid <- fread(here(wd, "nf_tmp_name2taxid.tsv.gz"), header=F) %>%
#   setnames(colnames(.), c("taxa", "taxid")) %>%
#   # sort by taxid in ascending order
#   arrange(taxid) %>%
#   # drop duplicates by taxa
#   distinct(taxa, .keep_all = TRUE)
# 
# df_anno_nr <- fread(here(wd, "targets_taxa.tsv"), header = F) %>% 
#   setnames(colnames(.), c("protid", "taxa")) %>%
#   left_join(name2taxid, by = "taxa") %>%
#   left_join(taxid2lineage, by = "taxid") %>%
#   distinct(protid, .keep_all = TRUE) %>%
#   # replace NA with "Unknown" in Kingdom, Phylum, Class, Order, Family, Genus, Species
#   mutate(Kingdom = ifelse(is.na(Kingdom), "Unknown", Kingdom),
#          Phylum = ifelse(is.na(Phylum), "Unknown", Phylum),
#          Class = ifelse(is.na(Class), "Unknown", Class),
#          Order = ifelse(is.na(Order), "Unknown", Order),
#          Family = ifelse(is.na(Family), "Unknown", Family),
#          Genus = ifelse(is.na(Genus), "Unknown", Genus),
#          Species = ifelse(is.na(Species), "Unknown", Species))
# 
# # write protein taxonomy annotation to file
# fwrite(df_anno_nr, path_target("annotation_nr.tsv"), sep = "\t", quote = F, row.names = F)
df_anno_nr <- fread(here(wd, "annotation_nr.tsv.gz"))
```

</div>

Annotations of IMGVR and others are stored in `df_anno_imgvr`, and saved
to file `annotation_imgvr.tsv`.

<div class="cell">

``` r
# df_anno_imgvr <- fread(here(wd, "HYDGfam_contig_info2.tsv")) %>%
#   setnames(colnames(.), make.names(colnames(.))) %>%
#   mutate(Host.taxonomy.prediction = ifelse(Host.prediction.method=="Isolate taxonomy", Host.taxonomy.prediction, "")) %>%
#   dplyr::select(c(UVIG, Ecosystem.classification, Topology, Taxonomic.classification, Host.taxonomy.prediction)) %>%
#   # split Taxonomic.classification to Domain, Kingdom, Phylum, Class, Order, Family, Genus based on ";"
#   tidyr::separate(Taxonomic.classification, c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>%
#   # split Host.taxonomytidyr::.prediction to host_Kingdom, host_Phylum, host_Class, host_Order, host_Family, host_Genus, host_Species based on ";"
#   tidyr::separate(Host.taxonomy.prediction, c("host_Kingdom", "host_Phylum", "host_Class", "host_Order", "host_Family", "host_Genus", "host_Species"), sep = ";") %>%
#   # split Ecosystem.classification to eco1, eco2, eco3, eco4 based on ";"
#   tidyr::separate(Ecosystem.classification, c("eco1", "eco2", "eco3", "eco4"), sep = ";") %>%
#   # replace all "^.*__" with ""
#   mutate_all(funs(str_replace_all(., "^.*__", "")))
# 
# # write protein taxonomy annotation to file
# fwrite(df_anno_imgvr, path_target("annotation_imgvr.tsv"), sep = "\t", quote = F, row.names = F)
df_anno_imgvr <- fread(here(wd, "annotation_imgvr.tsv.gz"))
```

</div>

### Tree visualization

Phylogenetic tree of each protein family was inferred using the maximum
likelihood method with the fasttree program (Price, Dehal, and Arkin
2010). The best-scoring ML tree was selected from 100 bootstrap
replicates. Taxonomy annotation of each protein was obtained from the
hit results and normalized using Taxonkit software (Shen and Ren 2021).
The tree was visualized using ggtree software (Xu et al. 2022).

TODO:

-   [ ] Plot tree for fam0253

-   [ ] Collapse tree to phylum level. That is, if nodes in the tree
    have same phylum, collapse them into one node. Because each tree
    contains thousands of nodes, we hope using this approach, we can
    have fewer nodes and make the tree easier to visualize.

-   [ ] Add Phylum name, environment source as a circular heatmap in the
    outside of the tree.

<div class="cell">

``` r
# create_tree(here(wd, "top_alkb_tree.nwk"), df_anno_imgvr, df_anno_nr)
```

</div>

Files were saved to
/Users/allen/github/rujinlong/phydgene/workflow/data/10_fig1-HYDG-tree.

TODO:

-   [ ] If the plot looks good, plot for all three trees

<div class="cell">

``` r
ftrees <- list.files(here(wd), pattern = "fam.*_tree.nwk", full.names = T) %>%
  # set names using famid
  setNames(stringr::str_extract(., "fam\\d+"))

# create tree visualization for each protein family using function `create_tree`
# very time-consuming
# lapply(ftrees, create_tree, df_anno_nr = df_anno_nr, df_anno_imgvr = df_anno_imgvr)
```

</div>

## Files written

These files have been written to the target directory,
`data/10_fig1-HYDG-tree`:

<div class="cell">

``` r
projthis::proj_dir_info(path_target())
```

<div class="cell-output-stdout">

    # A tibble: 11 × 4
       path                            type         size modification_time  
       <fs::path>                      <fct> <fs::bytes> <dttm>             
     1 IMGVR_HYDG_contigs.tsv          file       21.83K 2022-11-06 17:38:33
     2 annotation_imgvr.tsv            file      430.95K 2022-11-07 15:36:43
     3 annotation_nr.tsv               file       10.66M 2022-11-07 15:36:35
     4 fam0023_tree.tsv                file       10.45M 2022-11-06 15:57:37
     5 fam0023_tree2.pdf               file        4.28M 2022-11-06 15:29:36
     6 fam0026_tree2.pdf               file        1.58M 2022-11-06 15:30:12
     7 fam0253_tree2.pdf               file      767.79K 2022-11-06 15:30:25
     8 protein_taxonomy.tsv            file       10.66M 2022-11-06 22:22:50
     9 protein_taxonomy_no_kingdom.tsv file       10.93K 2022-11-06 15:25:49
    10 test_tree.pdf                   file      856.65K 2022-11-06 22:13:26
    11 top_alkb_tree.nwk.pdf           file       45.62K 2022-11-07 15:53:42

</div>

</div>

### References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-2010-FastTree_Price" class="csl-entry">

Price, Morgan N., Paramvir S. Dehal, and Adam P. Arkin. 2010. “FastTree
2 Approximately Maximum-Likelihood Trees for Large Alignments.” *PLOS
ONE* 5 (3): e9490. <https://doi.org/10.1371/journal.pone.0009490>.

</div>

<div id="ref-2021-TaxonKit_Shen" class="csl-entry">

Shen, Wei, and Hong Ren. 2021. “TaxonKit: A Practical and Efficient NCBI
Taxonomy Toolkit.” *Journal of Genetics and Genomics*, April.
<https://doi.org/10.1016/j.jgg.2021.03.006>.

</div>

<div id="ref-2022-Ggtree_Xu" class="csl-entry">

Xu, Shuangbin, Lin Li, Xiao Luo, Meijun Chen, Wenli Tang, Li Zhan, Zehan
Dai, Tommy T. Lam, Yi Guan, and Guangchuang Yu. 2022. “Ggtree: A
Serialized Data Object for Visualization of a Phylogenetic Tree and
Annotation Data.” *iMeta*, September. <https://doi.org/10.1002/imt2.56>.

</div>

</div>
