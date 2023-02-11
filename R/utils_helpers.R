#' Read amg_db file
#'
#' @param fin Path to amg_db file
#'
#' @return data.frame
#' @export
#'
#' @examples
#' @noRd
read_amgdb <- function(fin) {
  amgs <- data.table::fread(fin) %>%
    dplyr::filter(KO!="") %>%
    dplyr::rename(ko_id=KO, gene_func=gene) %>%
    dplyr::select(c(ko_id, gene_func)) %>%
    dplyr::mutate(amg=1)
  return(amgs)
}


#' Read the annotation file of the proteome
#'
#'
#' @param fp_anno_prot Path to the annotation file of the proteome
#' @param fp_amg Path to the annotation file of the AMGs
#'
#' @return
#' @export
#'
#' @examples
#' @noRd
read_anno_vprot <- function(fp_anno_prot, fp_amg) {
  amgs <- read_amgdb(fp_amg)
  anno_prot <- data.table::fread(fp_anno_prot) %>%
    dplyr::left_join(amgs, by = "ko_id")

  return(anno_prot)
}


#' Create a style for a bar in the vhydeg annotation table
#'
#' @param width with of the bar
#' @param fill color of the bar
#' @param height height of the bar
#' @param align alignment of the bar, either "left" or "right"
#' @param color color of the text, if NULL, the color of the text will be the same as the color of the bar
#'
#' @return
#' @export
#'
#' @examples
bar_style <- function(width = 1, fill = "#e6e6e6", height = "70%",
                      align = c("left", "right"), color = NULL) {
  align <- match.arg(align)
  if (align == "left") {
    position <- paste0(width * 100, "%")
    image <- sprintf("linear-gradient(90deg, %1$s %2$s, transparent %2$s)", fill, position)
  } else {
    position <- paste0(100 - width * 100, "%")
    image <- sprintf("linear-gradient(90deg, transparent %1$s, %2$s %1$s)", position, fill)
  }
  list(
    backgroundImage = image,
    backgroundSize = paste("100%", height),
    backgroundRepeat = "no-repeat",
    backgroundPosition = "center",
    color = color
  )
}


#' Filter tblastx results to only include the AMG hit and the corresponding flanking regions
#'
#' @param fin Path to the tblastx file
#' @param tarvir_id ID of the target virus (without the AMG hit)
#' @param hit_start Start position of AMG hit on the query virus
#' @param hit_end End position of the AMG hit on the query virus
#' @param flank_len Length of the flanking region
#' @param remove_seq_version Remove the version number of the sequence ID. Default is TRUE.
#' @param ... Other arguments
#'
#' @return data.frame
#' @export
#'
#' @examples
filter_tblastx <- function(fin, tarvir_id, hit_start, hit_end, flank_len, remove_seq_version=T, ...) {
  df_tblastx <- vpfkit::read_tblastx(fin, ...) %>%
    # order by start
    dplyr::arrange(.data$start)

  if (remove_seq_version) {
    df_tblastx$seq_id2 <- gsub("\\.[0-9]+$", "", df_tblastx$seq_id2)
    tarvir_id <- gsub("\\.[0-9]+$", "", tarvir_id)
  }

  tarbac_id <- df_tblastx %>%
    dplyr::filter(.data$seq_id2 != tarvir_id) %>%
    dplyr::filter(.data$start > hit_start & .data$end < hit_end) %>%
    dplyr::arrange(desc(.data$pident)) %>%
    head(1) %>%
    dplyr::pull("seq_id2")

  tar_ids <- c(tarbac_id, tarvir_id)
  df_tblastx_sel <- df_tblastx %>%
    dplyr::filter(.data$seq_id2 %in% tar_ids) %>%
    dplyr::filter(.data$start > hit_start - flank_len & .data$end < hit_end + flank_len)

  return(df_tblastx_sel)
}


#' Read bins of the metagenome
#'
#' @param fbin Path to the bin file
#' @param fbin_classify Path to the bin classification file
#'
#' @return data.frame
#' @export
#'
#' @examples
#' @noRd
read_bins <- function(fbin, fbin_classify) {
  bin_classify <- fread(fbin_classify) %>%
    column_to_rownames("binname") %>%
    setnames(colnames(.), paste0("bin_", colnames(.))) %>%
    rownames_to_column("binname")
  bins <- fread(fbin, header = FALSE, col.names = c("binname", "Contig")) %>%
    left_join(bin_classify, by = "binname")
  return(bins)
}


#' Parse pprmeta annotation results
#'
#' @param fin_pprmeta Path to the pprmeta file
#' @param fin_mapping Path to the mapping file
#'
#' @return data.frame
#' @export
#'
#' @examples
#' @noRd
parse_pprmeta <- function(fin_pprmeta, fin_mapping) {
  df_pprmeta <- fread(fin_pprmeta)
  df <- fread(fin_mapping, col.names = c("repid", "Header"), header = F) %>%
    inner_join(df_pprmeta, by = "Header") %>%
    dplyr::select(-Header)
  return(df)
}

#' Parse ganno annotation results
#'
#' @param fin Path to the annotation file
#' @param fin_mapping Path to the mapping file
#'
#' @return data.frame
#' @export
#'
#' @examples
#' @noRd
parse_ganno <- function(fin, fin_mapping) {
  df_ganno <- fread(fin)
  df <- fread(fin_mapping, col.names = c("query_id", "gene_id"), header = F) %>%
    inner_join(df_ganno, by = "query_id") %>%
    dplyr::select(-query_id) %>%
    mutate(ctgid = str_replace(gene_id, "_[0-9]*$", ""))
  return(df)
}




#' Read bakta (UniProt) annotation
#'
#' @param fpath_bakta Path to bakta annotation
#'
#' @return
#' @export
#'
#' @examples
read_bakta <- function(fpath_bakta) {
  df_bakta <- read.csv(fpath_bakta, sep = "\t") %>%
    dplyr::mutate(gene_id = ID, uniprot_gene=Gene, uniprot_product=Product) %>%
    dplyr::select(c(gene_id, uniprot_gene, uniprot_product))

  return(df_bakta)
}


#' Read InterProScan annotation
#'
#' @param fpath_ipr Path to InterProScan annotation
#'
#' @return
#' @export
#'
#' @examples
read_ipr <- function(fpath_ipr) {
  df_ipr <- read.csv(fpath_ipr, header = FALSE, sep = "\t")
  if (ncol(df_ipr) == 15) {
    colnames(df_ipr) <- c("gene_id", "md5", "length", "ipr_db_name", "ipr_db_acc", "ipr_db_desc", "start", "end", "ipr_evalue", "status", "date", "ipr_acc", "ipr_desc", "go_anno", "path_anno")
  } else if (ncol(df_ipr) == 13) {
    colnames(df_ipr) <- c("gene_id", "md5", "length", "ipr_db_name", "ipr_db_acc", "ipr_db_desc", "start", "end", "ipr_evalue", "status", "date", "ipr_acc", "ipr_desc")
  }
  df_ipr <- df_ipr %>%
    dplyr::select(c(gene_id, ipr_db_name, ipr_db_acc, ipr_db_desc, ipr_acc, ipr_desc, ipr_evalue)) %>%
    dplyr::distinct()

  return(df_ipr)
}


#' Read Phame gene2family annotation
#'
#' @param gene2family_raw gene2family raw data
#'
#' @return
#' @export
#'
#' @examples
get_gene2family <- function(gene2family_raw) {
  # split df_imgvr_families_sel into separate data frames for each gene_name, and store to a list named "families"
  families <- split(gene2family_raw, gene2family_raw$gene_name)

  # count the frequency of each family in each gene_name
  families_top1 <- lapply(families, function(x) {
    x %>%
      dplyr::select(family) %>%
      group_by(family) %>%
      summarise(count = n()) %>%
      mutate(gene_name = x$gene_name[1]) %>%
      arrange(desc(count)) %>%
      # select the top 1 family
      slice(1) %>%
      dplyr::select(family)
  })

  # combine the data frames in the list "families_top1" into a single data frame
  gene2family <- do.call(rbind, families_top1) %>%
    rownames_to_column("gene_name") %>%
    mutate(gene_family = str_c(gene_name, "__", family))

  return(gene2family)
}
