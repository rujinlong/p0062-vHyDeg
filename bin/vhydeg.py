#!/usr/bin/env python

import click
import pandas as pd

def read_canthyd_annotation(fin, skip_non_NdoB=True):
    df = pd.read_csv(fin, sep="\t")
    if skip_non_NdoB:
        df = df[df['dbid'] != 'non_NdoB_type']
        # df.reset_index(drop=True, inplace=True)
    df['genome_id'] = df['query_id'].str.replace(r"\|.*", "", regex=True)
    return df


def read_cutoffs(fin):
    df = pd.read_csv(fin, sep="\t")
    df.rename(columns={'Length': 'HMM_length'}, inplace=True)
    return df


def read_phrogs(fin_canthyd, fp_phrog_desc):
    # read phrog annotated by canthyd
    df = pd.read_csv(fin_canthyd, sep="\t")
    df["phrog_id"] = df[["query_id"]].replace(r"_[a-zA-Z].*", "", regex=True)
    df = df[["dbid", "bitscore", "evalue", "phrog_id"]].copy()
    df.rename(columns={"bitscore": "bitscore_vmax", "evalue": "phrog_evalue"}, inplace=True)
    # sort by bitscore
    df.sort_values(by=["bitscore_vmax"], ascending=False, inplace=True)
    # drop duplicates
    df.drop_duplicates(subset=["dbid"], keep="first", inplace=True)

    # read phrog description
    df_phrog_desc = pd.read_csv(fp_phrog_desc, sep="\t")
    df_phrog_desc["phrog_id"] = df_phrog_desc.apply(lambda x: "phrog_"+str(x["phrog"]), axis=1)
    df_phrog_desc.rename(columns={"annot": "phrog_annotation", "category": "phrog_category"}, inplace=True)
    df_phrog_desc = df_phrog_desc[["phrog_id", "phrog_annotation", "phrog_category"]].copy()

    df_phrog = df.merge(df_phrog_desc, on="phrog_id", how="left")
    return df_phrog


def read_genome_info(fin):
    df_genome_info = pd.read_csv(fin, sep="\t")
    df_genome_info.rename(columns={"UVIG": "genome_id", 
                                   "Length": "genome_length", 
                                   "Coordinates ('whole' if the UViG is the entire contig)": "Coordinates",
                                   "Ecosystem classification": "Ecosystem",
                                   "geNomad score": "geNomad_score",
                                   "Estimated completeness": "Completeness",
                                   "Estimated contamination": "Contamination",
                                   "MIUViG quality": "MIUViG_quality",
                                   "Taxonomic classification": "Taxa",
                                   "Taxonomic classification method": "Taxa_method",
                                   "Host taxonomy prediction": "Host_taxonomy",
                                   "Host prediction method": "Host_prediction_method",
                                   "Sequence origin (doi)": "Sequence_origin",
                                   "Gene content (total genes;cds;tRNA;geNomad marker)": "Gene_content"}, inplace=True)
    # remove Gene_content column
    df_genome_info.drop(columns=["Gene_content"], inplace=True)
    return df_genome_info


def merge_dfs(df_pvhydegs, df_cutoffs, df_vmax, df_genome_info):
    """
    bitscore > bitscore_vmax
    """
    df1 = df_pvhydegs.merge(df_cutoffs, on="dbid", how="left")
    df1['trusted_idx'] = df1['bitscore'] / df1['Trusted_Cutoff']
    df1['trusted_idx'] = df1['trusted_idx'].round(2)
    df1['noise_idx'] = df1['bitscore'] / df1['Noise_Cutoff']
    df1['noise_idx'] = df1['noise_idx'].round(2)

    df2 = df1.merge(df_vmax, on="dbid", how="left")
    df2.rename(columns={"dbid": "target_gene"}, inplace=True)
    df2 = df2.merge(df_genome_info, on="genome_id", how="left")
    df2 = df2[df2['bitscore'] > df2['bitscore_vmax']].copy()
    # rename columns 
    df2.rename(columns={"query_id": "gene_id", "target_gene": "gene_name", "hmm_description": "gene_description"}, inplace=True)
    return df2[['gene_id', 'genome_id', 'gene_name', 'gene_description', 'evalue', 'bitscore', 'trusted_idx', 'noise_idx', 'Eenzyme', 'Broad_Enzymatic_Group', 'Substrate', 'Hydrocarbon_Group', 'Respiration', 'Trusted_Cutoff', 'Noise_Cutoff', 'Seed_sequences', 'HMM_length', 'Closest_false_positive', 'Closely_related_genes_above_trusted_cutoff',  'bitscore_vmax', 'phrog_evalue', 'phrog_id', 'phrog_annotation', 'phrog_category', 'Taxon_oid', 'Scaffold_oid', 'Coordinates', 'Ecosystem', 'vOTU', 'genome_length', 'Topology', 'geNomad_score', 'Confidence', 'Completeness', 'Contamination', 'MIUViG_quality', 'Taxa', 'Taxa_method', 'Host_taxonomy', 'Host_prediction_method', 'Sequence_origin']]


# def filter_NRqry_proteins(df):
#     """
#     For each gene, select the top1 to search against NR database. Filter by,
#     1. vHYDEG_quality == "HQ"
#     2. Contamination == 0
#     
#     Update: manually select based on completeness of the contig
#     """
#     df2 = df[(df['vHYDEG_quality'] == "HQ") & (df['Contamination'] == 0) & (df['Topology'] != 'Provirus') & (df['Coordinates'] == 'whole')].copy()
# 
#     # sort by bitscore in descending order
#     df2.sort_values(by="bitscore", ascending=False, inplace=True)
# 
#     # remove duplicates based on target_gene
#     df2.drop_duplicates(subset="target_gene", inplace=True)
#     df2["NR_qry"] = "yes"
# 
#     df3 = df.merge(df2[["query_id", "NR_qry"]], on="query_id", how="left")
#     df3["NR_qry"] = df3["NR_qry"].fillna("no")
# 
#     return df3


@click.command()
@click.option("--fin_canthyd_anno", '-a', help="CANT-HYD annotation file")
@click.option("--fin_cutoffs", '-c', help="CANT-HYD cutoffs file")
@click.option("--fin_phrog", '-p', help="CANT-HYD annotation of PHROGs clustered protein families")
@click.option("--fin_phrog_desc", '-d', help="phrog description file")
@click.option("--fin_genome_info", '-g', help="genome info file")
@click.option("--genome_minlen", '-l', type=int, default=5000, help="minimum genome length")
@click.option("--idx", '-i', default="noise_idx", help="index to filter by [trusted_idx, noise_idx]")
@click.option("--idx_threshold", '-t', type=float, default=0.8, help="index threshold 0-1")
@click.option("--fout", '-o', help="output file name")
def main(fin_canthyd_anno, fin_cutoffs, fin_phrog, fin_phrog_desc, fin_genome_info, genome_minlen, idx, idx_threshold, fout):
    """
    Only keep rows that have,
    1. idx > idx_threshold
    2. bitscore > bitscore_vmax
    """
    df_pvhydegs = read_canthyd_annotation(fin_canthyd_anno, skip_non_NdoB=True)
    df_cutoffs = read_cutoffs(fin_cutoffs)
    df_vmax= read_phrogs(fin_phrog, fin_phrog_desc)
    df_genome_info = read_genome_info(fin_genome_info)

    df = merge_dfs(df_pvhydegs, df_cutoffs, df_vmax, df_genome_info)
    df["vHYDEG_quality"] = df.apply(lambda x: "HQ" if x[idx] > float(idx_threshold) else "MQ", axis=1)
    df = df[df['genome_length'] > genome_minlen].copy()
    # df2 = filter_NRqry_proteins(df)
    
    df2.to_csv(fout, sep="\t", index=False)
    print("file saved to: {}".format(fout))

if __name__ == '__main__':
    main()

