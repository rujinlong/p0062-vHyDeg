#!/usr/bin/env python

import click
from Bio import SeqIO
import os
import pandas as pd


@click.command()
@click.option("--fin", '-i', help="vHYDEGs.tsv file")
@click.option("--dirout", '-o', default="concat", help="Output directory")
def main(fin, dirout):
    vhydegs = pd.read_csv(fin, sep="\t")
    vhydegs["qryid"] = vhydegs["query_id"].apply(lambda x: x.replace("|", "__"))
    vhydegs = vhydegs[vhydegs["vHYDEG_quality"] == "HQ"][["qryid", "target_gene"]]

    # Create vHYDEGs dictionary
    vdict = {}
    for i in vhydegs.values:
        protid = i[0]
        gene_name = i[1]
        if not vdict.get(gene_name):
            vdict[gene_name] = [protid]
        else:
            vdict[gene_name].append(protid)

    for gene_name in vdict.keys():
        recs_sel = []
        recs_hit4clust = []
        df_all = []
        dirout2 = "{}_{}".format(dirout, gene_name)
        if not os.path.exists(dirout2):
            print(dirout2)
            os.mkdir(dirout2)
        for fid in vdict[gene_name]:
            dirrst = "rst_{}".format(fid)
            fin_hit4clust = os.path.join(dirrst, "{}_hit4cluster.faa".format(fid))
            fin_sel = os.path.join(dirrst, "{}_sel.faa".format(fid))
            fin_tsv = os.path.join(dirrst, "{}_blastp_nr2.tsv".format(fid))
            # check if fin_hit4clust exists
            if os.path.exists(fin_hit4clust):
                rec_hit4clust = list(SeqIO.parse(fin_hit4clust, "fasta"))
                rec_sel = list(SeqIO.parse(fin_sel, "fasta"))
                df = pd.read_csv(fin_tsv, sep="\t")
                recs_sel += rec_sel
                recs_hit4clust += rec_hit4clust
                df_all.append(df)

        recs_sel_nonredundant = []
        for rec in recs_sel:
            if rec.id not in [r.id for r in recs_sel_nonredundant]:
                recs_sel_nonredundant.append(rec)
        
        recs_hit4clust_nonredundant = []
        for rec in recs_hit4clust:
            if rec.id not in [r.id for r in recs_hit4clust_nonredundant]:
                recs_hit4clust_nonredundant.append(rec)
        if len(df_all) > 0:
            df_all = pd.concat(df_all)
            df_all.to_csv(os.path.join(dirout2, "{}_blastp_nr2.tsv".format(gene_name)), sep="\t", index=False)

        fout_sel = os.path.join(dirout2, "{}_sel.faa".format(gene_name))
        fout_hit4clust = os.path.join(dirout2, "{}_hit4cluster.faa".format(gene_name))
        SeqIO.write(recs_sel_nonredundant, fout_sel, "fasta")
        SeqIO.write(recs_hit4clust_nonredundant, fout_hit4clust, "fasta")


if __name__ == '__main__':
    main()