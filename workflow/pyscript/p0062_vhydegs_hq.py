#!/usr/bin/env python

import click
import os
import pandas as pd
from Bio import SeqIO


def extract_viral_nonviral_hits(fin_blast, topn):
    """
    For blast result, select top N non-viral hits and all viral hits as `hit_sel`.
    The rest hits are selected as `hit4cluster`.
    """
    df = pd.read_csv(fin_blast, sep="\t")
    df = df.sort_values(by="bitscore", ascending=False)

    # select virus kingdom
    df_virus = df[df["kingdom"] == "Viruses"]
    df_nonvirus = df[df["kingdom"] != "Viruses"]
    topn_nonvirus_hit = df_nonvirus["target_id"].tolist()[:topn]

    # check if there is a virus hit
    if df_virus.shape[0] > 0:
        virus_hit = df_virus["target_id"].tolist()
        hit_sel = virus_hit + topn_nonvirus_hit
    else:
        hit_sel = topn_nonvirus_hit

    hit4cluster = list(set(df["target_id"].tolist()) - set(hit_sel))

    return hit_sel, hit4cluster


@click.command()
@click.option("--fqry", '-i', help="input query protein file")
@click.option("--topn", '-n', default=10, type=int, help="top N hits to keep")
def main(fqry, topn):
    wd = os.getcwd()
    fid = ".".join(fqry.split(".")[:-1])
    fdir = os.path.join(wd, "rst_{}".format(fid))
    fblast = os.path.join(fdir, "{}_blastp_nr2.tsv".format(fid))
    fnr = os.path.join(fdir, "{}_blastp_nr.faa".format(fid))

    # Read qry FAA and all target FAAs
    rec_qry = SeqIO.read(fqry, format = "fasta")
    rec_target = list(SeqIO.parse(fnr, format = "fasta"))

    # For blast result, save top N non-viral hits and all viral hits as `rec_sel`.
    # The rest hits are saved as `rec_hit4cluster`.
    hit_sel, hit4cluster = extract_viral_nonviral_hits(fblast, topn)
    rec_sel = [rec_qry] + [x for x in rec_target if x.id in hit_sel]
    rec_hit4cluster = [x for x in rec_target if x.id in hit4cluster]

    fout_sel = os.path.join(fdir, "{}_sel.faa".format(fid))
    fout_hit4cluster = os.path.join(fdir, "{}_hit4cluster.faa".format(fid))

    SeqIO.write(rec_sel, fout_sel, format="fasta")
    SeqIO.write(rec_hit4cluster, fout_hit4cluster, format="fasta")


if __name__ == '__main__':
    main()