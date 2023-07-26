#!/usr/bin/env python

import click
import pandas as pd
import os
import shutil
from Bio import SeqIO

def filter_orf2fam(fin_orf2fam, nrqry_list):
    df_orf2fam = pd.read_csv(fin_orf2fam, sep="\t")
    df_orf2fam = df_orf2fam[df_orf2fam['orf'].isin(nrqry_list)]
    # df_orf2fam["fname"] = df_orf2fam.apply(lambda x: x["orf"].replace("|", "__"), axis=1)
    # df_orf2fam["fname2"] = df_orf2fam.apply(lambda x: x["fname"] + "__" + x["family"], axis=1)

    # fam2orf dict
    orf2fam = df_orf2fam.set_index('orf')['family'].to_dict()
    return orf2fam


@click.command()
@click.option("--fin_nrqry", '-i', help="file of nrqry list")
@click.option("--dirphams", '-d', help="directory of phams")
def main(fin_nrqry, dirphams):
    fin_orf2fam = os.path.join(dirphams, "orf2family.tsv")
    dir_fam = os.path.join(dirphams, "familiesFasta")
    dir_rst = os.path.join(dirphams, "familiesRename")
    # create directory dir_rst if not exist
    if not os.path.exists(dir_rst):
        os.mkdir(dir_rst)

    records = list(SeqIO.parse(fin_nrqry, "fasta"))
    nrqry_list = [x.id for x in records]
    records_dict = {x.id: x for x in records}

    orf2fam = filter_orf2fam(fin_orf2fam, nrqry_list)
    for nrqry_id in nrqry_list:
        fid = nrqry_id.replace("|", "__")
        famid = orf2fam.get(nrqry_id)
        if famid is None:
            print("no " + nrqry_id)
            fname_rst = os.path.join(dir_rst, fid + "_fam00.faa")
            SeqIO.write(records_dict[nrqry_id], fname_rst, "fasta")
        else:
            fname_rst = os.path.join(dir_rst, "{}_{}.faa".format(fid, famid))
            fname_fam = os.path.join(dir_fam, famid + ".fa")
            shutil.copyfile(fname_fam, fname_rst)
        print(fname_rst)
    


if __name__ == '__main__':
    main()