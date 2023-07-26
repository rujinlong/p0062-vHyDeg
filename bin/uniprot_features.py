#!/usr/bin/env python

import click
import json
import pandas as pd
import requests


def download_uniprot_features(uid):
    url = "https://www.ebi.ac.uk/proteins/api/features/" + uid
    response = requests.get(url)
    try:
        anno = response.json()
        return anno
    except:
        print("Error: ", response.status_code)


def extract_uniprot_features(uniprot_id, protein_id, kws=["BINDING", "DOMAIN"]):
    anno = download_uniprot_features(uniprot_id)
    features = []
    for feature in anno["features"]:
        if feature["type"] in kws:
            features.append([protein_id, uniprot_id, feature["begin"], feature["end"], feature["type"]])

    df = pd.DataFrame(features, columns=["prot_id", "uniprot_id", "begin", "end", "feature_type"])
    return df


@click.command()
@click.option("--protein_id", '-i', help="Query protein id")
@click.option("--uniprot_id", '-u', help="Target uniprot id")
def main(protein_id, uniprot_id):
    df = extract_uniprot_features(uniprot_id, protein_id)
    df.to_csv(f"{protein_id}_features_{uniprot_id}.tsv", index=False, sep="\t")

if __name__ == '__main__':
    main()