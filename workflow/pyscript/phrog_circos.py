import copy
import pandas as pd
from Bio import SeqIO, SeqFeature
import pycircos

def get_gbk(fgbk, sample_id=False):
    records = list(SeqIO.parse(fgbk, "genbank"))
    # only return the first record
    record = records[0]
    if sample_id:
        record.id = sample_id
        record.name = sample_id
    record.features.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(0, len(record)), type="source", strand=1))
    return record


def get_cds(record):
    # feature locations store in a list [feature_id, start, end, strand]
    feature_locations = []
    for feature in record.features:
        if feature.type == "CDS":
            feature_locations.append([feature.qualifiers["locus_tag"][0], feature.qualifiers["phrog"][0], feature.location.start.real, feature.location.end.real])

    # convert to dataframe
    df_cds = pd.DataFrame(feature_locations, columns=["feature_id", "phrog", "start", "end"])
    # width of each feature
    df_cds["width"] = df_cds["end"] - df_cds["start"] + 1
    return df_cds


def read_phrog_annot(fin_phrog):
    df_phrog = pd.read_csv(fin_phrog, sep="\t")
    # convert column "phrog" to str
    df_phrog["phrog"] = df_phrog["phrog"].astype(str)
    return df_phrog


def add_phrog_to_cds(df_cds, df_phrog):
    df_cds = df_cds.merge(df_phrog, on="phrog", how="left")
    # fill NaN in column "color" with #222223
    df_cds["color"] = df_cds["color"].fillna("#c9c9c9")
    # fill NaN in column "category" with "Unknown"
    df_cds["category"] = df_cds["category"].fillna("unknown function")
    return df_cds


def featureplot_with_bg(self, arc_id, source, feat, raxis_start, feat_width, facecolor_feat, facecolor_source="#ecfdfd"):
    self.featureplot(arc_id, source=source, raxis_range=(raxis_start, raxis_start + feat_width), facecolor=facecolor_source)
    self.featureplot(arc_id, source=feat,  raxis_range=(raxis_start, raxis_start + feat_width), facecolor=facecolor_feat)
    return self


def featureplot_gcskew(sample_id, garc, gcircle, raxis_start, feat_width, window_size=500):
    skews = garc.calc_nnskew(n1="G", n2="C", window_size=window_size)
    positive_skews = copy.deepcopy(skews)
    positive_skews[skews<0] = 0
    negative_skews = copy.deepcopy(skews)
    negative_skews[skews>=0] = 0
    gcircle.lineplot(sample_id, positive_skews, rlim=(min(skews), max(skews)), raxis_range=(raxis_start, raxis_start + feat_width), linecolor="#275ba5")
    gcircle.lineplot(sample_id, negative_skews, rlim=(min(skews), max(skews)), raxis_range=(raxis_start, raxis_start + feat_width), linecolor="#560061")
    return gcircle


def featureplot_gccontent(sample_id, garc, gcircle, raxis_start, feat_width, window_size=500):
    gc_content = garc.calc_nnratio(n1="G", n2="C", window_size=window_size)
    gcircle.lineplot(sample_id, gc_content, rlim=(min(gc_content),max(gc_content)), raxis_range=(raxis_start, raxis_start + feat_width), linecolor="#009CDE", linewidth=0.5)
    return gcircle


def extract_gbk_features(record):
    feat_source = []
    feat_CDS_plus  = []
    feat_CDS_minus = []
    feat_tRNA = []
    for feat in record.features:
        if feat.type == "source":
            feat_source.append(feat)
        if "CDS" in feat.type:
            if feat.strand == 1:
                feat_CDS_plus.append(feat)
            else:
                feat_CDS_minus.append(feat)
        if feat.type == "tRNA":
            feat_tRNA.append(feat)

    features = {"source": feat_source,
                "CDS_plus": feat_CDS_plus,
                "CDS_minus": feat_CDS_minus,
                "CDS": feat_CDS_plus + feat_CDS_minus,
                "tRNA": feat_tRNA}
    
    return features


def plot_phrog_circos(sample_id, record, fin_phrog, tick_interval, tick_label_size, window_size, protid_sulfur, raxmin=230, gap=10, figsize=8, arc_degree=5, dict_rbh=False, cdscolor="#E1F5FE", labelsize=10):
    feats = extract_gbk_features(record)
    df_cds = get_cds(record)
    df_phrog = read_phrog_annot(fin_phrog)
    df_cds_anno = add_phrog_to_cds(df_cds, df_phrog)

    # pycircos settings
    tick_max = len(record) + 1000
    rax = raxmin

    # pycircos plot
    garc = pycircos.Garc(arc_id=sample_id, record=record, interspace=0, linewidth=0, facecolor="black", raxis_range=(0,0), label="{}\n\n{} bp".format(sample_id, len(record)), label_visible=True, labelsize=labelsize)
    gcircle = pycircos.Gcircle(figsize=(figsize,figsize))
    gcircle.add_garc(garc)
    gcircle.set_garcs(arc_degree, 360-arc_degree)

    # plot gc skew + gc content
    gcircle = featureplot_gcskew(sample_id, garc, gcircle, rax, 70, window_size=window_size)
    rax += 70
    gcircle = featureplot_gccontent(sample_id, garc, gcircle, rax, 60, window_size=window_size)
    rax += 60

    # Plot tRNA
    gcircle = featureplot_with_bg(gcircle, sample_id, feats["source"], feats["tRNA"], rax, 40, "#AD1AAC")
    rax += 40 + gap

    # Plot CDS
    gcircle.featureplot(sample_id, source=feats["CDS_minus"],  raxis_range=(rax, rax+40), facecolor=cdscolor)
    gcircle.featureplot(sample_id, source=feats["CDS_plus"], raxis_range=(rax+40, rax+80), facecolor=cdscolor)
    if dict_rbh:
        # add sulfur gene annotation
        feats_sulfur_minus = [feat for feat in feats["CDS_minus"] if feat.qualifiers["locus_tag"][0] in protid_sulfur]
        if len(feats_sulfur_minus) > 0:
            gcircle.featureplot(sample_id, source=feats_sulfur_minus,  raxis_range=(rax, rax + 40), facecolor="#FFD600")
        feats_sulfur_plus = [feat for feat in feats["CDS_plus"] if feat.qualifiers["locus_tag"][0] in protid_sulfur]
        if len(feats_sulfur_plus) > 0:
            gcircle.featureplot(sample_id, source=feats_sulfur_plus, raxis_range=(rax + 40, rax + 80), facecolor="#FFD600")


        # add HYDG annotation
        hydg_locus_tag = [x["id_new"] for x in dict_rbh[record.id]]
        feats_hydg_minus = [feat for feat in feats["CDS_minus"] if feat.qualifiers["locus_tag"][0] in hydg_locus_tag]
        if len(feats_hydg_minus) > 0:
            gcircle.featureplot(sample_id, source=feats_hydg_minus,  raxis_range=(rax, rax + 40), facecolor="#D50000")
        feats_hydg_plus = [feat for feat in feats["CDS_plus"] if feat.qualifiers["locus_tag"][0] in hydg_locus_tag]
        if len(feats_hydg_plus) > 0:
            gcircle.featureplot(sample_id, source=feats_hydg_plus, raxis_range=(rax + 40, rax + 80), facecolor="#D50000")

    rax += 80 + gap

    # PHROGs
    gcircle.barplot(sample_id, data=[1]*len(df_cds_anno["start"]), positions=df_cds_anno["start"], width=df_cds_anno["width"], raxis_range=(rax, rax+70), facecolor=df_cds_anno["color"])
    rax += 70

    # add tickplot
    gcircle.featureplot(sample_id, source=feats["source"], raxis_range=(rax, rax+5), facecolor="#9E9E9E")
    rax += 5
    tlabels = [f"{tl / 1000:.0f} K" for tl in range(0, tick_max, tick_interval)]
    gcircle.tickplot(sample_id, raxis_range=(rax, rax+3), tickinterval=tick_interval, tickcolor="#9E9E9E", ticklabels=tlabels, ticklabelsize=tick_label_size, ticklabelorientation="horizontal")
    return gcircle