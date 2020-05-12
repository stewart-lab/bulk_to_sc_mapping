from pathlib import Path
import csv
import pandas as pd
import pdb
import pickle
from scipy.io import mmread

here = Path(__file__).parent
data = here / "../Data"
cell_types_map_file = data / "cell_types_map.csv"
raw_data = data / "Raw/atlas"
max_in_type = 2  # 00
max_genes = 5

import warnings

warnings.filterwarnings("error")


def get_metadata():
    df = pd.read_csv(
        raw_data / "meta.csv",
        usecols=[13],
        dtype={
            "cell": object,
            "barcode": object,
            "sample": int,
            "stage": object,
            "sequencing.batch": int,
            "theiler": object,
            "doub.density": float,
            "doublet": bool,
            "cluster": float,
            "cluster.sub": float,
            "cluster.stage": float,
            "cluster.theiler": float,
            "stripped": bool,
            "celltype": object,
            "colour": object,
            "umapX": float,
            "umapY": float,
            "haem_gephiX": float,
            "haem_gephiY": float,
            "haem_subclust": object,
            "endo_gephiX": float,
            "endo_gephiY": float,
            "endo_trajectoryName": object,
            "endo_trajectoryDPT": float,
            "endo_gutX": float,
            "endo_gutY": float,
            "endo_gutDPT": float,
            "endo_gutCluster": object,
        },
    ).dropna()
    df.index = df.index.rename("sample_num")
    return df


def get_samps(metadata):
    try:
        with open(cell_types_map_file, "r") as f:
            cell_types_map = {
                row[0]: [int(el) for el in row[1:]] for row in csv.reader(f)
            }
    except FileNotFoundError:
        cell_types_map = {}
        while len(metadata) > 0:
            row = metadata.sample(1)
            celltype = row["celltype"].values[0]
            row_indx = row.index[0]
            metadata = metadata.drop(row_indx)
            try:
                cell_types_map[celltype].append(row_indx)
            except KeyError:
                cell_types_map[celltype] = [row_indx]
            if len(cell_types_map[celltype]) == max_in_type:
                metadata = metadata[metadata["celltype"] != celltype]

        with open(cell_types_map_file, "w") as f:
            w = csv.writer(f)
            for celltype, rows in cell_types_map.items():
                w.writerow([celltype] + rows)

    return sorted([el for subl in cell_types_map.values() for el in subl])


def get_raw_counts(samps, metadata):
    try:
        with open("downsampled_df.pickle", "rb") as f:
            df = pickle.load(f)
    except FileNotFoundError:
        try:
            with open("raw_counts.pickle", "rb") as f:
                raw_counts = pickle.load(f)
        except FileNotFoundError:
            raw_counts = mmread(
                (raw_data / "raw_counts.mtx").__str__()
            ).tocsr()[:, samps]
            with open("raw_counts.pickle", "wb") as f:
                pickle.dump(raw_counts, f)
        df = pd.DataFrame(raw_counts.todense())
        df.columns = samps
        with open(raw_data / "genes.tsv", "r") as f:
            genes = [line.strip().split("\t")[1] for line in f]
        df.index = genes
        top_genes = (
            df.sum(axis=1).sort_values(ascending=False).head(max_genes).index
        )
        df = df.loc[top_genes, :]
        df.index =  df.index.rename('Genes')
        with open("downsampled_df.pickle", "wb") as f:
            pickle.dump(df, f)
    return df
    # load in size factors
    # Downsample result = dataframe.mul(series, axis=0)
    # norm df with them
    # save a cached copy
    # edit fn to look for cached copy first.
    # Downsample and save feature vector
    # Wait a moment, think that SingleCellExperiment might want raw counts, so:
    # Downsample and save size_factors vector


if __name__ == "__main__":
    metadata = get_metadata()
    samps = get_samps(metadata)
    raw_counts = get_raw_counts(samps, metadata)
    raw_counts.to_csv('raw_downsampled_counts.csv')
    metadata.loc[samps,:].to_csv('downsampled_celltypes.csv')
    pdb.set_trace()
