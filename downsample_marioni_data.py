from pathlib import Path
import csv
import pandas as pd
import pdb
from scipy.io import mmread

here = Path(__file__).parent
data = here / "../Data"
cell_types_map_file = data / "cell_types_map.csv"
raw_data = data / "Raw/atlas"
max_in_type = 200


def get_cell_types_map():
    try:
        with open(cell_types_map_file, 'r') as f:
            cell_types_map = {row[0]:[int(el) for el in row[1:]] for row in csv.reader(f)}
        return cell_types_map
    except FileNotFoundError:
        pass

    metadata = pd.read_csv(
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

    cell_types_map = {}

    while len(metadata) > 0:
        row = metadata.sample(1)
        celltype = row["celltype"].values[0]
        row_indx = row.index[0]
        metadata.drop(row_indx)
        try:
            cell_types_map[celltype].append(row_indx)
        except KeyError:
            cell_types_map[celltype] = [row_indx]
        if len(cell_types_map[celltype]) == max_in_type:
            metadata = metadata[metadata["celltype"] != celltype]

    with open(cell_types_map_file, 'w') as f:
        w = csv.writer(f)
        for celltype, rows in cell_types_map.items():
            w.writerow([celltype]+rows)

    return cell_types_map

def get_norm_expr(cell_types_map):
    rows = sorted([el for subl in cell_types_map.values() for el in subl])
    raw_counts = mmread((raw_data/"raw_counts.mtx").__str__()).tocsr()[:,rows]
    # <class 'scipy.sparse.coo.coo_matrix'>
    # https://stackoverflow.com/questions/7609108/slicing-sparse-scipy-matrix
    # 
    df = pd.SparseDataFrame(raw_counts).fillna(0)
    pdb.set_trace()
    return rows

if __name__ == "__main__":
    cell_types_map = get_cell_types_map()
    norm_expr = get_norm_expr(cell_types_map)

