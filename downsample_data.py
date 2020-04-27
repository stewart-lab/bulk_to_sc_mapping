from pathlib import Path
import pandas as pd
import pdb

here = Path(__file__).parent
data = here / "../Data"
raw_data = data / "Raw/atlas"
max_in_type = 200

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

[(el, len(cell_types_map[el])) for el in cell_types_map.keys()]

pdb.set_trace()
