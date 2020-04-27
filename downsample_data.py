from pathlib import Path
import pandas as pd

here = Path(__file__).parent
data = here / "../Data"
raw_data = data / "Raw/atlas"

metadata = pd.read_csv(raw_data / "meta.csv",
        dtype = {})
'''
            'cell':,
            barcode':,
            sample,
            stage,
            sequencing.batch,
            theiler,
            doub.density,
            doublet,
            cluster,
            cluster.sub,
            cluster.stage,
            cluster.theiler,
            stripped,
            celltype,
            colour,
            umapX,
            umapY,
            haem_gephiX,
            haem_gephiY,
            haem_subclust,
            endo_gephiX,
            endo_gephiY,
            endo_trajectoryName,
            endo_trajectoryDPT,
            endo_gutX,
            endo_gutY,
            endo_gutDPT,
            endo_gutCluster)
            '''
