from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py
from os.path import join
from anndata import AnnData
import scanpy as sc
from collections import defaultdict

def main():
    root = args[0]
    out_f = args[1]

    the_genes = None
    counts = None
    cells = []
    tumors = []

    print("Loading data...")
    df = pd.read_csv(
        join(root, 'merge10pts_counts.txt.gz'),
        sep='\t',
        index_col=0
    )
    df = df.transpose()
    df = df.loc[df['Tissue: 0=Normal, 1=Tumor'] == 1]
    print("done.")

    tumor_to_cells = defaultdict(lambda: [])
    for cell_id in df.index:
        tumor = cell_id.split('_')[0]
        tumor_to_cells[tumor].append(cell_id)
    print({k: len(v) for k,v in tumor_to_cells.items()})

    genes = [
        x.encode('utf-8')
        for x in df.columns[2:]
    ]

    for tumor_id, cells in tumor_to_cells.items():
        tumor = 'GSE144236.{}'.format(tumor_id)
        print("Processing {} cells in sample {}...".format(len(cells), tumor))
        df_sample = df.loc[cells][df.columns[2:]]
        print(df_sample.shape)
        cell_ids = [
            x.encode('utf-8')
            for x in cells
        ]
        X = np.array(df_sample, dtype=np.float32)
        ad = AnnData(
            X=X
        )
        sc.pp.normalize_total(ad, target_sum=1e6)
        sc.pp.log1p(ad)
        log1_tpm = ad.X
        with h5py.File(out_f, 'r+') as f:
            try:
                del f['per_tumor/{}'.format(tumor)]
            except KeyError:
                pass
            try:
                f.create_group('per_tumor')
            except:
                pass
            f['per_tumor'].create_group(tumor)
            f['per_tumor'][tumor].create_dataset(
                'log1_tpm',
                data=log1_tpm,
                compression="gzip"
            )
            f['per_tumor'][tumor].create_dataset(
                'cell',
                data=np.array(cell_ids),
                compression="gzip"
            )
            f['per_tumor'][tumor].create_dataset(
                'gene_name',
                data=np.array(genes),
                compression="gzip"
            )

if __name__ == '__main__':
    main()
