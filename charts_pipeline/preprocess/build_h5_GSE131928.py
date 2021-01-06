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
        join(root, 'GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz'),
        sep='\t',
        index_col=0
    )
    df = df.transpose()

    tumor_to_cells = defaultdict(lambda: [])
    for cell_id in df.index:
        raw_tum = cell_id.split('-')[0]
        if 'MGH' in raw_tum:
            tumor = raw_tum
        else:
            tumor = 'MGH{}'.format(raw_tum)
        tumor_to_cells[tumor].append(cell_id)
    tumor_to_cells = {
        tum: cells
        for tum, cells in tumor_to_cells.items()
        if 'CD45neg' not in tum
    }
    print({k: len(v) for k,v in tumor_to_cells.items()})
    print()

    # Aggregate tumor MGH105 cells
    MGH105_cells = []
    MGH105_tums = ['MGH105A', 'MGH105B', 'MGH105C', 'MGH105D']
    for tum in MGH105_tums:
        MGH105_cells += tumor_to_cells[tum]
    for tum in MGH105_tums:
        del tumor_to_cells[tum]
    tumor_to_cells['MGH105'] = MGH105_cells
        
    # Aggregate tumor MGH104 cells
    MGH104_cells = []
    MGH104_tums = ['MGH104', 'MGH104negP2', 'MGH104negP4', 'MGH104negP7']
    for tum in MGH104_tums:
        MGH104_cells += tumor_to_cells[tum]
    for tum in MGH104_tums:
        del tumor_to_cells[tum]
    tumor_to_cells['MGH104'] = MGH104_cells

    # Aggregate tumor MGH106 cells
    MGH106_cells = []
    MGH106_tums = ['MGH106', 'MGH106CD3posP1']
    for tum in MGH106_tums:
        MGH106_cells += tumor_to_cells[tum]
    for tum in MGH106_tums:
        del tumor_to_cells[tum]
    tumor_to_cells['MGH106'] = MGH106_cells

    print({k: len(v) for k,v in tumor_to_cells.items()})

    genes = [
        x.encode('utf-8')
        for x in df.columns
    ]   

    for tumor_id, cells in tumor_to_cells.items():
        tumor = 'GSE131928.{}'.format(tumor_id)
        print("Processing {} cells in sample {}...".format(len(cells), tumor))
        df_sample = df.loc[cells]
        cell_ids = [
            x.encode('utf-8')
            for x in cells
        ]
        X = np.array(df_sample, dtype=np.float32)
        #X = 10*(np.power(2, X)-1)
        ad = AnnData(
            X=X
        )
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












