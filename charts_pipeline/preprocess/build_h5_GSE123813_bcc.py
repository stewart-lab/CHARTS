from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py
from os.path import join
from anndata import AnnData
import scanpy as sc
from collections import Counter

def main():
    root = args[0]
    out_f = args[1]

    the_genes = None
    counts = None
    cells = []
    tumors = []

    print("Loading data...")
    df = pd.read_csv(
        join(root, 'bcc_scRNA_counts.txt.gz'),
        sep='\t',
        index_col=0
    )

    df = df.transpose()

    from collections import defaultdict
    cell_to_patient = {
        cell: cell.split('.')[1]
        for cell in df.index
    }
    patient_to_cells = defaultdict(lambda: [])
    for cell, patient in cell_to_patient.items():
        patient_to_cells[patient].append(cell)

    patient_to_pre_cells = defaultdict(lambda: [])
    patient_to_post_cells = defaultdict(lambda: [])
    for patient, cells in patient_to_cells.items():
        for cell in cells:
            is_pre = 'pre' in cell.split('.')[2]
            is_post = 'post' in cell.split('.')[2]
            if not is_pre and not is_post:
                print("Cell {} has neither pre/post".format(cell))
            if is_pre and is_post:
                print("Cell {} has BOTH pre/post".format(cell))
            
            if is_pre:
                patient_to_pre_cells[patient].append(cell)
            elif is_post:
                patient_to_post_cells[patient].append(cell)

    genes = [
        x.encode('utf-8')
        for x in df.columns
    ]

    out_f = '../../charts_db/charts.h5'
    for patient, cells in patient_to_pre_cells.items():
        del_tum = 'GSE123814.{}.pre_treatment'
        tumor = 'GSE123814.{}.pre_treatment'.format(patient)
        print("Processing {} cells in sample {}...".format(len(cells), tumor))
        df_sample = df.loc[cells]
        cell_ids = [
            x.encode('utf-8')
            for x in cells
        ]
        X = np.array(df_sample)
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
                del f['per_tumor/{}'.format(del_tum)]
            except:
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

    for patient, cells in patient_to_post_cells.items():
        tumor = 'GSE123814.{}.post_treatment'.format(patient)
        print("Processing {} cells in sample {}...".format(len(cells), tumor))
        df_sample = df.loc[cells]
        cell_ids = [
            x.encode('utf-8')
            for x in cells
        ]
        X = np.array(df_sample)
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
