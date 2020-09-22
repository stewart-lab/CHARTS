import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import sys
import os 
from os.path import join
import subprocess
from anndata import AnnData
import phate
from optparse import OptionParser
import json
import h5py

sys.path.append('..')

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-w", "--overwrite", action="store_true", help="Overwrite data if already present")
    parser.add_option("-c", "--cluster_id", help="Cluster ID to restrict plot")
    parser.add_option("-l", "--cluster_file", help="Files storing cluster assignments. Required if '-c' flat is used.")
    (options, args) = parser.parse_args()

    h5_f = args[0]
    overwrite = options.overwrite

    sc.settings.verbosity = 3
    sc.logging.print_versions()

    the_tumors = set()
    with h5py.File(h5_f, 'r') as f:
        the_tumors = f['per_tumor'].keys()
        the_tumors = sorted(the_tumors)
    print(the_tumors)

    for tumor in the_tumors:
        print("Running algorithms on tumor {}".format(tumor))
        with h5py.File(h5_f, 'r') as f:
            if not overwrite and 'umap_2'.format(tumor) in f['per_tumor/{}'.format(tumor)].keys():
                print('Detected tumor {} already in data. Skipping running dimensionality reduction.'.format(tumor))
                continue
            cells = [
                str(x)[2:-1]
                for x in f['per_tumor/{}/cell'.format(tumor)][:]
            ]
            expression = f['per_tumor/{}/log1_tpm'.format(tumor)][:]
        
        ad = AnnData(
            X=expression,
            obs=pd.DataFrame(data=cells, columns=['cell'])
        )

        try:
            X_phate_2 = run_phate(ad, 2)
            X_umap_2 = run_umap(ad, 2)
            X_phate_3 = run_phate(ad, 3)
            X_umap_3 = run_umap(ad, 3)

            with h5py.File(h5_f, 'r+') as f:
                try:
                    del f['per_tumor/{}/umap_2'.format(tumor)]
                except KeyError:
                    pass
                f.create_dataset(
                    'per_tumor/{}/umap_2'.format(tumor),
                    data=np.array(X_umap_2, dtype=np.float32),
                    compression="gzip"
                )
                try:
                    del f['per_tumor/{}/phate_2'.format(tumor)]
                except KeyError:
                    pass
                f.create_dataset(
                    'per_tumor/{}/phate_2'.format(tumor),
                    data=np.array(X_phate_2, dtype=np.float32),
                    compression="gzip"
                )
                try:
                    del f['per_tumor/{}/umap_3'.format(tumor)]
                except KeyError:
                    pass
                f.create_dataset(
                    'per_tumor/{}/umap_3'.format(tumor),
                    data=np.array(X_umap_3, dtype=np.float32),
                    compression="gzip"
                )
                try:
                    del f['per_tumor/{}/phate_3'.format(tumor)]
                except KeyError:
                    pass
                f.create_dataset(
                    'per_tumor/{}/phate_3'.format(tumor),
                    data=np.array(X_phate_3, dtype=np.float32),
                    compression="gzip"
                )
        except ValueError as e:
            print("Error processing tumor {}:".format(tumor))
            print(e)

def run_phate(ad, n_comps):
    print("Running PHATE with {} dimensions...".format(n_comps))
    phate_operator = phate.PHATE(
        n_jobs=-2,
        random_state=1,
        n_components=n_comps
    )
    X_phate = phate_operator.fit_transform(ad.X)
    return X_phate


def run_umap(ad, n_comps):
    print("Running UMAP with {} dimensions...".format(n_comps))
    sc.tl.pca(ad, n_comps=50)
    if ad.X.shape[0] < 300:
        n_neighbors = int(ad.X.shape[0] * 0.05)
    else:
        n_neighbors = 15
    sc.pp.neighbors(ad, n_neighbors=n_neighbors)
    sc.tl.umap(ad, n_components=n_comps)
    return ad.obsm['X_umap']    

if __name__ == '__main__':
    main()

