import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scanpy as sc
import sys
import os 
from os.path import join
import subprocess
from anndata import AnnData
from optparse import OptionParser
from collections import defaultdict
import h5py
import json

WINDOW_SIZE = 100

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-l", "--log_dir", help="Directory in which to write log files")
    parser.add_option("-w", "--overwrite", action='store_true', help="Overwrite data in the HDF5 file if there's a dataset already present")
    (options, args) = parser.parse_args()

    h5_f = args[0]
    tumor_sets_str = args[1]
    tumor_set_name = args[2]
    overwrite = options.overwrite
    out_dir = options.out_dir

    tumor_set = tumor_sets_str.split(',')
    print('Tumor set: ', tumor_set)

    # Map each chromosome to its genes and sort the genes
    # within each chromosome
    df_gene_loc = pd.read_csv('gene_locations.tsv', index_col=0, sep='\t')
    chrm_to_genes = defaultdict(lambda: set())
    for gene, chrm in zip(df_gene_loc.index, df_gene_loc['Chromosome/scaffold name']):
        chrm_to_genes[chrm].add(gene)
    chrm_to_ordered_genes = {}
    print("Sorting genes within each chromosome...")
    for chrm, genes in chrm_to_genes.items():
        locs = df_gene_loc.loc[genes]['Gene start (bp)']
        sorted_genes = sorted(
            [
                (gene, loc)
                for gene, loc in zip(genes, locs)
            ],
            key=lambda x: x[1]
        )
        sorted_genes = [x[0] for x in sorted_genes]
        chrm_to_ordered_genes[chrm] = sorted_genes
        
    # Remove all the contigs
    chrm_to_ordered_genes = {
        chrm: genes
        for chrm, genes in chrm_to_genes.items()
        if chrm.isnumeric() or chrm == 'X' or chrm == 'Y' or chrm == 'MT'
    }

    # Map each gene ID to it's gene Name
    with open('gene_id_to_name.json', 'r') as f:
        gene_id_to_name = json.load(f)
    gene_name_to_ids = defaultdict(lambda: [])
    for gene_id, gene_name in gene_id_to_name.items():
        gene_name_to_ids[gene_name].append(gene_id)

    # Determine whether to use gene names or ID's
    use_gene_id = False
    gene_ids_present = set()
    with h5py.File(h5_f, 'r') as f:
        for tumor in tumor_set:
            if 'gene_id' in f['per_tumor/{}'.format(tumor)].keys():
                gene_ids_present.add(True)
            else:
                gene_ids_present.add(False)
    if False in gene_ids_present:
        with open('gene_id_to_name.json', 'r') as f:
            gene_id_to_name = json.load(f)
        gene_name_to_ids = defaultdict(lambda: [])
        for gene_id, gene_name in gene_id_to_name.items():
            gene_name_to_ids[gene_name].append(gene_id)
    else:
        gene_id_to_name = {
            x: x
            for x in gene_id_to_name.keys()
        }

    tumor_dfs = []
    with h5py.File(h5_f, 'r') as f:
        for tumor in tumor_set:
            cells = [
                str(x)[2:-1]
                for x in f['per_tumor/{}/cell'.format(tumor)][:]
            ]
            if 'gene_id' in f['per_tumor/{}'.format(tumor)].keys():
                genes = [
                    str(x)[2:-1]
                    for x in f['per_tumor/{}/gene_id'.format(tumor)][:]
                ]
                if '.' in genes[0]:
                    print("Detected gene versions. Removing versions...")
                    genes = [
                        g.split('.')[0]
                        for g in genes
                    ]
            else:
                genes = [
                    str(x)[2:-1]
                    for x in f['per_tumor/{}/gene_name'.format(tumor)][:]
                ]
            expression = f['per_tumor/{}/log1_tpm'.format(tumor)][:]
            df = pd.DataFrame(
                data=expression,
                index=cells,
                columns=genes
            )
            tumor_dfs.append(df)

    # Select genes that are common to all tumors in the tumor set
    common_genes = None
    for df in tumor_dfs:
        if common_genes is None:
            common_genes = set(df.columns)
        else:
            common_genes &= set(df.columns)
    common_genes = sorted(common_genes)
    print("Found {} genes common to all tumors".format(len(common_genes)))
    tumor_dfs = [
        df[common_genes]
        for df in tumor_dfs
    ]

    # Compute the available genes that are common between the dataset
    # and the genomic coordinates file
    print("Filtering chromosome orderings by genes in dataset...")
    common_genes_s = set(common_genes)
    chrm_to_ordered_avail_genes = {
        chrm: [
            gene_id_to_name[gene]
            for gene in ordered_genes
            if gene in gene_id_to_name
            and gene_id_to_name[gene] in common_genes_s
        ]
        for chrm, ordered_genes in chrm_to_ordered_genes.items()
    }
    print("done")

    # Map gene to index
    gene_to_index = {
        gene: index
        for index, gene in enumerate(common_genes)
    }

    # Compute window-averages
    half_window = int(WINDOW_SIZE / 2)
    tumor_to_window_X = {}
    for tumor, df in zip(tumor_set, tumor_dfs):
        print("Computing window averages for tumor ", tumor)
        X = np.array(df)
        window_sums_X = []
        for chrm in sorted(chrm_to_ordered_avail_genes.keys()):
            print("Passing window over chromosome ", chrm)
            sorted_genes = chrm_to_ordered_avail_genes[chrm]
            for gene_i in np.arange(half_window, len(sorted_genes)-half_window):
                if (gene_i+1) % 500 == 0:
                    print("Chromosome {}, gene {}/{}".format(chrm, gene_i+1, len(sorted_genes)-2*half_window))
                window_indices = [
                    gene_to_index[gene]
                    for gene in sorted_genes[gene_i-half_window: gene_i+half_window]
                ]
                X_window = X[:,window_indices]
                X_window = np.exp(X_window)-1
                x_window = np.mean(X_window, axis=1)
                window_sums_X.append(x_window)
        window_sums_X = np.array(window_sums_X).T
        window_sums_X = np.log(window_sums_X + 1)
        print("Computed average-window-matrix for tumor {}. Shape is: {}".format(tumor, window_sums_X.shape))
        tumor_to_window_X[tumor] = window_sums_X

    # Cluster the tumors by their sliding window averages
    tumor_to_cells = {
        tumor: list(df.index)
        for tumor, df in zip(tumor_set, tumor_dfs)
    }
    ad = build_window_ad(tumor_to_window_X, tumor_to_cells) 
    print("Clustering tumors...")
    sc.pp.pca(ad)
    sc.pp.neighbors(ad)
    sc.tl.leiden(ad, resolution=4.0)

    if options.log_dir:
        sc.tl.umap(ad)
        fig, ax = plt.subplots(1,1,figsize=(8,6))
        sc.pl.umap(ad, color='tumor', size=8, ax=ax)
        plt.tight_layout()
        fig.savefig(
            join(
                options.log_dir, 
                '{}_malignant_score_clustering.umap_color_by_tumor.png'.format(tumor_set_name)
            ),
            format='png',
            dpi=150
            #bbox_inches='tight'
        )
        fig, ax = plt.subplots(1,1,figsize=(8,6))
        sc.pl.umap(ad, color='leiden', size=8, ax=ax)
        plt.tight_layout()
        fig.savefig(
            join(
                options.log_dir,
                '{}_malignant_score_clustering.umap_color_by_cluster.png'.format(tumor_set_name)
            ),
            format='png',
            dpi=150
            #bbox_inches='tight'
        )

        

    # Compute malignancy scores
    print("Computing malignancy scores...")
    cell_to_entropy = compute_malignancy_score(ad, tumor_to_cells)
    tumor_to_scores = {}
    for tumor, tumor_df in zip(tumor_set, tumor_dfs):
        tum_cells = tumor_to_cells[tumor]
        df = pd.DataFrame(
            data={
                'entropy': [
                    cell_to_entropy[cell]
                    for cell in tum_cells
                ]
            },
            index=tum_cells
        )
        df = df.loc[tumor_df.index]
        scores = np.array(df)
        scores = np.squeeze(scores)
        tumor_to_scores[tumor] = scores

    print(tumor_to_scores)

    with h5py.File(h5_f, 'r+') as f:
        for tumor in tumor_set:
            key = 'per_tumor/{}/malignancy_score'.format(tumor)
            scores = tumor_to_scores[tumor]
            if key not in f.keys():
                f.create_dataset(
                    key,
                    data=scores,
                    compression='gzip'
                )
            elif overwrite:
                try:
                    del f[key]
                except KeyError:
                    pass
                f.create_dataset(
                    key,
                    data=scores,
                    compression='gzip'
                )
            else:
                print("No data written for tumor {}. Dataset '{}' already present in HDF5 file.".format(
                    tumor,
                    key
                ))


def build_window_ad(tumor_to_window_X, tumor_to_cells):
    the_tumors = list(tumor_to_window_X.keys())
    X = np.concatenate([
        tumor_to_window_X[tumor]
        for tumor in the_tumors
    ])
    tumors = np.concatenate([
        np.full(len(tumor_to_window_X[tumor]), str(tumor))
        for tumor in the_tumors
    ])
    cells = np.concatenate([
        np.array(tumor_to_cells[tumor])
        for tumor in the_tumors
    ])
    ad = AnnData(
        X=X,
        obs=pd.DataFrame(
            data={
                'tumor': tumors
            },
            index=cells
        )
    )
    return ad


def compute_malignancy_score(ad, tumor_to_cells):
    clust_to_entropy = {}
    cell_to_entropy = {}
    for clust in sorted(set(ad.obs['leiden'])):
        clust_cells = ad.obs.loc[ad.obs['leiden'] == clust].index
        clust_tumor_to_frac = {
            tumor: len(set(tumor_to_cells[tumor]) & set(clust_cells)) / len(clust_cells)
            for tumor in tumor_to_cells.keys()
        }
        entropy = sum([
            -clust_tumor_to_frac[tumor] * np.log(clust_tumor_to_frac[tumor])
            for tumor in clust_tumor_to_frac
            if clust_tumor_to_frac[tumor] > 0
        ])
        print('Computed entropy of {} for cluster {}'.format(entropy, clust))
        for cell in clust_cells:
            cell_to_entropy[cell] = entropy
    return cell_to_entropy

if __name__ == '__main__':
    main()

