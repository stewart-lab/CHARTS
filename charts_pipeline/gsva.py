import matplotlib as mpl
#mpl.use('Agg')
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
from collections import defaultdict
import gseapy as gp
import h5py
import run_gsva

sys.path.append('..')

GENE_SETS = ['GO_Biological_Process_2018']
#GENE_SETS = ['GO_Molecular_Function_2018']

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-w", "--overwrite", action='store_true', help="Overwrite data in the HDF5 file if there's a dataset already present")
    (options, args) = parser.parse_args()

    h5_f = args[0]
    gene_set_f = args[1]
    collection_name = args[2]
    out_dir = options.out_dir

    the_tumors = set()
    with h5py.File(h5_f, 'r') as f:
        the_tumors = f['per_tumor'].keys()
        the_tumors = sorted(the_tumors)

    print(the_tumors)

    gene_set_to_genes = run_gsva._parse_gene_sets(gene_set_f)

    for tumor in the_tumors:
        print("Computing enrichment scores for tumor {}".format(tumor))
        with h5py.File(h5_f, 'r') as f:
            cells = [
                str(x)[2:-1]
                for x in f['per_tumor/{}/cell'.format(tumor)][:]
            ]
            genes = [
                str(x)[2:-1]
                for x in f['per_tumor/{}/gene_name'.format(tumor)][:]
            ]
            clusters = f['per_tumor/{}/leiden_res_4/cluster'.format(tumor)][:]
            expression = f['per_tumor/{}/log1_tpm'.format(tumor)][:]

        all_clusts = sorted(set(clusters))

        # Map each cluster to its cells
        clust_to_cells = defaultdict(lambda: [])
        for clust, cell in zip(clusters, cells):
            clust_to_cells[clust].append(cell)

        # Map each cell to its index
        cell_to_index = {
            cell: index
            for index, cell in enumerate(cells)
        }

        expression_clusts = []
        clusts = []
        for clust, clust_cells in clust_to_cells.items():
            print("Aggregating cluster {}...".format(clust))
            indices = [
                cell_to_index[cell]
                for cell in clust_cells
            ]
            X_clust = expression[indices,:]
            # Convert back to TPM
            X_clust = np.exp(X_clust) - 1
            x_clust = np.sum(X_clust, axis=0)
            s = sum(x_clust)
            x_clust = x_clust / sum(x_clust)
            expression_clusts.append(x_clust)
            clusts.append(clust)
        expression_clusts = np.array(expression_clusts)
        expression_clusts = np.log(expression_clusts+1)

        df = pd.DataFrame(
            data=expression_clusts,
            columns=genes,
            index=clusts
        )
        df = df.transpose()

        # Run GSVA
        try:
            df_gsva = run_gsva.run_GSVA(df, gene_set_to_genes)

            clust_to_scores = {
                int(clust[1:]): df_gsva[clust]
                for clust in df_gsva.columns
            }
            clust_scores = np.array([
                clust_to_scores[clust]
                for clust in all_clusts
            ], dtype=np.float32)
            #cell_scores = []
            #for cell, clust in zip(cells, clusters):
            #    cell_scores.append(clust_to_scores[clust])
            #cell_scores = np.array(cell_scores, dtype=np.float32)

            with h5py.File(h5_f, 'r+') as f:
                score_key = '{}_{}_gsva'.format(tumor, collection_name)
                set_key = '{}_{}_gene_set_name'.format(tumor, collection_name)

                ################## TODO REMOVE ###############################################
                try:
                    del f['per_tumor/{}/{}_gsva'.format(tumor, collection_name)]
                except KeyError:
                    pass
                try:
                    del f['per_tumor/{}/{}_gene_set_name'.format(tumor, collection_name)]
                except KeyError:
                    pass
                ##############################################################################

                try:
                    del f['per_tumor/{}/leiden_res_4/{}_gsva'.format(tumor, collection_name)]
                except KeyError:
                    pass
                try:
                    del f['per_tumor/{}/leiden_res_4/{}_gene_set_name'.format(tumor, collection_name)]
                except KeyError:
                    pass

                f['per_tumor/{}/leiden_res_4'.format(tumor)].create_dataset(
                    '{}_gsva'.format(collection_name),
                    data=clust_scores,
                    compression='gzip'
                )
                f['per_tumor/{}/leiden_res_4'.format(tumor)].create_dataset(
                    '{}_gene_set_name'.format(collection_name),
                    data=np.array([
                        x.encode('utf-8')
                        for x in df_gsva.index
                    ]),
                    compression='gzip'
                )
        except Exception as e:
            print("Error running GSVA on tumor {}...".format(tumor))
            print(e)

if __name__ == '__main__':
    main()

