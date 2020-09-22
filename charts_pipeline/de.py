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
    overwrite = options.overwrite
    out_dir = options.out_dir

    the_tumors = set()
    with h5py.File(h5_f, 'r') as f:
        the_tumors = f['per_tumor'].keys()
        the_tumors = sorted(the_tumors)
    print(the_tumors)

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
            clusters = f['per_tumor/{}/cluster'.format(tumor)][:]
            expression = f['per_tumor/{}/log1_tpm'.format(tumor)][:]

        # Convert expression back to TPM
        # for doing fold-change assessment
        expression = np.exp(expression)+1

        # Map each cluster to its cells
        clust_to_cells = defaultdict(lambda: [])
        for clust, cell in zip(clusters, cells):
            clust_to_cells[clust].append(cell)

        # Map each cell to its index
        cell_to_index = {
            cell: index
            for index, cell in enumerate(cells)
        }

        ad = AnnData(
            X=expression,
            var=pd.DataFrame(
                index=genes,
                data={'gene': genes}
            ),
            obs=pd.DataFrame(
                data={
                    'cluster': [str(x) for x in clusters]
                },
                index=cells
            )
        )

        # Run differential expression
        clust_to_de_df = run_de(ad)
        with h5py.File(h5_f, 'r+') as f:
            try:
                f.create_group('de')
            except:
                pass
            for clust, df in clust_to_de_df.items():
                tum_clust = '{}_{}'.format(tumor, clust)
                genes = np.array([
                    x.encode('utf-8')
                    for x in df.index
                ])
                if overwrite:
                    try:
                        del f['de/{}'.format(tumor_clust)]
                    except:
                        pass
                if tum_clust not in f['de'].keys():
                    print("writing data for tumor cluster {}".format(tum_clust))
                    f['de'].create_group(tum_clust)
                    f['de/{}'.format(tum_clust)].create_dataset('gene', data=genes)
                    f['de/{}'.format(tum_clust)].create_dataset('log_fc', data=np.array(df['log_fc']))
                    f['de/{}'.format(tum_clust)].create_dataset('pval_adj', data=np.array(df['pval_adj']))
                else:
                    print("Skipping adding cluster {} to database. Already present.".format(tum_clust))


def run_de(ad):

    try:
        sc.tl.rank_genes_groups(
            ad,
            groupby='cluster',
            method='wilcoxon',
            n_genes=len(ad.var)
        )
        clust_to_de_genes = {}
        for clust in sorted(set(ad.obs['cluster'])):
            gene_to_pval = {
                gene: pval
                for gene, pval in zip(
                    ad.uns['rank_genes_groups']['names'][clust],
                    ad.uns['rank_genes_groups']['pvals_adj'][clust]
                )
            }
            gene_to_logfold = {
                gene: lf
                for gene, lf in zip(
                    ad.uns['rank_genes_groups']['names'][clust],
                    ad.uns['rank_genes_groups']['logfoldchanges'][clust]
                )
            }
            de_genes = [
                gene
                for gene in gene_to_pval.keys()
                if gene_to_pval[gene] < 0.05
            ]
            df = pd.DataFrame(
                data={
                    'pval_adj': [
                        gene_to_pval[g]
                        for g in de_genes
                    ],
                    'log_fc': [
                        gene_to_logfold[g]
                        for g in de_genes
                    ]
                },
                index=de_genes
            )
            clust_to_de_genes[clust] = df
        return clust_to_de_genes
    except ZeroDivisionError:
        return {
            clust: pd.DataFrame({
                'pval_adj': [], 
                'log_fc': []
            }, index=[])
            for clust in sorted(set(ad.obs['cluster']))
        }

if __name__ == '__main__':
    main()

