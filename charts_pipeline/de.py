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

N_GENES = 50

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-w", "--overwrite", action='store_true', help="Overwrite data in the HDF5 file if there's a dataset already present")
    parser.add_option("-r", "--resolution", help="Clustering resolution")
    (options, args) = parser.parse_args()

    h5_f = args[0]
    overwrite = options.overwrite
    res = options.resolution

    the_tumors = set()
    with h5py.File(h5_f, 'r') as f:
        the_tumors = f['per_tumor'].keys()
        the_tumors = sorted(the_tumors)
    print(the_tumors)

    for tumor in the_tumors:

        # Delete the tumor's data if we're overwriting
        with h5py.File(h5_f, 'r+') as f:
            if overwrite:
                try:
                    del f['de/leiden_res_{}/{}'.format(res, tumor)]
                except Exception as e:
                    pass

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
            clusters = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res)][:]
            expression = f['per_tumor/{}/log1_tpm'.format(tumor)][:]

        # Convert expression back to TPM
        # for doing fold-change assessment
        #expression = np.exp(expression)+1

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

        # Map each gene to its index
        gene_to_index = {
            gene: index
            for index, gene in enumerate(genes)
        }

        # Run differential expression
        clust_to_de_df = run_de(ad)
        with h5py.File(h5_f, 'r+') as f:
            try:
                f.create_group('de/leiden_res_{}'.format(res))
            except:
                pass
            for clust, df in clust_to_de_df.items():
                #genes = np.array([
                #    x.encode('utf-8')
                #    for x in df.index
                #])

                # Get indices of most significant genes
                gene_indices = np.array([
                    gene_to_index[gene]
                    for gene in df.index
                ])
                if tumor not in f['de/leiden_res_{}'.format(res)].keys() \
                    or str(clust) not in f['de/leiden_res_{}/{}'.format(res, tumor)].keys():
                    print("Writing data for tumor {} cluster {}".format(tumor, clust))
                    try:
                        f['de/leiden_res_{}'.format(res)].create_group(tumor)
                    except Exception as e:
                        pass
                    f['de/leiden_res_{}/{}'.format(res, tumor)].create_group(clust)
                    f['de/leiden_res_{}/{}/{}'.format(res, tumor, clust)].create_dataset('gene_index', data=gene_indices, compression="gzip")
                    f['de/leiden_res_{}/{}/{}'.format(res, tumor, clust)].create_dataset('log_fc', data=np.array(df['log_fc']), compression="gzip")
                    #f['de/leiden_res_{}/{}'.format(res, tum_clust)].create_dataset('pval', data=np.array(df['pval']), compression="gzip")
                else:
                    print("Skipping adding tumor {}, cluster {} to database. Already present.".format(tumor, cluster))

def run_de(ad):
    try:
    #if True:
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
                    ad.uns['rank_genes_groups']['pvals'][clust]
                )
            }
            gene_to_logfold = {
                gene: lf
                for gene, lf in zip(
                    ad.uns['rank_genes_groups']['names'][clust],
                    ad.uns['rank_genes_groups']['logfoldchanges'][clust]
                )
            }
            #de_genes = [
            #    gene
            #    for gene in gene_to_pval.keys()
            #    if gene_to_pval[gene] < 0.05
            #]

            # Sort genes by their p-value
            de_genes = sorted(gene_to_pval.keys(), key=lambda x: gene_to_pval[x])
            de_genes = de_genes[:N_GENES]
            

            print('Top ten genes:')
            print([(g, gene_to_pval[g]) for g in de_genes[:10]])
            #print('Selected genes: ', de_genes[:100])

            df = pd.DataFrame(
                data={
                    'pval': [
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
                'pval': [], 
                'log_fc': []
            }, index=[])
            for clust in sorted(set(ad.obs['cluster']))
        }

if __name__ == '__main__':
    main()

