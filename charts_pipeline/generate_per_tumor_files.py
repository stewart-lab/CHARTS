import h5py
from optparse import OptionParser
import pandas as pd
import subprocess
import numpy as np
from os.path import join

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-l", "--log_dir", help="Directory in which to write log files")
    parser.add_option("-w", "--overwrite", action='store_true', help="Overwrite data in the HDF5 file if there's a dataset already present")
    (options, args) = parser.parse_args()

    h5_f = args[0]
    out_dir = options.out_dir

    res_low = '0.5'
    res_med = '1'
    res_high = '2' 

    the_tumors = set()
    with h5py.File(h5_f, 'r') as f:
        the_tumors = f['per_tumor'].keys()
        the_tumors = sorted(the_tumors)
    print(the_tumors)

    for tumor in the_tumors[:2]:
        # Create tumor's directory
        subprocess.run("mkdir {}".format(join(out_dir, tumor)), shell=True)

        print("Running algorithms on tumor {}".format(tumor))
        with h5py.File(h5_f, 'r') as f:
            cells = [
                str(x)[2:-1]
                for x in f['per_tumor/{}/cell'.format(tumor)][:]
            ]
            expression = f['per_tumor/{}/log1_tpm'.format(tumor)][:]

            umap_2 = f['per_tumor/{}/umap_2'.format(tumor)][:].T
            umap_3 = f['per_tumor/{}/umap_3'.format(tumor)][:].T
            phate_2 = f['per_tumor/{}/phate_2'.format(tumor)][:].T
            phate_3 = f['per_tumor/{}/phate_3'.format(tumor)][:].T
            print(phate_3.shape)
            malig_scores = f['per_tumor/{}/malignancy_score'.format(tumor)][:]
            print(malig_scores.shape)
            clusters_low = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res_low)][:]
            clusters_med = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res_med)][:]
            clusters_high = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res_high)][:]
        

            for res in ['0.5', '1', '2']:
                clusters =  [str(x) for x in range(len(set(f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res)][:])))]
                print("Clusters: ", clusters)
                cell_type_probs = f['per_tumor/{}/leiden_res_{}/cell_type_probability'.format(tumor, res)][:]
                cell_type_cols = [
                    '{} (probability)'.format(str(x)[2:-1])
                    for x in f['per_tumor/{}/leiden_res_{}/cell_type_probability_columns'.format(tumor, res)][:]
                ]
                pred_cell_type = [
                    str(x)[2:-1]
                    for x in f['per_tumor/{}/leiden_res_{}/predicted_cell_type'.format(tumor, res)][:]
                ]
                df_clust = pd.DataFrame(
                    data=cell_type_probs,
                    columns=cell_type_cols,
                    index=clusters
                )
                df_clust['Predicted Cell Type'] = pred_cell_type

                found_cancersea = True
                try:
                    cancersea = f['per_tumor/{}/leiden_res_{}/cancersea_gsva'.format(tumor, res)]
                    cancersea_gene_sets = [
                        str(x)[2:-1]
                        for x in f['per_tumor/{}/leiden_res_{}/cancersea_gene_set_name'.format(tumor, res)]
                    ]
                    df_clust_cancersea = pd.DataFrame(
                        data=cancersea,
                        columns=cancersea_gene_sets,
                        index=clusters
                    )
                    df_clust = df_clust.join(df_clust_cancersea)
                except KeyError:
                    print("No CancerSEA GSVA results found for tumor {}...".format(tumor))
                    found_cancersea = False

                found_hallmark = False
                try:
                    hallmark = f['per_tumor/{}/leiden_res_{}/hallmark_gsva'.format(tumor, res)]
                    hallmark_gene_sets = [
                        str(x)[2:-1]
                        for x in f['per_tumor/{}/leiden_res_{}/hallmark_gene_set_name'.format(tumor, res)]
                    ]
                    df_clust_hallmark = pd.DataFrame(
                        data=hallmark,
                        columns=hallmark_gene_sets,
                        index=clusters
                    )
                    df_clust = df_clust.join(df_clust_hallmark)
                except KeyError:
                    print("No CancerSEA GSVA results found for tumor {}...".format(tumor))
                    found_hallmark = False

                out_f = join(out_dir, tumor, 'per_cluster.res_{}.tsv'.format(res))
                print("Writing cluster information to file {}...".format(out_f))
                df_clust.to_csv(out_f, sep='\t')

            # Create per-cell information table
            df_cell = pd.DataFrame(
                data={
                    "UMAP1 (2D)": umap_2[0],
                    "UMAP2 (2D)": umap_2[1],
                    "UMAP1 (3D)": umap_3[0], 
                    "UMAP2 (3D)": umap_3[1], 
                    "UMAP3 (3D)": umap_3[2], 
                    "PHATE1 (2D)": phate_2[0],
                    "PHATE2 (2D)": phate_2[1],
                    "PHATE1 (3D)": phate_3[0],
                    "PHATE2 (3D)": phate_3[1],
                    "PHATE3 (3D)": phate_3[2],
                    "Cluster (Resolution = 0.5)": clusters_low,
                    "Cluster (Resolution = 1.0)": clusters_med,
                    "Cluster (Resolution = 2.0)": clusters_high,
                    "Malignancy Score (Entropy)": malig_scores
                },
                index=cells
            )
            df_cell.to_csv(join(out_dir, tumor, 'per_cell.tsv'), sep='\t')

            # Create expression matrix table
            expr = f['per_tumor/{}/log1_tpm'.format(tumor)][:]
            genes = [
                str(x)[2:-1]
                for x in f['per_tumor/{}/gene_name'.format(tumor)][:]
            ]
            df_expr = pd.DataFrame(
                data=expr,
                index=cells,
                columns=genes
            )
            df_expr.to_csv(join(out_dir, tumor, 'expression_log1_tpm.tsv'), sep='\t')

            cmd = "tar -C {chdir} -zcf {dirr}.tar.gz ./{tumor}".format(dirr=join(out_dir, tumor), tumor=tumor, chdir=out_dir)
            print(cmd)
            subprocess.run(cmd, shell=True)
            cmd = "rm -r {}".format(join(out_dir, tumor))
            print(cmd)
            subprocess.run(cmd, shell=True)

if __name__ == '__main__':
    main()
