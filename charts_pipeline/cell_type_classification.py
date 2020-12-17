import subprocess
import h5py
import sys
from anndata import AnnData
import scanpy as sc
from optparse import OptionParser
import pandas as pd
import os
from os.path import join
import dill
import numpy as np
import json

sys.path.append('../../../CellO')

import CellO

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-w", "--overwrite", action="store_true", help="Overwrite data in database")
    parser.add_option("-r", "--resolution", help="Clustering to classify")
    (options, args) = parser.parse_args()

    h5_f = args[0]
    model_dir = args[1]
    tumor_params_f = args[2]

    with open(tumor_params_f, 'r') as f:
        tumor_params = json.load(f)['per_tumor']
    
    resolution = options.resolution
    overwrite = options.overwrite

    sc.settings.verbosity = 3
    sc.logging.print_versions()

    the_tumors = set()
    with h5py.File(h5_f, 'r') as f:
        the_tumors = f['per_tumor'].keys()
        the_tumors = sorted(the_tumors)
        if not overwrite:
            the_tumors = [
                tumor
                for tumor in the_tumors 
                if 'predicted_cell_type' not in f['per_tumor/{}'.format(tumor)].keys()
            ]
    print(the_tumors)

    for tumor in the_tumors:
        if tumor == 'GSE103322.7':
            continue # TODO REMOVE!!!!
        print("Classifying cell types for tumor {}".format(tumor))
        with h5py.File(h5_f, 'r') as f:
            cells = [
                str(x)[2:-1]
                for x in f['per_tumor/{}/cell'.format(tumor)][:]
            ]
            if 'gene_id' in f['per_tumor/{}'.format(tumor)].keys():
                genes = [
                    str(x)[2:-1]
                    for x in f['per_tumor/{}/gene_id'.format(tumor)][:]
                ]
            else:
                genes = [
                    str(x)[2:-1]
                    for x in f['per_tumor/{}/gene_name'.format(tumor)][:]
                ]
            expression = f['per_tumor/{}/log1_tpm'.format(tumor)][:]
            clusters = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, resolution)][:]

        all_clusts = sorted(set(clusters))

        # Create expression AnnData object
        ad = AnnData(
            X=expression,
            obs=pd.DataFrame(data=cells, columns=['cell'], index=cells),
            var=pd.DataFrame(data=genes, columns=['gene'], index=genes)
        )
        ad.raw = ad

        # Map each cell to its cluster
        cell_to_cluster = {
            cell: clust
            for cell, clust in zip(cells, clusters)
        }

        mod = CellO.retreive_pretrained_model_from_local(ad, model_dir)
        if mod is None:
            mod = CellO.train_model(ad, algo='IR', log_dir=None)
            with open(join(model_dir, '{}_model.dill'.format(tumor)), 'wb') as f:
                dill.dump(mod, f)

        units = 'LOG1_TPM'
        
        # Get the anatomical entities to filter for the current tumor
        if tumor in tumor_params:
            filter_anatomical_terms = tumor_params[tumor]['exclude_anatomical_terms_from_cell_types']
        else:
            filter_anatomical_terms = None

        results_df, finalized_binary_results_df, ms_results_df = CellO.predict(
            ad,
            units,
            mod,
            assay='3_PRIME',
            algo='IR',
            cell_to_cluster=cell_to_cluster,
            remove_anatomical_subterms=filter_anatomical_terms
        )

        results_df.columns = [
            CellO.CELL_ONTOLOGY.id_to_term[x].name
            for x in results_df.columns
        ]
        finalized_binary_results_df.columns = [
            CellO.CELL_ONTOLOGY.id_to_term[x].name
            for x in finalized_binary_results_df.columns
        ]
        ms_results_df['most_specific_cell_type'] = [
            CellO.CELL_ONTOLOGY.id_to_term[x].name
            for x in ms_results_df['most_specific_cell_type']
        ]

        # Map each cluster to its predicted cell type
        clust_to_pred_cell_type = {
            cell_to_cluster[cell]: pred
            for cell, pred in zip(ms_results_df.index, ms_results_df['most_specific_cell_type'])
        }
        cell_types_predicted = np.array([
            clust_to_pred_cell_type[clust].encode('utf-8')
            for clust in all_clusts
        ])
        

        # Map each cluster to its cell type probabilities
        clust_to_cell_type_prob = {
            cell_to_cluster[cell]: np.array(row)
            for cell, (i, row) in zip(cells, results_df.iterrows())
        }
        cell_type_probs = np.array([
            clust_to_cell_type_prob[clust]
            for clust in all_clusts
        ])

        #cell_types_predicted = np.array([
        #    x.encode('utf-8')
        #    for x in ms_results_df['most_specific_cell_type']
        #])
        cell_type_prob_columns = np.array([
            x.encode('utf-8')
            for x in results_df.columns
        ])
        #cell_type_probs = np.array(results_df, dtype=np.float32)
        
        # This is a bit hacky, but we're going to pre-compute the text
        # that will be displayed upon each cell's hover-over event
        hover_texts = []
        for cell in finalized_binary_results_df.index:
            curr_pairs = []

            # TODO clean this up
            probs = results_df.loc[cell]
            cell_type_to_prob = {
                cell_type: prob
                for cell_type, prob in zip(results_df.columns, probs)
            }
            for cell_type, is_pred in zip(finalized_binary_results_df.columns, finalized_binary_results_df.loc[cell]):
                prob = cell_type_to_prob[cell_type]
                if (is_pred == 1 or prob > 0.5) and not cell_type in set(['cell', 'animal cell', 'native cell', 'eukaryotic cell']):
                    curr_pairs.append((
                        cell_type,
                        results_df.loc[cell][cell_type]
                    ))
            curr_pairs = sorted(curr_pairs, key=lambda x: x[1], reverse=True)
            txt = 'Cell Type (Probability)<br />-------------------------<br />'
            txt += '<br />'.join([
                '{} ({:.2f})'.format(x[0], x[1])
                for x in curr_pairs
            ])
            hover_texts.append(txt)

        # Map each cluster to its hover text
        clust_to_hover = {
            cell_to_cluster[cell]: txt
            for cell, txt in zip(cells, hover_texts)
        }
        hover_texts = np.array([
            clust_to_hover[clust].encode('utf-8')
            for clust in all_clusts
        ])

        with h5py.File(h5_f, 'r+') as f:
            try:
                del f['per_tumor/{}/leiden_res_{}/predicted_cell_type'.format(tumor, int(resolution))]
            except KeyError:
                pass
            f.create_dataset(
                'per_tumor/{}/leiden_res_{}/predicted_cell_type'.format(tumor, int(resolution)),
                data=cell_types_predicted,
                compression="gzip"
            )
            try:
                del f['per_tumor/{}/leiden_res_{}/cell_type_probability'.format(tumor, int(resolution))]
            except KeyError:
                pass
            f.create_dataset(
                'per_tumor/{}/leiden_res_{}/cell_type_probability'.format(tumor, int(resolution)),
                data=cell_type_probs,
                compression="gzip"
            )
            try:
                del f['per_tumor/{}/leiden_res_{}/cell_type_probability_columns'.format(tumor, int(resolution))]
            except KeyError:
                pass
            f.create_dataset(
                'per_tumor/{}/leiden_res_{}/cell_type_probability_columns'.format(tumor, int(resolution)),
                data=cell_type_prob_columns,
                compression="gzip"
            )
            try:
                del f['per_tumor/{}/leiden_res_{}/cell_type_hover_texts'.format(tumor, int(resolution))]
            except KeyError:
                pass
            f.create_dataset(
                'per_tumor/{}/leiden_res_{}/cell_type_hover_texts'.format(tumor, int(resolution)),
                data=hover_texts,
                compression="gzip"
            )

if __name__ == '__main__':
    main()
