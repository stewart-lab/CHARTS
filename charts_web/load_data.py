import pandas as pd
import h5py
import numpy as np
import json
from os.path import join
import pkg_resources as pr

resource_package = __name__


CONFIG_F = pr.resource_filename(
    resource_package,
    'config.json'
)

with open(CONFIG_F, 'r') as f:
    config = json.load(f)
    LOC = config['charts_db']
    DB_LOC = join(LOC, 'charts.h5')
    META_LOC = join(LOC, 'tumor_metadata.json')
    SUMMARY_LOC = join(LOC, 'data_summary.tsv')


def get_tumor_meta(tumor):
    with open(META_LOC, 'r') as f:
        j = json.load(f)
    return j[tumor]


def get_all_tumors():
    with h5py.File(DB_LOC, 'r') as f:
        tums = list(f['per_tumor'].keys())
        return tums


def load_tumor_gene_names(tumor):
    with h5py.File(DB_LOC, 'r') as f:
        return set([
            str(x)[2:-1]
            for x in f['{}_gene_name'.format(tumor)][:]
        ])

def load_genes_sorted(tumor):
    with h5py.File(DB_LOC, 'r') as f:
        return sorted([
            str(x)[2:-1]
            for x in f['per_tumor/{}/gene_name'.format(tumor)][:]
        ])


def hover_texts(tumor, res, cells):
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]
        clusters = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res)][:]
        # Map each cell to its cluster
        cell_to_clust = {
            cell: clust
            for cell, clust in zip(cells, clusters)
        }
        # Get each cluster's cell type 
        texts_per_clust = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/leiden_res_{}/cell_type_hover_texts'.format(tumor, res)]
        ]
        # Map each cell to its cell type
        texts = [
            texts_per_clust[cell_to_clust[cell]]
            for cell in cells
        ]
    df = pd.DataFrame(
        data={'hover_texts': texts},
        index=cells
    )
    df = df.loc[cells]
    return list(df['hover_texts'])     

def cell_type_probability_columns(tumor, res, min_prob=None):
    with h5py.File(DB_LOC, 'r') as f:
        cols = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/leiden_res_{}/cell_type_probability_columns'.format(tumor, res)]
        ]
        probs = f['per_tumor/{}/leiden_res_{}/cell_type_probability'.format(tumor, res)][:]
        mins = np.min(probs, axis=0) 
        df = pd.DataFrame(
            data={'probability': mins},
            index=cols
        )
        return set(df.loc[df['probability'] > min_prob].index)

def hallmark_gene_sets(tumor, res):
    with h5py.File(DB_LOC, 'r') as f:
        cols = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/leiden_res_{}/hallmark_gene_set_name'.format(tumor, res)]
        ]
    return sorted(set(cols))


def cancersea_gene_sets(tumor, res):
    with h5py.File(DB_LOC, 'r') as f:
        cols = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/leiden_res_{}/cancersea_gene_set_name'.format(tumor, res)]
        ]
    return sorted(set(cols))


def load_tumor_cell_type_classifications(tumor, res):
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]
        clusters = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res)][:]
        # Map each cell to its cluster
        cell_to_clust = {
            cell: clust
            for cell, clust in zip(cells, clusters)
        }
        # Get each cluster's cell type 
        cell_types_per_clust = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/leiden_res_{}/predicted_cell_type'.format(tumor, res)]
        ]
        # Map each cell to its cell type
        cell_types = [
            cell_types_per_clust[cell_to_clust[cell]]
            for cell in cells
        ]
    df = pd.DataFrame(
        data={
            'cell': cells,
            'color_by': cell_types
        }
    )
    df = df.set_index('cell')
    return df

def load_tumor_clusters(tumor, res):
    with h5py.File(DB_LOC, 'r') as f:
        return set([
            str(x)
            for x in f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res)][:]
        ])

def load_tumor_clusters_for_cells(tumor, res):
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]
        clusts = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res)][:]
    df = pd.DataFrame(
        data={
            'cell': cells,
            'color_by': clusts
        }
    )
    df = df.set_index('cell')
    return df


def load_tumor_gene(tumor, gene):
    gene = gene.lower().strip()
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]
        gene_names = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/gene_name'.format(tumor)][:]
        ]
        gene_name_to_index = {
            gene_name.lower(): index
            for index, gene_name in enumerate(gene_names)
        }
        if gene in gene_name_to_index:
            index = gene_name_to_index[gene]
            expressions = np.array(f['per_tumor/{}/log1_tpm'.format(tumor)][:,index])
        else:
            expressions = np.zeros(len(cells))
    df = pd.DataFrame(
        data={'color_by': expressions},
        index=cells
    )
    return df


def load_tumor_gene_w_cluster(tumor, res, gene):
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]
        clusts = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res)][:]
        gene_names = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/gene_name'.format(tumor)][:]
        ]
        gene_name_to_index = {
            gene_name: index
            for index, gene_name in enumerate(gene_names)
        }
        if gene in gene_name_to_index:
            index = gene_name_to_index[gene]
            expressions = np.array(f['per_tumor/{}/log1_tpm'.format(tumor)][:,index])
        else:
            expressions = np.zeros(len(cells))
    df = pd.DataFrame(
        data={'Expression log(TPM+1)': expressions, 'Cluster': clusts},
        index=cells
    )
    return df


def load_tumor_cell_type_probabilities(tumor, res, cell_type):
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]

        # Retrieve cell type probabilities for each cluster
        cell_types = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/leiden_res_{}/cell_type_probability_columns'.format(tumor, res)][:]
        ]
        cell_type_to_index = {
            cell_type: index
            for index, cell_type in enumerate(cell_types)
        }
        index = cell_type_to_index[cell_type]
        probs_per_clust = np.array(f['per_tumor/{}/leiden_res_{}/cell_type_probability'.format(tumor, res)][:,index])
      
        # Map cluster probabilities to each cell
        clusters = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res)][:]
        cell_to_clust = {
            cell: clust
            for cell, clust in zip(cells, clusters)
        }
        probabilities = [
            probs_per_clust[cell_to_clust[cell]]
            for cell in cells
        ]

    df = pd.DataFrame(
        data={'color_by': probabilities},
        index=cells
    )
    return df


def load_tumor_hallmark_enrichment(tumor, res, gene_set):
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]
        clusters = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res)][:]

        # Map each cell to its cluster
        cell_to_clust = {
            cell: clust
            for cell, clust in zip(cells, clusters)
        }

        # Get gene set enrichment scores
        gene_sets = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/leiden_res_{}/hallmark_gene_set_name'.format(tumor, res)][:]
        ]
        gene_set_to_index = {
            gene_set: index
            for index, gene_set in enumerate(gene_sets)
        }
        index = gene_set_to_index[gene_set]


        # Map scores to clusters
        clust_scores = np.array(f['per_tumor/{}/leiden_res_{}/hallmark_gsva'.format(tumor, res)][:,index])
        scores = np.array([
            clust_scores[cell_to_clust[cell]]
            for cell in cells
        ])

    df = pd.DataFrame(
        data={'color_by': scores},
        index=cells
    )
    return df


def load_tumor_cancersea_enrichment(tumor, res, gene_set):
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]
        clusters = f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, res)][:]

        # Map each cell to its cluster
        cell_to_clust = {
            cell: clust
            for cell, clust in zip(cells, clusters)
        }

        # Get gene set enrichment scores
        gene_sets = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/leiden_res_{}/cancersea_gene_set_name'.format(tumor, res)][:]
        ]
        gene_set_to_index = {
            gene_set: index
            for index, gene_set in enumerate(gene_sets)
        }
        index = gene_set_to_index[gene_set]

        
        # Map scores to clusters
        clust_scores = np.array(f['per_tumor/{}/leiden_res_{}/cancersea_gsva'.format(tumor, res)][:,index])
        scores = np.array([
            clust_scores[cell_to_clust[cell]]
            for cell in cells
        ])

    df = pd.DataFrame(
        data={'color_by': scores},
        index=cells
    )
    return df


def load_malignancy_score(tumor):
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]
        scores = np.array(f['per_tumor/{}/malignancy_score'.format(tumor)][:])
        scores *= -1
    df = pd.DataFrame(
        data={'color_by': scores},
        index=cells
    )
    return df


def load_color_by_real_value_mult_tumors(
        tum_1,
        tum_2,
        cells,
        data_suffix,
        col_suffix,
        feat
    ):
    """
    Load the dataframe used to color each point on a dimension-reduction
    plot using real-valued features such as gene expression or cell type
    probability.
    """
    with h5py.File(DB_LOC, 'r') as f:
        # Load data for first tumor
        cells_1 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_1)][:]
        ]
        cols_1 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/{}'.format(tum_1, col_suffix)][:]
        ]
        col_to_index_1 = {
            col: index
            for index, col in enumerate(cols_1)
        }
        if feat in col_to_index_1:
            index = col_to_index_1[feat]
            key = 'per_tumor/{}/{}'.format(tum_1, data_suffix)
            vals_1 = np.array(f[key][:,index])
        else:
            vals_1 = np.zeros(len(cells_1))
        df_1 = pd.DataFrame(
            data={'color_by': vals_1},
            index=cells_1
        )
        # Load data for second tumor
        cells_2 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_2)][:]
        ]
        cols_2 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/{}'.format(tum_2, col_suffix)][:]
        ]
        col_to_index_2 = {
            col: index
            for index, col in enumerate(cols_2)
        }
        if feat in col_to_index_2:
            index = col_to_index_2[feat]
            vals_2 = np.array(f['per_tumor/{}/{}'.format(tum_2, data_suffix)][:,index])
        else:
            vals_2 = np.zeros(len(cells_2))
        df_2 = pd.DataFrame(
            data={'color_by': vals_2},
            index=cells_2
        )
        # Join the two dataframes
        df = pd.concat([df_1, df_2])
        df = df.loc[cells]
        return df

def load_malignancy_score_mult_tumors(tum_1, tum_2, cells):
    with h5py.File(DB_LOC, 'r') as f:
        # Load data for first tumor
        cells_1 = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tum_1)][:]
        ]
        key = 'per_tumor/{}/malignancy_score'.format(tum_1)
        vals_1 = np.array(f[key][:])
        df_1 = pd.DataFrame(
            data={'color_by': vals_1},
            index=cells_1
        )
        # Load data for second tumor
        cells_2 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_2)][:]
        ]
        key = 'per_tumor/{}/malignancy_score'.format(tum_2)
        vals_2 = np.array(f[key][:])
        df_2 = pd.DataFrame(
            data={'color_by': vals_2},
            index=cells_2
        )
        # Join the two dataframes
        df = pd.concat([df_1, df_2])
        df = df.loc[cells]
        return df


def load_cell_type_classifications_mult_tumors(
        tum_1,
        tum_2,
        cells
    ):
    with h5py.File(DB_LOC, 'r') as f:
        cells_1 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_1)][:]
        ]
        cell_types_1 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/predicted_cell_type'.format(tum_1)][:]
        ]
        df_1 = pd.DataFrame(
            data={'color_by': cell_types_1},
            index=cells_1
        )

        cells_2 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_2)][:]
        ]
        cell_types_2 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/predicted_cell_type'.format(tum_2)][:]
        ]
        df_2 = pd.DataFrame(
            data={'color_by': cell_types_2},
            index=cells_2
        )

        df = pd.concat([df_1, df_2])
        df = df.loc[cells]
        return df


def load_hover_texts_mult_tumors(
        tum_1,
        tum_2,
        cells
    ):
    with h5py.File(DB_LOC, 'r') as f:
        cells_1 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_1)][:]
        ]
        texts_1 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell_type_hover_texts'.format(tum_1)]
        ]
        df_1 = pd.DataFrame(
            data={'hover_texts': texts_1},
            index=cells_1
        )

        cells_2 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_2)][:]
        ]
        texts_2 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell_type_hover_texts'.format(tum_2)]
        ]
        df_2 = pd.DataFrame(
            data={'hover_texts': texts_2},
            index=cells_2
        )

        df = pd.concat([df_1, df_2])
        df = df.loc[cells]
        return list(df['hover_texts'])


def load_clusters_mult_tumors(
        tum_1,
        tum_2,
        cells
    ):
    with h5py.File(DB_LOC, 'r') as f:
        cells_1 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_1)][:]
        ]
        clusters_1 = [
            '{}_{}'.format(tum_1, x)
            for x in f['per_tumor/{}/cluster'.format(tum_1)][:]
        ]
        df_1 = pd.DataFrame(
            data={'color_by': clusters_1},
            index=cells_1
        )

        cells_2 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_2)][:]
        ]
        clusters_2 = [
            '{}_{}'.format(tum_2, x)
            for x in f['per_tumor/{}/cluster'.format(tum_2)][:]
        ]
        df_2 = pd.DataFrame(
            data={'color_by': clusters_2},
            index=cells_2
        )

        df = pd.concat([df_1, df_2])
        df = df.loc[cells]
        return df


def load_tumors_for_cells_mult_tumors(tum_1, tum_2, cells):
    with h5py.File(DB_LOC, 'r') as f:
        cells_1 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_1)][:]
        ]
        tumors_1 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/tumor'.format(tum_1)][:]
        ]
        df_1 = pd.DataFrame(
            data={'color_by': tumors_1},
            index=cells_1
        )

        cells_2 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tum_2)][:]
        ]
        tumors_2 = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/tumor'.format(tum_2)][:]
        ]
        df_2 = pd.DataFrame(
            data={'color_by': tumors_2},
            index=cells_2
        )
        df = pd.concat([df_1, df_2])
        df = df.loc[cells]
        return df


def load_tumor_phate(tumor, num_dims):
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]
        X_umap = f['per_tumor/{}/phate_{}'.format(tumor, num_dims)][:]
    df = pd.DataFrame(
        data=X_umap,
        index=cells,
        columns=['PHATE {}'.format(x+1) for x in range(num_dims)]
    )
    return df


def load_tumor_umap(tumor, num_dims):
    with h5py.File(DB_LOC, 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['per_tumor/{}/cell'.format(tumor)][:]
        ]
        X_umap = f['per_tumor/{}/umap_{}'.format(tumor, num_dims)][:]
    df = pd.DataFrame(
        data=X_umap,
        index=cells,
        columns=['UMAP {}'.format(x+1) for x in range(num_dims)]
    )
    return df


def load_gsva_compare_cluster(collection_name):
    with h5py.File(DB_LOC, 'r') as f:
        cols = [
            ' '.join(str(x)[2:-1].split('_')[1:]).lower()
            for x in f['gsva_compare_cluster/{}/gene_set_name'.format(collection_name)][:]
        ]
        clusts = [
            str(x)[2:-1]
            for x in f['gsva_compare_cluster/{}/cluster'.format(collection_name)][:]
        ]
        scores = f['gsva_compare_cluster/{}/matrix'.format(collection_name)][:]
    
    df = pd.DataFrame(
        data=scores,
        index=clusts,
        columns=cols
    )
    return df


def load_de(res, tumor, cluster):
    with h5py.File(DB_LOC, 'r') as f:
        all_genes = np.array([
            str(x)[2:-1]
            for x in f['per_tumor/{}/gene_name'.format(tumor)][:]
        ])
        de_indices = f['de/leiden_res_{}/{}/{}/gene_index'.format(res, tumor, cluster)][:]
        genes = all_genes[de_indices]
        log_fc = f['de/leiden_res_{}/{}/{}/log_fc'.format(res, tumor, cluster)][:]
        sig_rank = list(np.arange(1,len(de_indices)+1))
    return genes, log_fc, sig_rank


def dataset_summary():
    return pd.read_csv(SUMMARY_LOC, sep='\t')
