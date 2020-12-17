import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

import load_data

def build_cell_type_probability_dropdown(tumor, res, html_id): 
    cell_types = load_data.cell_type_probability_columns(tumor, res, min_prob=0.00)
    options = [
        {'label': cell_type, 'value': cell_type}
        for cell_type in sorted(cell_types)
    ]
    return dcc.Dropdown(
        options=options,
        value='endothelial cell',
        id=html_id
    )


def build_cluster_res_dropdown(idd):
    return dcc.Dropdown(
        options=[
            {'label': '2', 'value': 2},
            {'label': '4', 'value': 4}
        ],
        value=4,
        id=idd
    )


def build_cell_type_probability_dropdown_mult_tumors(tum_1, tum_2, html_id):
    cell_types_1 = load_data.cell_type_probability_columns(tum_1, res, min_prob=0.00)
    cell_types_2 = load_data.cell_type_probability_columns(tum_2, res, min_prob=0.00)
    all_cell_types = cell_types_1 | cell_types_2
    options = [
        {'label': cell_type, 'value': cell_type}
        for cell_type in sorted(all_cell_types)
    ]
    return dcc.Dropdown(
        options=options,
        value='endothelial cell',
        id=html_id
    )

def build_hallmark_enrichment_dropdown(tumor, res, html_id):
    gene_sets = load_data.hallmark_gene_sets(tumor, res)
    options = [
        {'label': gene_set, 'value': gene_set}
        for gene_set in sorted(gene_sets)
    ]
    return dcc.Dropdown(
        options=options,
        value='HALLMARK_HYPOXIA',
        id=html_id
    )

#def build_hallmark_enrichment_dropdown_mult_tumors(tum_1, tum_2, html_id):
#    gene_sets_1 = set(load_data.hallmark_gene_sets(tum_1))
#    gene_sets_2 = set(load_data.hallmark_gene_sets(tum_2))
#    all_gene_sets = gene_sets_1 | gene_sets_2
#    options = [
#        {'label': gene_set, 'value': gene_set}
#        for gene_set in sorted(all_gene_sets)
#    ]
#    return dcc.Dropdown(
#        options=options,
#        value='HALLMARK_HYPOXIA',
#        id=html_id
#    )

def build_cancersea_enrichment_dropdown(tumor, res, html_id):
    gene_sets = load_data.cancersea_gene_sets(tumor, res)
    options = [
        {'label': gene_set, 'value': gene_set}
        for gene_set in sorted(gene_sets)
    ]
    return dcc.Dropdown(
        options=options,
        value='Hypoxia',
        id=html_id
    )

#def build_cancersea_enrichment_dropdown_mult_tumors(tum_1, tum_2, html_id):
#    gene_sets_1 = set(load_data.cancersea_gene_sets(tum_1))
#    gene_sets_2 = set(load_data.cancersea_gene_sets(tum_2))
#    all_gene_sets = gene_sets_1 | gene_sets_2
#    options = [
#        {'label': gene_set, 'value': gene_set}
#        for gene_set in sorted(all_gene_sets)
#    ]
#    return dcc.Dropdown(
#        options=options,
#        value='Hypoxia',
#        id=html_id
#    )


def build_tumor_dropdown(html_id, width=None):
    all_tumors = load_data.get_all_tumors()
    tum_to_type = {
        tum : load_data.get_tumor_meta(tum)['cancer_type_abbrev']
        for tum in all_tumors
    }
    options = [
        {
            'label': '{} ({})'.format(
                tum, tum_to_type[tum]   
            ), 
            'value': tum
        }
        for tum in sorted(all_tumors)
    ]
    if width:
        style={"width": width}
    else:
        style=None
    return dcc.Dropdown(
        options=options,
        value='GSE72056.89',
        id=html_id,
        style=style
    )
