import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

import load_data

def build_cell_type_probability_dropdown(tumor, html_id): 
    cell_types = load_data.cell_type_probability_columns(tumor, min_prob=0.00)
    options = [
        {'label': cell_type, 'value': cell_type}
        for cell_type in sorted(cell_types)
    ]
    return dcc.Dropdown(
        options=options,
        value='endothelial cell',
        id=html_id
    )

def build_cell_type_probability_dropdown_mult_tumors(tum_1, tum_2, html_id):
    cell_types_1 = load_data.cell_type_probability_columns(tum_1, min_prob=0.00)
    cell_types_2 = load_data.cell_type_probability_columns(tum_2, min_prob=0.00)
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

def build_hallmark_enrichment_dropdown(tumor, html_id):
    gene_sets = load_data.hallmark_gene_sets(tumor)
    options = [
        {'label': gene_set, 'value': gene_set}
        for gene_set in sorted(gene_sets)
    ]
    return dcc.Dropdown(
        options=options,
        value='HALLMARK_HYPOXIA',
        id=html_id
    )

def build_hallmark_enrichment_dropdown_mult_tumors(tum_1, tum_2, html_id):
    gene_sets_1 = set(load_data.hallmark_gene_sets(tum_1))
    gene_sets_2 = set(load_data.hallmark_gene_sets(tum_2))
    all_gene_sets = gene_sets_1 | gene_sets_2
    options = [
        {'label': gene_set, 'value': gene_set}
        for gene_set in sorted(all_gene_sets)
    ]
    return dcc.Dropdown(
        options=options,
        value='HALLMARK_HYPOXIA',
        id=html_id
    )

def build_cancersea_enrichment_dropdown(tumor, html_id):
    gene_sets = load_data.cancersea_gene_sets(tumor)
    options = [
        {'label': gene_set, 'value': gene_set}
        for gene_set in sorted(gene_sets)
    ]
    return dcc.Dropdown(
        options=options,
        value='Hypoxia',
        id=html_id
    )

def build_cancersea_enrichment_dropdown_mult_tumors(tum_1, tum_2, html_id):
    gene_sets_1 = set(load_data.cancersea_gene_sets(tum_1))
    gene_sets_2 = set(load_data.cancersea_gene_sets(tum_2))
    all_gene_sets = gene_sets_1 | gene_sets_2
    options = [
        {'label': gene_set, 'value': gene_set}
        for gene_set in sorted(all_gene_sets)
    ]
    return dcc.Dropdown(
        options=options,
        value='Hypoxia',
        id=html_id
    )


def build_tumor_dropdown(html_id, width=None):
    options=[
        {'label': 'PJ016 (GBM)', 'value': 'PJ016'},
        {'label': 'PJ017 (GBM)', 'value': 'PJ017'},
        {'label': 'PJ018 (GBM)', 'value': 'PJ018'},
        {'label': 'PJ025 (GBM)', 'value': 'PJ025'},
        {'label': 'PJ030 (LGG)', 'value': 'PJ030'},
        {'label': 'PJ032 (GBM)', 'value': 'PJ032'},
        {'label': 'PJ035 (GBM)', 'value': 'PJ035'},
        {'label': 'PJ048 (GBM)', 'value': 'PJ048'},
        {'label': 'LX653 (LUAD)', 'value': 'LX653'}, 
        {'label': 'LX661 (LUAD)', 'value': 'LX661'},
        {'label': 'LX675 (LUAD)', 'value': 'LX675'},
        {'label': 'LX676 (LUAD)', 'value': 'LX676'},
        {'label': 'LX679 (LUAD)', 'value': 'LX679'},
        {'label': 'LX680 (LUAD)', 'value': 'LX680'},
        {'label': 'LX682 (LUAD)', 'value': 'LX682'},
        {'label': 'LX684 (LUAD)', 'value': 'LX684'},
        {'label': 'GSE146026.1 (OV)', 'value': 'GSE146026.1'}, 
        {'label': 'GSE146026.2 (OV)', 'value': 'GSE146026.2'},
        {'label': 'GSE146026.3 (OV)', 'value': 'GSE146026.3'},
        {'label': 'GSE146026.4 (OV)', 'value': 'GSE146026.4'},
        {'label': 'GSE146026.5 (OV)', 'value': 'GSE146026.5'},
        {'label': 'GSE146026.6 (OV)', 'value': 'GSE146026.6'},
        {'label': 'GSE72056.58 (SKCM)', 'value': 'GSE72056.58'},
        {'label': 'GSE72056.59 (SKCM)', 'value': 'GSE72056.59'},
        {'label': 'GSE72056.60 (SKCM)', 'value': 'GSE72056.60'},
        {'label': 'GSE72056.65 (SKCM)', 'value': 'GSE72056.65'},
        {'label': 'GSE72056.67 (SKCM)', 'value': 'GSE72056.67'},
        {'label': 'GSE72056.71 (SKCM)', 'value': 'GSE72056.71'},
        {'label': 'GSE72056.72 (SKCM)', 'value': 'GSE72056.72'},
        {'label': 'GSE72056.74 (SKCM)', 'value': 'GSE72056.74'},
        {'label': 'GSE72056.75 (SKCM)', 'value': 'GSE72056.75'},
        {'label': 'GSE72056.78 (SKCM)', 'value': 'GSE72056.78'},
        {'label': 'GSE72056.79 (SKCM)', 'value': 'GSE72056.79'},
        {'label': 'GSE72056.80 (SKCM)', 'value': 'GSE72056.80'},
        {'label': 'GSE72056.81 (SKCM)', 'value': 'GSE72056.81'},
        {'label': 'GSE72056.82 (SKCM)', 'value': 'GSE72056.82'},
        {'label': 'GSE72056.84 (SKCM)', 'value': 'GSE72056.84'},
        {'label': 'GSE72056.88 (SKCM)', 'value': 'GSE72056.88'},
        {'label': 'GSE72056.89 (SKCM)', 'value': 'GSE72056.89'},
        {'label': 'GSE72056.94 (SKCM)', 'value': 'GSE72056.94'},
        {'label': 'GSE103322.5 (HNSC)', 'value': 'GSE103322.5'},
        {'label': 'GSE103322.6 (HNSC)', 'value': 'GSE103322.6'},
        {'label': 'GSE103322.7 (HNSC)', 'value': 'GSE103322.7'},
        {'label': 'GSE103322.8 (HNSC)', 'value': 'GSE103322.8'},
        {'label': 'GSE103322.10 (HNSC)', 'value': 'GSE103322.10'},
        {'label': 'GSE103322.12 (HNSC)', 'value': 'GSE103322.12'},
        {'label': 'GSE103322.13 (HNSC)', 'value': 'GSE103322.13'},
        {'label': 'GSE103322.16 (HNSC)', 'value': 'GSE103322.16'},
        {'label': 'GSE103322.17 (HNSC)', 'value': 'GSE103322.17'},
        {'label': 'GSE103322.18 (HNSC)', 'value': 'GSE103322.18'},
        {'label': 'GSE103322.20 (HNSC)', 'value': 'GSE103322.20'},
        {'label': 'GSE103322.22 (HNSC)', 'value': 'GSE103322.22'},
        {'label': 'GSE103322.23 (HNSC)', 'value': 'GSE103322.23'},
        {'label': 'GSE103322.24 (HNSC)', 'value': 'GSE103322.24'},
        {'label': 'GSE103322.25 (HNSC)', 'value': 'GSE103322.25'},
        {'label': 'GSE103322.26 (HNSC)', 'value': 'GSE103322.26'},
        {'label': 'GSE103322.28 (HNSC)', 'value': 'GSE103322.28'},
        {'label': 'GSE111672.A (PAAD)', 'value': 'GSE111672.A'},
        {'label': 'GSE111672.B (PAAD)', 'value': 'GSE111672.B'}
#        {'label': 'PJ016 + PJ018', 'value': 'PJ016&PJ018'},
#        {'label': 'PJ016 + PJ017', 'value': 'PJ016&PJ017'},
#        {'label': 'PJ016 + PJ025', 'value': 'PJ016&PJ025'},
#        {'label': 'PJ017 + PJ018', 'value': 'PJ017&PJ018'},
#        {'label': 'PJ017 + PJ025', 'value': 'PJ017&PJ025'},
#        {'label': 'PJ018 + PJ025', 'value': 'PJ018&PJ025'},
#        {'label': 'PJ016 + LX653', 'value': 'PJ016&LX653'},
#        {'label': 'PJ016 + LX676', 'value': 'PJ016&LX676'},
#        {'label': 'PJ016 + LX682', 'value': 'PJ016&LX682'},
#        {'label': 'PJ016 + GSE146026.5', 'value': 'PJ016&GSE146026.5'},
#        {'label': 'PJ017 + LX653', 'value': 'PJ017&LX653'},
#        {'label': 'PJ017 + LX676', 'value': 'PJ017&LX676'},
#        {'label': 'PJ017 + LX682', 'value': 'PJ017&LX682'},
#        {'label': 'PJ017 + GSE146026.5', 'value': 'PJ017&GSE146026.5'},
#        {'label': 'PJ018 + LX653', 'value': 'PJ018&LX653'},
#        {'label': 'PJ018 + LX676', 'value': 'PJ018&LX676'},
#        {'label': 'PJ018 + LX682', 'value': 'PJ018&LX682'},
#        {'label': 'PJ018 + GSE146026.5', 'value': 'PJ018&GSE146026.5'},
#        {'label': 'PJ025 + LX653', 'value': 'PJ025&LX653'},
#        {'label': 'PJ025 + LX676', 'value': 'PJ025&LX676'},
#        {'label': 'PJ025 + LX682', 'value': 'PJ025&LX682'},
#        {'label': 'PJ025 + GSE146026.5', 'value': 'PJ025&GSE146026.5'},
#        {'label': 'LX653 + LX676', 'value': 'LX653&LX676'},
#        {'label': 'LX653 + LX682', 'value': 'LX653&LX682'},
#        {'label': 'LX653 + GSE146026.5', 'value': 'LX653&GSE146026.5'},
#        {'label': 'LX676 + LX682', 'value': 'LX676&LX682'},
#        {'label': 'LX676 + GSE146026.5', 'value': 'LX676&GSE146026.5'}
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
