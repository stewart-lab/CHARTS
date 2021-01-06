#
#   Parse the counts data from the Gene Expression Omnibus for
#   accession GSE103224. Build an H5 file storing the counts.
#

from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py 
from os.path import join
from collections import defaultdict
import scanpy as sc
from anndata import AnnData

DATASET_TO_FILE = {
    'LX653': 'GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv',
    'LX661': 'GSM3516663_MSK_LX661_PRIMARY_TUMOUR_dense.csv',
    'LX675': 'GSM3516665_MSK_LX675_PRIMARY_TUMOUR_dense.csv',
    'LX676': 'GSM3516667_MSK_LX676_PRIMARY_TUMOUR_dense.csv',
    'LX679': 'GSM3516669_MSK_LX679_PRIMARY_TUMOUR_dense.csv', 
    'LX680': 'GSM3516670_MSK_LX680_PRIMARY_TUMOUR_dense.csv',
    'LX682': 'GSM3516672_MSK_LX682_PRIMARY_TUMOUR_dense.csv',
    'LX684': 'GSM3516674_MSK_LX684_PRIMARY_TUMOUR_dense.csv'
}
INTERSECTION = False

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    root = args[0]
    out_f = args[1]

    print('Finding genes common to all datasets...')
    counts = None
    cells = []
    datasets = []
    patients = []
    is_tumors = []
    sym_to_syns = defaultdict(lambda: set())
    the_genes = None

    with h5py.File(out_f, 'r+') as f:
        for dataset, ds_fname in DATASET_TO_FILE.items():
            ds_f = join(root, ds_fname)
            print('Loading {}...'.format(ds_f))
            df = pd.read_csv(
                ds_f,
                index_col=0
            )
            print(df.shape)

            # Cast to a numpy array
            curr_counts = np.array(df, dtype=np.float32)
            ad = AnnData(
                X=curr_counts
            )
            sc.pp.normalize_total(ad, target_sum=1e6)
            sc.pp.log1p(ad)
            curr_counts = ad.X

            cells = [
                '{}_{}'.format(dataset, i).encode('utf-8')
                for i in np.arange(curr_counts.shape[0])
            ]
            datasets = [
                dataset.encode('utf-8')
                for i in np.arange(curr_counts.shape[0])
            ]
            gene_names = [
                x.encode('utf-8')
                for x in df.columns
            ]
            try:
                del f['per_tumor/{}'.format(dataset)]
            except KeyError:
                pass
            try:
                f.create_group('per_tumor')
            except:
                pass
            f['per_tumor'].create_group(dataset)
            f['per_tumor'][dataset].create_dataset(
                'log1_tpm',
                data=curr_counts,
                compression="gzip"
            )
            f['per_tumor'][dataset].create_dataset(
                'cell',
                data=np.array(cells),
                compression="gzip"
            )
            f['per_tumor'][dataset].create_dataset(
                'tumor',
                data=np.array(datasets),
                compression="gzip"
            )
            f['per_tumor'][dataset].create_dataset(
                'gene_name',
                data=np.array(gene_names),
                compression="gzip"
            )

def _standard_gene(gene):
    try:
        return SYN_TO_SYM[gene]
    except KeyError:
        return gene

if __name__ == '__main__':
    main()
