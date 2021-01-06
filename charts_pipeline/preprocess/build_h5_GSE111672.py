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
    'GSE111672.A': 'GSE111672_PDAC-A-indrop-filtered-expMat.txt',
    'GSE111672.B': 'GSE111672_PDAC-B-indrop-filtered-expMat.txt'
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
                sep='\t',
                index_col=0
            )
            df = df.transpose()
            print(df)

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
