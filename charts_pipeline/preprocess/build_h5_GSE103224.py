#
#   Parse the counts data from the Gene Expression Omnibus for
#   accession GSE103224. Build an H5 file storing the counts.
#

from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py 
from os.path import join, exists
from anndata import AnnData
import scanpy as sc

TUMOR_TO_FILE = {
    'PJ016': 'GSM2758471_PJ016.filtered.matrix.txt.gz',
    'PJ018': 'GSM2758473_PJ018.filtered.matrix.txt.gz',
    'PJ030': 'GSM2758475_PJ030.filtered.matrix.txt.gz',
    'PJ035': 'GSM2758477_PJ035.filtered.matrix.txt.gz',
    'PJ017': 'GSM2758472_PJ017.filtered.matrix.txt.gz',
    'PJ025': 'GSM2758474_PJ025.filtered.matrix.txt.gz',
    'PJ032': 'GSM2758476_PJ032.filtered.matrix.txt.gz',
    'PJ048': 'GSM2940098_PJ048.filtered.matrix.txt.gz'
}

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    root = args[0]
    out_f = args[1]

    print('Finding genes common to all datasets...')
    the_genes = None
    counts = None
    cells = []
    tumors = []

    # Create HDF5 file if it doesn't exist
    if not exists(out_f):
        with h5py.File(out_f, 'w') as f:
            pass

    with h5py.File(out_f, 'r+') as f:
        for tumor, tumor_fname in TUMOR_TO_FILE.items():
            tumor_f = join(root, tumor_fname)
            print('Loading {}...'.format(tumor_f))
            df = pd.read_csv(
                tumor_f, 
                sep='\t', 
                header=None, 
                index_col=0
            )

            # Extract the genes
            genes = list(df.index)
            if the_genes is None:
                the_genes = list(genes)
                gene_names = list(df.iloc[:,0])
            assert frozenset(genes) == frozenset(the_genes)

            # Drop the gene-names column
            keep_cols = df.columns[1:]
            df = df[keep_cols]

            # Re-order rows according to the global gene 
            # ordering
            df = df.loc[the_genes]

            # Cast to a numpy array
            curr_counts = np.array(df, dtype=np.float32).T

            ad = AnnData(
                X=curr_counts
            )
            sc.pp.normalize_total(ad, target_sum=1e6)
            sc.pp.log1p(ad)
            curr_counts = ad.X

            cells = [
                '{}_{}'.format(tumor, i).encode('utf-8')
                for i in np.arange(curr_counts.shape[0])
            ]
            tumors = [
                tumor.encode('utf-8')
                for i in np.arange(curr_counts.shape[0])
            ]
            store_genes = [
                x.encode('utf-8')
                for x in the_genes
            ]
            store_gene_names = [
                x.encode('utf-8')
                for x in gene_names
            ]
            try:
                del f['per_tumor/{}'.format(tumor)]
            except KeyError:
                pass
            try:
                f.create_group('per_tumor')
            except:
                pass
            f['per_tumor'].create_group(tumor)
            f['per_tumor'][tumor].create_dataset(
                'log1_tpm', 
                data=curr_counts, 
                compression="gzip"
            )
            f['per_tumor'][tumor].create_dataset(
                'cell', 
                data=np.array(cells),
                compression="gzip"
            )
            f['per_tumor'][tumor].create_dataset(
                'tumor', 
                data=np.array(tumors),
                compression="gzip"
            )
            f['per_tumor'][tumor].create_dataset(
                'gene_id', 
                data=np.array(store_genes),
                compression="gzip"
            )
            f['per_tumor'][tumor].create_dataset(
                'gene_name', 
                data=np.array(store_gene_names),
                compression="gzip"
            )

if __name__ == '__main__':
    main()
