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
from anndata import AnnData
import scanpy as sc
from collections import defaultdict

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    root = args[0]
    out_f = args[1]

    the_genes = None
    counts = None
    cells = []
    tumors = []

    print("Loading data...")
    df = pd.read_csv(
        join(root, 'GSE118389_tpm_rsem.txt'), 
        sep='\t', 
        low_memory=False,
        index_col=0
    )
    df = df.transpose()
    print(df)

    tumor_to_cells = defaultdict(lambda: [])
    for cell_id in df.index:
        tumor_id = cell_id.split('_')[0]
        tumor_to_cells[tumor_id].append(cell_id)

    print("Found the following tumor ID's:")
    print(tumor_to_cells.keys())

    genes = [
        x.encode('utf-8')
        for x in df.columns
    ]

    for tumor, cells in tumor_to_cells.items():
        curr_df = df.loc[cells]
        tumor = 'GSE118389.{}'.format(tumor)
        print("{} cells in sample {}".format(len(cells), tumor))

        # Compute log1(TPM)
        X = np.array(curr_df, dtype=np.float32)
        X = np.log(X+1)
        print(X.shape)

        cells = [
            x.encode('utf-8')
            for x in curr_df.index
        ]
        tumors = [
            tumor.encode('utf-8')
            for i in np.arange(X.shape[0])
        ]

        print("Writing data to {}...".format(out_f))
        with h5py.File(out_f, 'r+') as f:
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
                data=X,
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
                'gene_name',
                data=np.array(genes),
                compression="gzip"
            )


if __name__ == '__main__':
    main()
