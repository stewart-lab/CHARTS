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
from collections import Counter

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
        join(root, 'GSE72056_melanoma_single_cell_revised_v2.txt'), 
        sep='\t', 
        header=None,
        low_memory=False,
        index_col=0
    )
    df = df.transpose()
    df = df.set_index('Cell')

    print(list(df.columns)[:10])
    keep_cols = df.columns[3:]
    print(len(keep_cols))
    genes = [
        x.encode('utf-8')
        for x in keep_cols
    ]


    gene_to_count = Counter(keep_cols)
    

    # Sum together genes with same symbol
    #df = df.T.groupby([i.split('.')[0] for i in df.T.index.values]).sum().T
    #print(df)

    dup_genes = set([
        gene
        for gene, count in gene_to_count.items()
        if count > 1
    ])
    print("Duplicate gene names:")
    print(dup_genes)
    for dup_gene in dup_genes:
        s = np.sum(
            np.array(df[dup_gene], dtype=np.float32),
            axis=1
        )
        df = df.drop([dup_gene], axis=1)
        df[dup_gene] = list(s)
    print(df.shape)
    print(df)

    print(np.array(df[keep_cols], dtype=np.float32))

    for sample in sorted(set(df['tumor'])):
        cells = df.loc[df['tumor'] == sample].index
        tumor = 'GSE72056.{}'.format(sample)
        print("{} cells in sample {}".format(len(cells), tumor))
        curr_df = df[keep_cols].loc[cells]

        # Convert to TPM
        X = np.array(curr_df, dtype=np.float32)
        X = np.power(2, X)-1
        # Compute log1(TPM)
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
