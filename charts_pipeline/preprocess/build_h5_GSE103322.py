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
from collections import Counter, defaultdict

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
        join(root, 'HNSCC_all_data.txt.gz'), 
        sep='\t', 
        low_memory=False,
        index_col=0
    )
    df = df.iloc[5:]
    df = df.transpose()
    #print(df.index)

    # Map each cell to its tumor
    cell_to_tumor = {}
    for cell in df.index:
        if 'HNSCC' in cell:
            tok = cell.split('_')[0]
            if tok == 'HNSCC':
                #tumor = 'unknown'
                continue
            else:
                tumor = int(cell.split('_')[0][5:])
        elif 'HN' in cell:
            tumor = int(cell.split('_')[0][2:])
        cell_to_tumor[cell] = tumor

    # Map each tumor to cells
    tumor_to_cells = defaultdict(lambda: [])
    for cell, tumor in cell_to_tumor.items():
        tumor_to_cells[tumor].append(cell)
    tumor_to_cells = dict(tumor_to_cells)

    print(tumor_to_cells.keys())
    print(df)

    #df = df.transpose()
    #df = df.set_index('Cell')

    genes = [
        x.encode('utf-8')[1:-1]
        for x in df.columns
    ]
    print(genes[:100])


    for tumor, cells in tumor_to_cells.items(): 
        tumor = 'GSE103322.{}'.format(tumor)
        print("{} cells in tumor {}".format(len(cells), tumor))
        curr_df = df.loc[cells]

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
        
        if len(cells) < 50:
            print('Skipping tumor {}. It contains {} cells, which is less than the 50 cutoof'.format(tumor, len(cells)))
            continue

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
