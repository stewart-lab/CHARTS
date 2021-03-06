import numpy as np
import pandas as pd
import scanpy as sc
import sys
import os 
from os.path import join
import subprocess
from anndata import AnnData
from optparse import OptionParser
from collections import defaultdict
import h5py

sys.path.append('..')

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-r", "--resolution", help="Resolution")
    parser.add_option("-w", "--overwrite", action="store_true", help="Overwrite existing clustering data in the database")
    (options, args) = parser.parse_args()

    h5_f = args[0]
    out_dir = options.out_dir
    overwrite = options.overwrite

    if options.resolution:
        resolution = float(options.resolution)
    else:
        resolution = 1.0
    
    sc.settings.verbosity = 3
    sc.logging.print_versions()
    
    the_tumors = set()
    with h5py.File(h5_f, 'r') as f:
        the_tumors = f['per_tumor'].keys()
        the_tumors = sorted(the_tumors)
    print(the_tumors)

    for tumor in the_tumors:
        # Determine whether to cluster this tumor or not
        if not overwrite:
            with h5py.File(h5_f, 'r') as f:
                if 'cluster' in f['per_tumor/{}'.format(tumor)].keys():
                    print('Found tumor {} in the database. Skipping clustering.'.format(tumor))
                    continue

        print("Clustering tumor {}".format(tumor))
        with h5py.File(h5_f, 'r') as f:
            cells = [
                str(x)[2:-1]
                for x in f['per_tumor/{}/cell'.format(tumor)][:]
            ]
            expression = f['per_tumor/{}/log1_tpm'.format(tumor)][:]

        print('Shape of matrix: ', expression.shape)

        #if expression.shape[0] > 10:
        ad = AnnData(
            X=expression, 
            obs=pd.DataFrame(data=cells, columns=['cell'], index=cells)
        )
        sc.pp.pca(ad, n_comps=min([50, expression.shape[0]-1]))
        sc.pp.neighbors(ad)

        if expression.shape[0] < 50:
            print("Number of cells is only {}. Overriding resolution to be 1.0".format(expression.shape[0]))
            use_resolution = 1.0
        else:
            use_resolution = resolution
        sc.tl.leiden(ad, resolution=use_resolution)

        clusters = ad.obs['leiden']
        assert tuple(ad.obs.index) == tuple(cells)
        
        with h5py.File(h5_f, 'r+') as f:
            try:
                del f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, int(resolution))]
            except KeyError:
                pass
            f.create_dataset(
                'per_tumor/{}/leiden_res_{}/cluster'.format(tumor, int(resolution)), 
                data=np.array(clusters, dtype=np.int32), 
                compression="gzip"
            )
        #else:
        #    with h5py.File(h5_f, 'r+') as f:
        #        try:
        #            del f['per_tumor/{}/leiden_res_{}/cluster'.format(tumor, int(resolution))]
        #        except KeyError:
        #            pass
        #        f.create_dataset(
        #            'per_tumor/{}/leiden_res_{}/cluster'.format(tumor, int(resolution)),
        #            data=np.array(
        #                [0 for i in range(expression.shape[0])], # Assign all cells to one cluster 
        #                dtype=np.int32
        #            ),
        #            compression="gzip"
        #        )

if __name__ == '__main__':
    main()

