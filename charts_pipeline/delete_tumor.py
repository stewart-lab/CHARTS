import sys
import h5py

h5_f = sys.argv[1]
tumor = sys.argv[2]

with h5py.File(h5_f, 'r+') as f:
    print(f['per_tumor'].keys())
    print()
    del f['per_tumor/{}'.format(tumor)]
    print(f['per_tumor'].keys())
