# /hwfssz5/ST_SUPERCELLS/P21Z10200N0171/USER/wangyijin/00.Software/miniconda/envs/scenicplus/bin/python
import argparse
import os, sys

import loompy as lp
import numpy as np
import scanpy as sc

output="./All.loom"
input="./myo_trace_genes_count.csv.gz"
def main():
    x=sc.read_csv(input)
    row_attrs = {"Gene": np.array(x.var_names),}
    col_attrs = {"CellID": np.array(x.obs_names)}
    lp.create(output, x.X.transpose(), row_attrs, col_attrs)

if __name__ == '__main__':
    main()
