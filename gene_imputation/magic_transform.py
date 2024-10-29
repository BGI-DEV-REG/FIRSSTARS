# /jdfssz1/ST_SUPERCELLS/P21Z10200N0171/USER/daixi/miniconda3/envs/r-reticulate/bin/python
# magic transform
import magic
import pandas as pd
import argparse
import pickle
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-M', '--model')
parser.add_argument('-G', '--gene')
parser.add_argument('-O', '--output')
parser.add_argument('-F', '--force', action='store_true')

args = parser.parse_args()
print(args)


if os.path.exists(args.output) and not args.force:
    print(f"文件 {args.output} 已存在。使用 --force 参数覆盖文件或选择其他文件名。")
    sys.exit(1)



with open(args.model, 'rb') as f:
    magic_operator = pickle.load(f)

X_magic = magic_operator.transform(genes = [args.gene])
X_magic.to_csv(args.output)

