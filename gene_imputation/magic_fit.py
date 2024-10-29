# /jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/02.ST/21.magic/01.script
# /jdfssz1/ST_SUPERCELLS/P21Z10200N0171/USER/daixi/miniconda3/envs/r-reticulate/bin/python
# magic fit
import magic
import pandas as pd
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input')
parser.add_argument('-O', '--output')

args = parser.parse_args()
print(args)



X = pd.read_csv(args.input, index_col=0)
magic_operator = magic.MAGIC() 
magic_operator.fit(X)


with open(args.output, 'wb') as f:
    pickle.dump(magic_operator, f)


