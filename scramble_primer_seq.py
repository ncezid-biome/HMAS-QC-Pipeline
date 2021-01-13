import pandas as pd
import random, sys

"""
This script will scramble the forward sequence of our oligo file, and generate 2 new files,
based on the sequence being scrambled is either of primers or barcodes
"""

# this need to be changed to the actual path
oligo_file = r"C:\CDC\zipped_inputs\amr2020.oligos"
tpes = ['primer','barcode']


for tpe in tpes:
    try:
        df = pd.read_csv(oligo_file, sep='\t', header=None, index_col=0, names=['forward', 'reverse', 'label'])
        fltr = (df.index == tpe)  #filter on the row
        df.loc[fltr,'forward'] = df.loc[fltr,'forward'].apply(lambda x: ''.join(random.sample(x, len(x))))
        df.to_csv(oligo_file + '_' + tpe + '_scramble', sep='\t', header=None)
    except FileNotFoundError:
        print(f"we can't locate the oligo file: {oligo_file}")
        sys.exit(1)

