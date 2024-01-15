#!/usr/bin/env python

import sys
import pandas as pd

F = sys.argv[1]

x = pd.read_parquet(F)
COLS = x.columns.to_list()
x[['#chrom', 'end']] = x.variant_id.str.split('_', expand=True).iloc[:,0:2]
x['start'] = x.end.astype(int) - 1
x[['#chrom', 'start', 'end'] + COLS].to_csv(sys.stdout, sep='\t', index=False)