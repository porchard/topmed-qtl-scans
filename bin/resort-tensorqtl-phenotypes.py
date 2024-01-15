#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd

# necessary to avoid triggering an error in sQTL mapping
PHENOTYPES = sys.argv[1]

tmp = pd.read_csv(PHENOTYPES, dtype=str, sep='\t')

tmp['gene'] = tmp.ID.str.split(':', expand=True)[4]
tmp['start'] = tmp.start.astype(int)

tmp.sort_values(['#chr', 'start', 'gene']).drop(columns='gene').to_csv(sys.stdout, sep='\t', index=False)