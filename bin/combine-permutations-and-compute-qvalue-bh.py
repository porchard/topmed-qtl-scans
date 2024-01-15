#!/usr/bin/env python

import sys
import pandas as pd
from statsmodels.stats.multitest import multipletests

PREFIX = sys.argv[1]
PERMUTATION_FILES = sys.argv[2:]

permutations = pd.concat([pd.read_csv(f, sep='\t') for f in PERMUTATION_FILES]).set_index('phenotype_id')
permutations['qval'] = multipletests(permutations.pval_beta, method='fdr_bh')[1]
out_file = f'{PREFIX}.cis_qtl.txt.gz'
permutations.to_csv(out_file, sep='\t', float_format='%.6g')
