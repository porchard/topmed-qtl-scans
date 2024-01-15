#!/usr/bin/env python
# coding: utf-8

import sys
from qtl import norm
import pandas as pd

COVARIATES = sys.argv[1]

covariates = pd.read_csv(COVARIATES, sep='\t', index_col=0, dtype=str)
samples = covariates.columns.to_list()

# for  norm.inverse_normal_transform, samples should be columns (as they already are in the covariates file)
phenotype_pcs = covariates[covariates.index.to_series().str.contains('phenotype_PC')].astype(float)
pcs = norm.inverse_normal_transform(phenotype_pcs)

pcs['#chr'] = 'chr0'
pcs['start'] = 1000
pcs['end'] = pcs.start + 1
pcs['gene_id'] = pcs.index.to_list()
pcs = pcs[['#chr', 'start', 'end', 'gene_id'] + samples]
pcs.to_csv(sys.stdout, sep='\t', index=False)
