#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd

#TOP = '/net/topmed10/working/porchard/rnaseq/data/tensorqtl-out/tensorqtl-freeze-2aRNA/trans-eqtl/maf005/trans-top/PBMC.trans_qtl.top.txt'
#PHENOTYPES = '/net/topmed10/working/porchard/rnaseq/work/tensorqtl-in/freeze-beta/results/tensorqtl-in/PBMC.tensorqtl-in.phenotypes.bed.gz'
#OUT = 

TOP, PHENOTYPES, OUT = sys.argv[1:]


top = pd.read_csv(TOP, sep='\t')
trans_eqtl = top[top.qvalue<=0.05]

if len(trans_eqtl) > 0:
    phenotypes = pd.read_csv(PHENOTYPES, sep='\t')
    phenotypes = phenotypes.drop(columns=['#chr', 'start', 'end'])
    phenotypes.columns = ['gene_id'] + phenotypes.columns.to_list()[1:]

    new_phenotypes = trans_eqtl[['variant_id', 'phenotype_id']].rename(columns={'phenotype_id': 'gene_id'})
    new_phenotypes[['#chr', 'end']] = new_phenotypes.variant_id.str.split('_', expand=True)[[0, 1]]
    new_phenotypes.end = new_phenotypes.end.astype(int)
    new_phenotypes['start'] = new_phenotypes.end - 1
    new_phenotypes = new_phenotypes[['#chr', 'start', 'end', 'gene_id']]
    new_phenotypes = new_phenotypes.merge(phenotypes)
    assert(len(new_phenotypes) == len(trans_eqtl))

    new_phenotypes = new_phenotypes.sort_values(['#chr', 'start'], ascending=True)
    new_phenotypes.rename(columns={'chr': '#chr'}).to_csv(OUT, sep='\t', index=False)
