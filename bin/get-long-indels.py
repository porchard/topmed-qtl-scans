#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import glob
import logging
import sys

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

#BIM_FILES = ['/net/topmed10/working/porchard/rnaseq/work/subset-topmed-bcf/freeze-alpha/results/plink-pass-filter/chr20.bim']
#BIM_FILES = glob.glob('/net/topmed10/working/porchard/rnaseq/work/subset-topmed-bcf/freeze-alpha/results/plink-pass-filter/*.bim')
BIM_FILES = sys.argv[1:]

long_variants = []
for f in BIM_FILES:
    logging.info(f'Scanning {f}')
    tmp = pd.read_csv(f, sep='\t', header=None, names=['chrom', 'variant', 'x', 'pos', 'allele_1', 'allele_2'])
    tmp['allele_1_len'] = tmp.allele_1.str.len()
    tmp['allele_2_len'] = tmp.allele_2.str.len()
    tmp['max_len'] = tmp[['allele_1_len', 'allele_2_len']].max(axis=1)
    long_variants.append(tmp[tmp.max_len>20])
long_variants = pd.concat(long_variants)

long_variants = long_variants[long_variants.max_len>=50]
long_variants[['variant']].to_csv(sys.stdout, index=False, header=None)