#!/usr/bin/env python

import pandas as pd
import sys

TOP, PHENOTYPES = sys.argv[1:]

top = pd.read_csv(TOP, sep='\t')
gene_to_variant = dict(zip(top.phenotype_id, top.variant_id))
phenotypes = pd.read_csv(PHENOTYPES, usecols=list(range(4)), sep='\t')
exon_to_gene = {exon: exon.split('_')[0] for exon in phenotypes.gene_id}
print('phenotype_id\tvariant_id')
for exon, gene in exon_to_gene.items():
    if gene in gene_to_variant:
        print('{}\t{}'.format(exon, gene_to_variant[gene]))