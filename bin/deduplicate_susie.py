#!/usr/bin/env python

import sys
import pandas as pd
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


def deduplicate_susie(cs_df, verbose=True):
    cs_full = cs_df.copy()

    # remove genes that have no variants that appear in > 1 credible set, as these are definitely not duplicates
    gene_variant_count = cs_full.groupby(['phenotype_id', 'variant_id']).size().rename('n').reset_index()
    gene_variant_count = gene_variant_count[gene_variant_count.n>1]
    cs_reduced = cs_full[cs_full.phenotype_id.isin(gene_variant_count.phenotype_id.to_list())].copy()

    if len(cs_reduced) == 0:
        return cs_full

    # now look for duplicates among the rest
    tmp = [[gene, cs_id, ','.join(df.sort_values('variant_id').variant_id) + '; ' + ','.join(df.sort_values('variant_id').pip.astype(str))] for (gene, cs_id), df in cs_reduced.groupby(['phenotype_id', 'cs_id'])]
    tmp = pd.DataFrame(tmp, columns=['phenotype_id', 'cs_id', 'info'])
    duplicates = []
    for (gene, info), df in tmp.groupby(['phenotype_id', 'info']):
        if len(df) > 1:
            duplicates.append(df)
    if len(duplicates) == 0:
        return cs_full
    for i in range(len(duplicates)):
        duplicates[i]['duplicate_set'] = (i + 1)
    duplicates = pd.concat(duplicates)

    cs_full = cs_full.merge(duplicates[['phenotype_id', 'cs_id', 'duplicate_set']], how='left')

    deduplicated = []
    if len(cs_full[cs_full.duplicate_set.isnull()]) > 0:
        deduplicated.append(cs_full[cs_full.duplicate_set.isnull()])
    for (phenotype_id, duplicate_set), df in cs_full[~cs_full.duplicate_set.isnull()].groupby(['phenotype_id', 'duplicate_set']):
        first_cs = df.cs_id.min()
        deduplicated.append(df[df.cs_id==first_cs])
    deduplicated = pd.concat(deduplicated)

    NUMBER_STARTING_CS = len(cs_df[['phenotype_id', 'cs_id']].drop_duplicates())
    NUMBER_ENDING_CS = len(deduplicated[['phenotype_id', 'cs_id']].drop_duplicates())

    if verbose:
        logging.info('Kept {:,} of {:,} CS'.format(NUMBER_ENDING_CS, NUMBER_STARTING_CS))
    
    return deduplicated[cs_df.columns.to_list()]


if __name__ == '__main__':
    SUSIE_CS = sys.argv[1]
    susie = pd.read_csv(SUSIE_CS, sep='\t', dtype=str)
    dedup = deduplicate_susie(susie)
    dedup.to_csv(sys.stdout, sep='\t', index=False)