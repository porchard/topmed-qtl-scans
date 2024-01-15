#!/usr/bin/env python3.8
# coding: utf-8

import pandas as pd
import pickle
import matplotlib.pyplot as plt
#import seaborn as sns
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--susie-pickle', nargs='+', required=True)
parser.add_argument('--susie-parquet', nargs='+', required=True)
parser.add_argument('--prefix', required=True)
parser.add_argument('--by-cluster', default=False, action='store_true', help='Collapse by cluster (default: by gene)')
args = parser.parse_args()

#SUSIE_PICKLE = glob.glob('/net/topmed10/working/porchard/rnaseq/data/tensorqtl-out/tensorqtl-freeze-2a/cis-sqtl/susie/maf001/Lung.*.10.susie.pickle')
#SUSIE_TXT = glob.glob('/net/topmed10/working/porchard/rnaseq/data/tensorqtl-out/tensorqtl-freeze-2a/cis-sqtl/susie/maf001/Lung.*.10.susie.txt.gz')
#NOMINAL_SIGNIF = '/net/topmed10/working/porchard/rnaseq/work/internal-freezes/freeze-2a/cis-sqtl/Lung.maf001.cis_qtl.signif_pairs.parquet'

SUSIE_PICKLE = args.susie_pickle
SUSIE_PARQUET = args.susie_parquet
PREFIX = args.prefix


# reduce credible sets: 
# among overlapping credible sets from different introns of the same gene, select the CS with the highest PIP.

# for each gene:
# identify overlapping credible sets
# for each set of overlapping sets, drop the credible set with the lowest max PIP
# will do this by:
# ranking credible sets by max pip
# then, in order of rank:
# check if CS overlaps any CS already seen
# if yes, skip to next CS
# if no, add it to the CS to keep


def _prune_cs_one_gene(df):
    cols = df.columns.to_list()
    df['unique_cs_id'] = df.phenotype_id + '___' + df.cs_id.astype(str)
    max_pip_per_cs = df.groupby('unique_cs_id').pip.max().reset_index().sort_values('pip', ascending=False)
    max_pip_per_cs['rnk'] = range(len(max_pip_per_cs))
    df = df.merge(max_pip_per_cs[['unique_cs_id', 'rnk']]).sort_values('rnk')

    included_snps = set()
    pruned = []
    for rnk, x in df.groupby('rnk', sort=True):
        if any(x.variant_id.isin(included_snps)):
            continue
        else:
            included_snps.update(x.variant_id.to_list())
            pruned.append(x)
    pruned = pd.concat(pruned)
    return pruned[cols]


def prune_cs(df, cluster_level=False):
    if cluster_level:
        pruned = pd.concat([_prune_cs_one_gene(x) for (g, c), x in df.groupby(['gene', 'cluster'])])
    else:    
        pruned = pd.concat([_prune_cs_one_gene(x) for g, x in df.groupby('gene')])
    return pruned



cs = pd.concat([pd.read_parquet(f) for f in SUSIE_PARQUET]) #######
cs['gene'] = cs.phenotype_id.str.split(':', expand=True)[4]
cs['cluster'] = cs.phenotype_id.str.split(':', expand=True)[3]
cs_pruned = prune_cs(cs) if not args.by_cluster else prune_cs(cs, cluster_level=True)

assert(cs.gene.nunique() == cs_pruned.gene.nunique())


# combine the SuSiE objects and filter to phenotypes with (pruned) credible sets
# note that many credible sets will remain within the SuSiE objects that will need to be filtered out from coloc results (e.g., a phenotype might have one credible set that survives pruning and one that does not)
susie = dict()
for f in SUSIE_PICKLE:
    with open(f, 'rb') as fh:
        susie.update(pickle.load(fh))



DROP_KEYS = [i for i in susie.keys() if i not in set(cs_pruned.phenotype_id.unique())]
for i in DROP_KEYS:
    del susie[i]

assert(len(susie) == cs_pruned.phenotype_id.nunique())

# output new pickle
# output pruned cs
cs_pruned[['phenotype_id', 'variant_id', 'pip', 'af', 'cs_id']].to_csv(f'{PREFIX}cs.txt', sep='\t', index=False)
cs[['phenotype_id', 'variant_id', 'pip', 'af', 'cs_id']].to_csv(f'{PREFIX}per-intron-cs.txt', sep='\t', index=False)

with open(f'{PREFIX}pickle', 'wb') as f:
    pickle.dump(susie, f)


#fig, axs = plt.subplots(figsize=(5*3, 3), ncols=3)
#
#ax = axs[0]
#tmp = cs.groupby('phenotype_id').cs_id.nunique().rename('n').reset_index()
#sns.histplot(x='n', data=tmp, discrete=True, ax=ax)
#ax.set_title('Per intron')
#ax.set_xlabel('Number of CS')
#
#ax = axs[1]
#tmp = cs[['gene', 'phenotype_id', 'cs_id']].drop_duplicates().groupby('gene').size().rename('n').reset_index()
#sns.histplot(x='n', data=tmp, discrete=True, ax=ax)
#ax.set_title('Per gene (pre-collapse)')
#ax.set_xlabel('Number of CS')
#
#ax = axs[2]
#tmp = cs_pruned[['gene', 'phenotype_id', 'cs_id']].drop_duplicates().groupby('gene').size().rename('n').reset_index()
#sns.histplot(x='n', data=tmp, discrete=True, ax=ax)
#ax.set_title('Per gene (post-collapse)')
#ax.set_xlabel('Number of CS')
#
#fig.tight_layout()
#fig.savefig(f'{PREFIX}number-cs.png', dpi=300, facecolor='white')
#fig.clf()
