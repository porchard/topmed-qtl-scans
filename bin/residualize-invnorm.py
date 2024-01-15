#!/usr/bin/env python3.8
# coding: utf-8

import torch
import pandas as pd
from tensorqtl import genotypeio, cis, trans, core
import logging
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--covariates', required=True)
parser.add_argument('--phenotypes', required=True)
parser.add_argument('--plink-prefix', required=True)
parser.add_argument('--prefix', required=True)
args = parser.parse_args()

# COVARIATES = '/net/topmed10/working/porchard/rnaseq/work/tensorqtl-in/freeze-beta/results/tensorqtl-in/T_cell.tensorqtl-in.30.covariates.tsv'
# PHENOTYPES = '/net/topmed10/working/porchard/rnaseq/work/tensorqtl-in/freeze-beta/results/tensorqtl-in/T_cell.tensorqtl-in.phenotypes.bed.gz'
# PLINK_PREFIX = '/net/topmed10/working/porchard/rnaseq/test/genotypes-for-residualization/T_cell.genotypes'

COVARIATES = args.covariates
PHENOTYPES = args.phenotypes
PLINK_PREFIX = args.plink_prefix
PREFIX = args.prefix

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')


logging.info('Reading covariates')
covariates = pd.read_csv(COVARIATES, sep='\t', index_col=0).T

logging.info('Reading phenotypes')
phenotypes = pd.read_csv(PHENOTYPES, sep='\t') # rows should be phenotypes, cols should be samples
phenotypes_for_residualization = phenotypes.iloc[:,3:].set_index('gene_id')
phenotype_columns = list(phenotypes_for_residualization.columns)
phenotype_index = list(phenotypes_for_residualization.index)

logging.info('Reading genotypes')
pr = genotypeio.PlinkReader(PLINK_PREFIX)
genotype_df = pr.load_genotypes()
genotype_df_for_residualization = genotype_df[covariates.index.to_list()]
genotype_columns = list(genotype_df_for_residualization.columns)
genotype_index = list(genotype_df_for_residualization.index)


assert list(phenotypes_for_residualization.columns) == list(covariates.index)
assert list(genotype_df_for_residualization.columns) == list(covariates.index)


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
residualizer = core.Residualizer(torch.tensor(covariates.values, dtype=torch.float32).to(device))

def residualize(df, residualizer, inverse_normal_transform, chunk_size = 10000):
    current_chunk_start = 0
    residuals = []
    while current_chunk_start <= len(df):
        current_chunk_end = min(chunk_size + current_chunk_start, len(df) + 1)
        r = residualizer.transform(torch.tensor(df.iloc[current_chunk_start:current_chunk_end,:].values, dtype=torch.float).to(device), inverse_normal_transform=inverse_normal_transform)
        residuals.append(pd.DataFrame(r.cpu().numpy(), columns=df.columns, index=df.index.to_list()[current_chunk_start:current_chunk_end]))
        current_chunk_start = current_chunk_start + chunk_size
    return pd.concat(residuals)

logging.info('Residualizing phenotypes')
residualized_phenotype_df = residualize(phenotypes_for_residualization, residualizer, inverse_normal_transform=True)

logging.info('Residualizing genotypes')
residualized_genotype_df = residualize(genotype_df_for_residualization, residualizer, inverse_normal_transform=False)

# output
logging.info('Outputting phenotypes')
residualized_phenotype_df.to_parquet(f'{PREFIX}residualized-phenotypes.parquet')

logging.info('Outputting genotypes')
residualized_genotype_df.to_parquet(f'{PREFIX}residualized-genotypes.parquet')
