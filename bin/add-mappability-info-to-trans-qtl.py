#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import logging
import pybedtools as bt
import sys
import argparse
import numpy as np
import re

# NOTE: use GENCODE gtf, which contains all genes in TOPMed GTF plus some (minus ERCCs)
# so might be some cross-mapping instances that aren't caught using TOPMed GTF

parser = argparse.ArgumentParser()
parser.add_argument('--trans')
parser.add_argument('--gtf', required=True, help='GTF file (used to infer TSS)')
parser.add_argument('--cross-mappability', required=True)
parser.add_argument('--cross-mappability-window', default=1000000, type=int, help='Window for removing crossmapping genes (default: 1000000)')
parser.add_argument('--gene-mappability', required=True)
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


# TRANS_TXT = '/net/topmed10/working/porchard/rnaseq/data/tensorqtl-out/tensorqtl-freeze-2aRNA/trans-eqtl/permutations/Whole_blood.trans_qtl_pairs.txt.gz'
# GTF = '/net/topmed10/working/porchard/rnaseq/work/crossmap-gencode-v30/gencode.v30.annotation.gtf'
# CROSS_MAPPABILITY = '/net/topmed10/working/porchard/rnaseq/work/crossmap-gencode-v30/results/cross_mappability/crossmap.txt'
# GENE_MAPPABILITY = '/net/topmed10/working/porchard/rnaseq/work/crossmap-gencode-v30/results/gene-mappability/gene_mappability/gene_mappability.txt'


TRANS_TXT = args.trans
GTF = args.gtf
CROSS_MAPPABILITY = args.cross_mappability
CROSS_MAPPABILITY_WINDOW = args.cross_mappability_window
GENE_MAPPABILITY = args.gene_mappability


def phenotype_id_to_gene_id(x):
    # first, try to match with a version
    ENSEMBL_RE_WITH_VERSION = 'ENSG\d+\.\d+'
    ENSEMBL_RE_WITHOUT_VERSION = 'ENSG\d+'
    with_version = re.search(ENSEMBL_RE_WITH_VERSION, x)
    without_version = re.search(ENSEMBL_RE_WITHOUT_VERSION, x)
    if with_version:
        return with_version.group(0)
    elif without_version:
        return without_version.group(0)
    else:
        raise ValueError(f'Not able to infer gene ID from {x}')



def variants_to_bed (variants):
    # variants is a list of variants in form chrom_pos_ref_alt
    # returns dataframe of chrom, start, end, variant
    if not isinstance(variants, (list, pd.Series, np.ndarray)):
        raise ValueError('variants must be a list, series, or array')
    df = pd.Series(variants).str.split('[_:]', expand=True, regex=True) if isinstance(variants, (np.ndarray, list)) else variants.str.split('[_:]', expand=True, regex=True)
    assert(len(df.columns) in [2, 4])
    df.columns = ['chrom', 'pos', 'ref', 'alt'] if len(df.columns) == 4 else ['chrom', 'pos']
    df['end'] = df.pos.astype(int)
    df['start'] = df.end - 1
    df['name'] = list(variants)
    return df[['chrom', 'start', 'end', 'name']]


def gtf_to_df(gtf):
    df = pd.read_csv(gtf, sep='\t', header=None, names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'], comment='#')
    return df


def parse_attributes(a, regex='[; ]*(.*?) "(.*?)"'):
    """
    Parse the attributes column of a (GENCODE/RefSeq) GTF file.

    Input:
    * a [str]: the attributes element (column 9 of the GTF file)
    * regex [str]: a regular expression that will be iteratively applied to the attribute string to capture attribute key, val pairs. Default should work for GENCODE/RefSeq
    """
    x = [m.groups() for m in re.finditer(regex, a)]
    return {key: val for key, val in x}


def gtf_to_tss(gtf, feature_id='gene_id'):
    """
    Given a GTF file, create a BED6-style DataFrame of TSS.

    Input:
    gtf: Path to GTF file
    feature_id: Attribute to use for labeling TSS (usually e.g. gene_id, gene_name, or transcript_id)

    Output:
    pandas DataFrame of TSS (chrom, start, end, feature_id, ., strand)
    """

    df = gtf_to_df(gtf)
    df = df[df.feature=='transcript']
    df['tss'] = np.where(df.strand == '+', df.start, df.end)
    df['tss_start'] = df.tss - 1  # BED indexes from 0
    df['tss_end'] = df.tss
    df['id'] = df.attributes.map(lambda x: parse_attributes(x)[feature_id])
    df['score'] = '.'
    return df[['chrom', 'tss_start', 'tss_end', 'id', 'score', 'strand']].rename(columns=lambda x: x.replace('tss_', ''))


# GTEx filtering:
# We filtered variants at MAF > 0.05 (within each tissue) and excluded any variant with mappability < 1, based on k-mer length 75
# Candidate trans-eGenes were restricted to protein-coding and lincRNA genes, as annotated in GENCODE v26. 
# Finally, we applied the hg38 cross-mapping filter as described in [16] with settings of k-mer length 75 for exons and 36 for UTRs, 
# applying this to the filtered set of variant-gene pairs with p-values below 10−5 to exclude any gene with mappability < 0.8 and 
# any variant-gene pair where the target eGene cross-maps with any gene within 1Mb of the variant.

gene_mappability = pd.read_csv(GENE_MAPPABILITY, sep='\t', header=None, names=['phenotype_id', 'mappability']).set_index('phenotype_id').mappability.to_dict()

cross_mappability = pd.read_csv(CROSS_MAPPABILITY, sep='\t', header=None, names=['gene_1', 'gene_2', 'n'])
cross_mappability = {gene: set(df.gene_2.unique()) for gene, df in cross_mappability.groupby('gene_1')}

trans = pd.read_csv(TRANS_TXT, sep='\t')

# get the genes w/in --cross-mappability-window of each variant
tss = gtf_to_tss(GTF)
tss = tss[['chrom', 'start', 'end', 'id']].drop_duplicates()
variant_bed = variants_to_bed(trans.variant_id.unique().tolist())
# break down by chrom for efficiency
CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX']
assert(all(variant_bed.chrom.isin(set(CHROMS))))
within_window = []
for chrom in CHROMS:
    within_window_chrom = bt.BedTool().from_dataframe(variant_bed[variant_bed.chrom==chrom]).sort().window(b=bt.BedTool().from_dataframe(tss[tss.chrom==chrom]).sort(), w=int(CROSS_MAPPABILITY_WINDOW)).to_dataframe()
    within_window.append(within_window_chrom[['name', 'thickEnd']].drop_duplicates())
within_window = pd.concat(within_window)



variant_to_nearby_genes = {variant: set() for variant in trans.variant_id.unique()}
for variant, gene in zip(within_window['name'], within_window['thickEnd']):
    variant_to_nearby_genes[variant].add(gene)

def gene_crossmaps_to_gene_near_variant(phenotype_id, variant_id):
    nearby_genes = variant_to_nearby_genes[variant_id]
    if phenotype_id in cross_mappability:
        for x in cross_mappability[phenotype_id]:
            if x in nearby_genes:
                return True
    for x in nearby_genes:
        if x in cross_mappability:
            if phenotype_id in cross_mappability[x]:
                return True
    return False

# Add biotype info
gtf_df = gtf_to_df(GTF)
gtf_df = gtf_df[gtf_df.feature=='gene']
gtf_df['gene_id'] = gtf_df.attributes.map(lambda x: parse_attributes(x)['gene_id'])
gtf_df['gene_type'] = gtf_df.attributes.map(lambda x: parse_attributes(x)['gene_type'])
gene_id_to_biotype = dict(zip(gtf_df.gene_id, gtf_df.gene_type))
# Add gene chromosome
gene_id_to_chrom = dict(zip(gtf_df.gene_id, gtf_df.chrom))


trans['gene_id'] = trans.phenotype_id.map(phenotype_id_to_gene_id)
trans['gene_mappability'] = trans.gene_id.map(gene_mappability)
trans['gene_crossmaps_to_gene_near_variant'] = [gene_crossmaps_to_gene_near_variant(gene_id, variant_id) for variant_id, gene_id in zip(trans.variant_id, trans.gene_id)]
trans['biotype'] = trans.gene_id.map(gene_id_to_biotype)
trans['phenotype_chr'] = trans.gene_id.map(gene_id_to_chrom)
trans = trans.drop(columns=['gene_id'])

trans.to_csv(sys.stdout, sep='\t', index=False)
