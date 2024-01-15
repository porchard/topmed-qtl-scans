#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import subprocess
import tempfile
import gzip
import re
import sys

SNP_MAPPABILITY, VARIANT_LIST = sys.argv[1:]


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


def bed_to_regions(bed_df):
    tmp = bed_df[['chrom', 'start', 'end']]
    return list(tmp.chrom.astype(str) + ':' + tmp.start.astype(str) + '-' + tmp.end.astype(str))


def parse_region(s):
    chrom, start, end = re.match('^(.*)[-:_](\d+)[-:_](\d+)$', s).groups()
    start, end = int(start), int(end)
    return (chrom, start, end)


def tabix(f, region, read_header=True):
    """
    f: file to tabix
    region: region, or list of regions, in format chr:start-end (can use :,-,_ as separators)
    read_header: use the last header line prefixed with '#' as column names in the returned dataframe
    """
    if not isinstance(region, list) and not isinstance(region, str):
        raise TypeError('region must be a list or a str')
    if not isinstance(f, str):
        raise TypeError('f must be a str (path to a file)')
    if not isinstance(read_header, bool):
        raise TypeError('read_header must be a bool')
    if isinstance(region, list):
        regions = []
        for r in region:
            chrom, start, end = parse_region(r)
            start, end = str(start), str(end)
            regions.append([chrom, start, end])
        with tempfile.NamedTemporaryFile() as tmpf:
            pd.DataFrame(regions).to_csv(tmpf.name, sep='\t', index=False, header=False)
            sp = subprocess.run(['tabix', '--regions', tmpf.name, f], capture_output=True, check=True)
    else:
        chrom, start, end = parse_region(region)
        sp = subprocess.run(['tabix', f, f'{chrom}:{start}-{end}'], capture_output=True, check=True)
    txt = sp.stdout.decode().split('\n')
    if txt[-1] == '':
        txt = txt[:-1]
    
    if len(txt) == 0:
        return None

    header = None
    if read_header:
        with gzip.open(f, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    header = line.lstrip('#').rstrip().split('\t')
                else:
                    break
    
    return pd.DataFrame([i.split('\t') for i in txt], columns=header)



variants = pd.read_csv(VARIANT_LIST, header=None, index_col=None)[0]

snp_mappability = tabix(SNP_MAPPABILITY, bed_to_regions(variants_to_bed(variants)))
snp_mappability.columns = ['chrom', 'start', 'end', 'mappability']
snp_mappability = {f'{chrom}_{end}': float(mappability) for chrom, end, mappability in zip(snp_mappability.chrom, snp_mappability.end, snp_mappability.mappability)}

for v in variants:
    chrom, pos, ref, alt = v.split('_')
    print('{}\t{}'.format(v, snp_mappability[f'{chrom}_{pos}']))
