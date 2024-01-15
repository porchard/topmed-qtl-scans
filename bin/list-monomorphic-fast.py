#!/usr/bin/env python
# coding: utf-8

import sys
import os
import tempfile
import subprocess
import pandas as pd
from pandas_plink import read_plink

PLINK_BED = sys.argv[1]

# assumes all individuals have no missing genotypes
# to be monomorphic, one of the following must be true
# * all inidividuals must be homozygous for same allele (maf = 0)
# * all individuals must be heterozygous (maf = 0.5)

# to avoid reading entire plink file, first find variants that meet these criteria
def plink_to_mafs(bfile, samples=None, variants=None, plink_command='plink'):
    """
    Get MAFs from plink binary files
    """
    assert(samples is None or isinstance(samples, list))
    assert(variants is None or isinstance(variants, list))
    
    abs_plink_path = os.path.abspath(bfile)
    starting_dir = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        args = [plink_command, '--freq', '--out', 'getmaf', '--bfile', abs_plink_path]
        if samples is not None:
            with open('keep-samples.txt', 'w') as fh:
                for s in samples:
                    fh.write('{}\t{}\n'.format(0, s))
            args += ['--keep', 'keep-samples.txt']
        if variants is not None:
            with open('keep-variants.txt', 'w') as fh:
                for s in variants:
                    fh.write('{}\n'.format(s))
            args += ['--extract', 'keep-variants.txt']
        sp = subprocess.run(args, capture_output=True, check=True)
        freqs = pd.read_csv('getmaf.frq', delim_whitespace=True)
        os.chdir(starting_dir)
    return freqs

mafs = plink_to_mafs(PLINK_BED.replace('.bed', ''))
variants_to_check = mafs[(mafs.MAF<=0.001) | (mafs.MAF>=0.499)].SNP.to_list()
    
if len(variants_to_check) > 0:
    (bim, fam, bed) = read_plink(PLINK_BED, verbose=False)
    m = bim['snp'].isin(variants_to_check).values
    bed = bed[m,:]
    bim = bim[m]
    df = pd.DataFrame(bed.compute(), index=bim.snp, columns=fam.iid)
    df = df.loc[variants_to_check]

    mono = (df.nunique(axis=1) == 1)

    for i in mono[mono].index:
        print(i)
