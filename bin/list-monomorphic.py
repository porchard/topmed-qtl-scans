#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
from pandas_plink import read_plink

PLINK_BED = sys.argv[1]

(bim, fam, bed) = read_plink(PLINK_BED, verbose=False)
df = pd.DataFrame(bed.compute(), index=bim.snp, columns=fam.iid)

mono = (df.nunique(axis=1) == 1)

for i in mono[mono].index:
    print(i)
