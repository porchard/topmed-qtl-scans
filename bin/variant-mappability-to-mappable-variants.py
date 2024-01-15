#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import sys

SNP_MAPPABILITY = sys.argv[1]

variants = pd.read_csv(SNP_MAPPABILITY, sep='\t', header=None, names=['variant_id', 'mappability'])
mappable_variants = variants[variants.mappability>=1.0]
mappable_variants = mappable_variants.variant_id.to_list()
for i in mappable_variants:
    print(i)
