#!/usr/bin/env python3.8

import sys
import pandas as pd
import logging

INPUT, OUTPUT = sys.argv[1:]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

tmp = pd.read_csv(INPUT, sep='\t', index_col=None)
tmp.to_parquet(OUTPUT)
