#!/usr/bin/env python3.8

import sys
import pandas as pd

OUTPUT_FILE = sys.argv[1]
INPUT_FILES = sys.argv[2:]

merged = pd.concat([pd.read_parquet(f) for f in INPUT_FILES])
merged.to_parquet(OUTPUT_FILE)