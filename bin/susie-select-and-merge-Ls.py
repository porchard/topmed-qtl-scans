#!/usr/bin/env python3
# coding: utf-8

import sys
import numpy as np
import pickle
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

SUSIE_FULL_RESULTS = sys.argv[2:]
OUT = sys.argv[1]


def get_number_cs(x, deduplicate=True):
    """
    x: SuSiE results
    """
    if x['sets']['cs'] is None:
        return 0
    else:
        if deduplicate:
            sets = set([tuple(v) for v in x['sets']['cs'].values()])
            return len(sets)
        else:
            return len(x['sets']['cs'])


susie = dict()
for f in SUSIE_FULL_RESULTS:
    logging.info(f'Loading file {f}')
    with open(f, 'rb') as fh:
        tmp = pickle.load(fh)
        # determine L based on max L observed
        L = max([v['lbf_variable'].shape[0] for v in tmp.values()])
        logging.info(f'File seems to represent L = {L}')
        for gene in tmp.keys():
            ## determine L
            #L = tmp[gene]['lbf_variable'].shape[0] # Note: L might be < the set L if the gene is being tested against fewer than L SNPs
            logging.info(f'Loading gene {gene}')
            if gene not in susie:
                susie[gene] = dict()
            assert(L not in susie[gene])
            susie[gene][L] = tmp[gene]

# now for each gene, select ideal L
susie_final = dict()

for gene in susie.keys():
    Ls = np.array(sorted(susie[gene].keys()))
    logging.info('Processing gene {} (has Ls {})'.format(gene, ', '.join(Ls.astype(str))))
    # determine if any converged. If none converged, keep the min L even though this will be ignored in most downstream analyses
    converged = np.array([susie[gene][L]['converged'] for L in Ls])
    #cs = [susie[gene][L]['sets']['cs'] for L in Ls]
    #number_cs = np.array([len(x) if x is not None else 0 for x in cs])
    number_cs = np.array([get_number_cs(susie[gene][L], deduplicate=False) for L in Ls])
    logging.info('Converged: {}'.format(', '.join(converged.astype(str))))
    logging.info('Number CS: {}'.format(', '.join(number_cs.astype(str))))
    if sum(converged) == 0:
        # If none converged, keep the min L even though this will be ignored in most downstream analyses
        selected_L = min(Ls)
        logging.info(f'None converged; keeping min L ({selected_L})')
        susie_final[gene] = susie[gene][selected_L]
        susie_final[gene]['L'] = selected_L
        continue
    # drop any Ls w/o convergence
    Ls = Ls[converged]
    number_cs = number_cs[converged]
    selected_L = min(Ls[Ls>=max(number_cs)])
    logging.info(f'Selected L = {selected_L}')
    susie_final[gene] = susie[gene][selected_L]
    susie_final[gene]['L'] = selected_L

# now output the new SuSiE pickled object
with open(OUT, 'wb') as f:
    pickle.dump(susie_final, f)
