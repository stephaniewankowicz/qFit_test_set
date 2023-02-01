#!/usr/bin/env python
#last edited: 2020-01-19
#last edited by: Stephanie Wankowicz
#git:

import pandas as pd
import os
import datetime
import argparse
import sys

def parse_log(log_holo, log_apo, holo, apo):
    rval_holo = pd.DataFrame()
    rval_apo = pd.DataFrame()
    log = open(log_holo, 'r')
    rval_holo.loc[1,'PDB'] = holo
    rval_apo.loc[2,'PDB'] = apo
    for line in log:
        if line.startswith('Final R-work'):
            rval_holo.loc[1,'Rwork'] = line.split('=')[1][1:6]
            rval_holo.loc[1,'Rfree'] = line.split('=')[2][1:6]
    log = open(log_apo, 'r')
    for line in log:
      if line.startswith('Final R-work'):
            rval_apo.loc[2,'Rwork'] = line.split('=')[1][1:6]
            rval_apo.loc[2,'Rfree'] = line.split('=')[2][1:6]
    rval_holo.to_csv(holo + '_prerefine_rvalues.csv', index=False)
    rval_apo.to_csv(apo + '_qFit_rvalues.csv', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('Holo_Log_File')
    parser.add_argument('Apo_Log_File')
    parser.add_argument('Holo')
    parser.add_argument('Apo')
    args = parser.parse_args()
    parse_log(args.Holo_Log_File, args.Apo_Log_File, args.Holo, args.Apo)
