#!/usr/bin/env python
"""
Script fo analyze precipitation forecasts and create verification plots 
Stephan Rasp
"""
# Imports
import argparse
import os
import numpy as np
from datetime import timedelta
from cosmo_utils.da import fdbkfile
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, make_timelist, \
    yyyymmddhhmmss
from scipy.stats import binned_statistic
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Define loading functions
def load_det(expid, date, t):
    topdir = datadir + expid + '/' + date + '/'
    gribfn = gribpref + t + precsuf
    detfn = topdir + 'det/' + gribfn
    detfobj = getfobj(detfn, fieldn = 'PREC_PERHOUR')
    return detfobj


# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--expid', metavar = 'expid', type=str, nargs = '+',
                    help = 'Experiment ID')
parser.add_argument('--date_start', metavar = 'date_start', type=str,
                    default = '20160525000000', help = 'Start date for date loop (yyyymmddhhmmss)')
parser.add_argument('--date_stop', metavar = 'date_stop', type=str,
                    default = '20160610000000', 
                    help = 'End date for date loop (yyyymmddhhmmss)')
parser.add_argument('--date_inc', metavar = 'date_inc', type=int, 
                    default = '24', 
                    help = 'Time increment between forecasts (h)')
parser.add_argument('--hint', metavar = 'hint', type=int, default =24,
                    help = 'Maximum forecast lead time')
args = parser.parse_args()

# Config for experiment
cdict = {'radar':'k',
             'REF':'navy',
             'REF_TL500':'darkgreen',
             'PSP_TL500':'orange',
            'DA_REF':'blue',
            'DA_REF_TL500':'cyan',
            'DA_PSP_TL500':'red',
            'DA_PSPv2_TL500':'magenta',
            'DA_PSP':'maroon',
             }

# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_start)
tend = yyyymmddhhmmss_strtotime(args.date_stop)
tint = timedelta(hours=args.date_inc)
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)



expid_str = ''
for ie, expid in enumerate(args.expid):
    print 'expid = ', expid
    expid_str += expid + '_'
    DATA_DIR='/e/uwork/extsrasp/cosmo_letkf/data_forecast/' + expid
    
    for t in timelist:
        print t
        for h in range(1, hint+1):
            
        
    
