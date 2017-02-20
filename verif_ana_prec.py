#!/usr/bin/env python
"""
Compute precipitation spread and rmse for analysis 
"""

# Imports
import argparse
import numpy as np
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cosmo_utils.helpers import yyyymmddhhmmss, yyyymmddhhmmss_strtotime, \
                                make_timelist, yymmddhhmm, ddhhmmss
from cosmo_utils.pyncdf import getfobj_ncdf
from cosmo_utils.pywgrib import getfobj_ens, getfobj
from cosmo_utils.diag import mean_spread_fieldobjlist
from datetime import timedelta

parser = argparse.ArgumentParser(description = 'Script to verify first guess precipitation forecast from assimilation cycling against radar observations.')
parser.add_argument('--expid', metavar = 'expid', type=str, nargs = '+',
                    help = 'Experiment ID')
parser.add_argument('--date_ana_start', metavar = 'date_ana_start', type=str,
                    help = 'Start date of verification (yyyymmddhhmmss)')
parser.add_argument('--date_ana_stop', metavar = 'date_ana_stop', type=str,
                    help = 'End date of verification (yyyymmddhhmmss)')
args = parser.parse_args()


plotstr = ('prec_' + args.date_ana_start + '_' + 
           args.date_ana_stop)

# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_ana_start)
tend = yyyymmddhhmmss_strtotime(args.date_ana_stop)
tint = timedelta(hours = 1)
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)

# Settings
radardir = '/e/uwork/extsrasp/radolan/'
radarpref = 'raa01-rw_10000-'
radarsufx = '-dwd---bin.nc'
precsuf = '_prec'
gribpref = 'lff'
nens = 40

cdict = {'radar':'k',
        'REF':'navy',
        'REF_TL500':'darkgreen',
        'PSP_TL500':'maroon',
        'DA_REF':'blue',
        'DA_REF_TL500':'cyan',
        'DA_PSP_TL500':'red',
        'DA_PSPv2_TL500':'magenta',
        'DA_PSP':'maroon',
        }


# Define loading functions (duplicate from plot_stamps.py)
def load_radar(t):
    dateobj = (yyyymmddhhmmss_strtotime(t))
    radardt = timedelta(minutes = 10)   # TODO Is this correct???
    radardateobj = dateobj - radardt
    radarfn = radardir + radarpref + yymmddhhmm(radardateobj) + radarsufx
    print 'radarfn = ', radarfn
    radarfobj = getfobj_ncdf(radarfn, fieldn = 'pr', dwdradar = True)
    return radarfobj

def load_det(DATA_DIR, t):
    detfn = DATA_DIR + gribpref + t + precsuf + '.det'
    print 'detfn = ', detfn
    detfobj = getfobj(detfn, fieldn = 'TOT_PREC_S')
    return detfobj

def load_ens(DATA_DIR, t):
    print 'Load ensemble'
    gribfn = gribpref + t + precsuf
    ensfobjlist = getfobj_ens(DATA_DIR, 'same', mems = nens, gribfn = gribfn,
                              fieldn = 'TOT_PREC_S', para = 4)
    return ensfobjlist

# Set up figure
fig, axarr = plt.subplots(1, 2, figsize = (10, 6))

# Loop over experiments
expid_str = ''
for ie, expid in enumerate(args.expid):
    print 'expid = ', expid
    expid_str += expid + '_'
    
    # Check if saved data is available
    savedir = '/e/uwork/extsrasp/save/' + expid + '/verif_ana_prec/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    savefn = savedir + plotstr + '.npy'
    print 'Try to load pre-saved data:', savefn
    if os.path.exists(savefn):
        print 'Found pre-saved data.'
        radarmean, detmean, detrmse, ensmean, ensspread = np.load(savefn)
    else:
        print 'Did not find pre-saved data, compute!'
        radarmean = []
        detmean = []
        detrmse = []
        ensmean = []
        ensspread = []
        # Loop over time
        for t in timelist:
            # Determine 3hrly storage time
            date_ana = yyyymmddhhmmss(t)
            print 'date_ana = ', date_ana
            date_fg = yyyymmddhhmmss(t - timedelta(hours = 1))
            date_store = yyyymmddhhmmss(t - timedelta(hours = t.hour % 3))

            # Directory where lff*_prec files are stored
            DATA_DIR = ('/e/uwork/extsrasp/cosmo_letkf/data/' + expid +
                        '/' + date_store + '/')

            # Load relevant files
            radarfobj = load_radar(date_ana)
            detfobj = load_det(DATA_DIR, date_fg)
            ensfobjlist = load_ens(DATA_DIR, date_fg)

            # Compute mean, spread and rmse
            mask = radarfobj.mask   # True where NO data!
            radarmean.append(np.mean(radarfobj.data[~mask]))
            detmean.append(np.mean(detfobj.data[~mask]))
            detrmse.append(np.sqrt(np.mean((radarfobj.data - 
                                            detfobj.data)[~mask]**2)))
            ensmeanfobj, ensspreadfobj = mean_spread_fieldobjlist(ensfobjlist)
            ensmean.append(np.mean(ensmeanfobj.data[~mask]))
            ensspread.append(np.mean(ensspreadfobj.data[~mask]))

        # End timeloop, save the data lists
        np.save(savefn, (radarmean, detmean, detrmse, ensmean, ensspread))
    
    # Plot
    axarr[0].plot(range(1, len(timelist)+1), radarmean, c = 'k', linewidth = 2)
    axarr[0].plot(range(1, len(timelist)+1), detmean, c = cdict[expid], linewidth = 2,
                  label = expid)
    axarr[0].plot(range(1, len(timelist)+1), ensmean, c = cdict[expid], linewidth = 2,
                  linestyle = '--')
    axarr[1].plot(range(1, len(timelist)+1), detrmse, c = cdict[expid], linewidth = 2)
    axarr[1].plot(range(1, len(timelist)+1), ensspread, c = cdict[expid], linewidth = 2,
                  linestyle = '--')

# End expid loop


# Finish the plots
axarr[0].set_xlabel('time from ' +  yyyymmddhhmmss(tstart - 
                                                 timedelta(hours = 1))+ ' [h]')
axarr[0].set_label('[mm/h]')
axarr[1].set_xlabel('time from ' +  yyyymmddhhmmss(tstart - 
                                                 timedelta(hours = 1))+ ' [h]')
axarr[1].set_label('[mm/h]')
axarr[0].set_title('domain mean precipitation\n det (solid), ens (dashed)')
axarr[1].set_title('det rmse (solid) and ens spread (dashed)')
axarr[0].legend(loc = 0, fontsize = 6)
plt.tight_layout(rect=[0, 0.0, 1, 0.95])

plotdir = '/e/uwork/extsrasp/plots/' + expid_str[:-1] + '/verif_ana_prec/'
if not os.path.exists(plotdir): os.makedirs(plotdir)
fig.suptitle(expid_str[:-1] + '  ' + plotstr)
print 'Save figure:', plotdir + plotstr
fig.savefig(plotdir + plotstr, dpi = 300)
plt.close('all')







