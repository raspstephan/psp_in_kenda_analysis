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
from scipy.stats import binned_statistic
from helpers import save_fig_and_log

from config import *

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--expid', metavar='expid', type=str, nargs='+',
                    help='Experiment ID')
parser.add_argument('--date_ana_start', metavar='date_ana_start', type=str,
                    help='Start date of verification (yyyymmddhhmmss)')
parser.add_argument('--date_ana_stop', metavar='date_ana_stop', type=str,
                    help='End date of verification (yyyymmddhhmmss)')
parser.add_argument('--composite', metavar='composite', type=str,
                    default='False',
                    help='If True diurnal composite is plotted for SYNOP')
args = parser.parse_args()

plotstr = ('prec_' + args.date_ana_start + '_' +
           args.date_ana_stop)
if args.composite == 'True':
    plotstr += '_composite'

# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_ana_start)
tend = yyyymmddhhmmss_strtotime(args.date_ana_stop)
tint = timedelta(hours=1)
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)

hourlist_plot = []
for i in range(args.hint + 1):
    hourlist_plot.append(str((tstart + timedelta(hours=i)).hour))

# Set up figure
fig, axarr = plt.subplots(1, 2, figsize=(10, 6))

# Loop over experiments
expid_str = ''
for ie, expid in enumerate(args.expid):
    print 'expid = ', expid
    expid_str += expid + '_'
    DATA_DIR = datadir + expid

    # Check if saved data is available
    savedir = savedir_base + expid + '/verif_ana_prec/'
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
        hourlist = []
        # Loop over time
        for t in timelist:
            # Determine 3hrly storage time
            date_ana = yyyymmddhhmmss(t)
            print 'date_ana = ', date_ana
            hourlist.append(t.hour)
            date_fg = yyyymmddhhmmss(t - timedelta(hours=1))
            date_store = yyyymmddhhmmss(t - timedelta(hours=t.hour % 3))

            # Directory where lff*_prec files are stored
            DATA_DIR = (datadir_da + expid +
                        '/' + date_store + '/')

            # Load relevant files
            radarfobj = load_radar(date_ana)
            detfobj = load_det_da(DATA_DIR, date_fg)
            ensfobjlist = load_ens_da(DATA_DIR, date_fg)

            # Compute mean, spread and rmse
            mask = radarfobj.mask  # True where NO data!
            radarmean.append(np.mean(radarfobj.data[~mask]))
            detmean.append(np.mean(detfobj.data[~mask]))
            detrmse.append(np.sqrt(np.mean((radarfobj.data -
                                            detfobj.data)[~mask] ** 2)))
            ensmeanfobj, ensspreadfobj = mean_spread_fieldobjlist(ensfobjlist)
            ensmean.append(np.mean(ensmeanfobj.data[~mask]))
            ensspread.append(np.mean(ensspreadfobj.data[~mask]))

        # Average if composite
        if args.composite == 'True':
            hour_bins = np.arange(-0.5, 24.5, 1)
            radarmean = binned_statistic(hourlist, radarmean,
                                         bins=hour_bins)[0]
            detmean = binned_statistic(hourlist, detmean,
                                       bins=hour_bins)[0]
            detrmse = binned_statistic(hourlist, detrmse,
                                       bins=hour_bins)[0]
            ensmean = binned_statistic(hourlist, ensmean,
                                       bins=hour_bins)[0]
            ensspread = binned_statistic(hourlist, ensspread,
                                         bins=hour_bins)[0]
        else:
            radarmean = np.array(radarmean)
            detmean = np.array(detmean)
            detrmse = np.array(detrmse)
            ensmean = np.array(ensmean)
            ensspread = np.array(ensspread)

        # End timeloop, save the data lists
        np.save(savefn, (radarmean, detmean, detrmse, ensmean, ensspread))

    # Plot
    axarr[0].plot(range(radarmean.shape[0]), radarmean, c='k', linewidth=2)
    axarr[0].plot(range(radarmean.shape[0]), detmean, c=cdict[expid],
                  linewidth=2,
                  label=expid)
    axarr[0].plot(range(radarmean.shape[0]), ensmean, c=cdict[expid],
                  linewidth=2,
                  linestyle='--')
    axarr[1].plot(range(radarmean.shape[0]), detrmse, c=cdict[expid],
                  linewidth=2)
    axarr[1].plot(range(radarmean.shape[0]), ensspread, c=cdict[expid],
                  linewidth=2,
                  linestyle='--')

# End expid loop


# Finish the plots
if args.composite == 'True':
    axarr[0].set_xlabel('time [UTC/h]')
    axarr[1].set_xlabel('time [UTC/h]')
else:
    axarr[0].set_xlabel('time from ' + yyyymmddhhmmss(tstart) + ' [h]')
    axarr[1].set_xlabel('time from ' + yyyymmddhhmmss(tstart) + ' [h]')
axarr[0].set_label('[mm/h]')
axarr[1].set_label('[mm/h]')
axarr[0].set_title('domain mean precipitation\n det (solid), ens (dashed)')
axarr[1].set_title('det rmse (solid) and ens spread (dashed)')
axarr[0].legend(loc=0, fontsize=6)
plt.tight_layout(rect=[0, 0.0, 1, 0.95])

plotdir = plotdir + expid_str[:-1] + '/verif_ana_prec/'
if not os.path.exists(plotdir): os.makedirs(plotdir)
fig.suptitle(expid_str[:-1] + '  ' + plotstr)
save_fig_and_log(fig, plotstr, plotdir)
plt.close('all')
