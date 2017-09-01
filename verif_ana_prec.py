#!/usr/bin/env python
"""
Compute precipitation spread and rmse for analysis 
"""

# Imports
import argparse
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from cosmo_utils.helpers import yyyymmddhhmmss, yyyymmddhhmmss_strtotime, \
    make_timelist, yymmddhhmm, ddhhmmss
from cosmo_utils.pyncdf import getfobj_ncdf
from cosmo_utils.pywgrib import getfobj_ens, getfobj
from cosmo_utils.diag import mean_spread_fieldobjlist
from datetime import timedelta
from scipy.stats import binned_statistic
from helpers import save_fig_and_log, load_radar, load_det_da, load_ens_da, \
    set_plot

from config import *

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--expid', metavar='expid', type=str, nargs='+',
                    help='Experiment ID')
parser.add_argument('--date_ana_start', metavar='date_ana_start', type=str,
                    help='Start date of verification (yyyymmddhhmmss)')
parser.add_argument('--date_ana_stop', metavar='date_ana_stop', type=str,
                    help='End date of verification (yyyymmddhhmmss)')
parser.add_argument('--composite',
                    dest='composite',
                    action='store_true',
                    help='Composite or individual plots.')
parser.set_defaults(composite=False)
args = parser.parse_args()

# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_ana_start)
tend = yyyymmddhhmmss_strtotime(args.date_ana_stop)
tint = timedelta(hours=1)
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)

# Set up figure
fig, axarr = plt.subplots(1, 2, figsize=(10, 6))

# Loop over experiments
expid_str = ''
exp_list = []
for ie, expid in enumerate(args.expid):
    print 'expid = ', expid
    expid_str += expid + '_'
    DATA_DIR = datadir + expid

    # Check if saved data is available
    savedir = savedir_base + expid + '/verif_ana_prec/'
    if not os.path.exists(savedir): os.makedirs(savedir)

    savefn = (savedir + expid + '_' + args.date_ana_start + '_' +
              args.date_ana_stop + '.npy')
    print 'Try to load pre-saved data:', savefn
    if os.path.exists(savefn):
        print 'Found pre-saved data.'
        # These objects are lists containing the hourly mean values
        # Corresponding to timelist. For each expid!
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

        print 'Save data: ', savefn
        np.save(savefn, (radarmean, detmean, detrmse, ensmean, ensspread))

    # Average if composite
    if args.composite:
        hourlist = [t.hour for t in timelist]
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

    exp_list.append([radarmean, detmean, detrmse, ensmean, ensspread])

# Now the plotting...
aspect = 0.75
if args.composite:
    x = np.unique(hourlist)
    plotstr = 'comp_' + args.date_ana_start + '_' + args.date_ana_stop
    fig1, ax1 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * aspect))
    fig2, ax2 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * aspect))
    for ie, expid in enumerate(args.expid):
        if ie == 0:
            ax1.plot(x, exp_list[ie][0],
                     c='k', linewidth=2, label='Radar')
        ax1.plot(x, exp_list[ie][1],
                 c=cdict[expid], linewidth=1.5, label=expid)
        ax1.plot(x, exp_list[ie][2],
                 c=cdict[expid], linewidth=1.5, linestyle='--')
        ax2.plot(x, exp_list[ie][3],
                 c=cdict[expid], linewidth=1.5, label=expid)
        ax2.plot(x, exp_list[ie][4],
                 c=cdict[expid], linewidth=1.5, linestyle='--')

    ax1.set_ylabel('Prec mean/rmse [mm/h]')
    ax2.set_ylabel('Prec ens rmse/spread [mm/h]')
    for ax in [ax1, ax2]:
        set_plot(ax, 'Det fc ' + args.date_ana_start[:-4] + '-' +
                 args.date_ana_stop[:-4], args, x)

    plotdir = plotdir + expid_str[:-1] + '/verif_ana_prec/'
    if not os.path.exists(plotdir): os.makedirs(plotdir)

    save_fig_and_log(fig1, 'det_' + plotstr, plotdir)
    save_fig_and_log(fig2, 'ens_' + plotstr, plotdir)
    plt.close('all')

else:
    x = range(24)
    # Get array indices
    index_start = 0
    index_stop = 24
    daylist = make_timelist(tstart, tend, timedelta(hours=24))[:-1]
    if tstart == tend:
        daylist = [tstart]
    for iday, day in enumerate(daylist):
        date = yyyymmddhhmmss(day)
        plotstr = date
        fig1, ax1 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * aspect))
        fig2, ax2 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * aspect))
        for ie, expid in enumerate(args.expid):
            if ie == 0:
                ax1.plot(x, exp_list[ie][0][index_start:index_stop],
                         c='k', linewidth=2, label='Radar')
            ax1.plot(x, exp_list[ie][1][index_start:index_stop],
                     c=cdict[expid], linewidth=1.5, label=expid)
            ax1.plot(x, exp_list[ie][2][index_start:index_stop],
                     c=cdict[expid], linewidth=1.5, linestyle='--')
            ax2.plot(x, exp_list[ie][3][index_start:index_stop],
                     c=cdict[expid], linewidth=1.5, label=expid)
            ax2.plot(x, exp_list[ie][4][index_start:index_stop],
                     c=cdict[expid], linewidth=1.5, linestyle='--')

        ax1.set_ylabel('Prec mean/rmse [mm/h]')
        ax2.set_ylabel('Prec ens rmse/spread [mm/h]')
        for ax in [ax1, ax2]:
            set_plot(ax, 'Det fc ' + date[:-4], args, x)

        pd = plotdir + expid_str[:-1] + '/verif_ana_prec/'
        if not os.path.exists(pd): os.makedirs(pd)

        save_fig_and_log(fig1, 'det_' + plotstr, pd)
        save_fig_and_log(fig2, 'ens_' + plotstr, pd)
        plt.close('all')

        index_start += 24
        index_stop += 24
#
# a=b
#     # Plot
#     axarr[0].plot(range(radarmean.shape[0]), radarmean, c='k', linewidth=2)
#     axarr[0].plot(range(radarmean.shape[0]), detmean, c=cdict[expid],
#                   linewidth=2,
#                   label=expid)
#     axarr[0].plot(range(radarmean.shape[0]), ensmean, c=cdict[expid],
#                   linewidth=2,
#                   linestyle='--')
#     axarr[1].plot(range(radarmean.shape[0]), detrmse, c=cdict[expid],
#                   linewidth=2)
#     axarr[1].plot(range(radarmean.shape[0]), ensspread, c=cdict[expid],
#                   linewidth=2,
#                   linestyle='--')
#
# # End expid loop
#
#
# # Finish the plots
# if args.composite == 'True':
#     axarr[0].set_xlabel('time [UTC/h]')
#     axarr[1].set_xlabel('time [UTC/h]')
# else:
#     axarr[0].set_xlabel('time from ' + yyyymmddhhmmss(tstart) + ' [h]')
#     axarr[1].set_xlabel('time from ' + yyyymmddhhmmss(tstart) + ' [h]')
# axarr[0].set_label('[mm/h]')
# axarr[1].set_label('[mm/h]')
# axarr[0].set_title('domain mean precipitation\n det (solid), ens (dashed)')
# axarr[1].set_title('det rmse (solid) and ens spread (dashed)')
# axarr[0].legend(loc=0, fontsize=6)
# plt.tight_layout(rect=[0, 0.0, 1, 0.95])
#
# plotdir = plotdir + expid_str[:-1] + '/verif_ana_prec/'
# if not os.path.exists(plotdir): os.makedirs(plotdir)
# fig.suptitle(expid_str[:-1] + '  ' + plotstr)
# save_fig_and_log(fig, plotstr, plotdir)
# plt.close('all')
