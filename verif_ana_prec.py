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
    set_plot, compute_ens_stats
from scipy.signal import convolve2d

from config import *

# Load radar tot mask
totmask = np.load('./radar_tot_mask.npy')
n = 21
kernel = np.ones((n, n)) / float((n * n))

# Arguments
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
parser.add_argument('--radar_thresh',
                    type=float,
                    default=100.,
                    help='Radar values above threshold will be set to nan.'
                         'Default = 100.')
parser.add_argument('--n_kernel',
                    type=int,
                    default=21,
                    help='Width of convolution kernel. Default=21')
parser.add_argument('--ens_norm_type',
                    type=int,
                    default=1,
                    help='Type of ensemble normalization. '
                         '0 (no normalization), 1 or 2.')
args = parser.parse_args()

assert args.n_kernel % 2 == 1, 'n_kernel must be odd'
kernel = (np.ones((args.n_kernel, args.n_kernel)) /
          float((args.n_kernel * args.n_kernel)))

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
              args.date_ana_stop + '_n_' + str(args.n_kernel) + '_norm_' +
              str(args.ens_norm_type) + '_.npy')
    print 'Try to load pre-saved data:', savefn
    if os.path.exists(savefn):
        print 'Found pre-saved data.'
        # These objects are lists containing the hourly mean values
        # Corresponding to timelist. For each expid!
        radarmean, detmean, detrmse, ensrmse, ensspread = np.load(savefn)
    else:
        print 'Did not find pre-saved data, compute!'
        radarmean = []
        detmean = []
        detrmse = []
        ensrmse = []
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
            try:
                radardata = load_radar(date_ana, return_array=True)
                detdata = load_det_da(DATA_DIR, date_fg, return_array=True)
                ensdatalist = load_ens_da(DATA_DIR, date_fg, return_array=True)
            except:
                print 'Some file was not found, write nans.'
                radardata = np.zeros((461, 421)) * np.nan
                detdata = np.zeros((461, 421)) * np.nan
                ensdatalist = [np.zeros((461, 421)) * np.nan for i in range(20)]

            # Thresholding and upscaling
            tmpmask = np.logical_or(totmask, radardata > args.radar_thresh)
            nanradar = radardata
            nanradar[tmpmask] = np.nan
            convradar = convolve2d(nanradar, kernel, mode='same')
            nandet = detdata
            nandet[tmpmask] = np.nan
            convdet = convolve2d(nandet, kernel, mode='same')
            convfieldlist = []
            for ensdata in ensdatalist:
                nanfield = ensdata
                nanfield[tmpmask] = np.nan
                convfield = convolve2d(nanfield, kernel,
                                       mode='same')
                convfieldlist.append(convfield)
            convfieldlist = np.array(convfieldlist)

            # Compute the statistics
            radarmean.append(np.nanmean(nanradar))
            detmean.append(np.nanmean(nandet))
            detrmse.append(np.sqrt(np.nanmean((convradar - convdet) ** 2)))
            s, r = compute_ens_stats(convradar, convfieldlist,
                                     args.ens_norm_type,
                                     norm_thresh=0.1)
            ensspread.append(s)
            ensrmse.append(r)

        print 'Save data: ', savefn
        np.save(savefn, (radarmean, detmean, detrmse, ensrmse, ensspread))

    # Average if composite
    if args.composite:
        hourlist = [t.hour for t in timelist]
        hour_bins = np.arange(-0.5, 24.5, 1)
        radarmean = binned_statistic(hourlist, radarmean,
                                     bins=hour_bins, statistic=np.nanmean)[0]
        detmean = binned_statistic(hourlist, detmean,
                                   bins=hour_bins, statistic=np.nanmean)[0]
        detrmse = binned_statistic(hourlist, detrmse,
                                   bins=hour_bins, statistic=np.nanmean)[0]
        ensrmse = binned_statistic(hourlist, ensrmse,
                                   bins=hour_bins, statistic=np.nanmean)[0]
        ensspread = binned_statistic(hourlist, ensspread,
                                     bins=hour_bins, statistic=np.nanmean)[0]

    exp_list.append([radarmean, detmean, detrmse, ensrmse, ensspread])

# Now the plotting...
aspect = 0.75
if args.composite:
    x = np.unique(hourlist)
    plotstr = ('comp_' + args.date_ana_start + '_' + args.date_ana_stop +
               '_n_' + str(args.n_kernel) + '_norm_' + str(args.ens_norm_type))
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
    ax2.set_ylabel('ens norm rmse/spread [mm/h]')
    for ax in [ax1, ax2]:
        set_plot(ax, args.date_ana_start[:-4] + '-' +
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
        plotstr = (date + '_n_' + str(args.n_kernel) + '_norm_' +
                   str(args.ens_norm_type))
        fig1, ax1 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * aspect))
        fig2, ax2 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * aspect))
        for ie, expid in enumerate(args.expid):
            if ie == 0:
                print len(exp_list[ie][0])
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
        ax2.set_ylabel('ens norm rmse/spread [mm/h]')
        for ax in [ax1, ax2]:
            set_plot(ax, date[:-4], args, x)

        pd = plotdir + expid_str[:-1] + '/verif_ana_prec/'
        if not os.path.exists(pd): os.makedirs(pd)

        save_fig_and_log(fig1, 'det_' + plotstr, pd)
        save_fig_and_log(fig2, 'ens_' + plotstr, pd)
        plt.close('all')

        index_start += 24
        index_stop += 24
