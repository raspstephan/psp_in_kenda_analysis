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
from helpers import *
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
parser.add_argument('--recompute',
                    dest='recompute',
                    action='store_true',
                    help='Recompute pre-processed files.')
parser.set_defaults(recompute=False)
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
    if os.path.exists(savefn) and not args.recompute:
        print 'Found pre-saved data.'
        # These objects are lists containing the hourly mean values
        # Corresponding to timelist. For each expid!
        radarmean, detmean, detrmse, fss01, fss10, ensrmse, ensrmv, ensbs, ensmean, \
            ensmean_std = np.load(savefn)
    else:
        print 'Did not find pre-saved data, compute!'
        radarmean = []
        detmean = []
        detrmse = []
        fss01 = []
        fss10 = []
        ensrmse = []
        ensrmv = []
        ensbs = []
        ensmean = []
        ensmean_std = []
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

            r = compute_det_stats(convradar, convdet, nanradar, nandet)
            detrmse.append(r[0])
            fss01.append(r[1])
            fss10.append(r[2])

            r = compute_ens_stats(convradar, convfieldlist,
                                  args.ens_norm_type,
                                  norm_thresh=0.1)
            ensrmse.append(r[0])
            ensrmv.append(r[1])
            ensbs.append(r[2])
            ensmean.append(r[3])
            ensmean_std.append(r[4])

        print 'Save data: ', savefn
        np.save(savefn, (radarmean, detmean, detrmse, fss01, fss10, ensrmse,
                         ensrmv, ensbs, ensmean, ensmean_std))

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
        fss01 = binned_statistic(hourlist, fss01,
                                   bins=hour_bins, statistic=np.nanmean)[0]
        fss10 = binned_statistic(hourlist, fss10,
                               bins=hour_bins, statistic=np.nanmean)[0]
        ensrmse = binned_statistic(hourlist, ensrmse,
                                   bins=hour_bins, statistic=np.nanmean)[0]
        ensrmv = binned_statistic(hourlist, ensrmv,
                                     bins=hour_bins, statistic=np.nanmean)[0]
        ensbs = binned_statistic(hourlist, ensbs,
                                     bins=hour_bins, statistic=np.nanmean)[0]
        ensmean = binned_statistic(hourlist, ensmean,
                                     bins=hour_bins, statistic=np.nanmean)[0]
        ensmean_std = binned_statistic(hourlist, ensmean_std,
                                   bins=hour_bins, statistic=np.nanmean)[0]

    exp_list.append([radarmean, detmean, detrmse, fss01, fss10, ensrmse, ensrmv,
                     ensbs, ensmean, ensmean_std])

exp_list = np.array(exp_list)

# Now the plotting...
if args.composite:
    x = range(24)
    plotstr = ('comp_' + args.date_ana_start + '_' + args.date_ana_stop +
               '_n_' + str(args.n_kernel) + '_norm_' + str(args.ens_norm_type))
    titlestr = args.date_ana_start[:-4] + '-' + args.date_ana_stop[:-4]

    radar = exp_list[:, 0]
    meanprec = exp_list[:, 1]
    rmse = exp_list[:, 2]
    fss01 = exp_list[:, 3]
    fss10 = exp_list[:, 4]
    det_fig = make_fig_fc_det_fss(args, x, radar, meanprec, fss01, fss10,
                                  titlestr)

    ensmean = exp_list[:, 8]
    ensmean_std = exp_list[:, 9]
    ensrmse = exp_list[:, 5]
    ensrmv = exp_list[:, 6]
    ensbs = exp_list[:, 7]
    ens_fig = make_fig_fc_ens(args, x, radar, ensmean, ensmean_std, ensrmse,
                              ensrmv, ensbs, titlestr)

    pd = plotdir + expid_str[:-1] + '/verif_ana_prec/'
    if not os.path.exists(pd): os.makedirs(pd)

    save_fig_and_log(det_fig, 'det_' + plotstr, pd)
    save_fig_and_log(ens_fig, 'ens_' + plotstr, pd)
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

        radar = exp_list[:, 0, index_start:index_stop]
        meanprec = exp_list[:, 1, index_start:index_stop]
        rmse = exp_list[:, 2, index_start:index_stop]
        fss01 = exp_list[:, 3, index_start:index_stop]
        fss10 = exp_list[:, 4, index_start:index_stop]
        det_fig = make_fig_fc_det_fss(args, x, radar, meanprec, fss01, fss10,
                                  date[:-4])

        ensmean = exp_list[:, 8, index_start:index_stop]
        ensmean_std = exp_list[:, 9, index_start:index_stop]
        ensrmse = exp_list[:, 5, index_start:index_stop]
        ensrmv = exp_list[:, 6, index_start:index_stop]
        ensbs = exp_list[:, 7, index_start:index_stop]

        ens_fig = make_fig_fc_ens(args, x, radar, ensmean, ensmean_std, ensrmse,
                                  ensrmv, ensbs, date[:-4])

        pd = plotdir + expid_str[:-1] + '/verif_ana_prec/'
        if not os.path.exists(pd): os.makedirs(pd)

        save_fig_and_log(det_fig, 'det_' + plotstr, pd)
        save_fig_and_log(ens_fig, 'ens_' + plotstr, pd)
        plt.close('all')

        index_start += 24
        index_stop += 24





