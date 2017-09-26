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
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, make_timelist, \
    yyyymmddhhmmss, ddhhmmss
from scipy.stats import binned_statistic
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
from cosmo_utils.scores.probab import FSS
from helpers import *
from config import *  # Import config file

# Load radar tot mask
totmask = np.load('../aux_data/radar_tot_mask.npy')

# Arguments
parser = argparse.ArgumentParser(description='Process input')
parser.add_argument('--expid', metavar='expid', type=str, nargs='+',
                    help='Experiment ID')
parser.add_argument('--date_start',
                    metavar='date_start',
                    type=str,
                    default='20160526000000',
                    help='Start date for date loop (yyyymmddhhmmss)')
parser.add_argument('--date_stop',
                    metavar='date_stop',
                    type=str,
                    default='20160609000000',
                    help='End date for date loop (yyyymmddhhmmss)')
parser.add_argument('--date_inc', metavar='date_inc', type=int,
                    default='24',
                    help='Time increment between forecasts (h)')
parser.add_argument('--hint', metavar='hint', type=int, default=24,
                    help='Maximum forecast lead time')
parser.add_argument('--ana', metavar='ana', type=str,
                    help='Type of analysis to be done [det or ens]')
parser.add_argument('--ens_norm_type',
                    type=int,
                    default=0,
                    help='Type of ensemble normalization. '
                         '0 (no normalization), 1 or 2.')
parser.add_argument('--n_kernel',
                    type=int,
                    default=21,
                    help='Width of convolution kernel. Default=21')
parser.add_argument('--composite',
                    dest='composite',
                    action='store_true',
                    help='Composite or individual plots.')
parser.set_defaults(composite=False)
parser.add_argument('--no_radar',
                    dest='no_radar',
                    action='store_true',
                    help='Do not plot radar curve')
parser.set_defaults(no_radar=False)
parser.add_argument('--radar_thresh',
                    type=float,
                    default=100.,
                    help='Radar values above threshold will be set to nan.'
                         'Default = 100.')
parser.add_argument('--recompute',
                    dest='recompute',
                    action='store_true',
                    help='Recompute pre-processed files.')
parser.set_defaults(recompute=False)
args = parser.parse_args()

assert args.ana in ['det', 'ens'], 'Wrong analysis!'

# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_start)
tend = yyyymmddhhmmss_strtotime(args.date_stop)
tint = timedelta(hours=args.date_inc)
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)

hourlist_plot = []
for i in range(args.hint + 1):
    hourlist_plot.append(str((tstart + timedelta(hours=i)).hour))

assert args.n_kernel % 2 == 1, 'n_kernel must be odd'
kernel = (np.ones((args.n_kernel, args.n_kernel)) /
          float((args.n_kernel * args.n_kernel)))


expid_str = ''
results_dict = {}
for ie, expid in enumerate(args.expid):
    print 'expid = ', expid
    expid_str += expid + '_'
    DATA_DIR = datadir + expid

    # Check if saved data is available
    savedir = savedir_base + expid + '/verif_fc_prec/'
    if not os.path.exists(savedir): os.makedirs(savedir)

    hourlist = []
    if args.ana == 'det':
        radar_list = []
        detmean_list = []
        detrmse_list = []
        fss01_list = []
        fss10_list = []
    else:
        radar_list = []
        ensrmse_list = []
        ensrmv_list = []
        ensbs_list = []
        ensmean_list = []
        ensmean_std_list = []

    for t in timelist:
        print t
        date = yyyymmddhhmmss(t)
        savestr = (args.ana + '_' + expid + '_' + date + '_n_' +
                   str(args.n_kernel) + '_norm_' + str(args.ens_norm_type))
        savefn = savedir + savestr + '.npy'
        print 'Try to load pre-saved data:', savefn
        if os.path.exists(savefn) and not args.recompute:
            print 'Found pre-saved data.'
            if args.ana == 'det':
                radar, detmean, detrmse, fss01, fss10 = np.load(savefn)
            if args.ana == 'ens':
                radar, ensrmse, ensrmv, ensbs, ensmean, ensmean_std = \
                    np.load(savefn)
        else:
            print 'Did not find pre-saved data, compute!'

            if args.ana == 'det':
                radar = []
                detmean = []
                detrmse = []
                fss01 = []
                fss10 = []
            else:
                radar = []
                ensrmse = []
                ensrmv = []
                ensbs = []
                ensmean = []
                ensmean_std = []

            # Load preloaded precipitation data
            if args.ana == 'det':
                preloaded_fn = (savedir_base + expid + '/prec_fields/det_' + 
                                date + '.npy')
                prec_fields = np.load(preloaded_fn)
            elif args.ana == 'ens':
                preloaded_fn = (savedir_base + expid + '/prec_fields/ens_' + 
                                date + '.npy')
                prec_fields = np.load(preloaded_fn)

            for h in range(1, args.hint + 1):
                print 'Hour', h
                hourlist.append(h)
                # Load data
                radarfobj = load_radar(date, ddhhmmss(timedelta(hours=h)))
                tmpmask = np.logical_or(totmask,
                                        radarfobj.data > args.radar_thresh)
                nanradar = radarfobj.data
                nanradar[tmpmask] = np.nan
                convradar = convolve2d(nanradar, kernel, mode='same')
                radar.append(np.mean(radarfobj.data[~tmpmask]))
                if args.ana == 'det':
                    # Det run
                    prec_hour = prec_fields[h - 1]

                    detmean.append(np.mean(prec_hour[~tmpmask]))

                    # Upscaled RMSE
                    nandet = prec_hour
                    nandet[tmpmask] = np.nan
                    convdet = convolve2d(nandet, kernel, mode='same')

                    r = compute_det_stats(convradar, convdet, nanradar, nandet)
                    detrmse.append(r[0])
                    fss01.append(r[1])
                    fss10.append(r[2])


                else:
                    ensfobjlist = load_ens(DATA_DIR, date,
                                           ddhhmmss(timedelta(hours=h)))
                    prec_hour = prec_fields[h -1]
                    convfieldlist = []
                    for p in prec_hour:
                        nanfield = p
                        nanfield[tmpmask] = np.nan
                        convfield = convolve2d(nanfield, kernel,
                                               mode='same')
                        convfieldlist.append(convfield)
                    convfieldlist = np.array(convfieldlist)
                    r = compute_ens_stats(convradar, convfieldlist,
                                          args.ens_norm_type,
                                          norm_thresh=0.1)
                    ensrmse.append(r[0])
                    ensrmv.append(r[1])
                    ensbs.append(r[2])
                    ensmean.append(r[3])
                    ensmean_std.append(r[4])

            # save the data
            if args.ana == 'det':
                print 'Save data: ', savefn
                np.save(savefn, (radar, detmean, detrmse, fss01, fss10))
            else:
                print 'Save data: ', savefn
                np.save(savefn, (radar, ensrmse, ensrmv, ensbs, ensmean,
                                 ensmean_std))

        # Append to list over days
        if args.ana == 'det':
            radar_list.append(radar)
            detmean_list.append(detmean)
            detrmse_list.append(detrmse)
            fss01_list.append(fss01)
            fss10_list.append(fss10)
        else:
            radar_list.append(radar)
            ensrmse_list.append(ensrmse)
            ensrmv_list.append(ensrmv)
            ensbs_list.append(ensbs)
            ensmean_list.append(ensmean)
            ensmean_std_list.append(ensmean_std)


    if args.composite:
        # Bin the data
        if args.ana == 'det':
            radarmean = np.mean(radar_list, axis=0)
            detmean = np.mean(detmean_list, axis=0)
            detrmse = np.mean(detrmse_list, axis=0)
            meanfss01 = np.mean(fss01_list, axis=0)
            meanfss10 = np.mean(fss10_list, axis=0)
            results_dict.update({
                expid : {
                'radar' : radarmean,
                'detmean' : detmean,
                'detrmse' : detrmse,
                'fss01' : meanfss01,
                'fss10' : meanfss10,
                }
            })
        if args.ana == 'ens':
            radarmean = np.mean(radar_list, axis=0)
            ensrmse = np.mean(ensrmse_list, axis=0)
            ensrmv = np.mean(ensrmv_list, axis=0)
            ensbs = np.mean(ensbs_list, axis=0)
            ensmean = np.mean(ensmean_list, axis=0)
            ensmean_std = np.mean(ensmean_std_list, axis=0)
            results_dict.update({
                expid : {
                'radar' : radarmean,
                'ensrmse' : ensrmse,
                'ensrmv' : ensrmv,
                'ensbs' : ensbs,
                'ensmean' : ensmean,
                'ensmean_std' : ensmean_std,
                }
            })

    else:
        if args.ana == 'det':
            results_dict.update({
                expid : {
                'radar' : np.asarray(radar_list),
                'detmean' : np.asarray(detmean_list),
                'detrmse' : np.asarray(detrmse_list),
                'fss01' : np.asarray(fss01_list),
                'fss10' : np.asarray(fss10_list),
                }
            })
        else:
            results_dict.update({
                expid : {
                'radar' : np.asarray(radar_list),
                'ensrmse' : np.asarray(ensrmse_list),
                'ensrmv' : np.asarray(ensrmv_list),
                'ensbs' : np.asarray(ensbs_list),
                'ensmean' : np.asarray(ensmean_list),
                'ensmean_std' : np.asarray(ensmean_std_list),
                }
            })


# Now plot!
# Set up figure
aspect = 0.75
x = range(1, args.hint + 1)
if args.composite:
    titlestr = args.date_start[:-4] + '-' + args.date_stop[:-4]
    plotstr = (args.ana + '_comp_' + args.date_start + '_' + args.date_stop + '_n_' +
                   str(args.n_kernel) + '_norm_' + str(args.ens_norm_type))
    pd = plotdir + expid_str[:-1] + '/verif_fc_prec/'
    if not os.path.exists(pd): os.makedirs(pd)

    if args.ana == 'det':
        det_fig = make_fig_fc_det_fss(args, x, results_dict,
                                  titlestr)

        save_fig_and_log(det_fig, plotstr, pd)

    else:

        ens_fig = make_fig_fc_ens(args, x, results_dict, titlestr)

        save_fig_and_log(ens_fig, plotstr, pd)

    plt.close('all')

else:
    for it, t in enumerate(timelist):
        date = yyyymmddhhmmss(t)
        plotstr = (args.ana + '_' + date + '_n_' +
                   str(args.n_kernel) + '_norm_' + str(args.ens_norm_type))
        pd = plotdir + expid_str[:-1] + '/verif_fc_prec/'
        if not os.path.exists(pd):
            os.makedirs(pd)

        if args.ana == 'det':

            det_fig = make_fig_fc_det_fss(args, x, results_dict, date[:-4], it)

            save_fig_and_log(det_fig, plotstr, pd)

        else:

            ens_fig = make_fig_fc_ens(args, x, results_dict, date[:-4], it)

            save_fig_and_log(ens_fig, plotstr, pd)

        plt.close('all')

