#!/usr/bin/env python
"""
Script fo analyze precipitation forecasts and create verification plots 
Stephan Rasp
"""
# Imports
import argparse
import numpy as np
from datetime import timedelta
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, make_timelist, \
    yyyymmddhhmmss, ddhhmmss
from scipy.stats import binned_statistic
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
from cosmo_utils.scores.probab import FSS
from helpers import save_fig_and_log, load_det, load_radar, load_ens, \
    strip_expid, set_plot
from config import *  # Import config file

# Load radar tot mask
totmask = np.load('./radar_tot_mask.npy')

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

n = 21
kernel = np.ones((n, n)) / float((n * n))


expid_str = ''
exp_list = []
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
        fss_list = []
    else:
        spread_list = []
        rmse_list = []

    for t in timelist:
        print t
        date = yyyymmddhhmmss(t)
        savestr = args.ana + '_' + expid + '_' + date
        savefn = savedir + savestr + '.npy'
        print 'Try to load pre-saved data:', savefn
        if os.path.exists(savefn):
            print 'Found pre-saved data.'
            if args.ana == 'det':
                radar, detmean, detrmse, fss = np.load(savefn)
            if args.ana == 'ens':
                spread, rmse = np.load(savefn)
        else:
            print 'Did not find pre-saved data, compute!'

            if args.ana == 'det':
                radar = []
                detmean = []
                detrmse = []
                fss = []
            else:
                spread = []
                rmse = []

            for h in range(1, args.hint + 1):
                print 'Hour', h
                hourlist.append(h)
                # Load data
                radarfobj = load_radar(date, ddhhmmss(timedelta(hours=h)))
                nanradar = radarfobj.data
                nanradar[totmask] = np.nan
                convradar = convolve2d(nanradar, kernel, mode='same')
                if args.ana == 'det':
                    # Det run
                    detfobj = load_det(DATA_DIR, date,
                                       ddhhmmss(timedelta(hours=h)))

                    radar.append(np.mean(radarfobj.data[~totmask]))
                    detmean.append(np.mean(detfobj.data[~totmask]))

                    # Upscaled RMSE
                    nandet = detfobj.data
                    nandet[totmask] = np.nan
                    convdet = convolve2d(nandet, kernel, mode='same')
                    detradarmean = 0.5 * (convradar + convdet)
                    detrmse.append(
                        np.sqrt(np.nanmean(((convradar - convdet) ** 2 /
                                            (detradarmean) ** 2)[
                                               detradarmean >= 0.1])))
                    fss.append(FSS(1, 21, nanradar, nandet,
                                   python_core=True))

                else:
                    ensfobjlist = load_ens(DATA_DIR, date,
                                           ddhhmmss(timedelta(hours=h)))
                    convfieldlist = []
                    for ensfobj in ensfobjlist:
                        nanfield = ensfobj.data
                        nanfield[totmask] = np.nan
                        convfield = convolve2d(nanfield, kernel,
                                               mode='same')
                        convfieldlist.append(convfield)
                    convfieldlist = np.array(convfieldlist)
                    meanfield = np.mean(convfieldlist, axis=0)
                    spread.append(np.nanmean((np.std(convfieldlist, axis=0) /
                                         meanfield)[meanfield >= 0.1]))
                    ensradarmean = 0.5 * (convradar + meanfield)
                    rmse.append(np.sqrt(
                        np.nanmean(((convradar - meanfield) ** 2 /
                                    (ensradarmean) ** 2)[
                                       ensradarmean >= 0.1])))

            # save the data
            if args.ana == 'det':
                print 'Save data: ', savefn
                np.save(savefn, (radar, detmean, detrmse, fss))
            else:
                print 'Save data: ', savefn
                np.save(savefn, (spread, rmse))

        # Append to list over days
        if args.ana == 'det':
            radar_list.append(radar)
            detmean_list.append(detmean)
            detrmse_list.append(detrmse)
            fss_list.append(fss)
        else:
            spread_list.append(spread)
            rmse_list.append(rmse)

    if args.composite:
        # Bin the data
        if args.ana == 'det':
            radarmean = np.mean(radar_list, axis=0)
            detmean = np.mean(detmean_list, axis=0)
            detrmse = np.mean(detrmse_list, axis=0)
            meanfss = np.mean(fss_list, axis=0)
            exp_list.append([radarmean, detmean, detrmse, meanfss])
        if args.ana == 'ens':
            ensspread = np.mean(spread_list, axis=0)
            ensrmse = np.mean(rmse_list, axis=0)
            exp_list.append([ensspread, ensrmse])

    else:
        if args.ana == 'det':
            exp_list.append([radar_list, detmean_list, detrmse_list,
                             fss_list])
        else:
            exp_list.append([spread_list, rmse_list])


# Now plot!
# Set up figure
aspect = 0.75
x = range(1, args.hint + 1)
if args.composite:
    plotstr = args.ana + '_comp_' + args.date_start + '_' + args.date_stop
    if args.no_radar:
        plotstr += '_no_radar'
    fig1, ax1 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * aspect))
    if args.ana == 'det':
        fig2, ax2 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * aspect))
    for ie, expid in enumerate(args.expid):
        if args.ana == 'det':
            if ie == 0 and not args.no_radar:
                ax1.plot(x, exp_list[ie][0],
                         c='k', linewidth=2, label='Radar')
            ax1.plot(x, exp_list[ie][1],
                     c=cdict[expid], linewidth=1.5, label=expid)
            ax2.plot(x, exp_list[ie][3],
                     c=cdict[expid], linewidth=1.5, label=expid)
        if args.ana == 'ens':
            ax1.plot(x, exp_list[ie][0],
                     c=cdict[expid], linewidth=1.5, linestyle='--')
            ax1.plot(x, exp_list[ie][1],
                     c=cdict[expid], linewidth=1.5, label=expid)

    # Finish the plots
    if args.ana == 'det':
        ax1.set_ylabel('Precipitation [mm/h]')
        ax2.set_ylabel('FSS [58.8km]')
        for ax in [ax1, ax2]:
            set_plot(ax, 'Det fc ' + args.date_start[:-4] + '-' +
                     args.date_stop[:-4], args, hourlist_plot)
    if args.ana == 'ens':
        ax1.set_ylabel('Norm spread/rmse')
        set_plot(ax1, 'Normalized spread (--) and RMSE (-)',
                 args, hourlist_plot)

    plotdir = plotdir + expid_str[:-1] + '/verif_fc_prec/'
    if not os.path.exists(plotdir): os.makedirs(plotdir)

    # Save figure and create log str
    if args.ana == 'det':
        save_fig_and_log(fig1, 'diprec_' + plotstr, plotdir)
        save_fig_and_log(fig2, 'fss_' + plotstr, plotdir)
    else:
        save_fig_and_log(fig1, 'ens_spread_' + plotstr, plotdir)

    plt.close('all')

else:
    for it, t in enumerate(timelist):
        date = yyyymmddhhmmss(t)
        plotstr = args.ana + '_' + date
        if args.no_radar:
            plotstr += '_no_radar'
        fig1, ax1 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * aspect))
        if args.ana == 'det':
            fig2, ax2 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * aspect))
        for ie, expid in enumerate(args.expid):
            if args.ana == 'det':
                if ie == 0 and not args.no_radar:
                    ax1.plot(x, exp_list[ie][0][it],
                             c='k', linewidth=2, label='Radar')
                ax1.plot(x, exp_list[ie][1][it],
                         c=cdict[expid], linewidth=1.5, label=expid)
                ax2.plot(x, exp_list[ie][3][it],
                         c=cdict[expid], linewidth=1.5, label=expid)
            if args.ana == 'ens':
                ax1.plot(x, exp_list[ie][0][it],
                         c=cdict[expid], linewidth=1.5, linestyle='--')
                ax1.plot(x, exp_list[ie][1][it],
                         c=cdict[expid], linewidth=1.5, label=expid)

        # Finish the plots
        if args.ana == 'det':
            ax1.set_ylabel('Precipitation [mm/h]')
            ax2.set_ylabel('FSS [58.8km]')
            for ax in [ax1, ax2]:
                set_plot(ax, 'Det fc ' + date[:-4], args, hourlist_plot)
        if args.ana == 'ens':
            ax1.set_ylabel('Norm spread/rmse')
            set_plot(ax1, 'Normalized spread (--) and RMSE (-)',
                     args, hourlist_plot)

        pd = plotdir + expid_str[:-1] + '/verif_fc_prec/'
        if not os.path.exists(pd):
            os.makedirs(pd)

        # Save figure and create log str
        if args.ana == 'det':
            save_fig_and_log(fig1, 'diprec_' + plotstr, pd)
            save_fig_and_log(fig2, 'fss_' + plotstr, pd)
        else:
            save_fig_and_log(fig1, 'ens_spread_' + plotstr, pd)

        plt.close('all')
