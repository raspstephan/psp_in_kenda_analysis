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
    strip_expid
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
                    help='Type of analysis to be done [di]')
args = parser.parse_args()

plotstr = (args.ana + '_' + args.date_start + '_' +
           args.date_stop)

# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_start)
tend = yyyymmddhhmmss_strtotime(args.date_stop)
tint = timedelta(hours=args.date_inc)
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)

hourlist = []
for i in range(args.hint + 1):
    hourlist.append(str((tstart + timedelta(hours=i)).hour))
print hourlist

n = 21
kernel = np.ones((n, n)) / float((n * n))

# Set up figure
fig1, ax1 = plt.subplots(1, 1, figsize=(6, 4))
if args.ana == 'det':
    fig2, ax2 = plt.subplots(1, 1, figsize=(6, 4))

expid_str = ''
for ie, expid in enumerate(args.expid):
    print 'expid = ', expid
    expid_str += expid + '_'
    DATA_DIR = datadir + expid

    # Check if saved data is available
    savedir = savedir_base + expid + '/verif_fc_prec/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    savefn = savedir + plotstr + '.npy'
    print 'Try to load pre-saved data:', savefn
    if os.path.exists(savefn):
        print 'Found pre-saved data.'
        if args.ana == 'det':
            radarmean, detmean, detrmse, meanfss = np.load(savefn)
        if args.ana == 'ens':
            ensspread, ensrmse = np.load(savefn)
    else:
        print 'Did not find pre-saved data, compute!'

        hourlist = []
        if args.ana == 'det':
            radarmean = []
            detmean = []
            detrmse = []
            fsslist = []

        if args.ana == 'ens':
            spreadlist = []
            rmselist = []

        for t in timelist:
            print t
            date = yyyymmddhhmmss(t)
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

                    radarmean.append(np.mean(radarfobj.data[~totmask]))
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
                    fss = FSS(1, 21, nanradar, nandet, python_core=True)
                    fsslist.append(fss)

                if args.ana == 'ens':
                    ensfobjlist = load_ens(DATA_DIR, date,
                                           ddhhmmss(timedelta(hours=h)))
                    convfieldlist = []
                    for ensfobj in ensfobjlist:
                        nanfield = ensfobj.data
                        nanfield[totmask] = np.nan
                        convfield = convolve2d(nanfield, kernel, mode='same')
                        convfieldlist.append(convfield)
                    convfieldlist = np.array(convfieldlist)
                    meanfield = np.mean(convfieldlist, axis=0)
                    spreadlist.append(
                        np.nanmean((np.std(convfieldlist, axis=0) /
                                    meanfield)[meanfield >= 0.1]))
                    ensradarmean = 0.5 * (convradar + meanfield)
                    rmselist.append(
                        np.sqrt(np.nanmean(((convradar - meanfield) ** 2 /
                                            (ensradarmean) ** 2)[
                                               ensradarmean >= 0.1])))



                    # End hourlist
        # End timelist

        # Bin the data 
        hour_bins = np.arange(0.5, 25.5, 1)
        if args.ana == 'det':
            radarmean = binned_statistic(hourlist, radarmean,
                                         bins=hour_bins)[0]
            detmean = binned_statistic(hourlist, detmean,
                                       bins=hour_bins)[0]
            detrmse = binned_statistic(hourlist, detrmse,
                                       bins=hour_bins)[0]
            meanfss = binned_statistic(hourlist, fsslist,
                                       bins=hour_bins)[0]
            np.save(savefn, (radarmean, detmean, detrmse, meanfss))

        if args.ana == 'ens':
            ensspread = binned_statistic(hourlist, spreadlist,
                                         bins=hour_bins)[0]
            ensrmse = binned_statistic(hourlist, rmselist,
                                       bins=hour_bins)[0]
            np.save(savefn, (ensspread, ensrmse))

    if args.ana == 'det':

        # Plot
        if ie == 0:
            ax1.plot(range(1, radarmean.shape[0] + 1), radarmean, c='k',
                     linewidth=3, label='Radar')
        ax1.plot(range(1, radarmean.shape[0] + 1), detmean, c=cdict[expid],
                 linewidth=2,
                 label=expid)

        # axarr[1].plot(range(radarmean.shape[0]), detrmse, c = cdict[expid],
        # linewidth = 2)
        ax2.plot(range(1, radarmean.shape[0] + 1), meanfss, c=cdict[expid],
                 linewidth=2, label=expid)
    if args.ana == 'ens':
        ax1.plot(range(1, ensspread.shape[0] + 1), ensspread, c=cdict[expid],
                 linewidth=2, linestyle='--')
        ax1.plot(range(1, ensspread.shape[0] + 1), ensrmse, c=cdict[expid],
                 linewidth=1.5, label=expid)

# End expid loop
# Finish the plots
if args.ana == 'det':
    ax1.set_ylabel('mean hourly precipitation [mm/h]')
    ax2.set_ylabel('FSS [Neighborhood size 58.8km]')
    for ax in [ax1, ax2]:
        plt.sca(ax)
        ax.set_xlabel('Time [UTC]')
        ax.legend(loc=0, fontsize=8, frameon=False)
        ymax = np.ceil(ax.get_ylim()[1] * 10) / 10.
        ax.set_ylim(0, ymax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_position(('outward', 10))
        ax.spines['bottom'].set_position(('outward', 10))
        ax.set_xticks(range(args.hint + 1)[::6])
        ax.set_xticklabels(hourlist[::6])
        ax.set_xlim(0, 24)
        ax.set_title('Deterministic forecasts ' + args.date_start[:-4] + '-' +
                     args.date_stop[:-4])
        plt.tight_layout()
if args.ana == 'ens':
    ax1.set_xlabel('Time [UTC]')
    # ax1.set_ylabel('Normalized spread / RMSE at 58.8km scale')
    ax1.set_ylabel('Normalized spread at 58.8km scale')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_position(('outward', 10))
    ax1.spines['bottom'].set_position(('outward', 10))
    ax.set_xticks(range(args.hint + 1)[::6])
    ax.set_xticklabels(hourlist[::6])
    ax1.legend(loc=0, fontsize=8, frameon=False)
    ax1.set_title('Normalized ensemble spread (--) and RMSE (-)')
    plt.tight_layout()

plotdir = plotdir + expid_str[:-1] + '/verif_fc_prec/'
if not os.path.exists(plotdir): os.makedirs(plotdir)

# Save figure and create log str
save_fig_and_log(fig1, 'diprec_' + plotstr, plotdir)
if args.ana == 'det':
    save_fig_and_log(fig2, 'fss_' + plotstr, plotdir)

plt.close('all')
