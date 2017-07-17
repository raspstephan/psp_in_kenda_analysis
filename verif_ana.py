#!/usr/bin/env python
"""
Verify ekf and fof files from analysis
"""
# Imports
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cosmo_utils.da import fdbkfile
from cosmo_utils.helpers import yyyymmddhhmmss, yyyymmddhhmmss_strtotime, \
                                make_timelist
from datetime import timedelta
from scipy.stats import binned_statistic

from config import *

# Arguments
parser = argparse.ArgumentParser(description = 'Script to verify analysis against SYNOP, TEMP and AIREP observations.')
parser.add_argument('--expid', metavar = 'expid', type=str, nargs = '+',
                    help = 'Experiment ID')
parser.add_argument('--date_ana_start', metavar = 'date_ana_start', type=str,
                    help = 'Start date of verification (yyyymmddhhmmss). In case of SYNOP these determine the start and end of the plot which contains verifications for every hour. In case of TEMP or AIREP these determine the start and end of the averaging period which is then displayed as one plot.')
parser.add_argument('--date_ana_stop', metavar = 'date_ana_stop', type=str,
                    help = 'End date of verification (yyyymmddhhmmss)')
parser.add_argument('--ver_int_start', metavar = 'ver_int_start', type=int, 
                    default = 0, 
                    help = 'Time of day (h/UTC) when verification starts, including')
parser.add_argument('--ver_int_stop', metavar = 'ver_int_stop', type=int, 
                    default = 23, 
                    help = 'Time of day (h/UTC) when verification stops, including')
parser.add_argument('--var', metavar = 'var', type=str, default = 'T2M',
                    help = 'Variable to be verified.')
parser.add_argument('--obs', metavar = 'obs', type=str, default = 'SYNOP',
                    help = 'Observation type: SYNOP, TEMP or AIREP')
parser.add_argument('--state', metavar = 'state', type=str, default = 'all',
                    help = 'Obs state: all (default), active or active_passive')
parser.add_argument('--composite', metavar = 'composite', type=str, 
                    default = 'False',
                    help = 'If True diurnal composite is plotted for SYNOP')
parser.add_argument('--delete_save', metavar = 'delete_save', type=str, 
                    default = 'False',
                    help = 'If True delete old save file')
args = parser.parse_args()

plotstr = (args.obs + '_' + args.var + '_' + args.date_ana_start + '_' + 
           args.date_ana_stop + '_' + str(args.ver_int_start).zfill(2) + 
           '_' + str(args.ver_int_stop).zfill(2) + '_' + args.state)
if args.composite == 'True':
    plotstr += '_composite'

tstart = yyyymmddhhmmss_strtotime(args.date_ana_start)
tend = yyyymmddhhmmss_strtotime(args.date_ana_stop)
tint = timedelta(hours = 1)
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)
    
# Filter out times which are not in ver_int 
new_timelist = []
for t in timelist:
    if t.hour >= args.ver_int_start and t.hour <= args.ver_int_stop:
        new_timelist.append(t)
timelist = new_timelist



# Set up figure
if args.obs == 'SYNOP':
    fig, ax = plt.subplots(1, 1, figsize = (7, 5))
    unitdict = {'T2M': 'K', 'RH2M': '%', 'PS': 'Pa'}

if args.obs in ['TEMP', 'AIREP']:
    fig, axarr = plt.subplots(1, 2, figsize = (10, 5))
    unitdict = {'T': 'K', 'RH': '%'}
    rmselimdict = {'T': (0,3), 'RH': (0, 40)}
    biaslimdict = {'T': (-3,3), 'RH': (-15,15)}

    if args.obs == 'TEMP':
        bin_edges = np.arange(0, 1025, 25) * 100. # Pa
        meanlev = (bin_edges[1:] + bin_edges[:-1])/2./100.  # hPa
    if args.obs == 'AIREP':
        bin_edges = np.arange(0, 12500, 250) # m
        meanlev = (bin_edges[1:] + bin_edges[:-1])/2.  # m



expid_str = ''
for ie, expid in enumerate(args.expid):
    print 'expid = ', expid
    expid_str += expid + '_'
    
    
    # Check if saved data is available
    savedir = savedir_base + expid + '/verif_ana/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    savefn = savedir + plotstr + '.npy'
    if args.delete_save == 'True':
        print 'delete old save file:', savefn
        os.system('rm ' + savefn)
    print 'Try to load pre-saved data:', savefn
    if os.path.exists(savefn):
        print 'Found pre-saved data.'
        mean_spread, mean_bias, rmse = np.load(savefn)
    else:
        print 'Did not find saved data'
        mean_spread = []
        mean_bias = []
        rmse = []
        spread = []
        bias = []
        lev = []
        hourlist = []
        for t in timelist:
            # Determine 3hrly storage time
            date_ana = yyyymmddhhmmss(t)
            print 'date_ana = ', date_ana
            date_fg = yyyymmddhhmmss(t - timedelta(hours = 1))
            date_store = yyyymmddhhmmss(t - timedelta(hours = t.hour % 3))

            # Directory where ekf's and fof's are stored
            FEED_DIR = (feeddir + expid +
                        '/' + date_store + '/')

            # Load files
            try:
                ekf = fdbkfile(FEED_DIR + 'ekf' + args.obs + '_' + date_ana + 
                               '.nc')
                fof = fdbkfile(FEED_DIR + 'fof_' + date_fg + '.nc')
            except:
                print 'Either ekf or fof file does not exsist. Skip over rest of loop.'
                print(FEED_DIR + 'ekf' + args.obs + '_' + date_ana + '.nc')
                print(FEED_DIR + 'fof_' + date_fg + '.nc')
                continue  # Skip over rest of loop

            ov_ekf = ekf.obs_veri(varno=args.var, obstype=args.obs, 
                                    veri = 'first guess ensemble spread',
                                    )
            ov_fof = fof.obs_veri(varno=args.var, obstype=args.obs, veri = 0)

            ens_spread = ov_ekf['veri_data']
            det_bias = ov_fof['veri_data'] - ov_fof['obs']  # neg means fc colder than obs
            print 'Total number of observations:', ov_ekf['state'].shape[0]
            if args.state == 'active_passive':
                active = [ov_ekf['state']<=5]   # Active and passive
            elif args.state == 'active':
                active = [ov_ekf['state']==1]
            else:
                active = [ov_ekf['state']>=0]   # all
            print 'Number of obs used for verification:', np.sum(active)
            if args.obs == 'TEMP':
                # Extend list
                lev_fof = ov_fof['level']
                lev_ekf = ov_ekf['level']
                
                # Make fof's and ekf's match
                match = np.zeros(lev_fof.shape, dtype = bool)
                l = 0
                for k in range(len(lev_fof)):
                    if l < len(lev_ekf):
                        if lev_fof[k] == lev_ekf[l]:
                            match[k] = True
                            l += 1
                assert np.array_equal(lev_fof[match], lev_ekf), 'Dims do not match'
                spread.extend(list(ens_spread[active]))
                bias.extend(list(det_bias[lev_fof >= 5000.][active]))
                lev.extend(list(lev_ekf[active]))

            elif args.obs == 'AIREP':
                lev_fof = ov_fof['level']
                lev_ekf = ov_ekf['level']
                assert lev_fof.shape[0] == lev_ekf.shape[0], 'Dims do not match'
                spread.extend(list(ens_spread[active]))
                bias.extend(list(det_bias[active]))
                lev.extend(list(lev_ekf[active]))

            elif args.obs == 'SYNOP':
                # Get mean values
                hourlist.append(t.hour)  # From 0 to 23
                mean_spread.append(np.mean(ens_spread[active]))
                mean_bias.append(np.mean(det_bias[active]))
                rmse.append(np.sqrt(np.mean(det_bias[active]**2)))

            else:
                raise Exception
        # End timeloop

        # Plot data for each expid
        if args.obs == 'SYNOP':
            mean_spread = np.array(mean_spread)
            mean_bias = np.array(mean_bias)
            rmse = np.array(rmse)
            if args.composite == 'True':
                hour_bins = np.arange(-0.5, 24.5, 1)
                mean_spread = binned_statistic(hourlist, mean_spread, 
                                               bins = hour_bins)[0]
                mean_bias = binned_statistic(hourlist, mean_bias, 
                                               bins = hour_bins)[0]
                rmse = binned_statistic(hourlist, rmse, 
                                               bins = hour_bins)[0]
                hourlist = binned_statistic(hourlist, hourlist, 
                                               bins = hour_bins)[0]
            np.save(savefn, (mean_spread, mean_bias, rmse))

        if args.obs in ['TEMP', 'AIREP']:
            if args.var == 'RH':
                spread = np.array(spread) * 100.   # %
                bias = np.array(bias) * 100.   # %
            # Height bin the data
            mean_spread = binned_statistic(lev, spread, bins = bin_edges)[0]
            mean_bias = binned_statistic(lev, bias, bins = bin_edges)[0]
            rmse = np.sqrt(binned_statistic(lev, np.array(bias)**2, 
                                        bins = bin_edges)[0])
            np.save(savefn, (mean_spread, mean_bias, rmse))
        
        # Plot data for each expid
    if args.obs == 'SYNOP':
        ax.plot(range(mean_spread.shape[0]), mean_spread, c = cdict[expid], 
                linewidth = 2, linestyle = '--')
        ax.plot(range(mean_spread.shape[0]), rmse, c = cdict[expid], 
                linewidth = 2, linestyle = '-', label = expid)
        ax.plot(range(mean_spread.shape[0]), mean_bias, c = cdict[expid], 
                linewidth = 2, linestyle = ':')

    if args.obs in ['TEMP', 'AIREP']:
        axarr[0].plot(rmse, meanlev, c = cdict[expid], linewidth = 2)
        axarr[0].plot(mean_spread, meanlev, c = cdict[expid], linewidth = 2,
                      linestyle = '--')
        axarr[1].plot(mean_bias, meanlev, c = cdict[expid], linewidth = 2,
                      label = expid)

# Finish the plots
if args.obs == 'SYNOP':
    if args.composite == 'True':
        ax.set_xlabel('Time [UTC/h]')
    else:
        ax.set_xlabel('time from ' +  yyyymmddhhmmss(tstart)+ ' [h]')
    ax.set_ylabel('[' + unitdict[args.var] +  ']')
    ax.axhline(0, c = 'gray', zorder = 0.1)
    ax.legend(loc = 0, fontsize = 6)
    ax.set_title(plotstr + '\n RMSE (solid), Spread (dashed), Bias (dotted)')
    plt.tight_layout()
    plotdir = plotdir + expid_str[:-1] + '/verif_ana/'
    if not os.path.exists(plotdir): os.makedirs(plotdir)
    print 'Save figure:', plotdir + plotstr
    fig.savefig(plotdir + plotstr, dpi = 300)
    plt.close('all')

if args.obs == 'TEMP':
    axarr[0].set_xlim(rmselimdict[args.var])
    axarr[0].set_ylim(0,1000)
    axarr[0].set_xlabel(args.var + ' RMSE [' + unitdict[args.var] +  ']')
    axarr[0].set_ylabel('Pressure [hPa]')
    axarr[0].set_title('RMSE (solid) / Spread (dashed) ' + args.var)
    axarr[0].invert_yaxis()
    axarr[1].plot([0, 0],[0,1000],c = 'gray')
    axarr[1].set_xlim(biaslimdict[args.var])
    axarr[1].set_ylim(0,1000)
    axarr[1].set_xlabel(args.var + ' BIAS [' + unitdict[args.var] +  ']')
    axarr[1].set_ylabel('Pressure [hPa]')
    axarr[1].set_title('Bias ' + args.var)
    axarr[1].invert_yaxis()
    axarr[1].legend(loc = 0, fontsize = 6)
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    plotdir = plotdir + expid_str[:-1] + '/verif_ana/'
    if not os.path.exists(plotdir): os.makedirs(plotdir)
    fig.suptitle(expid_str + '  ' + plotstr)
    print 'Save figure:', plotdir + plotstr
    fig.savefig(plotdir + plotstr, dpi = 300)
    plt.close('all')

if args.obs == 'AIREP':
    axarr[0].set_xlim(rmselimdict[args.var])
    axarr[0].set_ylim(0,12500)
    axarr[0].set_xlabel(args.var + ' RMSE [' + unitdict[args.var] +  ']')
    axarr[0].set_ylabel('Height [m]')
    axarr[0].set_title('RMSE (solid) / Spread (dashed) ' + args.var)
    axarr[1].plot([0, 0],[0,12500],c = 'gray')
    axarr[1].set_xlim(biaslimdict[args.var])
    axarr[1].set_ylim(0,12500)
    axarr[1].set_xlabel(args.var + ' BIAS [' + unitdict[args.var] +  ']')
    axarr[1].set_ylabel('Height [m]')
    axarr[1].set_title('Bias ' + args.var)
    axarr[1].legend(loc = 0, fontsize = 6)
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    plotdir = plotdir + expid_str[:-1] + '/verif_ana/'
    if not os.path.exists(plotdir): os.makedirs(plotdir)
    fig.suptitle(expid_str + '  ' + plotstr)
    print 'Save figure:', plotdir + plotstr
    fig.savefig(plotdir + plotstr, dpi = 300)
    plt.close('all')



