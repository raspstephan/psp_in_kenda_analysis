#!/usr/bin/env python
"""
Script fo analyze fof-Files and create verification plots 
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

# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--expid', metavar = 'expid', type=str, nargs = '+',
                    help = 'Experiment ID')
parser.add_argument('--ver_int_start', metavar = 'ver_int_start', type=int, 
                    default = 0, 
                    help = 'Start of verification interval (h from fc start)')
parser.add_argument('--ver_int_stop', metavar = 'ver_int_stop', type=int, 
                    default = 24, 
                    help = 'End of verification interval (h from fc start)')
parser.add_argument('--date_start', metavar = 'date_start', type=str,
                    default = '20160525000000', help = 'Start date for date loop (yyyymmddhhmmss)')
parser.add_argument('--date_stop', metavar = 'date_stop', type=str,
                    default = '20160610000000', 
                    help = 'End date for date loop (yyyymmddhhmmss)')
parser.add_argument('--date_inc', metavar = 'date_inc', type=int, 
                    default = '24', 
                    help = 'Time increment between forecasts (h)')
parser.add_argument('--var', metavar = 'var', type=str, default = 'T')
parser.add_argument('--obs', metavar = 'obs', type=str, default = 'TEMP')
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
    DATA_DIR='/e/uwork/extsrasp/cosmo_letkf/data_forecast/' + expid
    
    
    # Check if saved data is available ATTENTION NOT IMPLEMENTED
    #savedir = '/e/uwork/extsrasp/save/' + expid + '/verif_ana/'
    #if not os.path.exists(savedir): os.makedirs(savedir)
    #savefn = savedir + plotstr + '.npy'
    #print 'Try to load pre-saved data:', savefn
    #if os.path.exists(savefn):
        #print 'Found pre-saved data.'
        #mean_spread, mean_bias, rmse = np.load(savefn)
    #else:
    if True:
        obs_time = []   # in mins relative to obs_date
        obs = []
        verif = []
        bias = []
        obs_date = []   # datetime object
        obs_lev = []    # in Pa
        for t in timelist:
            print t
            # Load fof file
            foffn = (DATA_DIR + '/' + yyyymmddhhmmss(t) + '/det/fof_' +  
                    yyyymmddhhmmss(t) + '.nc')
            fof = fdbkfile(foffn)
            var_tab = fof.table_varno
            obs_tab = fof.table_obstype
            varno = var_tab[args.var]
            obsno = obs_tab[args.obs]
            if fof.fb_times()[-1]/60. > 24.:
                hint = 48
            
            # Get temp data, T
            ov = fof.obs_veri(varno=varno, obstype=obsno)
            hdr_inds = ov['hdr_inds']
            hdr_unique = np.unique(hdr_inds)  # Get unique obs labels 
            
            if args.obs in ['TEMP', 'AIREP']:
                obs_time.extend(ov['time'])
                obs.extend(ov['obs'])
                verif.extend(ov['veri_data'])
                bias.extend(ov['veri_data'] - ov['obs'])
                obs_lev.extend(ov['level'])

            elif args.obs == 'SYNOP':
                obs_time.extend(list(ov['time']))
                obs.extend(list(ov['obs']))
                verif.extend(list(ov['veri_data']))
                bias.extend(list(ov['veri_data'] - ov['obs']))

            else:
                raise Exception
        # End timeloop
    
    # Plot data for each expid 

    if args.obs in ['TEMP', 'AIREP']: 
        print 'Total number of observations:', len(obs)

        # Height bin the data
        bias = np.array(bias)
        obs_lev = np.array(obs_lev)
        obs_time = np.array(obs_time)
        obs_lev = np.array(obs_lev)
        mask = (obs_time >= args.ver_int_start*60.) & (obs_time <= 
                                                    args.ver_int_stop*60.)
        print 'Number of verif-observations:', np.sum(mask)

        mean_bias = binned_statistic(obs_lev[mask], bias[mask], 
                                     bins = bin_edges)[0]
        rmse = np.sqrt(binned_statistic(obs_lev[mask], np.array(bias[mask])**2, 
                                        bins = bin_edges)[0])
        if args.var == 'RH':   # convert to percent
            mean_bias *= 100.
            rmse *= 100.
        # Plot 
        axarr[0].plot(rmse, meanlev, c = cdict[expid], linewidth = 2)
        axarr[1].plot(mean_bias, meanlev, c = cdict[expid], linewidth = 2,
                      label = expid)
        
    if args.obs == 'SYNOP':
        print 'Total number of observations:', len(obs)

        bin_edges = np.arange(0, (args.hint+1)*60, 60)   # Hourly bins
        time_hist = np.histogram(obs_time, bins = bin_edges)
        mean_bias = binned_statistic(obs_time, bias, bins = bin_edges)[0]
        rmse = np.sqrt(binned_statistic(obs_time, np.array(bias)**2, 
                                        bins = bin_edges)[0])
        time_plot = bin_edges[:-1] / 60
        ax.plot(time_plot, mean_bias, linewidth = 2, c = cdict[expid], 
                label = expid, linestyle = ':')
        ax.plot(time_plot, rmse, linewidth = 2, c = cdict[expid], 
                linestyle = '-')
        
# Finish the plots
plotdir = '/e/uwork/extsrasp/plots/' + expid_str[:-1] + '/verif_fof/'
if not os.path.exists(plotdir): os.makedirs(plotdir)
if args.obs in ['TEMP', 'AIREP']:
    axarr[0].set_xlim(rmselimdict[args.var])
    axarr[0].set_ylim(0,np.max(meanlev))
    axarr[0].set_xlabel(args.var + ' RMSE [' + unitdict[args.var] +  ']')
    axarr[0].set_ylabel('Pressure [hPa]')
    axarr[0].set_title('RMSE ' + args.var)
    axarr[1].axvline(0 ,c = 'gray')
    axarr[1].set_xlim(biaslimdict[args.var])
    axarr[1].set_ylim(0,np.max(meanlev))
    axarr[1].set_xlabel(args.var + ' BIAS [' + unitdict[args.var] +  ']')
    axarr[1].set_ylabel('Pressure [hPa]')
    axarr[1].set_title('Bias ' + args.var)
    axarr[1].legend(loc = 0, fontsize = 8)
    if args.obs == 'TEMP':
        axarr[0].invert_yaxis()
        axarr[1].invert_yaxis()
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    
    plotstr = (args.obs + '_' + args.var + '_' + args.date_start + '_' + 
               args.date_stop + '_' + str(args.ver_int_start) + '_' + 
               str(args.ver_int_stop))
    fig.suptitle(plotstr)
    print 'Plotting:', plotdir + plotstr
    fig.savefig(plotdir + plotstr, dpi = 200)
    plt.close('all')


if args.obs == 'SYNOP':
    ax.set_xlabel('time')
    ax.set_ylabel(args.var + ' [' + unitdict[args.var] +  ']' )
    ax.legend(loc = 0, fontsize = 8)
    ax.axhline(y=0, zorder = 0.1, c = 'gray')
    plotstr = (args.obs + '_' + args.var + '_' + args.date_start + '_' + 
               args.date_stop)
    ax.set_title(plotstr + '\n RMSE (solid), Bias (dotted)')
    plt.tight_layout()
    print 'Plotting:', plotdir + plotstr
    fig.savefig(plotdir + plotstr, dpi = 200)
    plt.close('all')


        
