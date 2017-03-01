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
from matplotlib import gridspec

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

# General settings
if os.getcwd() == '/panfs/e/vol0/extsrasp/dwd_scripts':
    plotdir = '/e/uwork/extsrasp/plots/'
    datadir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/'
    radardir = '/e/uwork/extsrasp/radolan/'
    savedir_base = '/e/uwork/extsrasp/save/'
elif os.getcwd() == '/home/s/S.Rasp/repositories/dwd_scripts':
    datadir_cosmo = '/home/cosmo/stephan.rasp/dwd_data/data_forecast/'
    datadir_raid = '/home/data/raid_linux/stephan.rasp/dwd_data/data_forecast/'
    datadirdict = {'DA_REF': datadir_cosmo,
                   'DA_PSPv2': datadir_raid,
                   'DA_PSPv2_TL500': datadir_cosmo,
                   'DA_REF_TL500': datadir_raid,
                   }
    radardir = '/project/meteo/w2w/A6/radolan/netcdf_cosmo_de/'
    plotdir = '/home/s/S.Rasp/dwd_plots/plots/'
    savedir_base = '/home/cosmo/stephan.rasp/dwd_data/save/'
else: 
    raise Exception('Working directory not recognized:' + os.getcwd())


# Config for experiment
cdict = {'radar':'k',
             'REF':'navy',
             'REF_TL500':'darkgreen',
             'PSP_TL500':'orange',
            'DA_REF':'navy',
            'DA_REF_TL500':'cyan',
            'DA_PSP_TL500':'red',
            'DA_PSPv2_TL500':'fuchsia',
            'DA_PSPv2':'maroon',
             }

if args.obs in ['TEMP', 'AIREP']:
    plotstr = (args.obs + '_' + args.var + '_' + args.date_start + '_' + 
               args.date_stop + '_' + str(args.ver_int_start) + '_' + 
               str(args.ver_int_stop))
if args.obs == 'SYNOP':
    plotstr = (args.obs + '_' + args.var + '_' + args.date_start + '_' + 
               args.date_stop)


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
    #fig1, axarr1 = plt.subplots(1, 2, figsize = (5, 7))
    fig1 = plt.figure(figsize = (5,5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 4]) 
    axarr1 = [plt.subplot(gs[0]), plt.subplot(gs[1])]
    fig2 = plt.figure(figsize = (5,5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 4]) 
    axarr2 = [plt.subplot(gs[0]), plt.subplot(gs[1])]
    #fig2, axarr2 = plt.subplots(1, 2, figsize = (5, 7))
    unitdict = {'T': 'K', 'RH': '%'}
    rmselimdict = {'T': (0,3), 'RH': (0, 40)}
    biaslimdict = {'T': (-1.5,1.5), 'RH': (-15,15)}

    if args.obs == 'TEMP':
        binwidth = 25
        bin_edges = np.arange(200, 1000+binwidth, binwidth) * 100. # Pa
        meanlev = (bin_edges[1:] + bin_edges[:-1])/2./100.  # hPa
    if args.obs == 'AIREP':
        binwidth = 250
        bin_edges = np.arange(0, 10000 + binwidth, binwidth) # m
        meanlev = (bin_edges[1:] + bin_edges[:-1])/2.  # m


expid_str = ''
for ie, expid in enumerate(args.expid):
    print 'expid = ', expid
    expid_str += expid + '_'
    DATA_DIR = datadirdict[expid] + expid
    
    
    # Check if saved data is available
    savedir = savedir_base + expid + '/verif_fc/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    savefn = savedir + plotstr + '.npy'
    print 'Try to load pre-saved data:', savefn
    if os.path.exists(savefn):
        print 'Found pre-saved data.'
        mean_bias, rmse, count = np.load(savefn)
    else:
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
    
        # Save data 
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
            count = binned_statistic(obs_lev[mask], bias[mask], 
                                        bins = bin_edges, 
                                        statistic = 'count')[0]
            
            mean_bias = binned_statistic(obs_lev[mask], bias[mask], 
                                        bins = bin_edges)[0]
            rmse = np.sqrt(binned_statistic(obs_lev[mask], np.array(bias[mask])**2, 
                                            bins = bin_edges)[0])
            if args.var == 'RH':   # convert to percent
                mean_bias *= 100.
                rmse *= 100.
            np.save(savefn, (mean_bias, rmse, count))
            
        if args.obs == 'SYNOP':
            print 'Total number of observations:', len(obs)

            bin_edges = np.arange(0, (args.hint+1)*60, 60)   # Hourly bins
            time_hist = np.histogram(obs_time, bins = bin_edges)
            mean_bias = binned_statistic(obs_time, bias, bins = bin_edges)[0]
            rmse = np.sqrt(binned_statistic(obs_time, np.array(bias)**2, 
                                            bins = bin_edges)[0])
            np.save(savefn, (mean_bias, rmse, None))
    
    
    # Plot data for each expid 

    if args.obs in ['TEMP', 'AIREP']: 
        # Plot 
        axarr1[0].barh(meanlev, count, height = binwidth, color = 'gray', 
                       linewidth = 0)
        axarr2[0].barh(meanlev, count, height = binwidth, color = 'gray', 
                       linewidth = 0)
        axarr1[1].plot(rmse, meanlev, c = cdict[expid], linewidth = 2,
                 label = expid)
        axarr2[1].plot(mean_bias, meanlev, c = cdict[expid], linewidth = 2,
                      label = expid)
        
    if args.obs == 'SYNOP':
        time_plot = bin_edges[:-1] / 60
        ax.plot(time_plot, mean_bias, linewidth = 2, c = cdict[expid], 
                label = expid, linestyle = ':')
        ax.plot(time_plot, rmse, linewidth = 2, c = cdict[expid], 
                linestyle = '-')
        
# Finish the plots
plotdir = plotdir + expid_str[:-1] + '/verif_fof/'
if not os.path.exists(plotdir): os.makedirs(plotdir)
if args.obs in ['TEMP', 'AIREP']:
    
    
    axarr1[1].set_xlim(rmselimdict[args.var])
    axarr1[1].set_ylim(0,np.max(meanlev))
    axarr1[1].set_xlabel(args.var + ' RMSE [' + unitdict[args.var] +  ']')
    #axarr1[1].set_title('RMSE ' + args.var)
    axarr1[1].legend(loc = 0, fontsize = 8, frameon = False)
    
    axarr2[1].axvline(0 ,c = 'lightgray', zorder = 0.1)
    axarr2[1].set_xlim(biaslimdict[args.var])
    axarr2[1].set_ylim(0,np.max(meanlev))
    axarr2[1].set_xlabel(args.var + ' BIAS [' + unitdict[args.var] +  ']')
    #axarr2[1].set_title('Bias ' + args.var)
    axarr2[1].legend(loc = 0, fontsize = 8, frameon = False)
    
    axarr1[0].set_xlabel('# Obs')
    axarr2[0].set_xlabel('# Obs')
    #axarr1[1].set_title('RMSE')
    #axarr2[1].set_title('BIAS')
    
    for axarr in [axarr1, axarr2]:
        plt.sca(axarr[0])
        axarr[0].spines['top'].set_visible(False)
        axarr[0].spines['right'].set_visible(False)
        axarr[1].spines['right'].set_visible(False)
        axarr[1].spines['left'].set_visible(False)
        axarr[1].spines['top'].set_visible(False)
        axarr[0].spines['left'].set_position(('outward', 10))
        axarr[0].spines['bottom'].set_position(('outward', 10))
        axarr[1].spines['bottom'].set_position(('outward', 10))
        axarr[1].yaxis.set_ticks([])
        
        plt.tight_layout(rect=[0, 0.0, 1, 0.97])
    
    if args.obs == 'TEMP':
        axarr1[0].set_ylabel('Pressure [hPa]')
        axarr2[0].set_ylabel('Pressure [hPa]')
        axarr1[0].set_ylim(200, 1000)
        axarr1[1].set_ylim(200, 1000)
        axarr2[0].set_ylim(200, 1000)
        axarr2[1].set_ylim(200, 1000)
        axarr1[0].invert_yaxis()
        axarr1[1].invert_yaxis()
        axarr2[0].invert_yaxis()
        axarr2[1].invert_yaxis()
    if args.obs == 'AIREP':
        axarr1[0].set_ylabel('Height agu [m]')
        axarr2[0].set_ylabel('Height agu [m]')
        axarr1[0].set_ylim(0, 10000)
        axarr1[1].set_ylim(0, 10000)
        axarr2[0].set_ylim(0, 10000)
        axarr2[1].set_ylim(0, 10000)
        
    for axarr in [axarr1, axarr2]:
        plt.sca(axarr[0])
        plt.tight_layout(rect=[0, 0.0, 1, 0.97])
        
    
    
    fig1.suptitle(args.obs + ' ' + args.var + ' RMSE', fontsize = 10)
    fig2.suptitle(args.obs + ' ' + args.var + ' BIAS', fontsize = 10)
    print 'Plotting:', plotdir + 'rmse_' + plotstr + '.pdf'
    #fig1.savefig(plotdir + 'rmse_' + plotstr + '.pdf', format = 'pdf',
                 #transparent = True)
    #fig2.savefig(plotdir + 'bias_' + plotstr + '.pdf', format = 'pdf',
                 #transparent = True)
    fig1.savefig(plotdir + 'rmse_' + plotstr + '.png', format = 'png',
                 transparent = True, dpi = 300)
    fig2.savefig(plotdir + 'bias_' + plotstr + '.png', format = 'png',
                 transparent = True, dpi = 300)
    plt.close('all')


if args.obs == 'SYNOP':
    ax.set_xlabel('time')
    ax.set_ylabel(args.var + ' [' + unitdict[args.var] +  ']' )
    ax.legend(loc = 0, fontsize = 8)
    ax.axhline(y=0, zorder = 0.1, c = 'gray')

    ax.set_title(plotstr + '\n RMSE (solid), Bias (dotted)')
    plt.tight_layout()
    print 'Plotting:', plotdir + plotstr
    fig.savefig(plotdir + plotstr, dpi = 200)
    plt.close('all')


        
