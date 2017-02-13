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
                    default = 0, help = 'Start of verification interval (h from fc start)')
parser.add_argument('--ver_int_stop', metavar = 'ver_int_stop', type=int, 
                    default = 24, help = 'End of verification interval (h from fc start)')
parser.add_argument('--date_start', metavar = 'date_start', type=str,
                    default = '20160525000000', help = 'Start date for date loop (yyyymmddhhmmss)')
parser.add_argument('--date_stop', metavar = 'date_stop', type=str,
                    default = '20160610000000', help = 'End date for date loop (yyyymmddhhmmss)')
parser.add_argument('--date_inc', metavar = 'date_inc', type=ind, 
                    default = '24', help = 'Time increment between forecasts (h)')
parser.add_argument('--var', metavar = 'var', type=str, default = 'T')
parser.add_argument('--obs', metavar = 'obs', type=str, default = 'TEMP')
args = parser.parse_args()

# Config for experiment
DATA_DIR='/e/uwork/extsrasp/cosmo_letkf/data_forecast/' + args.expid
hint = 24.
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

# Config for verification
plotdir = '/e/uwork/extsrasp/plots/' + args.expid + '/verif/'
if not os.path.exists(plotdir): os.makedirs(plotdir)
savedir = '/e/uwork/extsrasp/save/' + args.expid + '/verif/'
if not os.path.exists(savedir): os.makedirs(savedir)


expid_str = ''
for ie, expid in enumerate(args.expid):
    print 'expid = ', expid
    expid_str += expid + '_'
    
    
    # Check if saved data is available
    savedir = '/e/uwork/extsrasp/save/' + expid + '/verif_ana/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    savefn = savedir + plotstr + '.npy'
    print 'Try to load pre-saved data:', savefn
    if os.path.exists(savefn):
        print 'Found pre-saved data.'
        mean_spread, mean_bias, rmse = np.load(savefn)
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
            
            if args.obs == 'TEMP':
                # Loop over obs 
                for obs_id in hdr_unique:
                    obs_date.append(t)
                    obs_time.append(ov['time'][hdr_inds == obs_id][0])
                    obs.append(ov['obs'][hdr_inds == obs_id])
                    verif.append(ov['veri_data'][hdr_inds == obs_id])
                    bias.append(ov['veri_data'][hdr_inds == obs_id] -
                                ov['obs'][hdr_inds == obs_id])
                    obs_lev.append(ov['level'][hdr_inds == obs_id])

            elif args.obs == 'SYNOP':
                obs_time.extend(list(ov['time']))
                obs.extend(list(ov['obs']))
                verif.extend(list(ov['veri_data']))
                bias.extend(list(ov['veri_data'] - ov['obs']))

            else:
                raise Exception

unitdict = {'T': 'K', 'RH': '%', 'T2M': 'K'} 
if args.obs == 'TEMP': 
    print 'Total number of observations:', len(obs)
    fig, ax = plt.subplots(1,1, figsize = (6, 4))
    hist_edges = np.arange(0, hint+1)
    tmp = ax.hist(np.aray(obs_time) / 60., bins = hist_edges, color = 'gray')
    xmax = np.max(tmp[0]) * 1.2
    ax.plot([args.ver_int_start, args.ver_int_start], [0, xmax],
            color = 'red')
    ax.plot([args.ver_int_stop, args.ver_int_stop], [0, xmax],
            color = 'red')
    ax.set_xlabel('Time UTC in h')
    ax.set_ylabel('Number of TEMP observations')
    plotstr = ('TEMP_timing_' + args.date_start + '_' + args.date_stop)
    ax.set_title(args.expid + '_' + plotstr)
    ax.set_ylim(0, xmax)
    plt.tight_layout()
    fig.savefig(plotdir + plotstr)
    plt.close('all')

    # Height bin the data
    bias = np.array(bias)
    obs_lev = np.array(obs_lev)
    obs_time = np.array(obs_time)
    mask = (obs_time >= args.ver_int_start*60.) & (obs_time <= 
                                                   args.ver_int_stop*60.)
    print 'Number of verif-observations:', np.sum(mask)

    # Flatten masked lists
    flatbias = [item for sublist in bias[mask] for item in sublist]
    flatlev = [item for sublist in obs_lev[mask] for item in sublist]
    bin_edges = np.arange(0, 1025, 25) * 100. # Pa

    mean_bias = binned_statistic(flatlev, flatbias, bins = bin_edges)[0]
    #print binned_statistic(flatlev, flatbias, bins = bin_edges, statistic = 'count')
    rmse = np.sqrt(binned_statistic(flatlev, np.array(flatbias)**2, 
                                    bins = bin_edges)[0])

    if args.var == 'RH':   # convert to percent
        mean_bias *= 100.
        rmse *= 100.

    # Plot 
    unitdict = {'T': 'K', 'RH': '%'}
    rmselimdict = {'T': (0,3), 'RH': (0, 40)}
    biaslimdict = {'T': (-3,3), 'RH': (-15,15)}
    fig, axarr = plt.subplots(1,2, figsize = (10,5))
    meanlev = (bin_edges[1:] + bin_edges[:-1])/2./100.  # hPa
    axarr[0].plot(rmse, meanlev, c = 'k', linewidth = 2)
    axarr[0].set_xlim(rmselimdict[args.var])
    axarr[0].set_ylim(0,1000)
    axarr[0].set_xlabel(args.var + ' RMSE [' + unitdict[args.var] +  ']')
    axarr[0].set_ylabel('Pressure [hPa]')
    axarr[0].set_title('RMSE ' + args.var)
    axarr[0].invert_yaxis()
    axarr[1].plot([0, 0],[0,1000],c = 'gray')
    axarr[1].plot(mean_bias, meanlev, c = 'k', linewidth = 2)
    axarr[1].set_xlim(biaslimdict[args.var])
    axarr[1].set_ylim(0,1000)
    axarr[1].set_xlabel(args.var + ' BIAS [' + unitdict[args.var] +  ']')
    axarr[1].set_ylabel('Pressure [hPa]')
    axarr[1].set_title('Bias ' + args.var)
    axarr[1].invert_yaxis()
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])

    plotstr = (args.obs + '_' + args.var + '_' + args.date_start + '_' + args.date_stop + '_'
               + str(args.ver_int_start) + '_' + str(args.ver_int_stop))
    fig.suptitle(args.expid + '  ' + plotstr)
    print 'Plotting:', plotdir + plotstr
    fig.savefig(plotdir + plotstr)
    plt.close('all')

    # Save array
    savefn = savedir + plotstr
    print 'Saving array as', savefn 
    np.save(savefn, (rmse, mean_bias, meanlev))

if args.obs == 'SYNOP':
    print 'Total number of observations:', len(obs)

    bin_edges = np.arange(0, (hint+1)*60, 60)   # Hourly bins
    time_hist = np.histogram(obs_time, bins = bin_edges)
    mean_bias = binned_statistic(obs_time, bias, bins = bin_edges)[0]
    rmse = np.sqrt(binned_statistic(obs_time, np.array(bias)**2, bins = bin_edges)[0])

    fig, ax = plt.subplots(1,1, figsize = (6,5))

    time_plot = bin_edges[:-1] / 60
    ax.plot(time_plot, mean_bias, linewidth = 2, label = 'bias', c = 'k')
    ax.plot(time_plot, rmse, linewidth = 2, label = 'rmse', c = 'k', linestyle = '--')
    ax.set_xlabel('time')
    ax.set_ylabel(args.var + ' [' + unitdict[args.var] +  ']' )
    ax.legend()
    plt.axhline(y=0, zorder = 0.1, c = 'gray')

    plt.tight_layout(rect=[0, 0.0, 1, 0.95])

    plotstr = (args.obs + '_' + args.var + '_' + args.date_start + '_' + args.date_stop)
    fig.suptitle(args.expid + '  ' + plotstr)
    print 'Plotting:', plotdir + plotstr
    fig.savefig(plotdir + plotstr)
    plt.close('all')

    # Save array
    savefn = savedir + plotstr
    print 'Saving array as', savefn 
    np.save(savefn, (rmse, mean_bias, time_plot))


        
