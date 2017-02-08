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

# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--expid', metavar = 'expid', type=str, nargs = '+')
parser.add_argument('--date_ana_start', metavar = 'date_ana_start', type=str)
parser.add_argument('--date_ana_stop', metavar = 'date_ana_stop', type=str)
parser.add_argument('--var', metavar = 'var', type=str, default = 'T2M')
parser.add_argument('--obs', metavar = 'obs', type=str, default = 'SYNOP')
args = parser.parse_args()

# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_ana_start)
tend = yyyymmddhhmmss_strtotime(args.date_ana_stop)
tint = timedelta(hours = 1)
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)

# Set up figure
if args.obs == 'SYNOP':
    fig, ax = plt.subplots(1, 1, figsize = (6, 5))

if args.obs in ['TEMP', 'AIREP']:
    fig, axarr = plt.subplots(1, 2, figsize = (10, 5))
    unitdict = {'T': 'K', 'RH': '%'}
    rmselimdict = {'T': (0,3), 'RH': (0, 40)}
    biaslimdict = {'T': (-3,3), 'RH': (-15,15)}
    cyc = ['b', 'g', 'r']

    if args.obs == 'TEMP':
        bin_edges = np.arange(0, 1025, 25) * 100. # Pa
        meanlev = (bin_edges[1:] + bin_edges[:-1])/2./100.  # hPa
    if args.obs == 'AIREP':
        bin_edges = np.arange(0, 12500, 250) # m
        meanlev = (bin_edges[1:] + bin_edges[:-1])/2.  # m



expid_str = ''
for ie, expid in enumerate(args.expid):
    expid_str += expid + '_'

    mean_spread = []
    mean_bias = []
    rmse = []

    spread = []
    bias = []
    lev = []

    for t in timelist:
        # Determine 3hrly storage time
        date_ana = yyyymmddhhmmss(t)
        print 'date_ana = ', date_ana
        date_fg = yyyymmddhhmmss(t - timedelta(hours = 1))
        date_store = yyyymmddhhmmss(t - timedelta(hours = t.hour % 3))

        # Directory where ekf's and fof's are stored
        FEED_DIR = ('/e/uwork/extsrasp/cosmo_letkf/feedback/' + expid +
                    '/' + date_store + '/')

        # Load files
        try:
            ekf = fdbkfile(FEED_DIR + 'ekf' + args.obs + '_' + date_ana + '.nc')
            fof = fdbkfile(FEED_DIR + 'fof_' + date_fg + '.nc')
        except:
            print 'Either ekf or fof file does not exsist. Skip over rest of loop.'
            print(FEED_DIR + 'ekf' + args.obs + '_' + date_ana + '.nc')
            print(FEED_DIR + 'fof_' + date_fg + '.nc')
            continue  # Skip over rest of loop

        ov_ekf = ekf.obs_veri(varno=args.var, obstype=args.obs, 
                                  veri = 'first guess ensemble spread')
        ov_fof = fof.obs_veri(varno=args.var, obstype=args.obs, veri = 0)

        ens_spread = ov_ekf['veri_data']
        det_bias = ov_fof['veri_data'] - ov_fof['obs']  # neg means fc colder than obs
        
        if args.obs == 'TEMP':
            # Extend list
            lev_fof = ov_fof['level']
            lev_ekf = ov_ekf['level']
            assert np.array_equal(lev_fof[lev_fof >= 5000.], lev_ekf), 'Dims do not match'
            spread.extend(list(ens_spread))
            bias.extend(list(det_bias[lev_fof >= 5000.]))
            lev.extend(list(lev_ekf))

        elif args.obs == 'AIREP':
            lev_fof = ov_fof['level']
            lev_ekf = ov_ekf['level']
            assert lev_fof.shape[0] == lev_ekf.shape[0], 'Dims do not match'
            spread.extend(list(ens_spread))
            bias.extend(list(det_bias))
            lev.extend(list(lev_ekf))

        elif args.obs == 'SYNOP':
            # Get mean values
            mean_spread.append(np.mean(ens_spread))
            mean_bias.append(np.mean(det_bias))
            rmse.append(np.sqrt(np.mean(det_bias**2)))

        else:
            raise Exception


    # Plot data for each expid
    if args.obs == 'SYNOP':
        ax.plot()

    if args.obs in ['TEMP', 'AIREP']:
        # Height bin the data
        mean_spread = binned_statistic(lev, spread, bins = bin_edges)[0]
        mean_bias = binned_statistic(lev, bias, bins = bin_edges)[0]
        rmse = np.sqrt(binned_statistic(lev, np.array(bias)**2, 
                                    bins = bin_edges)[0])
        axarr[0].plot(rmse, meanlev, c = cyc[ie], linewidth = 2)
        axarr[0].plot(mean_spread, meanlev, c = cyc[ie], linewidth = 2,
                      linestyle = '--')
        axarr[1].plot(mean_bias, meanlev, c = cyc[ie], linewidth = 2,
                      label = expid)

# Finish the plots
if args.obs == 'TEMP':
    axarr[0].set_xlim(rmselimdict[args.var])
    axarr[0].set_ylim(0,1000)
    axarr[0].set_xlabel(args.var + ' RMSE [' + unitdict[args.var] +  ']')
    axarr[0].set_ylabel('Pressure [hPa]')
    axarr[0].set_title('RMSE / Spread (dashed) ' + args.var)
    axarr[0].invert_yaxis()
    axarr[1].plot([0, 0],[0,1000],c = 'gray')
    axarr[1].set_xlim(biaslimdict[args.var])
    axarr[1].set_ylim(0,1000)
    axarr[1].set_xlabel(args.var + ' BIAS [' + unitdict[args.var] +  ']')
    axarr[1].set_ylabel('Pressure [hPa]')
    axarr[1].set_title('Bias ' + args.var)
    axarr[1].invert_yaxis()
    axarr[1].legend(loc = 1)
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])

    plotstr = (args.obs + '_' + args.var + '_' + args.date_ana_start + '_' + 
               args.date_ana_stop)
    plotdir = '/e/uwork/extsrasp/plots/' + expid_str[:-1] + '/verif_ana/'
    if not os.path.exists(plotdir): os.makedirs(plotdir)
    fig.suptitle(expid_str + '  ' + plotstr)
    print 'Save figure:', plotdir + plotstr
    fig.savefig(plotdir + plotstr)
    plt.close('all')

if args.obs == 'AIREP':
    axarr[0].set_xlim(rmselimdict[args.var])
    axarr[0].set_ylim(0,12500)
    axarr[0].set_xlabel(args.var + ' RMSE [' + unitdict[args.var] +  ']')
    axarr[0].set_ylabel('Pressure [hPa]')
    axarr[0].set_title('RMSE / Spread (dashed) ' + args.var)
    axarr[1].plot([0, 0],[0,12500],c = 'gray')
    axarr[1].set_xlim(biaslimdict[args.var])
    axarr[1].set_ylim(0,12500)
    axarr[1].set_xlabel(args.var + ' BIAS [' + unitdict[args.var] +  ']')
    axarr[1].set_ylabel('Pressure [hPa]')
    axarr[1].set_title('Bias ' + args.var)
    axarr[1].legend(loc = 1)
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])

    plotstr = (args.obs + '_' + args.var + '_' + args.date_ana_start + '_' + 
               args.date_ana_stop)
    plotdir = '/e/uwork/extsrasp/plots/' + expid_str[:-1] + '/verif_ana/'
    if not os.path.exists(plotdir): os.makedirs(plotdir)
    fig.suptitle(expid_str + '  ' + plotstr)
    print 'Save figure:', plotdir + plotstr
    fig.savefig(plotdir + plotstr)
    plt.close('all')



