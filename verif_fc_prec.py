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
    yyyymmddhhmmss, ddhhmmss, ddhhmmss_strtotime, yymmddhhmm
from cosmo_utils.pyncdf import getfobj_ncdf
from cosmo_utils.pywgrib import getfobj_ens, getfobj
from scipy.stats import binned_statistic
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
from cosmo_utils.scores.probab import FSS


# General settings
if os.getcwd() == '/panfs/e/vol0/extsrasp/dwd_scripts':
    plotdir = '/e/uwork/extsrasp/plots/'
    datadir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/'
    radardir = '/e/uwork/extsrasp/radolan/'
    savedir_base = '/e/uwork/extsrasp/save/'
elif os.getcwd() == '/home/s/S.Rasp/repositories/dwd_scripts':
    datadir = '/home/cosmo/stephan.rasp/dwd_data/data_forecast/'
    radardir = '/project/meteo/w2w/A6/radolan/netcdf_cosmo_de/'
    plotdir = '/home/s/S.Rasp/dwd_plots/plots/'
    savedir_base = '/home/cosmo/stephan.rasp/dwd_data/save/'
else: 
    raise Exception('Working directory not recognized:' + os.getcwd())

radarpref = 'raa01-rw_10000-'
radarsufx = '-dwd---bin.nc'
precsuf = '_15'
gribpref = 'lfff'
nens = 20
# Load radar tot mask

totmask = np.load('./radar_tot_mask.npy')

# Define loading functions
def load_det(datadir, date, t):
    topdir = datadir + '/' + date + '/'
    gribfn = gribpref + t + precsuf
    detfn = topdir + 'det/' + gribfn
    detfobj = getfobj(detfn, fieldn = 'PREC_PERHOUR')
    # Minus one hour
    #gribfnm1 = gribpref + ddhhmmss(ddhhmmss_strtotime(t) - 
                                 #timedelta(hours = 1)) + precsuf
    #detfnm1 = topdir + 'det/' + gribfnm1
    #detfobjm1 = getfobj(detfnm1, fieldn = 'TOT_PREC')
    #detfobj.data = detfobj.data - detfobjm1.data
    return detfobj

def load_radar(date, t):
    dateobj = (yyyymmddhhmmss_strtotime(date) + 
            ddhhmmss_strtotime(t))
    radardt = timedelta(minutes = 10)   # TODO Is this correct???
    radardateobj = dateobj - radardt
    radarfn = radardir + radarpref + yymmddhhmm(radardateobj) + radarsufx
    radarfobj = getfobj_ncdf(radarfn, fieldn = 'pr', dwdradar = True)
    return radarfobj

def load_ens(datadir, date, t):
    topdir = datadir + '/' + date + '/'
    gribfn = gribpref + t + precsuf
    ensfobjlist = getfobj_ens(topdir, 'sub', mems = nens, gribfn = gribfn, 
                              dir_prefix = 'ens', fieldn = 'PREC_PERHOUR', 
                              para = 4)
    return ensfobjlist


# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--expid', metavar = 'expid', type=str, nargs = '+',
                    help = 'Experiment ID')
parser.add_argument('--date_start', metavar = 'date_start', type=str,
                    default = '20160525000000', help = 'Start date for date loop (yyyymmddhhmmss)')
parser.add_argument('--date_stop', metavar = 'date_stop', type=str,
                    default = '20160610000000', 
                    help = 'End date for date loop (yyyymmddhhmmss)')
parser.add_argument('--date_inc', metavar = 'date_inc', type=int, 
                    default = '24', 
                    help = 'Time increment between forecasts (h)')
parser.add_argument('--hint', metavar = 'hint', type=int, default =24,
                    help = 'Maximum forecast lead time')
parser.add_argument('--ana', metavar = 'ana', type=str,
                    help = 'Type of analysis to be done [di]')
args = parser.parse_args()


plotstr = (args.ana + '_' + args.date_start + '_' + 
           args.date_stop)

# Config for experiment
cdict = {'radar':'k',
             'REF':'navy',
             'REF_TL500':'darkgreen',
             'PSP_TL500':'orange',
            'DA_REF':'navy',
            'DA_REF_ens':'blue',
            'DA_REF_TL500':'cyan',
            'DA_PSP_TL500':'red',
            'DA_PSPv2_TL500':'fuchsia',
            'DA_PSPv2_TL500_ens':'magenta',
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

n = 21
kernel = np.ones((n,n))/float((n*n))

# Set up figure
fig1, ax1 = plt.subplots(1, 1, figsize = (6, 4))
fig2, ax2 = plt.subplots(1, 1, figsize = (6, 4))

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
        
        hourlist =[]
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
            for h in range(1, args.hint+1):
                print 'Hour', h
                hourlist.append(h)
                # Load data 
                radarfobj = load_radar(date, ddhhmmss(timedelta(hours=h)))
                nanradar = radarfobj.data
                nanradar[totmask] = np.nan
                convradar = convolve2d(nanradar, kernel, mode='same')
                if args.ana == 'det':
                    # Det run 
                    detfobj = load_det(DATA_DIR, date, ddhhmmss(timedelta(hours=h)))
                    
                    radarmean.append(np.mean(radarfobj.data[~totmask]))
                    detmean.append(np.mean(detfobj.data[~totmask]))
                    
                    # Upscaled RMSE
                    nandet = detfobj.data
                    nandet[totmask] = np.nan 
                    convdet = convolve2d(nandet, kernel, mode='same')
                    detradarmean = 0.5*(convradar + convdet)
                    detrmse.append(np.sqrt(np.nanmean(((convradar - convdet)**2 / 
                                        (detradarmean)**2)[detradarmean >= 0.1])))
                    fss = FSS(1, 21, nanradar, nandet, python_core = True)
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
                    meanfield = np.mean(convfieldlist, axis = 0)
                    spreadlist.append(np.nanmean((np.std(convfieldlist, axis = 0)/
                                                meanfield)[meanfield >= 0.1]))
                    ensradarmean = 0.5*(convradar + meanfield)
                    rmselist.append(np.sqrt(np.nanmean(((convradar - meanfield)**2 / 
                                        (ensradarmean)**2)[ensradarmean >= 0.1])))
                    
                    
                    
            # End hourlist 
        # End timelist
        
        # Bin the data 
        hour_bins = np.arange(0.5, 25.5, 1)
        if args.ana == 'det':
            radarmean = binned_statistic(hourlist, radarmean, 
                                            bins = hour_bins)[0]
            detmean = binned_statistic(hourlist, detmean, 
                                            bins = hour_bins)[0]
            detrmse = binned_statistic(hourlist, detrmse, 
                                            bins = hour_bins)[0]
            meanfss = binned_statistic(hourlist, fsslist, 
                                            bins = hour_bins)[0]
            np.save(savefn, (radarmean, detmean, detrmse, meanfss))

        if args.ana == 'ens':
            ensspread = binned_statistic(hourlist, spreadlist, 
                                            bins = hour_bins)[0]
            ensrmse = binned_statistic(hourlist, rmselist, 
                                            bins = hour_bins)[0]
            np.save(savefn, (ensspread, ensrmse))


    
    if args.ana == 'det':

        # Plot
        if ie == 0:
            ax1.plot(range(radarmean.shape[0]), radarmean, c = 'k', 
                        linewidth = 3, label = 'Radar')
        ax1.plot(range(radarmean.shape[0]), detmean, c = cdict[expid], 
                    linewidth = 2,
                    label = expid)
        
        #axarr[1].plot(range(radarmean.shape[0]), detrmse, c = cdict[expid], 
                     #linewidth = 2)
        ax2.plot(range(radarmean.shape[0]), meanfss, c = cdict[expid], 
                    linewidth = 2, label = expid)
    if args.ana == 'ens':

        axarr[0].plot(range(ensspread.shape[0]), ensspread, c = cdict[expid], 
                    linewidth = 2, linestyle = '--')
        axarr[0].plot(range(ensspread.shape[0]), ensrmse, c = cdict[expid], 
                    linewidth = 2, label = expid)
    
    

# End expid loop
# Finish the plots
if args.ana == 'det':
    ax1.set_ylabel('mean hourly precipitation [mm/h]')
    ax2.set_ylabel('FSS [Neighborhood size 58.8km]')
    for ax in [ax1, ax2]:
        plt.sca(ax)
        ax.set_xlabel('Time [UTC and lead time in h]')
        ax.legend(loc = 0, fontsize = 8, frameon = False)
        ymax = np.ceil(ax.get_ylim()[1]*10)/10.
        ax.set_ylim(0, ymax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_position(('outward', 10))
        ax.spines['bottom'].set_position(('outward', 10))
        ax.set_xticks([0,6,12,18,24])
        ax.set_xlim(0, 24)
        ax.set_title('Deterministic forecasts ' + args.date_start)
        plt.tight_layout()
if args.ana == 'ens':
    axarr[0].set_xlabel('time [UTC/h]')
    
    
    
plotdir = plotdir + expid_str[:-1] + '/verif_fc_prec/'
if not os.path.exists(plotdir): os.makedirs(plotdir)

#print 'Save figure:', plotdir + plotstr + '.pdf'
fig1.savefig(plotdir + 'diprec_' + plotstr + '.pdf', format = 'pdf')
fig2.savefig(plotdir + 'fss_' + plotstr + '.pdf', format = 'pdf')

plt.close('all')

            
        
    
