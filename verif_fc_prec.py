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
    plotdir = '/home/s/S.Rasp/dwd_work/plots/'
    savedir_base = '/home/cosmo/stephan.rasp/dwd_data/save/'
else: 
    raise Exception('Working directory not recognized:' + os.getcwd())

radarpref = 'raa01-rw_10000-'
radarsufx = '-dwd---bin.nc'
precsuf = '_15'
gribpref = 'lfff'
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

n = 21
kernel = np.ones((n,n))/float((n*n))

# Set up figure
fig, axarr = plt.subplots(1, 2, figsize = (10, 6))

expid_str = ''
for ie, expid in enumerate(args.expid):
    print 'expid = ', expid
    expid_str += expid + '_'
    DATA_DIR = datadir + expid
    
    if args.ana == 'det':
        radarmean = []
        detmean = []
        hourlist =[]
        detrmse = []
        fsslist = []
    
    for t in timelist:
        print t
        date = yyyymmddhhmmss(t)
        for h in range(1, args.hint+1):
            print 'Hour', h
            hourlist.append(h)
            # Load data 
            if args.ana == 'det':
                # Det run 
                detfobj = load_det(DATA_DIR, date, ddhhmmss(timedelta(hours=h)))
                radarfobj = load_radar(date, ddhhmmss(timedelta(hours=h)))
                
                radarmean.append(np.mean(radarfobj.data[~totmask]))
                detmean.append(np.mean(detfobj.data[~totmask]))
                
                # Upscaled RMSE
                nanradar = radarfobj.data
                nandet = detfobj.data
                nanradar[totmask] = np.nan
                nandet[totmask] = np.nan 
                convradar = convolve2d(nanradar, kernel, mode='same')
                convdet = convolve2d(nandet, kernel, mode='same')
                detrmse.append(np.sqrt(np.nanmean((convradar - convdet)**2 / 
                                       (0.5*(convradar + convdet))**2)))
                fss = FSS(1, 21, nanradar, nandet, python_core = True)
                fsslist.append(fss)
        # End hourlist 
    # End timelist
    
    # Bin the data 
    hour_bins = np.arange(0.5, 25.5, 1)
    radarmean = binned_statistic(hourlist, radarmean, 
                                    bins = hour_bins)[0]
    detmean = binned_statistic(hourlist, detmean, 
                                    bins = hour_bins)[0]
    detrmse = binned_statistic(hourlist, detrmse, 
                                    bins = hour_bins)[0]
    meanfss = binned_statistic(hourlist, fsslist, 
                                    bins = hour_bins)[0]
    
    # Plot
    axarr[0].plot(range(radarmean.shape[0]), radarmean, c = 'k', linewidth = 2)
    axarr[0].plot(range(radarmean.shape[0]), detmean, c = cdict[expid], 
                  linewidth = 2,
                  label = expid)
    axarr[1].plot(range(radarmean.shape[0]), detrmse, c = cdict[expid], 
                  linewidth = 2)
    axarr[1].plot(range(radarmean.shape[0]), meanfss, c = cdict[expid], 
                  linewidth = 2, linestyle = '--')

# End expid loop
# Finish the plots

axarr[0].set_xlabel('time [UTC/h]')
axarr[1].set_xlabel('time [UTC/h]')
axarr[0].set_label('[mm/h]')
axarr[1].set_label('FSS')
axarr[0].set_title('domain mean precipitation')
axarr[1].set_title('FSS n = ' + str(n*2.8) + 'km')
axarr[0].legend(loc = 0, fontsize = 6)
plt.tight_layout(rect=[0, 0.0, 1, 0.95])

plotdir = plotdir + expid_str[:-1] + '/verif_fc_prec/'
if not os.path.exists(plotdir): os.makedirs(plotdir)
fig.suptitle(expid_str[:-1] + '  ' + plotstr)
print 'Save figure:', plotdir + plotstr
fig.savefig(plotdir + plotstr, dpi = 300)
plt.close('all')

            
        
    
