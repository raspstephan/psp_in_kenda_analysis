"""
Plotting precipitation stamps for my KENDA runs
Stephan Rasp

"""


# Imports
import argparse
import matplotlib
matplotlib.use("Agg")
from cosmo_utils.pywgrib import getfobj_ens, getfobj, fieldobj
from cosmo_utils.pyncdf import getfobj_ncdf, getfobj_ncdf_timeseries
from cosmo_utils.plot import ax_contourf
from cosmo_utils.diag import mean_spread_fieldobjlist, get_totmask
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, ddhhmmss_strtotime, \
    yymmddhhmm, ddhhmmss, make_timelist, yyyymmddhhmmss
from cosmo_utils.scores.probab import neighborhood_avg
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import timedelta
from scipy.signal import convolve2d
from helpers import save_fig_and_log, load_det, load_radar

from config import *

# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
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
                    default=24,
                    help='Time increment between forecasts (h)')
parser.add_argument('--hint', metavar='hint', type=int, default=24,
                    help='Maximum forecast lead time')
args = parser.parse_args()



#plotdir = '/e/uwork/extsrasp/plots/'
#datadir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/'# + args.dateid[0]   # TODO Date is hardcoded
precsuf = '_15'
gribpref = 'lfff'
nens = 20


#radardir = '/e/uwork/extsrasp/radolan/'
radarpref = 'raa01-rw_10000-'
radarsufx = '-dwd---bin.nc'

cmPrec = ((1    , 1     , 1    ), 
          (0    , 0.627 , 1    ),
          (0.137, 0.235 , 0.98 ),
          (0.392, 0     , 0.627),
          (0.784, 0     , 0.627))
           #(0.1  , 0.1   , 0.784),
levelsPrec = [0, 0.3, 1, 3, 10, 30]
levelsPrec_smooth = [0, 0.03, 0.1, 0.3, 1, 3]
cmtauc = ("#DAFF47","#EDB400","#E16F56","#BC2D82","#7F008D","#001889")
cmcape = ("#FFFFFF","#FFFFB1","#FFFFAA","#FFFDA3","#FFF99E","#FFF399","#FFED94","#FFE68E","#FFDE89","#FFD584","#FFCB7F","#FFC17A","#FDB675","#F6AA6F","#EE9D6A","#E59064","#DC825E","#D17357","#C56451","#B8534A","#AA4243","#9C2E3C","#8E1334","#7F002D")
levelstauc = [0, 1, 3, 6, 10, 20, 100]
levelscape = np.linspace(0, 2500, 25)
levelsprob = np.linspace(0,1,21)
cmprob = [(1,1,1,1)] + [plt.cm.plasma_r(i) for i in np.linspace(0.1, 1, 19)]

# Define loading functions
def load_radar_ts(time):
    radardt = timedelta(minutes = 10)
    dateobj_start = (yyyymmddhhmmss_strtotime(args.date[0]) + 
                     ddhhmmss_strtotime(time[0]) - radardt)
    dateobj_end = (yyyymmddhhmmss_strtotime(args.date[0]) + 
                   ddhhmmss_strtotime(time[1])- radardt)
    tinc = timedelta(hours = 1)
    radarts = getfobj_ncdf_timeseries(radardir + radarpref, dateobj_start, 
                                      dateobj_end, tinc, 
                                      refdate = yyyymmddhhmmss_strtotime(args.date[0]), 
                                      ncdffn_sufx = radarsufx, fieldn = 'pr',
                                      abs_datestr='yymmddhhmm',
                                      dwdradar = True)
    return radarts


# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_start)
tend = yyyymmddhhmmss_strtotime(args.date_stop)
tint = timedelta(hours=args.date_inc)
if tstart == tend:
    datelist = [tstart]
else:
    datelist = make_timelist(tstart, tend, tint)

hourlist = []
for h in range(1, args.hint + 1):
    hourlist

for idate, date in enumerate(datelist):
    print date
    datestr = yyyymmddhhmmss(date)
    for h in range(1, args.hint + 1):
        hstr = ddhhmmss(timedelta(hours=h))
        detfobjlist = []
        expid_str = ''
        for exp in args.expid:
            expid_str += exp + '_'
            DATA_DIR = datadir + exp
            detfobjlist.append(load_det(DATA_DIR, datestr, hstr))
        radarfobj = load_radar(datestr, hstr)

        plotfobjlist = [radarfobj] + detfobjlist
        titlelist = ['Radar']
        for e in args.expid:
            titlelist.append(e)

        nrows = int(np.ceil(len(plotfobjlist) / 3.))
        fig, axlist = plt.subplots(1, len(plotfobjlist),
                                   figsize=(4 * len(plotfobjlist), 5))
        # axlist = np.ravel(np.transpose(axmat))

        for i, fobj in enumerate(plotfobjlist):
            plt.sca(axlist[i])
            if i == 0:
                mask = fobj.data > 100
            else:
                mask = None
            cf, tmp = ax_contourf(axlist[i], fobj,
                                  Basemap_drawrivers=False,
                                  npars=0, nmers=0,
                                  colors=cmPrec,
                                  pllevels=levelsPrec,
                                  sp_title=titlelist[i],
                                  extend='max', mask=mask)
            if i < 0:
                cb = fig.colorbar(cf, orientation='horizontal',
                                  fraction=0.05, pad=0.0)
        plotstr = datestr + '_' + hstr
        fig.suptitle(plotstr, fontsize=18)

        plt.tight_layout(rect=[0, 0.0, 1, 0.95])
        plotdir = plotdir + expid_str[:-1] + '/prec_stamps/'
        if not os.path.exists(plotdir): os.makedirs(plotdir)
        save_fig_and_log(fig, 'det_stamps_' + plotstr, plotdir)

        
        
        
        
        
        
        
