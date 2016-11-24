"""
Plotting precipitation stamps for my KENDA runs
Stephan Rasp

"""


# Imports
import argparse
import matplotlib
matplotlib.use("Agg")
from cosmo_utils.pywgrib import getfobj_ens, getfobj, fieldobj
from cosmo_utils.pyncdf import getfobj_ncdf
from cosmo_utils.plot import ax_contourf
from cosmo_utils.diag import mean_spread_fieldobjlist
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, ddhhmmss_strtotime, yymmddhhmm, ddhhmmss
from cosmo_utils.scores.probab import neighborhood_avg
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import timedelta

# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--expid', metavar = 'expid', type=str, nargs = '+')
parser.add_argument('--date', metavar = 'date', type=str, nargs = '+')
parser.add_argument('--time', metavar = 'time', type=str, nargs = '+')
parser.add_argument('--plot', metavar = 'plot', type=str, nargs = '+')
args = parser.parse_args()

# General settings
plotdir = '/e/uwork/extsrasp/plots/'
datadir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/20160606_00_12_'   # TODO Date is hardcoded
precsuf = '_15'
gribpref = 'lfff'
nens = 20

radardir = '/e/uwork/extsrasp/radolan/'
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
cmtauc = ("#FFFFFF", "#DFE3BA","#F1D275","#E8A55E","#C06348","#7F002D")
cmcape = ("#FFFFFF","#FFFFB1","#FFFFAA","#FFFDA3","#FFF99E","#FFF399","#FFED94","#FFE68E","#FFDE89","#FFD584","#FFCB7F","#FFC17A","#FDB675","#F6AA6F","#EE9D6A","#E59064","#DC825E","#D17357","#C56451","#B8534A","#AA4243","#9C2E3C","#8E1334","#7F002D")
levelstauc = [0, 1, 3, 6, 10, 20, 100]
levelscape = np.linspace(0, 2500, 25)

# TODO Only works for one date and expid, introduce loops later

# Time loop
if len(args.time) == 1:
    timelist = [ddhhmmss(timedelta(hours=int(args.time[0])))]
else:
    tstart = int(args.time[0])
    tend = int(args.time[1])
    timelist = range(tstart, tend + 1)
    for i in range(len(timelist)):
        timelist[i] = ddhhmmss(timedelta(hours=timelist[i]))
print 'timelist', timelist
for t in timelist:
    # Load the data
    topdir = datadir + args.expid[0] + '/' + args.date[0] + '/'
    gribfn = gribpref + t + precsuf
    ensfobjlist = getfobj_ens(topdir, 'sub', mems = nens, gribfn = gribfn, 
                            dir_prefix = 'ens', fieldn = 'PREC_PERHOUR', 
                            para = 2)
    # Calculate ensemble mean and spread
    meanfobj, spreadfobj = mean_spread_fieldobjlist(ensfobjlist)


    if 'fc_obs_stamps' in args.plot:
        # Load deterministic data
        detfn = topdir + 'det/' + gribfn
        detfobj = getfobj(detfn, fieldn = 'PREC_PERHOUR')
        
        # Calculate tau_c  TODO Once cape and prec in same file, implement as derived field
        capefn = gribpref + t + '_60'
        capefobjlist = getfobj_ens(topdir, 'sub', mems = nens, gribfn = capefn, 
                                dir_prefix = 'ens', fieldn = 'CAPE_ML_S', 
                                para = 2)
        meancapefobj, tmp = mean_spread_fieldobjlist(capefobjlist)
        # This is now taken from derive_ncdf
        ## Coarse graining
        #n = np.ceil(60./2.8) // 2 * 2 + 1  # Closes odd neighborhood size
        #prec_field = neighborhood_avg(detfobj.data, n, boundary = True)
        #cape_field = neighborhood_avg(capefobj.data, n, boundary = True)

        # Gaussian filter
        sig = 60./2.8/2.   # Sigma for Gaussian filtering 60 km
        prec_field = gaussian_filter(meanfobj.data, sig)
        cape_field = gaussian_filter(meancapefobj.data, sig)

        # Calculate tau_c
        tau_c = 0.5*(49.58/3600.)*cape_field/prec_field   # Factor from F.Heinlein
        tau_c[prec_field < 0.02] = np.nan   # Apply low precipitation threshold
        
        # Write new fobj
        tau_c_fobj = fieldobj(\
            data = tau_c,
            fieldn = 'tau_c',
            fieldn_long = 'Convective adjustment timescale',
            unit = 'hours',
            levs_inp = detfobj.levs,
            rlats = detfobj.rlats[:,0],
            rlons = detfobj.rlons[0,:],
            polelat = detfobj.polelat,
            polelon = detfobj.polelon,
            gribfn = detfn,
            )

        
        # Load radar data
        dateobj = (yyyymmddhhmmss_strtotime(args.date[0]) + 
                ddhhmmss_strtotime(t))
        radardt = timedelta(minutes = 70)   # TODO Is this correct???
        radardateobj = dateobj - radardt
        radarfn = radardir + radarpref + yymmddhhmm(radardateobj) + radarsufx
        radarfobj = getfobj_ncdf(radarfn, fieldn = 'pr', dwdradar = True)
        
    # Plot what is to be plotted
    if 'ens_stamps' in args.plot:
        print 'Plotting ens_stamps'
        plotdirsub = (plotdir + '20160606_00_12_' + args.expid[0] + '/' + 
                    args.date[0] + '/ens_stamps/')
        if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
        
        fig, axmat = plt.subplots(4, 5, figsize = (20,20))
        axlist = np.ravel(axmat)
        
        for i, (ax, fobj) in enumerate(zip(axlist, ensfobjlist)):
            plt.sca(ax)
            
            cf, tmp = ax_contourf(ax, fobj, Basemap_drawrivers = False, 
                                  npars = 0, 
                                  nmers = 0, colors = cmPrec, 
                                  pllevels = levelsPrec,
                                  sp_title = 'mem ' + str(i+1))
        
        titlestr = args.expid[0] + ' ' + args.date[0] + ' + ' + t
        fig.suptitle(titlestr, fontsize = 18)
        
        plt.tight_layout(rect=[0, 0.0, 1, 0.98])
        plotstr = str(t)
        plt.savefig(plotdirsub + plotstr, dpi = 300)
        plt.close('all')
        
        
    if 'fc_obs_stamps' in args.plot:
        print 'Plotting fc_obs_stamps'
        plotdirsub = (plotdir + '20160606_00_12_' + args.expid[0] + '/' + 
                    args.date[0] + '/fc_obs_stamps/')
        if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
        
        
        
        fig, axmat = plt.subplots(2, 3, figsize = (12,10))
        axlist = np.ravel(axmat)
        
        tmpfobjlist = [detfobj, meancapefobj, tau_c_fobj, radarfobj, 
                    meanfobj, spreadfobj]
        titlelist = ['det prec [mm/h]', 'ens mean CAPE [J/kg]', 
                     'ens mean TAU_C [h]', 'radar prec [mm/h]', 
                     'ens mean prec [mm/h]',  'ens spread prec [mm/h]']
        colorslist = [cmPrec, cmcape, cmtauc, cmPrec, cmPrec, cmPrec]
        levelslist = [levelsPrec, levelscape, levelstauc, levelsPrec, 
                      levelsPrec,levelsPrec]
        
        for i, ax in enumerate(axlist):
            plt.sca(ax)
            
            cf, tmp = ax_contourf(ax, tmpfobjlist[i],
                                  Basemap_drawrivers = False, 
                                  npars = 0, nmers = 0, 
                                  colors = colorslist[i], 
                                  pllevels = levelslist[i],
                                  sp_title = titlelist[i],
                                  extend = 'max')
            if i < 3:
                cb = fig.colorbar(cf, orientation = 'horizontal', 
                                fraction = 0.05, pad = 0.0)
        titlestr = args.expid[0] + ' ' + args.date[0] + ' + ' + t
        fig.suptitle(titlestr, fontsize = 18)
        
        plt.tight_layout(rect=[0, 0.0, 1, 0.95])
        plotstr = str(t)
        print 'Saving as ', plotdirsub + plotstr
        plt.savefig(plotdirsub + plotstr, dpi = 300)
        plt.close('all')
