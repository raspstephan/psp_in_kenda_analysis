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
parser.add_argument('--dateid', metavar = 'dateid', type=str, nargs = '+')
parser.add_argument('--date', metavar = 'date', type=str, nargs = '+')
parser.add_argument('--time', metavar = 'time', type=str, nargs = '+')
parser.add_argument('--plot', metavar = 'plot', type=str, nargs = '+')
parser.add_argument('--plotens', metavar = 'plotens', type=str, default = 'False')
args = parser.parse_args()

# General settings
plotdir = '/e/uwork/extsrasp/plots/'
datadir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/' + args.dateid[0]   # TODO Date is hardcoded
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

# Define loading functions
def load_det(expid, t):
    topdir = datadir + expid + '/' + args.date[0] + '/'
    gribfn = gribpref + t + precsuf
    detfn = topdir + 'det/' + gribfn
    detfobj = getfobj(detfn, fieldn = 'PREC_PERHOUR')
    return detfobj

def load_radar(t):
    dateobj = (yyyymmddhhmmss_strtotime(args.date[0]) + 
            ddhhmmss_strtotime(t))
    radardt = timedelta(minutes = 70)   # TODO Is this correct???
    radardateobj = dateobj - radardt
    radarfn = radardir + radarpref + yymmddhhmm(radardateobj) + radarsufx
    radarfobj = getfobj_ncdf(radarfn, fieldn = 'pr', dwdradar = True)
    return radarfobj

def load_radar_ts(time):
    radardt = timedelta(minutes = 70)
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

def load_ens(expid, t):
    topdir = datadir + expid + '/' + args.date[0] + '/'
    gribfn = gribpref + t + precsuf
    ensfobjlist = getfobj_ens(topdir, 'sub', mems = nens, gribfn = gribfn, 
                              dir_prefix = 'ens', fieldn = 'PREC_PERHOUR', 
                              para = 2)
    return ensfobjlist

def load_det_cape(expid, t):
    topdir = datadir + expid + '/' + args.date[0] + '/'
    if expid in ['ref', 'std'] and args.dateid == '20160606_00_12_':
        gribfn = gribpref + t + '_60'
    else:
        gribfn = gribpref + t + '_15'
    detfn = topdir + 'det/' + gribfn
    detfobj = getfobj(detfn, fieldn = 'CAPE_ML_S')
    return detfobj







# Time loop
if len(args.time) == 1:
    timelist = [ddhhmmss(timedelta(hours=int(args.time[0])))]
else:
    tstart = int(args.time[0])
    tend = int(args.time[1])
    tplot = range(tstart, tend + 1)
    timelist = []
    for i in range(len(tplot)):
        timelist.append(ddhhmmss(timedelta(hours=tplot[i])))
print 'timelist', timelist

if 'prec_time' in args.plot:
    radarts = load_radar_ts([timelist[0], timelist[-1]])
    totmask = get_totmask(radarts)
    
savelist = []
savelist_ens = []
for it, t in enumerate(timelist):
    # Load the data
    #topdir = datadir + args.expid[0] + '/' + args.date[0] + '/'
    #gribfn = gribpref + t + precsuf
    #ensfobjlist = getfobj_ens(topdir, 'sub', mems = nens, gribfn = gribfn, 
                            #dir_prefix = 'ens', fieldn = 'PREC_PERHOUR', 
                            #para = 2)
    ## Calculate ensemble mean and spread
    #meanfobj, spreadfobj = mean_spread_fieldobjlist(ensfobjlist)


    #if 'fc_obs_stamps' in args.plot:
        
        ## Calculate tau_c  TODO Once cape and prec in same file, implement as derived field
        #capefn = gribpref + t + '_60'
        #capefobjlist = getfobj_ens(topdir, 'sub', mems = nens, gribfn = capefn, 
                                #dir_prefix = 'ens', fieldn = 'CAPE_ML_S', 
                                #para = 2)
        #meancapefobj, tmp = mean_spread_fieldobjlist(capefobjlist)
        ## This is now taken from derive_ncdf
        ### Coarse graining
        ##n = np.ceil(60./2.8) // 2 * 2 + 1  # Closes odd neighborhood size
        ##prec_field = neighborhood_avg(detfobj.data, n, boundary = True)
        ##cape_field = neighborhood_avg(capefobj.data, n, boundary = True)

        ## Gaussian filter
        #sig = 60./2.8/2.   # Sigma for Gaussian filtering 60 km
        #prec_field = gaussian_filter(meanfobj.data, sig)
        #cape_field = gaussian_filter(meancapefobj.data, sig)

        ## Calculate tau_c
        #tau_c = 0.5*(49.58/3600.)*cape_field/prec_field   # Factor from F.Heinlein
        #tau_c[prec_field < 0.02] = np.nan   # Apply low precipitation threshold
        
        ## Write new fobj
        #tau_c_fobj = fieldobj(\
            #data = tau_c,
            #fieldn = 'tau_c',
            #fieldn_long = 'Convective adjustment timescale',
            #unit = 'hours',
            #levs_inp = detfobj.levs,
            #rlats = detfobj.rlats[:,0],
            #rlons = detfobj.rlons[0,:],
            #polelat = detfobj.polelat,
            #polelon = detfobj.polelon,
            #gribfn = detfn,
            #)

        
        
        
    # Plot what is to be plotted
    if 'ens_stamps' in args.plot:
        print 'Plotting ens_stamps'
        plotdirsub = (plotdir + args.dateid[0] + args.expid[0] + '/' + 
                    args.date[0] + '/ens_stamps/')
        if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
        
        # Load data
        ensfobjlist = load_ens(args.expid[0], t)
        
        
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
        # TODO Not updated yet
        plotdirsub = (plotdir + args.dateid[0] + args.expid[0] + '/' + 
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
        
    if 'prec_comp' in args.plot:
        print 'Plotting prec_comp'
        exptag = ''
        detfobjlist = []
        for exp in args.expid:
            exptag += exp + '_'
            detfobjlist.append(load_det(exp, t))
        exptag = exptag [:-1]
        capefobj = load_det_cape(args.expid[0], t)    
        radarfobj = load_radar(t)
        
        plotdirsub = (plotdir + args.dateid[0] + exptag + '/' + 
                    args.date[0] + '/prec_comp/')
        if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
        
        fig, axmat = plt.subplots(2, 2, figsize = (8, 10))
        axlist = np.ravel(axmat)
        
        plotfobjlist = [radarfobj, capefobj] + detfobjlist
        titlelist = ['radar prec [mm/h]', 'ens mean CAPE [J/kg]',
                     args.expid[0] + ' det prec [mm/h]',
                     args.expid[1] + ' det prec [mm/h]']
        colorslist = [cmPrec, cmcape, cmPrec, cmPrec]
        levelslist = [levelsPrec, levelscape, levelsPrec, levelsPrec]
        
        for i, ax in enumerate(axlist):
            plt.sca(ax)
            cf, tmp = ax_contourf(ax, plotfobjlist[i],
                                  Basemap_drawrivers = False, 
                                  npars = 0, nmers = 0, 
                                  colors = colorslist[i], 
                                  pllevels = levelslist[i],
                                  sp_title = titlelist[i],
                                  extend = 'max')
            if i < 2:
                cb = fig.colorbar(cf, orientation = 'horizontal', 
                                fraction = 0.05, pad = 0.0)
        titlestr = args.date[0] + ' + ' + t
        fig.suptitle(titlestr, fontsize = 18)
        
        plt.tight_layout(rect=[0, 0.0, 1, 0.95])
        plotstr = str(t)
        print 'Saving as ', plotdirsub + plotstr
        plt.savefig(plotdirsub + plotstr, dpi = 300)
        plt.close('all')
        
    if 'prec_time' in args.plot:
        print 'Plotting prec_time'
        radarfobj = radarts[it]
        means = [np.mean(radarfobj.data[~totmask])]
        exptag = ''
        tmplist = []
        for exp in args.expid:
            exptag += exp + '_'
            detfobj = load_det(exp, t)
            means += [np.mean(detfobj.data[~totmask])]
            if args.plotens == 'True':
                ensfobjlist = load_ens(exp, t)
                tmplist2 = []
                for ensfobj in ensfobjlist:
                    tmplist2.append(np.mean(ensfobj.data[~totmask]))
                tmplist.append(tmplist2)
        if args.plotens == 'True':
            savelist_ens.append(tmplist)
        exptag = exptag [:-1]
        savelist.append(means)

if 'prec_time' in args.plot:
    plotdirsub = (plotdir + args.dateid[0] + exptag + '/' + 
                    args.date[0] + '/prec_time/')
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    savemat = np.array(savelist)
    savemat_ens = np.array(savelist_ens)
    clist = ['k', 'lime', 'orangered']
    cycg = [plt.cm.Greens(i) for i in np.linspace(0.1, 0.9, nens)]
    cycr = [plt.cm.Reds(i) for i in np.linspace(0.1, 0.9, nens)]
    cyclist = [cycg, cycr]
    labelslist = ['radar'] + args.expid
    fig, ax = plt.subplots(1, 1, figsize = (6, 4))
    for i in range(savemat.shape[1]):
        ax.plot(tplot, savemat[:,i], c = clist[i], label = labelslist[i],
                linewidth = 3)
    if args.plotens == 'True':
        for i in range(savemat_ens.shape[1]):
            for j in range(savemat_ens.shape[2]):
                ax.plot(tplot, savemat_ens[:,i,j], c = cyclist[i][j],
                        zorder = 0.1)
    
    ax.legend(loc = 2)
    ax.set_xlabel('time [UTC]')
    ax.set_ylabel('di prec [mm/h]')
    ax.set_title(args.date[0])
    plt.tight_layout()
    
    plotstr = str(args.time[0]) + '_' + str(args.time[1])
    print 'Saving as ', plotdirsub + plotstr
    plt.savefig(plotdirsub + plotstr, dpi = 300)
    plt.close('all')
        
        
        
        
        
        
        
