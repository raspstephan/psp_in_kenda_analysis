"""
Plotting precipitation stamps for my KENDA runs
Stephan Rasp

"""


# Imports
import argparse
from cosmo_utils.pywgrib import getfobj_ens
from cosmo_utils.plot import ax_contourf
import matplotlib.pyplot as plt
import numpy as np
import os

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
cmPrec = ((1    , 1     , 1    ), 
          (0    , 0.627 , 1    ),
          (0.137, 0.235 , 0.98 ),
          (0.392, 0     , 0.627),
          (0.784, 0     , 0.627))
           #(0.1  , 0.1   , 0.784),
levelsPrec = [0, 1, 3, 10, 30, 100.]

# TODO Only works for one date and expid, introduce loops later
# Load the data
topdir = datadir + args.expid[0] + '/' + args.date[0] + '/'
gribfn = gribpref + args.time[0] + precsuf
ensfobjlist = getfobj_ens(topdir, 'sub', mems = nens, gribfn = gribfn, 
                          dir_prefix = 'ens', fieldn = 'TOT_PREC_S', para = 2)

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
        
        cf, tmp = ax_contourf(ax, fobj, Basemap_drawrivers = False, npars = 0, 
                              nmers = 0, colors = cmPrec, pllevels = levelsPrec,
                              sp_title = 'mem ' + str(i+1))
        
    titlestr = args.expid[0] + ' ' + args.date[0] + ' + ' + args.time[0]
    fig.suptitle(titlestr, fontsize = 18)
    
    plt.tight_layout(rect=[0, 0.0, 1, 0.98])
    plotstr = str(args.time[0])
    plt.savefig(plotdirsub + plotstr, dpi = 300)
    plt.close('all')
    
    

