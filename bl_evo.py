"""
Plot BL evolution
"""

# Imports
import argparse
import matplotlib
matplotlib.use("Agg")
from cosmo_utils.pywgrib import getfobj_ens, getfobj, fieldobj, getfield
from cosmo_utils.plot import ax_contourf
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, ddhhmmss_strtotime, yymmddhhmm, ddhhmmss
import matplotlib.pyplot as plt
import numpy as np
from datetime import timedelta
import os
from copy import deepcopy

# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--expid', metavar = 'expid', type=str, nargs = '+')
parser.add_argument('--date', metavar = 'date', type=str)
parser.add_argument('--time', metavar = 'time', type=str)
parser.add_argument('--box_coo', metavar = 'box_coo', type=int, nargs = '+')
args = parser.parse_args()



# Settings
exptag = ''
for exp in args.expid:
    exptag += exp + '_'
exptag = exptag[:-1]
plotdir = '/e/uwork/extsrasp/plots/' + exptag + '/bl_evo/'
if not os.path.exists(plotdir): os.makedirs(plotdir)
datadir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/'
precsuf = '_15'
tsuf = '_60'
gribpref = 'lfff'
t = ddhhmmss(timedelta(hours = int(args.time)))

cmPrec = ((1    , 1     , 1    ), 
          (0    , 0.627 , 1    ),
          (0.137, 0.235 , 0.98 ),
          (0.392, 0     , 0.627),
          (0.784, 0     , 0.627))
           #(0.1  , 0.1   , 0.784),
levelsPrec = [0, 0.3, 1, 3, 10, 30]
lev_max = 30

# Create box
ie1, ie2, je1, je2 = args.box_coo # first indices (ie*) are N-S
mask = np.zeros((461, 421))
mask[ie1:ie2, je1:je2] = 1

# Load HH
# hhfn = datadir + args.expid[0] + '/' + args.date + '/det/' + gribpref + ddhhmmss(timedelta(hours = 0)) + tsuf
# hhfobj = getfobj(hhfn, fieldn = 'HH')


##############################

ncols = len(args.expid)
fig1, axmat1 = plt.subplots(3, ncols, figsize = (4*ncols, 14))
if ncols == 1:
    axmat1 = [axmat1]

fig2, axmat2 = plt.subplots(2, 3, figsize = (12, 10))
cyc = ['b', 'g', 'r']
for i, exp in enumerate(args.expid):

    # Load precipitation 
    precfn = datadir + exp + '/' + args.date + '/det/' + gribpref + t  + precsuf
    precfobj = getfobj(precfn, fieldn = 'PREC_PERHOUR')

    maskfobj = deepcopy(precfobj)
    maskfobj.data = mask
    # Plot precipitation plus box
    plt.sca(axmat1[0,i])
    cf, tmp = ax_contourf(axmat1[0,i], precfobj,
                          Basemap_drawrivers = False, 
                          npars = 0, nmers = 0, 
                          colors = cmPrec, 
                          pllevels = levelsPrec,
                          sp_title = exp,
                          extend = 'max',
                          fieldobj_opt = maskfobj,
                          pllevels_opt = [0.5],
                          )
    cb = fig1.colorbar(cf, orientation = 'vertical', 
                            fraction = 0.05, pad = 0.0)

    # Load T or Theta
    tfn = datadir + exp + '/' + args.date + '/det/' + gribpref + t  + tsuf
    tfobj = getfobj(tfn, fieldn = 'T')
    tfield = tfobj.data
    pfield = getfield(tfn, fieldn = 'PS')
    wfobj = getfobj(tfn, fieldn = 'W')

    thetafield = tfield * (1000.e2 / pfield)**(2./7.)

    # Average for box and create vertical profile
    tprof = np.mean(np.mean(tfield[:, ie1:ie2, je1:je2], axis = 2), axis = 1)
    thetaprof = np.mean(np.mean(thetafield[:, ie1:ie2, je1:je2], axis = 2), 
                        axis = 1)
    if i == 0:
        refthetafield = thetafield
        reftheta = thetaprof
        refexp = exp
        refwfobj = wfobj
        
    difftheta = thetaprof - reftheta
    
    if i ==1:
        diffthetafield = thetafield - refthetafield
        diffthetafobj = deepcopy(tfobj)
        diffthetafobj.data = diffthetafield
        diffwfobj = deepcopy(wfobj)
        diffwfobj.data = wfobj.data - refwfobj.data
        tl500wfobj = wfobj
        
    if i ==2:

        pspwfobj = wfobj

    # Plot
    plt.sca(axmat2[0,0])
    axmat2[0,0].plot(tprof[lev_max:], range(tprof.shape[0])[lev_max:], 
                   label = exp, c = cyc[i])
    axmat2[0,1].plot(thetaprof[lev_max:], range(tprof.shape[0])[lev_max:], 
                   label = exp, c = cyc[i])
    axmat2[0,2].plot(difftheta[lev_max:], range(tprof.shape[0])[lev_max:], 
                   label = refexp + '-' + exp, c = cyc[i])

for il, lev in enumerate([49,45,41]):
    plt.sca(axmat1[1,il])
    cf, tmp = ax_contourf(axmat1[1,il], diffthetafobj,
                            Basemap_drawrivers = False, 
                            npars = 0, nmers = 0, lev = lev,
                            sp_title = 'REF_TL500 - REF theta ' + str(lev),
                            extend = 'both',
                            fieldobj_opt = maskfobj,
                            pllevels_opt = [0.5],
                            ji0 = (ie1, je1), ji1 = (ie2, je2),
                            pllevels = np.arange(-1,1.1,.1)
                            )
    cb = fig1.colorbar(cf, orientation = 'vertical', 
                            fraction = 0.05, pad = 0.0)
    if il > 0:
        plt.sca(axmat1[2,il])
        cf, tmp = ax_contourf(axmat1[2,il], diffwfobj,
                                Basemap_drawrivers = False, 
                                npars = 0, nmers = 0, lev = lev,
                                sp_title = 'REF_TL500 - REF w ' + str(lev),
                                extend = 'both',
                                fieldobj_opt = maskfobj,
                                pllevels_opt = [0.5],
                                ji0 = (ie1, je1), ji1 = (ie2, je2),
                                pllevels = np.arange(-1,1.1,.1)
                                )
        cb = fig1.colorbar(cf, orientation = 'vertical', 
                                fraction = 0.05, pad = 0.0)
        
plt.sca(axmat1[2,0])
cf, tmp = ax_contourf(axmat1[2,0], wfobj,
                        Basemap_drawrivers = False, 
                        npars = 0, nmers = 0, lev = 41,
                        sp_title = 'REF w ' + str(41),
                        extend = 'both',
                        fieldobj_opt = maskfobj,
                        pllevels_opt = [0.5],
                        ji0 = (ie1, je1), ji1 = (ie2, je2),
                        pllevels = np.arange(-1,1.1,.1)
                        )
cb = fig1.colorbar(cf, orientation = 'vertical', 
                    fraction = 0.05, pad = 0.0)


plt.tight_layout(rect=[0, 0.0, 1, 1])
plotstr1 = args.date + '_' + str(t) + '_prec'
print 'Saving as ', plotdir + plotstr1
fig1.savefig(plotdir + plotstr1, dpi = 150)


###
refw = refwfobj.data[:, ie1:ie2, je1:je2]
refwarray = []
tl500w = tl500wfobj.data[:, ie1:ie2, je1:je2]
tl500warray = []
psp500w = pspwfobj.data[:, ie1:ie2, je1:je2]
psp500warray = []
for z in range(refw.shape[0])[lev_max:]:
    refwarray.append(np.ravel(refw[z,:,:]))
    tl500warray.append(np.ravel(tl500w[z,:,:]))
    psp500warray.append(np.ravel(psp500w[z,:,:]))
    
plt.sca(axmat2[1,0])
axmat2[1,0].boxplot(refwarray, vert = False, labels = range(refw.shape[0])[lev_max:])
#axmat2[1,0].set_yticks(range(refw.shape[0])[lev_max:])
axmat2[1,0].set_ylabel('level')
axmat2[1,0].invert_yaxis()
axmat2[1,0].set_title('REF W')
axmat2[1,0].set_xlabel('w [m/s]')
axmat2[1,0].set_xlim((-2,2))
axmat2[1,1].boxplot(tl500warray, vert = False, labels = range(refw.shape[0])[lev_max:])
#axmat2[1,1].set_yticks(range(refw.shape[0])[lev_max:])
axmat2[1,1].set_ylabel('level')
axmat2[1,1].invert_yaxis()
axmat2[1,1].set_title('REF_TL500 W')
axmat2[1,1].set_xlabel('w [m/s]')
axmat2[1,1].set_xlim((-2,2))
axmat2[1,2].boxplot(psp500warray, vert = False, labels = range(refw.shape[0])[lev_max:])
#axmat2[1,1].set_yticks(range(refw.shape[0])[lev_max:])
axmat2[1,2].set_ylabel('level')
axmat2[1,2].invert_yaxis()
axmat2[1,2].set_title('PSP_TL500 W')
axmat2[1,2].set_xlabel('w [m/s]')
axmat2[1,2].set_xlim((-2,2))


axmat2[0,0].legend(loc = 1, fontsize = 6)
axmat2[0,0].set_xlabel('T [K]')
axmat2[0,0].set_ylabel('level')
axmat2[0,0].invert_yaxis()
axmat2[0,1].set_xlabel('Theta [K]')
axmat2[0,1].set_ylabel('level')
axmat2[0,1].invert_yaxis()
axmat2[0,2].legend(loc = 1, fontsize = 6)
axmat2[0,2].set_xlabel('Difference Theta [K]')
axmat2[0,2].set_ylabel('level')
axmat2[0,2].invert_yaxis()
plt.tight_layout(rect=[0, 0.0, 1, 1])
plotstr2 = args.date + '_' + str(t) + '_prof'
print 'Saving as ', plotdir + plotstr2
fig2.savefig(plotdir + plotstr2, dpi = 150)


plt.close('all')


