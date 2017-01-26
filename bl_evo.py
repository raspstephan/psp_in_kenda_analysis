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

box_coo_str = ''
for b in args.box_coo:
    box_coo_str += str(b) + '_'
box_coo_str = box_coo_str[:-1]

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
hhfn = datadir + args.expid[0] + '/' + args.date + '/det/' + gribpref + ddhhmmss(timedelta(hours = 0)) + 'c' + tsuf
hhfield = getfield(hhfn, fieldn = 'HH')
hhfield = (hhfield[1:] + hhfield[:-1]) / 2.
meanhh = np.mean(np.mean(hhfield[:, ie1:ie2, je1:je2], axis = 2), axis = 1)
del hhfield


##############################

ncols = len(args.expid)
fig1, axmat1 = plt.subplots(3, ncols, figsize = (4.5*ncols, 14))
if ncols == 1:
    axmat1 = [axmat1]

fig2, axmat2 = plt.subplots(3, 3, figsize = (12, 13))
cyc = ['b', 'g', 'r']
for i, exp in enumerate(args.expid):
    print 'Exp loop:', i
    
    # Load precipitation 
    print 'Load precipitation'
    precfn = datadir + exp + '/' + args.date + '/det/' + gribpref + t  + precsuf
    precfobj = getfobj(precfn, fieldn = 'PREC_PERHOUR')

    maskfobj = deepcopy(precfobj)
    maskfobj.data = mask
    
    # Plot precipitation plus box
    print 'Plot precipitation'
    plt.sca(axmat1[0,i])
    cf, tmp = ax_contourf(axmat1[0,i], precfobj,
                          Basemap_drawrivers = False, 
                          npars = 0, nmers = 0, 
                          colors = cmPrec, 
                          pllevels = levelsPrec,
                          sp_title = 'PREC [mm/h] ' + exp,
                          extend = 'max',
                          fieldobj_opt = maskfobj,
                          pllevels_opt = [0.5],
                          )
    cb = fig1.colorbar(cf, orientation = 'vertical', 
                            fraction = 0.05, pad = 0.0)

    # Load T or Theta
    print 'Load fields'
    tfn = datadir + exp + '/' + args.date + '/det/' + gribpref + t  + tsuf
    tfobj = getfobj(tfn, fieldn = 'T')
    tfield = tfobj.data
    qvfield = getfield(tfn, fieldn = 'QV')
    pfield = getfield(tfn, fieldn = 'PS')
    wfobj = getfobj(tfn, fieldn = 'W')
    #tpertfield = getfield(tfn, fieldn = 'var230')
    #qvpertfield = getfield(tfn, fieldn = 'var232')
    
    thetafield = tfield * (1000.e2 / pfield)**(2./7.)
    
    thetagrad = thetafield[49, ie1:ie2, je1:je2] - thetafield[40, ie1:ie2, je1:je2]
    #Rd = 287.
    #Rv = 461.
    #satfwa = 1.0007
    #satfwb = 3.46e-8
    #satewa = 611.21
    #satewb = 17.502
    #satewc = 32.18
    #rddrv = Rd / Rv
    
    #f = satfwa + satfwb * pfield
    #esl = f * satewa * np.exp( satewb*(tfield-273.15)/(tfield-satewc) )
    #qvsatfield =  rddrv * esl / (pfield-(1.0-rddrv)*esl)
    #rhfield = qvfield / qvsatfield * 100 # %

    # Average for box and create vertical profile
    print 'Average fields'
    tprof = np.mean(np.mean(tfield[:, ie1:ie2, je1:je2], axis = 2), axis = 1)
    #pprof = np.mean(np.mean(pfield[:, ie1:ie2, je1:je2], axis = 2), axis = 1)
    #rhprof = np.mean(np.mean(rhfield[:, ie1:ie2, je1:je2], axis = 2), axis = 1)
    qvprof = np.mean(np.mean(qvfield[:, ie1:ie2, je1:je2], axis = 2), axis = 1)
    thetaprof = np.mean(np.mean(thetafield[:, ie1:ie2, je1:je2], axis = 2), 
                        axis = 1)
    if i == 0:
        refthetafield = thetafield
        reftheta = thetaprof
        refexp = exp
        refwfobj = wfobj
        refthetagrad = thetagrad
        
    difftheta = thetaprof - reftheta
    
    if i ==1:
        diffthetafield = thetafield - refthetafield
        diffthetafobj = deepcopy(tfobj)
        diffthetafobj.data = diffthetafield
        diffwfobj = deepcopy(wfobj)
        diffwfobj.data = wfobj.data - refwfobj.data
        tl500wfobj = wfobj
        tl500thetagrad = thetagrad
        
    if i ==2:
        tstdfield = np.sqrt(getfield(tfn, fieldn = 'var245'))
        qvstdfield = np.sqrt(getfield(tfn, fieldn = 'var247'))
        tkefield = getfield(tfn, fieldn = 'TKE')
        tstdprof = np.mean(np.mean(tstdfield[:, ie1:ie2, je1:je2], axis = 2), 
                           axis = 1)
        qvstdprof = np.mean(np.nanmean(qvstdfield[:, ie1:ie2, je1:je2], axis = 2), 
                           axis = 1)
        tkeprof = np.mean(np.mean(tkefield[:, ie1:ie2, je1:je2], axis = 2), 
                           axis = 1)
        plt.sca(axmat2[0,0])
        axmat2[2,0].plot(tstdprof[lev_max:], range(tstdprof.shape[0])[lev_max:])
        axmat2[2,1].plot(qvstdprof[lev_max:]*1000., range(tstdprof.shape[0])[lev_max:])
        axmat2[2,2].plot(tkeprof[lev_max:], range(tstdprof.shape[0])[lev_max:])
        
        pspwfobj = wfobj
        pspthetagrad = thetagrad

    # Plot
    print 'Plot profiles'
    plt.sca(axmat2[0,0])
    axmat2[0,0].plot(qvprof[lev_max:]*1000., range(tprof.shape[0])[lev_max:], 
                   label = exp, c = cyc[i])
    axmat2[0,1].plot(thetaprof[lev_max:], range(tprof.shape[0])[lev_max:], 
                   label = exp, c = cyc[i])
    axmat2[0,2].plot(difftheta[lev_max:], range(tprof.shape[0])[lev_max:], 
                   label = exp + '-' + refexp, c = cyc[i])


for il, lev in enumerate([49,45,41]):
    plt.sca(axmat1[1,il])
    cf, tmp = ax_contourf(axmat1[1,il], diffthetafobj,
                            Basemap_drawrivers = False, 
                            npars = 0, nmers = 0, lev = lev,
                            sp_title = 'THETA [K] (REF_TL500 - REF) lvl ' + str(lev+1),
                            extend = 'both',
                            fieldobj_opt = maskfobj,
                            pllevels_opt = [0.5],
                            ji0 = (ie1, je1), ji1 = (ie2, je2),
                            pllevels = np.arange(-1,1.1,.1)
                            )
    cb = fig1.colorbar(cf, orientation = 'vertical', 
                       fraction=0.043, pad = 0.0)
    #if il > 0:
        #plt.sca(axmat1[2,il])
        #cf, tmp = ax_contourf(axmat1[2,il], diffwfobj,
                                #Basemap_drawrivers = False, 
                                #npars = 0, nmers = 0, lev = lev,
                                #sp_title = 'REF_TL500 - REF w ' + str(lev),
                                #extend = 'both',
                                #fieldobj_opt = maskfobj,
                                #pllevels_opt = [0.5],
                                #ji0 = (ie1, je1), ji1 = (ie2, je2),
                                #pllevels = np.arange(-1,1.1,.1)
                                #)
        #cb = fig1.colorbar(cf, orientation = 'vertical', 
                           #fraction=0.046, pad = 0.0)
        
plt.sca(axmat1[2,0])
#cf, tmp = ax_contourf(axmat1[2,0], wfobj,
                        #Basemap_drawrivers = False, 
                        #npars = 0, nmers = 0, lev = 41,
                        #sp_title = 'REF w ' + str(41),
                        #extend = 'both',
                        #fieldobj_opt = maskfobj,
                        #pllevels_opt = [0.5],
                        #ji0 = (ie1, je1), ji1 = (ie2, je2),
                        #pllevels = np.arange(-1,1.1,.1)
                        #)
#cb = fig1.colorbar(cf, orientation = 'vertical', 
                    #fraction = 0.05, pad = 0.0)
axmat1[2,0].hist(np.ravel(refthetagrad), range = (-1, 2), bins = 20)
axmat1[2,1].hist(np.ravel(tl500thetagrad), range = (-1, 2), bins = 20)
axmat1[2,2].hist(np.ravel(pspthetagrad), range = (-1, 2), bins = 20)
axmat1[2,0].set_xlabel('THETA grad [K]')
axmat1[2,0].set_ylabel('Number')
axmat1[2,0].set_title('THETA lvl 50 - lvl 41 (REF)')
axmat1[2,1].set_xlabel('THETA grad [K]')
axmat1[2,1].set_ylabel('Number')
axmat1[2,1].set_title('THETA lvl 50 - lvl 41 (REF_TL500)')
axmat1[2,2].set_xlabel('THETA grad [K]')
axmat1[2,2].set_ylabel('Number')
axmat1[2,2].set_title('THETA lvl 50 - lvl 41 (PSP_TL500)')


fig1.suptitle(args.date + ' + ' +  str(t) + 'h', fontsize = 14)
plt.tight_layout(rect=[0, 0.0, 1, 0.98])
plotstr1 = args.date + '_' + str(t) + '_' + box_coo_str + '_prec'
print 'Saving as ', plotdir + plotstr1
fig1.savefig(plotdir + plotstr1, dpi = 150)


###
refw = refwfobj.data[:, ie1:ie2, je1:je2]
refwarray = []
tl500w = tl500wfobj.data[:, ie1:ie2, je1:je2]
tl500warray = []
psp500w = pspwfobj.data[:, ie1:ie2, je1:je2]
psp500warray = []
psp500w = pspwfobj.data[:, ie1:ie2, je1:je2]
psp500warray = []
for z in range(refw.shape[0])[lev_max:]:
    refwarray.append(np.ravel(refw[z,:,:]))
    tl500warray.append(np.ravel(tl500w[z,:,:]))
    psp500warray.append(np.ravel(psp500w[z,:,:]))
    
plt.sca(axmat2[1,0])

# Plot w box plots
warr = [refwarray, tl500warray, psp500warray]
for ip in range(3):
    axmat2[1,ip].boxplot(warr[ip], vert = False, 
                         labels = range(refw.shape[0])[lev_max:])
    axmat2[1,ip].set_ylabel('level')
    axmat2[1,ip].invert_yaxis()
    axmat2[1,ip].set_title('W distribution ' + args.expid[ip])
    axmat2[1,ip].set_xlabel('w [m/s]')
    axmat2[1,ip].set_xlim((-2,2))
    
    # Settings for top row
    axmat2[0,ip].set_ylabel('level')
    axmat2[0,ip].invert_yaxis()
    
    # Settings for bottom row
    axmat2[2,ip].set_ylabel('level')
    axmat2[2,ip].invert_yaxis()
    

# individual settings top row
axmat2[0,0].legend(loc = 1, fontsize = 6)
axmat2[0,2].legend(loc = 1, fontsize = 6)
axmat2[0,0].set_xlabel('QV [g/kg]')
axmat2[0,1].set_xlabel('Theta [K]')
axmat2[0,2].set_xlabel('Difference Theta [K]')
axmat2[0,0].set_title('QV mean profiles')
axmat2[0,1].set_title('THETA mean profiles')
axmat2[0,2].set_title('THETA difference from REF')

# individual settings top row
axmat2[2,0].set_xlabel('std(T) [K]')
axmat2[2,1].set_xlabel('std(QV) [g/kg]')
axmat2[2,2].set_xlabel('TKE [J/kg]')
axmat2[2,0].set_xlim((0,0.5))
axmat2[2,1].set_xlim((0,0.5))
axmat2[2,2].set_xlim((0,2))
axmat2[2,0].set_title('Mean BL variability of T (PSP)')
axmat2[2,1].set_title('Mean BL variability of QV (PSP)')
axmat2[2,2].set_title('Mean TKE (PSP)')


fig2.suptitle(args.date + ' + ' + str(t) + 'h', fontsize = 14)
plt.tight_layout(rect=[0, 0.0, 1, 0.95])
plotstr2 = args.date + '_' + str(t) + '_' + box_coo_str + '_prof'
print 'Saving as ', plotdir + plotstr2
fig2.savefig(plotdir + plotstr2, dpi = 150)


plt.close('all')


