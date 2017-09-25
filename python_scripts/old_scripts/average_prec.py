
# Imports
import os
import numpy as np
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, make_timelist, \
    yyyymmddhhmmss
from datetime import timedelta

# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--expid', metavar = 'expid', type=str, nargs = '+')
parser.add_argument('--date_ini', metavar = 'date_ini', type=str,
                    default = '20160525000000')
parser.add_argument('--date_end', metavar = 'date_end', type=str,
                    default = '20160610000000')
parser.add_argument('--date', metavar = 'date', type=str, nargs = '+')
parser.add_argument('--time', metavar = 'time', type=str, nargs = '+')
args = parser.parse_args()

exptag = ''
for exp in args.expid:
    exptag += exp + '_'
exptag = exptag[:-1]

savedir = '/e/uwork/extsrasp/save/' + exptag + '/prec_time/'
plotdir = '/e/uwork/extsrasp/plots/' + exptag + '/prec_time/'
if not os.path.exists(plotdir): os.makedirs(plotdir)


# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_ini)
tend = yyyymmddhhmmss_strtotime(args.date_end)
tint = timedelta(days=1)
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)

matlist = []
for t in timelist:
    savefn = savedir + yyyymmddhhmmss(t) + '_' + str(args.time[0]) + '_' + str(args.time[1])
    tplot, savemat, labelslist = np.load(savefn + '.npy')
    matlist.append(savemat)
    
cdict = {'radar':'k',
        'REF':'b',
        'REF_TL500':'green',
        'PSP_TL500':'r',
        #'sig1':'green',
        #'time20':'orange',
        #'time10':'orange',
        #'const3':'purple',
        #'const1':'purple',
        #'nolowest':'gray',
        #'ref_tl500':'lightblue',
        #'std2_tl500':'darkblue',
        }
lwdict = {'radar':2,
        'REF':2,
        'REF_TL500':2,
        'PSP_TL500':2,
        #'sig1':1.5,
        #'time20':1.5,
        #'time10':1.5,
        #'const3':1.5,
        #'const1':1.5,
        #'nolowest':1.5,
        #'ref_tl500':1.5,
        #'std2_tl500':1.5,
        }
lsdict = {'radar':'-',
        'REF':'-',
        'REF_TL500':'-',
        'PSP_TL500':'-',
        #'sig1':'--',
        #'time20':'-',
        #'time10':'--',
        #'const3':'-',
        #'const1':'--',
        #'nolowest':'-',
        #'ref_tl500':'-',
        #'std2_tl500':'-',
        }

savemat = np.mean(matlist, axis = 0)
fig, ax = plt.subplots(1, 1, figsize = (6, 4))
for i in range(savemat.shape[1]):
    l = labelslist[i]
    ax.plot(tplot, savemat[:,i], c = cdict[l], label = l,
            linewidth = lwdict[l], linestyle = lsdict[l])
ax.legend(loc = 2, fontsize = 6)
ax.set_xlabel('time [UTC]')
ax.set_ylabel('di prec [mm/h]')
ax.set_title(args.date_ini + '_' + args.date_end)
plt.tight_layout()

plotstr = (args.date_ini + '_' + args.date_end + '_' + str(args.time[0]) + 
           '_' + str(args.time[1]))
print 'Saving as ', plotdir + plotstr
plt.savefig(plotdir + plotstr, dpi = 150)
plt.close('all')


    
