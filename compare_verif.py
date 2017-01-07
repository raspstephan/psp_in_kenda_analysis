
# Imports
import os
import numpy as np
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--expid', metavar = 'expid', type=str, nargs = '+')
parser.add_argument('--ver_start_min', metavar = 'ver_start_min', type=int)
parser.add_argument('--ver_end_min', metavar = 'ver_end_min', type=int)
parser.add_argument('--date_ini', metavar = 'date_ini', type=str,
                    default = '20160525000000')
parser.add_argument('--date_end', metavar = 'date_end', type=str,
                    default = '20160610000000')
args = parser.parse_args()

plotstr = ('TEMP_T_' + args.date_ini + '_' + args.date_end + '_'
           + str(args.ver_start_min) + '_' + str(args.ver_end_min))

# Plot 
cyc = ['b', 'g', 'r']
fig, axarr = plt.subplots(1,2, figsize = (10,5))
expid_str = ''
for i, expid in enumerate(args.expid):
	expid_str += expid + '_'
	# Load data
	savedir = '/e/uwork/extsrasp/save/' + expid + '/verif/'
	savefn = savedir + plotstr
	print 'Loading', savefn 
	rmse, mean_bias, meanlev =  np.load(savefn + '.npy')

	axarr[0].plot(rmse, meanlev, c = cyc[i], linewidth = 2, label = expid)
	axarr[1].plot(mean_bias, meanlev, c = cyc[i], linewidth = 2, label = expid)
expid_str = expid_str[:-1]

axarr[0].set_xlim(0,3)
axarr[0].set_ylim(0,1100)
axarr[0].set_xlabel('T RMSE [K]')
axarr[0].set_ylabel('Pressure [hPa]')
axarr[0].set_title('RMSE Temperature')
axarr[0].invert_yaxis()
axarr[1].plot([0, 0],[0,1100],c = 'gray')
axarr[1].set_xlim(-3,3)
axarr[1].set_ylim(0,1100)
axarr[1].set_xlabel('T bias [K]')
axarr[1].set_ylabel('Pressure [hPa]')
axarr[1].set_title('Bias Temperature')
axarr[1].invert_yaxis()
axarr[1].legend(loc = 1)
plt.tight_layout(rect=[0, 0.0, 1, 0.95])

plotdir = '/e/uwork/extsrasp/plots/' + expid_str + '/verif/'
if not os.path.exists(plotdir): os.makedirs(plotdir)
fig.suptitle(expid_str + '  ' + plotstr)
fig.savefig(plotdir + plotstr)
plt.close('all')





