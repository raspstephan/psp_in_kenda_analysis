"""
Script fo analyze fof-Files and create verification plots 
Stephan Rasp
"""
# Imports
import argparse
import os
import numpy as np
from datetime import timedelta
from cosmo_utils.da import fdbkfile
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, make_timelist, \
    yyyymmddhhmmss
from scipy.stats import binned_statistic
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--expid', metavar = 'expid', type=str)
parser.add_argument('--ver_start_min', metavar = 'ver_start_min', type=int)
parser.add_argument('--ver_end_min', metavar = 'ver_end_min', type=int)
parser.add_argument('--date_ini', metavar = 'date_ini', type=str,
                    default = '20160525000000')
parser.add_argument('--date_end', metavar = 'date_end', type=str,
                    default = '20160610000000')
parser.add_argument('--var', metavar = 'var', type=str, default = 'T')
args = parser.parse_args()


# Config for experiment
FCINT='86400'  # forecast start interval in seconds (24h)
DATA_DIR='/e/uwork/extsrasp/cosmo_letkf/data_forecast/' + args.expid

# Config for verification
plotdir = '/e/uwork/extsrasp/plots/' + args.expid + '/verif/'
if not os.path.exists(plotdir): os.makedirs(plotdir)
savedir = '/e/uwork/extsrasp/save/' + args.expid + '/verif/'
if not os.path.exists(savedir): os.makedirs(savedir)

# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_ini)
tend = yyyymmddhhmmss_strtotime(args.date_end)
tint = timedelta(seconds=int(FCINT))
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)


obs_time = []   # in mins relative to obs_date
obs = []
verif = []
bias = []
obs_date = []   # datetime object
obs_lev = []    # in Pa
for t in timelist:
    print t
    # Load fof file
    foffn = (DATA_DIR + '/' + yyyymmddhhmmss(t) + '/det/fof_' +  
             yyyymmddhhmmss(t) + '.nc')
    fof = fdbkfile(foffn)
    var_tab = fof.table_varno
    varno = var_tab[args.var]
    
    # Get temp data, T
    ov = fof.obs_veri(varno=varno, obstype=5)
    hdr_inds = ov['hdr_inds']
    hdr_unique = np.unique(hdr_inds)  # Get unique obs labels 
    
    # Loop over obs 
    for obs_id in hdr_unique:
        obs_date.append(t)
        obs_time.append(ov['time'][hdr_inds == obs_id][0])
        obs.append(ov['obs'][hdr_inds == obs_id])
        verif.append(ov['veri_data'][hdr_inds == obs_id])
        bias.append(ov['veri_data'][hdr_inds == obs_id] -
                    ov['obs'][hdr_inds == obs_id])
        obs_lev.append(ov['level'][hdr_inds == obs_id])
    
    
print 'Total number of observations:', len(obs)
fig, ax = plt.subplots(1,1, figsize = (6, 4))
hist_edges = np.arange(0, 60*25, 60)
tmp = ax.hist(obs_time, bins = hist_edges, color = 'gray')
xmax = np.max(tmp[0]) * 1.2
ax.plot([args.ver_start_min, args.ver_start_min], [0, xmax],
        color = 'red')
ax.plot([args.ver_end_min, args.ver_end_min], [0, xmax],
        color = 'red')
ax.set_xlabel('Time UTC in minutes')
ax.set_ylabel('Number of TEMP observations')
plotstr = ('TEMP_timing_' + args.date_ini + '_' + args.date_end)
ax.set_title(args.expid + '_' + plotstr)
ax.set_ylim(0, xmax)
plt.tight_layout()
fig.savefig(plotdir + plotstr)
plt.close('all')

# Height bin the data
bias = np.array(bias)
obs_lev = np.array(obs_lev)
obs_time = np.array(obs_time)
mask = (obs_time >= args.ver_start_min) & (obs_time <= args.ver_end_min)
print 'Number of verif-observations:', np.sum(mask)

# Flatten masked lists
flatbias = [item for sublist in bias[mask] for item in sublist]
flatlev = [item for sublist in obs_lev[mask] for item in sublist]
bin_edges = np.arange(0, 1000, 50) * 100. # Pa

mean_bias = binned_statistic(flatlev, flatbias, bins = bin_edges)[0]
rmse = np.sqrt(binned_statistic(flatlev, np.array(flatbias)**2, 
                                bins = bin_edges)[0])

if args.var == 'RH':   # convert to percent
    mean_bias *= 100.
    rmse *= 100.

# Plot 
unitdict = {'T': 'K', 'RH': '%'}
rmselimdict = {'T': (0,3), 'RH': (0, 40)}
biaslimdict = {'T': (-3,3), 'RH': (-15,15)}
fig, axarr = plt.subplots(1,2, figsize = (10,5))
meanlev = (bin_edges[1:] + bin_edges[:-1])/2./100.  # hPa
axarr[0].plot(rmse, meanlev, c = 'k', linewidth = 2)
axarr[0].set_xlim(rmselimdict[args.var])
axarr[0].set_ylim(0,1000)
axarr[0].set_xlabel(args.var + ' RMSE [' + unitdict[args.var] +  ']')
axarr[0].set_ylabel('Pressure [hPa]')
axarr[0].set_title('RMSE ' + args.var)
axarr[0].invert_yaxis()
axarr[1].plot([0, 0],[0,1000],c = 'gray')
axarr[1].plot(mean_bias, meanlev, c = 'k', linewidth = 2)
axarr[1].set_xlim(biaslimdict[args.var])
axarr[1].set_ylim(0,1000)
axarr[1].set_xlabel(args.var + ' BIAS [' + unitdict[args.var] +  ']')
axarr[1].set_ylabel('Pressure [hPa]')
axarr[1].set_title('Bias ' + args.var)
axarr[1].invert_yaxis()
plt.tight_layout(rect=[0, 0.0, 1, 0.95])

plotstr = ('TEMP_' + args.var + '_' + args.date_ini + '_' + args.date_end + '_'
           + str(args.ver_start_min) + '_' + str(args.ver_end_min))
fig.suptitle(args.expid + '  ' + plotstr)
print 'Plotting:', plotdir + plotstr
fig.savefig(plotdir + plotstr)
plt.close('all')

# Save array
savefn = savedir + plotstr
print 'Saving array as', savefn 
np.save(savefn, (rmse, mean_bias, meanlev))




        
