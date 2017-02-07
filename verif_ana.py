"""
Verify ekf and fof files from analysis
"""
# Imports
import argparse
import numpy as np
from cosmo_utils.da import fdbkfile
from cosmo_utils.helpers import yyyymmddhhmmss, yyyymmddhhmmss_strtotime, \
                                make_timelist
from datetime import timedelta

# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--expid', metavar = 'expid', type=str)
parser.add_argument('--date_ana_start', metavar = 'date_ana_start', type=str)
parser.add_argument('--date_ana_stop', metavar = 'date_ana_stop', type=str)
parser.add_argument('--var', metavar = 'var', type=str, default = 'T2M')
parser.add_argument('--obs', metavar = 'obs', type=str, default = 'SYNOP')
args = parser.parse_args()


# Loop over time
tstart = yyyymmddhhmmss_strtotime(args.date_ana_start)
tend = yyyymmddhhmmss_strtotime(args.date_ana_stop)
tint = timedelta(hours = 1)
if tstart == tend:
    timelist = [tstart]
else:
    timelist = make_timelist(tstart, tend, tint)

mean_spread = []
mean_bias = []
rmse = []

for t in timelist:
    # Determine 3hrly storage time
    date_ana = yyyymmddhhmmss(t)
    print 'date_ana = ', date_ana
    date_fg = yyyymmddhhmmss(t - timedelta(hours = 1))
    date_store = yyyymmddhhmmss(t - timedelta(hours = t.hour % 3))

    # Directory where ekf's and fof's are stored
    FEED_DIR = ('/e/uwork/extsrasp/cosmo_letkf/feedback/' + args.expid +
                '/' + date_store + '/')

    # Load files
    try:
        ekf = fdbkfile(FEED_DIR + 'ekf' + args.obs + '_' + date_ana + '.nc')
        fof = fdbkfile(FEED_DIR + 'fof_' + date_fg + '.nc')
    except:
        print 'Either ekf or fof file does not exsist. Skip over rest of loop.'
        print(FEED_DIR + 'ekf' + args.obs + '_' + date_ana + '.nc')
        print(FEED_DIR + 'fof_' + date_fg + '.nc')
        continue  # Skip over rest of loop

    ov_ekf = ekf.obs_veri(varno=args.var, obstype=args.obs, 
                              veri = 'first guess ensemble spread')
    ov_fof = fof.obs_veri(varno=args.var, obstype=args.obs, veri = 0)

    if args.obs == 'TEMP':
        # NOTE: ekf and fof files have different hdr_inds, but data is in same order
        # NO something is not right here ...
        hdr_inds_ekf = ov_ekf['hdr_inds']
        hdr_inds_fof = ov_ekf['hdr_inds']
        hdr_unique = np.unique(hdr_inds)

    elif args.obs == 'SYNOP':
        ens_spread = ov_ekf['veri_data']
        det_bias = ov_fof['veri_data'] - ov_fof['obs']  # neg means fc colder than obs

        # Get mean values
        mean_spread.append(np.mean(ens_spread))
        mean_bias.append(np.mean(det_bias))
        rmse.append(np.sqrt(np.mean(det_bias**2)))

    else:
        raise Exception

print mean_spread, mean_bias, rmse


