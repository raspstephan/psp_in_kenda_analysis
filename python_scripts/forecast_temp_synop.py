"""Script to verify forecasts against TEMP and SYNOP observations.

Can do the following plots:
- deterministic rmse and bias
- ensemble rmse and spread
for
- TEMP/AIREP for a apecific time window
- SYNOP as a function of time
"""

from argparse import ArgumentParser
import sys
import numpy as np
import os
import config
import helpers as h
sys.path.append('/home/s/S.Rasp/repositories/kendapy')
from ekf import Ekf
from scipy.stats import binned_statistic


def compute_synop(inargs, exp_id, date):
    """Compute synop verification for one exp_id and date.
    
    Args:
        inargs: Command line arguments
        exp_id: Experiment ID string
        date: datetime object

    Returns:
        results: det --> [rmse, bias [time]]
                 ens --> [rmse, std [time]]
    """

    date_str = h.dt_to_yyyymmddhhmmss(date)

    if inargs.fc_type == 'det':
        fof_fn = (config.datadir + exp_id + date_str + '/det/fof_' + date_str +
                  '.nc')
        # Load ekf object and get data
        fof = Ekf(fof_fn)
        fkw = {'obstype': 'SYNOP', 'varname': inargs.var, 'state': inargs.state}
        obs = Ekf.obs(**fkw)
        fc = Ekf.fg(**fkw)
        time = Ekf.obs(param='time', **fkw)

        # Compute statistics
        diff = obs - fc

        # Bin statistic according to time
        bin_edges = np.arange(0, (args.hint + 1) * 60, 60)  # Hourly bins
        bias = binned_statistic(time, diff, bins=bin_edges)[0]
        rmse = np.sqrt(binned_statistic(time, np.array(diff) ** 2,
                                        bins=bin_edges)[0])
        results = [rmse, bias]
    elif inargs.fc_type == 'ens':
        # Loop over ensemble members
        ens_list = []
        for ie in range(20):   # ATTENTION: Hard coded ensemble size
            fof_fn = (config.datadir + exp_id + date_str + '/' +
                      str(ie+1).zfill(3) + '/fof_' + date_str + '.nc')
            # Load ekf object and get data
            fof = Ekf(fof_fn)
            fkw = {'obstype': 'SYNOP', 'varname': inargs.var,
                   'state': inargs.state}
            obs = Ekf.obs(**fkw)
            fc = Ekf.fg(**fkw)
            time = Ekf.obs(param='time', **fkw)
            ens_list.append(fc)

        # compute statistics
        ens_mean = np.mean(ens_list, axis=0)
        ens_spread = np.std(ens_list, axis=0, ddof=1)

        # Bin statistic according to time
        bin_edges = np.arange(0, (args.hint + 1) * 60, 60)  # Hourly bins
        spread = binned_statistic(time, ens_spread, bins=bin_edges)[0]
        rmse = np.sqrt(binned_statistic(time, np.array(obs - ens_mean) ** 2,
                                        bins=bin_edges)[0])
        results = [rmse, spread]
    else:
        raise TypeError, 'Wrong fc_type.'

    return results


def get_synop_for_one_day(inargs, exp_id, date):
    """Computes Synop verification for one day
    
    Args:
        inargs: Command line arguments
        exp_id: Experiment ID string
        date: datetime object

    Returns:
        results: Numpy array with dimensions [time, metric_dim]

    """

    # Create savestr
    save_fn = (config.savedir_base + exp_id + '/fof_' + inargs.obs + '_' +
               inargs.var + '_' + inargs.fc_type + '_' +
               h.dt_to_yyyymmddhhmmss(date) + '.npy')
    # Check if save_fn exists or recompute
    print 'Check if pre-computed file exists: %s' % save_fn
    if not os.path.exists(save_fn) or inargs.recompute:
        print 'Nope. Compute!'
        r = compute_synop(inargs, exp_id, date)
        np.save(save_fn, r)
    else:
        print 'Yep. Load!'
        r = np.load(save_fn)
    return r


def get_verification_for_all_dates(inargs, exp_id):
    """
    
    Args:
        inargs: Command line arguments
        exp_id: Experiment ID string 

    Returns:
        results: List
                 [rmse, bias] for det
                 [rmse, spread] for ens
                 of Numpy arrays with dimensions 
                 [date, time] for SYNOP
                 [date, height] for TEMP/AIREP
    """
    date_list = []
    # Loop over dates and collect data
    for date in h.make_timelist(inargs.date_start, inargs.date_stop,
                                inargs.hours_inc):
        if inargs.obs == 'SYNOP':
            date_list.append(get_synop_for_one_day(inargs, exp_id, date))
        elif inargs.obs in ['TEMP', 'AIREP']:
            date_list.append(get_temp_airep_for_one_day(inargs, exp_id, date))
        else:
            raise TypeError, 'Wrong obs type.'

    # TODO: Some dimension transform
    return something


# Main program
def main(inargs):
    """Main program

    Args:
        inargs: Command line arguments
    """
    # Compute verification
    results_list = []
    for e in inargs.exp_id:
        results_list.append(get_verification_for_all_dates(inargs, e))



if __name__ == '__main__':
    # Get command line arguments
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--exp_id', metavar='exp_id', type=str, nargs='+',
                        help='Experiment ID')
    parser.add_argument('--date_start',
                        type=str,
                        default='20160526000000',
                        help='Start date for date loop (yyyymmddhhmmss)')
    parser.add_argument('--date_stop',
                        type=str,
                        default='20160609000000',
                        help='End date for date loop (yyyymmddhhmmss)')
    parser.add_argument('--obs',
                        type=str,
                        help='SYNOP, TEMP or AIREP')
    parser.add_argument('--var',
                        type=str,
                        help='Variable name.')
    parser.add_argument('--state',
                        type=str,
                        default='all',
                        help='Observation state.')
    parser.add_argument('--fc_type',
                        type=str,
                        default='det',
                        help='det or ens')
    parser.add_argument('--recompute',
                        dest='recompute',
                        action='store_true',
                        help='Recompute pre-processed files.')
    parser.set_defaults(recompute=False)
    parser.add_argument('--composite',
                        dest='composite',
                        action='store_true',
                        help='Composite over all days')
    parser.set_defaults(composite=False)

    args = parser.parse_args()

    main(args)
