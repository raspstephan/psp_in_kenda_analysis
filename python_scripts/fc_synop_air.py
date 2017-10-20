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
import pdb


def compute_det_synop(diff, time):
    """Compute deterministic SYNOP statistics

    Args:
        diff: Numpy array with difference obs - fc for each data point
        time: Corresponding time array

    Returns:
        [rmse, bias]
    """
    # Bin statistic according to time
    print 'SYNOP total obs: %i' % diff.shape[0]
    bin_edges = np.arange(60, (24 + 1) * 60, 60)  # Hourly bins
    bias = binned_statistic(time, diff, bins=bin_edges)[0]
    rmse = np.sqrt(binned_statistic(time, np.array(diff) ** 2,
                                    bins=bin_edges)[0])
    return [rmse, bias]


def compute_det_air(diff, time, plevel, verif_time, obs_type):
    """Compute deterministic upper air statistics

    Args:
        diff: Numpy array with difference fc - obs for each data point
        time: Corresponding time array
        plevel: Pressure for obs
        verif_time: list with start and stop verif time iin hours
        obs_type: TEMP or AIREP

    Returns:
        [rmse, bias]
    """
    # Cut out requested time
    mask = (time >= verif_time[0] * 60.) & (time <= verif_time[1] * 60.)
    plevel = plevel[mask]
    diff = diff[mask]
    print 'Obs: %i/%i' % (np.sum(mask), mask.shape[0])

    # Bin statistic according to pressure level
    bin_edges = config.temp_bin_edges if obs_type == 'TEMP' \
        else config.airep_bin_edges
    bias = binned_statistic(plevel, diff, bins=bin_edges)[0]
    rmse = np.sqrt(binned_statistic(plevel, np.array(diff) ** 2,
                                    bins=bin_edges)[0])
    return [rmse, bias]


def compute_ens_air(ens_mean, ens_spread, obs, time, plevel, verif_time,
                    obs_type):
    """Compute ensemble upper air statistics

    Args:
        diff: Numpy array with difference obs - fc for each data point
        time: Corresponding time array
        plevel: Pressure for obs
        verif_time: list with start and stop verif time iin hours
        obs_type: TEMP or AIREP

    Returns:
        [rmse, spread]
    """
    # Cut out requested time
    mask = (time >= verif_time[0] * 60.) & (time <= verif_time[1] * 60.)
    plevel = plevel[mask]
    ens_mean = ens_mean[mask]
    ens_spread = ens_spread[mask]
    obs = obs[mask]
    print 'Obs: %i/%i' % (np.sum(mask), mask.shape[0])

    # compute statistics
    diff = ens_mean - obs

    # Bin statistic according to pressure level
    bin_edges = config.temp_bin_edges if obs_type == 'TEMP' \
        else config.airep_bin_edges
    spread = binned_statistic(plevel, ens_spread, bins=bin_edges)[0]
    rmse = np.sqrt(binned_statistic(plevel, np.array(diff) ** 2,
                                    bins=bin_edges)[0])
    return [rmse, spread]


def compute_verif(inargs, exp_id, date):
    """Compute verification for one exp_id and date.
    
    Args:
        inargs: Command line arguments
        exp_id: Experiment ID string
        date: datetime object

    Returns:
        results: det --> [rmse, bias]
                 ens --> [rmse, std]
                 with [time] for SYNOP
                      [height] for TEMP/AIREP
    """

    date_str = h.dt_to_yyyymmddhhmmss(date)
    fkw = {'obstype': inargs.obs, 'varname': inargs.var, 'state': inargs.state}

    if inargs.fc_type == 'det':
        fof_fn = (config.datadir + exp_id + '/' + date_str + '/det/fof_' +
                  date_str + '.nc')
        # Load ekf object and get data
        fof = Ekf(fof_fn, suppress_warnings=True)
        obs = fof.obs(**fkw)
        fc = fof.fg(**fkw)
        time = fof.obs(param='time', **fkw)
        plevel = fof.obs(param='level', **fkw)

        # Compute statistics
        diff = fc - obs

        if inargs.obs == 'SYNOP':
            results = compute_det_synop(diff, time)
        elif inargs.obs in ['TEMP', 'AIREP']:
            results = compute_det_air(diff, time, plevel, inargs.air_verif_time,
                                      inargs.obs)
        else:
            raise TypeError, 'Wrong obs type.'
    elif inargs.fc_type == 'ens':
        # Loop over ensemble members
        ens_list = []
        for ie in range(20):   # ATTENTION: Hard coded ensemble size
            fof_fn = (config.datadir + exp_id + '/' + date_str + '/ens' +
                      str(ie+1).zfill(3) + '/fof_' + date_str + '.nc')
            print fof_fn
            # Load ekf object and get data
            fof = Ekf(fof_fn, suppress_warnings=True)
            obs = fof.obs(**fkw)
            fc = fof.fg(**fkw)
            time = fof.obs(param='time', **fkw)
            plevel = fof.obs(param='level', **fkw)
            ens_list.append(fc)

        # compute statistics
        ens_mean = np.mean(ens_list, axis=0)
        ens_spread = np.std(ens_list, axis=0, ddof=1)

        if inargs.obs == 'SYNOP':
            print 'SYNOP data missing for ensembles...'
            results = compute_ens_synop(diff, time)
        elif inargs.obs in ['TEMP', 'AIREP']:
            results = compute_ens_air(ens_mean, ens_spread, obs, time, plevel,
                                      inargs.air_verif_time, inargs.obs)
        else:
            raise TypeError, 'Wrong obs type.'

    else:
        raise TypeError, 'Wrong fc_type.'

    return results


def get_verification_for_one_day(inargs, exp_id, date):
    """Computes verification for one day
    
    Args:
        inargs: Command line arguments
        exp_id: Experiment ID string
        date: datetime object

    Returns:
        results: det --> [rmse, bias ]
                 ens --> [rmse, std]
                 with [time] for SYNOP
                      [height] for TEMP/AIREP

    """

    # Create savestr
    verif_str = '_%ih-%ih' % tuple(inargs.air_verif_time) \
        if inargs.air_verif_time is not None else ''
    save_fn = (config.savedir_base + exp_id + '/fof_' + inargs.obs + '_' +
               inargs.var + '_' + inargs.fc_type + '_' +
               h.dt_to_yyyymmddhhmmss(date) + verif_str + '.npy')
    # Check if save_fn exists or recompute
    print 'Check if pre-computed file exists: %s' % save_fn
    if not os.path.exists(save_fn) or inargs.recompute:
        print 'Nope. Compute!'
        r = compute_verif(inargs, exp_id, date)
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
        results: Numpy array with dimensions
                 [metric, date, time] for SYNOP
                 [metric, date, height] for TEMP/AIREP
                 where metric is
                 [rmse, bias] for det
                 [rmse, spread] for ens
    """
    date_list = []
    # Loop over dates and collect data
    for date in h.make_timelist(inargs.date_start, inargs.date_stop,
                                inargs.hours_inc):
            date_list.append(get_verification_for_one_day(inargs, exp_id, date))

    # Change dimensions from [date, metric, time] to [metric, date, time]
    r = np.array(date_list)
    r = np.rollaxis(r, 1, 0)
    return r


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

    # Get ylabel
    if inargs.var == 'T2M': ylabel = 'T2M [K]'
    elif inargs.var == 'PS': ylabel = 'PS [Pa]'
    elif inargs.var == 'RH2M': ylabel = 'RH2M [%]'
    elif inargs.var == 'T': ylabel = 'T [K]'
    elif inargs.var in ['U10M', 'V10M']: ylabel = inargs.var + ' [m/s]'
    else: raise TypeError, 'Var not implemented.'

    ylabel = ylabel + ' rmse(-) and bias(--)' if inargs.fc_type == 'det' \
        else ylabel + ' rmse(-) and spread(--)'

    # Plot results
    if inargs.composite:
        plot_list = [np.nanmean(r, axis=1) for r in results_list]
        title_str = args.date_start[:-4] + '-' + args.date_stop[:-4]
        if inargs.obs == 'SYNOP':
            fig = h.plot_synop(plot_list, inargs.exp_id, title_str, ylabel)
        else:
            title_str += '_%ih-%ih' % tuple(inargs.air_verif_time)
            fig = h.plot_air(plot_list, inargs.exp_id, title_str, ylabel,
                             inargs.obs)

    else:
        for idate, date in enumerate(h.make_timelist(inargs.date_start,
                                                     inargs.date_stop,
                                                     inargs.hours_inc)):
            plot_list = [r[:, idate] for r in results_list]
            title_str = h.dt_to_yyyymmddhhmmss(date)[:-4]
            if inargs.obs == 'SYNOP':
                fig = h.plot_synop(plot_list, inargs.exp_id, title_str, ylabel)
            else:
                title_str += '_%ih-%ih' % tuple(inargs.air_verif_time)
                fig = h.plot_air(plot_list, inargs.exp_id, title_str, ylabel,
                                 inargs.obs)

    exp_id_str = '_'.join(inargs.exp_id)
    plot_dir = config.plotdir + '/' + exp_id_str + '/forecast_synop_air/'
    plot_str = (inargs.obs + '_' + inargs.var + '_' + inargs.fc_type + '_' +
                title_str)
    h.save_fig_and_log(fig, plot_str, plot_dir)


if __name__ == '__main__':
    # Get command line arguments
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--exp_id', metavar='exp_id', type=str, nargs='+',
                        help='Experiment ID')
    parser.add_argument('--date_start',
                        type=str,
                        default='20160527000000',
                        help='Start date for date loop (yyyymmddhhmmss)')
    parser.add_argument('--date_stop',
                        type=str,
                        default='20160609000000',
                        help='End date for date loop (yyyymmddhhmmss)')
    parser.add_argument('--hours_inc',
                        type=int,
                        default=24,
                        help='Time increment between forecasts (h)')
    parser.add_argument('--obs',
                        type=str,
                        help='SYNOP, TEMP or AIREP')
    parser.add_argument('--var',
                        type=str,
                        help='Variable name.')
    parser.add_argument('--air_verif_time',
                        type=int,
                        nargs='+',
                        help='Verification start and stop time in hours.')
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
