"""
Compute and plot metrics from KENDA forecast files.

"""

# Imports
from argparse import ArgumentParser
import numpy as np
import os
import config
import helpers as h



def compute_metric(inargs, exp_id, date):
    """Actually computes the metric for one exp_id and date.
    NaNs are handled on a daily basis
    
    Args:
        inargs: 
        exp_id: 
        date: 

    Returns:
        metric: Numpy array with dimensions [time, metric_dim]
    """

    # Load presaved forecast data
    fc_data = ???
    # Load presaved radar data
    radar_data = ???

    # Handle NaNs
    radar_data, fc_data = h.handle_nans(radar_data, fc_data)

    # Pass data to
    if inargs.metric == 'rmse':
        m = h.compute_rmse(radar_data, fc_data)
    else:
        raise ValueError('Metric %s does not exist.' % inargs.metric)

    return m




def get_metric_for_one_day(inargs, exp_id, date):
    """Returns the requested metric for one day and one exp_id.
    Saves the data as a numpy array and loads the pre-computed files if it
    exists.
    
    Args:
        inargs: Command line arguments
        exp_id: Experiment ID string
        date: datetime object

    Returns:
        results: Numpy array with dimensions [time, metric_dim]
    """

    # Create savestr
    save_fn = (config.savedir_base + exp_id + '/' + inargs.metric + '_' +
               h.dt_to_yyyymmddhhmmss(date))

    # Check if save_fn exists or recompute
    print 'Check if pre-computed file exists: %s' % save_fn
    if not os.path.exists(save_fn) or inargs.recompute:
        print 'Nope. Compute!'
        r = compute_metric(inargs, exp_id, date)
        np.save(save_fn, r)
    else:
        print 'Yep. Load!'
        r = np.load(save_fn)
    return r


def get_metric_for_all_dates(inargs, exp_id):
    """Get the requested metric for all dates
    
    Args:
        inargs: Command line arguments
        exp_id: Experiment ID string

    Returns:
        results: Numpy array with dimensions [date, time, metric_dim]
    """
    date_list = []
    # Loop over dates and collect data
    for date in h.make_timelist(inargs.date_start, inargs.date_stop,
                                inargs.hours_inc):
        date_list.append(get_metric_for_one_day(args, exp_id, date))
    return np.array(date_list)


# Main program
def main(inargs):
    """
    Parameters:
        
    """
    # Loop over exp_ids and collect data in list
    results_list = []
    for e in inargs.exp_id:
        results_list.append(get_metric_for_all_dates(inargs, e))

    if inargs.composite:
        plot_list = [np.nanmean(r, axis=0) for r in results_list]
        plot_panel(plot_list)

    else:
        for idate, date in enumerate(???):
            plot_list = [r[idate] for r in results_list]
            plot_panel(plot_list)


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
    parser.add_argument('--hours_inc',
                        type=int,
                        default='24',
                        help='Time increment between forecasts (h)')
    parser.add_argument('--metric',
                        type=str,
                        help='Metric to be computed and plotted.')
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