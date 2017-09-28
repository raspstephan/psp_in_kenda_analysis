"""
Compute and plot metrics from KENDA forecast files.

"""

# Imports
from argparse import ArgumentParser
import numpy as np
import os
import config
import helpers as h


# Intermediate functions
def compute_metric(inargs, exp_id, date):
    """Actually computes the metric for one exp_id and date.
    NaNs are handled on a daily basis
    
    Args:
        inargs: Command line arguments
        exp_id: Experiment ID string
        date: datetime object

    Returns:
        metric: Numpy array with dimensions [time, metric_dim]
    """
    date_str = h.dt_to_yyyymmddhhmmss(date)

    # Load presaved forecast data
    if exp_id == 'radar':   # Make an exception for radar as exp_id
        radar_fn = (config.savedir_base + 'radar/prec_fields/' + 'radar_' +
                    date_str + '.npy')
        fc_data = np.load(radar_fn)
    else:
        fc_fn = (config.savedir_base + exp_id + '/prec_fields/' +
                 config.metric_dict[inargs.metric.split('-')[0]]['det_or_ens'] +
                 '_' + date_str + '.npy')
        fc_data = np.load(fc_fn)

    if config.metric_dict[inargs.metric.split('-')[0]]['use_radar']:
        # Load presaved radar data
        radar_fn = (config.savedir_base + 'radar/prec_fields/' + 'radar_' +
                    date_str + '.npy')
        radar_data = np.load(radar_fn)
        radar_data, fc_data = h.handle_nans(radar_data, fc_data,
                                            inargs.radar_thresh)
    # Pass data to computation functions
    if inargs.metric == 'det_mean_prec':
        m = h.compute_det_mean_prec(fc_data)
    elif inargs.metric == 'det_rmse':
        m = h.compute_det_rmse(radar_data, fc_data)
    elif 'det_fss' in inargs.metric:
        # Parse
        _, fss_thresh, fss_size = inargs.metric.split('-')
        fss_thresh = float(fss_thresh) / 10.
        fss_size = int(fss_size)
        # Update dictionary
        config.metric_dict[inargs.metric.split('-')[0]]['ylabel'] = \
            'FSS ' + str(fss_thresh) + 'mm/h ' + str(fss_size * 2.8) + 'km'
        m = h.compute_det_fss(radar_data, fc_data, fss_thresh, fss_size)
    elif inargs.metric == 'det_prec_hist':
        m = h.compute_det_prec_hist(fc_data)
    elif inargs.metric == 'ens_crps':
        m = h.compute_ens_crps(radar_data, fc_data)
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
               h.dt_to_yyyymmddhhmmss(date) + '.npy')

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


# Plotting functions
def plot_panel(inargs, plot_list, title_str):
    """
    
    Args:
        inargs: Command line arguments
        plot_list: List of metrics [exp_id][time, metric_dim]
    """

    if config.metric_dict[inargs.metric.split('-')[0]]['plot_type'] == 'line':
        fig = h.plot_line(plot_list, inargs.exp_id, inargs.metric.split('-')[0],
                          title_str)
    elif config.metric_dict[inargs.metric.split('-')[0]]['plot_type'] == 'hist':
        fig = h.plot_hist(plot_list, inargs.exp_id, inargs.metric.split('-')[0],
                          title_str, inargs.hist_normalize)
    else:
        raise ValueError('Plot type %s does not exist.' %
                         config.metric_dict[inargs.metric.split('-')[0]]
                         ['plot_type'])
    exp_id_str = '_'.join([e for e in inargs.exp_id if not e == 'radar'])
    plot_dir = config.plotdir + '/' + exp_id_str + '/forecast_metrics/'
    plot_str = inargs.metric + '_' + title_str
    h.save_fig_and_log(fig, plot_str, plot_dir)


# Main program
def main(inargs):
    """Main program
    
    Args:
        inargs: Command line arguments
    """
    # Loop over exp_ids and collect data in list
    results_list = []
    for e in inargs.exp_id:
        results_list.append(get_metric_for_all_dates(inargs, e))

    if inargs.composite:
        plot_list = [np.nanmean(r, axis=0) for r in results_list]
        title_str = args.date_start[:-4] + '-' + args.date_stop[:-4]
        plot_panel(inargs, plot_list, title_str)

    else:
        for idate, date in enumerate(h.make_timelist(inargs.date_start,
                                                     inargs.date_stop,
                                                     inargs.hours_inc)):
            plot_list = [r[idate] for r in results_list]
            title_str = h.dt_to_yyyymmddhhmmss(date)[:-4]
            plot_panel(inargs, plot_list, title_str)


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
    parser.add_argument('--radar_thresh',
                        type=float,
                        default=100.,
                        help='Radar values above threshold will be set to nan.'
                             'Default = 100.')
    parser.add_argument('--hist_normalize',
                        dest='hist_normalize',
                        action='store_true',
                        help='Normalize histogram by all precipitation points')
    parser.set_defaults(hist_normalize=False)

    args = parser.parse_args()

    main(args)
