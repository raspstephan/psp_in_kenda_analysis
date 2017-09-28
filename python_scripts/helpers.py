"""
Helper functions for my KENDA scripts

"""
import sys
from datetime import datetime, timedelta
from subprocess import check_output
from git import Repo
from cosmo_utils.pywgrib import getfobj_ens, getfobj
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, ddhhmmss_strtotime, \
    yymmddhhmm
from config import *  # Import config file
from cosmo_utils.pyncdf import getfobj_ncdf
import numpy as np
import matplotlib.pyplot as plt
from cosmo_utils.scores.probab import FSS
sys.path.append('/home/s/S.Rasp/repositories/enstools/')
from enstools.scores import crps_sample


def save_fig_and_log(fig, fig_name, plot_dir):
    """
    Save the given figure along with a log file
    """

    # Step 1: save figure
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    print('Saving figure: %s' % (plot_dir + '/' + fig_name + '.pdf'))
    fig.savefig(plot_dir + '/' + fig_name + '.pdf')

    # Step 2: Create and save log file
    time_stamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    pwd = check_output(['pwd']).rstrip()  # Need to remove trailing /n
    git_dir = pwd.rsplit('/', 1)[0]
    git_hash = Repo(git_dir).heads[0].commit
    exe_str = ' '.join(sys.argv)
    s = check_output(['conda', 'env', 'list'])
    for l in s.split('\n'):
        if '*' in l:
            py_env = l
    assert 'py_env' in locals(), 'No active conda environemnt found.'

    log_str = ("""
Time: %s\n
Executed command:\n
python %s\n
In directory: %s\n
Git hash: %s\n
Anaconda environment: %s\n
    """ % (time_stamp, exe_str, pwd, str(git_hash), py_env))

    logf = open(plot_dir + '/' + fig_name + '.log', 'w+')
    logf.write(log_str)
    logf.close()


radarpref = 'raa01-rw_10000-'
radarsufx = '-dwd---bin.nc'
precsuf = '_15'
precsuf_da = '_prec'
gribpref = 'lfff'
gribpref_da = 'lff'
nens = 20
nens_da = 40


# Define loading functions
def load_det(datadir, date, t, return_array=False):
    topdir = datadir + '/' + date + '/'
    gribfn = gribpref + t + precsuf
    detfn = topdir + 'det/' + gribfn
    detfobj = getfobj(detfn, fieldn='PREC_PERHOUR')
    # Minus one hour
    # gribfnm1 = gribpref + ddhhmmss(ddhhmmss_strtotime(t) -
    # timedelta(hours = 1)) + precsuf
    # detfnm1 = topdir + 'det/' + gribfnm1
    # detfobjm1 = getfobj(detfnm1, fieldn = 'TOT_PREC')
    # detfobj.data = detfobj.data - detfobjm1.data
    if return_array:
        return detfobj.data
    else:
        return detfobj


def load_radar(date, t='00000000', return_array=False):
    dateobj = (yyyymmddhhmmss_strtotime(date) + ddhhmmss_strtotime(t))
    radardt = timedelta(minutes=10)  # TODO Is this correct???
    radardateobj = dateobj - radardt
    radarfn = radardir + radarpref + yymmddhhmm(radardateobj) + radarsufx
    radarfobj = getfobj_ncdf(radarfn, fieldn='pr', dwdradar=True)
    if return_array:
        return radarfobj.data
    else:
        return radarfobj


def load_ens(datadir, date, t, return_array=False):
    topdir = datadir + '/' + date + '/'
    gribfn = gribpref + t + precsuf
    ensfobjlist = getfobj_ens(topdir, 'sub', mems=nens, gribfn=gribfn,
                              dir_prefix='ens', fieldn='PREC_PERHOUR',
                              para=4)
    if return_array:
        return [fobj.data for fobj in ensfobjlist]
    else:
        return ensfobjlist


def load_det_da(datadir, t, return_array=False):
    detfn = datadir + gribpref_da + t + precsuf_da + '.det'
    detfobj = getfobj(detfn, fieldn='TOT_PREC_S')
    if return_array:
        return detfobj.data
    else:
        return detfobj

def load_ens_da(datadir, t, return_array=False):
    gribfn = gribpref_da + t + precsuf_da
    ensfobjlist = getfobj_ens(datadir, 'same', mems=nens, gribfn=gribfn,
                              fieldn='TOT_PREC_S', para=4)
    if return_array:
        return [fobj.data for fobj in ensfobjlist]
    else:
        return ensfobjlist


def strip_expid(expid):
    return expid.replace('DA_', '').replace('_ens', '').replace('v2', ''). \
        replace('_2JUN', '')


def set_plot(ax, title, args, hourlist_plot, adjust=True):
    plt.sca(ax)
    ax.set_xlabel('Time [UTC]')
    ax.legend(loc=0, fontsize=8, frameon=False)
    ymax = np.ceil(ax.get_ylim()[1] * 10) / 10.
    ax.set_ylim(0, ymax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    ax.spines['bottom'].set_position(('outward', 3))
    try:
        ax.set_xticks(range(args.hint + 1)[::6])
        ax.set_xticklabels(hourlist_plot[::6])
        ax.set_xlim(0, 24)
    except AttributeError:   # Must be DA plots
        ax.set_xticks(range(24)[::6])
        ax.set_xlim(0, 24)
    ax.set_title(title)
    if adjust:
        plt.subplots_adjust(bottom=0.18, left=0.18, right=0.97)


def set_panel(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    ax.spines['bottom'].set_position(('outward', 3))
    ax.set_xticks(range(24)[::6])
    ax.set_xlim(0, 24)
    ax.set_xlabel('Time [UTC]')


def compute_det_stats(nanradar, nandet, convradar, convdet):
    """
    Compute RMSE and FSS for determinsitc forecast
    """
    rmse = np.sqrt(np.nanmean((convradar - convdet) ** 2))
    fss01 = FSS(0.1, 21, nanradar, nandet, python_core=True)
    fss10 = FSS(1.0, 21, nanradar, nandet, python_core=True)

    return rmse, fss01, fss10


def compute_prec_hist(nanfield, bin_edges):
    """
    Compute precipitation histograms.
    """
    a = np.ravel(nanfield)
    a = a[np.isfinite(a)]
    hist = np.histogram(a, bin_edges)[0]
    return hist


def compute_ens_stats(convradar, convfieldlist, ens_norm_type,
                      norm_thresh=0.1, bs_threshold=1.):
    """
    convfieldlist dimensions = (ens, x, y)
    Compute spread and rmse of ensemble with given normalization type.
    0: no normalization
    1: grid point normalization
    2: domain normalization
    """
    meanfield = np.mean(convfieldlist, axis=0)
    ensradarmean = 0.5 * (convradar + meanfield)
    if ens_norm_type == 0:  # No normalization
        spread = np.nanmean(np.std(convfieldlist, axis=0, ddof=1))
        rmv = np.sqrt(np.nanmean(np.var(convfieldlist, axis=0, ddof=1)))
        rmse = np.sqrt(np.nanmean((convradar - meanfield) ** 2))
        mean = np.nanmean(meanfield)
        mean_std = np.std(np.nanmean(convfieldlist, axis=(1, 2)), ddof=1)
        bs = compute_bs(convradar, convfieldlist, threshold=bs_threshold)

        results = (rmse, rmv, bs, mean, mean_std)
    elif ens_norm_type == 1:   # Grid point normalization
        spread = np.nanmean((np.std(convfieldlist, axis=0, ddof=1) /
                             meanfield)[meanfield >= norm_thresh])
        rmse = np.sqrt(np.nanmean(((convradar - meanfield) ** 2 /
                                   ensradarmean ** 2)[ensradarmean >= 0.1]))
        results = (spread, rmse)

    elif ens_norm_type == 2:   # Domain normalization
        spread = (np.nanmean(np.std(convfieldlist, axis=0, ddof=1)
                             [meanfield >= 0.1]) /
                  np.nanmean(meanfield[meanfield >= 0.1]))
        rmse = (np.sqrt(np.nanmean(((convradar - meanfield) ** 2)
                                    [ensradarmean >= 0.1])) /
                np.nanmean(ensradarmean[ensradarmean >= 0.1]))
        results = (spread, rmse)

    else:
        raise Exception, 'Wrong ens_norm_type'

    return results


def compute_crps(obs, enslist):
    """
    Compute the sample CRPS gripointwise and then average
    """
    obs_flat = np.ravel(obs)
    mask = np.isfinite(obs_flat)
    obs_flat = obs_flat[mask]
    ens_flat = np.array([np.ravel(e)[mask] for e in enslist])
    crps = crps_sample(obs_flat, ens_flat)
    return np.mean(crps)


def compute_bs(obs, enslist, threshold):
    """
    Compute the Brier score for a given threshold.
    """
    prob_fc = np.sum((enslist > threshold), axis=0) / enslist.shape[0]
    bin_obs = obs > threshold
    bs = np.nanmean((prob_fc - bin_obs) ** 2)
    return bs


def make_fig_fc_ens(args, x, rd, title, it=None):
    """
    Plot ensemble panels: 
    1. Mean precipitation with ensemble error bars
    2. RMSE and RMV
    3. BS
    """

    aspect = 0.4
    # Set up figure
    fig, axes = plt.subplots(1, 3, figsize=(pw, aspect * pw))

    # Loop over experiments
    for ie, expid in enumerate(args.expid):
        if ie == 0:
            if it is None:
                it = np.ones(rd[expid]['radar'].shape[0], dtype=bool)
            axes[0].plot(x, rd[expid]['radar'][it], c='k', lw=2, label='Radar')
        # Plot mean with error bars in first panel
        axes[0].errorbar(x, rd[expid]['ensmean'][it], yerr=rd[expid]['ensmean_std'], c=cdict[expid],
                         label=expid)

        # Plot RMSE and RMV in second plot
        axes[1].plot(x, rd[expid]['ensrmse'][it], lw=1.5, c=cdict[expid])
        axes[1].plot(x, rd[expid]['ensrmv'][it], lw=1.5, c=cdict[expid], ls='--')

        # Plot Brier Score in thrid panel
        axes[2].plot(x, rd[expid]['ensbs'][it], lw=1.5, c=cdict[expid])

    # Define labels
    axes[0].set_title('Mean precip pm std')
    axes[0].set_ylabel('[mm/h]')
    axes[1].set_title('RMSE and RMV')
    axes[1].set_ylabel('[mm/h]')
    axes[2].set_title('Brier Score')

    # Adjust spines and x-label
    for ax in axes:
        set_panel(ax)

    axes[0].legend(loc=0, fontsize=5, frameon=False)
    fig.suptitle(title, fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, 0.95))

    return fig


def make_fig_fc_det_rmse(args, x, rd, title, it=None):
    """
    Plot deterministic panels.
    1. Mean precipitation
    2. RMSE
    3. FSS
    """
    aspect = 0.4
    # Set up figure
    fig, axes = plt.subplots(1, 3, figsize=(pw, aspect * pw))

    # Loop over experiments
    for ie, expid in enumerate(args.expid):
        if ie == 0:
            if it is None:
                it = np.ones(rd[expid]['radar'].shape[0], dtype=bool)
            axes[0].plot(x, rd[expid]['radar'][it], c='k', lw=2, label='Radar')

        axes[0].plot(x, rd[expid]['detmean'][it], c=cdict[expid], label=expid)

        axes[1].plot(x, rd[expid]['detrmse'][it], c=cdict[expid])

        axes[2].plot(x, rd[expid]['fss10'][it], c=cdict[expid])

    # Define labels
    axes[0].set_title('Mean precip')
    axes[0].set_ylabel('[mm/h]')
    axes[1].set_title('RMSE')
    axes[1].set_ylabel('[mm/h]')
    axes[2].set_title('FSS')
    axes[2].set_ylabel('[1 mm/h, 60 km]')

    # Adjust spines and x-label
    for ax in axes:
        set_panel(ax)

    axes[0].legend(loc=0, fontsize=5, frameon=False)
    fig.suptitle(title, fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, 0.95))

    return fig


def make_fig_fc_det_fss(args, x, rd, title, it=None):
    """
    Plot deterministic panels.
    1. Mean precipitation
    2. FSS 0.1
    3. FSS 1.0
    """
    aspect = 0.4
    # Set up figure
    fig, axes = plt.subplots(1, 3, figsize=(pw, aspect * pw))

    # Loop over experiments
    for ie, expid in enumerate(args.expid):
        if ie == 0: 
            if it is None:
                it = np.ones(rd[expid]['radar'].shape[0], dtype=bool)
            axes[0].plot(x, rd[expid]['radar'][it], c='k', lw=2, label='Radar')

        axes[0].plot(x, rd[expid]['detmean'][it], c=cdict[expid], label=expid)

        axes[1].plot(x, rd[expid]['fss01'][it], c=cdict[expid])

        axes[2].plot(x, rd[expid]['fss10'][it], c=cdict[expid])

    # Define labels
    axes[0].set_title('Mean precip')
    axes[0].set_ylabel('[mm/h]')
    axes[1].set_title('FSS')
    axes[1].set_ylabel('[0.1 mm/h, 60 km]')
    axes[2].set_title('FSS')
    axes[2].set_ylabel('[1 mm/h, 60 km]')

    # Adjust spines and x-label
    for ax in axes:
        set_panel(ax)

    axes[0].legend(loc=0, fontsize=5, frameon=False)
    fig.suptitle(title, fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, 0.95))

    return fig


def make_fig_hist(args, rd, title, it=None):
    """
    Plot precipitation histogram of different experiments.
    """

    x = np.arange(bin_edges.shape[0] - 1)
    fig, ax = plt.subplots(1, 1, figsize=(pw/2, pw/2))

    for ie, expid in enumerate(args.expid):
        if ie == 0: 
            if it is None:
                it = np.ones(rd[expid]['radar_hist'].shape[0], dtype=bool)
            ax.bar(x, rd[expid]['radar_hist'][it], 0.1, color='k', 
                   label='Radar')

        ax.bar(x + 0.1*(ie+1), rd[expid]['hist'][it], 0.1, 
               color=cdict[expid], label=expid)
    ax.set_xticks(range(bin_edges.shape[0]))
    ax.set_xticklabels(['%.1f'%i for i in list(bin_edges)])
    ax.set_yscale('log')
    ax.set_xlabel('mm/h')
    ax.legend(loc=0, fontsize=5)

    plt.title(title)
    plt.tight_layout()

    return fig


def make_timelist(date_start, date_stop, hours_inc):
    """
    
    Args:
        date_start: yyyymmddhhmmss time string
        date_stop: yyyymmddhhmmss time string
        hour_inc: increment in h

    Returns:
        timelist: List with datetime objects
    """
    dt_start = yyyymmddhhmmss_to_dt(date_start)
    dt_stop = yyyymmddhhmmss_to_dt(date_stop)
    td_inc = timedelta(hours=hours_inc)

    timelist = []
    t = dt_start
    while t <= dt_stop:
        timelist.append(t)
        t += td_inc
    return timelist


def yyyymmddhhmmss_to_dt(yyyymmddhhmmss):
    f = '%Y%m%d%H%M%S'
    return datetime.strptime(yyyymmddhhmmss, f)


def dt_to_yyyymmddhhmmss(dt):
    f = '%Y%m%d%H%M%S'
    return datetime.strftime(dt, f)


def handle_nans(radar_data, fc_data, radar_thresh):
    """Handle NaNs on a daily basis.
    
    Args:
        radar_data: Radar data array with dimensions [time, x, y]
        fc_data: Forecast data array with dimensions [time, x, y] or 
                 [time, ens, x, y]
        radar_thresh: Threshold, NaNs above.

    Returns:
        radar_data, fc_data: Same arrays with NaNs 
    """
    mask = np.max(radar_data, axis=0) > radar_thresh
    radar_data[:, mask] = np.nan
    if fc_data.ndim == 3:
        fc_data[:, mask] = np.nan
    else:
        fc_data[:, :, mask] = np.nan
    return radar_data, fc_data


# New compute_*metric* functions
# These work with input of dimension [hour, x, y] or [hour, ens, x, y]
def compute_det_rmse(radar_data, fc_data):
    """Compute deterministic rmse
    
    Args:
        radar_data: Radar data
        fc_data: forecast data

    Returns:
        rmse: Numpy array with dimensions [hour]
    """
    rmse = np.sqrt(np.nanmean((radar_data - fc_data) ** 2, axis=(1, 2)))
    return rmse


def compute_det_mean_prec(data):
    """Compute deterministic mean precipitation

    Args:
        data: data

    Returns:
        mean: Numpy array with dimensions [hour]
    """
    mean = np.nanmean(data, axis=(1, 2))
    return mean


# Panel plotting functions
def plot_line(plot_list, exp_ids, metric, title):
    """
    
    Args:
        plot_list: List of metrics [exp_id][time, metric_dim]

    Returns:
        fig: Figure object
    """

    fig, ax = plt.subplots(1, 1, figsize=(0.5 * pw, 0.5 * pw))

    x = np.arange(1, 25)
    for ie, e in enumerate(exp_ids):
        ax.plot(x, plot_list[ie], c=cdict[e], label=e)

    ax.set_xlabel('Forecast lead time [h]')
    ax.set_ylabel(metric_dict[metric]['ylabel'])
    ax.set_title(title)
    ax.legend(loc=0, fontsize=6)

    set_panel(ax)

    plt.tight_layout()

    return fig




