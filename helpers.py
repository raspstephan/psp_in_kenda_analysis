"""
Helper functions for my KENDA scripts

"""
import sys
from datetime import datetime
from subprocess import check_output
from git import Repo
from cosmo_utils.pywgrib import getfobj_ens, getfobj
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, ddhhmmss_strtotime, \
    yymmddhhmm
from datetime import timedelta
from config import *  # Import config file
from cosmo_utils.pyncdf import getfobj_ncdf
import numpy as np
import matplotlib.pyplot as plt
from cosmo_utils.scores.probab import FSS


def save_fig_and_log(fig, fig_name, plot_dir):
    """
    Save the given figure along with a log file
    """

    # Step 1: save figure
    print('Saving figure: %s' % (plot_dir + '/' + fig_name + '.pdf'))
    fig.savefig(plot_dir + '/' + fig_name + '.pdf')

    # Step 2: Create and save log file
    time_stamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    pwd = check_output(['pwd']).rstrip()  # Need to remove trailing /n
    git_hash = Repo(pwd).heads[0].commit
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


def compute_bs(obs, enslist, threshold):
    """
    Compute the Brier score for a given threshold.
    """
    prob_fc = np.sum((enslist > threshold), axis=0) / enslist.shape[0]
    bin_obs = obs > threshold
    bs = np.nanmean((prob_fc - bin_obs) ** 2)
    return bs


def make_fig_fc_ens(args, x, radar, ensmean, ensmean_std, rmse, rmv, bs, title):
    """
    Plot ensemble panels: 
    1. Mean precipitation with ensemble error bars
    2. RMSE and RMV
    3. BS
    """

    aspect = 0.4
    # Set up figure
    fig, axes = plt.subplots(1, 3, figsize=(pw, aspect * pw))

    # Plot radar in first panel
    axes[0].plot(x, radar[0], c='k', lw=2, label='Radar')

    # Loop over experiments
    for ie, expid in enumerate(args.expid):

        # Plot mean with error bars in first panel
        axes[0].errorbar(x, ensmean[ie], yerr=ensmean_std[ie], c=cdict[expid],
                         label=expid)

        # Plot RMSE and RMV in second plot
        axes[1].plot(x, rmse[ie], lw=1.5, c=cdict[expid])
        axes[1].plot(x, rmv[ie], lw=1.5, c=cdict[expid], ls='--')

        # Plot Brier Score in thrid panel
        axes[2].plot(x, bs[ie], lw=1.5, c=cdict[expid])

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


def make_fig_fc_det_rmse(args, x, radar, meanprec, rmse, fss, title):
    """
    Plot deterministic panels.
    1. Mean precipitation
    2. RMSE
    3. FSS
    """
    aspect = 0.4
    # Set up figure
    fig, axes = plt.subplots(1, 3, figsize=(pw, aspect * pw))

    # Plot radar in first panel
    axes[0].plot(x, radar[0], c='k', lw=2, label='Radar')

    # Loop over experiments
    for ie, expid in enumerate(args.expid):

        # Mean precipitation in first panel
        axes[0].plot(x, meanprec[ie], c=cdict[expid], label=expid)

        # RMSE in second panel
        axes[1].plot(x, rmse[ie], c=cdict[expid])

        # FSS in third panel
        axes[2].plot(x, fss[ie], c=cdict[expid])

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


def make_fig_fc_det_fss(args, x, radar, meanprec, fss01, fss10, title):
    """
    Plot deterministic panels.
    1. Mean precipitation
    2. FSS 0.1
    3. FSS 1.0
    """
    aspect = 0.4
    # Set up figure
    fig, axes = plt.subplots(1, 3, figsize=(pw, aspect * pw))

    # Plot radar in first panel
    axes[0].plot(x, radar[0], c='k', lw=2, label='Radar')

    # Loop over experiments
    for ie, expid in enumerate(args.expid):

        # Mean precipitation in first panel
        axes[0].plot(x, meanprec[ie], c=cdict[expid], label=expid)

        # RMSE in second panel
        axes[1].plot(x, fss01[ie], c=cdict[expid])

        # FSS in third panel
        axes[2].plot(x, fss10[ie], c=cdict[expid])

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