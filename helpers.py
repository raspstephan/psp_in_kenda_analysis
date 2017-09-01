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
def load_det(datadir, date, t):
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
    return detfobj


def load_radar(date, t='00000000'):
    dateobj = (yyyymmddhhmmss_strtotime(date) + ddhhmmss_strtotime(t))
    radardt = timedelta(minutes=10)  # TODO Is this correct???
    radardateobj = dateobj - radardt
    radarfn = radardir + radarpref + yymmddhhmm(radardateobj) + radarsufx
    radarfobj = getfobj_ncdf(radarfn, fieldn='pr', dwdradar=True)
    return radarfobj


def load_ens(datadir, date, t):
    topdir = datadir + '/' + date + '/'
    gribfn = gribpref + t + precsuf
    ensfobjlist = getfobj_ens(topdir, 'sub', mems=nens, gribfn=gribfn,
                              dir_prefix='ens', fieldn='PREC_PERHOUR',
                              para=4)
    return ensfobjlist


def load_det_da(datadir, t):
    detfn = datadir + gribpref_da + t + precsuf_da + '.det'
    detfobj = getfobj(detfn, fieldn='TOT_PREC_S')
    return detfobj


def load_ens_da(datadir, t):
    gribfn = gribpref_da + t + precsuf_da
    ensfobjlist = getfobj_ens(datadir, 'same', mems=nens, gribfn=gribfn,
                              fieldn='TOT_PREC_S', para=4)
    return ensfobjlist


def strip_expid(expid):
    return expid.replace('DA_', '').replace('_ens', '').replace('v2', ''). \
        replace('_2JUN', '')


def set_plot(ax, title, args, hourlist_plot):
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
    plt.subplots_adjust(bottom=0.18, left=0.18, right=0.97)
