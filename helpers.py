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
from config import *   # Import config file
from cosmo_utils.pyncdf import getfobj_ncdf


def save_fig_and_log(fig, fig_name, plot_dir):
    """
    Save the given figure along with a log file
    """

    # Step 1: save figure
    print('Saving figure: %s' % (plot_dir + '/' + fig_name + '.pdf'))
    fig.savefig(plot_dir + '/' + fig_name + '.pdf')

    # Step 2: Create and save log file
    time_stamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    pwd = check_output(['pwd']).rstrip()   # Need to remove trailing /n
    git_hash = Repo(pwd).heads[0].commit
    exe_str = ' '.join(sys.argv)

    log_str = ("""
Time: %s\n
Executed command:\n
python %s\n
In directory: %s\n
Git hash: %s\n
    """ % (time_stamp, exe_str, pwd, str(git_hash)))

    logf = open(plot_dir + '/' + fig_name + '.log', 'w+')
    logf.write(log_str)
    logf.close()


radarpref = 'raa01-rw_10000-'
radarsufx = '-dwd---bin.nc'
precsuf = '_15'
gribpref = 'lfff'
nens = 20


# Define loading functions
def load_det(datadir, date, t):
    topdir = datadir + '/' + date + '/'
    gribfn = gribpref + t + precsuf
    detfn = topdir + 'det/' + gribfn
    detfobj = getfobj(detfn, fieldn = 'PREC_PERHOUR')
    # Minus one hour
    #gribfnm1 = gribpref + ddhhmmss(ddhhmmss_strtotime(t) -
                                 #timedelta(hours = 1)) + precsuf
    #detfnm1 = topdir + 'det/' + gribfnm1
    #detfobjm1 = getfobj(detfnm1, fieldn = 'TOT_PREC')
    #detfobj.data = detfobj.data - detfobjm1.data
    return detfobj


def load_radar(date, t):
    dateobj = (yyyymmddhhmmss_strtotime(date) + ddhhmmss_strtotime(t))
    radardt = timedelta(minutes = 10)   # TODO Is this correct???
    radardateobj = dateobj - radardt
    radarfn = radardir + radarpref + yymmddhhmm(radardateobj) + radarsufx
    radarfobj = getfobj_ncdf(radarfn, fieldn = 'pr', dwdradar = True)
    return radarfobj


def load_ens(datadir, date, t):
    topdir = datadir + '/' + date + '/'
    gribfn = gribpref + t + precsuf
    ensfobjlist = getfobj_ens(topdir, 'sub', mems = nens, gribfn = gribfn,
                              dir_prefix = 'ens', fieldn = 'PREC_PERHOUR',
                              para = 4)
    return ensfobjlist


def strip_expid(expid):
    return expid.replace('DA_', '').replace('_ens', '').replace('v2', '').\
        replace('_2JUN', '')
