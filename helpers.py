"""
Helper functions for my KENDA scripts

"""
import sys
from datetime import datetime
from subprocess import check_output
from git import Repo


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