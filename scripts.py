# Imports
import argparse
import os
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, make_timelist, \
    yyyymmddhhmmss
from datetime import timedelta

# Arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--tstart', metavar = 'tstart', type=str, 
                    default = '20160525000000')
parser.add_argument('--tend', metavar = 'tend', type=str, 
                    default = '20160610000000')
parser.add_argument('--ver_start', metavar = 'ver_start', type=str, 
                    default = '600')
parser.add_argument('--ver_end', metavar = 'ver_end', type=str, 
                    default = '840')
parser.add_argument('--var', metavar = 'var', type=str, default = 'T')
parser.add_argument('--obs', metavar = 'obs', type=str, default = 'TEMP')
parser.add_argument('--hint', metavar = 'hint', type=str, default = '24')
parser.add_argument('--plot', metavar = 'plot', type=str, nargs = '+')
args = parser.parse_args()



tstart = yyyymmddhhmmss_strtotime('20160525000000')
tend = yyyymmddhhmmss_strtotime('20160610000000')
tint = timedelta(days = 1)
timelist = make_timelist(tstart, tend, tint)
ver_start = '600'
ver_end = '840'
var = 'T2M'
obs = 'SYNOP'
hint = '48'

for t in timelist:
    if 'prec_time' in args.plot:
        print('python plot_stamps.py --date ' + yyyymmddhhmmss(t) + 
            ' --time 1 '  + hint + ' --plot prec_time --expid REF REF_TL500 PSP_TL500')
        os.system('python plot_stamps.py --date ' + yyyymmddhhmmss(t) + 
                ' --time 1 '  + hint + ' --plot prec_time --expid REF REF_TL500 PSP_TL500')
    if 'prec_comp' in args.plot:
        print('python plot_stamps.py --date ' + yyyymmddhhmmss(t) + 
            ' --time 1 '  + hint + ' --plot prec_comp --expid REF REF_TL500 PSP_TL500')
        os.system('python plot_stamps.py --date ' + yyyymmddhhmmss(t) + 
                ' --time 1 '  + hint + ' --plot prec_comp --expid REF REF_TL500 PSP_TL500')
    
    if 'verif' in args.plot:
        # T verif for each day individually
        for e in ['REF', 'REF_TL500', 'PSP_TL500']:
            print('python verif_fof.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
                yyyymmddhhmmss(t) + ' --date_end ' + yyyymmddhhmmss(t) + ' --var ' + var +  
                ' --expid ' + e + ' --obs ' + obs)
            os.system('python verif_fof.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
                yyyymmddhhmmss(t) + ' --date_end ' + yyyymmddhhmmss(t) + ' --var ' + var +  
                ' --expid ' + e + ' --obs ' + obs)
        
        # Comparison plot
        allexp = 'REF REF_TL500 PSP_TL500'
        print('python compare_verif.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
                yyyymmddhhmmss(t) + ' --date_end ' + yyyymmddhhmmss(t) + ' --var ' + var +  
                ' --expid ' + allexp + ' --obs ' + obs)
        os.system('python compare_verif.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
                yyyymmddhhmmss(t) + ' --date_end ' + yyyymmddhhmmss(t) + ' --var ' + var +  
                ' --expid ' + allexp + ' --obs ' + obs)

# Average prec 
if 'prec_time' in args.plot:
    allexp = 'REF REF_TL500 PSP_TL500'
    print('python average_prec.py --date_ini '+ yyyymmddhhmmss(tstart) + 
        ' --date_end '+ yyyymmddhhmmss(tend) +' --time 1 '  + hint + ' --expid ' + allexp)
    os.system('python average_prec.py --date_ini '+ yyyymmddhhmmss(tstart) + 
        ' --date_end '+ yyyymmddhhmmss(tend) +' --time 1 '  + hint + ' --expid ' + allexp)


if 'verif' in args.plot:
    # T verif for entire period
    for e in ['REF', 'REF_TL500', 'PSP_TL500']:
        print('python verif_fof.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
            yyyymmddhhmmss(tstart) + ' --date_end ' + yyyymmddhhmmss(tend) + ' --var ' + var +  
            ' --expid ' + e + ' --obs ' + obs)
        os.system('python verif_fof.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
            yyyymmddhhmmss(tstart) + ' --date_end ' + yyyymmddhhmmss(tend) + ' --var ' + var +  
            ' --expid ' + e + ' --obs ' + obs)
            
    allexp = 'REF REF_TL500 PSP_TL500'
    print('python compare_verif.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
            yyyymmddhhmmss(tstart) + ' --date_end ' + yyyymmddhhmmss(tend) + ' --var ' + var +  
            ' --expid ' + allexp + ' --obs ' + obs)
    os.system('python compare_verif.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
            yyyymmddhhmmss(tstart) + ' --date_end ' + yyyymmddhhmmss(tend) + ' --var ' + var +  
            ' --expid ' + allexp + ' --obs ' + obs)
