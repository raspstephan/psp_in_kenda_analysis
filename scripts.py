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
parser.add_argument('--ind_days', metavar = 'ind_days', type=str, default = 'True')
args = parser.parse_args()



tstart = yyyymmddhhmmss_strtotime(args.tstart)
tend = yyyymmddhhmmss_strtotime(args.tend)
tint = timedelta(days = 1)
timelist = make_timelist(tstart, tend, tint)
ver_start = args.ver_start
ver_end = args.ver_end
var = args.var
obs = args.obs
hint = args.hint

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
    
    if 'verif' in args.plot and args.ind_days == 'True':
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
    
if 'bl_evo' in args.plot:
    for t in range(12,24):
        print('python bl_evo.py --expid DA_REF DA_PSP DA_REF_TL500 DA_PSP_TL500 DA_PSPv2_TL500 --date 20160604000000 --time ' + str(t) + ' --box_coo 300 330 200 230')
        os.system('python bl_evo.py --expid DA_REF DA_PSP DA_REF_TL500 DA_PSP_TL500 DA_PSPv2_TL500 --date 20160604000000 --time ' + str(t) + ' --box_coo 300 330 200 230 &')
