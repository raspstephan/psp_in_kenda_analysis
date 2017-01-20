import os
from cosmo_utils.helpers import yyyymmddhhmmss_strtotime, make_timelist, \
    yyyymmddhhmmss
from datetime import timedelta

tstart = yyyymmddhhmmss_strtotime('20160525120000')
tend = yyyymmddhhmmss_strtotime('20160609120000')
tint = timedelta(days = 1)
timelist = make_timelist(tstart, tend, tint)
ver_start = '600'
ver_end = '840'
var = 'T2M'
obs = 'SYNOP'

for t in timelist:
    print('python plot_stamps.py --date ' + yyyymmddhhmmss(t) + 
          ' --time 1 24 --plot prec_time --expid REF REF_TL500 PSP_TL500')
    os.system('python plot_stamps.py --date ' + yyyymmddhhmmss(t) + 
              ' --time 1 24 --plot prec_time --expid REF REF_TL500 PSP_TL500')
    print('python plot_stamps.py --date ' + yyyymmddhhmmss(t) + 
          ' --time 1 24 --plot prec_comp --expid REF REF_TL500 PSP_TL500')
    os.system('python plot_stamps.py --date ' + yyyymmddhhmmss(t) + 
              ' --time 1 24 --plot prec_comp --expid REF REF_TL500 PSP_TL500')
    
    # # T verif for each day individually
    # for e in ['REF', 'REF_TL500', 'PSP_TL500']:
    #     print('python verif_fof.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
    #         yyyymmddhhmmss(t) + ' --date_end ' + yyyymmddhhmmss(t) + ' --var ' + var +  
    #         ' --expid ' + e + ' --obs ' + obs)
    #     os.system('python verif_fof.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
    #         yyyymmddhhmmss(t) + ' --date_end ' + yyyymmddhhmmss(t) + ' --var ' + var +  
    #         ' --expid ' + e + ' --obs ' + obs)
    
    # # Comparison plot
    # allexp = 'REF REF_TL500 PSP_TL500'
    # print('python compare_verif.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
    #         yyyymmddhhmmss(t) + ' --date_end ' + yyyymmddhhmmss(t) + ' --var ' + var +  
    #         ' --expid ' + allexp + ' --obs ' + obs)
    # os.system('python compare_verif.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
    #         yyyymmddhhmmss(t) + ' --date_end ' + yyyymmddhhmmss(t) + ' --var ' + var +  
    #         ' --expid ' + allexp + ' --obs ' + obs)

# Average prec 
allexp = 'REF REF_TL500 PSP_TL500'
print('python average_prec.py --date_ini '+ yyyymmddhhmmss(tstart) + 
      ' --date_end '+ yyyymmddhhmmss(tend) +' --time 1 24 --expid ' + allexp)
os.system('python average_prec.py --date_ini '+ yyyymmddhhmmss(tstart) + 
      ' --date_end '+ yyyymmddhhmmss(tend) +' --time 1 24 --expid ' + allexp)



# # T verif for entire period
# for e in ['REF', 'REF_TL500', 'PSP_TL500']:
#     print('python verif_fof.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
#         yyyymmddhhmmss(tstart) + ' --date_end ' + yyyymmddhhmmss(tend) + ' --var ' + var +  
#         ' --expid ' + e + ' --obs ' + obs)
#     os.system('python verif_fof.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
#         yyyymmddhhmmss(tstart) + ' --date_end ' + yyyymmddhhmmss(tend) + ' --var ' + var +  
#         ' --expid ' + e + ' --obs ' + obs)
        
# allexp = 'REF REF_TL500 PSP_TL500'
# print('python compare_verif.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
#         yyyymmddhhmmss(tstart) + ' --date_end ' + yyyymmddhhmmss(tend) + ' --var ' + var +  
#         ' --expid ' + allexp + ' --obs ' + obs)
# os.system('python compare_verif.py --ver_start_min ' + ver_start + ' --ver_end_min ' + ver_end + ' --date_ini ' +
#         yyyymmddhhmmss(tstart) + ' --date_end ' + yyyymmddhhmmss(tend) + ' --var ' + var +  
#         ' --expid ' + allexp + ' --obs ' + obs)
