import os
datelist = ['20160526000000','20160527000000','20160528000000','20160529000000',
             '20160530000000','20160531000000','20160601000000','20160602000000',
             '20160603000000','20160604000000','20160605000000','20160606000000',
             '20160607000000','20160608000000','20160609000000']
for date in datelist:
    os.system('python plot_stamps.py --expid DA_REF DA_PSPv2 DA_REF_TL500 DA_PSPv2_TL500 --date ' +
              date + ' --time 1 24 --plot prec_comp')

