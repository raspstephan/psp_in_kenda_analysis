import os
datelist = ['20160526000000','20160527000000','20160528000000','20160529000000',
             '20160530000000','20160531000000','20160601000000','20160602000000',
             '20160603000000','20160604000000','20160605000000','20160606000000',
             '20160607000000','20160608000000','20160609000000']
datelist12 = []
for date in datelist[:-1]:
    date12 = str(int(date) + 120000)
    datelist12.append(date12)
    os.system('python verif_fc_prec.py --expid DA_REF DA_PSPv2 DA_REF_TL500 DA_PSPv2_TL500 --date_start ' +
              date12 + ' --date_stop ' + date12 +
              ' --date_inc 24 --hint 24 --ana det')

os.system('python verif_fc_prec.py --expid DA_REF DA_PSPv2 DA_REF_TL500 DA_PSPv2_TL500 --date_start ' +
              datelist12[0] + ' --date_stop ' + datelist12[-1] +
              ' --date_inc 24 --hint 24 --ana det')