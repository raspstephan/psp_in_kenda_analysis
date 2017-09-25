import os
datelist = ['20160526000000','20160527000000','20160528000000','20160529000000',
             '20160530000000','20160531000000','20160601000000','20160602000000',
             '20160603000000','20160604000000',]
datelist12 = []
for date in datelist[:-1]:
    os.system(
        'python verif_ana.py --expid DA_REF DA_PSPv2 DA_REF_TL500 DA_PSPv2_TL500 --date_ana_start ' +
        date + ' --date_ana_stop ' + date +
        ' --date_inc 24 --hint 24 --ana det')
    date12 = str(int(date) + 120000)
    datelist12.append(date12)
    os.system('python verif_ana.py --expid DA_REF DA_PSPv2 DA_REF_TL500 DA_PSPv2_TL500 --date_ana_start ' +
              date12 + ' --date_ana_stop ' + date12 +
              ' --date_inc 24 --hint 24 --ana det')

os.system('python verif_ana.py --expid DA_REF DA_PSPv2 DA_REF_TL500 DA_PSPv2_TL500 --date_start ' +
              datelist12[0] + ' --date_stop ' + datelist12[-1] +
              ' --date_inc 24 --hint 24 --ana det')