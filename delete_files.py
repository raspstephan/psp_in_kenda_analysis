"""
Script to delete COSMO output files to save space

"""
import os

tdir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/'

dates = ['20160529_00_','20160605_00_','20160606_00_12_']
exps = ['ref','std','const1','const3','nolowest','sig1','sig3','time10','time20']
starts = ['000000', '120000']
nens = 20

for d in dates:
    for e in exps:
        for s in starts:
            for i in range(1,nens+1):
                if d+e not in ['20160606_00_12_ref', '20160606_00_12_std']:
                    edir = (tdir + d + e + '/' + d[:8] + s + '/ens' + 
                            str(i).zfill(3))
                    print 'Deleting in dir:', edir
                    os.system('rm ' + edir + '/*_60')
            
