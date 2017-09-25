import datetime as dt
import os

dstart = '20160526000000'
dstop = '20160609000000'
td_int = dt.timedelta(hours=12)

dformat = '%Y%m%d%H%M%S'

timelist = [dstart]
dobj = dt.datetime.strptime(dstart, dformat)
dstr = dstart
while dstr != dstop:
    dobj += td_int
    dstr = dt.datetime.strftime(dobj, dformat)
    timelist.append(dstr)

base_dir = '/project/meteo/w2w/A6/S.Rasp/kenda_psp_data/data/DA_REF/'
copy_dir = 'det_anai'

for t in timelist:
    s = 'cp %s%s/laf%s.det %s%s/%s/' % (base_dir, t, t, base_dir, copy_dir, t)
    print(s)
    os.system('mkdir %s%s/%s' % (base_dir, copy_dir, t))
    os.system(s)
