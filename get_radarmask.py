from cosmo_utils.pyncdf import getfobj_ncdf_timeseries
radarpref = '/project/meteo/w2w/A6/radolan/netcdf_cosmo_de/raa01-rw_10000-'
radarsufx = '-dwd---bin.nc'
from datetime import timedelta
from cosmo_utils.helpers import yyyymmddhh_strtotime
tstart = '2016052601'
tend = '2016061000'

dtradar = timedelta(minutes = 10)

radarts = getfobj_ncdf_timeseries(radarpref, yyyymmddhh_strtotime(tstart)-dtradar, 
                                  yyyymmddhh_strtotime(tend)-dtradar, 
                                  timedelta(minutes=60), 
                                  reftime=yyyymmddhh_strtotime(tstart), ncdffn_sufx=radarsufx, 
                                  fieldn = 'pr', abs_datestr='yymmddhhmm',dwdradar = True)

from cosmo_utils.diag import get_totmask
radarmask = get_totmask(radarts)

import matplotlib.pyplot as plt
plt.imshow(radarmask)
plt.savefig('./radarmask')

import numpy as np
np.save('./radar_tot_mask.npy',radarmask)
history

