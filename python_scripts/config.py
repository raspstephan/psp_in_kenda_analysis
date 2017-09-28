import os
import numpy as np

# General settings
if os.getcwd() == '/panfs/e/vol0/extsrasp/psp_in_kenda_analysis/python_scripts':
    plotdir = '/e/uwork/extsrasp/plots/'
    datadir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/'
    radardir = '/e/uwork/extsrasp/radolan/'
    savedir_base = '/e/uwork/extsrasp/save/'
elif os.getcwd() == '/home/s/S.Rasp/repositories/psp_in_kenda_analysis/python_scripts':
    datadir = '/project/meteo/w2w/A6/S.Rasp/kenda_psp_data/data_forecast/'
    datadir_da = '/project/meteo/w2w/A6/S.Rasp/kenda_psp_data/data/'
    datadir_cosmo2 = '/project/meteo/cosmo/stephan.rasp/dwd_data/data/'
    datadir_raid2 = '/project/meteo/data/raid_linux/stephan.rasp/dwd_data/data/'

    radardir = '/project/meteo/w2w/A6/radolan/netcdf_cosmo_de/'
    plotdir = '/home/s/S.Rasp/repositories/psp_in_kenda_analysis/figures/'
    savedir_base = '/project/meteo/w2w/A6/S.Rasp/kenda_psp_data/save/'
    feeddir = '/project/meteo/w2w/A6/S.Rasp/kenda_psp_data/feedback/'
else:
    raise Exception('Working directory not recognized:' + os.getcwd())

# Config for experiment
cdict = {'radar': 'k',
         'REF': 'navy',
         'REF_TL500': 'darkgreen',
         'PSP_TL500': 'orange',
         'DA_REF': 'navy',
         'DA_REF_ens': 'navy',
         'DA_REF_TL500': 'cyan',
         'DA_REF_TL500_ens': 'cyan',
         'DA_PSP_TL500': 'red',
         'DA_PSPv2_TL500': 'fuchsia',
         'DA_PSPv2_TL500_ens': 'fuchsia',
         'DA_PSPv2': 'red',
         'DA_PSPv2_ens': 'red',
         'noDA_PSPv2': 'maroon',
         'noDA_PSPv2_ens': 'maroon',
         }

pw = 7.87

bin_edges = np.append(0, np.logspace(-1, 2, 10))


metric_dict = {
    'det_rmse': {
        'var': 'prec',
        'full_name': 'Deterministic RMSE',
        'use_radar': True,
        'det_or_ens': 'det',
        'plot_type': 'line',
        'ylabel': 'RMSE [mm/h]',
    },
    'det_mean_prec': {
        'var': 'prec',
        'full_name': 'Deterministic mean precipitation',
        'use_radar': True,
        'det_or_ens': 'det',
        'plot_type': 'line',
        'ylabel': 'Precipitation [mm/h]',
    },
    'det_mean_cape': {
        'var': 'cape',
        'full_name': 'Deterministic mean CAPE',
        'use_radar': False,
        'det_or_ens': 'det',
        'plot_type': 'line',
        'ylabel': 'CAPE [J/kg]',
    },
    'det_fss': {
        'var': 'prec',
        'full_name': 'Deterministic FSS',
        'use_radar': False,
        'det_or_ens': 'det',
        'plot_type': 'line',
        'ylabel': 'FSS',
    },
    'det_prec_hist': {
        'var': 'prec',
        'full_name': 'Deterministic precipitation histogram',
        'use_radar': True,
        'det_or_ens': 'det',
        'plot_type': 'hist',
        'ylabel': 'Frequency',
    },
    'ens_crps': {
        'var': 'prec',
        'full_name': 'Ensemble CRPS',
        'use_radar': True,
        'det_or_ens': 'ens',
        'plot_type': 'line',
        'ylabel': 'CRPS',
    },
}