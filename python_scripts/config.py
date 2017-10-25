import os
import numpy as np

# General settings
if os.getcwd() == '/panfs/e/vol0/extsrasp/psp_in_kenda_analysis/python_scripts':
    plotdir = '/e/uwork/extsrasp/plots/'
    datadir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/'
    radardir = '/e/uwork/extsrasp/radolan/'
    savedir_base = '/e/uwork/extsrasp/save/'
elif 'home/s/S.Rasp/' in os.getcwd():
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
         'DA_REF': '#377eb8',
         'DA_REF_ens': '#377eb8',
         'DA_REF_TL500': '#4daf4a',
         'DA_REF_TL500_ens': '#4daf4a',
         'DA_PSPv2_TL500': '#ff7f00',
         'DA_PSPv2_TL500_ens': '#ff7f00',
         'DA_PSPv2': '#e41a1c',
         'DA_PSPv2_ens': '#e41a1c',
         'noDA_PSPv2': '#984ea3',
         'noDA_PSPv2_ens': '#984ea3',
         }

# cdict = {'radar': 'k',
#          'REF': 'navy',
#          'REF_TL500': 'darkgreen',
#          'PSP_TL500': 'orange',
#          'DA_REF': '#80b1d3',
#          'DA_REF_ens': '#80b1d3',
#          'DA_REF_TL500': '#8dd3c7',
#          'DA_REF_TL500_ens': '#8dd3c7',
#          'DA_PSPv2_TL500': '#ffffb3',
#          'DA_PSPv2_TL500_ens': '#ffffb3',
#          'DA_PSPv2': '#fb8072',
#          'DA_PSPv2_ens': '#fb8072',
#          'noDA_PSPv2': '#bebada',
#          'noDA_PSPv2_ens': '#bebada',
#          }

pw = 7.87

bin_edges = np.append(0, np.logspace(-1, 2, 10))
temp_bin_edges = np.arange(200, 1050, 50) * 100.
airep_bin_edges = np.arange(0, 10000 + 250., 250.)  # m

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
    'det_mean_cin': {
        'var': 'cin',
        'full_name': 'Deterministic mean CIN',
        'use_radar': False,
        'det_or_ens': 'det',
        'plot_type': 'line',
        'ylabel': 'CIN [J/kg]',
    },
    'det_median_prec': {
        'var': 'prec',
        'full_name': 'Deterministic median precipitation',
        'use_radar': True,
        'det_or_ens': 'det',
        'plot_type': 'line',
        'ylabel': 'Precipitation [mm/h]',
    },
    'det_median_cape': {
        'var': 'cape',
        'full_name': 'Deterministic median CAPE',
        'use_radar': False,
        'det_or_ens': 'det',
        'plot_type': 'line',
        'ylabel': 'CAPE [J/kg]',
    },
    'det_median_cin': {
        'var': 'cin',
        'full_name': 'Deterministic median CIN',
        'use_radar': False,
        'det_or_ens': 'det',
        'plot_type': 'line',
        'ylabel': 'CIN [J/kg]',
    },
    'det_sal': {
        'var': 'prec',
        'full_name': 'Deterministic SAL',
        'use_radar': True,
        'det_or_ens': 'det',
        'plot_type': 'sal',
        'ylabel': 'S (-) A(--) L(:) %s mm/h',
    },
    'det_fss': {
        'var': 'prec',
        'full_name': 'Deterministic FSS',
        'use_radar': True,
        'det_or_ens': 'det',
        'plot_type': 'line',
        'ylabel': 'FSS %s mm/h %s pts',
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
    'ens_rmse': {
        'var': 'prec',
        'full_name': 'RMSE of ensemble mean',
        'use_radar': True,
        'det_or_ens': 'ens',
        'plot_type': 'line',
        'ylabel': 'RMSE [mm/h]',
    },
    'ens_rmv': {
        'var': 'prec',
        'full_name': 'RMV of ensemble',
        'use_radar': False,
        'det_or_ens': 'ens',
        'plot_type': 'line',
        'ylabel': 'RMV [mm/h]',
    },
    'ens_bs': {
        'var': 'prec',
        'full_name': 'Ensemble Brier Score',
        'use_radar': True,
        'det_or_ens': 'ens',
        'plot_type': 'line',
        'ylabel': 'BS %s mm/h %s pts',
    },
}