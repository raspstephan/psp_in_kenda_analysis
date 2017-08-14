import os

# General settings
if os.getcwd() == '/panfs/e/vol0/extsrasp/dwd_scripts':
    plotdir = '/e/uwork/extsrasp/plots/'
    datadir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/'
    radardir = '/e/uwork/extsrasp/radolan/'
    savedir_base = '/e/uwork/extsrasp/save/'
elif os.getcwd() == '/home/s/S.Rasp/repositories/dwd_scripts':
    datadir = '/project/meteo/w2w/A6/S.Rasp/kenda_psp_data/data_forecast/'
    datadir_cosmo2 = '/project/meteo/cosmo/stephan.rasp/dwd_data/data/'
    datadir_raid2 = '/project/meteo/data/raid_linux/stephan.rasp/dwd_data/data/'

    radardir = '/project/meteo/w2w/A6/radolan/netcdf_cosmo_de/'
    plotdir = '/home/s/S.Rasp/repositories/dwd_scripts/figures/'
    savedir_base = '/project/meteo/w2w/A6/S.Rasp/kenda_psp_data/save/'
    feeddir = '/project/meteo/w2w/A6/S.Rasp/kenda_psp_data/feedback/'
else:
    raise Exception('Working directory not recognized:' + os.getcwd())

# Config for experiment
cdict = {'radar': 'k',
         'REF': 'navy',
         'REF_2JUN_ens': 'navy',
         'REF_TL500': 'darkgreen',
         'PSP_TL500': 'orange',
         'DA_REF': 'navy',
         'DA_REF_2JUN': 'navy',
         'DA_REF_ens': 'navy',
         'DA_REF_TL500': 'cyan',
         'DA_REF_TL500_2JUN': 'cyan',
         'DA_REF_TL500_2JUN_ens': 'cyan',
         'DA_REF_TL500_ens': 'cyan',
         'DA_PSP_TL500': 'red',
         'DA_PSPv2_TL500': 'fuchsia',
         'DA_PSPv2_TL500_2JUN': 'fuchsia',
         'DA_PSPv2_TL500_2JUN_ens': 'fuchsia',
         'DA_PSPv2_TL500_ens': 'fuchsia',
         'DA_PSPv2': 'maroon',
         'DA_PSPv2_2JUN': 'maroon',
         'DA_PSPv2_ens': 'maroon',
         'DA_PSPv2_2JUN_ens': 'maroon',
         }
