import os
# General settings
if os.getcwd() == '/panfs/e/vol0/extsrasp/dwd_scripts':
    plotdir = '/e/uwork/extsrasp/plots/'
    datadir = '/e/uwork/extsrasp/cosmo_letkf/data_forecast/'
    radardir = '/e/uwork/extsrasp/radolan/'
    savedir_base = '/e/uwork/extsrasp/save/'
elif os.getcwd() == '/home/s/S.Rasp/repositories/dwd_scripts':
    datadir_cosmo = '/project/meteo/cosmo/stephan.rasp/dwd_data/data_forecast/'
    datadir_raid = '/project/meteo/data/raid_linux/stephan.rasp/dwd_data/data_forecast/'
    datadir_cosmo2 = '/project/meteo/cosmo/stephan.rasp/dwd_data/data/'
    datadir_raid2 = '/project/meteo/data/raid_linux/stephan.rasp/dwd_data/data/'
    datadirdict = {'REF': datadir_cosmo,
                   'REF_2JUN_ens': datadir_raid,
                   'DA_REF': datadir_cosmo,
                   'DA_REF_2JUN': datadir_cosmo,
                   'DA_PSPv2': datadir_raid,
                   'DA_PSPv2_2JUN': datadir_raid,
                   'DA_PSPv2_ens': datadir_cosmo,
                   'DA_PSPv2_2JUN_ens': datadir_cosmo,
                   'DA_PSPv2_TL500': datadir_cosmo,
                   'DA_PSPv2_TL500_2JUN': datadir_cosmo,
                   'DA_REF_TL500': datadir_raid,
                   'DA_REF_TL500_2JUN': datadir_cosmo,
                   'DA_REF_TL500_2JUN_ens': datadir_raid,
                   'DA_REF_ens': datadir_raid,
                   'DA_PSPv2_ens': datadir_raid,
                   'DA_PSPv2_TL500_ens': datadir_raid,
                   'DA_PSPv2_TL500_2JUN_ens': datadir_cosmo,
                   'DA_REF_TL500_ens': datadir_raid,
                   }
    datadirdict2 = {'DA_REF': datadir_cosmo2,
                   'DA_PSPv2': datadir_raid2,
                   'DA_PSPv2_TL500': datadir_cosmo2,
                   'DA_REF_TL500': datadir_raid2,
                   }
    radardir = '/project/meteo/w2w/A6/radolan/netcdf_cosmo_de/'
    plotdir = '/home/s/S.Rasp/dwd_plots/plots/'
    savedir_base = '/project/meteo/cosmo/stephan.rasp/dwd_data/save/'
    feeddir = '/project/meteo/cosmo/stephan.rasp/dwd_data/feedback/'
else: 
    raise Exception('Working directory not recognized:' + os.getcwd())


# Config for experiment
cdict = {'radar':'k',
             'REF':'navy',
             'REF_2JUN_ens': 'navy',
             'REF_TL500':'darkgreen',
             'PSP_TL500':'orange',
            'DA_REF':'navy',
            'DA_REF_2JUN':'navy',
            'DA_REF_ens':'navy',
            'DA_REF_TL500':'cyan',
            'DA_REF_TL500_2JUN':'cyan',
            'DA_REF_TL500_2JUN_ens':'cyan',
            'DA_REF_TL500_ens':'cyan',
            'DA_PSP_TL500':'red',
            'DA_PSPv2_TL500':'fuchsia',
            'DA_PSPv2_TL500_2JUN':'fuchsia',
            'DA_PSPv2_TL500_2JUN_ens':'fuchsia',
            'DA_PSPv2_TL500_ens':'fuchsia',
            'DA_PSPv2':'maroon',
            'DA_PSPv2_2JUN':'maroon',
            'DA_PSPv2_ens':'maroon',
            'DA_PSPv2_2JUN_ens':'maroon',
            'DA_PSPv2_ens':'maroon',
             }

def strip_expid(expid):
    return expid.replace('DA_', '').replace('_ens', '').replace('v2', '').replace('_2JUN', '')
