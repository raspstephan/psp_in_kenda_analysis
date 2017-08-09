# Scripts for analyzing my KENDA runs

The data is stored in the directory: `/project/meteo/w2w/A6/S.Rasp/kenda_psp_data/`

There might however be missing data.

4 experiments: `DA_REF`, `DA_PSPv2`, `DA_REF_TL500`, `DA_PSPv2_TL500`

PSPv2 has the following settings: blpert_sigma=2.5, blpert_const=2.0, blpert_fixedtime=600.

If TL500, blpert_const=0.6.

The DA cycling is run from 26 May 00UTC to 9 June 00UTC with ensemble analyses written every 12 hours. 

Deterministic forecasts are run every 12 hours with a forecast lead time of 24 hours. 

Ensemble forecasts with 20 members are run on 29 May 00UTC and from 4 June 00UTC to 6 June 12 UTC.

The DA runscripts are saved in `bacy_6000.04`, while the forecast runscripts are saved in `bacy_hreich`.