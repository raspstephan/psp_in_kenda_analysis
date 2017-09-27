# Lab Notebook for PSP in KENDA Analysis

## August 25

### Next Steps from meeting
- Separate analysis by period: strong and weak forcing. 
- Additional metrics
    - Precipitation histograms
    - SAL
    - CRPS

### Preload precipitation files
Preload precipitation files to make analysis a lot faster. 
Preloading done in ipynb, currently running for forecasts, det and ens.
Tomorrow I need to adjust the precipitation forecast analysis script to read the preloaded files. 

## August 26
Adjusted verif_fc_prec.py to use preloaded precipitation files. Added basic precipitation histograms and started with crps (data analysis working, but no plotting yet.)

Note: I mixed up convolved and not-convolved fields for FSS and RMSE for my previous plots. 

The preloading of the noDA_PSPv2 fields is not finished. This is because some files are missing (12UTC + 12h mem 16 26 May)

Continue tomorrow:
- Finish CRPS implementation and compute all fields. 
- Check DA_REF det run and see whether SYNOP data is there
- Start splitting basic metrics and putting plots in latex file!

## August 27