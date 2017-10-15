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

ToDo:
- finish CRPS
    - Plot
    - Run for all days
- Fix histograms
    - Rescale
- Look at CAPE/CIN

I did not actually manage to do any of these things but I created a good framework for my analysis.

Another ToDo is the missing values for the noDA run!


## August 28

ToDo:
- Radar plots

CRPS for 00UTC running in screen.

I think I need some sort of upscaling. Add as additional parameter.


## October 15

Goal: Find out about ensemble spread in precipitation!

1. Try Brier Score for 0.1 and 1.0 mm/h --> Running in screen for upscale 1 and 31