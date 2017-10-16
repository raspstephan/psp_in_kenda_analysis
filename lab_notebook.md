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

1. Try Brier Score for 0.1 and 1.0 mm/h --> Running in screen for upscale 1 and 31. Done
2. Run CRPS again with upscaling --> In screen. Done

Check out status of runs!

- Restarted noDA_PSPv2_ens with correct values --> Should be done end of the week!
- Restart DA_REF det --> Probably doesn't work. Check tomorrow.

ToDo:
- Scores (FSS) for 12UTC init. Is noDA still better than DA_PSP? --> In screen. Done
- Precip accum for 12 UTC --> Done
- FSS for 11 grid points
- Forecast T (+PS) verification --> Do that now!
- Analysis precipitation (ens!) verification
- Start coherent paper.tex --> started


## October 16

Runs:
- noDA still running
- REF restarted, should work this time!

Put results (BS, CRPS, FSS12) in paper.tex --> 
- FSS 12UTC: Is noDA still better than DA_PSP? No big difference. Does this indicate that using PSP in DA... not sure.
- BS/CRPS. There is some suggestions that at the time of triggering there is an icrease in spread, and therefore better reliability.
- 12UTC, initial spin up peak. Also for strong forcing, too much precip. contradicts 00UTC too little precip... 

Forecast SYNOP / TEMP verification
- apparently SYNOP is missing for ensemble... crap!
- Det Synop done, tomorrow TEMP/AIREP, also for ensemble.


