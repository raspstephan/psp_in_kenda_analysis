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


## October 17 

ToDo: 
- Check DWD runs
    - noDA at 060512
    - Da_REF done! --> Copy! --> **SYNOP NOT IN FOF** --> why? is the original data not there?
- Continue with SYNOP/TEMP verification, note down missing data, make plan for runs!

### Why are SYNOP obs not in fof's?
- Synop files have wrong name: cdfin_synop_1.nc instead of cdfin_synop, cdfin_synop.2 --> Why?. There are two synop files in each tarred obs file. 
    - It's handled differently in run_cosmo, but maybe wrong??
    - Tried fixing it --> I think it worked!

Tomorrow:
- Check DWD runs, copy, start new!
- Finish ensemble airep
    - check which data is missing
- Start analysis precipitation verif--> Done

## October 18

Check which verifications are missing.
- SYNOP for det noDA_PSPv2 --> Done
- SYNOP for all forecast ensembles

Tomorrow:
- Redo ensemble fc precipitation with noDA --> Preprocessing done --> Running in screen
- Do fg precipitation analysis and put in paper.tex. --> Done
- Do Synop and air analysis for forecasts and put in paper.tex --> Done

## October 19

- Check DWD runs --> Synch last missing run. I should have all forecast data now except ensemble SYNOP.
    - Need to preprocess as well --> Done
- Bug in DA preprocessing, check tomorrow when nb is done! --> TMR

Tmr:
- include new ensemble plots --> Done
- Think about story for G! --> Done


## October 20

Debug:
- last time step for fg forecasts --> Found bug, rerunning fg preprocessing
- first time step for SYNOP (some of them...)
    - SYNOP pbs are not the same for all runs! (Probably one file missing or obs changed...)
    - But why is first bin off: Hardly any observations... Leave out? Yup!
    - Rerunning det fc's with new SYNOP --> Running

Next time:
- Check DWD runs and start next forecast

## October 23

Runs:
- DA_PSPv2 det didn't run --> restart

Tasks:
- See if fg bug is fixed! --> Yes. Recomputing all scores!
- Start paper!

Still to do:
- time axis label
- colors --> Doing that now!
- All days di prec
- Stamps
- Upscaled BS does not make any sense! (1)
- FSS for different scales (5 15) --> Done. put in paper --> Done
- FSS for larger threshold 5mm/h for 15 pts --> Running --> Done!

## October 24

- Check DWD runs --> Again, didn't work. Why? Wrong data directory --> Restart
- Keep writing!

(1) I think I made it somewhat analog to part of the FSS, recomputeing for 0.1, 1.0 and 15/31 and fg 15, probably should also do 5 mm, running for fg


Tmr: 
- Check DWD run, restart
- Include new BS results. --> Not sure I am doing this right...


## October 25

- DWD runs
    - not quite sure what is going on with DA_REF and DA_REF_new...
    - Started DA_REF_TL500
- Ensemble skill
    - Make a notebook