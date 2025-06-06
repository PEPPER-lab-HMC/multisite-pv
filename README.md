# multisite-pv
PV curve calculations for SRER, ONAQ, and CdM.

## Folders

data_clean
.csv's of ONAQ (sage_pv_curve.csv), CdM (juniper_pv_curve.csv), and SRER (pv_comb_202403308) 
shrubs. Each file contains weight and WP measurements to create PV curves.

source
Functions for creating .csv's for PV curves. pv - function2 generates the .csv's and plots
each point as they are collected, c.loc - function and corner - function fit models to the
PV curves.

scripts -> mod-1

scripts -> mod-2

## Scripts

01-test-function
Using pv - function2 to create SRER PV curves. Outputted files here are not included in repo.

02-clean-data-20231030
Manually generates PV curves for the SRER data (not using source functions). Incorporates
dry weight measurements. Creates pv_comb_20231030.csv (stored in data_clean).

02-clean-data-20240308
Manually generates PV curves for the SRER data (not using source functions). Incorporates
dry weight measurements. Creates pv_comb_20240308.csv (stored in data_clean).

03-clean-ONAQ
Manually generates PV curves for the ONAQ data. Incorporates dry weight measurements.
Creates sage_pv_curve.csv (stored in data_clean).

03-clean-CdM
Manually generates PV curves for the CdM data. Incorporates dry weight measurements. Creates
juniper_pv_curve.csv (stored in data_clean).

04-combine-sites
Joins together all clean data into one .csv. Creates multisite_pv_data.csv (stored in
data_clean). 

