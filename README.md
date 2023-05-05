# ASA2023-OAEs
Analysis code for OAE data for ASA 2023 project

## Initial analysis 
Input the subject ID and ear to be analyzed in either analysis script in the form: 
`SFanalysis('Q421', 'R')' or 'DPanalysis('S343', 'L')' 
If the subject ID starts with an "S" as is in the convention in the SNAPlab, you'll also want to be sure you have the FPL calibration file associated with that subject and ear. Values with be converted to EPL. 
Figures will be plotted of the response and the variables will be saved in a structure called 'res'. 

'pullData.m' is designed to run the analysis script on many subjects and ears at once. If a data file is not found, it will be skipped. 

'CompareSFandDP.m' has code for averaging responses across half-octave frequency bands and looking at the data for each exposure groups. DP and SF amplitudes are also plotted against each other using this code. Resulting data is stored in a table which can be read into R for LME modeling/statistics. 

## Data Collection
Code for data collection is stored in two separate repositories DPOAEswept and SFOAEswept. 
