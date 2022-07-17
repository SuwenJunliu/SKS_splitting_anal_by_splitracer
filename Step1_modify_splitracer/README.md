# Step.1 Modify the SplitRacer

## What is the SplitRacer

`SplitRacer: MATLAB Code and GUI for Semiautomated Analysis and Interpretation of Teleseismic Shear‐Wave Splitting 
https://doi.org/10.1785/0220160191`

SplitRacer is a tool for analysis the Shear‐Wave, especially the SKS phase. The Chinese version of introduction is available at `https://doi.org/10.16738/j.cnki.issn.1003-3238.201802006, SplitRacer:用于远震剪切波分裂半自动分析和解释的MATLAB代码和图形用户界面`

## Why we choose and modify the SplitRacer 

There are many useful tools for Shear‐Wave Splitting such as SplitLab (https://github.com/IPGP/splitlab). The reasons why I choose SplitRacer are below:

1. SplitRacer is based on Matlab, which is easy for coding and debugging.
2. SplitLab is running on old versions of MATLAB, especially compared to SplitRacer. I met a lot problems while debugging the SplitLab.
3. In 2021, Frederik Link et al. update the SplitRacer, and add the auto processing function, which makes it possible for analyzing the dense array. 

`An automatized XKS-splitting procedure for large data sets: Extension package for SplitRacer and application to the USArray https://doi.org/10.1016/j.cageo.2021.104961`

Unfortunately, the format of seismic data using by SplitRacer is miniSeed. In my case, our data format is SAC. So I modify the serveral matlab code file. There are many tiny modifications to make it suitable for our real stituation.

## What is the difference?

Here are the list for code file

For non-auto code:
- pre-processing/get_channel_info.m
- pre-processing/sub_window_new.m
- pre-processing/prep_data.m
- pre-processing/re_sample.m
- pre-processing/snr.m

For auto code:
- auto/qc_phases_auto.m
- auto/extract_adv_auto_sac.m
- auto/prep_data_auto.m
- auto/get_channel_info_auto.m

We also add some code file

- pre-processing/extract_adv_sac.m
- pre-processing/extract_sac.m
- auto/extract_adv_auto_sac.m

## How to use those code

Just subsistude those code file.


## Why I modify those code?

### pre-processing/get_channel_info.m

Sometimes our sacfile will contain some space in description, so we need to add the separator `|`. Same reason for file `auto/get_channel_info_auto.m`

### pre-processing/sub_window_new.m

Add the earthquake info (origin time) to the figure, which make it easier for selecting the event by output figure. 

### pre-processing/prep_data.m

Main program for pre-processing. Substitude the MiniSeed reading function (extract_adv/extract) with SAC reading function. Same reason for file `auto/prep_data_auto.m`


### pre-processing/re_sample.m

For "view waveforms" processing, it will show the one-hour length of each waveform. However, it doesn't last for so long in our data. The waveform will contain a linear incresement (or decline) by interpolation. So we change the resample length to 30 min.

### pre-processing/snr.m

In Splitracer, the signal noise ratio is calculated as follows. It will read the waveform with the length of 100s, then the center of the waveform will be considered as the arrival time of SKS phase (or other splitting phase). The signal window is [center, center+20], and the noise window is [center-15, center]. Sometimes the SKS phase will appear earlier then the therotical arrival time, so the SNR will significant small. To avoid this issue, we set the "center" a little bit "earlier" (around 5s).

### auto/qc_phases_auto.m

Adjust the parameters of STA/LTA for low SNR data.

### auto/extract_adv_auto_sac.m

The sacfile reading component, same reason for file `pre-processing/extract_adv_sac.m` `pre-processing/extract_sac.m` `auto/extract_adv_auto_sac.m`

