# Data sources

**Strains**
HKK796 : slo-1::GFP
NM1968 : slo-1(js379)
HKK1165: slo-1::GFP panneuronal degredation of SLO-1, rgef-1 promoter
VG902 : slo-1::GFP ciliated neuron degradation of SLO-1, osm-6 promoter
VG903 : slo-1::GFP mechanosensory neuron degradation of SLO-1, mec-18 promoter
N2 : wildtype

## Procedure - TWR/TAR/Initial
1. use python to get sets of MWT paths : get_mwtpath_gfp_control.ipynb
2. import paths using matlab and run through TWR, TAR, Initial : run_std_twr_tar_initial.m
3. get raw data from TWR and TAR : get_Dance_rawdata2csv.m
4. calculate TAR last 3 taps : TAR_last3taps.py


## Procedure - integral (depgreciated)
1. convert .mat data into csv data (rawdata.csv and mwtdb.csv) : `convert_shanespark_raw2csv.m`
2. combine all raw data in one csv (rawdata_all.csv, mwtdb_all.csv) : `rawdata_csv.py`
3. do analysis on integral: integral_diff.ipynb

