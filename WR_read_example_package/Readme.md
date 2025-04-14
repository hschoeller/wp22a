# WR_read_examples.ipynb
## Jupyter Notebook demonstrating reading WR data (IWR and LC) and doing simple analysis
### Notebook description

This notebook intends to simplify the use of the 7 year-round weather regimes for the North-Atlantic European region introduced by Grams et al. (2017) [doi:10.1038/nclimate3338](https://doi.org/10.1038/nclimate3338) and updated on ERA5 as described by Hauser et al. (2024) [doi:10.5194/wcd-5-633-2024](https://doi.org/10.5194/wcd-5-633-2024). Christian thanks Sera and Dominik for help with coding this ipynb.

**Part 1** uses fct_wrera_db (V Nov 2020) and fct_wrlcera (V Sep 2019) of Dominik Büeler to read the WR data as explained below. The former reads the simpler file containing IWR time series and the categorical maxIWR and LCattr attributions. The latter generates LC objects.

**Part 2** generates an nice time series plots of IWR, with the active life cycles marked, and a marker for the life cycle attribution.

**Part 3** computes a frequency climatology in a given period and plots the frequency absolute as well as as an anomaly.


author contact: christian.grams@gmx.de

based on templates by: seraphine.hauser@gmx.de and dominik.bueeler@meteoswiss.ch

date: 3 February 2025



# USE OF THIS DATA AND NOTEBOOK:

Please note that at the moment I still try to keep track who is using which version of the data, because I am working on detailed documentation paper. The paper will describe the methodology as well as key characteristics and trends in WR life cycles. With the paper I will make the data freely available. In the mean time I want to be able to update users individually, to ensure they use the correct data. So far the methodology is described in the #LSDPatKIT teams's various papers and initially in Grams et al. 2017  [doi:10.1038/nclimate3338](https://doi.org/10.1038/nclimate3338) (methodology section). Please cite this paper in the meantime and Hauser et al. 2024 [doi:10.5194/wcd-5-633-2024](https://doi.org/10.5194/wcd-5-633-2024) for the slight modifications with the update on ERA5. Until the definition paper is out, please let me know before you share the data beyond your project. 




# Readme for ERA5 year-round 7 regime life cycles. 
## Configuration 

- 1979-2019 reference climatology for Z500 anomaly computation and IWR scaling 
- 10d low-pass filter 
- period 19790111_00 – 20191231_21 used for EOF clustering
- 7 leading EOFs explain 74,4% of variance

### Weather regime indices in correct order

**Index: **           1 6 7 2 4 5 3 0 

**Abbreviation: **    AT ZO ScTr AR EuBL ScBL GL no  

### Data directory

``DATADIR=./ec.era5/eof/update_1950/latwgt/eof_AUTO_Z500_N161_-80E40E30N90N_year/norminput/``

### 90d running mean frequencies and normalization weights

`cd <DATADIR>`

``cat Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_proj_local_nrmfct.txt``


```
    mean and std proj. (Michel and Riviere, 2011) 19790111_00 to 20191231_18
    AT ZO ScTr AR EuBL ScBL GL 
    mean  -0.0023841809 -0.0078344550 -0.0015250760  0.0021150061  0.0042730025  0.0089835953  0.0049633454
    stdv   0.1711293906  0.2818541229  0.2004271448  0.2492269725  0.2292326987  0.2332117558  0.3415428102
```

## Changes to ERA-Interim definition

- normalization weights using the 'lat-weighted spatial mean (in EOF domain) of the grid-point based 30d temporal stddev of Z0500'. This improves summer regime identification.
- EOF-clustering performed for 6h data (insensitive to 6h, 12, 48h time interval) 
- IWR (projection) computed for 3h data
- Z500 normweight for 3,9,15,21 UTC uses normweight of 0,6,12,18 UTC
- Standard abbreviations in filenames (no longer use abbreviations ZOWE, ZOEA, BL)




## Weather regime attribution and IWR files

### Regime Attribution file 

This file is read in by calling `fct_wrera_db.py`.

**ERA5: **`$DATADIR/Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local.txt`

```
WR index after Michel and Riviere (2011) of filtered data N161 low pass > 10days, Z0@500, normed: intersection times, OVERLAP LCs allowed
--------------------------------------------------------
Cluster Class Index 1-4:  AT AR GL EuBL ScBL ZO ScTr
--------------------------------------------------------
time in h since 19790101_00| YYYYMMDD_HH| cluster class index| max WR index| lifecycle WR index 
--------------------------------------------------------

     240 19790111_00    7    7    0
     243 19790111_03    0    7    0
     246 19790111_06    7    7    0
     249 19790111_09    0    7    0
     252 19790111_12    7    7    0
     255 19790111_15    0    7    0
```

The file contains 5 columns, the 5th is most relevant:

1. hour since 19790101_00
2. date in yyyymmdd_hh
3. eof attribution -> timestep contributing to cluster XY based on EOFs (only EOF clustering period 11.1.1979-31.12.2019 only for 0 6 12 18 UTC times, other times are "0")
4. max WR index (Michel and Rivière, 2011).  -> based on weather regime projection in physical space
5. LC attribution based on WR index (Michel and Rivière, 2011). -> THIS IS THE LIFECYCLE attribution INCLUDING NO REGIME. USE THIS.

***Regime indices and regime order: *** 

The indices 0-7 refer to the following regimes: NAME (ABBREVIATION IN FILES for ERA5)

0-7: no AT AR GL EuBL ScBL ZO ScTr

0=no regime [only 5th column],
1=Atlantic Trough (AT),
2=Atlantic Ridge (AR),
3=Greenland Blocking (GL),
4=European Blocking (EUBL),
5=Scandinavian Blocking (ScBLL),
6=Zonal regime(ZO=NAO+),
7=Scandinavian Trough(ScTr),

As in this Notebook, please write your code so that the regimes appear sorted in the following way (this will group related regimes). **Please use this order when plotting panels / grouping regimes as well.** Reason: related regimes are grouped next to each other, and cyclonic / blocked are grouped together. In all our papers we use this scheme.

Index re-order:             1 6 7 2 4 5 3 0

Abbreviation:                AT ZO ScTr AR EuBL ScBL GL no

If possible also use the following color codes to colors regimes in plots. Reason again is that the mixed colors reflect relations between regimes, and this is the color scheme we use in all our studies (better intercomparison of studies). An RGB Table and python matplotlib names used for Bueeler et al. 2021 (https://doi.org/10.1002/qj.4178) and Osman et al. 2023 (https://doi.org/10.1002/qj.4512) are listed below and used in this notebook.
You might want to read this file in python.
`wr_metadata_rgb_new_colors.txt`

```
    Long name		        Name	Color		RGB
    Atlantic trough		    AT	    indigo		(75,   0,130)
    Zonal	 		        ZO	    red		    (255,  0,  0)
    Scandinavian trough	    ScTr	darkorange	(255,140,  0)
    Atlantic ridge		    AR	    gold		(255,215,  0)
    European blocking	    EuBL	yellowgreen	(154,205, 50)
    Scandinavian blocking	ScBL	darkgreen	(  0,100,  0)
    Greenland blocking	    GL	    blue		(  0,  0,255)
    No regime 		        no	    grey		(128,128,128) 
```

* * *

### Regime projection vector 

**ERA5: **`$DATADIR/Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_proj_local.txt`

```
WR index after Michel and Riviere (2011) of filtered data N161 low pass > 10days, Z0@500 normed: intersection times
--------------------------------------------------------
Cluster Class Index 1-4:  AT AR GL EuBL ScBL ZO ScTr
--------------------------------------------------------
time in h since 19790101_00, YYYYMMDD_HH, cluster class index, max WR index, WR index for regimes: AT ZO ScTr AR EuBL ScBL GL 
--------------------------------------------------------

     240 19790111_00    7    7   0.86587346   0.75074303   1.35784054  -0.45948577  -1.18579698  -1.07059395  -0.34373909
     243 19790111_03    0    7   0.85441613   0.72267860   1.31580353  -0.44945070  -1.16703987  -1.03862977  -0.32952571
     246 19790111_06    7    7   0.83927369   0.69352823   1.27102554  -0.43822104  -1.14390016  -1.00381625  -0.31613618
     249 19790111_09    0    7   0.82096398   0.66378951   1.22433603  -0.42610571  -1.11709785  -0.96681082  -0.30378893
     252 19790111_12    7    7   0.79906243   0.63312817   1.17505550  -0.41290590  -1.08599937  -0.92708772  -0.29231820
     255 19790111_15    0    7   0.77413946   0.60201460   1.12398231  -0.39892274  -1.05132735  -0.88528168  -0.28191039
     
```

This file contains the projection of instantaneous normalized Z500 anomalies in the WR patterns (7 dimensional WR projection vector) following Michel and Rivière 2011. It better describes a current flow situation in terms of ressemblance to each of the 7 cluster mean EOF patterns. It is used to objectively identify the regime life cycles.

The files contains 11 columns

1. hour since 19790101_00
2. date in yyyymmdd_hh
3. eof attribution -> timestep contributing to cluster XY based on EOFs (only EOF clustering period 11.1.1979-31.12.2019 only for 0 6 12 18 UTC times, other times are "0")
4. max WR index (Michel and Rivière, 2011).  -> based on weather regime projection in physical space
5. to 11. WR index following Michel and Rivière, 2011 **correctly ordered**: AT ZO ScTr AR EuBL ScBL GL

Part 2 generates an example time series plot for 2021 which shows each WR time series and the identified active life cycles (in bold). On the bottom the "dominant active life cycle" is marked. The latter corresponds to the unambigous "LC attribution" (5th column in `Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local.txt`) and is the active LC with maximum projection at that time (if two or more LC coexist). 



* * *
## WR lifecycle files

These files are read in by calling `fct_wrlcera.py`.

The files `$DATADIR/Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local_cluster[0..7].txt` contain the objective life cycle definition.

The header for regimes 1-7 look as follows, for no regime (index 0) it is reduced and defined differently:

```
WR lifecycle based on projections of filtered data N161 low pass > 10days, Z0@500 normed: intersection times, OVERLAP LCs allowed
Projection P after Michel and Riviere (2011): cos-wgt ano*classmean and WR index: I=(P(t)-avg(P))/std(P(t)
--------------------------------------------------------
 LIFECYCLE for members in Cluster 1
AT
 clsfd EOF : 7958/59860 (13.2944%)
 total mxI : 18275/126376 (14.4608%)
--------------------------------------------------------
number    onset     sat start      mx       sat end     decay        dcfr  dcto dctoID dctoDATE    onfr onto onfrID onfromDATE trfr trfrID trfromDATE trto trtoID trtoDATE 
--------------------------------------------------------

     0 19790518_06 19790525_15 19790526_18 19790527_18 19790530_03     AT   EuBL    1 19790531_06  none    AT -999 19790518_03    AT    0 19790519_03    AT    0 19790529_21
     1 19791209_06 19791211_12 19791212_21 19791213_21 19791215_18     AT   none -999 19791215_21    ZO    ZO    1 19791209_03    ZO    1 19791209_12  none -999 19791215_21
     2 19801020_06 19801020_06 19801023_21 19801027_03 19801027_03     AT   none -999 19801027_06    GL    GL   11 19801020_03    GL   11 19801021_21  none -999 19801027_06
```

### key parameters defining the life cycle

- **number**: exclusive ID of the lifecycle
- **onset**: onset date
- *sat start*: begin of saturation stage date (not used)
- **mx**: maximum stage date
- *sat end*: end of saturation stage date (not used)
- **decay**: decay date

### Transitions within a time window and based on life cycle (use these)

- *dcfr*: regime type of active dominant life cycle at *decay* (not used)
- **dcto**: within 4 days after the decay (dc, dc+96h), type of first active dominant lifecycle or none
- **dctoID**: within 4 days after the decay(dc, dc+96h), ID of first active dominant lifecycle or -999 for none
- **dctoDATE**: within 4 days after the decay(dc, dc+96h), date when the other active dominant lifecycle is identified for the first time. For none this is *dt* after the *decay*
- **onfr**: up to 4 days prior to onset (on-96h,on), type of first (backward looking) active dominant lifecycle or none
- *onto*: regime type of active dominant lifecycle at *onset* (not used)
- **onfrID**: up to 4 days prior to onset (on-96h,on), ID of first active dominant lifecycle or -999 for none
- **onfrDATE**: up to 4 days prior to onset (on-96h,on), date when the other active dominant lifecycle is identified for the first time. For none this is *dt* before the *onset*

### Immediate transitions based on dominant active life cycle (not used)

- *trfr*: regime type of active dominant life cycle when the current life cycle becomes dominant for the first time. It can be the life cycle itself (in the case that max projection is reached for the first time but another projection was larger without contributing to a LC (not persistent enough)).
- *trfrID*: ID of the *trfr* LC
- *trfrDATE*: date when the *trfr* life cycle was dominant for the last time (this is *-dt*h (one time step) before the considered LC becomes dominant for the first time)
- *trto*: regime type of active dominant life cycle when the current life cycle does no longer have the strongest projection for the first time. It can be the life cycle itself (in the case that another regime index is higher but this regime does not become an active life cycles (not persistent enough)).
- *trtoID*: ID of the *trto* LC
- *trtoDATE*: date when the *trto* life cycle is dominant for the first time


### For the no regime two simplified files are contained:

`Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local_cluster0_times.txt`: simply all time steps with no regime

`Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local_cluster0.txt`: a simplified life cycle file indicating the begin and end of a no regime period.

