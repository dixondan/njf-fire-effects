Overview
--------

This repository is associated with the [open access paper](https://www.sciencedirect.com/science/article/pii/S0048969723034514, "Fire effects on eucalypt flowering")
 "Fire reduces eucalypt forest flowering phenology at the landscape-scale" in the journal *Science of the Total Environmnet* by Dan Dixon, Nik Callow, John Duncan, Sam Setterfield, and Natasha Paul (2023). 

<p align="center">
  <img src="graphabs.png" />
</p>

Included are the following data:
--------

  - Gridded flowering rasters for 2018, 2019 and 2020 found
  - Fire histry rasters including time since fire, severity and unique fire identifiers
  - CSV files to recreate estimates in Figure 3
 
And code to reproduce the analysis:
--------
   -  Sampling treatment and control forest pairs in Google Earth Engine (process/1_sample_from_ee.py)
   -  Organize those paris to run Levels and SFD models (process/2_organize.py)
   -  R code to run SFD/Levels model (process/3_make_fig.R)

Cite
--------


Questions
--------
Dan J. Dixon

Email: 1dandixon@gmail.com
Cooperative Research Centre for Honeybee Products: https://www.crchoneybeeproducts.com/
