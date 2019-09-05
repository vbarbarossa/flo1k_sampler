# flo1k_sampler
R scripts to get streamflow time series of FLO1K at user-specified locations (inclusive of reallocation routine on the FLO1K river network)

If you end up using this script, please ackowledge the FLO1K paper:
Barbarossa et al. 2018. FLO1K, global maps of mean, maximum and minimum annual streamflow at 1 km resolution from 1960 through 2015. https://www.nature.com/articles/sdata201852

To run the main script a basin knowledge of R is required. The scripts have been tested on R version 3.4.2.

Usage:
- make sure the libraries "raster", "sf" and "foreach" are installed
- modify the user settings in the main script "sample_flo1k.R"
- details on input data requirements are also specified in the same "sample_FLO1K.R" script
- output is stored in the "out" folder
