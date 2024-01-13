
## *In preparation* R-code for survival analysis using combined dead and live encounter data of buzzards *Buteo buteo*

<img src="plots/Fig_1.png" style="display: block; margin: auto;" />

*Due to a 100 MB size limit, the data file `models.barker.adj.RDS` can
not be uploaded to GitHub in the moment*

-   **01_R\_code.Rmd** contains the code to reproduce all analyses and
    figures presented in the manuscript.
-   **02_R\_code.Rmd** Supplementary analysis of the
    capture-mark-resighting data using the burnham models.
-   **03_R\_code.Rmd** Supplementary analysis of the
    capture-mark-resighting data using the CJS models.

## **`data`**

This folder contains raw data and cached results of the
capture-mark-resighting analysis. The most important raw data files are
listed below.

-   **ringing_data.csv**: Individual-specific data of 3463 common
    buzzards. Note, a subset of 2723 individuals is used in
    capture-mark-resighting models.

-   **resighting_data.csv**: Data on 2940 resighting events, of which
    2793 correspond to individuals represented in
    capture-mark-resighting models.

-   **nao.txt** contains North Atlantic Oscillation (NAO) indices
    derived from
    <https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-pc-based>

-   **climate.data.RDS**: Climate data (air temperature, precipitation
    levels and days with snow) from eight stations of the German Weather
    Service (DWD).

## **`R`**

-   **custom_functions.R**: R scripts implementing several data
    wrangling steps, including the creation of capture histories for
    capture-mark-resighting models, bootstrapping procedures etc.

## **plots**

Contains all figures that are presented in the manuscript or
supplementary files.
