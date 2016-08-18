# <a name="top"></a> GT16\_JGRA Code and Input Data Supplement

---

This software package for regional climatological analysis using NOAA NCEI GHCN-Daily datasets is a component of the author's Ph.D. dissertation research at the University of Wisconsin–Madison and accompanies the paper

> Garcia, M., and P.A. Townsend (in review): "Recent climatological trends and potential influences on forest phenology around western Lake Superior, USA." Submitted to the *Journal of Geophysical Research – Atmospheres* on 5 April 2016, revised 20 August 2016.

That paper manuscript will be included in this package upon acceptance for publication. The corresponding author is Matthew Garcia (<matt.e.garcia@gmail.com>). Prof. Townsend is the author's dissertation advisor.

All scripts and modules (software) included here are licensed for free and fair use under the [Gnu GPL, version 3](https://www.gnu.org/licenses/) (see [LICENSE](./LICENSE_GnuGPLv3.txt) for details). Copyright on the original work (except where noted, especially community contributions) is retained by the author. See also the accompanying [DISCLAIMER](./DISCLAIMER.txt) and [CITATION](./CITATION.txt) pertaining to your use of this software.

All python scripts in this package were written for compatibility with python v2.7. We assume no responsibility for differences in functionality (and resulting errors) if you use these scripts with python v3.x. We will, however, attempt to work into ongoing upgrades those fixes compatible with python v3.x that do not break the python v2.7 functionality, so please still report those "bugs" as you find them.

---

## Package Documentation

[Acknowledgements](#acknowledgements)  
[Setup](#setup)  
[Modules](#modules)  
[Datasets and Docs](#datasets)  
[Tools](#tools)  
[Processing Steps](#processing)  
— [Basic gridded climatological analysis](#basic)  
— [Notes regarding our use of ParallelPython](#pp)  
— [Notes regarding stop–start capability](#stopstart)  
— Usage instructions for serial and checkpointed alternative script versions (coming soon)   
— [Further analysis using EPA Level-IV ecoregions](#further)  
— [Notes regarding our ecoregion clustering procedure](#cluster)

---

## <a name="acknowledgements"></a> Acknowledgements

Support for this work was provided by the USDA Forest Service–Northern Research Station, USDA NIFA McIntire–Stennis funding to the University of Wisconsin–Madison (project WIS01621), UW–Madison College of Agricultural and Life Sciences, NASA Biodiversity and Ecological Forecasting Project NNX14AC36G, the author’s 2015-2016 Wisconsin Space Grant Consortium Research Fellowship, and the author's Spring 2016 UW–Madison Department of Forest & Wildlife Ecology Research Fellowship. 

Portions of this work were developed using the computing resources and assistance of the UW–Madison Center For High Throughput Computing (CHTC) in the Department of Computer Sciences. The CHTC is supported by UW–Madison, the Advanced Computing Initiative, the Wisconsin Alumni Research Foundation, the Wisconsin Institutes for Discovery, and the National Science Foundation, and is an active member of the Open Science Grid, which is supported by the National Science Foundation and the US Department of Energy's Office of Science. 

Supporting entities assume neither responsibility for, nor ownership of, the intellectual property in this work.

[Back to top](#top)

---

## <a name="setup"></a> Setup

After you have cloned, forked, or otherwise downloaded this software package, in the main directory run 

`$ python setup.py`

which will

1. verify scripts, modules, tools, sample data, documents, and auxiliary files
2. verify availability of python dependencies used by the various scripts (listed below)
3. uncompress certain large example data files
4. build directory structure to receive script processing outputs
5. copy source scripts, modules, and tools into main package directory, then archive the originals 

[Back to top](#top)

---

## <a name="modules"></a> Modules

The **setup.py** script checks for several user-installed python packages that are called from the modules and scripts listed below:

* **os** is used in just a few scripts for file and directory management
* **sys** is used in all scripts for stdout messages and exit commands
* **datetime** is used at the beginning and end of every script to help the user gauge processing time
* **glob** is used for directory contents listings
* **numpy** is the essential mathematics module used in almost all scripts
* **pandas** is used for input meteorological data organization, cleaning, and sorting (NOTE: pandas v0.17 or above is required for the sorting syntax used here)
* **h5py** is used for almost all I/O operations (except csv and text files)
* **pp** is [ParallelPython](http://www.parallelpython.com/) for multi-processor application, used during the grid generation phase of analysis; serial versions of these scripts may also be developed
* **gdal** is used for map operations along with...
* **osgeo.osr** for conversion of map coordinates (input station locations and output plot boundaries)
* **scipy.interpolate** is used in several of the available spatial interpolation methods
* **scipy.ndimage** is used during generation of output map plots
* **scipy.stats** is used for time series correlation, linear regression, and significance testing
* **matplotlib.pyplot** is used for generation of all output plots and maps

All of these modules are invoked and/or aliased using

`import <module>` or `import <module> as <alias>`

near the head of the calling script, so subroutines are called according to the 

`<module>.<subroutine>` or `<alias>.<subroutine>`

convention. Note that ParallelPython requires unaliased function calls in the parallelized subroutines, so some scripts have mixed usage.

In addition to those, this package contains several original python modules:

* **Date\_Convert.py** converts between calendar date and day-of-year
* **Interpolation.py** contains several spatial interpolation methods using **numpy** and **scipy**
* **Plots.py** contains all plotting and mapping interaction with **matplotlib**
* **Read\_Header\_Files.py** is for use with ArcGIS-style header files that accompany binary datasets
* **Stats.py** contains routines providing a specific collection of statistics, including trends and *p*-values, using **numpy** and **scipy**
* **Teleconnections.py** is used to read and process NCEP and similarly-formatted climate teleconnection index datasets 
* **UTM\_Geo\_Convert.py** converts between lat/lon (geographic) and UTM coordinate systems using **gdal** and **osgeo.osr**

Except for the **Interpolation** module, all of these modules are invoked using

`from <module> import *`

near the head of the calling script, so subroutines are called directly. **Interpolation.py** is invoked using

`import Interpolation`

at the head of the calling script, and subroutines are called according to the 

`Interpolation.<subroutine>`

convention.

[Back to top](#top)

---

## <a name="datasets"></a> Datasets and Docs

### Example input datasets

Most of these example or sample input datasets are specific to our own selected study area and period. Note that the operating space for interpolations and masking (including ecoregion-based analyses) is in rectilinear coordinates with distances in meters. For our study that is UTM zone 15N.

* **NCEI\_WLS\_19830101-20151031.csv** obtained by on-line order from the [NOAA NCEI](http://www.ncei.noaa.gov/) web interface according to the [GHCN-Daily](http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn) instructions document listed below
* **NLCD\_2011\_WLS\_UTM15N.bil** is a binary data layer that as reprojected and clipped to our desired study area from the widely available [USGS NLCD 2011](http://www.mrlc.gov/nlcd2011.php) product using ArcGIS; the actual content of the NLCD data layer is not used in the analyses below, but a simplified classification based on this NLCD layer is the basemap shown in Figure 1 of our paper
* **NLCD\_2011\_WLS\_UTM15N.hdr** is the grid-defining header file to accompany the binary data layer; note that the UTM zone is given in the file name, but not in the file contents
* **EPA\_L4\_Ecoregions\_WLS\_polygonIDs.txt** is a text file indicating the correspondence of data layer polygon identifiers (integers) with ecoregion common IDs (alphanumeric) in our study area
* **EPA\_L4\_Ecoregions\_WLS\_UTM15N.bil** is a binary data layer that was gridded, reprojected, and clipped to our desired study area (matching the NLCD layer above) from the [EPA Level-IV ecoregion polygon product](http://archive.epa.gov/wed/ecoregions/web/html/level_iii_iv-2.html) using ArcGIS
* **EPA\_L4\_Ecoregions\_WLS\_UTM15N.hdr** is the grid-defining header file to accompany the binary data layer; note that the UTM zone is given in the file name, but not in the file contents
* **NCEP\_CPC\_AMO\_indices.csv** is the recent time series of monthly mean Atlantic Multidecadal Oscillation index values from the [NOAA Earth System Research Laboratory website](http://www.esrl.noaa.gov/psd/data/timeseries/AMO/)
* **NCEP\_CPC\_AO\_indices.csv** is the recent time series of monthly mean Arctic Oscillation index values from the [NOAA/NCEP Climate Prediction Center website](http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.shtml)
* **NCEP\_CPC\_ENSO\_indices.csv** is the recent time series of monthly mean El Niño–Southern Oscillation SST anomaly values from the [NOAA/NCEP Climate Prediction Center website](http://www.cpc.ncep.noaa.gov/products/precip/CWlink/MJO/enso.shtml)
* **NCEP\_CPC\_NAO\_indices.csv** is the recent time series of monthly mean North Atlantic Oscillation index values from the [NOAA/NCEP Climate Prediction Center website](http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/teleconnections.shtml)
* **NCEP\_CPC\_PDO\_indices.csv** is the recent time series of monthly mean Pacific Decadal Oscillation index values from the [NOAA National Centers for Environmental Information website](https://www.ncdc.noaa.gov/teleconnections/pdo/)
* **NCEP\_CPC\_PNA\_indices.csv** is the recent time series of monthly mean Pacific–North America pattern index values from the [NOAA/NCEP Climate Prediction Center website](http://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/pna.shtml)
* **NSIDC\_MIFL\_Superior\_Ice.csv** is the recent time series of annual dates for ice-on and ice-off conditions at Bayfield WI from the [NOAA National Snow and Ice Data Center website](http://nsidc.org/) and the [Madeline Island Ferry Line website](http://madferry.com/)
* **Query\_locations\_dates\_sample.csv** is a sample collection of query locations (borrowed from a colleague's actual study locations) for demonstration of the **query\_NCEI\_grids.py** tool listed below

### GHCN-Daily documents

* **How\_to\_get\_NCEI\_GHCND\_data.txt** describes in detail using the [NCEI](http://www.ncei.noaa.gov/) map-based interface to order [GHCN-Daily](http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn) datasets; see below for a note on what period you should request, based on your actual desired analysis period
* **NCEI\_GHCND\_documentation.pdf** is a copy of the [GHCN-Daily](http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn) dataset documentation available through the [NCEI](http://www.ncei.noaa.gov/) website

[Back to top](#top)

---

## <a name="processing"></a> Processing Steps

### <a name="basic"></a> Basic gridded climatological analysis

These procedures were developed to work with [GHCN-Daily](http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn) station temperature and precipitation datasets in CSV format obtained from the [NOAA National Centers for Environmental Information (NCEI)](http://www.ncei.noaa.gov/). The parameters of the GHCN-Daily dataset are given in the included [NCEI GHCN-Daily documentation (pdf)](./docs/NCEI_GHCND_documentation.pdf). Additional references regarding GHCN-Daily QA/QC procedures are listed in our paper. For instructions on using of the NCEI website map interface to make data requests, and especially the variables required for this analysis, see [How to get NCEI GHCND data](./docs/How_to_get_NCEI_GHCND_data.txt). Other datasets available from NOAA/NCEI and other sources (e.g. WUnderground) might also work, as long as the temperature and precipitation observations are accompanied by data quality flags. However, in our experience, interpolated fields using raw data without quality indicators are quite messy. To proceed without QA/QC indicators, we would need to develop a workaround for the first script in this sequence.

Your daily meteorological station dataset should start *at least* 6 months (preferably a full calendar year) prior to the start of your desired climatological analysis period, because of the way that these climatological derivatives are accumulated over as much as one year, and because some of the critical seasonal indicators (CD, GDD) are reset on alternating cycles. This is the reason that our own dataset begins on 1 Jan 1983 but our analysis period actually begins in 1984.

Note that the user is expected to proceed through these scripts in numerical sequence, though some tolerances are included in case you need to go back to change and re-run a script. For those scripts that use ParallelPython, some specific notes regarding execution are provided below. For some of these, serial versions are also provided in this code package. 

1. **process\_NCEI\_00.py**  
<u>Function</u>: QA/QC ('cleaning') of daily meteorological station data in NOAA/NCEI GHCN-Daily datasets  
<u>Usage</u>: `python process_NCEI_00.py NCEI_WLS_19830101-20151031.csv ./data`  
<u>Input</u>: 1 station meteorological data file from NOAA/NCEI in '.csv' format (the script knows that this file will be in the 'data' subdirectory)  
<u>Output</u>: 1 '.csv' file with the cleaned version of the input dataset; 1 '.csv' file with an accounting of errors cleaned, listed by station and variable; 1 '.h5' file with preliminary metadata (all of these will be generated in the 'data' subdirectory)  
<u>Methods</u>: Most QA/QC decisions are based on missing values (location and/or observation) and reported data flags; trace precipitation (reported as 0) is adjusted to a value of 0.1 mm; for any flag indicating a possibly erroneous observation, that value is set to indicate missing data; currently checks for outliers in temperature observations (see script for details) but does not yet check for outliers in precipitation observations.  
<u>Notes</u>: There are some known issues with cooperative precipitation reports that are not yet corrected, such as a 'date shift' problem. We're thinking about them, though...

2. **process\_NCEI\_01.py**  
<u>Function</u>: Extraction of daily meteorological station data from cleaned NOAA/NCEI dataset  
<u>Usage</u>: `python process_NCEI_01.py NCEI_WLS_19830101-20151031 ./data`  
<u>Input</u>: 2 output files from **process\_NCEI\_00.py** in '.csv' and '.h5' formats (same root file name, in 'data' subdirectory)  
<u>Output</u>: Updated input '.h5' file with sorted meteorological data (no new files)  

3. **process\_NCEI\_02.py**  
<u>Function</u>: Gridded interpolation of daily Prcp/Tmax/Tmin station data via user's method of choice, and calculation of daily Tavg field  
<u>Usage</u>: `python process_NCEI_02.py NLCD_2011_WLS_UTM15N NCEI_WLS_19830101-20151031 ./grids 480 RBF 1`  
where 'NLCD\_2011\_WLS\_UTM15N' is the root of the header file name that defines the study area (grid geographic location and extent) in the 'data' subdirectory, './grids' is the path to the desired output location, '480' is the desired interpolation output grid resolution (in meters), 'RBF' is the desired spatial interpolation method (see Notes below), and '1' is the default flag value for daily map output graphics  
<u>Input</u>: 1 NLCD (or other) binary grid header file in text format with '.hdr' extension (in 'data' subdirectory); 1 output file from **process\_NCEI\_01.py** in '.h5' format (in 'data' subdirectory)   
<u>Output</u>: Daily '.h5' files with original meteorological data and four gridded fields (1 new '.h5' file per day, in 'grids' subdirectory); corresponding daily mapped variables (4 new '.png' files, in 'images' subdirectory, if requested)   
<u>Methods</u>: The operating space for interpolations is in rectilinear coordinates (UTM, with distances in meters); there are 4 spatial interpolation methods in **Interpolation.py** currently available to the user:  
— RBF: radial basis functions, a **scipy.interpolate** built-in multiquadric method (the fastest and least computationally expensive method, according to our tests)  
— CSP: cubic splines via **griddata**, a **scipy.interpolate** built-in method  
— BSP: bivariate cubic B-splines, a **scipy.interpolate** built-in method  
— IDW: inverse-distance-squared (for temperature) and -cubed (for precipitation)   
<u>Notes</u>: This script uses ParallelPython (see notes below) but a serial version is also available; if a study area NLCD or other grid is not available, the user can spoof the required header file (see script for details)   
<u>To Do</u>: Instructions for use of the serial and version of this script and its preprocessing 'helper' script will be provided soon

4. **process\_NCEI\_03.py**  
<u>Function</u>: Temporal accumulation of various gridded daily meteorological data values as climatological derivatives   
<u>Usage</u>: `python process_NCEI_03.py NCEI_WLS_19830101-20151031 ./grids`  
<u>Input</u>: 1 output file from **process\_NCEI\_01.py** in '.h5' format (in 'data' subdirectory); daily output files from **process\_NCEI\_02.py** in '.h5' format (in 'grids' subdirectory)   
<u>Output</u>: Copied daily input files with new accumulation grids (1 new '.h5' file per day, in 'grids' subdirectory)  
<u>Notes</u>: The following 15 climatological indicators are calculated:   
— FD: freezing days for Tmin, Tavg, and Tmax (counted from 1 Jul)   
— CD: chilling days using Tavg relative to 5°C (counted from 1 Jul)   
— CDD: chilling degree days using Tavg relative to 5°C (accumulated from 1 Jul)   
— GDD: growing degree days using Tavg relative to 5°C (accumulated from 1 Jan)   
— GDD\_base0: growing degree days using Tavg relative to 0°C (accumulated from 1 Jan)   
— Tx\_90d: 90-day mean and variance of Tmin, Tavg, and Tmax   
— P\_Xd: 30-day, 90-day, 180-day, and 365-day precipitation totals   
— P\_90d\_X: 90-day precipitation days at thresholds of 0, 10, and 25 mm   
<u>Notes</u>: This script uses ParallelPython (see notes below) but a serial version is also available; daily map output graphics of the calculated indicators are not yet available but would be easy to implement (though they would inflate your 'images' subdirectory tremendously)   
<u>To Do</u>: Instructions for use of the serial and checkpointed versions of this script will be provided soon

5. **process\_NCEI\_04.py**  
<u>Function</u>: Aggregation of climatological derivatives for specific dates and time periods over the desired analysis period. A total of 89 gridded climatological indicators are derived or calculated on an annual basis, and 109 climatological statistics grids (each with 7 indicators) are derived or calculated for the entire designated study period. Examples: CD values at VEQ (vernal equinox) and other seasonal boundaries, finding CD between VEQ and SSOL (summer solstice), finding beginning and end of CD plateau, calculating length (in both days and GDD) of CD plateau, and calculating statistics on all of these (mean/stdev/trend/*p*-value over analysis period)  
<u>Usage</u>: `python process_NCEI_04.py 1984 2013 ./grids`  
where the beginning and ending years of the analysis period are given  
<u>Input</u>: Daily output files from **process\_NCEI\_03.py** in '.h5' format (in 'grids' subdirectory)  
<u>Output</u>: 1 new '.h5' file with aggregated grid datacubes and statistics grids (in 'analyses' subdirectory)  
<u>Notes</u>: This script uses ParallelPython (see notes below) but a serial version is also available in two parts (for more economical memory usage) as **process\_NCEI\_04a.py** and **process\_NCEI\_04b.py** to be executed in sequence   
<u>To Do</u>: Specific instructions for use of the serial 2-part version of this script will be provided soon

6. **process\_NCEI\_05.py**  
<u>Function</u>: Summarize several seasonal temperature and precipitation and cold/warm season statistics on study area grids, grid-mean annual time series (and their statistics), and cross-indicator/season correlations and their *p*-values, including autumn-to-winter comparisons (to close the annual cycle of potential season-to-season correlations)  
<u>Usage</u>: `python process_NCEI_05.py 1984 2013 ./analyses`  
where the beginning and ending years of the analysis period are given  
<u>Input</u>: 1 output file from **process\_NCEI\_04.py** in '.h5' format (in 'analyses' subdirectory)  
<u>Output</u>: Several '.csv' files with aggregated seasonal statistics (in 'analyses' subdirectory)  

7. **process\_NCEI\_06.py**  
<u>Function</u>: Summarize key meteorological variables on a daily basis  
<u>Usage</u>: `python process_NCEI_06.py 1984 2013 ./grids`  
where the beginning and ending years of the analysis period are given  
<u>Input</u>: Daily output files from **process\_NCEI\_03.py** in '.h5' format (in 'grids' subdirectory)  
<u>Output</u>: 7 new '.csv' files and 1 new '.h5' file with daily time series and aggregated mean daily variable values (in 'analyses' subdirectory)  

8. **process\_NCEI\_07.py**  
<u>Function</u>: Plot key meteorological variables on a daily/seasonal/annual/decadal basis  
<u>Usage</u>: `python process_NCEI_07.py 1984 2013 ./analyses`  
where the beginning and ending years of the analysis period are given  
<u>Input</u>: Output files from **process\_NCEI\_05.py** in '.csv' format and **process\_NCEI\_06.py** in '.h5' format (in 'analyses' subdirectory)  
<u>Output</u>: Graphical analyses in '.png' format (in 'images/analyses' subdirectory)  

[Back to Processing Steps](#processing) — [Back to top](#top)

### <a name="pp"></a> Notes regarding our use of ParallelPython (pp)

[ParallelPython](http://www.parallelpython.com/) is used in three of the scripts above: **process\_NCEI\_02.py**, **process\_NCEI\_03.py**, and **process\_NCEI\_04.py**. In all of these, the code is written for execution on a single-node multi-processor (SMP) system. On a linux system you can run `$ nproc` or `$ lscpu` to find out how many processors/cores are available. Where the number of **pp** processors is requested in a script, that number does *not* include the processor on which the main script is already runnning. For example, if you have one 4-processor system, in SMP mode you can only request 3 processors for **pp** usage. If you have a cluster (multiple systems, each with multiple processors) you can change the **pp** specifications to take fuller advantage of your cluster capabilities, but the limitation then becomes the contruction of the loops involved—they would need to be rewritten entirely for better use of cluster parallelization, and then still only **process\_NCEI\_02.py** and **process\_NCEI\_04.py** would likely benefit from that particular change. 

Depending on the duration of your study period, the size of your study area (thus grid size), and the speed of your system/cluster CPUs, **process\_NCEI\_02.py** and **process\_NCEI\_03.py** can take a long time to run. For our study (524 x 611 grid, ~33 data years) on a 3.0GHz SMP system, **process\_NCEI\_02.py** using 3 processors (plus one for the master process) took ~30.5 hours, **process\_NCEI\_03.py** using 8 processors (plus one for the master process) took ~48 hours, and **process\_NCEI\_04.py** using 6 processors (plus one for the master process) took less than 5 hours.

Implementing **pp** is not entirely straightforward. Using the example of **process\_NCEI\_03.py** (simplified here for brevity) the essentials in a given script are:

```
import numpy as np  # note numpy can be aliased here
import pp

...

def datacube_sum(<input arguments>):
    ...  # statements in which numpy is NOT aliased
    return <output arguments>

def datacube_mean(<input arguments>):
    ...  # statements in which numpy is NOT aliased
    return <output arguments>

...

ncpus = 8
ppservers0 = ()  # empty tuple indicates SMP usage on local node
job_server0 = pp.Server(ncpus, ppservers=ppservers0)
dc_sum = pp.Template(job_server0, datacube_sum, depfuncs=(), modules=("numpy",))
dc_mean = pp.Template(job_server0, datacube_mean, depfuncs=(), modules=("numpy",))

...

for date in dates:
    ...
    jobs_dc_sum = []  # to filled with <datacube_sum input arguments> and then replaced by <datacube_sum output arguments>
    jobs_dc_mean = []  # to be filled with <datacube_mean input arguments> and then replaced by <datacube_mean output arguments>
    ...
    jobs_dc_sum.append(dc_sum.submit(<datacube_sum input arguments>))
    jobs_dc_mean.append(dc_mean.submit(<datacube_mean input arguments>))
    ...
    for job in jobs_dc_sum:
        returns = job()  # retrieves <datacube_sum output arguments> as a list
        a = returns[0]
        b = returns[1]
        ...
    for job in jobs_dc_mean:
        returns = job()  # retrieves <datacube_mean output arguments> as a list
        c = returns[0]
        d = returns[1]
        ...
    ...
...     
        
job_server0.destroy()  # releases those requested processors cleanly

```
You can specify multiple systems in the `ppservers` declaration, and multiple job servers as separate and unique `job_server` declarations (see the [ParallelPython](http://www.parallelpython.com/) site for details). You can specify any number of `pp.Template` to keep your parallelization clear and your function calls and responses from any confusion. Since we're sending `<input arguments>` to, and receiving `<output arguments>` from, the same `jobs_x` list in separate steps, it is good practice for these arguments to identify uniquely the data chunk being operated on in the parallel environment. The returned information is then an unambiguous item in the `jobs_x` list and can be properly assigned to your waiting variables.

Note that we do not execute any I/O with external files within a parallelized routine. Most of our I/O deals with HDF-5 files, and there *is* a parallel-capable version of the HDF-5 library available, but it's not common even on clusters. If you figure out how to use it, especially on an SMP system, please let us know.

As mentioned above, refactoring internal loops to take better advantage of a larger system/cluster would help with **process\_NCEI\_02.py** because each variable grid and date is calculated independently of all others. Instead of sending each variable on a given day to a different processor (thus 3 processors requested), clearing the returned grids, and then doing it again for the next day, we might queue up the thousands of days for distribution among any number of processors, speeding up the the overall process considerably. If we rewrite the loops in **process\_NCEI\_02.py** we will certainly post and explain that update here.

The script **process\_NCEI\_04.py** is well-parallelized, but not necessariy as far as it could possibly go. As currently written, parallelization operates at the level of individual variables on individual grid points to obtain linear regression statistics over the study period. Certainly, with 17 variables over potentially millions of grid points, more processors would help this part of the process move faster. Another loop, to obtain seasonal values in individual years, could be parallelized by year but probably would not add much speed to the process. As above, if we rewrite the loops in **process\_NCEI\_04.py** for better performance on parallel systems, we will certainly post and explain that update here.

The other script, **process\_NCEI\_03.py**, is a special case here: each day's calculations of accumulated climatological variable values depend on those of the previous days. We cannot therefore parallelize over the dates in the study period, but we have parallelized the calculations of numerous varibles within a single date, like in **process\_NCEI\_02.py**. While we may rewrite the loops there for better use of more parallel CPUs, the same would be difficult for **process\_NCEI\_03.py**. In the meantime, we are working on checkpoint-based operations to ease the memory footprint and computation time of this script (see below), and along the way may find a better way to parallelize methods as well.

Finally, watch your RAM usage (directly correlated with grid size) as well as resource availability for any other users of the system/cluster. If your process crashes because of RAM usage or a resource conflict with another system user, the notes regarding stop–start capability below may be helpful.

[Back to Processing Steps](#processing) — [Back to top](#top)

### <a name="stopstart"></a> Notes regarding stop–start capability

<u>This section has not yet been updated to reflect the availability and use of alternative serial and checkpointed versions of the scripts mentioned here.</u>   

If your **process\_NCEI\_02.py** execution stops because of server or other issues (other than a program bug) you can simply edit the loop control in the script to pick up where it left off. First, you need a count of how many **grids/\*\_1.h5** files were generated before the crash. In practice, it's likely be best to delete the last couple of dates generated in case of data corruption, and then take that file count. Say you have *n* files completed. Then in **process\_NCEI\_02.py** look for the line  
`for date in dates:`  
and edit it to start at the *n*th date:  
`for date in dates[n:]:`  
When you start the script running again (using the same command as the previous time) you should find that it picks up with the next date in the sequence. If not, stop the process, adjust your *n* accordingly, and try again.

If your **process\_NCEI\_03.py** execution stops, the only option currently available is to start it over from the beginning, because the temporal accumulation datacubes are not saved outside of the execution memory. We made that decision during development based on I/O speed vs processing time, and because a dedicated SMP system was available to us. However, with access to larger systems, we are currently working on an alternative version using saved variables at designated checkpoints (e.g. yearly) to minimize data loss and reprocessing time in case of a crash, and to optimize execution and file management on systems with limited disk space or processing times. We will post that new script version here with updated instructions as soon as it's completed and tested. 

Yes, it's all very complicated, but then you're a scientist. If you're using this package rather than something that ends with '.exe', you also likely have some knowledge of Python (or other programming languages), or you may have access to someone who can help you figure it out. We have confidence in you. If all else fails, email us.

[Back to Processing Steps](#processing) — [Back to top](#top)

### <a name="further"></a> Further analysis using EPA Level-IV ecoregions  

If you have a map of the [EPA Level-IV ecoregions](http://archive.epa.gov/wed/ecoregions/web/html/level_iii_iv-2.html) for your study area, such as the map in our paper corresponding to the sample binary file listed above (derived via ArcGIS from the EPA-provided shapefiles for our study region), or some other polygon-based division of your study area (e.g. states, counties), the following scripts provide further analyses where the climatological derivatives calculated above are aggregated over individual and clustered polygons/regions. Correspondence of study area time series with climatological teleconnections, using data from the [NOAA/NCEP Climate Prediction Center (CPC)](http://www.cpc.ncep.noaa.gov/) and the [NOAA National Snow and Ice Data Center (NSIDC)](https://nsidc.org/) (because of our own study area's proximity to Lake Superior), are also addressed in these further analyses.

1. **process\_NCEI\_08.py**  
<u>Function</u>: Generate land masks for analyses and map graphics, based on given map of EPA Level-IV ecoregions  
<u>Usage</u>: `python process_NCEI_08.py EPA_L4_Ecoregions_WLS_UTM15N ./data`  
where 'EPA\_L4\_Ecoregions\_WLS\_UTM15N' is the root of the EPA Level-IV ecoregion binary and header files in the 'data' subdirectory  
<u>Input</u>: 1 '.bil' and 1 '.hdr' file with the ecoregion map for the study area (same root file name, in 'data' subdirectory)  
<u>Output</u>: 1 file **clipped\_ecoregions.h5** (in 'data' subdirectory)  

2. **process\_NCEI\_09.py**  
<u>Function</u>: Generate maps of numerous climatological derivatives on annual and summary bases  
<u>Usage</u>: `python process_NCEI_09.py 1984 2013 ./analyses 0 1`  
where the beginning and ending years of the analysis period are given, '0' is the default flag value for annual maps, and '1' is the default flag value for overall study period summary maps  
<u>Input</u>: Aggregated grid datacubes and statistics grids in '.h5' file from **process\_NCEI\_04.py** (in 'analyses' subdirectory); ecoregion maps and grid information from **process\_NCEI\_08.py** (in 'data' subdirectory)  
<u>Output</u>: Many many map graphics in '.png' format (in 'analyses/annual\_maps' and 'analyses/summary\_maps' subdirectories)   
<u>Q</u>: *Would you say this script generates a plethora of maps?*  
<u>A</u>: Oh yes, you will have a plethora.

3. **process\_NCEI\_10.py**  
<u>Function</u>: Calculate climatological statistics over time on grid-wide and ecoregion areas  
<u>Usage</u>: `python process_NCEI_10.py ecoregion_polygonIDs.txt 1984 2013 ./analyses`  
where **ecoregion\_polygonIDs.txt** indicates the correspondence between ecoregion designations and their polygon ID values in **EPA\_L4\_Ecoregions\_WLS\_UTM15N.bil** (both in the 'data' subdirectory), and the beginning and ending years of the analysis period are given  
<u>Input</u>: Aggregated grid datacubes and statistics grids in '.h5' file from **process\_NCEI\_04.py** (in 'analyses' subdirectory); ecoregion maps and grid information from **process\_NCEI\_08.py** (in 'data' subdirectory)  
<u>Output</u>: Aggregated full-grid and ecoregion-based statistics in '.h5' and '.csv' files (in 'analyses' subdirectory); maps of individual ecoregions in '.png' format (in 'analyses/ecoregion\_maps' subdirectory)  

4. **process\_NCEI\_11.py**  
<u>Function</u>: Summarize climatological statistics for individual ecoregions and perform various statistical tests on time series across ecoregions for each variable; calculate dissimilarity (analogous to Euclidean distance in 4D space) for each variable and ecoregion pair for use in **process\_NCEI\_12.py** clustering procedure  
<u>Usage</u>: `python process_NCEI_11.py ecoregion_polygonIDs.txt 1984 2013 ./analyses`  
where **ecoregion\_polygonIDs.txt** is as for **process\_NCEI\_10.py**, and the beginning and ending years of the analysis period are given  
<u>Input</u>: Ecoregion-based climatological time series and statistics in '.h5' files from **process\_NCEI\_10.py** (in 'analyses' subdirectory)  
<u>Output</u>: Ecoregion-based time series correlation statistics in '.h5' files (in 'analyses' subdirectory)  
<u>Notes</u>: Dissimilarity is calculated as
<math fontsize="12pt" display="inline">
  <mrow>
    <mi>D</mi>
    <mo>=</mo>
    <mn>[(1</mn>
    <mo>-</mo>
    <mi>R</mi>
    <msup>
      <mn>)</mn>
      <mn>2</mn>
    </msup>
    <mo>+</mo>
    <msup>
      <msub>
        <mi>R</mi>
        <mi>s</mi>
      </msub>
      <mn>2</mn>
    </msup>
    <mo>+</mo>
    <mn>(1</mn>
    <mo>-</mo>
    <msub>
      <mi>T</mi>
      <mi>s</mi>
    </msub>
    <msup>
      <mn>)</mn>
      <mn>2</mn>
    </msup>
    <mo>+</mo>
    <mn>(1</mn>
    <mo>-</mo>
    <msub>
      <mi>L</mi>
      <mi>s</mi>
    </msub>
    <msup>
      <mn>)</mn>
      <mn>2</mn>
    </msup>
    <msup>
      <mn>]</mn>
      <mn>1/2</mn>
    </msup>
  </mrow>
</math>
as decribed in our paper using the Pearson correlation coefficient 
<math fontsize="12pt" display="inline">
  <mrow>
    <mi>R</mi>
  </mrow>
</math>
and the *p*-values (statistical significance) from the least-squares linear regression
<math fontsize="12pt" display="inline">
  <mrow>
    <mn>(</mn>
    <msub>
      <mi>R</mi>
      <mi>s</mi>
    </msub>
    <mn>),</mn>
  </mrow>
</math>
Student's t-test for differing means
<math fontsize="12pt" display="inline">
  <mrow>
    <mn>(</mn>
    <msub>
      <mi>T</mi>
      <mi>s</mi>
    </msub>
    <mn>),</mn>
  </mrow>
</math>
and Levine's test for differing variances
<math fontsize="12pt" display="inline">
  <mrow>
    <mn>(</mn>
    <msub>
      <mi>L</mi>
      <mi>s</mi>
    </msub>
    <mn>).</mn>
  </mrow>
</math>
Several other formulations were examined and are still included on commented lines in the script file, should the user wish to try another calculation.

5. **process\_NCEI\_12.py**  
<u>Function</u>: Use statistical test results to cluster ecoregions using a time-series similarity measure that was defined in **process\_NCEI\_11.py**, then examine the clusters for subsets, overlaps, opportunities to merge/separate, singleton ecoregions, etc.  
<u>Usage</u>: `python process_NCEI_12.py 1984 2013 ./analyses`  
where the beginning and ending years of the analysis period are given  
<u>Input</u>: Ecoregion time series statistical tests in '.h5' files from **process\_NCEI\_11.py** (in 'analyses' subdirectory)  
<u>Output</u>: Ecoregion clusters in '.txt' file (in 'analyses' subdirectory)  

6. **process\_NCEI\_13.py**  
<u>Function</u>: Calculate climatological statistics over time on ecoregion clusters (similar to **process\_NCEI\_10.py**)   
<u>Usage</u>: `python process_NCEI_13.py 1984 2013 ./analyses`  
where the beginning and ending years of the analysis period are given  
<u>Input</u>: Aggregated grid datacubes and statistics grids in '.h5' file from **process\_NCEI\_04.py** (in 'analyses' subdirectory); ecoregion maps and grid information from **process\_NCEI\_08.py** (in 'data' subdirectory); aggregated ecoregion-based time series in '.h5' file from **process\_NCEI\_10.py** (in 'analyses' subdirectory); ecoregion clusters in '.txt' file from **process\_NCEI\_12.py** (in 'analyses' subdirectory)  
<u>Output</u>: Aggregated ecoregion cluster-based statistics in '.h5' files (in 'analyses' subdirectory); maps of ecoregion clusters in '.png' format (in 'analyses/cluster\_maps' subdirectory)   

7. **process\_NCEI\_14.py**  
<u>Function</u>: Summarize climatological statistics for ecoregion clusters and perform various statistical tests on time series across clusters for each variable (similar to **process\_NCEI\_11.py**)  
<u>Usage</u>: `python process_NCEI_14.py 1984 2013 ./analyses`  
where the beginning and ending years of the analysis period are given  
<u>Input</u>: Ecoregion cluster time series and stats in '.h5' files from **process\_NCEI\_13.py** (in 'analyses' subdirectory)   
<u>Output</u>: Aggregated ecoregion cluster-based statistics in '.h5' and '.csv' files (in 'analyses' subdirectory)   

8. **process\_NCEI\_15.py**  
<u>Function</u>: Cross-correlation of key dates and statistics from variable time series by ecoregion cluster, along with several climatological (teleconnection) indices  
<u>Usage</u>: `python process_NCEI_15.py 1984 2013 ./analyses`  
where the beginning and ending years of the analysis period are given  
<u>Input</u>: Full-grid time series in '.h5' files from **process\_NCEI\_10.py** (in 'analyses' subdirectory); ecoregion cluster time series and stats in '.h5' files from **process\_NCEI\_13.py** (in 'analyses' subdirectory); climatological teleconnection time series from various sources in '.csv' files  (in 'data' subdirectory)  
<u>Output</u>: Ecoregion cluster-based time series correlation statistics in '.h5' files (in 'analyses' subdirectory)  

[Back to Processing Steps](#processing) — [Back to top](#top)

### <a name="cluster"></a> Notes regarding our ecoregion clustering procedure

In order to reproduce *exactly* the ecoregion cluster map and statistical results shown in our paper, it will be necessary to iterate on **process\_NCEI\_13.py** and **process\_NCEI\_14.py**. At the end of **process\_NCEI\_14.py** we found that several clusters were still quite similar to each other, with calculated dissimilarity
<math fontsize="12pt" display="inline">
  <mrow>
    <mn>(</mn>
    <mi>D</mi>
    <mn>)</mn>
  </mrow>
</math>
values more than 2 standard deviations below the mean inter-cluster dissimilarity value. We manually edited the ecoregion clusters '.txt' file (in the 'analyses' subdirectory) that was produced by **process\_NCEI\_12.py** to combine only those cluster pairs with 
<math fontsize="12pt" display="inline">
  <mrow>
    <mi>D</mi>
    <mo><</mo>
    <mn>&mu;</mn>
    <mo>-</mo>
    <mn>2</mn>
    <mo>&InvisibleTimes;</mo>
    <mi>&sigma;</mi>
  </mrow>
</math>
and then ran this new ecoregion clusters '.txt' file through **process\_NCEI\_13.py** and **process\_NCEI\_14.py** a second time. It was the product of this second iteration that was then passed on to **process\_NCEI\_15.py** for final cross-correlation analysis.

[Back to Processing Steps](#processing) — [Back to top](#top)

---

## <a name="tools"></a> Tools

1. **query\_NCEI\_grids.py**  
<u>Function</u>: Queries climatological derivatives grids for desired locations and dates   
<u>Usage</u>: `python query_NCEI_grids.py Query_locations_dates_sample.csv`  
<u>Input</u>: Input file of query locations and dates in '.csv' format; daily output files from **process\_NCEI\_03.py** in '.h5' format (in 'grids' subdirectory)  
<u>Output</u>: Copied input '.csv' file with new columns for various climatological derivative values (in 'data' subdirectory)  
<u>Notes</u>: Query locations are expected in decimal latitude and longitude pairs; query dates are expected in mm/dd/yy format

2. R script to obtain GHCN-Daily data from NCEI via REST API  
(contributed by UW–Madison Ph.D. student W. Beckett Hills)   
\*\*COMING SOON\*\*

[Back to top](#top)

---

README markdown generated using [MacDown](http://macdown.uranusjr.com/)
