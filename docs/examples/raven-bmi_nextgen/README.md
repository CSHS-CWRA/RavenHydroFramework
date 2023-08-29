# Raven BMI - NextGen

The examples assume you have a *shared library* version of Raven (```libraven.so``` in Linux, ```libraven.dll``` in Windows) and have an *executable* version of NOAA's NextGen (command ```ngen``` working).

All examples have correponding versions for being run as *Standalone*.


## Alouette5_60hours

Model for the two first days of 2016 at *hourly* temporal resolution (hypotetical disagegation of the daily data used in example ```Alouette4_60days```).

Single basing segmented into 23 HRUs.

**Before running:** the files ```../ng_Alouette5_realization.json``` and ```Alouette5.yaml``` need to be edited and every ```{PATH}``` needs to be replaced by the respective absolute path.

In a bash terminal (Linux):

```bash
$ mkdir -p Alouette5_60hours/output
$ cd Alouette5_60hours/output
$ ngen "../ng_catchment.geojson" "cat-3" "../ng_nexus.geojson" "nex-2" "../ng_Alouette5_realization.json"
```

At the end of the simulation, both the output files produced by NextGen (in this case: ```cat-3.csv``` and ```nex-2.csv```) and the output files produced by Raven will be at the working directory (```Alouette5_60hours/output```).

Implemented and tested using the commit [bf9e14b of NextGen](https://github.com/NOAA-OWP/ngen/commit/bf9e14be42dbe7a48b70c13c383e003a75c50229).