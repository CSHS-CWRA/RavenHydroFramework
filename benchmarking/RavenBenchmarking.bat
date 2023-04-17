REM Benchmarking Batch File For Raven version evaluation
echo off
Setlocal EnableDelayedExpansion
cls

REM Location of Working directory (no end slash)
SET workingdir=%~d0\James\Software\Raven\BenchmarkProtocol

REM version name 
SET ver_name=vTest

echo Benchmarking Raven version %ver_name%...

REM Location of Raven executable
SET ravexe=%workingdir%\_Executables\%ver_name%\Raven.exe

if not exist %ravexe% (
  echo raven file executable %ravexe% doesn't exist. BENCHMARKING FAILED.
  pause
  exit
) 

if not exist %workingdir% (
  echo working directory %workingdir% doesn't exist. BENCHMARKING FAILED.
  pause
  exit
)

mkdir %workingdir%\out_%ver_name%\Alouette\
chdir %workingdir%\_InputFiles\Alouette\
%ravexe% Alouette_ws -o %workingdir%\out_%ver_name%\Alouette\

mkdir %workingdir%\out_%ver_name%\LaJoie\
chdir %workingdir%\_InputFiles\LaJoie\
%ravexe% La_Joie_ws -o %workingdir%\out_%ver_name%\LaJoie\

mkdir %workingdir%\out_%ver_name%\Revelstoke\
chdir %workingdir%\_InputFiles\Revelstoke\
%ravexe% Revelstoke_ws -o %workingdir%\out_%ver_name%\Revelstoke\

mkdir %workingdir%\out_%ver_name%\Williston_Finlay\
chdir %workingdir%\_InputFiles\Williston_Finlay\
%ravexe% Williston_Finlay_ws -o %workingdir%\out_%ver_name%\Williston_Finlay\
 
mkdir %workingdir%\out_%ver_name%\Alouette2\
chdir %workingdir%\_InputFiles\Alouette2\
%ravexe% Alouette2 -o %workingdir%\out_%ver_name%\Alouette2\

mkdir %workingdir%\out_%ver_name%\Irondequoit\
chdir %workingdir%\_InputFiles\Irondequoit\
%ravexe% Irondequoit -o %workingdir%\out_%ver_name%\Irondequoit\

mkdir %workingdir%\out_%ver_name%\Nith\
chdir %workingdir%\_InputFiles\Nith\
%ravexe% Nith -o %workingdir%\out_%ver_name%\Nith\

mkdir %workingdir%\out_%ver_name%\York\
chdir %workingdir%\_InputFiles\York\
%ravexe% York_gridded_m_daily_i_daily -o %workingdir%\out_%ver_name%\York\

mkdir %workingdir%\out_%ver_name%\LOTW\
chdir %workingdir%\_InputFiles\LOTW\
%ravexe% LOWRL -o %workingdir%\out_%ver_name%\LOTW\

mkdir %workingdir%\out_%ver_name%\Salmon_GR4J\
chdir %workingdir%\_InputFiles\Salmon_GR4J\
%ravexe% raven-gr4j-salmon -o %workingdir%\out_%ver_name%\Salmon_GR4J\

mkdir %workingdir%\out_%ver_name%\Salmon_HBV\
chdir %workingdir%\_InputFiles\Salmon_HBV\
%ravexe% raven-hbv-salmon -o %workingdir%\out_%ver_name%\Salmon_HBV\

mkdir %workingdir%\out_%ver_name%\Salmon_HMETS\
chdir %workingdir%\_InputFiles\Salmon_HMETS\
%ravexe% raven-hmets-salmon -o %workingdir%\out_%ver_name%\Salmon_HMETS\

mkdir %workingdir%\out_%ver_name%\Salmon_MOHYSE\
chdir %workingdir%\_InputFiles\Salmon_MOHYSE\
%ravexe% raven-mohyse-salmon -o %workingdir%\out_%ver_name%\Salmon_MOHYSE\

mkdir %workingdir%\out_%ver_name%\Salmon_Mixed\
chdir %workingdir%\_InputFiles\Salmon_Mixed\
%ravexe% raven_weighted_processes -o %workingdir%\out_%ver_name%\Salmon_Mixed\

mkdir %workingdir%\out_%ver_name%\Salmon_HYMOD\
chdir %workingdir%\_InputFiles\Salmon_HYMOD\
%ravexe% raven-hymod-salmon -o %workingdir%\out_%ver_name%\Salmon_HYMOD\

mkdir %workingdir%\out_%ver_name%\Salmon_SACSMA\
chdir %workingdir%\_InputFiles\Salmon_SACSMA\
%ravexe% raven-sacsma-salmon -o %workingdir%\out_%ver_name%\Salmon_SACSMA\

mkdir %workingdir%\out_%ver_name%\Lievre\
chdir %workingdir%\_InputFiles\Lievre\
%ravexe% lievre -o %workingdir%\out_%ver_name%\Lievre\

mkdir %workingdir%\out_%ver_name%\York_nc1\
chdir %workingdir%\_InputFiles\York_nc\
%ravexe% York -o %workingdir%\out_%ver_name%\York_nc1\

mkdir %workingdir%\out_%ver_name%\York_nc2\
chdir %workingdir%\_InputFiles\York_nc\
%ravexe% York_gridded_m_daily_i_daily -o %workingdir%\out_%ver_name%\York_nc2\

mkdir %workingdir%\out_%ver_name%\York_nc3\
chdir %workingdir%\_InputFiles\York_nc\
%ravexe% York_gridded_m_daily_i_subdaily -o %workingdir%\out_%ver_name%\York_nc3\

mkdir %workingdir%\out_%ver_name%\York_nc4\
chdir %workingdir%\_InputFiles\York_nc\
%ravexe% York_gridded_m_subdaily_i_daily -o %workingdir%\out_%ver_name%\York_nc4\

mkdir %workingdir%\out_%ver_name%\York_nc5\
chdir %workingdir%\_InputFiles\York_nc\
%ravexe% York_gridded_m_subdaily_i_subdaily -o %workingdir%\out_%ver_name%\York_nc5\

mkdir %workingdir%\out_%ver_name%\York_nc6\
chdir %workingdir%\_InputFiles\York_nc\
%ravexe% York_nongridded_m_subdaily_i_daily -o %workingdir%\out_%ver_name%\York_nc6\

mkdir %workingdir%\out_%ver_name%\York_nc6\
chdir %workingdir%\_InputFiles\York_nc\
%ravexe% York_nongridded_m_subdaily_i_daily -o %workingdir%\out_%ver_name%\York_nc6\

mkdir %workingdir%\out_%ver_name%\ElbowHistoric\
chdir %workingdir%\_InputFiles\ElbowHistoric\
%ravexe% elbow -o %workingdir%\out_%ver_name%\ElbowHistoric\

mkdir %workingdir%\out_%ver_name%\ElbowHistoric_csv\
chdir %workingdir%\_InputFiles\ElbowHistoric_csv\
%ravexe% elbow -o %workingdir%\out_%ver_name%\ElbowHistoric_csv\

mkdir %workingdir%\out_%ver_name%\ElbowForecast\
chdir %workingdir%\_InputFiles\ElbowForecast\
%ravexe% elbow -o %workingdir%\out_%ver_name%\ElbowForecast\

mkdir %workingdir%\out_%ver_name%\Elbow-forecastrun\
chdir %workingdir%\_InputFiles\Elbow-forecastrun\work\
%ravexe% elbow -o %workingdir%\out_%ver_name%\Elbow-forecastrun\

mkdir %workingdir%\out_%ver_name%\Elbow-historicrun\
chdir %workingdir%\_InputFiles\Elbow-historicrun\work\
%ravexe% elbow -o %workingdir%\out_%ver_name%\Elbow-historicrun\

mkdir %workingdir%\out_%ver_name%\02EB011\
chdir %workingdir%\_InputFiles\02EB011\
%ravexe% raven-blended-02EB011 -o %workingdir%\out_%ver_name%\02EB011\


"C:\Program Files\Beyond Compare 4\BCompare.exe" %workingdir%\out_%ver_name%\ %workingdir%\out_v406\
 
echo -------------------------------------------
echo ...BENCHMARKING DONE.
echo -------------------------------------------
pause
