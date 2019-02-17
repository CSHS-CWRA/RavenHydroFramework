REM Benchmarking Batch File For Raven version evaluation
echo off
Setlocal EnableDelayedExpansion
cls

REM Location of Working directory (no end slash)
SET workingdir=C:\James\Software\Raven\BenchmarkProtocol

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

echo -------------------------------------------
echo ...BENCHMARKING DONE.
echo -------------------------------------------
pause
