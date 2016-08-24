REM Benchmarking Batch File For Raven version evaluation
echo off
Setlocal EnableDelayedExpansion
cls
echo benchmarking raven...

REM Location of Working directory (no end slash)
SET workingdir=C:\James\Software\Raven\BenchmarkProtocol

REM version name 
SET ver_name=v619

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

echo -------------------------------------------
echo ...BENCHMARKING DONE.
echo -------------------------------------------
pause
