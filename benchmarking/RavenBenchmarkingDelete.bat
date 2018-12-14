REM Benchmarking Batch File For Raven version evaluation
echo off
Setlocal EnableDelayedExpansion
cls
echo deleting raven benchmark files...

REM Location of Working directory (no end slash)
SET workingdir=C:\James\Software\Raven\BenchmarkProtocol

REM version name 
SET ver_name=v142

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

del %workingdir%\out_%ver_name%\Alouette\*.*

del %workingdir%\out_%ver_name%\LaJoie\*.*

del %workingdir%\out_%ver_name%\Revelstoke\*.*

del %workingdir%\out_%ver_name%\Williston_Finlay\*.*
 
del %workingdir%\out_%ver_name%\Alouette2\*.*

del %workingdir%\out_%ver_name%\Irondequoit\*.*

del %workingdir%\out_%ver_name%\Nith\*.*

echo -------------------------------------------
echo ...BENCHMARK CLEANUP DONE.
echo -------------------------------------------
pause
