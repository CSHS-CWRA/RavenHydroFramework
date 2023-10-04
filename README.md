# RavenHydroFramework

The code repository for the Raven Hydrological Modelling Framework developed at the University of Waterloo.

Release versions, tutorials, documentation, and more distributed at [the Raven website](http://raven.uwaterloo.ca/Main.html).

Intended for use with Visual Studio Community Edition 2021, but also provided with Windows/linux/unix/MacOS g++ makefile and CMake configuration file.

Note unconventional two-space tabbing conventions.
If you would like to work with the active development of Raven's core, please do so on a branch and coordinate commits to the trunk with the Raven development team.

Please contact us while you're at it - we love to have people helping out.

## Building Raven

Raven is built using CMake. A suggested sequence of commands to build Raven executable file is:

```bash
mkdir build
cd build
cmake [OPTIONS] ../
make
```

The `cmake` command can be configured with the following optional arguments (```[OPTIONS]``` above):

* `-DCOMPILE_LIB` to build Raven as a library (default: `OFF`)
* `-DCOMPILE_EXE` to build Raven as an executable (default: `ON`)

So that the ```cmake``` command to build Raven as solely a dynamic library becomes:

```bash
cmake -DCOMPILE_LIB=ON -DCOMPILE_EXE=OFF ../
```

# Refactoring plan

1. [```DONE```] Have a ```ExitGracefully``` method of the ```"CModel"``` class (keeping the public function);
   1. Requirement: Must be able to exit finely from anywhere in the code;
   2. *DONE*: Create method;
   3. *DONE*: Make the public function not use the static ```optStruct``` variable (uses the one from global CModel instead);
   4. Commit: [699b25e](https://github.com/adlzanchetta/RavenHydroFramework/commit/699b25ec47dd05c8df11007f8bdc8678d00979e2).
2. [```DONE```] Make ```static optStruct Options;``` a class member of the "CModel" class;
    1. **TODO**: Document changes done;
    2. ```Options``` usually becomes ```pModel->GetOptStruct()```;
    3. Commit: [c0b81b4](https://github.com/adlzanchetta/RavenHydroFramework/commit/c0b81b406ab0807f74bab6001d1ccdef752a26c1).
3. [```DONE```] Make ```CGlobalParams``` a class member of the "CModel" class;
    1. E.g., ```CGlobalParams::GetParams()``` becomes:
        1. ```this->_pGlobalParams->GetParams()``` (within ```CModel```);
        2. ```_pModel->GetGlobalParams()->GetParams()``` (outside ```CModel```).
     2. Commit: [9f6396c](https://github.com/adlzanchetta/RavenHydroFramework/commit/9f6396cee470f7ed2cc4ef32e49253814dc8abe5).
4. [```DONE```] Make ```CLandUseClass```  non-static and a class member of ```CModel```;
   1. At: ```LandUseClass.cpp```, ```SoilAndLandClasses.h``` (53 "static" words);
   2. attributes ```pAllLUClasses``` and ```NumLUClasses``` must go to ```CModel```;
   3. methods must be non-static;
   4. Commit: [d8ecaae](https://github.com/adlzanchetta/RavenHydroFramework/commit/d8ecaaeca934eebf2b7fbaa2e22501b356c9c48e).
5. [```DONE```] Make ```CSoilClass``` non-static and a class member of ```CModel```;
   1. At: ```CSoilClass.cpp```, ```SoilAndLandClasses.h``` (46 "static" words);
   2. attributes must be non-static;
   3. methods must be non-static;
6. [```DONE```] Make ```CVegetationClass``` non-static and a class member of ```CModel```;
   1. At: ```VegetationClass.cpp```, ```SoilAndLandClasses.h``` (39 "static" words);
   2. methods must be non-static;
7. [```DONE```] Make ```CTerrainClass``` non-static and a class member of ```CModel```;
   1. At: ```TerrainClass.cpp```, ```SoilAndLandClasses.h``` (33 "static" words);
   2. methods must be non-static;
8. [```TODO```] Make ```CSoilProfile``` non-static and a class member of ```CModel```;
   1. At: ```SoilProfile.cpp```, ```SoilAndLandClasses.h``` (31 "static" words);
   2. methods must be non-static;
9. [```TODO```] Make ```CChannelXSect``` non-static and a class member of ```CModel```;
   1. attributes must be non-static;
   2. methods must be non-static;
10. [```TODO```] Make ```CModel``` a non-static class;
    1. TODO - detail
    2. TODO - detail


## Notes:

- ```CHydroProcessABC```: has a ```CModelABC *pModel```;
