cmake_minimum_required(VERSION 3.20)

# Setup Project
PROJECT(Raven CXX Fortran)

# Create MFUSG Library
add_library(MFUSGLib
  MFUSGLib/cln2basu1.f
  MFUSGLib/cln2props1.f
  MFUSGLib/density.f
  MFUSGLib/disu2gncb1.f
  MFUSGLib/disu2gncn1.f
  MFUSGLib/Extern.f90
  MFUSGLib/fpath.f90
  MFUSGLib/glo2basu1.f
  MFUSGLib/glo2btnu1.f
  MFUSGLib/glo2sgbu1.f
  MFUSGLib/glo2sms-u1.f
  MFUSGLib/gmodules.f
  MFUSGLib/gwf2basu1.f
  MFUSGLib/gwf2bcf-lpf-u1.f
  MFUSGLib/gwf2chd7u1.f
  MFUSGLib/gwf2dpf1u1.f
  MFUSGLib/gwf2drn7u1.f
  MFUSGLib/gwf2drt8u.f
  MFUSGLib/gwf2ets8u1.f
  MFUSGLib/gwf2evt8u1.f
  MFUSGLib/gwf2fhb7u1.f
  MFUSGLib/gwf2gag7u1.f
  MFUSGLib/gwf2ghb7u1.f
  MFUSGLib/gwf2hfb7u1.f
  MFUSGLib/gwf2lak7u1.f
  MFUSGLib/gwf2pbju1.f
  MFUSGLib/gwf2QRT8u.f
  MFUSGLib/gwf2rch8u1.f
  MFUSGLib/gwf2riv7u1.f
  MFUSGLib/gwf2sfr7u1.f
  MFUSGLib/gwf2str7u1.f
  MFUSGLib/gwf2sub7u1.f
  MFUSGLib/gwf2wel7u1.f
  MFUSGLib/gwt2bndsu1.f
  MFUSGLib/gwt2dptu1.f
  MFUSGLib/lak_gag_sfr_modules.f
  MFUSGLib/mfusg.f
  MFUSGLib/MFUSGLib.vfproj
  MFUSGLib/parutl7.f
  MFUSGLib/pcgu7.f
  MFUSGLib/rvn2pbju1.f90
  MFUSGLib/sparse.f
  MFUSGLib/tvmu1.f
  MFUSGLib/tvmu2.f
  MFUSGLib/utl7u1.f
  MFUSGLib/xmd.f
  MFUSGLib/xmdlib_2.f
)

# Create Cpp Executable
include_directories(Raven)
add_executable(Raven
  Raven/Abstraction.cpp
  Raven/Advection.cpp
  Raven/Advection.h
  Raven/Albedo.h
  Raven/Assimilate.cpp
  Raven/Baseflow.cpp
  Raven/CapillaryRise.cpp
  Raven/ChannelXSect.cpp
  Raven/ChannelXSect.h
  Raven/CommonFunctions.cpp
  Raven/Convolution.cpp
  Raven/Convolution.h
  Raven/CropGrowth.cpp
  Raven/CropGrowth.h
  Raven/CustomOutput.cpp
  Raven/CustomOutput.h
  Raven/Decay.cpp
  Raven/Decay.h
  Raven/DepressionProcesses.cpp
  Raven/DepressionProcesses.h
  Raven/Diagnostics.cpp
  Raven/Diagnostics.h
  Raven/Drain.cpp
  Raven/EnergyTransport.cpp
  Raven/Evaporation.cpp
  Raven/Flush.cpp
  Raven/ForcingGrid.cpp
  Raven/ForcingGrid.h
  Raven/Forcings.cpp
  Raven/Forcings.h
  Raven/Gauge.cpp
  Raven/Gauge.h
  Raven/GlacierProcesses.cpp
  Raven/GlacierProcesses.h
  Raven/GlobalParams.cpp
  Raven/GlobalParams.h
  Raven/GreenAmpt.cpp
  Raven/GroundwaterModel.cpp
  Raven/GroundwaterModel.h
  Raven/GWRiverConnection.cpp
  Raven/GWRiverConnection.h
  Raven/GWSWProcessABC.cpp
  Raven/GWSWProcesses.h
  Raven/HeatConduction.cpp
  Raven/HeatConduction.h
  Raven/HRUGroups.cpp
  Raven/HydroProcessABC.cpp
  Raven/HydroProcessABC.h
  Raven/HydroUnits.cpp
  Raven/HydroUnits.h
  Raven/Infiltration.cpp
  Raven/Infiltration.h
  Raven/Interflow.cpp
  Raven/IrregularTimeSeries.cpp
  Raven/IrregularTimeSeries.h
  Raven/LandUseClass.cpp
  Raven/LatAdvection.cpp
  Raven/LatAdvection.h
  Raven/LateralExchangeABC.cpp
  Raven/LateralExchangeABC.h
  Raven/LatFlush.cpp
  Raven/MassRouting.cpp
  Raven/MFUSGpp.h
  Raven/Model.cpp
  Raven/Model.h
  Raven/ModelABC.h
  Raven/ModelEnsemble.cpp
  Raven/ModelEnsemble.h
  Raven/ModelForcingGrids.cpp
  Raven/ModelInitialize.cpp
  Raven/ModelParamCheck.cpp
  Raven/OpenWaterEvap.cpp
  Raven/OpenWaterEvap.h
  Raven/OrographicCorrections.cpp
  Raven/ParseEnsembleFile.cpp
  Raven/ParseGWFile.cpp
  Raven/ParseHRUFile.cpp
  Raven/ParseInitialConditionFile.cpp
  Raven/ParseInput.cpp
  Raven/ParseLib.cpp
  Raven/ParseLib.h
  Raven/ParseLiveFile.cpp
  Raven/ParsePropertyFile.cpp
  Raven/ParseTimeSeriesFile.cpp
  Raven/PartitionPrecip.cpp
  Raven/Percolation.cpp
  Raven/PotentialMelt.cpp
  Raven/PrairieSnow.cpp
  Raven/PrairieSnow.h
  Raven/Precipitation.h
  Raven/ProcessGroup.cpp
  Raven/ProcessGroup.h
  Raven/Properties.h
  Raven/PropertyClassABC.h
  Raven/Radiation.cpp
  Raven/Radiation.h
  Raven/Raven.rc
  Raven/Raven.sdf
  Raven/Raven.sln
  Raven/Raven.vcxproj
  Raven/Raven.vcxproj.filters
  Raven/Raven.vcxproj.user
  Raven/RavenIcon.ico
  Raven/RavenInclude.h
  Raven/RavenMain.cpp
  Raven/Recharge.cpp
  Raven/Reservoir.cpp
  Raven/Reservoir.h
  Raven/resource.h
  Raven/RiverReach.cpp
  Raven/RiverReach.h
  Raven/SnowAlbedo.cpp
  Raven/SnowBalance.cpp
  Raven/SnowMeltRefreeze.cpp
  Raven/SnowMovers.h
  Raven/SnowParams.cpp
  Raven/SnowTempEvolve.cpp
  Raven/SoilAndLandClasses.h
  Raven/SoilClass.cpp
  Raven/SoilEvaporation.cpp
  Raven/SoilParams.cpp
  Raven/SoilProfile.cpp
  Raven/SoilProfile.h
  Raven/SoilWaterMovers.h
  Raven/SoilWaterMovers2.h
  Raven/Solvers.cpp
  Raven/StandardOutput.cpp
  Raven/StateVariables.cpp
  Raven/StateVariables.h
  Raven/SubBasin.cpp
  Raven/SubBasin.h
  Raven/Sublimation.cpp
  Raven/SurfaceEnergyExchange.cpp
  Raven/SurfaceEnergyExchange.h
  Raven/TerrainClass.cpp
  Raven/TimeSeries.cpp
  Raven/TimeSeries.h
  Raven/TimeSeriesABC.cpp
  Raven/TimeSeriesABC.h
  Raven/Transformation.cpp
  Raven/TransientParam.cpp
  Raven/TransientParam.h
  Raven/Transport.cpp
  Raven/Transport.h
  Raven/UnitTesting.cpp
  Raven/UnitTesting.h
  Raven/UpdateForcings.cpp
  Raven/UTM_to_LatLong.cpp
  Raven/VegetationClass.cpp
  Raven/VegetationMovers.cpp
  Raven/VegetationMovers.h
  Raven/VegetationParams.cpp
)

source_group("Source Files\\Hydrological Processes\\Lake\\Glacier\\Wetland Processes" FILES 
  Raven/Abstraction.cpp
  Raven/DepressionProcesses.cpp
  Raven/GlacierProcesses.cpp
  Raven/OpenWaterEvap.cpp
)
source_group("Source Files\\Transport Processes" FILES 
  Raven/Advection.cpp
  Raven/Decay.cpp
  Raven/LatAdvection.cpp
  Raven/MassRouting.cpp
  Raven/Transformation.cpp
  Raven/Transport.cpp
)
source_group("Source Files\\_Driver" FILES 
  Raven/Assimilate.cpp
  Raven/Model.cpp
  Raven/ModelEnsemble.cpp
  Raven/ModelForcingGrids.cpp
  Raven/ModelInitialize.cpp
  Raven/ModelParamCheck.cpp
  Raven/RavenMain.cpp
  Raven/Solvers.cpp
)
source_group("Source Files\\Hydrological Processes\\Soil Water Processes" FILES 
  Raven/Baseflow.cpp
  Raven/CapillaryRise.cpp
  Raven/GreenAmpt.cpp
  Raven/Infiltration.cpp
  Raven/Interflow.cpp
  Raven/Percolation.cpp
  Raven/SoilEvaporation.cpp
)
source_group("Source Files\\Profiles\\XSects\\Reservoirs" FILES 
  Raven/ChannelXSect.cpp
  Raven/Reservoir.cpp
  Raven/SoilProfile.cpp
)
source_group("Source Files" FILES 
  Raven/CommonFunctions.cpp
  Raven/RiverReach.cpp
  Raven/StateVariables.cpp
  Raven/UnitTesting.cpp
  Raven/UTM_to_LatLong.cpp
)
source_group("Source Files\\Hydrological Processes" FILES 
  Raven/Convolution.cpp
  Raven/Flush.cpp
  Raven/HydroProcessABC.cpp
  Raven/LateralExchangeABC.cpp
  Raven/LatFlush.cpp
  Raven/PartitionPrecip.cpp
  Raven/ProcessGroup.cpp
)
source_group("Source Files\\Hydrological Processes\\Vegetation Processes" FILES 
  Raven/CropGrowth.cpp
  Raven/VegetationMovers.cpp
)
source_group("Source Files\\_Driver\\Output" FILES 
  Raven/CustomOutput.cpp
  Raven/Diagnostics.cpp
  Raven/StandardOutput.cpp
)
source_group("Source Files\\Groundwater\\GWSW Processes" FILES 
  Raven/Drain.cpp
  Raven/GWSWProcessABC.cpp
  Raven/Recharge.cpp
)
source_group("Source Files\\Heat&Energy" FILES 
  Raven/EnergyTransport.cpp
  Raven/HeatConduction.cpp
  Raven/SurfaceEnergyExchange.cpp
)
source_group("Source Files\\Forcing Functions" FILES 
  Raven/Evaporation.cpp
  Raven/Forcings.cpp
  Raven/OrographicCorrections.cpp
  Raven/PotentialMelt.cpp
  Raven/Radiation.cpp
  Raven/UpdateForcings.cpp
)
source_group("Source Files\\Forcing Functions\\Gauge\\Time Series\\ForcingGrid" FILES 
  Raven/ForcingGrid.cpp
  Raven/Gauge.cpp
  Raven/IrregularTimeSeries.cpp
  Raven/TimeSeries.cpp
  Raven/TimeSeriesABC.cpp
)
source_group("Source Files\\Parameters\\Classes" FILES 
  Raven/GlobalParams.cpp
  Raven/LandUseClass.cpp
  Raven/SoilClass.cpp
  Raven/TerrainClass.cpp
  Raven/TransientParam.cpp
  Raven/VegetationClass.cpp
)
source_group("Source Files\\Groundwater" FILES 
  Raven/GroundwaterModel.cpp
  Raven/GWRiverConnection.cpp
)
source_group("Source Files\\HRUs and Subbasins" FILES 
  Raven/HRUGroups.cpp
  Raven/HydroUnits.cpp
  Raven/SubBasin.cpp
)
source_group("Source Files\\_Driver\\Input Parsing" FILES 
  Raven/ParseEnsembleFile.cpp
  Raven/ParseGWFile.cpp
  Raven/ParseHRUFile.cpp
  Raven/ParseInitialConditionFile.cpp
  Raven/ParseInput.cpp
  Raven/ParseLib.cpp
  Raven/ParseLiveFile.cpp
  Raven/ParsePropertyFile.cpp
  Raven/ParseTimeSeriesFile.cpp
)
source_group("Source Files\\Hydrological Processes\\Snow Water Processes" FILES 
  Raven/PrairieSnow.cpp
  Raven/SnowAlbedo.cpp
  Raven/SnowBalance.cpp
  Raven/SnowMeltRefreeze.cpp
  Raven/SnowTempEvolve.cpp
  Raven/Sublimation.cpp
)
source_group("Source Files\\Hydrological Routines" FILES 
  Raven/SnowParams.cpp
  Raven/SoilParams.cpp
  Raven/VegetationParams.cpp
)
source_group("Header Files\\Transport" FILES 
  Raven/Advection.h
  Raven/Decay.h
  Raven/LatAdvection.h
  Raven/Transport.h
)
source_group("Header Files\\Hydrological Processes" FILES 
  Raven/Albedo.h
  Raven/Convolution.h
  Raven/CropGrowth.h
  Raven/DepressionProcesses.h
  Raven/GlacierProcesses.h
  Raven/HydroProcessABC.h
  Raven/Infiltration.h
  Raven/LateralExchangeABC.h
  Raven/OpenWaterEvap.h
  Raven/PrairieSnow.h
  Raven/Precipitation.h
  Raven/ProcessGroup.h
  Raven/SnowMovers.h
  Raven/SoilWaterMovers.h
  Raven/VegetationMovers.h
)
source_group("Header Files\\Profiles\\XSects\\Reservoirs" FILES 
  Raven/ChannelXSect.h
  Raven/Reservoir.h
  Raven/SoilProfile.h
)
source_group("Header Files\\Input\\Output" FILES 
  Raven/CustomOutput.h
  Raven/Diagnostics.h
  Raven/ParseLib.h
)
source_group("Header Files\\Forcing Functions" FILES 
  Raven/ForcingGrid.h
  Raven/Forcings.h
  Raven/Radiation.h
)
source_group("Header Files\\Forcing Functions\\Gauge\\Time Series" FILES 
  Raven/Gauge.h
  Raven/IrregularTimeSeries.h
  Raven/TimeSeries.h
  Raven/TimeSeriesABC.h
)
source_group("Header Files\\Parameters\\Classes" FILES 
  Raven/GlobalParams.h
  Raven/SoilAndLandClasses.h
  Raven/TransientParam.h
)
source_group("Header Files\\Groundwater" FILES 
  Raven/GroundwaterModel.h
  Raven/GWRiverConnection.h
  Raven/GWSWProcesses.h
  Raven/MFUSGpp.h
)
source_group("Header Files" FILES 
  Raven/HeatConduction.h
  Raven/Properties.h
  Raven/PropertyClassABC.h
  Raven/RavenInclude.h
  Raven/RiverReach.h
  Raven/SoilWaterMovers2.h
  Raven/StateVariables.h
  Raven/SurfaceEnergyExchange.h
  Raven/UnitTesting.h
)
source_group("Header Files\\Subbasins\\HRUs" FILES 
  Raven/HydroUnits.h
  Raven/SubBasin.h
)
source_group("Header Files\\Model" FILES 
  Raven/Model.h
  Raven/ModelABC.h
  Raven/ModelEnsemble.h
)
source_group("Resource Files" FILES 
  Raven/resource.h
  Raven/Raven.rc
)
target_link_libraries(Raven MFUSGLib)
set_target_properties(Raven PROPERTIES LINKER_LANGUAGE CXX)
