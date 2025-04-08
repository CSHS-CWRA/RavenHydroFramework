/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2025 the Raven Development Team

  Includes CModel routines for writing output headers and contents:
    CModel::CloseOutputStreams()
    CModel::WriteOutputFileHeaders()
    CModel::WriteMinorOutput()
    CModel::WriteMajorOutput()
    CModel::SummarizeToScreen()
    CModel::RunDiagnostics()
    Ensim output routines
    NetCDF output routines
  ----------------------------------------------------------------*/
#include "Model.h"
#include "StateVariables.h"

#if defined(_WIN32)
#include <direct.h>
#elif defined(__linux__)
#include <sys/stat.h>
#elif defined(__unix__)
#include <sys/stat.h>
#elif defined(__APPLE__)
#include <sys/stat.h>
#endif
int  NetCDFAddMetadata  (const int fileid,const int time_dimid,                  string shortname,string longname,string units);
int  NetCDFAddMetadata2D(const int fileid,const int time_dimid,int nbasins_dimid,string shortname,string longname,string units);
void WriteNetCDFGlobalAttributes(const int out_ncid,const optStruct &Options,const string descript);
void AddSingleValueToNetCDF     (const int out_ncid,const string &label,const size_t time_index,const double &value);
void WriteNetCDFBasinList       (const int ncid,const int varid,const int varid_name,const CModel* pModel,bool is_res,const optStruct &Options);
//////////////////////////////////////////////////////////////////
/// \brief returns true if specified observation time series is the flow series for subbasin SBID
/// \param pObs [in] observation time series
/// \param SBID [in] subbasin ID
//
bool IsContinuousFlowObs(const CTimeSeriesABC *pObs,long long SBID)
{
 // clears up  terribly ugly repeated if statements
  if (pObs==NULL)                                   { return false; }
  if (pObs->GetLocID() != SBID)                     { return false; }
  if (pObs->GetType() != CTimeSeriesABC::TS_REGULAR){ return false; }
  return (!strcmp(pObs->GetName().c_str(),"HYDROGRAPH")); //name ="HYDROGRAPH"
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if specified observation time series is the flow series for subbasin SBID
/// \param pObs [in] observation time series
/// \param SBID [in] subbasin ID
//
bool IsContinuousLevelObs(const CTimeSeriesABC *pObs,long long SBID)
{
 // clears up  terribly ugly repeated if statements
  if (pObs==NULL)                                   { return false; }
  if (pObs->GetLocID() != SBID)                     { return false; }
  if (pObs->GetType() != CTimeSeriesABC::TS_REGULAR){ return false; }
  return (!strcmp(pObs->GetName().c_str(),"WATER_LEVEL")); //name ="WATER_LEVEL"
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if specified observation time series is the flow series for subbasin SBID
/// \param pObs [in] observation time series
/// \param SBID [in] subbasin ID
//
bool IsContinuousConcObs(const CTimeSeriesABC *pObs,const long long SBID, const int c)
{
 // clears up  terribly ugly repeated if statements
  if (pObs==NULL)                                   { return false; }
  if (pObs->GetLocID() != SBID)                     { return false; }
  if (pObs->GetType() != CTimeSeriesABC::TS_REGULAR){ return false; }
  if (pObs->GetConstitInd() != c                   ){ return false; }
  return ((!strcmp(pObs->GetName().c_str(),"STREAM_CONCENTRATION")) ||
          (!strcmp(pObs->GetName().c_str(),"STREAM_TEMPERATURE"  ))); //name ="STREAM_CONCENTRATION" or "STREAM_TEMPERATURE"
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if specified observation time series is the reservoir stage series for subbasin SBID
/// \param pObs [in] observation time series
/// \param SBID [in] subbasin ID
//
bool IsContinuousStageObs(CTimeSeriesABC *pObs,long long SBID)
{
 // clears up  terribly ugly repeated if statements
  if (pObs==NULL)                                   { return false; }
  if (pObs->GetLocID() != SBID)                     { return false; }
  if (pObs->GetType() != CTimeSeriesABC::TS_REGULAR){ return false; }
  return (!strcmp(pObs->GetName().c_str(),"RESERVOIR_STAGE")); //name ="RESERVOIR_STAGE"
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if specified observation time series is the reservoir inflow series for subbasin SBID
/// \param pObs [in] observation time series
/// \param SBID [in] subbasin ID
//
bool IsContinuousInflowObs(CTimeSeriesABC *pObs, long long SBID)
{
  if (pObs==NULL)                                   { return false; }
  if (pObs->GetLocID() != SBID)                     { return false; }
  if (pObs->GetType() != CTimeSeriesABC::TS_REGULAR){ return false; }
  return (!strcmp(pObs->GetName().c_str(),"RESERVOIR_INFLOW")); //name ="RESERVOIR_INFLOW"
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if specified observation time series is the reservoir inflow series for subbasin SBID
/// \param pObs [in] observation time series
/// \param SBID [in] subbasin ID
//
bool IsContinuousNetInflowObs(CTimeSeriesABC *pObs, long long SBID)
{
  if (pObs==NULL)                                   { return false; }
  if (pObs->GetLocID() != SBID)                     { return false; }
  if (pObs->GetType() != CTimeSeriesABC::TS_REGULAR){ return false; }
  return (!strcmp(pObs->GetName().c_str(),"RESERVOIR_NETINFLOW")); //name ="RESERVOIR_NETINFLOW"
}


//////////////////////////////////////////////////////////////////
/// \brief Adds output directory & prefix to base file name
/// \param filebase [in] base filename, with extension, no directory information
/// \param &Options [in] Global model options information
//
string FilenamePrepare(string filebase, const optStruct &Options)
{
  string fn;
  if (Options.run_name==""){fn=Options.output_dir+filebase;}
  else                     {fn=Options.output_dir+Options.run_name+"_"+filebase;}
  return fn;
}

//////////////////////////////////////////////////////////////////
/// \brief Closes output file streams
/// \details after end of simulation from Main() or in ExitGracefully; All file streams are opened in WriteOutputFileHeaders() routine
//
void CModel::CloseOutputStreams()
{
  for (int c=0;c<_nCustomOutputs;c++){
    _pCustomOutputs[c]->CloseFiles(*_pOptStruct);
  }
  _pTransModel->CloseOutputFiles();
  if ( _STORAGE.is_open()){ _STORAGE.close();}
  if (   _HYDRO.is_open()){   _HYDRO.close();}
  if (_FORCINGS.is_open()){_FORCINGS.close();}
  if (_RESSTAGE.is_open()){_RESSTAGE.close();}
  if ( _DEMANDS.is_open()){ _DEMANDS.close();}
  if (  _LEVELS.is_open()){  _LEVELS.close();}

  if (_pDO!=NULL){_pDO->CloseOutputStreams();}

#ifdef _RVNETCDF_

  int    retval;      // error value for NetCDF routines
  if (_HYDRO_ncid != -9)    {retval = nc_close(_HYDRO_ncid);    HandleNetCDFErrors(retval); }
  _HYDRO_ncid    = -9;
  if (_STORAGE_ncid != -9)  {retval = nc_close(_STORAGE_ncid);  HandleNetCDFErrors(retval); }
  _STORAGE_ncid  = -9;
  if (_FORCINGS_ncid != -9) {retval = nc_close(_FORCINGS_ncid); HandleNetCDFErrors(retval); }
  _FORCINGS_ncid = -9;
  if(_RESSTAGE_ncid != -9)  {retval = nc_close(_RESSTAGE_ncid); HandleNetCDFErrors(retval); }
  _RESSTAGE_ncid = -9;
  if(_RESMB_ncid != -9)     {retval = nc_close(_RESMB_ncid);    HandleNetCDFErrors(retval); }
  _RESMB_ncid = -9;

#endif   // end compilation if NetCDF library is available
}


//////////////////////////////////////////////////////////////////
/// \brief Write output file headers
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CModel::WriteOutputFileHeaders(const optStruct &Options)
{
  int i,j,p;
  string tmpFilename;

  if(Options.noisy) { cout<<"  Writing Output File Headers..."<<endl; }

  if (Options.output_format==OUTPUT_STANDARD)
  {

    //WatershedStorage.csv
    //--------------------------------------------------------------
    if (Options.write_watershed_storage)
    {
      tmpFilename=FilenamePrepare("WatershedStorage.csv",Options);
      _STORAGE.open(tmpFilename.c_str());
      if (_STORAGE.fail()){
        ExitGracefully(("CModel::WriteOutputFileHeaders: unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
      }

      int iCumPrecip=GetStateVarIndex(ATMOS_PRECIP);
      _STORAGE<<"time [d],date,hour,rainfall [mm/day],snowfall [mm/d SWE],Channel Storage [mm],Reservoir Storage [mm],Rivulet Storage [mm]";
      for (i=0;i<GetNumStateVars();i++){
        if ((CStateVariable::IsWaterStorage(_aStateVarType[i]) && (i!=iCumPrecip))){
            _STORAGE << "," << CStateVariable::GetStateVarLongName(_aStateVarType[i],_aStateVarLayer[i],_pTransModel) << " [mm]";
        }
      }
      _STORAGE<<", Total [mm], Cum. Inputs [mm], Cum. Outflow [mm], MB Error [mm]"<<endl;
    }

    //Hydrographs.csv
    //--------------------------------------------------------------
    tmpFilename=FilenamePrepare("Hydrographs.csv",Options);
    _HYDRO.open(tmpFilename.c_str());
    if (_HYDRO.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }

    CSubBasin *pSB;
    _HYDRO<<"time,date,hour";
    _HYDRO<<",precip [mm/day]";
    for (p=0;p<_nSubBasins;p++)
    {
      pSB=_pSubBasins[p];
      if (pSB->IsGauged() && pSB->IsEnabled())
      {
        if (pSB->GetName()=="")      {_HYDRO<<",ID="<<pSB->GetID()  <<" [m3/s]";}
        else                         {_HYDRO<<","   <<pSB->GetName()<<" [m3/s]";}

        for (i = 0; i < _nObservedTS; i++){
          if (IsContinuousFlowObs(_pObservedTS[i],pSB->GetID()))
          {
            if (pSB->GetName()=="")  {_HYDRO<<",ID="<<pSB->GetID()  <<" (observed) [m3/s]";}
            else                     {_HYDRO<<","   <<pSB->GetName()<<" (observed) [m3/s]";}
          }
        }
        if (Options.write_localflow) {
          if (pSB->GetName()=="")    {_HYDRO<<",ID="<<pSB->GetID()  <<" (local) [m3/s]";}
          else                       {_HYDRO<<","   <<pSB->GetName()<<" (local) [m3/s]";}
        }
        if (pSB->GetReservoir() != NULL)
        {
          if (pSB->GetName()=="")    {_HYDRO<<",ID="<<pSB->GetID()  <<" (res. inflow) [m3/s]";}
          else                       {_HYDRO<<","   <<pSB->GetName()<<" (res. inflow) [m3/s]";}
          if (Options.write_netresinflow){
            if (pSB->GetName()=="")    {_HYDRO<<",ID="<<pSB->GetID()  <<" (res. net inflow) [m3/s]";}
            else                       {_HYDRO<<","   <<pSB->GetName()<<" (res. net inflow) [m3/s]";}
          }

          for(i = 0; i < _nObservedTS; i++){
            if(IsContinuousInflowObs(_pObservedTS[i],pSB->GetID()))
            {
              if (pSB->GetName()==""){_HYDRO<<",ID="<<pSB->GetID()  <<" (obs. res. inflow) [m3/s]";}
              else                   {_HYDRO<<","   <<pSB->GetName()<<" (obs. res. inflow) [m3/s]";}
            }
          }
        }
      }
    }
    _HYDRO<<endl;

    //WaterLevels.csv
    //--------------------------------------------------------------
    if (Options.write_waterlevels){
      tmpFilename=FilenamePrepare("WaterLevels.csv",Options);
      _LEVELS.open(tmpFilename.c_str());
      if (_LEVELS.fail()){
        ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
      }

      CSubBasin *pSB;
      _LEVELS<<"time,date,hour";
      _LEVELS<<",precip [mm/day]";
      for (p=0;p<_nSubBasins;p++)
      {
        pSB=_pSubBasins[p];
        if (pSB->IsGauged() && pSB->IsEnabled())
        {
          if (pSB->GetName()=="")      {_LEVELS<<",ID="<<pSB->GetID()  <<" [m]";}
          else                         {_LEVELS<<","   <<pSB->GetName()<<" [m]";}

          for (i = 0; i < _nObservedTS; i++){
            if (IsContinuousLevelObs(_pObservedTS[i],pSB->GetID()))
            {
              if (pSB->GetName()=="")  {_LEVELS<<",ID="<<pSB->GetID()  <<" (observed) [m]";}
              else                     {_LEVELS<<","   <<pSB->GetName()<<" (observed) [m]";}
            }
          }
        }
      }
      _LEVELS<<endl;
    }

    //ReservoirStages.csv
    //--------------------------------------------------------------
    if(Options.write_reservoir)
    {
      tmpFilename=FilenamePrepare("ReservoirStages.csv",Options);
      _RESSTAGE.open(tmpFilename.c_str());
      if(_RESSTAGE.fail()) {
        ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
      }

      _RESSTAGE<<"time,date,hour";
      _RESSTAGE<<",precip [mm/day]";
      for(p=0;p<_nSubBasins;p++)
      {
        if((_pSubBasins[p]->IsGauged()) && (_pSubBasins[p]->IsEnabled()) && (_pSubBasins[p]->GetReservoir()!=NULL)) {
          if(_pSubBasins[p]->GetName()=="") { _RESSTAGE<<",ID="<<_pSubBasins[p]->GetID()  <<" "; }
          else                              { _RESSTAGE<<","   <<_pSubBasins[p]->GetName()<<" "; }


          for(i = 0; i < _nObservedTS; i++) {
            if(IsContinuousStageObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
            {
              if(_pSubBasins[p]->GetName()=="") { _RESSTAGE<<",ID="<<_pSubBasins[p]->GetID()  <<" (observed) [m]"; }
              else                              { _RESSTAGE<<","   <<_pSubBasins[p]->GetName()<<" (observed) [m]"; }
            }
          }
        }
      }
      _RESSTAGE<<endl;
    }

    //ReservoirMassBalance.csv
    //--------------------------------------------------------------
    if (Options.write_reservoirMB)
    {
      ofstream RES_MB;
      string name;
      tmpFilename=FilenamePrepare("ReservoirMassBalance.csv",Options);
      RES_MB.open(tmpFilename.c_str());
      if (RES_MB.fail()){
        ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
      }
      RES_MB<<"time,date,hour";
      RES_MB<<",precip [mm/day]";
      for(p=0;p<_nSubBasins;p++){
        if((_pSubBasins[p]->IsGauged())  && (_pSubBasins[p]->IsEnabled())  && (_pSubBasins[p]->GetReservoir()!=NULL)) {

          if(_pSubBasins[p]->GetName()==""){ name=to_string(_pSubBasins[p]->GetID())+"="+to_string(_pSubBasins[p]->GetID()); }
          else                             { name=_pSubBasins[p]->GetName(); }
          RES_MB<<","   <<name<<" stage [m]";
          RES_MB<<","   <<name<<" area [m2]";
          RES_MB<<","   <<name<<" inflow [m3]";
          RES_MB<<","   <<name<<" outflow [m3]"; //from main outlet
          for (int i = 0; i < _pSubBasins[p]->GetReservoir()->GetNumControlStructures(); i++) {
            string name2=_pSubBasins[p]->GetReservoir()->GetControlName(i);
            RES_MB<<","   <<name<<" ctrl outflow "<<name2<<" [m3]";
            RES_MB<<","   <<name<<" ctrl regime "<<name2;
          }
          RES_MB<<","   <<name<<" precip [m3]";
          RES_MB<<","   <<name<<" evap [m3]";
          RES_MB<<","   <<name<<" seepage [m3]";
          RES_MB<<","   <<name<<" volume [m3]";
          RES_MB<<","   <<name<<" losses [m3]";
          RES_MB<<","   <<name<<" MB error [m3]";
          RES_MB<<","   <<name<<" constraint";
        }
      }
      RES_MB<<endl;
      RES_MB.close();
    }

    //ForcingFunctions.csv
    //--------------------------------------------------------------
    if(Options.write_forcings)
    {
      tmpFilename=FilenamePrepare("ForcingFunctions.csv",Options);
      _FORCINGS.open(tmpFilename.c_str());
      if (_FORCINGS.fail()){
        ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
      }
      _FORCINGS<<"time [d], date, hour, day_angle,";
      _FORCINGS<<" rain [mm/d], snow [mm/d], temp [C], temp_daily_min [C], temp_daily_max [C],temp_daily_ave [C],temp_monthly_min [C],temp_monthly_max [C],";
      _FORCINGS<<" air_density [kg/m3], air_pressure [KPa], rel_humidity [-],";
      _FORCINGS<<" cloud_cover [-],";
      _FORCINGS<<" ET_radiation [MJ/m2/d], SW_radiation [MJ/m2/d], net_SW_radiation [MJ/m2/d], LW_incoming [MJ/m2/d], net_LW_radiation [MJ/m2/d], wind_speed [m/s],";
      _FORCINGS<<" PET [mm/d], OW_PET [mm/d],";
      _FORCINGS<<" daily_correction [-], potential_melt [mm/d]";
      _FORCINGS<<endl;
    }
  }
  else if (Options.output_format==OUTPUT_ENSIM)
  {
    WriteEnsimStandardHeaders(Options);
  }
  else if (Options.output_format==OUTPUT_NETCDF)
  {
    WriteNetcdfStandardHeaders(Options);  // creates NetCDF files, writes dimensions and creates variables (without writing actual values)
  }

  //Demands.csv
  //--------------------------------------------------------------
  CSubBasin *pSB;
  if((Options.write_demandfile) && (Options.output_format!=OUTPUT_NONE))
  {
    tmpFilename=FilenamePrepare("Demands.csv",Options);
    _DEMANDS.open(tmpFilename.c_str());
    if(_DEMANDS.fail()) {
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }

    _DEMANDS<<"time,date,hour";
    for(p=0;p<_nSubBasins;p++) {
      pSB=_pSubBasins[p];
      if((pSB->IsEnabled()) && (pSB->IsGauged()) && (pSB->HasIrrigationDemand()))
      {
        string name;
        if(pSB->GetName()=="") { name="ID="+to_string(pSB->GetID()); }
        else                              { name=pSB->GetName(); }
        _DEMANDS<<","<<name<<" [m3/s]";
        _DEMANDS<<","<<name<<" (min.) [m3/s]";
        _DEMANDS<<","<<name<<" (total dem.) [m3/s]";
        _DEMANDS<<","<<name<<" (total del.) [m3/s]";
        _DEMANDS<<","<<name<<" (total ret.) [m3/s]";
        for (int ii=0;ii<pSB->GetNumWaterDemands(); ii++){
          name=pSB->GetWaterDemandObj(ii)->GetName();
          _DEMANDS<<","<<name<<" (demand) [m3/s]";
          _DEMANDS<<","<<name<<" (delivery) [m3/s]";
          if (pSB->GetWaterDemandObj(ii)->HasReturnFlow()){
          _DEMANDS<<","<<name<<" (return) [m3/s]";
          }
        }

        if (pSB->GetReservoir()!=NULL){
          for (int ii=0;ii<pSB->GetReservoir()->GetNumWaterDemands(); ii++){
            name=pSB->GetReservoir()->GetWaterDemandObj(ii)->GetName();
            _DEMANDS<<","<<name<<" (res. demand) [m3/s]";
            _DEMANDS<<","<<name<<" (res. delivery) [m3/s]";
            if (pSB->GetReservoir()->GetWaterDemandObj(ii)->HasReturnFlow()){
            _DEMANDS<<","<<name<<" (res. return) [m3/s]";
            }
          }
        }

      }
    }
    _DEMANDS<<endl;
  }

  //WatershedMassEnergyBalance.csv
  //--------------------------------------------------------------
  if (Options.write_mass_bal)
  {
    ofstream MB;
    tmpFilename=FilenamePrepare("WatershedMassEnergyBalance.csv",Options);
    MB.open(tmpFilename.c_str());
    if (MB.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    MB<<"time,date,hour";
    for (j=0;j<_nProcesses;j++){
      for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
        MB<<","<<GetProcessName(_pProcesses[j]->GetProcessType());
        MB<<"["<<CStateVariable::GetStateVarUnits(_aStateVarType[_pProcesses[j]->GetFromIndices()[q]])<<"]";
      }
    }
    MB<<endl;
    MB<<",,from:";
    for (j=0;j<_nProcesses;j++){
      for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
        sv_type typ=GetStateVarType (_pProcesses[j]->GetFromIndices()[q]);
        int     ind=GetStateVarLayer(_pProcesses[j]->GetFromIndices()[q]);
        MB << "," <<  _pStateVar->SVTypeToString(typ, ind);
      }
    }
    MB<<endl;
    MB<<",,to:";
    for (j=0;j<_nProcesses;j++){
      for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
        sv_type typ=GetStateVarType (_pProcesses[j]->GetToIndices()[q]);
        int     ind=GetStateVarLayer(_pProcesses[j]->GetToIndices()[q]);
        MB << "," << _pStateVar->SVTypeToString(typ, ind);
      }
    }
    MB<<endl;
    MB.close();
  }

  //HRUGroup_MassEnergyBalance.csv
  //--------------------------------------------------------------
  if (Options.write_group_mb!=DOESNT_EXIST)
  {
    int kk=Options.write_group_mb;
    ofstream HGMB;
    tmpFilename=FilenamePrepare(_pHRUGroups[kk]->GetName()+"_MassEnergyBalance.csv",Options);

    HGMB.open(tmpFilename.c_str());
    if (HGMB.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    HGMB<<"time,date,hour";
    for (j=0;j<_nProcesses;j++){
      for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
        HGMB<<","<<GetProcessName(_pProcesses[j]->GetProcessType());
        HGMB<<"["<<CStateVariable::GetStateVarUnits(_aStateVarType[_pProcesses[j]->GetFromIndices()[q]])<<"]";
      }
    }
    HGMB<<endl;
    HGMB<<",,from:";
    for (j=0;j<_nProcesses;j++){
      for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
        sv_type typ=GetStateVarType (_pProcesses[j]->GetFromIndices()[q]);
        int     ind=GetStateVarLayer(_pProcesses[j]->GetFromIndices()[q]);
        HGMB << "," << _pStateVar->SVTypeToString(typ, ind);
      }
    }
    HGMB<<endl;
    HGMB<<",,to:";
    for (j=0;j<_nProcesses;j++){
      for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
        sv_type typ=GetStateVarType (_pProcesses[j]->GetToIndices()[q]);
        int     ind=GetStateVarLayer(_pProcesses[j]->GetToIndices()[q]);
        HGMB << "," << _pStateVar->SVTypeToString(typ,ind);
      }
    }
    HGMB<<endl;
    HGMB.close();
  }

  //ExhaustiveMassBalance.csv
  //--------------------------------------------------------------
  if (Options.write_exhaustiveMB)
  {
    ofstream MB;
    tmpFilename=FilenamePrepare("ExhaustiveMassBalance.csv",Options);
    MB.open(tmpFilename.c_str());
    if (MB.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    MB<<"time,date,hour";
    bool first;
    for (i=0;i<_nStateVars;i++){
      if (CStateVariable::IsWaterStorage(_aStateVarType[i]))
      {
        MB << "," << _pStateVar->SVTypeToString(_aStateVarType[i],_aStateVarLayer[i]);
        first=true;
        for (j=0;j<_nProcesses;j++){
          for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
            if (_pProcesses[j]->GetFromIndices()[q]==i){
              if (!first){MB<<",";}first=false;
            }
          }
        }
        for (j=0;j<_nProcesses;j++){
          for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
            if (_pProcesses[j]->GetToIndices()[q]==i){
              if (!first){MB<<",";}first=false;
            }
          }
        }
        for(j=0;j<_nProcesses;j++){
          for(int q=0;q<_pProcesses[j]->GetNumLatConnections();q++){
            CLateralExchangeProcessABC *pProc=static_cast<CLateralExchangeProcessABC *>(_pProcesses[j]);
            if(pProc->GetLateralToIndices()[q]==i){
              if(!first){ MB<<","; }first=false;break;
            }
            if (pProc->GetLateralFromIndices()[q]==i){
              if (!first){MB<<",";}first=false;break;
            }
          }
        }
        MB<<",,,";//cum, stor, error
      }
    }
    MB<<endl;
    MB<<",,";//time,date,hour
    for (i=0;i<_nStateVars;i++){
      if (CStateVariable::IsWaterStorage(_aStateVarType[i]))
      {
        for (j=0;j<_nProcesses;j++){
          for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
            if (_pProcesses[j]->GetFromIndices()[q]==i){MB<<","<<GetProcessName(_pProcesses[j]->GetProcessType());}
          }
        }
        for (j=0;j<_nProcesses;j++){
          for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
            if (_pProcesses[j]->GetToIndices()[q]==i){MB<<","<<GetProcessName(_pProcesses[j]->GetProcessType());}
          }
        }
        for(j=0;j<_nProcesses;j++){
          for(int q=0;q<_pProcesses[j]->GetNumLatConnections();q++){
            CLateralExchangeProcessABC *pProc=static_cast<CLateralExchangeProcessABC *>(_pProcesses[j]);
            if (pProc->GetLateralToIndices  ()[q]==i){MB<<","<<GetProcessName(_pProcesses[j]->GetProcessType());break;}
            if (pProc->GetLateralFromIndices()[q]==i){MB<<","<<GetProcessName(_pProcesses[j]->GetProcessType());break;}
          }
        }
        MB<<",cumulative,storage,error";
      }
    }
    MB<<endl;
    MB.close();
  }

  // HRU Storage files
  //--------------------------------------------------------------
  if (_pOutputGroup!=NULL){
    for (int kk=0; kk<_pOutputGroup->GetNumHRUs();kk++)
    {
      ofstream HRUSTOR;
      tmpFilename="HRUStorage_"+to_string(_pOutputGroup->GetHRU(kk)->GetHRUID())+".csv";
      tmpFilename=FilenamePrepare(tmpFilename,Options);
      HRUSTOR.open(tmpFilename.c_str());
      if (HRUSTOR.fail()){
        ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
      }
      int iCumPrecip=GetStateVarIndex(ATMOS_PRECIP);
      HRUSTOR<<"time,date,hour,rainfall [mm/day],snowfall [mm/d SWE]";
      for (i=0;i<GetNumStateVars();i++){
        if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iCumPrecip)){
            HRUSTOR << "," << CStateVariable::GetStateVarLongName(_aStateVarType[i],_aStateVarLayer[i],_pTransModel) << " [mm]";
        }
      }
      HRUSTOR<<", Total [mm]"<<endl;
      HRUSTOR.close();
    }
  }

  // Assimilation adjustments
  //--------------------------------------------------------------
  if (Options.assimilate_flow){
    ofstream ASSIM;
    tmpFilename=FilenamePrepare("AssimilationAdjustments_BySubbasin.csv",Options);
    ASSIM.open(tmpFilename.c_str());
    if (ASSIM.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    ASSIM<<"time,date";
    
    for (int p = 0; p < _nSubBasins; p++) {
      if((_pSubBasins[p]->IsGauged()) && (_pSubBasins[p]->IsEnabled())){
        ASSIM<<","<<_pSubBasins[p]->GetID()<<" "; 
      }
    }
    ASSIM<<endl;
    for (int p = 0; p < _nSubBasins; p++) {
      if((_pSubBasins[p]->IsGauged()) && (_pSubBasins[p]->IsEnabled())){
        ASSIM<<", "; 
      }
    }
    ASSIM<<endl;
    ASSIM.close();
  }

  // Custom output files
  //--------------------------------------------------------------
  for (int c=0;c<_nCustomOutputs;c++)
  {
    _pCustomOutputs[c]->WriteFileHeader(Options);
  }

  // Transport output files
  //--------------------------------------------------------------
  _pTransModel->WriteOutputFileHeaders(Options);

  // Management output files
  //--------------------------------------------------------------
  if (Options.management_optimization) {
    _pDO->WriteOutputFileHeaders(Options);
  }

  //raven_debug.csv
  //--------------------------------------------------------------
  if (Options.debug_mode)
  {
    ofstream DEBUG;
    tmpFilename=FilenamePrepare("raven_debug.csv",Options);
    DEBUG.open(tmpFilename.c_str());
    if (DEBUG.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    DEBUG<<"time[d],date,hour,debug1,debug2,debug3,debug4,debug5,debug6,debug7,debug8,debug9,debug10"<<endl;
    DEBUG.close();
  }

  //opens and closes diagnostics.csv so that this warning doesn't show up at end of simulation
  //--------------------------------------------------------------
  if ((_nObservedTS>0) && (_nDiagnostics>0))
  {
    ofstream DIAG;
    tmpFilename=FilenamePrepare("Diagnostics.csv",Options);
    DIAG.open(tmpFilename.c_str());
    if(DIAG.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    DIAG.close();
  }

}

//////////////////////////////////////////////////////////////////
/// \brief Writes minor output to file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time *at the end of* the pertinent time step
//
void CModel::WriteMinorOutput(const optStruct &Options,const time_struct &tt)
{

  int     i,iCumPrecip,k;
  double  output_int = 0.0;
  double  mod_final = 0.0;
  double  S,currentWater;
  string  thisdate;
  string  thishour;
  bool    silent=true; //for debugging
  bool    quiet=true;  //for debugging
  double  t;

  CSubBasin* pSB;

  string tmpFilename;

  if ((tt.model_time==0) && (Options.suppressICs)){return;}

  if ((_pEnsemble != NULL) && (_pEnsemble->DontWriteOutput())) { return; } //specific to EnKF

  //converts the 'write every x timesteps' into a 'write at time y' value
  output_int = Options.output_interval * Options.timestep;
  mod_final  = ffmod(tt.model_time,output_int);

  iCumPrecip=GetStateVarIndex(ATMOS_PRECIP);

  if(fabs(mod_final) <= 0.5*Options.timestep)  //checks to see if sufficiently close to timestep
                                               //(this should account for any roundoff error in timestep calcs)
  {
    thisdate=tt.date_string;                   //refers to date and time at END of time step
    thishour=DecDaysToHours(tt.julian_day);
    t       =tt.model_time;

    time_struct prev;
    JulianConvert(t-Options.timestep,Options.julian_start_day,Options.julian_start_year,Options.calendar,prev); //get start of time step, prev

    double usetime=tt.model_time;
    string usedate=thisdate;
    string usehour=thishour;
    if(Options.period_starting){
      usedate=prev.date_string;
      usehour=DecDaysToHours(prev.julian_day);
      usetime=tt.model_time-Options.timestep;
    }

    // Console output
    //----------------------------------------------------------------
    if ((quiet) && (!Options.silent) && (tt.day_of_month==1) && ((tt.julian_day)-floor(tt.julian_day+TIME_CORRECTION)<Options.timestep/2))
    {
      cout<<thisdate <<endl;
    }
    if(!silent)
    {
      cout <<thisdate<<" "<<thishour<<":";
      if (t!=0){cout <<" | P: "<< setw(6)<<setiosflags(ios::fixed) << setprecision(2)<<GetAveragePrecip();}
      else     {cout <<" | P: ------";}
    }

    //Write current state of water storage in system to WatershedStorage.csv (ALWAYS DONE if not switched OFF)
    //----------------------------------------------------------------
    if (Options.output_format==OUTPUT_STANDARD)
    {
      if (Options.write_watershed_storage)
      {
        double snowfall      =GetAverageSnowfall();
        double precip        =GetAveragePrecip();
        double channel_stor  =GetTotalChannelStorage();
        double reservoir_stor=GetTotalReservoirStorage();
        double rivulet_stor  =GetTotalRivuletStorage();

        _STORAGE<<tt.model_time <<","<<thisdate<<","<<thishour; //instantaneous, so thishour rather than usehour used.

        if (t!=0){_STORAGE<<","<<precip-snowfall<<","<<snowfall;}//precip
        else     {_STORAGE<<",---,---";}
        _STORAGE<<","<<channel_stor<<","<<reservoir_stor<<","<<rivulet_stor;

        currentWater=0.0;
        for (i=0;i<GetNumStateVars();i++)
        {
          if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iCumPrecip))
          {
            S=GetAvgStateVar(i);
            if (!silent){cout<<"  |"<< setw(6)<<setiosflags(ios::fixed) << setprecision(2)<<S;}
            _STORAGE<<","<<FormatDouble(S);
            currentWater+=S;
          }
        }
        currentWater+=channel_stor+rivulet_stor+reservoir_stor;
        if(t==0){
          // \todo [fix]: this fixes a mass balance bug in reservoir simulations, but there is certainly a more proper way to do it
          // JRC: I think somehow this is being double counted in the delta V calculations in the first timestep
          for(int p=0;p<_nSubBasins;p++){
            if(_pSubBasins[p]->GetReservoir()!=NULL){
              currentWater+=_pSubBasins[p]->GetIntegratedReservoirInflow(Options.timestep)/2.0/_WatershedArea*MM_PER_METER/M2_PER_KM2;
              currentWater-=_pSubBasins[p]->GetIntegratedOutflow        (Options.timestep)/2.0/_WatershedArea*MM_PER_METER/M2_PER_KM2;
            }
            //currentWater-=_pSubBasins[p]->GetIntegratedSpecInflow(0,Options.timestep)/2.0/_WatershedArea*MM_PER_METER/M2_PER_KM2;
          }
        }

        _STORAGE<<","<<currentWater<<","<<_CumulInput<<","<<_CumulOutput<<","<<FormatDouble((currentWater-_initWater)+(_CumulOutput-_CumulInput));
        _STORAGE<<endl;
      }

      //Write hydrographs for gauged watersheds to Hydrographs.csv (ALWAYS DONE)
      //----------------------------------------------------------------
      if (Options.ave_hydrograph)
      {
        _HYDRO<<usetime<<","<<usedate<<","<<usehour;
        if(t!=0) { _HYDRO<<","<<GetAveragePrecip(); }//watershed-wide precip + irrigation
        else     { _HYDRO<<",---";                  }

        for (int p=0;p<_nSubBasins;p++)
        {
          pSB=_pSubBasins[p];
          if (pSB->IsGauged()  && (pSB->IsEnabled()))
          {
            _HYDRO<<","<<pSB->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);

            for (i = 0; i < _nObservedTS; i++)
            {
              if (IsContinuousFlowObs(_pObservedTS[i],pSB->GetID()))
              {
                double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep); //time shift handled in CTimeSeries::Parse
                if ((val != RAV_BLANK_DATA) && (tt.model_time>0)){ _HYDRO << "," << val; }
                else                                             { _HYDRO << ",";       }
              }
            }
            if (Options.write_localflow){
              _HYDRO<<","<<pSB->GetIntegratedLocalOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
            }
            if (pSB->GetReservoir() != NULL){
              _HYDRO<<","<<pSB->GetIntegratedReservoirInflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
                if (Options.write_netresinflow){
                  double Qin=pSB->GetIntegratedReservoirInflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
                  double P  =pSB->GetReservoir()->GetReservoirPrecipGains(Options.timestep)/(Options.timestep*SEC_PER_DAY);
                  double E  =pSB->GetReservoir()->GetReservoirEvapLosses (Options.timestep)/(Options.timestep*SEC_PER_DAY);
                  double GW =pSB->GetReservoir()->GetReservoirGWLosses   (Options.timestep)/(Options.timestep*SEC_PER_DAY);
                  _HYDRO<<","<<Qin+P-E-GW;
                }
              for(i = 0; i < _nObservedTS; i++)
              {
                if(IsContinuousInflowObs(_pObservedTS[i],pSB->GetID()))
                {
                  double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep); //time shift handled in CTimeSeries::Parse
                  if((val != RAV_BLANK_DATA) && (tt.model_time>0)) { _HYDRO << "," << val; }
                  else                                             { _HYDRO << ","; }
                }
              }
            }
          }
        }
        _HYDRO<<endl;
      }
      else //point value hydrograph or t==0
      {
        if((Options.period_starting) && (t==0)){}//don't write anything at time zero
        else{
          _HYDRO<<t<<","<<thisdate<<","<<thishour;
          if(t!=0){ _HYDRO<<","<<GetAveragePrecip(); }//watershed-wide precip+irrigation
          else    { _HYDRO<<",---";                  }
          for(int p=0;p<_nSubBasins;p++)
          {
            pSB=_pSubBasins[p];
            if(pSB->IsGauged()  && (pSB->IsEnabled()))
            {
              _HYDRO<<","<<pSB->GetOutflowRate();

              for(i = 0; i < _nObservedTS; i++){
                if(IsContinuousFlowObs(_pObservedTS[i],pSB->GetID()))
                {
                  double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep);
                  if((val != RAV_BLANK_DATA) && (tt.model_time>0)){ _HYDRO << "," << val; }
                  else                                            { _HYDRO << ","; }
                }
              }
              if (Options.write_localflow){
                _HYDRO<<","<<pSB->GetLocalOutflowRate();
              }
              if(pSB->GetReservoir() != NULL){
                _HYDRO<<","<<pSB->GetReservoirInflow();
                if (Options.write_netresinflow){
                  double Qin=pSB->GetReservoirInflow();
                  double P  =pSB->GetReservoir()->GetReservoirPrecipGains(Options.timestep)/(Options.timestep*SEC_PER_DAY);
                  double E  =pSB->GetReservoir()->GetReservoirEvapLosses (Options.timestep)/(Options.timestep*SEC_PER_DAY);
                  double GW =pSB->GetReservoir()->GetReservoirGWLosses   (Options.timestep)/(Options.timestep*SEC_PER_DAY);
                  _HYDRO<<","<<Qin+P-E-GW;
                }
                for(i = 0; i < _nObservedTS; i++)
                {
                  if(IsContinuousInflowObs(_pObservedTS[i],pSB->GetID()))
                  {
                    double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep); //time shift handled in CTimeSeries::Parse
                    if((val != RAV_BLANK_DATA) && (tt.model_time>0)) { _HYDRO << "," << val; }
                    else                                             { _HYDRO << ",";        }
                  }
                }
              }
            }
          }
          _HYDRO<<endl;
        }
      }
    }
    else if (Options.output_format==OUTPUT_ENSIM)
    {
      WriteEnsimMinorOutput(Options,tt);
    }
    else if (Options.output_format==OUTPUT_NETCDF)
    {
      WriteNetcdfMinorOutput(Options,tt);
    }

    //Write cumulative mass balance info to HRUGroup_MassEnergyBalance.csv
    //----------------------------------------------------------------
    if (Options.write_group_mb!=DOESNT_EXIST)
    {
      if((Options.period_starting) && (t==0)){}//don't write anything at time zero
      else{
        double sum;
        int kk=Options.write_group_mb;
        ofstream HGMB;
        tmpFilename=FilenamePrepare(_pHRUGroups[kk]->GetName()+"_MassEnergyBalance.csv",Options);
        HGMB.open(tmpFilename.c_str(),ios::app);
        HGMB<<usetime<<","<<usedate<<","<<usehour;
        double areasum=0.0;
        for(k = 0; k < _nHydroUnits; k++){
          if(_pHRUGroups[kk]->IsInGroup(k)){
            areasum+=_pHydroUnits[k]->GetArea();
          }
        }
        for(int js=0;js<_nTotalConnections;js++)
        {
          sum=0.0;
          for(k = 0; k < _nHydroUnits; k++){
            if(_pHRUGroups[kk]->IsInGroup(k)){
              sum += _aCumulativeBal[k][js] * _pHydroUnits[k]->GetArea();
            }
          }
          HGMB<<","<<sum/areasum;
        }
        HGMB<<endl;
        HGMB.close();
      }
    }

    //Write cumulative mass balance info to WatershedMassEnergyBalance.csv
    //----------------------------------------------------------------
    if (Options.write_mass_bal)
    {
      if((Options.period_starting) && (t==0)){}//don't write anything at time zero
      else{
        double sum;
        ofstream MB;
        tmpFilename=FilenamePrepare("WatershedMassEnergyBalance.csv",Options);
        MB.open(tmpFilename.c_str(),ios::app);

        MB<<usetime<<","<<usedate<<","<<usehour;
        for(int js=0;js<_nTotalConnections;js++)
        {
          sum=0.0;
          for(k=0;k<_nHydroUnits;k++){
            if(_pHydroUnits[k]->IsEnabled())
            {
              sum+=_aCumulativeBal[k][js]*_pHydroUnits[k]->GetArea();
            }
          }
          MB<<","<<sum/_WatershedArea;
        }
        MB<<endl;
        MB.close();
      }
    }

    //WaterLevels.csv
    //Write hydrographs for gauged watersheds (ALWAYS DONE)
    //----------------------------------------------------------------
    if (Options.write_waterlevels)
    {
      if((Options.period_starting) && (t==0)){}//don't write anything at time zero
      else{
        _LEVELS<<t<<","<<thisdate<<","<<thishour;
        if(t!=0){ _LEVELS<<","<<GetAveragePrecip(); }//watershed-wide precip+irrigation
        else    { _LEVELS<<",---";                  }
        for(int p=0;p<_nSubBasins;p++)
        {
          pSB=_pSubBasins[p];
          if(pSB->IsGauged()  && (pSB->IsEnabled()))
          {
            _LEVELS<<","<<pSB->GetWaterLevel();

            for(i = 0; i < _nObservedTS; i++){
              if(IsContinuousLevelObs(_pObservedTS[i],pSB->GetID()))
              {
                double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep);
                if((val != RAV_BLANK_DATA) && (tt.model_time>0)){ _LEVELS << "," << val; }
                else                                            { _LEVELS << ","; }
              }
            }
          }
        }
        _LEVELS<<endl;
      }
    }



    //ReservoirStages.csv
    //--------------------------------------------------------------
    if ((Options.write_reservoir) && (Options.output_format==OUTPUT_STANDARD))
    {
      int nn=(int)((tt.model_time+TIME_CORRECTION)/Options.timestep);//current timestep index

      if((Options.period_starting) && (t==0)){}//don't write anything at time zero
      else{
	    _RESSTAGE<< t<<","<<thisdate<<","<<thishour<<","<<GetAveragePrecip();
	    for (int p=0;p<_nSubBasins;p++)
        {
          pSB=_pSubBasins[p];
	        if ((pSB->IsGauged())  && (pSB->IsEnabled()) && (pSB->GetReservoir()!=NULL))
          {
	          _RESSTAGE<<","<<pSB->GetReservoir()->GetResStage();

	          for (i = 0; i < _nObservedTS; i++){
	            if (IsContinuousStageObs(_pObservedTS[i],pSB->GetID()))
	            {
                  double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep);
	              if ((val != RAV_BLANK_DATA) && (tt.model_time>0)){ _RESSTAGE << "," << val; }
	              else                                             { _RESSTAGE << ",";        }
	            }
            }
          }
        }
	      _RESSTAGE<<endl;
			}
    }

    // Demands.csv
    //----------------------------------------------------------------
    double irr,Qd,Qr,Q,eF;
    double demsum,delsum,retsum;
    if((Options.write_demandfile) && (Options.output_format!=OUTPUT_NONE))
    {
      if((Options.period_starting) && (t==0)) {}//don't write anything at time zero
      else {
        _DEMANDS<< t<<","<<thisdate<<","<<thishour;
        for(int p=0;p<_nSubBasins;p++)
        {
          pSB=_pSubBasins[p];
          if((pSB->IsEnabled()) && (pSB->IsGauged()) && (pSB->HasIrrigationDemand()))
          {
            Q   =pSB->GetOutflowRate     (); //AFTER irrigation removed -INSTANTANEOUS FLOW
            eF  =pSB->GetEnviroMinFlow   (tt.model_time);
            _DEMANDS<<","<<Q<<","<<eF;

            demsum=delsum=retsum=0.0;
            for (int ii=0;ii<pSB->GetNumWaterDemands();ii++)
            {
              demsum+=pSB->GetWaterDemand(ii);
              delsum+=pSB->GetDemandDelivery(ii);
              retsum+=pSB->GetReturnFlow(ii);
            }
            if (pSB->GetReservoir()!=NULL){
              for (int ii=0;ii<pSB->GetReservoir()->GetNumWaterDemands(); ii++){
                demsum =pSB->GetReservoir()->GetWaterDemand(ii);
                delsum =pSB->GetReservoir()->GetDemandDelivery(ii);
                retsum =pSB->GetReservoir()->GetReturnFlow(ii);
              }
            }
            _DEMANDS<<","<<demsum<<","<<delsum<<","<<retsum;

            for (int ii=0;ii<pSB->GetNumWaterDemands();ii++)
            {
              irr =pSB->GetWaterDemand(ii);
              Qd  =pSB->GetDemandDelivery(ii);
              Qr  =pSB->GetReturnFlow(ii);
              //double unmet=max(irr-Qd,0.0);
              _DEMANDS<<","<<irr<<","<<Qd;
              if (pSB->GetWaterDemandObj(ii)->HasReturnFlow()){
                _DEMANDS<<","<<Qr;
              }
            }
            if (pSB->GetReservoir()!=NULL){
              for (int ii=0;ii<pSB->GetReservoir()->GetNumWaterDemands(); ii++){
                irr =pSB->GetReservoir()->GetWaterDemand(ii);
                Qd  =pSB->GetReservoir()->GetDemandDelivery(ii);
                Qr  =pSB->GetReservoir()->GetReturnFlow(ii);
                //double unmet=max(irr-Qd,0.0);
                _DEMANDS<<","<<irr<<","<<Qd;
                if (pSB->GetReservoir()->GetWaterDemandObj(ii)->HasReturnFlow()){
                  _DEMANDS<<","<<Qr;
                }
              }
            }
          }
        }
        _DEMANDS<<endl;
      }
    }

    //ReservoirMassBalance.csv
    //----------------------------------------------------------------
    if((Options.write_reservoirMB) && (Options.output_format==OUTPUT_STANDARD))
    {
      if((Options.period_starting) && (t==0)){}//don't write anything at time zero
      else{
        ofstream RES_MB;
        tmpFilename=FilenamePrepare("ReservoirMassBalance.csv",Options);
        RES_MB.open(tmpFilename.c_str(),ios::app);
        if(RES_MB.fail()){
          ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
        }

        RES_MB<< usetime<<","<<usedate<<","<<usehour<<","<<GetAveragePrecip();
        double in,out,loss,stor,oldstor,precip,evap,seepage,stage,area;
        for(int p=0;p<_nSubBasins;p++)
        {
          pSB=_pSubBasins[p];
          if((pSB->IsGauged()) &&  (pSB->IsEnabled()) && (pSB->GetReservoir()!=NULL))
          {
            string constraint_str;

            stage         =pSB->GetReservoir()->GetResStage();//m
            area          =pSB->GetReservoir()->GetSurfaceArea();//m2
            in            =pSB->GetIntegratedReservoirInflow(Options.timestep);//m3
            out           =pSB->GetIntegratedOutflow        (Options.timestep);//m3
            RES_MB<<","<<stage<<","<<area;
            RES_MB<<","<<in<<","<<out;
            for (int i = 0; i < pSB->GetReservoir()->GetNumControlStructures(); i++) {
              double Qc=pSB->GetReservoir()->GetIntegratedControlOutflow(i,Options.timestep);
              RES_MB<<","<<Qc;//m3
              RES_MB<<","<<pSB->GetReservoir()->GetRegimeName(i,prev);//evaluated at start of timestep
              if (pSB->GetReservoir()->GetControlFlowTarget(i)!=pSB->GetDownstreamID()){out+=Qc;}
            }
            stor          =pSB->GetReservoir()->GetStorage             ();//m3
            oldstor       =pSB->GetReservoir()->GetOldStorage          ();//m3
            loss          =pSB->GetReservoir()->GetReservoirLosses     (Options.timestep);//m3 = GW+ET
            precip        =pSB->GetReservoir()->GetReservoirPrecipGains(Options.timestep);//m3
            evap          =pSB->GetReservoir()->GetReservoirEvapLosses (Options.timestep);//m3
            seepage       =pSB->GetReservoir()->GetReservoirGWLosses   (Options.timestep);//m3
            constraint_str=pSB->GetReservoir()->GetCurrentConstraint   ();
            if(tt.model_time==0.0){ in=0.0; }
            RES_MB<<","<<precip<<","<<evap<<","<<seepage;
            RES_MB<<","<<stor<<","<<loss<<","<<in-out-loss+precip-(stor-oldstor)<<","<<constraint_str;
          }
        }
        RES_MB<<endl;
        RES_MB.close();
      }
    }


    // ExhaustiveMassBalance.csv
    //--------------------------------------------------------------
    if (Options.write_exhaustiveMB)
    {
      if((Options.period_starting) && (t==0)){}//don't write anything at time zero
      else{
        int j,js,q;
        double cumsum;
        double sum;

        ofstream MB;
        tmpFilename=FilenamePrepare("ExhaustiveMassBalance.csv",Options);
        MB.open(tmpFilename.c_str(),ios::app);

        MB<<usetime<<","<<usedate<<","<<usehour;
        for(i=0;i<_nStateVars;i++)
        {
          if(CStateVariable::IsWaterStorage(_aStateVarType[i]))
          {
            cumsum=0.0;
            js=0;
            for(j=0;j<_nProcesses;j++){
              for(q=0;q<_pProcesses[j]->GetNumConnections();q++){
                if(_pProcesses[j]->GetFromIndices()[q]==i)
                {
                  sum=0.0;
                  for(k=0;k<_nHydroUnits;k++){
                    if(_pHydroUnits[k]->IsEnabled())
                    {
                      sum+=_aCumulativeBal[k][js]*_pHydroUnits[k]->GetArea();
                    }
                  }
                  MB<<","<<-sum/_WatershedArea;
                  cumsum-=sum/_WatershedArea;
                }
                js++;
              }
            }
            js=0;
            for(j=0;j<_nProcesses;j++){
              for(q=0;q<_pProcesses[j]->GetNumConnections();q++){
                if(_pProcesses[j]->GetToIndices()[q]==i)
                {
                  sum=0.0;
                  for(k=0;k<_nHydroUnits;k++){
                    if(_pHydroUnits[k]->IsEnabled())
                    {
                      sum+=_aCumulativeBal[k][js]*_pHydroUnits[k]->GetArea();
                    }
                  }
                  MB<<","<<sum/_WatershedArea;
                  cumsum+=sum/_WatershedArea;
                }
                js++;
              }
            }
            js=0;
            bool found;
            for(j=0;j<_nProcesses;j++){
              sum=0;
              found=false;
              for(q=0;q<_pProcesses[j]->GetNumLatConnections();q++){
                CLateralExchangeProcessABC *pProc=static_cast<CLateralExchangeProcessABC *>(_pProcesses[j]);
                if (pProc->GetLateralToIndices()[q]==i){
                  sum+=_aCumulativeLatBal[js];found=true;
                }
                if (pProc->GetLateralFromIndices()[q]==i){
                  sum-=_aCumulativeLatBal[js];found=true;
                }
                js++;
              }
              if((_pProcesses[j]->GetNumLatConnections()>0) && (found==true)){
                MB<<","<<sum/_WatershedArea;
                cumsum+=sum/_WatershedArea;
              }
            }

            //Cumulative, storage, error
            double Initial_i=0.0; //< \todo [bug] need to evaluate and store initial storage actross watershed!!!
            MB<<","<<cumsum<<","<<GetAvgStateVar(i)<<","<<cumsum-GetAvgStateVar(i)-Initial_i;
          }
        }
        MB<<endl;
        MB.close();
      }
    }

    // ForcingFunctions.csv
    //----------------------------------------------------------------
    if ((Options.write_forcings) && (Options.output_format==OUTPUT_STANDARD))
    {
      if(t==0){}//don't write anything at time zero
      else{
        force_struct *pFave;
        force_struct faveStruct = GetAverageForcings();
        pFave = &faveStruct;
        _FORCINGS<<usetime<<","<<usedate<<","<<usehour<<",";
        _FORCINGS<<pFave->day_angle<<",";
        _FORCINGS<<pFave->precip*(1-pFave->snow_frac) <<",";
        _FORCINGS<<pFave->precip*(pFave->snow_frac) <<",";
        _FORCINGS<<pFave->temp_ave<<",";
        _FORCINGS<<pFave->temp_daily_min<<",";
        _FORCINGS<<pFave->temp_daily_max<<",";
        _FORCINGS<<pFave->temp_daily_ave<<",";
        _FORCINGS<<pFave->temp_month_min<<",";
        _FORCINGS<<pFave->temp_month_max<<",";
        _FORCINGS<<pFave->air_dens<<",";
        _FORCINGS<<pFave->air_pres<<",";
        _FORCINGS<<pFave->rel_humidity<<",";
        _FORCINGS<<pFave->cloud_cover<<",";
        _FORCINGS<<pFave->ET_radia<<",";
        _FORCINGS<<pFave->SW_radia<<",";
        _FORCINGS<<pFave->SW_radia_net<<",";
        //_FORCINGS<<pFave->SW_radia_subcan<<",";
        //_FORCINGS<<pFave->SW_subcan_net<<",";
        _FORCINGS<<pFave->LW_incoming<<",";
        _FORCINGS<<pFave->LW_radia_net<<",";
        _FORCINGS<<pFave->wind_vel<<",";
        _FORCINGS<<pFave->PET<<",";
        _FORCINGS<<pFave->OW_PET<<",";
        _FORCINGS<<pFave->subdaily_corr<<",";
        _FORCINGS<<pFave->potential_melt;
        _FORCINGS<<endl;
      }
    }

    // Transport output files
    //--------------------------------------------------------------
    _pTransModel->WriteMinorOutput(Options,tt);

    if (Options.management_optimization) {
      _pDO->WriteMinorOutput(Options,tt);
    }

    // raven_debug.csv
    //--------------------------------------------------------------
    if (Options.debug_mode)
    {
      ofstream DEBUG;
      tmpFilename=FilenamePrepare("raven_debug.csv",Options);
      DEBUG.open(tmpFilename.c_str(),ios::app);
      DEBUG<<t<<","<<thisdate<<","<<thishour;
      for(i=0;i<10;i++){DEBUG<<","<<g_debug_vars[i];}
      DEBUG<<endl;
      DEBUG.close();
    }

    // HRU storage output
    //--------------------------------------------------------------
    if (_pOutputGroup!=NULL)
    {
      for (int kk=0;kk<_pOutputGroup->GetNumHRUs();kk++)
      {
        ofstream HRUSTOR;
        tmpFilename="HRUStorage_"+to_string(_pOutputGroup->GetHRU(kk)->GetHRUID())+".csv";
        tmpFilename=FilenamePrepare(tmpFilename,Options);
        HRUSTOR.open(tmpFilename.c_str(),ios::app);

        const force_struct *F=_pOutputGroup->GetHRU(kk)->GetForcingFunctions();

        HRUSTOR<<tt.model_time <<","<<thisdate<<","<<thishour;//instantaneous -no period starting correction

        if (t!=0){HRUSTOR<<","<<F->precip*(1-F->snow_frac)<<","<<F->precip*(F->snow_frac);}//precip
        else     {HRUSTOR<<",---,---";}

        currentWater=0;
        for (i=0;i<GetNumStateVars();i++)
        {
          if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iCumPrecip))
          {
            S=_pOutputGroup->GetHRU(kk)->GetStateVarValue(i);
            HRUSTOR<<","<<S;
            currentWater+=S;
          }
        }
        HRUSTOR<<","<<currentWater;
        HRUSTOR<<endl;
        HRUSTOR.close();
      }
    }

    // Assimilation adjustments
    //--------------------------------------------------------------
    if (Options.assimilate_flow){
      ofstream ASSIM;
      tmpFilename=FilenamePrepare("AssimilationAdjustments_BySubbasin.csv",Options);
      ASSIM.open(tmpFilename.c_str(),ios::app);

      ASSIM<<usetime<<","<<usedate;
      for (int p = 0; p < _nSubBasins; p++) {
        if((_pSubBasins[p]->IsGauged()) && (_pSubBasins[p]->IsEnabled())){
          if (Options.assim_method==DA_ECCC){ASSIM<<","<<_aDAQadjust[p]<<" ";}
          else                              {ASSIM<<","<<_aDAscale  [p]<<" ";}
        }
      }
      ASSIM<<endl;
      ASSIM.close();
    }
  } // end of write output interval if statement

  // Custom output files
  //--------------------------------------------------------------
  for (int c=0;c<_nCustomOutputs;c++)
  {
    _pCustomOutputs[c]->WriteCustomOutput(tt,Options);
  }

  // Write major output, if necessary
  //--------------------------------------------------------------
  if ((_nOutputTimes>0) && (_currOutputTimeInd<_nOutputTimes) && (tt.model_time>_aOutputTimes[_currOutputTimeInd]-0.5*Options.timestep))
  {
    string thishour=DecDaysToHours(tt.julian_day);
    _currOutputTimeInd++;
    tmpFilename="state_"+tt.date_string.substr(0,4)+tt.date_string.substr(5,2)+tt.date_string.substr(8,2)+"_"+thishour.substr(0,2)+thishour.substr(3,2);
    WriteMajorOutput(Options,tt,tmpFilename,false);
  }

}


//////////////////////////////////////////////////////////////////
/// \brief Replaces the WriteOutputFileHeaders function by not requiring the Options structure as an argument
/// \param &tt [in] Local (model) time *at the end of* the pertinent time step
/// \param solfile [in] Name of the solution file to be written
/// \param final [in] Whether this is the final solution file to be written
//
void CModel::WriteMajorOutput(const time_struct &tt, string solfile, bool final) const
{
  int i,k;
  string tmpFilename;
  const optStruct* Options = this->_pOptStruct;  // just to make the code more readable

  if (Options->output_format==OUTPUT_NONE){return;} //:SuppressOutput is on

  // WRITE {RunName}_solution.rvc - final state variables file
  ofstream RVC;
  tmpFilename=FilenamePrepare(solfile+".rvc", *_pOptStruct);
  RVC.open(tmpFilename.c_str());
  if (RVC.fail()){
    WriteWarning(("CModel::WriteMajorOutput: Unable to open output file "+tmpFilename+" for writing.").c_str(),
                  Options->noisy);
  }
  RVC<<":TimeStamp "<<tt.date_string<<" "<<DecDaysToHours(tt.julian_day)<<endl;

  //Header--------------------------
  //write in blocks of 80 state variables
  int mini,maxi;
  int M=80;
  for (int j=0; j<ceil(GetNumStateVars()/(double)(M)); j++){
    mini=j*M;
    maxi=min(GetNumStateVars(),(j+1)*M);
    RVC<<":HRUStateVariableTable"<<endl;
    RVC<<"  :Attributes,";
    for (i=mini;i<maxi;i++)
    {
      RVC << _pStateVar->SVTypeToString(_aStateVarType[i], _aStateVarLayer[i]);
      if (i!=GetNumStateVars()-1){RVC<<",";}
    }
    RVC<<endl;
    RVC<<"  :Units,";
    for (i=mini;i<maxi;i++)
    {
      RVC<<CStateVariable::GetStateVarUnits(_aStateVarType[i]);
      if (i!=GetNumStateVars()-1){RVC<<",";}
    }
    RVC<<endl;
    //Data----------------------------
    for (k=0;k<_nHydroUnits;k++)
    {
      RVC<<std::fixed; RVC.precision(5);
      RVC<<"  "<<_pHydroUnits[k]->GetHRUID()<<",";
      for (i=mini;i<maxi;i++)
      {
        RVC<<_pHydroUnits[k]->GetStateVarValue(i);
        if (i!=GetNumStateVars()-1){RVC<<",";}
      }
      RVC<<endl;
    }
    RVC<<":EndHRUStateVariableTable"<<endl;
  }
  //By basin------------------------
  RVC<<":BasinStateVariables"<<endl;
  for (int p=0;p<_nSubBasins;p++){
    RVC<<"  :BasinIndex "<<_pSubBasins[p]->GetID()<<",";
    _pSubBasins[p]->WriteToSolutionFile(RVC);
  }
  RVC<<":EndBasinStateVariables"<<endl;

  _pTransModel->WriteMajorOutput(RVC);

  RVC.close();

  // SubbasinProperties.csv
  //--------------------------------------------------------------
  if (Options->write_basinfile){
    ofstream BAS;
    bool has_obs;
    tmpFilename=FilenamePrepare("SubbasinProperties.csv",*_pOptStruct);
    BAS.open(tmpFilename.c_str());
    if(BAS.fail()) {
      WriteWarning(("CModel::WriteMinorOutput: Unable to open output file "+tmpFilename+" for writing.").c_str(),Options->noisy);
    }
    BAS<<"ID,Qref[m3/s],reach_length[m],area[km2],drainage_area[km2],t_conc[d],t_peak[d],gamma_sh,gamma_sc[1/d],celerity[m/s],diffusivity[m2/s],has observed flow,N,UH[0],UH[1],UH[2],..."<<endl;
    for(int pp=0;pp<_nSubBasins;pp++) {
      has_obs=false;
      BAS<<_pSubBasins[pp]->GetID()<<",  "<<_pSubBasins[pp]->GetReferenceFlow();
      BAS<<","<<_pSubBasins[pp]->GetReachLength();
      BAS<<","<<_pSubBasins[pp]->GetBasinArea();
      BAS<<","<<_pSubBasins[pp]->GetDrainageArea();
      BAS<<","<<_pSubBasins[pp]->GetBasinProperties("TIME_CONC");
      BAS<<","<<_pSubBasins[pp]->GetBasinProperties("TIME_TO_PEAK");
      BAS<<","<<_pSubBasins[pp]->GetBasinProperties("GAMMA_SHAPE");
      BAS<<","<<_pSubBasins[pp]->GetBasinProperties("GAMMA_SCALE");
      BAS<<","<<_pSubBasins[pp]->GetBasinProperties("CELERITY");
      BAS<<","<<_pSubBasins[pp]->GetBasinProperties("DIFFUSIVITY");
      //Has flow observations
      for (i = 0; i < _nObservedTS; i++){
        if (IsContinuousFlowObs(_pObservedTS[i],_pSubBasins[pp]->GetID())) {has_obs=true;}
      }
      if (has_obs){BAS<<", TRUE";}
      else        {BAS<<",FALSE";}

      BAS<<","<<_pSubBasins[pp]->GetLatHistorySize();
      for (int i = 0; i < _pSubBasins[pp]->GetLatHistorySize(); i++) {
        BAS<<","<<_pSubBasins[pp]->GetUnitHydrograph()[i];
      }
      BAS<<endl;
    }
    BAS.close();
  }

  // rating_curves.csv
  //--------------------------------------------------------------
  if(Options->write_channels){
    WriteRatingCurves(*Options);
  }
}


//////////////////////////////////////////////////////////////////
/// \brief Writes major output to file at the end of simulation
/// \details Writes:
/// - Solution file of all state variables; and
/// - Autogenerated parameters
///
/// \param &Options [in] Global model options information
//
void CModel::WriteMajorOutput(const optStruct &Options, const time_struct &tt, string solfile, bool final) const
{
  WriteMajorOutput(tt,solfile,final);
}
//////////////////////////////////////////////////////////////////
/// \brief Writes simple output to file
/// \note  written at start of time step before external script is read, after forcings processed.
/// \param &Options      [in] Global model options information
/// \param &tt           [in] time structure
//
void CModel::WriteSimpleOutput(const optStruct &Options, const time_struct &tt)
{
  if (Options.write_simpleout)
  {
    string thisdate=tt.date_string;                   //refers to date and time at START of time step
    string thishour=DecDaysToHours(tt.julian_day);
    double thistime=tt.model_time;

    force_struct *pFave;
    force_struct faveStruct = GetAverageForcings();
    pFave = &faveStruct;

    ofstream _SIMPLE;
    string tmpFilename=FilenamePrepare("simple_output.csv",Options);
    _SIMPLE.open(tmpFilename.c_str());
    if (_SIMPLE.fail()){
      ExitGracefully(("CModel::WriteSimpleOutput: unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }

    _SIMPLE<<thistime<<", "<<thisdate<<", "<<thishour<<", ";
    _SIMPLE<<pFave->precip<<", ";
    _SIMPLE<<pFave->temp_daily_min<<", ";
    _SIMPLE<<pFave->temp_daily_max<<", ";
    _SIMPLE<<pFave->SW_radia<<", ";
    _SIMPLE<<pFave->day_length<<", ";
    _SIMPLE<<endl;

    _SIMPLE.close();
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Writes progress file in JSON format (mainly for PAVICS runs)
///        Looks like:
///               {
///               	"% progress": 65,
///               	"seconds remaining": 123
///               }
///
/// \note  Does not account for initialization (reading) and final writing of model outputs. Only pure modeling time.
///
/// \param &Options      [in] Global model options information
/// \param &elapsed_time [in] elapsed time  (computational time markers)
/// \param elapsed_steps [in] elapsed number of simulation steps to perform (to determine % progress)
/// \param total_steps   [in]   total number of simulation steps to perform (to determine % progress)
//
void CModel::WriteProgressOutput(const optStruct &Options, clock_t elapsed_time, int elapsed_steps, int total_steps)
{
  if (Options.pavics)
  {
    ofstream PROGRESS;
    PROGRESS.open((Options.main_output_dir+"Raven_progress.txt").c_str());
    if (PROGRESS.fail()){
      PROGRESS.close();
      ExitGracefully("ParseInput:: Unable to open Raven_progress.txt. Bad output directory specified?",RUNTIME_ERR);
    }

    float total_time = (float(total_steps) * float(elapsed_time) / float(elapsed_steps)) / CLOCKS_PER_SEC;
    if (Options.benchmarking){ total_time =float(elapsed_time) / CLOCKS_PER_SEC;}
    PROGRESS<<"{"<<endl;
    PROGRESS<<"       \"% progress\": "       << int( float(elapsed_steps) * 100.0 / float(total_steps) ) <<","<< endl;
    PROGRESS<<"       \"seconds remaining\": "<< total_time - float(elapsed_time) / CLOCKS_PER_SEC <<endl;
    PROGRESS<<"}"<<endl;

    PROGRESS.close();
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Writes model summary information to screen
/// \param &Options [in] Global model options information
//
void CModel::SummarizeToScreen  (const optStruct &Options) const
{
  int rescount=0;
  for (int p = 0; p < _nSubBasins; p++){
    if (_pSubBasins[p]->GetReservoir() != NULL){rescount++;}
  }
  int disablecount=0;
  double allarea=0.0;
  for(int k=0;k<_nHydroUnits; k++){
    if(!_pHydroUnits[k]->IsEnabled()){disablecount++;}
    allarea+=_pHydroUnits[k]->GetArea();
  }
  int SBdisablecount=0;
  for(int p=0;p<_nSubBasins; p++){
    if(!_pSubBasins[p]->IsEnabled()){SBdisablecount++;}
  }
  time_struct tt,tt2;
  double day;
  int year;
  JulianConvert(0,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);
  AddTime(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.calendar,day,year);
  JulianConvert(0,day,year,Options.calendar,tt2);

  if(!Options.silent)
  {
    cout <<"==MODEL SUMMARY======================================="<<endl;
    cout <<"       Model Run: "<<Options.run_name    <<endl;
    cout <<"      Start time: "<<tt.date_string      <<endl;
    cout <<"        End time: "<<tt2.date_string     <<" (duration="<<Options.duration<<" days)"<<endl;
    cout <<"    rvi filename: "<<Options.rvi_filename<<endl;
    cout <<"Output Directory: "<<Options.main_output_dir  <<endl;
    cout <<"     # SubBasins: "<<GetNumSubBasins()   << " ("<< rescount << " reservoirs) ("<<SBdisablecount<<" disabled)"<<endl;
    cout <<"          # HRUs: "<<GetNumHRUs()        << " ("<<disablecount<<" disabled)"<<endl;
    cout <<"        # Gauges: "<<GetNumGauges()      <<endl;
    cout <<"#State Variables: "<<GetNumStateVars()   <<endl;
    for (int i=0;i<GetNumStateVars();i++){
      //don't write if convolution storage or advection storage?
      cout << "                - ";
      cout << CStateVariable::GetStateVarLongName(_aStateVarType[i],
                                                  _aStateVarLayer[i],
                                                  _pTransModel) << " (";
      cout << _pStateVar->SVTypeToString(_aStateVarType[i],_aStateVarLayer[i]) << ")" << endl;
    }
    cout <<"     # Processes: "<<GetNumProcesses()   <<endl;
    for (int j=0;j<GetNumProcesses();j++)
    {
      cout<<"                - ";
      cout<<GetProcessName(GetProcessType(j))<<endl;
    }
    cout <<"    #Connections: "<<_nTotalConnections          <<endl;
    cout <<"#Lat.Connections: "<<_nTotalLatConnections       <<endl;
    cout <<"        Duration: "<<Options.duration            <<" d"<<endl;
    cout <<"       Time step: "<<Options.timestep            <<" d ("<<(Options.timestep*MIN_PER_DAY)<<" min)"<<endl;
    cout <<"  Watershed Area: "<<_WatershedArea              <<" km2 (simulated) of "<<allarea<<" km2"<<endl;
    cout <<"======================================================"<<endl;
    cout <<endl;

  }
}

//////////////////////////////////////////////////////////////////
/// \brief run model diagnostics (at end of simulation)
///
/// \param &Options [in] global model options
//
void CModel::RunDiagnostics(const optStruct &Options)
{
  if((_nObservedTS==0) || (_nDiagnostics==0)) { return; }

  ofstream DIAG;
  string tmpFilename;
  tmpFilename=FilenamePrepare("Diagnostics.csv",Options);
  DIAG.open(tmpFilename.c_str());
  if(DIAG.fail()) {
    ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
  }
  //header
  DIAG<<"observed_data_series,filename,";
  for(int j=0; j<_nDiagnostics;j++) {
    DIAG<<_pDiagnostics[j]->GetName()<<",";
  }
  DIAG<<endl;
  //body
  bool skip;
  for (int d=0; d<_nDiagPeriods; d++){
    double starttime   = _pDiagPeriods[d]->GetStartTime();
    double endtime     = _pDiagPeriods[d]->GetEndTime();
    comparison compare = _pDiagPeriods[d]->GetComparison();
    double     thresh  = _pDiagPeriods[d]->GetThreshold();
    for(int i=0;i<_nObservedTS;i++)
    {
      skip=false;
      string datatype=_pObservedTS[i]->GetName();
      if((datatype=="HYDROGRAPH"          ) || (datatype=="RESERVOIR_STAGE")     || (datatype=="LAKE_AREA")   ||
         (datatype=="RESERVOIR_INFLOW"    ) || (datatype=="RESERVOIR_NETINFLOW") || (datatype=="WATER_LEVEL") ||
         (datatype=="STREAM_CONCENTRATION") || (datatype=="STREAM_TEMPERATURE"))
      {
        CSubBasin *pBasin=GetSubBasinByID(_pObservedTS[i]->GetLocID());
        if ((pBasin==NULL) || (!pBasin->IsEnabled())){skip=true;}
      }

      if (!skip)
      {

        DIAG<<_pObservedTS[i]->GetName()<<"_"<<_pDiagPeriods[d]->GetName()<<"["<<_pObservedTS[i]->GetLocID()<<"],"<<_pObservedTS[i]->GetSourceFile() <<",";//append to end of name for backward compatibility
        for(int j=0; j<_nDiagnostics;j++) {
          DIAG<<_pDiagnostics[j]->CalculateDiagnostic(_pModeledTS[i],_pObservedTS[i],_pObsWeightTS[i],starttime,endtime,compare,thresh,Options)<<",";
        }
        DIAG<<endl;
      }
    }
    //compute aggregate diagnostics
    for (int ii = 0; ii < _nAggDiagnostics; ii++)
    {
      int            kk  =_pAggDiagnostics[ii]->kk;
      string agg_datatype=_pAggDiagnostics[ii]->datatype;

      string name="";//"Average_HYDROGRAPH_GroupX"
      switch(_pAggDiagnostics[ii]->aggtype)
      {
        case AGG_AVERAGE:               name+="Average"; break;
        case AGG_MAXIMUM:               name+="Maximum"; break;
        case AGG_MINIMUM:               name+="Minimum"; break;
        case AGG_MEDIAN:                name+="Median";  break;
      }
      //name+=AggStatToString(_pAggDiagnostics[ii]->aggtype); //replace above with this
      name+="_"+agg_datatype;
      name+="_"+_pDiagPeriods[d]->GetName();
      if (kk!=DOESNT_EXIST){
        if (agg_datatype=="HYDROGRAPH"){name+="_"+_pSBGroups[kk]->GetName();} // \todo[funct]- handle more datatypes
      }

      DIAG<<name<<",[multiple],";
      for(int j=0; j<_nDiagnostics;j++) {
        DIAG<<CalculateAggDiagnostic(ii,j,starttime,endtime,compare,thresh,Options)<<",";
      }
      DIAG<<endl;
    }
  }

  DIAG.close();

  //reset for ensemble mode
  for(int i=0;i<_nObservedTS;i++)
  {
    _aObsIndex[i]=0;
  }

  // Post-simulation QA/QC reporting
  //----------------------------------------------------------------------------
  if (_pDO!=NULL){_pDO->Closure(Options); }

  for (int p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->GetReservoir()!=NULL){
      int N=_pSubBasins[p]->GetReservoir()->GetNumDryTimesteps();
      if (N>0)
      {
        string warn="CModel::RunDiagnostics: reservoir in basin "+to_string(_pSubBasins[p]->GetID())+ " went dry "+to_string(N) + " time steps during the course of the simulation";
        WriteWarning(warn,Options.noisy);
      }
    }
  }
}


//////////////////////////////////////////////////////////////////
/// \brief run model diagnostics (at end of simulation)
///
/// \param &Options [in] global model options
//
double CModel::CalculateAggDiagnostic(const int ii, const int j, const double &starttime, const double &endtime, const comparison compare,const double &thresh, const optStruct &Options)
{
  bool skip;
  double val;
  double stat=0;
  int N=0;
  agg_stat       type=_pAggDiagnostics[ii]->aggtype;
  int            kk  =_pAggDiagnostics[ii]->kk;
  string agg_datatype=_pAggDiagnostics[ii]->datatype;

  double *data=new double [_nObservedTS]; //maximum number of diagnostics that could be considered

  if     (type==AGG_AVERAGE) { stat=0; }
  else if(type==AGG_MAXIMUM) { stat=-ALMOST_INF; }
  else if(type==AGG_MINIMUM) { stat=ALMOST_INF; }
  else if(type==AGG_MEDIAN ) { stat=0; }
  for(int i=0;i<_nObservedTS;i++)
  {
    data[i]=0;
    skip=false;
    string datatype=_pObservedTS[i]->GetName();
    if (datatype!=agg_datatype){skip=true;}
    else if((datatype=="HYDROGRAPH"          ) || (datatype=="RESERVOIR_STAGE")     || (datatype=="LAKE_AREA")   ||
            (datatype=="RESERVOIR_INFLOW"    ) || (datatype=="RESERVOIR_NETINFLOW") || (datatype=="WATER_LEVEL") ||
            (datatype=="STREAM_CONCENTRATION") || (datatype=="STREAM_TEMPERATURE")) //subbasin-linked metrics
    {
      CSubBasin *pBasin=GetSubBasinByID(_pObservedTS[i]->GetLocID());
      if ((pBasin==NULL) || (!pBasin->IsEnabled())){skip=true;}
      if ((pBasin!=NULL) && (kk!=DOESNT_EXIST) && (!IsInSubBasinGroup(pBasin->GetID(),_pSBGroups[kk]->GetName()))){skip=true;}
    }
    else{ //HRU-linked
      int global_k=GetHRUByID(_pObservedTS[i]->GetLocID())->GetGlobalIndex();
      if ((kk!=DOESNT_EXIST) && (!IsInHRUGroup(global_k,_pHRUGroups[kk]->GetName()))){skip=true;}
    }

    if (!skip)
    {
      val=_pDiagnostics[j]->CalculateDiagnostic(_pModeledTS[i],_pObservedTS[i],_pObsWeightTS[i],starttime,endtime,compare, thresh,Options);

      if      (type==AGG_AVERAGE){stat+=val; N++;}
      else if (type==AGG_MAXIMUM){upperswap(stat,val);N++;}
      else if (type==AGG_MINIMUM){lowerswap(stat,val);N++;}
      else if (type==AGG_MEDIAN ){data[N]=val; N++;}
    }
  }
  int n;
  if(type==AGG_MEDIAN){quickSort(data,0,N-1);n=(int)rvn_floor((double)(N)/2.0+0.01);}

  if (N==0){return -ALMOST_INF;}
  if       (type==AGG_AVERAGE){stat/=N;}
  else if ((type==AGG_MEDIAN ) && (N%2==1)){stat=data[n];} //odd N
  else if ((type==AGG_MEDIAN ) && (N%2==0)){stat=0.5*(data[n]+data[n-1]);} //even N
  delete [] data;
  return stat;
}

//////////////////////////////////////////////////////////////////
/// \brief Writes output headers for WatershedStorage.tb0 and Hydrographs.tb0
///
/// \param &Options [in] global model options
//
void CModel::WriteEnsimStandardHeaders(const optStruct &Options)
{
  int i;
  time_struct tt, tt2;

  JulianConvert(0.0,              Options.julian_start_day, Options.julian_start_year, Options.calendar, tt);//start of the timestep
  JulianConvert(Options.timestep, Options.julian_start_day, Options.julian_start_year, Options.calendar, tt2);//end of the timestep

  //WatershedStorage.tb0
  //--------------------------------------------------------------
  int iCumPrecip = GetStateVarIndex(ATMOS_PRECIP);
  string tmpFilename;

  if (Options.write_watershed_storage)
  {
    tmpFilename = FilenamePrepare("WatershedStorage.tb0", Options);

    _STORAGE.open(tmpFilename.c_str());
    if (_STORAGE.fail()){
      ExitGracefully(("CModel::WriteEnsimStandardHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    _STORAGE << "#########################################################################" << endl;
    _STORAGE << ":FileType tb0 ASCII EnSim 1.0" << endl;
    _STORAGE << "#" << endl;
    _STORAGE << ":Application   Raven" << endl;
    if(!Options.benchmarking){
      _STORAGE << ":Version       " << Options.version << endl;
      _STORAGE << ":CreationDate  " << GetCurrentMachineTime() << endl;
    }
    _STORAGE << "#" << endl;
    _STORAGE << "#------------------------------------------------------------------------" << endl;
    _STORAGE << "#" << endl;
    _STORAGE << ":RunName       " << Options.run_name << endl;
    _STORAGE << ":Format         Instantaneous" << endl;
    _STORAGE << "#" << endl;

    if (Options.suppressICs){
      _STORAGE << ":StartTime " << tt2.date_string << " " << DecDaysToHours(tt2.julian_day) << endl;
    }
    else{
      _STORAGE << ":StartTime " << tt.date_string << " " << DecDaysToHours(tt.julian_day) << endl;
    }
    if (Options.timestep != 1.0){ _STORAGE << ":DeltaT " << DecDaysToHours(Options.timestep) << endl; }
    else                        { _STORAGE << ":DeltaT 24:00:00.00" << endl; }
    _STORAGE << "#" << endl;

    _STORAGE<<":ColumnMetaData"<<endl;
    _STORAGE<<"  :ColumnName rainfall snowfall \"Channel storage\" \"Rivulet storage\"";
    for (i=0;i<GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iCumPrecip)){
  	_STORAGE<<" \""<<CStateVariable::GetStateVarLongName(_aStateVarType[i],
                                                         _aStateVarLayer[i],
                                                         _pTransModel)<<"\"";}}
    _STORAGE<<" \"Total storage\" \"Cum. precip\" \"Cum. outflow\" \"MB error\""<<endl;

    _STORAGE<<"  :ColumnUnits mm/d mm/d mm mm ";
    for (i=0;i<GetNumStateVars();i++){
    if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iCumPrecip)){_STORAGE<<" mm";}}
    _STORAGE<<" mm mm mm mm"<<endl;

    _STORAGE<<"  :ColumnType float float float float";
    for (i=0;i<GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iCumPrecip)){_STORAGE<<" float";}}
    _STORAGE<<" float float float float"<<endl;

    _STORAGE << "  :ColumnFormat -1 -1 0 0";
    for (i = 0; i < GetNumStateVars(); i++){
      if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i != iCumPrecip)){
	    _STORAGE << " 0";
      }
    }
    _STORAGE << " 0 0 0 0" << endl;

    _STORAGE << ":EndColumnMetaData" << endl;

    _STORAGE << ":EndHeader" << endl;
  }

  //Hydrographs.tb0
  //--------------------------------------------------------------
  tmpFilename = FilenamePrepare("Hydrographs.tb0", Options);
  _HYDRO.open(tmpFilename.c_str());
  if (_HYDRO.fail()){
    ExitGracefully(("CModel::WriteEnsimStandardHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
  }
  _HYDRO << "#########################################################################" << endl;
  _HYDRO << ":FileType tb0 ASCII EnSim 1.0" << endl;
  _HYDRO << "#" << endl;
  _HYDRO << ":Application   Raven" << endl;
  if(!Options.benchmarking){
    _HYDRO << ":Version       " << Options.version << endl;
    _HYDRO << ":CreationDate  " << GetCurrentMachineTime() << endl;
  }
  _HYDRO << "#" << endl;
  _HYDRO << "#------------------------------------------------------------------------" << endl;
  _HYDRO << "#" << endl;
  _HYDRO << ":RunName       " << Options.run_name << endl;
  _HYDRO << "#" << endl;

  if (Options.ave_hydrograph){
    _HYDRO << ":Format         PeriodEnding" << endl;
  }
  else{
    _HYDRO << ":Format         Instantaneous" << endl;
  }
  if (((Options.period_ending) && (Options.ave_hydrograph)) || (Options.suppressICs )){
    _HYDRO << ":StartTime " << tt2.date_string << " " << DecDaysToHours(tt2.julian_day) << endl;
  }
  else{
    _HYDRO << ":StartTime " << tt.date_string << " " << DecDaysToHours(tt.julian_day) << endl;
  }

  if (Options.timestep!=1.0){_HYDRO<<":DeltaT " <<DecDaysToHours(Options.timestep)<<endl;}
  else                      {_HYDRO<<":DeltaT 24:00:00.00"  <<endl;}
  _HYDRO<<"#"<<endl;

  double val = 0; //snapshot hydrograph
  double val2=1;
  if (Options.ave_hydrograph){ val = 1; } //continuous hydrograph
  if (Options.period_ending) { val *= -1; val2*=-1;}//period ending

  _HYDRO<<":ColumnMetaData"<<endl;
  _HYDRO<<"  :ColumnName precip";
  for (int p=0;p<_nSubBasins;p++){_HYDRO<<" Q_"<<_pSubBasins[p]->GetID();}_HYDRO<<endl;
  _HYDRO<<"  :ColumnUnits mm/d";
  for (int p=0;p<_nSubBasins;p++){_HYDRO<<" m3/s";}_HYDRO<<endl;
  _HYDRO<<"  :ColumnType float";
  for (int p=0;p<_nSubBasins;p++){_HYDRO<<" float";}_HYDRO<<endl;
  _HYDRO<<"  :ColumnFormat "<<val2;
  for (int p=0;p<_nSubBasins;p++){_HYDRO<<" "<<val;}_HYDRO<<endl;
  _HYDRO<<":EndColumnMetaData"<<endl;
  _HYDRO<<":EndHeader"<<endl;
}

//////////////////////////////////////////////////////////////////
/// \brief Writes output headers for WatershedStorage.nc and Hydrographs.nc
///
/// \param &Options [in] global model options
/// original code developed by J. Mai
//
void CModel::WriteNetcdfStandardHeaders(const optStruct &Options)
{

#ifdef _RVNETCDF_

  time_struct tt;                                    // start time structure
  const int   ndims1 = 1;
  const int   ndims2 = 2;
  int         dimids1[ndims1];                       // array which will contain all dimension ids for a variable
  int         ncid, varid_pre;                       // When we create netCDF variables and dimensions, we get back an ID for each one.
  int         time_dimid, varid_time;                // dimension ID (holds number of time steps) and variable ID (holds time values) for time
  int         nSim, nbasins_dimid, varid_bsim,varid_bsim2;       
  //                                                 // # of sub-basins with simulated outflow, dimension ID, and
  //                                                 // variable to write basin IDs for simulated outflows
  int         varid_qsim;                            // variable ID for simulated outflows
  int         varid_qobs;                            // variable ID for observed outflows
  int         varid_qin;                             // variable ID for observed inflows
  int         varid_qloc;
  int         varid_qnet_in;                         // variable ID for reservoir net inflows

  int         retval;                                // error value for NetCDF routines
  string      tmpFilename;
  int         p;                                     // loop over all sub-basins
  string      tmp,tmp2,tmp3;

  // initialize all potential file IDs with -9 == "not existing and hence not opened"
  _HYDRO_ncid    = -9;   // output file ID for Hydrographs.nc         (-9 --> not opened)
  _STORAGE_ncid  = -9;   // output file ID for WatershedStorage.nc    (-9 --> not opened)
  _FORCINGS_ncid = -9;   // output file ID for ForcingFunctions.nc    (-9 --> not opened)
  _RESSTAGE_ncid = -9;   // output file ID for ReservoirStages.nc    (-9 --> not opened)
  _RESMB_ncid    = -9;

  //converts start day into "hours since YYYY-MM-DD HH:MM:SS"  (model start time)
  char  starttime[200]; // start time string in format 'hours since YYY-MM-DD HH:MM:SS'
  JulianConvert( 0.0,Options.julian_start_day, Options.julian_start_year, Options.calendar, tt);
  strcpy(starttime, "hours since ") ;
  strcat(starttime, tt.date_string.c_str()) ;
  strcat(starttime, " ");
  strcat(starttime, DecDaysToHours(tt.julian_day,true).c_str());
  if(Options.time_zone!=0) { strcat(starttime,TimeZoneToString(Options.time_zone).c_str()); }

  //====================================================================
  //  Hydrographs.nc
  //====================================================================
  // Create the file.
  tmpFilename = FilenamePrepare("Hydrographs.nc", Options);
  retval = nc_create(tmpFilename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);  HandleNetCDFErrors(retval);
  _HYDRO_ncid = ncid;

  // ----------------------------------------------------------
  // global attributes
  // ----------------------------------------------------------
  WriteNetCDFGlobalAttributes(_HYDRO_ncid,Options,"Standard Output");

  // ----------------------------------------------------------
  // time
  // ----------------------------------------------------------
  // (a) Define the DIMENSIONS. NetCDF will hand back an ID
  retval = nc_def_dim(_HYDRO_ncid, "time", NC_UNLIMITED, &time_dimid);  HandleNetCDFErrors(retval);

  /// Define the time variable. Assign units attributes to the netCDF VARIABLES.
  dimids1[0] = time_dimid;
  retval = nc_def_var(_HYDRO_ncid, "time", NC_DOUBLE, ndims1,dimids1, &varid_time); HandleNetCDFErrors(retval);
  retval = nc_put_att_text(_HYDRO_ncid, varid_time, "units"   ,      strlen(starttime)  , starttime);   HandleNetCDFErrors(retval);
  retval = nc_put_att_text(_HYDRO_ncid, varid_time, "calendar",      strlen("gregorian"), "gregorian"); HandleNetCDFErrors(retval);
  retval = nc_put_att_text(_HYDRO_ncid, varid_time, "standard_name", strlen("time"),      "time");      HandleNetCDFErrors(retval);

  // define precipitation variable
  varid_pre= NetCDFAddMetadata(_HYDRO_ncid, time_dimid,"precip","Precipitation","mm d**-1");

  // ----------------------------------------------------------
  // simulated/observed outflows
  // ----------------------------------------------------------
  // (a) count number of simulated outflows "nSim"
  nSim = 0;
  for (p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->IsGauged()  && (_pSubBasins[p]->IsEnabled())){nSim++;}
  }

  if (nSim > 0)
  {
    // (b) create dimension "nbasins"
    retval = nc_def_dim(_HYDRO_ncid, "nbasins", nSim, &nbasins_dimid);                             HandleNetCDFErrors(retval);

    // (c) create variable  and set attributes for"basin_name"
    dimids1[0] = nbasins_dimid;
    retval = nc_def_var(_HYDRO_ncid, "basin_name", NC_STRING, ndims1, dimids1, &varid_bsim);       HandleNetCDFErrors(retval);
    tmp ="ID of sub-basins with simulated outflows";
    tmp2="timeseries_id";
    tmp3="1";
    retval = nc_put_att_text(_HYDRO_ncid, varid_bsim, "long_name",  tmp.length(), tmp.c_str());    HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_HYDRO_ncid, varid_bsim, "cf_role"  , tmp2.length(),tmp2.c_str());    HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_HYDRO_ncid, varid_bsim, "units"    , tmp3.length(),tmp3.c_str());    HandleNetCDFErrors(retval);
    
    // (d) create variable  and set attributes for"basin_fullname"
    dimids1[0] = nbasins_dimid;
    retval = nc_def_var(_HYDRO_ncid, "basin_fullname", NC_STRING, ndims1, dimids1, &varid_bsim2);       HandleNetCDFErrors(retval);
    tmp ="Name of sub-basins with simulated outflows";
    tmp2="timeseries_id";
    tmp3="1";
    retval = nc_put_att_text(_HYDRO_ncid, varid_bsim2, "long_name",  tmp.length(), tmp.c_str());    HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_HYDRO_ncid, varid_bsim2, "cf_role"  , tmp2.length(),tmp2.c_str());    HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_HYDRO_ncid, varid_bsim2, "units"    , tmp3.length(),tmp3.c_str());    HandleNetCDFErrors(retval);

    varid_qsim= NetCDFAddMetadata2D(_HYDRO_ncid, time_dimid,nbasins_dimid,"q_sim","Simulated outflows","m**3 s**-1");
    varid_qobs= NetCDFAddMetadata2D(_HYDRO_ncid, time_dimid,nbasins_dimid,"q_obs","Observed outflows" ,"m**3 s**-1");
    varid_qin = NetCDFAddMetadata2D(_HYDRO_ncid, time_dimid,nbasins_dimid,"q_in" ,"Simulated reservoir inflows"  ,"m**3 s**-1");
    if (Options.write_localflow){
    varid_qloc= NetCDFAddMetadata2D(_HYDRO_ncid, time_dimid,nbasins_dimid,"q_loc" ,"Local inflow contribution"  ,"m**3 s**-1");
    }
    if (Options.write_netresinflow){
    varid_qnet_in=NetCDFAddMetadata2D(_HYDRO_ncid, time_dimid,nbasins_dimid,"qnet_in" ,"Simulated reservoir net inflows"  ,"m**3 s**-1");
    }
  }// end if nSim>0

  // End define mode. This tells netCDF we are done defining metadata.
  retval = nc_enddef(_HYDRO_ncid);  HandleNetCDFErrors(retval);

  // (a) write gauged basin names/IDs to variable "basin_name"
  if (nSim>0){
    WriteNetCDFBasinList(_HYDRO_ncid,varid_bsim,varid_bsim2,this,false,Options);
  }

  //====================================================================
  //  ReservoirStages.nc
  //====================================================================
  if(Options.write_reservoir)
  {
    // Create the file.
    tmpFilename = FilenamePrepare("ReservoirStages.nc",Options);
    retval = nc_create(tmpFilename.c_str(),NC_CLOBBER|NC_NETCDF4,&ncid);  HandleNetCDFErrors(retval);
    _RESSTAGE_ncid = ncid;

    // ----------------------------------------------------------
    // global attributes
    // ----------------------------------------------------------
    WriteNetCDFGlobalAttributes(_RESSTAGE_ncid,Options,"Standard Output");

    // ----------------------------------------------------------
    // time
    // ----------------------------------------------------------
    // (a) Define the DIMENSIONS. NetCDF will hand back an ID
    retval = nc_def_dim(_RESSTAGE_ncid,"time",NC_UNLIMITED,&time_dimid);  HandleNetCDFErrors(retval);

    /// Define the time variable. Assign units attributes to the netCDF VARIABLES.
    dimids1[0] = time_dimid;
    retval = nc_def_var     (_RESSTAGE_ncid,"time",NC_DOUBLE,ndims1,dimids1,&varid_time);           HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_RESSTAGE_ncid,varid_time,"units",strlen(starttime),starttime);        HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_RESSTAGE_ncid,varid_time,"calendar",strlen("gregorian"),"gregorian"); HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_RESSTAGE_ncid,varid_time,"standard_name",strlen("time"),"time");      HandleNetCDFErrors(retval);

    // define precipitation variable
    varid_pre= NetCDFAddMetadata(_RESSTAGE_ncid,time_dimid,"precip","Precipitation","mm d**-1");

    // ----------------------------------------------------------
    // simulated/observed stages
    // ----------------------------------------------------------
    // (a) count number of simulated reservoir stages "nSim"
    nSim = 0;
    for(p=0;p<_nSubBasins;p++) {
      if((_pSubBasins[p]->IsGauged())  && (_pSubBasins[p]->IsEnabled())  && (_pSubBasins[p]->GetReservoir()!=NULL)) { nSim++; }
    }

    if(nSim > 0)
    {
      // (b) create dimension "nbasins"
      retval = nc_def_dim(_RESSTAGE_ncid,"nbasins",nSim,&nbasins_dimid);                             HandleNetCDFErrors(retval);

      // (c) create variable  and set attributes for"basin_name"
      dimids1[0] = nbasins_dimid;
      retval = nc_def_var(_RESSTAGE_ncid,"basin_name",NC_STRING,ndims1,dimids1,&varid_bsim);         HandleNetCDFErrors(retval);
      tmp ="subbasin ID of reservoirs with simulated stage";
      tmp2="timeseries_id";
      tmp3="1";
      retval = nc_put_att_text(_RESSTAGE_ncid,varid_bsim,"long_name", tmp.length(), tmp.c_str());    HandleNetCDFErrors(retval);
      retval = nc_put_att_text(_RESSTAGE_ncid,varid_bsim,"cf_role",  tmp2.length(),tmp2.c_str());    HandleNetCDFErrors(retval);
      retval = nc_put_att_text(_RESSTAGE_ncid,varid_bsim,"units",    tmp3.length(),tmp3.c_str());    HandleNetCDFErrors(retval);

      // (d) create variable  and set attributes for"basin_fullname"
      dimids1[0] = nbasins_dimid;
      retval = nc_def_var(_RESSTAGE_ncid, "basin_fullname", NC_STRING, ndims1, dimids1, &varid_bsim2);       HandleNetCDFErrors(retval);
      tmp ="subbasin name of reservoirs with simulated stage";
      tmp2="timeseries_id";
      tmp3="1";
      retval = nc_put_att_text(_RESSTAGE_ncid, varid_bsim2, "long_name",  tmp.length(), tmp.c_str());    HandleNetCDFErrors(retval);
      retval = nc_put_att_text(_RESSTAGE_ncid, varid_bsim2, "cf_role"  , tmp2.length(),tmp2.c_str());    HandleNetCDFErrors(retval);
      retval = nc_put_att_text(_RESSTAGE_ncid, varid_bsim2, "units"    , tmp3.length(),tmp3.c_str());    HandleNetCDFErrors(retval);


      varid_qsim= NetCDFAddMetadata2D(_RESSTAGE_ncid,time_dimid,nbasins_dimid,"h_sim","Simulated stage","m");
      varid_qobs= NetCDFAddMetadata2D(_RESSTAGE_ncid,time_dimid,nbasins_dimid,"h_obs","Observed stage","m");

    }// end if nSim>0

    // End define mode. This tells netCDF we are done defining metadata.
    retval = nc_enddef(_RESSTAGE_ncid);  HandleNetCDFErrors(retval);

    // write values to NetCDF
    // (a) write gauged reservoir basin names/IDs to variable "basin_name"
    if (nSim>0){
      WriteNetCDFBasinList(_RESSTAGE_ncid,varid_bsim,varid_bsim2,this,true,Options);
    }
  }

  //====================================================================
  //  WatershedStorage.nc
  //====================================================================
  if (Options.write_watershed_storage)
  {
    tmpFilename = FilenamePrepare("WatershedStorage.nc", Options);
    retval = nc_create(tmpFilename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);  HandleNetCDFErrors(retval);
    _STORAGE_ncid = ncid;

    // ----------------------------------------------------------
    // global attributes
    // ----------------------------------------------------------
    WriteNetCDFGlobalAttributes(_STORAGE_ncid,Options,"Standard Output");

    // ----------------------------------------------------------
    // time vector
    // ----------------------------------------------------------
    // Define the DIMENSIONS. NetCDF will hand back an ID
    retval = nc_def_dim(_STORAGE_ncid, "time", NC_UNLIMITED, &time_dimid);  HandleNetCDFErrors(retval);

    /// Define the time variable.
    dimids1[0] = time_dimid;
    retval = nc_def_var(_STORAGE_ncid, "time", NC_DOUBLE, ndims1,dimids1, &varid_time); HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_STORAGE_ncid, varid_time, "units"   , strlen(starttime)  , starttime);   HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_STORAGE_ncid, varid_time, "calendar", strlen("gregorian"), "gregorian"); HandleNetCDFErrors(retval);

    // ----------------------------------------------------------
    // precipitation / channel storage / state vars / MB diagnostics
    // ----------------------------------------------------------
    int varid;
    varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"rainfall","rainfall","mm d**-1");
    varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"snowfall","snowfall","mm d**-1");
    varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"channel_storage","Channel Storage","mm");
    varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"reservoir_storage","Reservoir Storage","mm");
    varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"rivulet_storage","Rivulet Storage","mm");

    int iCumPrecip=GetStateVarIndex(ATMOS_PRECIP);
    for(int i=0;i<_nStateVars;i++){
      if((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iCumPrecip)){
	      string name = GetStateVarInfo()->GetStateVarLongName(_aStateVarType[i],
                                                              _aStateVarLayer[i],
                                                              GetTransportModel());
	      varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,name,name,"mm");
      }
    }
    varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"total","total water storage","mm");
    varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"cum_input","cumulative water input","mm");
    varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"cum_outflow","cumulative water output","mm");
    varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"MB_error","mass balance error","mm");

    // End define mode. This tells netCDF we are done defining metadata.
    retval = nc_enddef(_STORAGE_ncid);  HandleNetCDFErrors(retval);
  }

  //====================================================================
  //  ForcingFunctions.nc
  //====================================================================
  if(Options.write_forcings)
  {
    tmpFilename = FilenamePrepare("ForcingFunctions.nc",Options);
    retval = nc_create(tmpFilename.c_str(),NC_CLOBBER|NC_NETCDF4,&ncid);  HandleNetCDFErrors(retval);
    _FORCINGS_ncid = ncid;

    WriteNetCDFGlobalAttributes(_FORCINGS_ncid,Options,"Standard Output");

    // ----------------------------------------------------------
    // time vector
    // ----------------------------------------------------------
    retval = nc_def_dim(_FORCINGS_ncid,"time",NC_UNLIMITED,&time_dimid);  HandleNetCDFErrors(retval);
    dimids1[0] = time_dimid;
    retval = nc_def_var(_FORCINGS_ncid,"time",NC_DOUBLE,ndims1,dimids1,&varid_time); HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_FORCINGS_ncid,varid_time,"units",strlen(starttime),starttime);   HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_FORCINGS_ncid,varid_time,"calendar",strlen("gregorian"),"gregorian"); HandleNetCDFErrors(retval);

    // ----------------------------------------------------------
    int varid;
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"rainfall","rainfall","mm d**-1");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"snowfall","snowfall","mm d**-1");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"temp","temp","C");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"temp_daily_min","temp_daily_min","C");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"temp_daily_max","temp_daily_max","C");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"temp_daily_ave","temp_daily_ave","C");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"air_density","air density","kg m**-3");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"air_pressure","air pressure","kPa");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"relative_humidity","relative humidity","");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"cloud_cover","cloud cover","");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"ET_radiation","ET radiation","MJ m**-2 d**-1");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"SW_radiation","SW radiation","MJ m**-2 d**-1");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"net_SW_radiation","net SW radiation","MJ m**-2 d**-1");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"LW_incoming","LW incoming","MJ m**-2 d**-1");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"LW_radiation","LW radiation","MJ m**-2 d**-1");
    //varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"SW_radia_subcan","SW subcanopy","MJ m**-2 d**-1");
    //varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"SW_subcan_net","net SW subcanopy","MJ m**-2 d**-1");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"wind_velocity","wind velocity","m s**-1");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"PET","PET","mm d**-1");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"OW_PET","OW PET","mm d**-1");
    varid= NetCDFAddMetadata(_FORCINGS_ncid,time_dimid,"potential_melt","potential melt","mm d**-1");

    // End define mode. This tells netCDF we are done defining metadata.
    retval = nc_enddef(_FORCINGS_ncid);  HandleNetCDFErrors(retval);
  }

  //====================================================================
  //  ReservoirMassBalance.nc
  //====================================================================
  if(Options.write_reservoirMB)
  {
    tmpFilename = FilenamePrepare("ReservoirMassBalance.nc",Options);
    retval = nc_create(tmpFilename.c_str(),NC_CLOBBER|NC_NETCDF4,&ncid);  HandleNetCDFErrors(retval);
    _RESMB_ncid    = ncid;

    WriteNetCDFGlobalAttributes(_RESMB_ncid,Options,"Standard Output");

    // ----------------------------------------------------------
    // time vector
    // ----------------------------------------------------------
    retval = nc_def_dim(_RESMB_ncid,"time",NC_UNLIMITED,&time_dimid);  HandleNetCDFErrors(retval);
    dimids1[0] = time_dimid;
    retval = nc_def_var(_RESMB_ncid,"time",NC_DOUBLE,ndims1,dimids1,&varid_time);                HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_RESMB_ncid,varid_time,"units",strlen(starttime),starttime);        HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_RESMB_ncid,varid_time,"calendar",strlen("gregorian"),"gregorian"); HandleNetCDFErrors(retval);

    // ----------------------------------------------------------
    int varid;
    varid= NetCDFAddMetadata(_RESMB_ncid,time_dimid,"precip","Precipitation","mm d**-1");

    // (a) count number of simulated reservoirs "nSim"
    nSim = 0;
    for(p=0;p<_nSubBasins;p++) {
      if((_pSubBasins[p]->IsGauged())  && (_pSubBasins[p]->IsEnabled())  && (_pSubBasins[p]->GetReservoir()!=NULL)) { nSim++; }
    }

    if(nSim > 0)
    {
      // (b) create dimension "nbasins"
      retval = nc_def_dim(_RESMB_ncid,"nbasins",nSim,&nbasins_dimid);                             HandleNetCDFErrors(retval);

      // (c) create variable  and set attributes for"basin_name"
      dimids1[0] = nbasins_dimid;
      retval = nc_def_var(_RESMB_ncid,"basin_name",NC_STRING,ndims1,dimids1,&varid_bsim);         HandleNetCDFErrors(retval);
      tmp ="subbasin ID of reservoirs with simulated stage";
      tmp2="timeseries_id";
      tmp3="1";
      retval = nc_put_att_text(_RESMB_ncid,varid_bsim,"long_name", tmp.length(), tmp.c_str());    HandleNetCDFErrors(retval);
      retval = nc_put_att_text(_RESMB_ncid,varid_bsim,"cf_role",  tmp2.length(),tmp2.c_str());    HandleNetCDFErrors(retval);
      retval = nc_put_att_text(_RESMB_ncid,varid_bsim,"units",    tmp3.length(),tmp3.c_str());    HandleNetCDFErrors(retval);

      // (d) create variable  and set attributes for"basin_fullname"
      dimids1[0] = nbasins_dimid;
      retval = nc_def_var(_RESMB_ncid, "basin_fullname", NC_STRING, ndims1, dimids1, &varid_bsim2);       HandleNetCDFErrors(retval);
      tmp ="subbasin name of reservoirs with simulated stage";
      tmp2="timeseries_id";
      tmp3="1";
      retval = nc_put_att_text(_RESMB_ncid, varid_bsim2, "long_name",  tmp.length(), tmp.c_str());    HandleNetCDFErrors(retval);
      retval = nc_put_att_text(_RESMB_ncid, varid_bsim2, "cf_role"  , tmp2.length(),tmp2.c_str());    HandleNetCDFErrors(retval);
      retval = nc_put_att_text(_RESMB_ncid, varid_bsim2, "units"    , tmp3.length(),tmp3.c_str());    HandleNetCDFErrors(retval);

      varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"stage",    "stage",  "m");
      varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"area",     "area",  "m2");
      varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"inflow",   "inflow", "m3");
      varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"outflow",  "outflow","m3");
      varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"precip_m3","precip", "m3");
      varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"evap",     "evap",   "m3");
      varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"seepage",  "seepage","m3");
      varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"volume",   "volume", "m3");
      varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"losses",   "losses", "m3");
      varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"MB_error","Mass balance error","m3");
      //varid= NetCDFAddMetadata2D(_RESMB_ncid,time_dimid,nbasins_dimid,"constraint","constraint","");

      //NetCDF DOES NOT Report control structure outflow details (see .csv)

    }// end if nSim>0

    // End define mode. This tells netCDF we are done defining metadata.
    retval = nc_enddef(_RESMB_ncid);  HandleNetCDFErrors(retval);

    // write values to NetCDF
    // (a) write gauged reservoir basin names/IDs to variable "basin_name"
    if (nSim>0){
      WriteNetCDFBasinList(_RESMB_ncid,varid_bsim,varid_bsim2,this,true,Options);
    }
  }
#endif   // end compilation if NetCDF library is available
}


//////////////////////////////////////////////////////////////////
/// \brief Writes minor output data to WatershedStorage.tb0 and Hydrographs.tb0
///
/// \param &Options [in] global model options
/// \param &tt [in] current time structure
/// \todo [reorg] merge with WriteMinorOutput - too much complex code repetition here when only difference is (1) delimeter and (2) timestep info included in the .csv file
//
void  CModel::WriteEnsimMinorOutput (const optStruct &Options,
                                     const time_struct &tt)
{
  double currentWater;
  double S;
  int i;
  int iCumPrecip=GetStateVarIndex(ATMOS_PRECIP);

  double snowfall      =GetAverageSnowfall();
  double precip        =GetAveragePrecip();
  double channel_stor  =GetTotalChannelStorage();
  double reservoir_stor=GetTotalReservoirStorage();
  double rivulet_stor  =GetTotalRivuletStorage();

  if ((tt.model_time==0) && (Options.suppressICs==true) && (Options.period_ending)){return;}

  if ((_pEnsemble != NULL) && (_pEnsemble->DontWriteOutput())){return;}

  //----------------------------------------------------------------
  // write watershed state variables  (WatershedStorage.tb0)
  if (Options.write_watershed_storage)
  {
    if (tt.model_time!=0){_STORAGE<<" "<<precip-snowfall<<" "<<snowfall;}//precip
    else                 {_STORAGE<<" 0.0 0.0";}
    _STORAGE<<" "<<channel_stor+reservoir_stor<<" "<<rivulet_stor;

    currentWater=0.0;
    for (i=0;i<GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) &&  (i!=iCumPrecip)){
	      S=GetAvgStateVar(i);_STORAGE<<" "<<FormatDouble(S);currentWater+=S;
      }
    }
    currentWater+=channel_stor+rivulet_stor;
    if(tt.model_time==0){
      // \todo [fix]: this fixes a mass balance bug in reservoir simulations, but there is certainly a more proper way to do it
      // JRC: I think somehow this is being double counted in the delta V calculations in the first timestep
      for(int p=0;p<_nSubBasins;p++){
	      if(_pSubBasins[p]->GetReservoir()!=NULL){
	        currentWater+=_pSubBasins[p]->GetIntegratedReservoirInflow(Options.timestep)/2.0/_WatershedArea*MM_PER_METER/M2_PER_KM2;
	        currentWater-=_pSubBasins[p]->GetIntegratedOutflow        (Options.timestep)/2.0/_WatershedArea*MM_PER_METER/M2_PER_KM2;
	      }
      }
    }
    _STORAGE<<" "<<currentWater<<" "<<_CumulInput<<" "<<_CumulOutput<<" "<<FormatDouble((currentWater-_initWater)+(_CumulOutput-_CumulInput));
    _STORAGE<<endl;
  }

  //----------------------------------------------------------------
  //Write hydrographs for gauged watersheds (ALWAYS DONE) (Hydrographs.tb0)
  if (Options.ave_hydrograph)
  {
    if(tt.model_time!=0) { _HYDRO<<" "<<GetAveragePrecip(); }
    else                 { _HYDRO<<" 0.0"; }
    for (int p=0;p<_nSubBasins;p++){
      if (_pSubBasins[p]->IsGauged()  && (_pSubBasins[p]->IsEnabled()))
      {
        _HYDRO<<" "<<_pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
      }
    }
    _HYDRO<<endl;
  }
  else
  {
    if (tt.model_time!=0){_HYDRO<<" "<<GetAveragePrecip();}
    else                 {_HYDRO<<" 0.0";}
    for (int p=0;p<_nSubBasins;p++){
      if (_pSubBasins[p]->IsGauged() && (_pSubBasins[p]->IsEnabled()))
      {
        _HYDRO<<" "<<_pSubBasins[p]->GetOutflowRate();
      }
    }
    _HYDRO<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Writes minor output data to WatershedStorage.nc and Hydrographs.nc
///
/// \param &Options [in] global model options
/// \param &tt      [in] current time structure
//
void  CModel::WriteNetcdfMinorOutput ( const optStruct   &Options,
                                       const time_struct &tt)
{
#ifdef _RVNETCDF_

  int    retval;                // error value for NetCDF routines
  int    time_id;               // variable id in NetCDF for time
  size_t time_index[1], count1[1];  // determines where and how much will be written to NetCDF; 1D variable (pre, time)
  size_t start2[2], count2[2];  // determines where and how much will be written to NetCDF; 2D variable (qsim, qobs, qin)
  double current_time[1];       // current time in hours since start time
  double current_prec[1];       // precipitation of current time step
  size_t time_ind2;
  current_time[0] = tt.model_time*HR_PER_DAY;
  current_time[0]=RoundToNearestMinute(current_time[0]);

  time_index [0] = int(rvn_round(tt.model_time/Options.timestep));   // element of NetCDF array that will be written
  time_ind2       =int(rvn_round(tt.model_time/Options.timestep));
  count1[0] = 1;                                                 // writes exactly one time step

  //====================================================================
  //  Hydrographs.nc
  //====================================================================
  int    precip_id;             // variable id in NetCDF for precipitation
  int    qsim_id;               // variable id in NetCDF for simulated outflow
  int    qobs_id;               // variable id in NetCDF for observed outflow
  int    qin_id;                // variable id in NetCDF for observed inflow
  int    qin_net_id;            // variable id in NetCDF for simulated net res inflow

  // (a) count how many values need to be written for q_obs, q_sim, q_in
  int iSim, nSim; // current and total # of sub-basins with simulated outflows
  nSim = 0;
  for (int p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->IsGauged()  && (_pSubBasins[p]->IsEnabled())){nSim++;}
  }

  // (b) allocate memory if necessary
  double *outflow_obs=NULL;   // q_obs
  double *outflow_sim=NULL;   // q_sim
  double *inflow_obs =NULL;   // q_in
  double *outflow_loc=NULL;   // q_loc
  double *inflow_net =NULL;   // q_net_res
  if(nSim>0){
    outflow_sim=new double[nSim];
    outflow_obs=new double[nSim];
    inflow_obs =new double[nSim];
    outflow_loc=new double[nSim];
    inflow_net =new double[nSim];
  }

  // (c) obtain data
  iSim = 0;
  current_prec[0] = NETCDF_BLANK_VALUE;
  if (Options.ave_hydrograph)
  {
    if(tt.model_time != 0.0) { current_prec[0] = GetAveragePrecip(); } //watershed-wide precip+irrigation
    else                     { current_prec[0] = NETCDF_BLANK_VALUE; } // was originally '---'
    for (int p=0;p<_nSubBasins;p++)
    {
      if (_pSubBasins[p]->IsGauged() && (_pSubBasins[p]->IsEnabled())){
        outflow_sim[iSim] = _pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
        outflow_obs[iSim] = NETCDF_BLANK_VALUE;
        for (int i = 0; i < _nObservedTS; i++){
          if (IsContinuousFlowObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
          {
            double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep); //time shift handled in CTimeSeries::Parse
            if ((val != RAV_BLANK_DATA) && (tt.model_time>0)){ outflow_obs[iSim] = val;    }
          }
        }
        outflow_loc[iSim] =_pSubBasins[p]->GetIntegratedLocalOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);

        inflow_obs[iSim]=NETCDF_BLANK_VALUE;
        inflow_net[iSim]=NETCDF_BLANK_VALUE;
        if (_pSubBasins[p]->GetReservoir() != NULL){
          inflow_obs[iSim] = _pSubBasins[p]->GetIntegratedReservoirInflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
          if (Options.write_netresinflow){
              double Qin=inflow_obs[iSim];
              double P  =_pSubBasins[p]->GetReservoir()->GetReservoirPrecipGains(Options.timestep)/(Options.timestep*SEC_PER_DAY);
              double E  =_pSubBasins[p]->GetReservoir()->GetReservoirEvapLosses (Options.timestep)/(Options.timestep*SEC_PER_DAY);
              double GW =_pSubBasins[p]->GetReservoir()->GetReservoirGWLosses   (Options.timestep)/(Options.timestep*SEC_PER_DAY);
              inflow_net[iSim]=Qin+P-E-GW;
          }
        }
        iSim++;
      }
    }
  }
  else {  // point-value hydrograph
    if (tt.model_time != 0.0){current_prec[0] = GetAveragePrecip();} //watershed-wide precip+irrigation
    else                     {current_prec[0] = NETCDF_BLANK_VALUE;} // was originally '---'
    for (int p=0;p<_nSubBasins;p++)
    {
      if (_pSubBasins[p]->IsGauged() && (_pSubBasins[p]->IsEnabled())){
        outflow_sim[iSim] = _pSubBasins[p]->GetOutflowRate();
        outflow_obs[iSim] = NETCDF_BLANK_VALUE;
        for (int i = 0; i < _nObservedTS; i++){
          if (IsContinuousFlowObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
          {
            double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep);
            if ((val != RAV_BLANK_DATA) && (tt.model_time>0)){ outflow_obs[iSim] = val;    }
          }
        }
        outflow_loc[iSim] =_pSubBasins[p]->GetLocalOutflowRate();

        inflow_obs[iSim]=NETCDF_BLANK_VALUE;
        inflow_net[iSim]=NETCDF_BLANK_VALUE;
        if (_pSubBasins[p]->GetReservoir() != NULL){
          inflow_obs[iSim] = _pSubBasins[p]->GetReservoirInflow();
          if (Options.write_netresinflow){
              double Qin=inflow_obs[iSim];
              double P  =_pSubBasins[p]->GetReservoir()->GetReservoirPrecipGains(Options.timestep)/(Options.timestep*SEC_PER_DAY);
              double E  =_pSubBasins[p]->GetReservoir()->GetReservoirEvapLosses (Options.timestep)/(Options.timestep*SEC_PER_DAY);
              double GW =_pSubBasins[p]->GetReservoir()->GetReservoirGWLosses   (Options.timestep)/(Options.timestep*SEC_PER_DAY);
              inflow_net[iSim]=Qin+P-E-GW;
          }
        }
        iSim++;
      }
    }
  }

  // write new time step
  retval = nc_inq_varid      (_HYDRO_ncid, "time",   &time_id);                            HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_HYDRO_ncid, time_id, time_index, count1, &current_time[0]); HandleNetCDFErrors(retval);

  // write precipitation values
  retval = nc_inq_varid      (_HYDRO_ncid, "precip", &precip_id);                          HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_HYDRO_ncid,precip_id,time_index,count1,&current_prec[0]);   HandleNetCDFErrors(retval);

  // write simulated outflow/obs outflow/obs inflow values
  if (nSim > 0){
    start2[0] = int(rvn_round(tt.model_time/Options.timestep)); // element of NetCDF array that will be written
    start2[1] = 0;                                              // element of NetCDF array that will be written
    count2[0] = 1;      // writes exactly one time step
    count2[1] = nSim;   // writes exactly nSim elements
    retval = nc_inq_varid (_HYDRO_ncid, "q_sim", &qsim_id);                             HandleNetCDFErrors(retval);
    retval = nc_inq_varid (_HYDRO_ncid, "q_obs", &qobs_id);                             HandleNetCDFErrors(retval);
    retval = nc_inq_varid (_HYDRO_ncid, "q_in",  &qin_id);                              HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_HYDRO_ncid, qsim_id, start2, count2, &outflow_sim[0]); HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_HYDRO_ncid, qobs_id, start2, count2, &outflow_obs[0]); HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_HYDRO_ncid, qin_id,  start2, count2, &inflow_obs[0]);  HandleNetCDFErrors(retval);
    if (Options.write_localflow){
      int qloc_id;
      retval = nc_inq_varid (_HYDRO_ncid, "q_loc",  &qloc_id);                            HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_HYDRO_ncid, qloc_id, start2, count2, &outflow_loc[0]); HandleNetCDFErrors(retval);
    }
    if (Options.write_netresinflow){
      retval = nc_inq_varid (_HYDRO_ncid, "qnet_in",  &qin_net_id);                     HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_HYDRO_ncid, qin_net_id,  start2, count2, &inflow_net[0]);  HandleNetCDFErrors(retval);
    }
  }

  delete[] outflow_obs;
  delete[] outflow_sim;
  delete[] inflow_obs;
  delete[] outflow_loc;

  //====================================================================
  //  ReservoirStages.nc
  //====================================================================
  if(Options.write_reservoir)
  {
    int    precip_id;             // variable id in NetCDF for precipitation
    int    hsim_id;               // variable id in NetCDF for simulated stage
    int    hobs_id;               // variable id in NetCDF for observed stage

    // (a) count how many values need to be written for h_obs, h_sim
    int iSim,nSim; // current and total # of sub-basins with simulated stages
    nSim = 0;
    for(int p=0;p<_nSubBasins;p++) {
      if(_pSubBasins[p]->IsGauged()  && (_pSubBasins[p]->IsEnabled()) && (_pSubBasins[p]->GetReservoir()!=NULL)) { nSim++; }
    }

    // (b) allocate memory if necessary
    double* stage_obs=NULL;   // h_obs
    double* stage_sim=NULL;   // h_sim
    if(nSim>0) {
      stage_sim=new double[nSim];
      stage_obs=new double[nSim];
    }

    // (c) obtain data
    iSim = 0;
    current_prec[0] = NETCDF_BLANK_VALUE;
    if(tt.model_time != 0.0) { current_prec[0] = GetAveragePrecip(); } //watershed-wide precip+irrigation

    for(int p=0;p<_nSubBasins;p++)
    {
      if(_pSubBasins[p]->IsGauged() && (_pSubBasins[p]->IsEnabled()) && (_pSubBasins[p]->GetReservoir()!=NULL))
      {
        stage_sim[iSim] = _pSubBasins[p]->GetReservoir()->GetResStage();
        stage_obs[iSim] = NETCDF_BLANK_VALUE;
        for(int i = 0; i < _nObservedTS; i++) {
          if(IsContinuousStageObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
          {
            double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep); //time shift handled in CTimeSeries::Parse
            if((val != RAV_BLANK_DATA) && (tt.model_time>0)) { stage_obs[iSim] = val; }
          }
        }
        iSim++;
      }
    }

    // write new time step
    retval = nc_inq_varid      (_RESSTAGE_ncid,"time",&time_id);                              HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_RESSTAGE_ncid,time_id,time_index,count1,&current_time[0]);   HandleNetCDFErrors(retval);

    // write precipitation values
    retval = nc_inq_varid      (_RESSTAGE_ncid,"precip",&precip_id);                          HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_RESSTAGE_ncid,precip_id,time_index,count1,&current_prec[0]); HandleNetCDFErrors(retval);

    // write simulated stage/obs stage values
    if(nSim > 0) {
      start2[0] = int(rvn_round(tt.model_time/Options.timestep)); // element of NetCDF array that will be written
      start2[1] = 0;                                              // element of NetCDF array that will be written
      count2[0] = 1;      // writes exactly one time step
      count2[1] = nSim;   // writes exactly nSim elements
      retval = nc_inq_varid(_RESSTAGE_ncid,"h_sim",&hsim_id);                          HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_RESSTAGE_ncid,"h_obs",&hobs_id);                          HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESSTAGE_ncid,hsim_id,start2,count2,&stage_sim[0]); HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESSTAGE_ncid,hobs_id,start2,count2,&stage_obs[0]); HandleNetCDFErrors(retval);
    }

    delete[] stage_obs;
    delete[] stage_sim;
  }

  //====================================================================
  //  WatershedStorage.nc
  //====================================================================
  if (Options.write_watershed_storage)
  {
    double snowfall      =GetAverageSnowfall();
    double precip        =GetAveragePrecip();
    double channel_stor  =GetTotalChannelStorage();
    double reservoir_stor=GetTotalReservoirStorage();
    double rivulet_stor  =GetTotalRivuletStorage();

    // write new time step
    retval = nc_inq_varid (_STORAGE_ncid, "time",   &time_id);                                 HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_STORAGE_ncid, time_id, time_index, count1, &current_time[0]); HandleNetCDFErrors(retval);

    if(tt.model_time!=0){
      AddSingleValueToNetCDF(_STORAGE_ncid,"rainfall"       ,time_ind2,precip-snowfall);
      AddSingleValueToNetCDF(_STORAGE_ncid,"snowfall"       ,time_ind2,snowfall);
    }
    AddSingleValueToNetCDF(_STORAGE_ncid,"channel_storage"  ,time_ind2,channel_stor);
    AddSingleValueToNetCDF(_STORAGE_ncid,"reservoir_storage",time_ind2,reservoir_stor);
    AddSingleValueToNetCDF(_STORAGE_ncid,"rivulet_storage"  ,time_ind2,rivulet_stor);

    double currentWater=0.0;
    double S;string short_name;
    int iCumPrecip=GetStateVarIndex(ATMOS_PRECIP);
    for (int i=0;i<GetNumStateVars();i++)
    {
      if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iCumPrecip))
	  {
	    S=FormatDouble(GetAvgStateVar(i));
	    short_name = GetStateVarInfo()->GetStateVarLongName(_aStateVarType[i],
                                                          _aStateVarLayer[i],
                                                          GetTransportModel());
	    AddSingleValueToNetCDF(_STORAGE_ncid, short_name.c_str(),time_ind2,S);
	    currentWater+=S;
	  }
    }

    currentWater+=channel_stor+rivulet_stor+reservoir_stor;
    if(tt.model_time==0){
      // \todo [fix]: this fixes a mass balance bug in reservoir simulations, but there is certainly a more proper way to do it
      // JRC: I think somehow this is being double counted in the delta V calculations in the first timestep
      for(int p=0;p<_nSubBasins;p++){
	      if(_pSubBasins[p]->GetReservoir()!=NULL){
	        currentWater+=_pSubBasins[p]->GetIntegratedReservoirInflow(Options.timestep)/2.0/_WatershedArea*MM_PER_METER/M2_PER_KM2;
	        currentWater-=_pSubBasins[p]->GetIntegratedOutflow        (Options.timestep)/2.0/_WatershedArea*MM_PER_METER/M2_PER_KM2;
	      }
      }
    }
    AddSingleValueToNetCDF(_STORAGE_ncid,"total"        ,time_ind2,currentWater);
    AddSingleValueToNetCDF(_STORAGE_ncid,"cum_input"    ,time_ind2,_CumulInput);
    AddSingleValueToNetCDF(_STORAGE_ncid,"cum_outflow"  ,time_ind2,_CumulOutput);
    AddSingleValueToNetCDF(_STORAGE_ncid,"MB_error"     ,time_ind2,FormatDouble((currentWater-_initWater)+(_CumulOutput-_CumulInput)));
  }

  //====================================================================
  //  ForcingFunctions.nc
  //====================================================================
  if(Options.write_forcings)
  {
    force_struct *pFave;
    force_struct faveStruct = GetAverageForcings();
    pFave = &faveStruct;

    // write new time step
    retval = nc_inq_varid      (_FORCINGS_ncid, "time",   &time_id);                            HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_FORCINGS_ncid, time_id, time_index, count1, &current_time[0]); HandleNetCDFErrors(retval);

    // write data
    AddSingleValueToNetCDF(_FORCINGS_ncid,"rainfall"         ,time_ind2,pFave->precip*(1.0-pFave->snow_frac));
    AddSingleValueToNetCDF(_FORCINGS_ncid,"snowfall"         ,time_ind2,pFave->precip*(    pFave->snow_frac));
    AddSingleValueToNetCDF(_FORCINGS_ncid,"temp"             ,time_ind2,pFave->temp_ave);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"temp_daily_min"   ,time_ind2,pFave->temp_daily_min);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"temp_daily_max"   ,time_ind2,pFave->temp_daily_max);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"temp_daily_ave"   ,time_ind2,pFave->temp_daily_ave);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"air_density"      ,time_ind2,pFave->air_dens);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"air_pressure"     ,time_ind2,pFave->air_pres);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"relative_humidity",time_ind2,pFave->rel_humidity);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"cloud_cover"      ,time_ind2,pFave->cloud_cover);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"ET_radiation"     ,time_ind2,pFave->ET_radia);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"SW_radiation"     ,time_ind2,pFave->SW_radia);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"net_SW_radiation" ,time_ind2,pFave->SW_radia_net);
    //AddSingleValueToNetCDF(_FORCINGS_ncid,"SW_radia_subcan"  ,time_ind2,pFave->SW_rad_subcanopy);
    //AddSingleValueToNetCDF(_FORCINGS_ncid,"SW_subcan_net"    ,time_ind2,pFave->SW_subcan_net);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"LW_incoming"      ,time_ind2,pFave->LW_incoming);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"LW_radiation"     ,time_ind2,pFave->cloud_cover);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"wind_velocity"    ,time_ind2,pFave->wind_vel);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"PET"              ,time_ind2,pFave->PET);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"OW_PET"           ,time_ind2,pFave->OW_PET);
    AddSingleValueToNetCDF(_FORCINGS_ncid,"potential_melt"   ,time_ind2,pFave->potential_melt);
  }

  //====================================================================
  //  ReservoirMassBalance.nc
  //====================================================================
  if(Options.write_reservoirMB)
  {

    int    this_id;

    // (a) count how many values need to be written
    int iSim,nSim; // current and total # of sub-basins with simulated stages
    nSim = 0;
    for(int p=0;p<_nSubBasins;p++) {
      if(_pSubBasins[p]->IsGauged()  && (_pSubBasins[p]->IsEnabled()) && (_pSubBasins[p]->GetReservoir()!=NULL)) { nSim++; }
    }

    // (b) allocate memory if necessary
    double* stage=NULL;
    double* area=NULL;
    double* inflow=NULL;
    double* outflow=NULL;
    double* evap=NULL;
    double* seepage=NULL;
    double *stor=NULL;
    double* MBerr=NULL;
    double* losses=NULL;
    double* precip=NULL;
    if(nSim>0) {
      stage  =new double [nSim];
      area   =new double [nSim];
      inflow =new double [nSim];
      outflow=new double [nSim];
      evap   =new double [nSim];
      seepage=new double [nSim];
      stor   =new double [nSim];
      MBerr  =new double [nSim];
      losses =new double [nSim];
      precip =new double [nSim];
    }

    // (c) obtain data
    iSim = 0;
    current_prec[0] = NETCDF_BLANK_VALUE;
    if(tt.model_time != 0.0) { current_prec[0] = GetAveragePrecip(); } //watershed-wide precip+irrigation
    double oldstor;
    for(int p=0;p<_nSubBasins;p++)
    {
      CSubBasin *pSB=_pSubBasins[p];
      if(pSB->IsGauged() && (pSB->IsEnabled()) && (pSB->GetReservoir()!=NULL))
      {
        stage  [iSim]=pSB->GetReservoir()->GetResStage();
        area   [iSim]=pSB->GetReservoir()->GetSurfaceArea();
        inflow [iSim]=pSB->GetIntegratedReservoirInflow(Options.timestep);//m3
        outflow[iSim]=pSB->GetIntegratedOutflow        (Options.timestep);//m3
        evap   [iSim]=pSB->GetReservoir()->GetReservoirEvapLosses (Options.timestep);//m3
        seepage[iSim]=pSB->GetReservoir()->GetReservoirGWLosses   (Options.timestep);//m3
        stor   [iSim]=pSB->GetReservoir()->GetStorage             ();//m3
        oldstor      =pSB->GetReservoir()->GetOldStorage          ();//m3
        losses [iSim]=pSB->GetReservoir()->GetReservoirLosses     (Options.timestep);//m3 = GW+ET
        precip [iSim]=pSB->GetReservoir()->GetReservoirPrecipGains(Options.timestep);//m3
        MBerr  [iSim]=inflow[iSim]-outflow[iSim]-losses[iSim]+precip[iSim]-(stor[iSim]-oldstor);
        //constraint_str[iSim]=pSB->GetReservoir()->GetCurrentConstraint   ();
        if(tt.model_time==0.0){ inflow[iSim]=0.0; }

        iSim++;
      }
    }

    // write new time step
    retval = nc_inq_varid      (_RESMB_ncid,"time",&time_id);                              HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_RESMB_ncid,time_id,time_index,count1,&current_time[0]);   HandleNetCDFErrors(retval);

    // write precipitation values
    retval = nc_inq_varid      (_RESMB_ncid,"precip",&this_id);                          HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_RESMB_ncid,this_id,time_index,count1,&current_prec[0]); HandleNetCDFErrors(retval);

    // write simulated stage/obs stage values
    if(nSim > 0) {
      start2[0] = int(rvn_round(tt.model_time/Options.timestep)); // element of NetCDF array that will be written
      start2[1] = 0;                                              // element of NetCDF array that will be written
      count2[0] = 1;      // writes exactly one time step
      count2[1] = nSim;   // writes exactly nSim elements
      retval = nc_inq_varid(_RESMB_ncid,"stage",&this_id);                          HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&stage[0]);     HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_RESMB_ncid,"area",&this_id);                           HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&area[0]);      HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_RESMB_ncid,"inflow",&this_id);                         HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&inflow[0]);    HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_RESMB_ncid,"outflow",&this_id);                        HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&outflow[0]);   HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_RESMB_ncid,"evap",&this_id);                           HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&evap[0]);      HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_RESMB_ncid,"seepage",&this_id);                        HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&seepage[0]);   HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_RESMB_ncid,"volume",&this_id);                         HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&stor[0]);      HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_RESMB_ncid,"MB_error",&this_id);                       HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&MBerr[0]);     HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_RESMB_ncid,"losses",&this_id);                         HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&losses[0]);    HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_RESMB_ncid,"precip_m3",&this_id);                      HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&precip[0]);    HandleNetCDFErrors(retval);
      //retval = nc_inq_varid(_RESMB_ncid,"constraint",&this_id);                      HandleNetCDFErrors(retval);
      //retval = nc_put_vara_double(_RESMB_ncid,this_id,start2,count2,&constraint[0]); HandleNetCDFErrors(retval);
    }
    delete [] stage;
    delete [] area;
    delete [] inflow;
    delete [] outflow;
    delete [] evap;
    delete [] seepage;
    delete [] stor;
    delete [] MBerr;
    delete [] losses;
    delete [] precip;
  }
#endif
}


//////////////////////////////////////////////////////////////////
/// \brief creates specified output directory, if needed
///
/// \param &Options [in] global model options
//
void PrepareOutputdirectory(const optStruct &Options)
{
  if (Options.output_dir!="")
  {
#if defined(_WIN32)
    _mkdir(Options.output_dir.c_str());
#elif defined(__linux__)
    mkdir(Options.output_dir.c_str(), 0777);
#elif defined(__APPLE__)
    mkdir(Options.output_dir.c_str(),0777);
#elif defined(__unix__)
    mkdir(Options.output_dir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
  }
  g_output_directory=Options.main_output_dir;//necessary evil
}

//////////////////////////////////////////////////////////////////
/// \brief returns directory path given filename
///
/// \param fname [in] filename, e.g., C:\temp\thisfile.txt returns c:\temp
//
string GetDirectoryName(const string &fname)
{
  size_t pos = fname.find_last_of("\\/");
  if (std::string::npos == pos){ return ""; }
  else                         { return fname.substr(0, pos);}
}

//////////////////////////////////////////////////////////////////
/// \brief returns file extension given filename
///
/// \param filename [in] , e.g., C:\temp\thisfile.txt returns txt
//
string GetFileExtension(string filename)
{
  if (filename==""){return ""; }
  return filename.substr(filename.find_last_of(".")+ 1);
}

//////////////////////////////////////////////////////////////////
/// \brief returns directory path given filename and relative path
///
/// \param filename [in] filename, e.g., C:/temp/thisfile.txt returns c:/temp
/// \param relfile [in] filename of reference file
/// e.g.,
///       absolute path of reference file is adopted
///       if filename = something.txt         and relfile= c:/temp/myfile.rvi,  returns c:/temp/something.txt
///
///       relative path of reference file is adopted
///       if filename = something.txt         and relfile= ../dir/myfile.rvi,   returns ../dir/something.txt
///
///       if path of reference file is same as file, then nothing changes
///       if filename = ../temp/something.txt and relfile= ../temp/myfile.rvi,  returns ../temp/something.txt
///
///       if absolute paths of file is given, nothing changes
///       if filename = c:/temp/something.txt and relfile= ../temp/myfile.rvi,  returns c:/temp/something.txt
//
string CorrectForRelativePath(const string filename,const string relfile)
{
  string filedir = GetDirectoryName(relfile); //if a relative path name, e.g., "/path/model.rvt", only returns e.g., "/path"

  if (StringToUppercase(filename).find(StringToUppercase(filedir)) == string::npos)//checks to see if absolute dir already included in redirect filename
  {
    string firstchar  = filename.substr(0, 1);   // if '/' --> absolute path on UNIX systems
    string secondchar = filename.substr(1, 1);   // if ':' --> absolute path on WINDOWS system

    if ( (firstchar.compare("/") != 0) && (secondchar.compare(":") != 0) ){
      // cout << "This is not an absolute filename!  --> " << filename << endl;
      //+"//"
      //cout << "StandardOutput: corrected filename: " << filedir + "//" + filename << endl;
      return filedir + "//" + filename;
    }
  }
  // cout << "StandardOutput: corrected filename: " << filename << endl;
  return filename;
}

//////////////////////////////////////////////////////////////////
/// \brief adds metadata of attribute to NetCDF file
/// \param fileid [in] NetCDF output file id
/// \param time_dimid [in], identifier of time attribute
/// \param shortname [in] attribute short name
/// \param longname [in] attribute long name
/// \param units [in] attribute units as string
//
int NetCDFAddMetadata(const int fileid,const int time_dimid,string shortname,string longname,string units)
{
 int varid(0);
#ifdef _RVNETCDF_
  int retval;
  int dimids[1];
  dimids[0] = time_dimid;

  static double fill_val[] = {NETCDF_BLANK_VALUE};
  static double miss_val[] = {NETCDF_BLANK_VALUE};

  // (a) create variable precipitation
  retval = nc_def_var(fileid,shortname.c_str(),NC_DOUBLE,1,dimids,&varid); HandleNetCDFErrors(retval);

  // (b) add attributes to variable
  retval = nc_put_att_text  (fileid,varid,"units",units.length(),units.c_str());              HandleNetCDFErrors(retval);
  retval = nc_put_att_text  (fileid,varid,"long_name",longname.length(),longname.c_str());    HandleNetCDFErrors(retval);
  retval = nc_put_att_double(fileid,varid,"_FillValue",NC_DOUBLE,1,fill_val);                 HandleNetCDFErrors(retval);
  retval = nc_put_att_double(fileid,varid,"missing_value",NC_DOUBLE,1,miss_val);              HandleNetCDFErrors(retval);
#endif
  return varid;
}
//////////////////////////////////////////////////////////////////
/// \brief adds metadata of attribute to NetCDF file for 2D data (time,nbasins)
/// \param fileid [in] NetCDF output file id
/// \param time_dimid [in], identifier of time attribute
/// \param nbasins_dimid [in], identifier of # basins attribute
/// \param shortname [in] attribute short name
/// \param longname [in] attribute long name
/// \param units [in] attribute units as string
//
int NetCDFAddMetadata2D(const int fileid,const int time_dimid,int nbasins_dimid,string shortname,string longname,string units)
{
  int    varid(0);
#ifdef _RVNETCDF_
  int    retval;
  int    dimids2[2];
  string tmp;

  static double fill_val[] = {NETCDF_BLANK_VALUE};
  static double miss_val[] = {NETCDF_BLANK_VALUE};

  dimids2[0] = time_dimid;
  dimids2[1] = nbasins_dimid;

  // (a) create variable
  retval = nc_def_var(fileid,shortname.c_str(),NC_DOUBLE,2,dimids2,&varid); HandleNetCDFErrors(retval);

  tmp = "basin_name";

  // (b) add attributes to variable
  retval = nc_put_att_text(  fileid,varid,"units",         units.length(),   units.c_str());    HandleNetCDFErrors(retval);
  retval = nc_put_att_text(  fileid,varid,"long_name",     longname.length(),longname.c_str()); HandleNetCDFErrors(retval);
  retval = nc_put_att_double(fileid,varid,"_FillValue",    NC_DOUBLE,1,      fill_val);         HandleNetCDFErrors(retval);
  retval = nc_put_att_double(fileid,varid,"missing_value", NC_DOUBLE,1,      miss_val);         HandleNetCDFErrors(retval);
  retval = nc_put_att_text(  fileid,varid,"coordinates",   tmp.length(),     tmp.c_str());      HandleNetCDFErrors(retval);

#endif
  return varid;
}
//////////////////////////////////////////////////////////////////
/// \brief writes global attributes to netcdf file identified with out_ncid. Must be timeSeries featureType.
/// \param out_ncid [in] NetCDF file identifier
/// \param Options [in] model Options structure
/// \param descript [in] contents of NetCDF description attribute
//
void WriteNetCDFGlobalAttributes(const int out_ncid,const optStruct &Options,const string descript)
{
  ExitGracefullyIf(out_ncid==-9,"WriteNetCDFGlobalAttributes: netCDF file not open",RUNTIME_ERR);

#ifdef _RVNETCDF_
  string att,val;
  int retval(0);
  retval = nc_put_att_text(out_ncid, NC_GLOBAL, "Conventions", strlen("CF-1.6"),          "CF-1.6");           HandleNetCDFErrors(retval);
  retval = nc_put_att_text(out_ncid, NC_GLOBAL, "featureType", strlen("timeSeries"),      "timeSeries");       HandleNetCDFErrors(retval);
  retval = nc_put_att_text(out_ncid, NC_GLOBAL, "history",     strlen("Created by Raven"),"Created by Raven"); HandleNetCDFErrors(retval);
  retval = nc_put_att_text(out_ncid, NC_GLOBAL, "description", strlen(descript.c_str()),  descript.c_str());  HandleNetCDFErrors(retval);

  for(int i=0;i<Options.nNetCDFattribs;i++){
    att=Options.aNetCDFattribs[i].attribute;
    val=Options.aNetCDFattribs[i].value;
    retval = nc_put_att_text(out_ncid,NC_GLOBAL,att.c_str(),strlen(val.c_str()),val.c_str());
    HandleNetCDFErrors(retval);
  }
#endif
}
//////////////////////////////////////////////////////////////////
/// \brief adds single value of attribute 'shortname' linked to time time_index to NetCDF file
/// \param out_ncid [in] NetCDF file identifier
/// \param shortname [in] short name of attribute (e.g., 'rainfall')
/// \param time_index [in] index of current time step
/// \param value [in] value of attribute at current time step
//
void AddSingleValueToNetCDF(const int out_ncid,const string &shortname,const size_t time_index,const double &value)
{
#ifdef _RVNETCDF_
  static size_t count1[1]; //static for speed of execution
  static size_t time_ind[1];
  static double val[1];
  time_ind[0]=time_index;
  count1  [0]=1;
  val     [0]=value;
  int var_id(0);
  int retval(0);
  retval = nc_inq_varid      (out_ncid,shortname.c_str(),&var_id);      HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(out_ncid,var_id,time_ind,count1,&val[0]); HandleNetCDFErrors(retval);
#endif
}
//////////////////////////////////////////////////////////////////
/// \brief writes list of Basin IDs to NetCDF file
/// used for hydrographs.nc, reservoirstages.nc, pollutographs.nc
/// \param ncid [in] NetCDF file identifier
/// \param varid [in] netCDF ID of existing basin ID attribute
/// \param varid_name [in] netCDF ID of existing basin name attribute
/// \param is_res [in] true if this is a reservoir file and only reservoir basins should be included
//
void WriteNetCDFBasinList(const int ncid,const int varid,const int varid_name,const CModel* pModel,bool is_res,const optStruct& Options)
{
#ifdef _RVNETCDF_
  int ibasin = 0;
  int retval;
  size_t      start[1],count[1];                    // determines where and how much will be written to NetCDF
  char* current_basin_name[1];                      // current name or ID  of basin
  current_basin_name[0]=new char[200];
  for(int p=0;p<pModel->GetNumSubBasins();p++) {
    if(pModel->GetSubBasin(p)->IsGauged()  && (pModel->GetSubBasin(p)->IsEnabled())) {
      if (!( (is_res) && (pModel->GetSubBasin(p)->GetReservoir()==NULL))){
        string bID,bname;
        bID   = to_string(pModel->GetSubBasin(p)->GetID()); 
        bname = pModel->GetSubBasin(p)->GetName(); 
        start[0] = ibasin;
        count[0] = 1;
        strcpy(current_basin_name[0],bID.c_str());
        retval = nc_put_vara_string(ncid,varid,start,count,(const char**)current_basin_name);  HandleNetCDFErrors(retval);
        strcpy(current_basin_name[0],bname.c_str());
        retval = nc_put_vara_string(ncid,varid_name,start,count,(const char**)current_basin_name);  HandleNetCDFErrors(retval);
        ibasin++;
      }
    }
  }
  delete[] current_basin_name[0];
#endif

}


//JRC \todo[clean] - find a best place to put this. Might eventually require separate file?
//////////////////////////////////////////////////////////////////
/// \brief return calibration objective function
/// \notes right now only supports hydrograph goodness of fit metrics at one subbasin
/// \param &calib_SBID [in] target subbasin ID
/// \param &calib_Obj [in] calibration objective diagnostics (e.g., NSE)
//
double CModel::GetObjFuncVal(long long calib_SBID,diag_type calib_Obj, const string calib_period) const
{
  double starttime=0.0;
  double endtime  =0.0;
  comparison compare=COMPARE_GREATERTHAN;
  double thresh=-ALMOST_INF;
  double objval;
  int ii=DOESNT_EXIST; // observation index
  int jj=DOESNT_EXIST; // diagnostic measure

  //- grab diagnostic information -----------------------------------
  for(int d=0;d<_nDiagPeriods;d++) {
    if(_pDiagPeriods[d]->GetName()==calib_period) {
      starttime=_pDiagPeriods[d]->GetStartTime();
      endtime  =_pDiagPeriods[d]->GetEndTime();
      compare  =_pDiagPeriods[d]->GetComparison();
      thresh   =_pDiagPeriods[d]->GetThreshold();
    }
  }

  for(int i=0;i<_nObservedTS;i++)
  {
    if((_pObservedTS[i]->GetName()=="HYDROGRAPH") && (_pObservedTS[i]->GetLocID()==calib_SBID)) { ii=i; }
  }
  for(int j=0; j<_nDiagnostics;j++) {
    if(_pDiagnostics[j]->GetType()==calib_Obj) { jj=j; }
  }
  ExitGracefullyIf(ii==DOESNT_EXIST,"GetObjFuncVal: unable to find calibration target time series (hydrograph in basin :CalibrationSBID)",BAD_DATA);
  ExitGracefullyIf(jj==DOESNT_EXIST,"GetObjFuncVal: unable to find calibration target diagnostic ",BAD_DATA);
  ExitGracefullyIf(endtime==0.0,    "GetObjFuncVal: unable to find calibration period with this name ",BAD_DATA);

  //- Calculate objective function -----------------------------------
  objval=_pDiagnostics[jj]->CalculateDiagnostic(_pModeledTS[ii],_pObservedTS[ii],_pObsWeightTS[ii],starttime,endtime,compare,thresh,*_pOptStruct);

  if((calib_Obj==DIAG_NASH_SUTCLIFFE) ||
     (calib_Obj==DIAG_KLING_GUPTA))
  {
    objval*=-1;
  }
  return objval;
}
