/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team

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
#endif
int NetCDFAddMetadata(const int fileid,const int time_dimid,string shortname,string longname,string units);
int NetCDFAddMetadata2D(const int fileid,const int time_dimid,int nbasins_dimid,string shortname,string longname,string units);
//////////////////////////////////////////////////////////////////
/// \brief returns true if specified observation time series is the flow series for subbasin SBID
/// \param pObs [in] observation time series
/// \param SBID [in] subbasin ID
//
bool IsContinuousFlowObs(CTimeSeriesABC *pObs,long SBID)
{
 // clears up  terribly ugly repeated if statements
  if(pObs==NULL){return false;}
  if (s_to_l(pObs->GetTag().c_str()) != SBID){ return false; }//SBID is correct
  if(pObs->GetType() != CTimeSeriesABC::ts_regular){ return false; }
  return (!strcmp(pObs->GetName().c_str(),"HYDROGRAPH")); //name ="HYDROGRAPH"      
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if specified observation time series is the reservoir stage series for subbasin SBID
/// \param pObs [in] observation time series
/// \param SBID [in] subbasin ID
//
bool IsContinuousStageObs(CTimeSeriesABC *pObs,long SBID)
{
 // clears up  terribly ugly repeated if statements
  if(pObs==NULL){return false;}
  return (
    (!strcmp(pObs->GetName().c_str(),"RESERVOIR_STAGE")) &&
    (s_to_l(pObs->GetTag().c_str()) == SBID) &&
    (pObs->GetType() == CTimeSeriesABC::ts_regular)
    );
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if specified observation time series is the reservoir inflow series for subbasin SBID
/// \param pObs [in] observation time series
/// \param SBID [in] subbasin ID
//
bool IsContinuousInflowObs(CTimeSeriesABC *pObs, long SBID)
{
	// clears up  terribly ugly repeated if statements
	if (pObs == NULL) { return false; }
	return (
		(!strcmp(pObs->GetName().c_str(), "RESERVOIR_INFLOW")) &&
		(s_to_l(pObs->GetTag().c_str()) == SBID) &&
		(pObs->GetType() == CTimeSeriesABC::ts_regular)
		);
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if specified observation time series is the reservoir inflow series for subbasin SBID
/// \param pObs [in] observation time series
/// \param SBID [in] subbasin ID
//
bool IsContinuousNetInflowObs(CTimeSeriesABC *pObs, long SBID)
{
	// clears up  terribly ugly repeated if statements
	if (pObs == NULL) { return false; }
	return (
		(!strcmp(pObs->GetName().c_str(), "RESERVOIR_NETINFLOW")) &&
		(s_to_l(pObs->GetTag().c_str()) == SBID) &&
		(pObs->GetType() == CTimeSeriesABC::ts_regular)
		);
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
    _pCustomOutputs[c]->CloseFiles();
  }
  _pTransModel->CloseOutputFiles();
  if ( _STORAGE.is_open()){ _STORAGE.close();}
  if (   _HYDRO.is_open()){   _HYDRO.close();}
  if (_FORCINGS.is_open()){_FORCINGS.close();}
  if (_RESSTAGE.is_open()){_RESSTAGE.close();}

#ifdef _RVNETCDF_
  
  /* close netcdfs */
  int    retval;      // error value for NetCDF routines
  if (_HYDRO_ncid != -9)    {retval = nc_close(_HYDRO_ncid);    HandleNetCDFErrors(retval); }
  _HYDRO_ncid    = -9;
  if (_STORAGE_ncid != -9)  {retval = nc_close(_STORAGE_ncid);  HandleNetCDFErrors(retval); }
  _STORAGE_ncid  = -9;
  if (_FORCINGS_ncid != -9) {retval = nc_close(_FORCINGS_ncid); HandleNetCDFErrors(retval); }
  _FORCINGS_ncid = -9;
  
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

  if (Options.output_format==OUTPUT_STANDARD)
  {

    //WatershedStorage.csv
    //--------------------------------------------------------------
    tmpFilename=FilenamePrepare("WatershedStorage.csv",Options);
    _STORAGE.open(tmpFilename.c_str());
    if (_STORAGE.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }

    int iAtmPrecip=GetStateVarIndex(ATMOS_PRECIP);
    _STORAGE<<"time [d],date,hour,rainfall [mm/day],snowfall [mm/d SWE],Channel Storage [mm],Reservoir Storage [mm],Rivulet Storage [mm]";
    for (i=0;i<GetNumStateVars();i++){
      if (CStateVariable::IsWaterStorage(_aStateVarType[i])){
        if (i!=iAtmPrecip){
          _STORAGE<<","<<CStateVariable::GetStateVarLongName(_aStateVarType[i],_aStateVarLayer[i])<<" [mm]";
          //_STORAGE<<","<<CStateVariable::SVTypeToString(_aStateVarType[i],_aStateVarLayer[i])<<" [mm]";
        }
      }
    }
    _STORAGE<<", Total [mm], Cum. Inputs [mm], Cum. Outflow [mm], MB Error [mm]"<<endl;

    //Hydrographs.csv
    //--------------------------------------------------------------
    tmpFilename=FilenamePrepare("Hydrographs.csv",Options);
    _HYDRO.open(tmpFilename.c_str());
    if (_HYDRO.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }

    _HYDRO<<"time,date,hour";
    _HYDRO<<",precip [mm/day]";
    for (p=0;p<_nSubBasins;p++){
      if (_pSubBasins[p]->IsGauged() && _pSubBasins[p]->IsEnabled()){
        string name;
        if (_pSubBasins[p]->GetName()==""){_HYDRO<<",ID="<<_pSubBasins[p]->GetID()  <<" [m3/s]";}
        else                              {_HYDRO<<","   <<_pSubBasins[p]->GetName()<<" [m3/s]";}
        //if (Options.print_obs_hydro)
        {
          for (int i = 0; i < _nObservedTS; i++){
            if (IsContinuousFlowObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
            {
              if (_pSubBasins[p]->GetName()==""){_HYDRO<<",ID="<<_pSubBasins[p]->GetID()  <<" (observed) [m3/s]";}
              else                              {_HYDRO<<","   <<_pSubBasins[p]->GetName()<<" (observed) [m3/s]";}
            }
          }
        }

        if (_pSubBasins[p]->GetReservoir() != NULL){
          if (_pSubBasins[p]->GetName()==""){_HYDRO<<",ID="<<_pSubBasins[p]->GetID()  <<" (res. inflow) [m3/s]";}
          else                              {_HYDRO<<","   <<_pSubBasins[p]->GetName()<<" (res. inflow) [m3/s]";}
        }
      }
    }
    _HYDRO<<endl;
  }
  else if (Options.output_format==OUTPUT_ENSIM)
  {
    WriteEnsimStandardHeaders(Options);
  }
  else if (Options.output_format==OUTPUT_NETCDF)
  {
    WriteNetcdfStandardHeaders(Options);  // creates NetCDF files, writes dimensions and creates variables (without writing actual values)   
  }

  //WatershedEnergyStorage.csv
  //--------------------------------------------------------------
  if (Options.write_energy)
  {
    ofstream EN_STORAGE;
    tmpFilename=FilenamePrepare("WatershedEnergyStorage.csv",Options);
    EN_STORAGE.open(tmpFilename.c_str());
    if (EN_STORAGE.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }

    EN_STORAGE<<"time[d],date,hour,temp[C],net incoming [MJ/m2/d]";
    for (i=0;i<GetNumStateVars();i++){
      if (CStateVariable::IsEnergyStorage(_aStateVarType[i])){
        EN_STORAGE<<","<<CStateVariable::SVTypeToString(_aStateVarType[i],_aStateVarLayer[i])<<" [MJ/m2]";
      }
    }
    EN_STORAGE<<", Total [MJ/m2], Cum. In [MJ/m2], Cum. Out [MJ/m2], EB Error [MJ/m2]"<<endl;
    EN_STORAGE.close();
  }

  //ReservoirStages.csv
  //--------------------------------------------------------------
  if((Options.write_reservoir) && (Options.output_format!=OUTPUT_NONE))
  {
    tmpFilename=FilenamePrepare("ReservoirStages.csv",Options);
    _RESSTAGE.open(tmpFilename.c_str());
    if(_RESSTAGE.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
  
    _RESSTAGE<<"time,date,hour";
    _RESSTAGE<<",precip [mm/day]";
    for(p=0;p<_nSubBasins;p++){
      if((_pSubBasins[p]->IsGauged()) && (_pSubBasins[p]->IsEnabled()) && (_pSubBasins[p]->GetReservoir()!=NULL)) {
        string name;
        if(_pSubBasins[p]->GetName()==""){ _RESSTAGE<<",ID="<<_pSubBasins[p]->GetID()  <<" "; }
        else                             { _RESSTAGE<<","   <<_pSubBasins[p]->GetName()<<" "; }
      }
      //if (Options.print_obs_hydro)
      {
        for(int i = 0; i < _nObservedTS; i++){
          if(IsContinuousStageObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
          {
            if(_pSubBasins[p]->GetName()==""){ _RESSTAGE<<",ID="<<_pSubBasins[p]->GetID()  <<" (observed) [m3/s]"; }
            else                             { _RESSTAGE<<","   <<_pSubBasins[p]->GetName()<<" (observed) [m3/s]"; }
          }
        }
      }
    }
    _RESSTAGE<<endl;
  }

  //ReservoirMassBalance.csv
  //--------------------------------------------------------------
  if ((Options.write_reservoirMB) && (Options.output_format!=OUTPUT_NONE))
  {
    ofstream RES_MB;
    string name;
    cout<<"RESERVOIR"<<endl;
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
        else                              { name=_pSubBasins[p]->GetName(); }
        RES_MB<<","   <<name<<" inflow [m3]";
        RES_MB<<","   <<name<<" outflow [m3]";
        RES_MB<<","   <<name<<" volume [m3]";
        RES_MB<<","   <<name<<" losses [m3]";
        RES_MB<<","   <<name<<" MB error [m3]";
      }
    }
    RES_MB<<endl;
    RES_MB.close();
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
    MB<<"time [d],date,hour";
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
        MB<<","<<CStateVariable::SVTypeToString(typ,ind);
      }
    }
    MB<<endl;
    MB<<",,to:";
    for (j=0;j<_nProcesses;j++){
      for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
        sv_type typ=GetStateVarType (_pProcesses[j]->GetToIndices()[q]);
        int     ind=GetStateVarLayer(_pProcesses[j]->GetToIndices()[q]);
        MB<<","<<CStateVariable::SVTypeToString(typ,ind);
      }
    }
    MB<<endl;
    MB.close();
  }

  //WatershedMassEnergyBalance.csv
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
    HGMB<<"time [d],date,hour";
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
        HGMB<<","<<CStateVariable::SVTypeToString(typ,ind);
      }
    }
    HGMB<<endl;
    HGMB<<",,to:";
    for (j=0;j<_nProcesses;j++){
      for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
        sv_type typ=GetStateVarType (_pProcesses[j]->GetToIndices()[q]);
        int     ind=GetStateVarLayer(_pProcesses[j]->GetToIndices()[q]);
        HGMB<<","<<CStateVariable::SVTypeToString(typ,ind);
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
    MB<<"time[d],date,hour";
    bool first;
    for (int i=0;i<_nStateVars;i++){
      if (CStateVariable::IsWaterStorage(_aStateVarType[i]))
      {
        MB<<","<<CStateVariable::SVTypeToString(_aStateVarType[i],_aStateVarLayer[i]);
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
    for (int i=0;i<_nStateVars;i++){
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
            if (pProc->GetLateralToIndices()[q]==i){MB<<","<<GetProcessName(_pProcesses[j]->GetProcessType());break;}
            if (pProc->GetLateralFromIndices()[q]==i){MB<<","<<GetProcessName(_pProcesses[j]->GetProcessType());break;}
          }
        }
        MB<<",cumulative,storage,error";
      }
    }
    MB<<endl;
    MB.close();
  }

  //ForcingFunctions.csv
  //--------------------------------------------------------------
  if (Options.write_forcings)
  {
    tmpFilename=FilenamePrepare("ForcingFunctions.csv",Options);
    _FORCINGS.open(tmpFilename.c_str());
    if (_FORCINGS.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    _FORCINGS<<"time [d],date,hour,day_angle,";
    _FORCINGS<<" rain [mm/d], snow [mm/d], temp [C], temp_daily_min [C], temp_daily_max [C],temp_daily_ave [C],temp_monthly_min [C],temp_monthly_max [C],";
    _FORCINGS<<" air dens. [kg/m3], air pres. [KPa], rel hum. [-],";
    _FORCINGS<<" cloud cover [-],";
    _FORCINGS<<" ET radiation [MJ/m2/d], SW radiation [MJ/m2/d], net SW radiation [MJ/m2/d], LW radiation [MJ/m2/d], wind vel. [m/s],";
    _FORCINGS<<" PET [mm/d], OW PET [mm/d],";
    _FORCINGS<<" daily correction [-], potential melt [mm/d]";
    _FORCINGS<<endl;
  }

  // HRU Storage files
  //--------------------------------------------------------------
  if (_pOutputGroup!=NULL){
    for (int kk=0; kk<_pOutputGroup->GetNumHRUs();kk++)
    {
      ofstream HRUSTOR;
      tmpFilename="HRUStorage_"+to_string(_pOutputGroup->GetHRU(kk)->GetID())+".csv";
      tmpFilename=FilenamePrepare(tmpFilename,Options);
      HRUSTOR.open(tmpFilename.c_str());
      if (HRUSTOR.fail()){
        ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
      }
      int iAtmPrecip=GetStateVarIndex(ATMOS_PRECIP);
      HRUSTOR<<"time [d],date,hour,rainfall [mm/day],snowfall [mm/d SWE]";
      for (i=0;i<GetNumStateVars();i++){
        if (CStateVariable::IsWaterStorage(_aStateVarType[i])){
          if (i!=iAtmPrecip){
            HRUSTOR<<","<<CStateVariable::GetStateVarLongName(_aStateVarType[i],_aStateVarLayer[i])<<" [mm]";
            //HRUSTOR<<","<<CStateVariable::SVTypeToString(_aStateVarType[i],_aStateVarLayer[i])<<" [mm]";
          }
        }
      }
      HRUSTOR<<", Total [mm]"<<endl;
    }
  }

  // Custom output files
  //--------------------------------------------------------------
  for (int c=0;c<_nCustomOutputs;c++)
  {
    _pCustomOutputs[c]->WriteFileHeader(Options);
  }

  // Transport output files
  //--------------------------------------------------------------
  if (Options.output_format==OUTPUT_STANDARD)
  {
    _pTransModel->WriteOutputFileHeaders(Options);
  }
  else if (Options.output_format==OUTPUT_ENSIM)
  {
    _pTransModel->WriteEnsimOutputFileHeaders(Options);
  }

  //raven_debug.csv
  //--------------------------------------------------------------
  if (Options.debug_mode)
  {
    ofstream DEBUG;
    DEBUG.open("raven_debug.csv");
    if (DEBUG.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    DEBUG<<"time[d],date,hour,debug1,debug2,debug3,debug4,debug5"<<endl;
    DEBUG.close();
  }

  //opens and closes diagnostics.csv so that this warning doesn't show up at end of simulation
  //--------------------------------------------------------------
  if ((_nObservedTS>0) && (_nDiagnostics>0))
  {
    ofstream DIAG;
    string tmpFilename;
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
  bool    silent=true;
  bool    quiet=true;
  double  t;

  string tmpFilename;
  
  if ((tt.model_time==0) && (Options.suppressICs)){return;}

  //converts the 'write every x timesteps' into a 'write at time y' value
  output_int = Options.output_interval * Options.timestep;
  mod_final  = ffmod(tt.model_time,output_int);
  //cout<<"time: "<<setprecision(12)<<t<<"  timestep: "<<setprecision(12)<<Options.timestep<<"  output at: "<<setprecision(12)<<output_int<<"  answer: "<<setprecision(12)<<test_mod<<"  mod: "<<setprecision(12)<<mod_final<<endl;

  iCumPrecip=GetStateVarIndex(ATMOS_PRECIP);

  if(fabs(mod_final) <= 0.5*Options.timestep)  //checks to see if sufficiently close to timestep
                                               //(this should account for any roundoff error in timestep calcs)
  {
    thisdate=tt.date_string;                   //refers to date and time at END of time step
    thishour=DecDaysToHours(tt.julian_day);
    t       =tt.model_time;

    time_struct prev;
    JulianConvert(t-Options.timestep,Options.julian_start_day,Options.julian_start_year,prev); //get start of time step, prev

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
    if ((quiet) && (!Options.silent) && (tt.day_of_month==1) && ((tt.julian_day)-floor(tt.julian_day)<Options.timestep/2))
    {
      cout<<thisdate <<endl;
    }
    if(!silent)
    {
      cout <<thisdate<<" "<<thishour<<":";
      if (t!=0){cout <<" | P: "<< setw(6)<<setiosflags(ios::fixed) << setprecision(2)<<GetAveragePrecip();}
      else     {cout <<" | P: ------";}
    }


    //Write current state of water storage in system to WatershedStorage.csv (ALWAYS DONE)
    //----------------------------------------------------------------
    if (Options.output_format==OUTPUT_STANDARD)
    {
      double snowfall      =GetAverageSnowfall();
      double precip        =GetAveragePrecip();
      double channel_stor  =GetTotalChannelStorage();
      double reservoir_stor=GetTotalReservoirStorage();
      double rivulet_stor  =GetTotalRivuletStorage();

      _STORAGE<<tt.model_time <<","<<thisdate<<","<<thishour;

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
        }
      }

      _STORAGE<<","<<currentWater<<","<<_CumulInput<<","<<_CumulOutput<<","<<FormatDouble((currentWater-_initWater)+(_CumulOutput-_CumulInput));
      _STORAGE<<endl;

      //Write hydrographs for gauged watersheds (ALWAYS DONE)
      //----------------------------------------------------------------
      if ((Options.ave_hydrograph) && (t!=0.0))
      {
        
        _HYDRO<<usetime<<","<<usedate<<","<<usehour<<","<<GetAveragePrecip();
        for (int p=0;p<_nSubBasins;p++){
          if (_pSubBasins[p]->IsGauged()  && (_pSubBasins[p]->IsEnabled()))
          {
            _HYDRO<<","<<_pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);

            //if (Options.print_obs_hydro)
            {
              for (int i = 0; i < _nObservedTS; i++)
              {
                if (IsContinuousFlowObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
                {
                  //int nn=int(floor((tt.model_time+TIME_CORRECTION)/ Options.timestep));
                  //double val=_pObservedTS[i]->GetSampledValue(nn); //fails for interval>timestep
                  double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep); //time shift handled in CTimeSeries::Parse
                  if ((val != RAV_BLANK_DATA) && (tt.model_time>0)){ _HYDRO << "," << val; }
                  else                                             { _HYDRO << ",";       }
                }
              }
            }
            if (_pSubBasins[p]->GetReservoir() != NULL){
              _HYDRO<<","<<_pSubBasins[p]->GetIntegratedReservoirInflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
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
          if(t!=0){ _HYDRO<<","<<GetAveragePrecip(); }//watershed-wide precip
          else     { _HYDRO<<",---"; }
          for(int p=0;p<_nSubBasins;p++){
            if(_pSubBasins[p]->IsGauged()  && (_pSubBasins[p]->IsEnabled()))
            {
              _HYDRO<<","<<_pSubBasins[p]->GetOutflowRate();

              //if (Options.print_obs_hydro)
              {
                for(int i = 0; i < _nObservedTS; i++){
                  if(IsContinuousFlowObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
                  {
                    //int nn=int(floor((tt.model_time+TIME_CORRECTION)/ Options.timestep));
                    //double val=_pObservedTS[i]->GetSampledValue(nn); //fails for interval>timestep
                    double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep);
                    if((val != RAV_BLANK_DATA) && (tt.model_time>0)){ _HYDRO << "," << val; }
                    else                                                         { _HYDRO << ","; }
                  }
                }
              }
              if(_pSubBasins[p]->GetReservoir() != NULL){
                _HYDRO<<","<<_pSubBasins[p]->GetReservoirInflow();
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

    //ReservoirStages.csv
    //--------------------------------------------------------------
    if ((Options.write_reservoir) && (Options.output_format!=OUTPUT_NONE))
    {
			if((Options.period_starting) && (t==0)){}//don't write anything at time zero
      else{
	      _RESSTAGE<< t<<","<<thisdate<<","<<thishour<<","<<GetAveragePrecip();
	      for (int p=0;p<_nSubBasins;p++){
	        if ((_pSubBasins[p]->IsGauged())  && (_pSubBasins[p]->IsEnabled()) && (_pSubBasins[p]->GetReservoir()!=NULL)) {
	          _RESSTAGE<<","<<_pSubBasins[p]->GetReservoir()->GetResStage();
	        }
	        //if (Options.print_obs_hydro)
	        {
	          for (int i = 0; i < _nObservedTS; i++){
	            if (IsContinuousStageObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
	            {
	              //int nn=int(floor((tt.model_time+TIME_CORRECTION)/ Options.timestep));
	              double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep);
	              if ((val != RAV_BLANK_DATA) && (tt.model_time>0)){ _RESSTAGE << "," << val; }
	              else                                             { _RESSTAGE << ",";       }
	            }
	          }
	        }

	      }
	      _RESSTAGE<<endl;
			}
    }

    //ReservoirMassBalance.csv
    //----------------------------------------------------------------
    if((Options.write_reservoirMB) && (Options.output_format!=OUTPUT_NONE))
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
        double in,out,loss,stor,oldstor;
        for(int p=0;p<_nSubBasins;p++){
          if((_pSubBasins[p]->IsGauged()) &&  (_pSubBasins[p]->IsEnabled()) && (_pSubBasins[p]->GetReservoir()!=NULL)) {

            string name;
            if(_pSubBasins[p]->GetName()==""){ name=to_string(_pSubBasins[p]->GetID())+"="+to_string(_pSubBasins[p]->GetID()); }
            else                             { name=_pSubBasins[p]->GetName(); }

            in     =_pSubBasins[p]->GetIntegratedReservoirInflow(Options.timestep);//m3
            out    =_pSubBasins[p]->GetIntegratedOutflow(Options.timestep);//m3
            stor   =_pSubBasins[p]->GetReservoir()->GetStorage();//m3
            oldstor=_pSubBasins[p]->GetReservoir()->GetOldStorage();//m3
            loss   =_pSubBasins[p]->GetReservoir()->GetReservoirLosses(Options.timestep);//m3
            if(tt.model_time==0.0){ in=0.0; }
            RES_MB<<","<<in<<","<<out<<","<<stor<<","<<loss<<","<<in-out-loss-(stor-oldstor);
          }
        }
        RES_MB<<endl;
        RES_MB.close();
      }
    }


    //Write current state of energy storage in system to WatershedEnergyStorage.csv
    //----------------------------------------------------------------
    double sum=0.0;
    if (Options.write_energy)
    {
      force_struct F=GetAverageForcings();

      ofstream EN_STORAGE;
      tmpFilename=FilenamePrepare("WatershedEnergyStorage.csv",Options);
      EN_STORAGE.open(tmpFilename.c_str(),ios::app);

      EN_STORAGE<<t<<","<<thisdate<<","<<thishour;
      EN_STORAGE<<","<<F.temp_ave;
      EN_STORAGE<<",TMP_DEBUG";
      //  STOR<<","<<GetAverageNetRadiation();//TMP DEBUG
      for (i=0;i<GetNumStateVars();i++)
      {
        if (CStateVariable::IsEnergyStorage(_aStateVarType[i]))
        {
          S=GetAvgStateVar(i);
          if (!silent){cout<<"  |"<< setw(6)<<setiosflags(ios::fixed) << setprecision(2)<<S;}
          EN_STORAGE<<","<<S;
          sum+=S;
        }
      }
      EN_STORAGE<<","<<sum<<","<<_CumEnergyGain<<","<<_CumEnergyLoss<<","<<_CumEnergyGain-sum-_CumEnergyLoss;
      EN_STORAGE<<endl;
      EN_STORAGE.close();
    }

    if (!silent){cout<<endl;}

    //ExhaustiveMassBalance.csv
    //--------------------------------------------------------------
    if (Options.write_exhaustiveMB)
    {
      if((Options.period_starting) && (t==0)){}//don't write anything at time zero
      else{
        int j,js,q;
        double cumsum;

        ofstream MB;
        tmpFilename=FilenamePrepare("ExhaustiveMassBalance.csv",Options);
        MB.open(tmpFilename.c_str(),ios::app);

        MB<<usetime<<","<<usedate<<","<<usehour;
        for(int i=0;i<_nStateVars;i++)
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
              for(int q=0;q<_pProcesses[j]->GetNumLatConnections();q++){
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

    //Write wshed-averaged forcing functions to ForcingFunctions.csv
    //----------------------------------------------------------------
    if (Options.write_forcings)
    {
      if((Options.period_starting) && (t==0)){}//don't write anything at time zero
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
        _FORCINGS<<pFave->LW_radia<<",";
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
    if (Options.output_format==OUTPUT_STANDARD)
    {
      _pTransModel->WriteMinorOutput(Options,tt);
    }
    else if (Options.output_format==OUTPUT_ENSIM)
    {
      _pTransModel->WriteEnsimMinorOutput(Options,tt);
    }

    // raven_debug.csv
    //--------------------------------------------------------------
    if (Options.debug_mode)
    {
      ofstream DEBUG;
      DEBUG.open("raven_debug.csv",ios::app);

      DEBUG<<t<<","<<thisdate<<","<<thishour;
      for(int i=0;i<5;i++){DEBUG<<","<<g_debug_vars[i];}
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
        tmpFilename="HRUStorage_"+to_string(_pOutputGroup->GetHRU(kk)->GetID())+".csv";
        tmpFilename=FilenamePrepare(tmpFilename,Options);
        HRUSTOR.open(tmpFilename.c_str(),ios::app);

        const force_struct *F=_pOutputGroup->GetHRU(kk)->GetForcingFunctions();

        HRUSTOR<<tt.model_time <<","<<thisdate<<","<<thishour;

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
  } // end of write output interval if statement

  // Custom output files
  //--------------------------------------------------------------
  for (int c=0;c<_nCustomOutputs;c++)
  {
    _pCustomOutputs[c]->WriteCustomOutput(tt,Options);
  }

  // Write major output, if necessary
  //--------------------------------------------------------------
  if ((_nOutputTimes>0) && (tt.model_time>_aOutputTimes[_currOutputTimeInd]-0.5*Options.timestep))
  {
    _currOutputTimeInd++;
    tmpFilename="state_"+tt.date_string;
    WriteMajorOutput(tmpFilename,Options,tt,false);
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
void CModel::WriteMajorOutput(string solfile, const optStruct &Options, const time_struct &tt, bool final) const
{
  int i,k;
  string tmpFilename;

  // WRITE {RunName}_solution.rvc - final state variables file
  ofstream OUT;
  tmpFilename=FilenamePrepare(solfile+".rvc",Options);
  OUT.open(tmpFilename.c_str());
  if (OUT.fail()){
    WriteWarning(("CModel::WriteMajorOutput: Unable to open output file "+tmpFilename+" for writing.").c_str(),Options.noisy);
  }
  OUT<<":TimeStamp "<<tt.date_string<<" "<<DecDaysToHours(tt.julian_day)<<endl;

  //Header--------------------------
  OUT<<":HRUStateVariableTable"<<endl;
  OUT<<"  :Attributes,";
  for (i=0;i<GetNumStateVars();i++)
  {
    OUT<<CStateVariable::SVTypeToString(_aStateVarType[i],_aStateVarLayer[i]);
    if (i!=GetNumStateVars()-1){OUT<<",";}
  }
  OUT<<endl;
  OUT<<"  :Units,";
  for (i=0;i<GetNumStateVars();i++)
  {
    OUT<<CStateVariable::GetStateVarUnits(_aStateVarType[i]);
    if (i!=GetNumStateVars()-1){OUT<<",";}
  }
  OUT<<endl;
  //Data----------------------------
  for (k=0;k<_nHydroUnits;k++)
  {
    OUT<<"  "<<_pHydroUnits[k]->GetID()<<",";
    for (i=0;i<GetNumStateVars();i++)
    {
      OUT<<_pHydroUnits[k]->GetStateVarValue(i);
      if (i!=GetNumStateVars()-1){OUT<<",";}
    }
    OUT<<endl;
  }
  OUT<<":EndHRUStateVariableTable"<<endl;

  //By basin------------------------
  OUT<<":BasinStateVariables"<<endl;
  for (int p=0;p<_nSubBasins;p++){
    OUT<<"  :BasinIndex "<<_pSubBasins[p]->GetID()<<",";
    _pSubBasins[p]->WriteToSolutionFile(OUT);
  }
  OUT<<":EndBasinStateVariables"<<endl;
  OUT.close();

  if(Options.write_channels){
    CChannelXSect::WriteRatingCurves();
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
  if (Options.pavics){

    ofstream PROGRESS((Options.output_dir+"Raven_progress.txt").c_str());
    if (PROGRESS.fail()){
      ExitGracefully("ParseInput:: Unable to open Raven_progress.txt. Bad output directory specified?",RUNTIME_ERR);
    }
    
    float total_time = (float(total_steps) * float(elapsed_time) / float(elapsed_steps)) / CLOCKS_PER_SEC;
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
  if(!Options.silent){
    cout <<"==MODEL SUMMARY======================================="<<endl;
    cout <<"       Model Run: "<<Options.run_name    <<endl;
    cout <<"    rvi filename: "<<Options.rvi_filename<<endl;
    cout <<"Output Directory: "<<Options.output_dir  <<endl;
    cout <<"     # SubBasins: "<<GetNumSubBasins()   << " ("<< rescount << " reservoirs) ("<<SBdisablecount<<" disabled)"<<endl;
    cout <<"          # HRUs: "<<GetNumHRUs()        << " ("<<disablecount<<" disabled)"<<endl;
    cout <<"        # Gauges: "<<GetNumGauges()      <<endl;
    cout <<"#State Variables: "<<GetNumStateVars()   <<endl;
    for (int i=0;i<GetNumStateVars();i++){
      //don't write if convolution storage or advection storage?
      cout<<"                - ";
      cout<<CStateVariable::GetStateVarLongName(_aStateVarType[i],_aStateVarLayer[i])<<" (";
      cout<<CStateVariable::SVTypeToString     (_aStateVarType[i],_aStateVarLayer[i])<<")"<<endl;
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
    cout <<"       Time step: "<<Options.timestep            <<" d"<<endl;
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
void CModel::RunDiagnostics (const optStruct &Options)
{
  if ((_nObservedTS==0) || (_nDiagnostics==0)) {return;}

  ofstream DIAG;
  string tmpFilename;
  tmpFilename=FilenamePrepare("Diagnostics.csv",Options);
  DIAG.open(tmpFilename.c_str());
  if (DIAG.fail()){
    ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
  }
  //header
  DIAG<<"observed data series,filename,";
  for (int j=0; j<_nDiagnostics;j++){
    DIAG<<_pDiagnostics[j]->GetName()<<",";
  }
  DIAG<<endl;
  //body
  for (int i=0;i<_nObservedTS;i++)
  {
    DIAG<<_pObservedTS[i]->GetName()<<","<<_pObservedTS[i]->GetSourceFile() <<",";
    for (int j=0; j<_nDiagnostics;j++){
      DIAG<<_pDiagnostics[j]->CalculateDiagnostic(_pModeledTS[i],_pObservedTS[i],_pObsWeightTS[i],Options)<<",";
    }
    DIAG<<endl;
  }
  DIAG.close();
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

  JulianConvert(0.0, Options.julian_start_day, Options.julian_start_year, tt);//start of the timestep
  JulianConvert(Options.timestep, Options.julian_start_day, Options.julian_start_year, tt2);//end of the timestep

  //WatershedStorage.tb0
  //--------------------------------------------------------------
  int iAtmPrecip = GetStateVarIndex(ATMOS_PRECIP);
  string tmpFilename;

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
    _STORAGE << ":CreationDate  " << GetCurrentTime() << endl;
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
    if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iAtmPrecip)){
      _STORAGE<<" \""<<CStateVariable::GetStateVarLongName(_aStateVarType[i],_aStateVarLayer[i])<<"\"";}}
  _STORAGE<<" \"Total storage\" \"Cum. precip\" \"Cum. outflow\" \"MB error\""<<endl;

  _STORAGE<<"  :ColumnUnits mm/d mm/d mm mm ";
  for (i=0;i<GetNumStateVars();i++){
    if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iAtmPrecip)){
      _STORAGE<<" mm";}}
  _STORAGE<<" mm mm mm mm"<<endl;

  _STORAGE<<"  :ColumnType float float float float";
  for (i=0;i<GetNumStateVars();i++){
    if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iAtmPrecip)){
      _STORAGE<<" float";}}
  _STORAGE<<" float float float float"<<endl;

  _STORAGE << "  :ColumnFormat -1 -1 0 0";
  for (i = 0; i < GetNumStateVars(); i++){
    if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i != iAtmPrecip)){
      _STORAGE << " 0";
    }
  }
  _STORAGE << " 0 0 0 0" << endl;

  _STORAGE << ":EndColumnMetaData" << endl;

  _STORAGE << ":EndHeader" << endl;

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
    _HYDRO << ":CreationDate  " << GetCurrentTime() << endl;
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
  int         nSim, nbasins_dimid, varid_bsim;       // # of sub-basins with simulated outflow, dimension ID, and
  //                                                 // variable to write basin IDs for simulated outflows
  int         varid_qsim;                            // variable ID for simulated outflows
  int         varid_qobs;                            // variable ID for observed outflows
  int         varid_qin;                             // variable ID for observed inflows

  int         retval;                                // error value for NetCDF routines
  string      tmpFilename;
  int         ibasin, p;                             // loop over all sub-basins
  size_t      start[1], count[1];                    // determines where and how much will be written to NetCDF
  const char *current_basin_name[1];                 // current time in days since start time

  string      tmp,tmp2;
  static double fill_val[] = {NETCDF_BLANK_VALUE};
  static double miss_val[] = {NETCDF_BLANK_VALUE}; 

  // initialize all potential file IDs with -9 == "not existing and hence not opened" 
  _HYDRO_ncid    = -9;   // output file ID for Hydrographs.nc         (-9 --> not opened)
  _STORAGE_ncid  = -9;   // output file ID for WatershedStorage.nc    (-9 --> not opened)
  _FORCINGS_ncid = -9;   // output file ID for ForcingFunctions.nc    (-9 --> not opened)

  //converts start day into "days since YYYY-MM-DD HH:MM:SS"  (model start time)
  char  starttime[200]; // start time string in format 'days since YYY-MM-DD HH:MM:SS'
  JulianConvert( 0.0,Options.julian_start_day, Options.julian_start_year, tt);
  strcpy(starttime, "days since ") ;
  strcat(starttime, tt.date_string.c_str()) ;
  strcat(starttime, " 00:00:00");

  //====================================================================
  //  Hydrographs.nc
  //====================================================================
  // Create the file. 
  tmpFilename = FilenamePrepare("Hydrographs.nc", Options);
  retval = nc_create(tmpFilename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);  HandleNetCDFErrors(retval);
  _HYDRO_ncid = ncid;

  // ---------------------------------------------------------- 
  // time                                                       
  // ---------------------------------------------------------- 
  // (a) Define the DIMENSIONS. NetCDF will hand back an ID 
  retval = nc_def_dim(_HYDRO_ncid, "time", NC_UNLIMITED, &time_dimid);  HandleNetCDFErrors(retval);

  /// Define the time variable. Assign units attributes to the netCDF VARIABLES. 
  dimids1[0] = time_dimid;
  retval = nc_def_var(_HYDRO_ncid, "time", NC_DOUBLE, ndims1,dimids1, &varid_time); HandleNetCDFErrors(retval);
  retval = nc_put_att_text(_HYDRO_ncid, varid_time, "units"   , strlen(starttime)  , starttime);   HandleNetCDFErrors(retval);
  retval = nc_put_att_text(_HYDRO_ncid, varid_time, "calendar", strlen("gregorian"), "gregorian"); HandleNetCDFErrors(retval);

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
    tmp ="Name/ID of sub-basins with simulated outflows";
    tmp2="timeseries_ID";
    retval = nc_put_att_text(_HYDRO_ncid, varid_bsim, "long_name",  tmp.length(), tmp.c_str());    HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_HYDRO_ncid, varid_bsim, "cf_role"  , tmp2.length(),tmp2.c_str());    HandleNetCDFErrors(retval);
    
    varid_qsim= NetCDFAddMetadata2D(_HYDRO_ncid, time_dimid,nbasins_dimid,"q_sim","Simulated outflows","m**3 s**-1");
    varid_qobs= NetCDFAddMetadata2D(_HYDRO_ncid, time_dimid,nbasins_dimid,"q_obs","Observed outflows" ,"m**3 s**-1");
    varid_qin = NetCDFAddMetadata2D(_HYDRO_ncid, time_dimid,nbasins_dimid,"q_in" ,"Observed inflows"  ,"m**3 s**-1");
  }// end if nSim>0

  // End define mode. This tells netCDF we are done defining metadata. 
  retval = nc_enddef(_HYDRO_ncid);  HandleNetCDFErrors(retval);
    
  // write values to NetCDF 
  // (a) write gauged basin names/IDs to variable "basin_name" 
  ibasin = 0;
  for(p=0;p<_nSubBasins;p++){
    if(_pSubBasins[p]->IsGauged()  && (_pSubBasins[p]->IsEnabled())){
      if(_pSubBasins[p]->GetName()==""){ current_basin_name[0] = (to_string(_pSubBasins[p]->GetID())).c_str(); }
      else                             { current_basin_name[0] = (_pSubBasins[p]->GetName()).c_str(); }
      // write sub-basin name
      start[0] = ibasin;
      count[0] = 1;
      retval = nc_put_vara_string(_HYDRO_ncid,varid_bsim,start,count,&current_basin_name[0]);  HandleNetCDFErrors(retval);
      ibasin++;
    }
  }
  //====================================================================
  //  WatershedStorage.nc
  //====================================================================
  tmpFilename = FilenamePrepare("WatershedStorage.nc", Options);
  retval = nc_create(tmpFilename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);  HandleNetCDFErrors(retval);
  _STORAGE_ncid = ncid;

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
  varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"Channel Storage","Channel Storage","mm");
  varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"Reservoir Storage","Reservoir Storage","mm");
  varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"Rivulet Storage","Rivulet Storage","mm");

  int iAtmPrecip=GetStateVarIndex(ATMOS_PRECIP);
  for(int i=0;i<_nStateVars;i++){
    if((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iAtmPrecip)){
      string name =CStateVariable::GetStateVarLongName(_aStateVarType[i],_aStateVarLayer[i]);
      varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,name,name,"mm");
    }
  }
  varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"Total","total water storage","mm");
  varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"Cum. Input","cumulative water input","mm");
  varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"Cum. Outflow","cumulative water output","mm");
  varid= NetCDFAddMetadata(_STORAGE_ncid, time_dimid,"MB Error","mass balance error","mm");

  // End define mode. This tells netCDF we are done defining metadata. 
  retval = nc_enddef(_STORAGE_ncid);  HandleNetCDFErrors(retval);
    
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

  double snowfall    =GetAverageSnowfall();
  double precip      =GetAveragePrecip();
  double channel_stor=GetTotalChannelStorage();
  double reservoir_stor=GetTotalReservoirStorage();
  double rivulet_stor=GetTotalRivuletStorage();

  if ((tt.model_time==0) && (Options.suppressICs==true) && (Options.period_ending)){return;}

  //----------------------------------------------------------------
  // write watershed state variables  (WatershedStorage.tb0)
  if (tt.model_time!=0){_STORAGE<<" "<<precip-snowfall<<" "<<snowfall;}//precip
  else                 {_STORAGE<<" 0.0 0.0";}
  _STORAGE<<" "<<channel_stor+reservoir_stor<<" "<<rivulet_stor;
  //_STORAGE<<" "<<channel_stor<<" "<<reservoir_stor<<" "<<rivulet_stor;  // \todo [update] - backwards incompatible

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

  //----------------------------------------------------------------
  //Write hydrographs for gauged watersheds (ALWAYS DONE) (Hydrographs.tb0)
  if ((Options.ave_hydrograph) && (tt.model_time!=0))
  {
    _HYDRO<<" "<<GetAveragePrecip();
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
  size_t start1[1], count1[1];  // determines where and how much will be written to NetCDF; 1D variable (pre, time)
  size_t start2[2], count2[2];  // determines where and how much will be written to NetCDF; 2D variable (qsim, qobs, qin)
  double current_time[1];       // current time in days since start time
  double current_prec[1];       // precipitation of current time step

  current_time[0] = tt.model_time;

  //====================================================================
  //  Hydrographs.nc
  //====================================================================
  int    precip_id;             // variable id in NetCDF for precipitation
  int    qsim_id;               // variable id in NetCDF for simulated outflow
  int    qobs_id;               // variable id in NetCDF for observed outflow
  int    qin_id;                // variable id in NetCDF for observed inflow

  // (a) count how many values need to be written for q_obs, q_sim, q_in 
  int iSim, nSim; // current and total # of sub-basins with simulated outflows
  nSim = 0;
  for (int p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->IsGauged()  && (_pSubBasins[p]->IsEnabled())){nSim++;}
  }      

  // (b) allocate memory if necessary 
  double *outflow_obs=NULL;   // q_obs
  double *outflow_sim=NULL;   // q_sim
  double *inflow_obs =NULL;    // q_in
  if(nSim>0){
    outflow_sim=new double[nSim];
    outflow_obs=new double[nSim];
    inflow_obs =new double[nSim];
  }

  // (c) obtain data 
  iSim = 0;
  current_prec[0] = NETCDF_BLANK_VALUE;
  if ((Options.ave_hydrograph) && (current_time[0] != 0.0))
  {
    current_prec[0] = GetAveragePrecip();
    for (int p=0;p<_nSubBasins;p++)
    {
      if (_pSubBasins[p]->IsGauged() && (_pSubBasins[p]->IsEnabled())){
        outflow_sim[iSim] = _pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
        outflow_obs[iSim] = NETCDF_BLANK_VALUE;
        for (int i = 0; i < _nObservedTS; i++){
          if (IsContinuousFlowObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
          {
            double val = _pObservedTS[i]->GetAvgValue(current_time[0],Options.timestep); //time shift handled in CTimeSeries::Parse
            if ((val != RAV_BLANK_DATA) && (current_time[0]>0)){ outflow_obs[iSim] = val;    }
          }
        }
        inflow_obs[iSim] =NETCDF_BLANK_VALUE;
        if (_pSubBasins[p]->GetReservoir() != NULL){
          inflow_obs[iSim] = _pSubBasins[p]->GetIntegratedReservoirInflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
        }
        iSim++;
      }
    }      
  }
  else {  // point-value hydrograph
    if (current_time[0] != 0.0){current_prec[0] = GetAveragePrecip();} //watershed-wide precip
    else                       {current_prec[0] = NETCDF_BLANK_VALUE;} // was originally '---'
    for (int p=0;p<_nSubBasins;p++)
    {
      if (_pSubBasins[p]->IsGauged() && (_pSubBasins[p]->IsEnabled())){
        outflow_sim[iSim] = _pSubBasins[p]->GetOutflowRate();
        outflow_obs[iSim] = NETCDF_BLANK_VALUE;
        for (int i = 0; i < _nObservedTS; i++){
          if (IsContinuousFlowObs(_pObservedTS[i],_pSubBasins[p]->GetID()))
          {
            double val = _pObservedTS[i]->GetAvgValue(current_time[0],Options.timestep);
            if ((val != RAV_BLANK_DATA) && (current_time[0]>0)){ outflow_obs[iSim] = val;    }
          }
        }
        inflow_obs[iSim] =NETCDF_BLANK_VALUE;
        if (_pSubBasins[p]->GetReservoir() != NULL){
          inflow_obs[iSim] = _pSubBasins[p]->GetReservoirInflow();
        }
        iSim++;
      }
    }
  }

  // (d) write values to NetCDF (one time step) 
  start1[0] = int(round(current_time[0]/Options.timestep));   // element of NetCDF array that will be written
  count1[0] = 1;                                              // writes exactly one time step

  // write new time step 
  retval = nc_inq_varid (_HYDRO_ncid, "time",   &time_id);                             HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_HYDRO_ncid, time_id, start1, count1, &current_time[0]); HandleNetCDFErrors(retval);
  
  // write precipitation values 
  retval = nc_inq_varid (_HYDRO_ncid, "precip", &precip_id);                           HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_HYDRO_ncid,precip_id,start1,count1,&current_prec[0]);   HandleNetCDFErrors(retval);

  // write simulated outflow/obs outflow/obs inflow values 
  if (nSim > 0){
    start2[0] = int(round(current_time[0]/Options.timestep));   // element of NetCDF array that will be written
    start2[1] = 0;                                              // element of NetCDF array that will be written
    count2[0] = 1;      // writes exactly one time step
    count2[1] = nSim;   // writes exactly nSim elements
    retval = nc_inq_varid (_HYDRO_ncid, "q_sim", &qsim_id);                             HandleNetCDFErrors(retval);
    retval = nc_inq_varid (_HYDRO_ncid, "q_obs", &qobs_id);                             HandleNetCDFErrors(retval);
    retval = nc_inq_varid (_HYDRO_ncid, "q_in",  &qin_id);                              HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_HYDRO_ncid, qsim_id, start2, count2, &outflow_sim[0]); HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_HYDRO_ncid, qobs_id, start2, count2, &outflow_obs[0]); HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_HYDRO_ncid, qin_id, start2, count2, &inflow_obs[0]);   HandleNetCDFErrors(retval);
  }
  
  delete[] outflow_obs;
  delete[] outflow_sim;
  delete[] inflow_obs;

  //====================================================================
  //  WatershedStorage.nc
  //====================================================================
  double snowfall      =GetAverageSnowfall();
  double precip        =GetAveragePrecip();
  double channel_stor  =GetTotalChannelStorage();
  double reservoir_stor=GetTotalReservoirStorage();
  double rivulet_stor  =GetTotalRivuletStorage();
  double currentWater;
  int var_id;
  double value[1];

  // write new time step 
  retval = nc_inq_varid (_STORAGE_ncid, "time",   &time_id);                             HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_STORAGE_ncid, time_id, start1, count1, &current_time[0]); HandleNetCDFErrors(retval);
  

  if(tt.model_time!=0){
    value[0]=precip-snowfall;
    retval = nc_inq_varid (_STORAGE_ncid, "rainfall",   &var_id);                  HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_STORAGE_ncid, var_id, start1, count1, &value[0]); HandleNetCDFErrors(retval);
    value[0]=snowfall;
    retval = nc_inq_varid (_STORAGE_ncid, "snowfall",   &var_id);                  HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_STORAGE_ncid, var_id, start1, count1, &value[0]); HandleNetCDFErrors(retval);
  }
  value[0]=channel_stor;
  retval = nc_inq_varid (_STORAGE_ncid, "Channel Storage",   &var_id);           HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_STORAGE_ncid, var_id, start1, count1, &value[0]); HandleNetCDFErrors(retval);
  value[0]=reservoir_stor;
  retval = nc_inq_varid (_STORAGE_ncid, "Reservoir Storage",   &var_id);         HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_STORAGE_ncid, var_id, start1, count1, &value[0]); HandleNetCDFErrors(retval);
  value[0]=rivulet_stor;
  retval = nc_inq_varid (_STORAGE_ncid, "Rivulet Storage",   &var_id);           HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_STORAGE_ncid, var_id, start1, count1, &value[0]); HandleNetCDFErrors(retval);

  currentWater=0.0;
  double S;string short_name;
  int iAtmPrecip=GetStateVarIndex(ATMOS_PRECIP);
  for (int i=0;i<GetNumStateVars();i++)
  {
    if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) && (i!=iAtmPrecip))
    {
      S=FormatDouble(GetAvgStateVar(i));
      value[0]=S;
      short_name=CStateVariable::GetStateVarLongName(_aStateVarType[i],_aStateVarLayer[i]);
      retval = nc_inq_varid (_STORAGE_ncid, short_name.c_str(),   &var_id);          HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_STORAGE_ncid, var_id, start1, count1, &value[0]); HandleNetCDFErrors(retval);
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
  value[0]=currentWater;
  retval = nc_inq_varid (_STORAGE_ncid, "Total",   &var_id);                      HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_STORAGE_ncid, var_id, start1, count1, &value[0]);  HandleNetCDFErrors(retval);
  value[0]=_CumulInput;
  retval = nc_inq_varid (_STORAGE_ncid, "Cum. Input",   &var_id);                 HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_STORAGE_ncid, var_id, start1, count1, &value[0]);  HandleNetCDFErrors(retval);
  value[0]=_CumulOutput;
  retval = nc_inq_varid (_STORAGE_ncid, "Cum. Outflow",   &var_id);               HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_STORAGE_ncid, var_id, start1, count1, &value[0]);  HandleNetCDFErrors(retval);
  value[0]=FormatDouble((currentWater-_initWater)+(_CumulOutput-_CumulInput));
  retval = nc_inq_varid (_STORAGE_ncid, "MB Error",   &var_id);                   HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_STORAGE_ncid, var_id, start1, count1, &value[0]);  HandleNetCDFErrors(retval);

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
#endif
  }
  g_output_directory=Options.output_dir;//necessary evil
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
/// \brief returns directory path given filename and relative path
///
/// \param filename [in] filename, e.g., C:/temp/thisfile.txt returns c:/temp
/// \param relfile [in] filename of reference file 
/// e.g.,
///       absolute path of reference file is adopted 
///       if filename = something.txt         and relfile= c:/temp/myfile.rvi,  returns c:/temp/something.txt
///
///       relative path of reference file is adopted
///       if filename = something.txt         and relfile= ../dir/myfile.rvi,   returns ../dir/myfile.rvi
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
  
  if (StringToUppercase(filename).find(StringToUppercase(filedir)) == string::npos){ //checks to see if absolute dir already included in redirect filename
    
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
int NetCDFAddMetadata2D(const int fileid,const int time_dimid,int nbasins_dimid,string shortname,string longname,string units)
{
  int varid(0);
#ifdef _RVNETCDF_
  int retval;
  int dimids2[2];

  static double fill_val[] = {NETCDF_BLANK_VALUE};
  static double miss_val[] = {NETCDF_BLANK_VALUE}; 
  
  dimids2[0] = time_dimid;
  dimids2[1] = nbasins_dimid;

  // (a) create variable 
  retval = nc_def_var(fileid,shortname.c_str(),NC_DOUBLE,2,dimids2,&varid); HandleNetCDFErrors(retval);

  // (b) add attributes to variable
  retval = nc_put_att_text  (fileid,varid,"units",units.length(),units.c_str());              HandleNetCDFErrors(retval);
  retval = nc_put_att_text  (fileid,varid,"long_name",longname.length(),longname.c_str());    HandleNetCDFErrors(retval);
  retval = nc_put_att_double(fileid,varid,"_FillValue",NC_DOUBLE,1,fill_val);                 HandleNetCDFErrors(retval);
  retval = nc_put_att_double(fileid,varid,"missing_value",NC_DOUBLE,1,miss_val);              HandleNetCDFErrors(retval);
#endif
  return varid;
}
