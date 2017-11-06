/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"
#include "StateVariables.h"

#if defined(_WIN32)
#include <direct.h>
#elif defined(__linux__)
#include <sys/stat.h>
#endif

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

#ifdef _RVNETCDF_
  
  /* close netcdfs */
  int    retval;      // error value for NetCDF routines
  if (_HYDRO_ncid != -9)    { if ((retval = nc_close(_HYDRO_ncid)))    { nc_error_exit(retval); }; }
  _HYDRO_ncid    = -9;
  if (_STORAGE_ncid != -9)  { if ((retval = nc_close(_STORAGE_ncid)))  { nc_error_exit(retval); }; }
  _STORAGE_ncid  = -9;
  if (_FORCINGS_ncid != -9) { if ((retval = nc_close(_FORCINGS_ncid))) { nc_error_exit(retval); }; }
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
      if (_pSubBasins[p]->IsGauged()){
        string name;
        if (_pSubBasins[p]->GetName()==""){_HYDRO<<",ID="<<_pSubBasins[p]->GetID()  <<" [m3/s]";}
        else                              {_HYDRO<<","   <<_pSubBasins[p]->GetName()<<" [m3/s]";}
        //if (Options.print_obs_hydro)
        {
          for (int i = 0; i < _nObservedTS; i++){
            if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "HYDROGRAPH")) &&
                (s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
                (_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular))
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
    ofstream RES_STAGE;
    tmpFilename=FilenamePrepare("ReservoirStages.csv",Options);
    RES_STAGE.open(tmpFilename.c_str());
    if (RES_STAGE.fail()){
      ExitGracefully(("CModel::WriteOutputFileHeaders: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }

    RES_STAGE<<"time,date,hour";
    RES_STAGE<<",precip [mm/day]";
    for (p=0;p<_nSubBasins;p++){
      if ((_pSubBasins[p]->IsGauged()) && (_pSubBasins[p]->GetReservoir()!=NULL)) {
        string name;
        if (_pSubBasins[p]->GetName()==""){RES_STAGE<<",ID="<<_pSubBasins[p]->GetID()  <<" [m]";}
        else                              {RES_STAGE<<","   <<_pSubBasins[p]->GetName()<<" [m]";}
      }
      //if (Options.print_obs_hydro)
      {
        for (int i = 0; i < _nObservedTS; i++){
          if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "RESERVOIR_STAGE")) &&
              (s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
              (_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular))
          {
            if (_pSubBasins[p]->GetName()==""){RES_STAGE<<",ID="<<_pSubBasins[p]->GetID()  <<" (observed) [m3/s]";}
            else                              {RES_STAGE<<","   <<_pSubBasins[p]->GetName()<<" (observed) [m3/s]";}
          }
        }
      }
    }
    RES_STAGE<<endl;
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
              if (!first){MB<<",";}
              first=false;
            }
          }
        }
        for (j=0;j<_nProcesses;j++){
          for (int q=0;q<_pProcesses[j]->GetNumConnections();q++){
            if (_pProcesses[j]->GetToIndices()[q]==i){
              if (!first){MB<<",";}
              first=false;
            }
          }
        }
        MB<<",,,";
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
    _FORCINGS<<" ET radiation [W/m2], SW radiation [W/m2], LW radiation [W/m2], wind vel. [m/s],";
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

    // Console output
    //----------------------------------------------------------------
    if ((quiet) && (!Options.silent) && (tt.day_of_month==1) && ((tt.julian_day)-floor(tt.julian_day)<Options.timestep))
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
          _STORAGE<<","<<S;
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

      _STORAGE<<","<<currentWater<<","<<_CumulInput<<","<<_CumulOutput<<","<<(currentWater-_initWater)+(_CumulOutput-_CumulInput);
      _STORAGE<<endl;

      //Write hydrographs for gauged watersheds (ALWAYS DONE)
      //----------------------------------------------------------------
      if ((Options.ave_hydrograph) && (t!=0.0))
      {
        _HYDRO<<t<<","<<thisdate<<","<<thishour<<","<<GetAveragePrecip();
        for (int p=0;p<_nSubBasins;p++){
          if (_pSubBasins[p]->IsGauged())
          {
            _HYDRO<<","<<_pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);

            //if (Options.print_obs_hydro)
            {
              for (int i = 0; i < _nObservedTS; i++)
              {
                if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "HYDROGRAPH")) &&
                    (s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
                    (_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular))
                {

                  //int nn=int(floor((tt.model_time+TIME_CORRECTION)/ Options.timestep));
                  //double val=_pObservedTS[i]->GetSampledValue(nn); //fails for interval>timestep
                  double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep); //time shift handled in CTimeSeries::Parse
                  if ((val != CTimeSeriesABC::BLANK_DATA) && (tt.model_time>0)){ _HYDRO << "," << val; }
                  else                                                         { _HYDRO << ", ";       }
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
      else //point value hydrograph
      {
        _HYDRO<<t<<","<<thisdate<<","<<thishour;
        if (t!=0){_HYDRO<<","<<GetAveragePrecip();}//watershed-wide precip
        else     {_HYDRO<<",---";}
        for (int p=0;p<_nSubBasins;p++){
          if (_pSubBasins[p]->IsGauged())
          {
            _HYDRO<<","<<_pSubBasins[p]->GetOutflowRate();

            //if (Options.print_obs_hydro)
            {
              for (int i = 0; i < _nObservedTS; i++){
                if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "HYDROGRAPH")) &&
                    (s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
                    (_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular))
                {
                  //int nn=int(floor((tt.model_time+TIME_CORRECTION)/ Options.timestep));
                  //double val=_pObservedTS[i]->GetSampledValue(nn); //fails for interval>timestep
                  double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep);
                  if ((val != CTimeSeriesABC::BLANK_DATA) && (tt.model_time>0)){ _HYDRO << "," << val; }
                  else                                                         { _HYDRO << ", ";       }
                }
              }
            }
            if (_pSubBasins[p]->GetReservoir() != NULL){
              _HYDRO<<","<<_pSubBasins[p]->GetReservoirInflow();
            }
          }
        }
        _HYDRO<<endl;
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
      double sum;
      int kk=Options.write_group_mb;
      ofstream HGMB;
      tmpFilename=FilenamePrepare(_pHRUGroups[kk]->GetName()+"_MassEnergyBalance.csv",Options);
      HGMB.open(tmpFilename.c_str(),ios::app);
      HGMB<<t<<","<<thisdate<<","<<thishour;
      double areasum=0.0;
      for(k = 0; k < _nHydroUnits; k++){
        if(_pHRUGroups[kk]->IsInGroup(k)){
          areasum+=_pHydroUnits[k]->GetArea();
        }
      }
      for (int js=0;js<_nTotalConnections;js++)
      {
        sum=0.0;
        for (k = 0; k < _nHydroUnits; k++){
          if (_pHRUGroups[kk]->IsInGroup(k)){
            sum += _aCumulativeBal[k][js] * _pHydroUnits[k]->GetArea();
          }
        }
        HGMB<<","<<sum/areasum;
      }
      HGMB<<endl;
      HGMB.close();
    }

    //Write cumulative mass balance info to WatershedMassEnergyBalance.csv
    //----------------------------------------------------------------
    if (Options.write_mass_bal)
    {
      double sum;
      ofstream MB;
      tmpFilename=FilenamePrepare("WatershedMassEnergyBalance.csv",Options);
      MB.open(tmpFilename.c_str(),ios::app);

      MB<<t<<","<<thisdate<<","<<thishour;
      for (int js=0;js<_nTotalConnections;js++)
      {
        sum=0.0;
        for (k=0;k<_nHydroUnits;k++){
          sum+=_aCumulativeBal[k][js]*_pHydroUnits[k]->GetArea();
        }
        MB<<","<<sum/_WatershedArea;
      }
      MB<<endl;
      MB.close();
    }

    //ReservoirStages.csv
    //--------------------------------------------------------------
    if ((Options.write_reservoir) && (Options.output_format!=OUTPUT_NONE))
    {
      ofstream RES_STAGE;
      tmpFilename=FilenamePrepare("ReservoirStages.csv",Options);
      RES_STAGE.open(tmpFilename.c_str(),ios::app);

      RES_STAGE<< t<<","<<thisdate<<","<<thishour<<","<<GetAveragePrecip();
      for (int p=0;p<_nSubBasins;p++){
        if ((_pSubBasins[p]->IsGauged()) && (_pSubBasins[p]->GetReservoir()!=NULL)) {
          RES_STAGE<<","<<_pSubBasins[p]->GetReservoir()->GetStage();
        }
        //if (Options.print_obs_hydro)
        {
          for (int i = 0; i < _nObservedTS; i++){
            if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "RESERVOIR_STAGE")) &&
                (s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
                (_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular))
            {
              //int nn=int(floor((tt.model_time+TIME_CORRECTION)/ Options.timestep));
              double val = _pObservedTS[i]->GetAvgValue(tt.model_time,Options.timestep);
              if ((val != CTimeSeriesABC::BLANK_DATA) && (tt.model_time>0)){ RES_STAGE << "," << val; }
              else                                                         { RES_STAGE << ", ";       }
            }
          }
        }

      }
      RES_STAGE<<endl;
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

      int j,js,q;
      double cumsum;

      ofstream MB;
      tmpFilename=FilenamePrepare("ExhaustiveMassBalance.csv",Options);
      MB.open(tmpFilename.c_str(),ios::app);

      MB<<t<<","<<thisdate<<","<<thishour;
      for (int i=0;i<_nStateVars;i++)
      {
        if (CStateVariable::IsWaterStorage(_aStateVarType[i]))
        {
          cumsum=0.0;
          js=0;
          for (j=0;j<_nProcesses;j++){
            for (q=0;q<_pProcesses[j]->GetNumConnections();q++){
              if (_pProcesses[j]->GetFromIndices()[q]==i)
              {
                sum=0.0;
                for (k=0;k<_nHydroUnits;k++){
                  sum+=_aCumulativeBal[k][js]*_pHydroUnits[k]->GetArea();
                }
                MB<<","<<-sum/_WatershedArea;
                cumsum-=sum/_WatershedArea;
              }
              js++;
            }
          }
          js=0;
          for (j=0;j<_nProcesses;j++){
            for (q=0;q<_pProcesses[j]->GetNumConnections();q++){
              if (_pProcesses[j]->GetToIndices()[q]==i)
              {
                sum=0.0;
                for (k=0;k<_nHydroUnits;k++){
                  sum+=_aCumulativeBal[k][js]*_pHydroUnits[k]->GetArea();
                }
                MB<<","<<sum/_WatershedArea;
                cumsum+=sum/_WatershedArea;
              }
              js++;
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

    //Write wshed-averaged forcing functions to ForcingFunctions.csv
    //----------------------------------------------------------------
    if (Options.write_forcings)
    {
      force_struct *pFave;
      force_struct faveStruct = GetAverageForcings();
      pFave = &faveStruct;
      _FORCINGS<<t<<","<<thisdate<<","<<thishour<<","<<pFave->day_angle<<",";
      _FORCINGS<<pFave->precip*(1-pFave->snow_frac) <<","<<pFave->precip*(pFave->snow_frac) <<",";
      _FORCINGS<<pFave->temp_ave<<","<<pFave->temp_daily_min<<","<<pFave->temp_daily_max<<","<<pFave->temp_daily_ave<<","<<pFave->temp_month_min<<","<<pFave->temp_month_max<<",";
      _FORCINGS<<pFave->air_dens<<","<<pFave->air_pres<<","<<pFave->rel_humidity<<",";
      _FORCINGS<<pFave->cloud_cover<<",";
      _FORCINGS<<pFave->ET_radia*MJ_PER_D_TO_WATT<<","<<pFave->SW_radia*MJ_PER_D_TO_WATT<<","<<pFave->LW_radia*MJ_PER_D_TO_WATT<<","<<pFave->wind_vel<<",";
      _FORCINGS<<pFave->PET<<","<<pFave->OW_PET<<",";
      _FORCINGS<<pFave->subdaily_corr<<","<<pFave->potential_melt;
      _FORCINGS<<endl;
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

  // parameters.csv - parameters file
  if ((Options.write_params) && (final==true))
  {
    ofstream PARAMS;
    tmpFilename=FilenamePrepare("parameters.csv",Options);
    PARAMS.open(tmpFilename.c_str());
    if (PARAMS.fail()){
      ExitGracefully(("CModel::WriteMajorOutput: Unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    CGlobalParams   ::WriteParamsToFile(PARAMS);
    CSoilClass      ::WriteParamsToFile(PARAMS);
    CVegetationClass::WriteParamsToFile(PARAMS);
    CLandUseClass   ::WriteParamsToFile(PARAMS);
    CTerrainClass   ::WriteParamsToFile(PARAMS);
    PARAMS.close();
  }

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
  if (!Options.silent){
    cout <<"==MODEL SUMMARY======================================="<<endl;
    cout <<"       Model Run: "<<Options.run_name    <<endl;
    cout <<"     # SubBasins: "<<GetNumSubBasins()   << " ("<< rescount << " reservoirs)"<<endl;
    cout <<"          # HRUs: "<<GetNumHRUs()        <<endl;
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
    cout <<"  Watershed Area: "<<_WatershedArea              <<" km2"<<endl;
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
  _STORAGE << ":Version       " << Options.version << endl;
  _STORAGE << ":CreationDate  " << GetCurrentTime() << endl;
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
  _HYDRO << ":Version       " << Options.version << endl;
  _HYDRO << ":CreationDate  " << GetCurrentTime() << endl;
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
//
void CModel::WriteNetcdfStandardHeaders(const optStruct &Options)
{

#ifdef _RVNETCDF_

  time_struct tt;                                    // start time structure
  char        starttime[200];                        // start time string in format 'days since YYY-MM-DD HH:MM:SS'
  int         ndims1 = 1;
  int         ndims2 = 2;
  int         dimids1[ndims1];                       // array which will contain all dimension ids for a variable
  int         dimids2[ndims2];                       // array which will contain all dimension ids for a variable
  int         ncid, varid_pre;                       // When we create netCDF variables and dimensions, we get back an ID for each one.
  int         time_dimid, varid_time;                // dimension ID (holds number of time steps) and variable ID (holds time values) for time
  int         Nsim, nbasin_sim_dimid, varid_bsim;    // # of sub-basins with simulated outflow, dimension ID, and
  //                                                 // variable to write basin IDs for simulated outflows
  int         Nobs, nbasin_obs_dimid, varid_bobs;    // # of sub-basins with observed outflow, dimension ID, and
  //                                                 //  variable to write basin IDs for observed outflows
  int         Nres, nbasin_res_dimid, varid_bres;    // # of sub-basins with observed (reservoir) inflows, dimension ID, and
  //                                                 //  variable to write basin IDs for observed (reservoir) inflows
  int         varid_qsim;                            // variable ID for simulated outflows
  int         varid_qobs;                            // variable ID for observed outflows
  int         varid_qin;                             // variable ID for observed inflows

  int         retval;                                // error value for NetCDF routines
  string      tmpFilename;
  int         ibasin, p;                             // loop over all sub-basins
  size_t      start[1], count[1];                    // determines where and how much will be written to NetCDF
  const char *current_basin_name[1];                 // current time in days since start time

  /* initialize all potential file IDs with -9 == "not existing and hence not opened" */
  _HYDRO_ncid    = -9;   // output file ID for Hydrographs.nc         (-9 --> not opened)
  _STORAGE_ncid  = -9;   // output file ID for WatershedStorage.nc    (-9 --> not opened)
  _FORCINGS_ncid = -9;   // output file ID for ForcingFunctions.nc    (-9 --> not opened)
  
  tmpFilename = FilenamePrepare("Hydrographs.nc", Options);

  /* Create the file. The NC_CLOBBER parameter tells netCDF to
   * overwrite this file, if it already exists. NC_NETCDF4 tells 
   * that NetCDF version 4 file is created. */
  if ((retval = nc_create(tmpFilename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid))){
    nc_error_exit(retval);
  }
  _HYDRO_ncid = ncid;

  /* ---------------------------------------------------------- */
  /* time                                                       */
  /* ---------------------------------------------------------- */
  /* (a) Define the DIMENSIONS. NetCDF will hand back an ID for each. */
  if ((retval = nc_def_dim(_HYDRO_ncid, "time", NC_UNLIMITED, &time_dimid))){
    nc_error_exit(retval);
  }

  /* (b) The dimids array is used to pass the IDs of the dimensions of
   *     the variable. */
  dimids1[0] = time_dimid;
  
  /* (c) Define the VARIABLES. The type of the variable can be: 
   * 	 NC_BYTE   	   signed 1 byte integer, 
   * 	 NC_CHAR   	   ISO/ASCII character, 
   * 	 NC_SHORT  	   signed 2 byte integer, 
   * 	 NC_INT    	   signed 4-byte integer,  
   * 	 NC_FLOAT  	   single precision floating point number,
   * 	 NC_DOUBLE 	   double precision floating point number,
   * 	 nc_UBYTE  	   unsigned 1 byte int,
   * 	 nc_USHORT 	   unsigned 2-byte int,
   * 	 nc_UINT   	   unsigned 4-byte int,
   * 	 nc_INT64  	   signed 8-byte int,
   * 	 nc_UINT64 	   unsigned 8-byte int, or
   * 	 nc_STRING 	   string */
  if ((retval = nc_def_var(_HYDRO_ncid, "time", NC_DOUBLE, ndims1, 
			   dimids1, &varid_time))){
    nc_error_exit(retval);
  }

  /* (d) Assign units attributes to the netCDF VARIABLES. */
  /*     --> converts start day into "days since YYYY-MM-DD HH:MM:SS" */
  JulianConvert( 0.0,Options.julian_start_day, Options.julian_start_year, tt);
  strcpy(starttime, "days since ") ;
  strcat(starttime, tt.date_string.c_str()) ;
  strcat(starttime, " 00:00:00");
  if ((retval = nc_put_att_text(_HYDRO_ncid, varid_time, "units", 
				strlen(starttime), starttime))){
    nc_error_exit(retval);
  }
  if ((retval = nc_put_att_text(_HYDRO_ncid, varid_time, "calendar", 
				strlen("gregorian"), "gregorian"))){
    nc_error_exit(retval);
  }

  /* ---------------------------------------------------------- */
  /* precipitation                                              */
  /* ---------------------------------------------------------- */
  /* (a) create variable precipitation */
  if ((retval = nc_def_var(_HYDRO_ncid, "precip", NC_DOUBLE, ndims1, 
			   dimids1, &varid_pre))){
    nc_error_exit(retval);
  }
  /* (b) add attributes to variable */
  if ((retval = nc_put_att_text(_HYDRO_ncid, varid_pre, "units", 
				strlen("mm d**-1"), "mm d**-1"))){
    nc_error_exit(retval);
  }
  if ((retval = nc_put_att_text(_HYDRO_ncid, varid_pre, "long_name", 
				strlen("Precipitation"), "Precipitation"))){
    nc_error_exit(retval);
  }
  static double fill_val[] = {-9999.0}; /* attribute vals */
  if ((retval = nc_put_att_double(_HYDRO_ncid, varid_pre, "_FillValue", NC_DOUBLE,
				  1, fill_val))){
    nc_error_exit(retval);
  }

  static double miss_val[] = {-9999.0}; /* attribute vals */
  if ((retval = nc_put_att_double(_HYDRO_ncid, varid_pre, "missing_value", NC_DOUBLE,
				  1, miss_val))){
    nc_error_exit(retval);
  }

  /* ---------------------------------------------------------- */
  /* simulated outflows                                         */
  /* ---------------------------------------------------------- */
  /* (a) count number of simulated outflows "Nsim"              */
  Nsim = 0;
  for (p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->IsGauged()){
      Nsim = Nsim + 1;
    }
  }
  /* (b) create dimension "nbasin_sim"                          */
  if (Nsim > 0){
    if ((retval = nc_def_dim(_HYDRO_ncid, "nbasin_sim", Nsim, &nbasin_sim_dimid))){
      nc_error_exit(retval);
    }
  }
  /* (c) create variable "basin_name_sim"                       */
  if (Nsim > 0){
    dimids1[0] = nbasin_sim_dimid;
    if ((retval = nc_def_var(_HYDRO_ncid, "basin_name_sim", NC_STRING, ndims1, 
			     dimids1, &varid_bsim))){
      nc_error_exit(retval);
    }
  }
  /* (d) set some attributes to variable "basin_name_sim"       */
  if (Nsim > 0){
    if ((retval = nc_put_att_text(_HYDRO_ncid, varid_bsim, "long_name", 
				  strlen("Name/ID of sub-basins with simulated outflows"),
				  "Name/ID of sub-basins with simulated outflows"))){
      nc_error_exit(retval);
    }
  }
  /* (e) create variable "q_sim" which will contain at the end simulated outflows; shape=(tt,Nsim) */
  if (Nsim > 0){
    dimids2[0] = time_dimid;
    dimids2[1] = nbasin_sim_dimid;
    if ((retval = nc_def_var(_HYDRO_ncid, "q_sim", NC_DOUBLE, ndims2, 
			     dimids2, &varid_qsim))){
      nc_error_exit(retval);
    }
  }
  /* (f) set some attributes to variable "q_sim"       */
  if (Nsim > 0){
    if ((retval = nc_put_att_text(_HYDRO_ncid, varid_qsim, "long_name", 
				  strlen("Simulated outflows"),
				  "Simulated outflows"))){
      nc_error_exit(retval);
    }
    if ((retval = nc_put_att_text(_HYDRO_ncid, varid_qsim, "units", 
				  strlen("m**3 s**-1"),
				  "m**3 s**-1"))){
      nc_error_exit(retval);
    }
    static double fill_val[] = {-9999.0}; /* attribute vals */
    if ((retval = nc_put_att_double(_HYDRO_ncid, varid_qsim, "_FillValue", NC_DOUBLE,
				    1, fill_val))){
      nc_error_exit(retval);
    }
    
    static double miss_val[] = {-9999.0}; /* attribute vals */
    if ((retval = nc_put_att_double(_HYDRO_ncid, varid_qsim, "missing_value", NC_DOUBLE,
				    1, miss_val))){
      nc_error_exit(retval);
    }
  }
  
  /* ---------------------------------------------------------- */
  /* observed outflows                                          */
  /* ---------------------------------------------------------- */
  /* (a) count number of simulated outflows "Nobs  "            */
  Nobs = 0;
  for (p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->IsGauged()){
        {
          for (int i = 0; i < _nObservedTS; i++){
            if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "HYDROGRAPH")) &&
                (s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
                (_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular))
            {
              Nobs = Nobs + 1;
            }
          }
        }
      }
  }
  /* (b) create dimension "nbasin_obs"                          */
  if (Nobs > 0){
    if ((retval = nc_def_dim(_HYDRO_ncid, "nbasin_obs", Nobs, &nbasin_obs_dimid))){
      nc_error_exit(retval);
    }
  }
  /* (c) create variable "basin_name_obs"                       */
  if (Nobs > 0){
    dimids1[0] = nbasin_obs_dimid;
    if ((retval = nc_def_var(_HYDRO_ncid, "basin_name_obs", NC_STRING, ndims1, 
			     dimids1, &varid_bobs))){
      nc_error_exit(retval);
    }
  }
  /* (d) set some attributes to variable "basin_name_obs"       */
  if (Nobs > 0){
    if ((retval = nc_put_att_text(_HYDRO_ncid, varid_bobs, "long_name", 
				  strlen("Name/ID of sub-basins with observed outflows"),
				  "Name/ID of sub-basins with observed outflows"))){
      nc_error_exit(retval);
    }
  }
  /* (e) create variable "q_obs" which will contain at the end observed outflows; shape=(tt,Nobs) */
  if (Nobs > 0){
    dimids2[0] = time_dimid;
    dimids2[1] = nbasin_obs_dimid;
    if ((retval = nc_def_var(_HYDRO_ncid, "q_obs", NC_DOUBLE, ndims2, 
			     dimids2, &varid_qobs))){
      nc_error_exit(retval);
    }
  }
  /* (f) set some attributes to variable "q_obs"       */
  if (Nobs > 0){
    if ((retval = nc_put_att_text(_HYDRO_ncid, varid_qobs, "long_name", 
				  strlen("Observed outflows"),
				  "Observed outflows"))){
      nc_error_exit(retval);
    }
    if ((retval = nc_put_att_text(_HYDRO_ncid, varid_qobs, "units", 
				  strlen("m**3 s**-1"),
				  "m**3 s**-1"))){
      nc_error_exit(retval);
    }
    static double fill_val[] = {-9999.0}; /* attribute vals */
    if ((retval = nc_put_att_double(_HYDRO_ncid, varid_qobs, "_FillValue", NC_DOUBLE,
				    1, fill_val))){
      nc_error_exit(retval);
    }
    
    static double miss_val[] = {-9999.0}; /* attribute vals */
    if ((retval = nc_put_att_double(_HYDRO_ncid, varid_qobs, "missing_value", NC_DOUBLE,
				    1, miss_val))){
      nc_error_exit(retval);
    }
  }

  /* ---------------------------------------------------------- */
  /* inflows                                                    */
  /* ---------------------------------------------------------- */
  /* (a) count number of simulated outflows "Nres  "            */
  Nres = 0;
  for (p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->IsGauged()){
      if (_pSubBasins[p]->GetReservoir() != NULL){
	Nres = Nres + 1;
      }
    }
  }
  /* (b) create dimension "nbasin_res"                          */
  if (Nres > 0){
    if ((retval = nc_def_dim(_HYDRO_ncid, "nbasin_res", Nres, &nbasin_res_dimid))){
      nc_error_exit(retval);
    }
  }
  /* (c) create variable "basin_name_res"                       */
  if (Nres > 0){
    dimids1[0] = nbasin_res_dimid;
    if ((retval = nc_def_var(_HYDRO_ncid, "basin_name_res", NC_STRING, ndims1, 
			     dimids1, &varid_bres))){
      nc_error_exit(retval);
    }
  }
  /* (d) set some attributes to variable "basin_name_res"       */
  if (Nres > 0){
    if ((retval = nc_put_att_text(_HYDRO_ncid, varid_bres, "long_name", 
				  strlen("Name/ID of sub-basins with inflows"),
				  "Name/ID of sub-basins with inflows"))){
      nc_error_exit(retval);
    }
  }
  /* (e) create variable "q_in" which will contain at the end observed inflows; shape=(tt,Nres) */
  if (Nres > 0){
    dimids2[0] = time_dimid;
    dimids2[1] = nbasin_res_dimid;
    if ((retval = nc_def_var(_HYDRO_ncid, "q_in", NC_DOUBLE, ndims2, 
			     dimids2, &varid_qin))){
      nc_error_exit(retval);
    }
  }
  /* (f) set some attributes to variable "q_in"       */
  if (Nres > 0){
    if ((retval = nc_put_att_text(_HYDRO_ncid, varid_qin, "long_name", 
				  strlen("Observed outflows"),
				  "Observed outflows"))){
      nc_error_exit(retval);
    }
    if ((retval = nc_put_att_text(_HYDRO_ncid, varid_qin, "units", 
				  strlen("m**3 s**-1"),
				  "m**3 s**-1"))){
      nc_error_exit(retval);
    }
    static double fill_val[] = {-9999.0}; /* attribute vals */
    if ((retval = nc_put_att_double(_HYDRO_ncid, varid_qin, "_FillValue", NC_DOUBLE,
				    1, fill_val))){
      nc_error_exit(retval);
    }
    
    static double miss_val[] = {-9999.0}; /* attribute vals */
    if ((retval = nc_put_att_double(_HYDRO_ncid, varid_qin, "missing_value", NC_DOUBLE,
				    1, miss_val))){
      nc_error_exit(retval);
    }
  }

  /* End define mode. This tells netCDF we are done defining
   * metadata. */
  if ((retval = nc_enddef(_HYDRO_ncid)))
    nc_error_exit(retval);

  /* write values to NetCDF */
  /* (a) write basin names/IDs to variable "basin_name_sim"   */
  ibasin = 0;
  for (p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->IsGauged()){
      if (_pSubBasins[p]->GetName()==""){current_basin_name[0] = (to_string(_pSubBasins[p]->GetID())).c_str();}
      else                              {current_basin_name[0] = (_pSubBasins[p]->GetName()).c_str();}
      // write sub-basin name
      start[0] = ibasin;
      count[0] = 1;
      if ((retval = nc_put_vara_string(_HYDRO_ncid, varid_bsim, start, count, &current_basin_name[0]))){
	nc_error_exit(retval);
      }
      ibasin = ibasin + 1;
    }
  }
  /* (b) write basin names/IDs to variable "basin_name_obs"   */
  ibasin = 0;
  for (p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->IsGauged()){
        {
          for (int i = 0; i < _nObservedTS; i++){
            if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "HYDROGRAPH")) &&
                (s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
                (_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular))
            {
              if (_pSubBasins[p]->GetName()==""){current_basin_name[0] = (to_string(_pSubBasins[p]->GetID())).c_str();}
              else                              {current_basin_name[0] = (_pSubBasins[p]->GetName()).c_str();}
	      // write sub-basin name
	      start[0] = ibasin;
	      count[0] = 1;
	      if ((retval = nc_put_vara_string(_HYDRO_ncid, varid_bobs, start, count, &current_basin_name[0]))){
		nc_error_exit(retval);
	      }
	      ibasin = ibasin + 1;
            }
          }
        }
      }
  }
  /* (c) write basin names/IDs to variable "basin_name_res"   */
  ibasin = 0;
  for (p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->IsGauged()){
      if (_pSubBasins[p]->GetReservoir() != NULL){
	if (_pSubBasins[p]->GetName()==""){current_basin_name[0] = (to_string(_pSubBasins[p]->GetID())).c_str();}
	else                              {current_basin_name[0] = (_pSubBasins[p]->GetName()).c_str();}
	// write sub-basin name
	start[0] = ibasin;
	count[0] = 1;
	if ((retval = nc_put_vara_string(_HYDRO_ncid, varid_bres, start, count, &current_basin_name[0]))){
	  nc_error_exit(retval);
	}
	ibasin = ibasin + 1;
      }
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

  double snowfall    =GetAverageSnowfall();
  double precip      =GetAveragePrecip();
  double channel_stor=GetTotalChannelStorage();
  double reservoir_stor=GetTotalReservoirStorage();
  double rivulet_stor=GetTotalRivuletStorage();

  if ((tt.model_time==0) && (Options.suppressICs==true) && (Options.period_ending)){return;}

  //----------------------------------------------------------------
  // write watershed state variables
  if (tt.model_time!=0){_STORAGE<<" "<<precip-snowfall<<" "<<snowfall;}//precip
  else                 {_STORAGE<<" 0.0 0.0";}
  _STORAGE<<" "<<channel_stor+reservoir_stor<<" "<<rivulet_stor;
  //_STORAGE<<" "<<channel_stor<<" "<<reservoir_stor<<" "<<rivulet_stor;  // \todo [update] - backwards incompatible

  currentWater=0.0;
  for (i=0;i<GetNumStateVars();i++){
    if ((CStateVariable::IsWaterStorage(_aStateVarType[i])) &&  (i!=iCumPrecip)){
      S=GetAvgStateVar(i);_STORAGE<<" "<<S;currentWater+=S;
    }
  }
  currentWater+=channel_stor+rivulet_stor;

  _STORAGE<<" "<<currentWater<<" "<<_CumulInput<<" "<<_CumulOutput<<" "<<(currentWater-_initWater)+(_CumulOutput-_CumulInput);
  _STORAGE<<endl;

  //----------------------------------------------------------------
  //Write hydrographs for gauged watersheds (ALWAYS DONE)
  if ((Options.ave_hydrograph) && (tt.model_time!=0))
  {
    _HYDRO<<" "<<GetAveragePrecip();
    for (int p=0;p<_nSubBasins;p++){
      if (_pSubBasins[p]->IsGauged())
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
      if (_pSubBasins[p]->IsGauged())
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
  
  int    retval;      	      	// error value for NetCDF routines
  int    time_id;     	      	// variable id in NetCDF for time
  int    precip_id;   	      	// variable id in NetCDF for precipitation
  int    qsim_id;  	      	// variable id in NetCDF for simulated outflow
  int    qobs_id;  	      	// variable id in NetCDF for observed outflow
  int    qin_id;  	      	// variable id in NetCDF for observed inflow
  size_t start1[1], count1[1];  // determines where and how much will be written to NetCDF; 1D variable (pre, time)
  size_t start2[2], count2[2];  // determines where and how much will be written to NetCDF; 2D variable (qsim, qobs, qin)
  double current_time[1];       // current time in days since start time
  double current_prec[1];       // precipitation of current time step

  int Isim, Nsim; // current and total # of sub-basins with simulated outflows
  int Iobs, Nobs; // current and total # of sub-basins with observed outflows
  int Ires, Nres; // current and total # of sub-basins with observed inflows

  /* (a) count how many values need to be written for q_obs, q_sim, q_in */
  Nsim = 0;
  Nobs = 0;
  Nres = 0;
  if ((Options.ave_hydrograph) && (current_time[0] != 0.0)){
    for (int p=0;p<_nSubBasins;p++){
      if (_pSubBasins[p]->IsGauged()){
	Nsim = Nsim+1;
	for (int i = 0; i < _nObservedTS; i++){
	  if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "HYDROGRAPH")) &&
	      (s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
	      (_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular)){
	    Nobs = Nobs+1;
	  }
	}
	if (_pSubBasins[p]->GetReservoir() != NULL){ Nres = Nres+1;}
      }
    }      
  }
  else {  // point-value hydrograph
    for (int p=0;p<_nSubBasins;p++){
      if (_pSubBasins[p]->IsGauged()){
	Nsim = Nsim+1;
	for (int i = 0; i < _nObservedTS; i++){
	  if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "HYDROGRAPH")) &&
	      (s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
	      (_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular)){
	    Nobs = Nobs+1;
	  }
	}
	if (_pSubBasins[p]->GetReservoir() != NULL){ Nres = Nres+1;}
      }
    }
  }

  /* (b) allocate memory if necessary */
  double outflow_obs[Nobs];   // q_obs
  double outflow_sim[Nsim];   // q_sim
  double inflow_obs[Nres];    // q_in

  /* (c) read data */
  Isim = 0;
  Iobs = 0;
  Ires = 0;
  current_time[0] = tt.model_time;
  current_prec[0] = -9999.0;
  if ((Options.ave_hydrograph) && (current_time[0] != 0.0)){
    current_prec[0] = GetAveragePrecip();
    
    for (int p=0;p<_nSubBasins;p++){
      if (_pSubBasins[p]->IsGauged()){
	outflow_sim[Isim] = _pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
	Isim = Isim+1;
	
	//if (Options.print_obs_hydro)
	{
	  for (int i = 0; i < _nObservedTS; i++){
	    if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "HYDROGRAPH")) &&
		(s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
		(_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular)){
	      double val = _pObservedTS[i]->GetAvgValue(current_time[0],Options.timestep); //time shift handled in CTimeSeries::Parse
	      if ((val != CTimeSeriesABC::BLANK_DATA) && (current_time[0]>0)){ outflow_obs[Iobs] = val;    }
	      else                                                           { outflow_obs[Iobs] = -9999.; }
	      Iobs = Iobs+1;
	    }
	  }
	}
	if (_pSubBasins[p]->GetReservoir() != NULL){
	  inflow_obs[Ires] = _pSubBasins[p]->GetIntegratedReservoirInflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
	  Ires = Ires+1;
	}
      }
    }      
  }
  else {  // point-value hydrograph
    if (current_time[0] != 0.0){current_prec[0] = GetAveragePrecip();} //watershed-wide precip
    else                       {current_prec[0] = -9999.;}   // was originally '---'
    for (int p=0;p<_nSubBasins;p++){
      if (_pSubBasins[p]->IsGauged()){
	outflow_sim[Isim] = _pSubBasins[p]->GetOutflowRate();
	Isim = Isim+1;
	
	//if (Options.print_obs_hydro)
	{
	  for (int i = 0; i < _nObservedTS; i++){
	    if ((!strcmp(_pObservedTS[i]->GetName().c_str(), "HYDROGRAPH")) &&
		(s_to_l(_pObservedTS[i]->GetTag().c_str()) == _pSubBasins[p]->GetID()) &&
		(_pObservedTS[i]->GetType() == CTimeSeriesABC::ts_regular)){
	      double val = _pObservedTS[i]->GetAvgValue(current_time[0],Options.timestep);
	      if ((val != CTimeSeriesABC::BLANK_DATA) && (current_time[0]>0)){ outflow_obs[Iobs] = val;    }
	      else                                                           { outflow_obs[Iobs] = -9999.; }
	      Iobs = Iobs+1;
	    }
	  }
	}
	if (_pSubBasins[p]->GetReservoir() != NULL){
	  inflow_obs[Ires] = _pSubBasins[p]->GetReservoirInflow();
	  Ires = Ires+1;
	}
      }
    }
  }

  /* (d) write values to NetCDF (one time step) */
  start1[0] = round(current_time[0]/Options.timestep);   // element of NetCDF array that will be written
  count1[0] = 1;                                         // writes exactly one time step
  start2[0] = round(current_time[0]/Options.timestep);   // element of NetCDF array that will be written
  start2[1] = 0;                                         // element of NetCDF array that will be written
  count2[0] = 1;                                         // writes exactly one time step

  /* (d1) inquire variable IDs */
  if ((retval = nc_inq_varid (_HYDRO_ncid, "time",   &time_id))){
    nc_error_exit(retval);
  }
  if ((retval = nc_inq_varid (_HYDRO_ncid, "precip", &precip_id))){
    nc_error_exit(retval);
  }
  if (Nsim > 0){
    if ((retval = nc_inq_varid (_HYDRO_ncid, "q_sim", &qsim_id))){
      nc_error_exit(retval);}
  }
  if (Nobs > 0){
    if ((retval = nc_inq_varid (_HYDRO_ncid, "q_obs", &qobs_id))){
      nc_error_exit(retval);}
  }
  if (Nres > 0){
    if ((retval = nc_inq_varid (_HYDRO_ncid, "q_in",  &qin_id))){
      nc_error_exit(retval);}
  }

  /* (d2) write new time step */
  if ((retval = nc_put_vara_double(_HYDRO_ncid, time_id, start1, count1, &current_time[0]))){
    nc_error_exit(retval);
  }
  
  /* (d3) write precipitation values */
  if ((retval = nc_put_vara_double(_HYDRO_ncid, precip_id, start1, count1, &current_prec[0]))){
    nc_error_exit(retval);
  }

  /* (d4) write simulated outflow values */
  if (Nsim > 0){
    count2[1] = Nsim;   // writes exactly Nsim elements
    if ((retval = nc_put_vara_double(_HYDRO_ncid, qsim_id, start2, count2, &outflow_sim[0]))){
      nc_error_exit(retval);
    }
  }
  /* (d5) write observed  outflow values */
  if (Nobs > 0){
    count2[1] = Nobs;   // writes exactly Nobs elements
    if ((retval = nc_put_vara_double(_HYDRO_ncid, qobs_id, start2, count2, &outflow_obs[0]))){
      nc_error_exit(retval);
    }
  }
  /* (d6) write observed  inflow  values */
  if (Nres > 0){
    count2[1] = Nres;   // writes exactly Nres elements
    if ((retval = nc_put_vara_double(_HYDRO_ncid, qin_id, start2, count2, &inflow_obs[0]))){
      nc_error_exit(retval);
    }
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

///////////////////////////////////////////////////////////////////
/// \brief NetCDF error handling
/// \return Error string and exit code
//
void nc_error_exit(int error_code){

#ifdef _RVNETCDF_
  printf("Error: %s\n", nc_strerror(error_code));
  ExitGracefullyIf(error_code,
                   "NetCDF: NetCDF error occured",BAD_DATA);
#endif
}
