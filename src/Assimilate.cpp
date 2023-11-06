/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"

/*****************************************************************
   Model Streamflow Assimilation Routines
*****************************************************************/
bool IsContinuousFlowObs2(const CTimeSeriesABC *pObs,long SBID)
{
 // clears up  terribly ugly repeated if statements
  if(pObs==NULL)                                   { return false; }
  if(pObs->GetLocID() != SBID)                     { return false; }
  if(pObs->GetType() != CTimeSeriesABC::TS_REGULAR){ return false; }
  return (!strcmp(pObs->GetName().c_str(),"HYDROGRAPH")); //name ="HYDROGRAPH"
}
/////////////////////////////////////////////////////////////////
/// \brief initializes memory for assimilation and sets up connections to observation data
/// \param Options [in] current model options structure
//
void CModel::InitializeDataAssimilation(const optStruct &Options)
{
  // Initialize Streamflow assimilation
  if(Options.assimilate_flow)
  {
    _aDAscale    =new double[_nSubBasins];
    _aDAlength   =new double[_nSubBasins];
    _aDAtimesince=new double[_nSubBasins];
    _aDAoverride =new bool  [_nSubBasins];
    _aDAobsQ     =new double[_nSubBasins];
    _aDAlast     =new double[_nSubBasins];
    for(int p=0; p<_nSubBasins; p++) {
      _aDAscale    [p]=1.0;
      _aDAlength   [p]=0.0;
      _aDAtimesince[p]=0.0;
      _aDAoverride [p]=false;
      _aDAobsQ     [p]=0.0;
      _aDAlast     [p]=1.0;
    }
    int count=0;
    for(int p=0; p<_nSubBasins; p++) {
      if(_pSubBasins[p]->UseInFlowAssimilation())
      {
        count++;
      }
    }
    if(count==0) {
      ExitGracefully("InitializeDataAssimlation: :AssimilateStreamflow command was used in .rvi file, but no observations were enabled for assimilation using :AssimilateStreamflow command in .rvt file",BAD_DATA_WARN);
    }
  }

  //Connect Lake assimilation observations
  if(Options.assimilate_stage)
  {
    for(int i=0;i<_nObservedTS;i++) {
      if((!strcmp(_pObservedTS[i]->GetName().c_str(),"RESERVOIR_STAGE")) && (_pObservedTS[i]->GetType() == CTimeSeriesABC::TS_REGULAR)) {
        CSubBasin *pBasin=GetSubBasinByID(_pObservedTS[i]->GetLocID());
        CReservoir *pRes=pBasin->GetReservoir();
        if(pRes!=NULL) {
          pRes->TurnOnAssimilation(_pObservedTS[i]);
        }
        else {
          ExitGracefully("Reservoir stage observations assigned to basin without lake or reservoir. Data cannot be assimilated",BAD_DATA_WARN);
        }
      }
    }
    WriteAdvisory("Lake stage data assimilation will lead to mass balance error estimates in the WatershedStorage output file. This is a natural side effect of assimilation.",Options.noisy);
  }

}
/////////////////////////////////////////////////////////////////
/// \brief Overrides flows with observations at gauges, or
/// \param p [in] global subbasin index
/// \param Options [in] current model options structure
/// \param tt [in] current model time structure
//
void CModel::AssimilationOverride(const int p,const optStruct& Options,const time_struct& tt)
{
  if(!Options.assimilate_flow)                                        { return; }

  //Get current, up-to-date flow and calculate scaling factor
  //---------------------------------------------------------------------
  if(_aDAoverride[p])
  {
    double Qobs,Qmod;
    double alpha = this->_pGlobalParams->GetParams()->assimilation_fact;

    Qobs = _aDAobsQ[p];
    //Option A: mean flow
    //Qmod = _pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
    //Option B: instantaneous flow
    Qmod = _pSubBasins[p]->GetOutflowRate();
    if(Qmod>PRETTY_SMALL) {
      _aDAscale[p]=1.0+alpha*((Qobs-Qmod)/Qmod);
    }
    else {
      _aDAscale[p]=1.0;
    }
    _aDAlast[p]=1.0;//no need to scale previous
  }

  // actually scale flows
  //---------------------------------------------------------------------
  double mass_added=_pSubBasins[p]->ScaleAllFlows(_aDAscale[p]/_aDAlast[p],_aDAoverride[p],Options.timestep,tt.model_time);
  if(mass_added>0.0){_CumulInput +=mass_added/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;}
  else              {_CumulOutput-=mass_added/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;}

}
/////////////////////////////////////////////////////////////////
/// \brief Calculates updated DA scaling coefficients for all subbasins for forthcoming timestep
/// \note most scale coefficients will be 1.0  except:
/// (1) gauges with valid observation data, which scale to override OR
/// (2) basins at or upstream of missing observation data
/// actual scaling performed in CModel::AssimilationOverride() during routing
/// \param tt [in] current model time structure
/// \param Options [in] current model options structure
//
void CModel::PrepareAssimilation(const optStruct &Options,const time_struct &tt)
{
  if (!Options.assimilate_flow)                                      {return;}
  if(tt.model_time<(Options.assimilation_start+Options.timestep/2.0)){return;}//assimilates all data after assimilation date

  int    p,pdown;
  double Qobs;
  double t_observationsOFF=ALMOST_INF;//Only used for debugging - keep as ALMOST_INF otherwise

  int nn=(int)((tt.model_time+TIME_CORRECTION)/Options.timestep);//current timestep index

  for(p=0; p<_nSubBasins; p++) {
    _aDAlast[p]=_aDAscale[p];
  }

  for(int pp=_nSubBasins-1; pp>=0; pp--)//downstream to upstream
  {
    p=GetOrderedSubBasinIndex(pp);

    bool ObsExists=false; //observation available in THIS basin
    // observations in this basin, determine scaling variables based upon blank/not blank
    //----------------------------------------------------------------
    if(_pSubBasins[p]->UseInFlowAssimilation())
    {
      for(int i=0; i<_nObservedTS; i++) //determine whether flow observation is available
      {
        if(IsContinuousFlowObs2(_pObservedTS[i],_pSubBasins[p]->GetID()))//flow observation is available and linked to this subbasin
        {
          Qobs = _pObservedTS[i]->GetSampledValue(nn+1);//+1: correction for period ending storage of hydrograph

          //bool fakeblank=((tt.model_time>30) && (tt.model_time<40)) || ((tt.model_time>45) && (tt.model_time<47));//TMP DEBUG
          //if (fakeblank){Qobs=RAV_BLANK_DATA;}

          if((Qobs!=RAV_BLANK_DATA) && (tt.model_time<t_observationsOFF))
          {
            //_aDAscale[p] calculated live in AssimilationOverride when up-to-date modelled flow available
            _aDAlength   [p]=0.0;
            _aDAtimesince[p]=0.0;
            _aDAoverride [p]=true;
            _aDAobsQ     [p]=Qobs;
          }
          else
          { //found a blank or zero flow value
            _aDAscale    [p]=_aDAscale[p];
            _aDAtimesince[p]+=Options.timestep;
            _aDAlength   [p]=0.0;
            _aDAoverride [p]=false;
            _aDAobsQ     [p]=0.0;
          }
          ObsExists=true;
          break; //avoids duplicate observations
        }
      }
    }
    // no observations in this basin, get scaling from downstream
    //----------------------------------------------------------------
    pdown=GetDownstreamBasin(p);
    if(ObsExists==false) { //observations may be downstream, propagate scaling upstream
      //if ((pdown!=DOESNT_EXIST) && (!_aDAoverride[pdown])){ //alternate - allow information to pass through reservoirs
      if((pdown!=DOESNT_EXIST) && (_pSubBasins[p]->GetReservoir()==NULL) && (!_aDAoverride[pdown])) {
        _aDAscale      [p]= _aDAscale    [pdown];
        _aDAlength     [p]+=_pSubBasins  [pdown]->GetReachLength();
        _aDAtimesince  [p]= _aDAtimesince[pdown];
        _aDAoverride   [p]=false;
      }
      else{ //Nothing downstream or reservoir present in this basin, no assimilation
        _aDAscale    [p]=1.0;
        _aDAlength   [p]=0.0;
        _aDAtimesince[p]=0.0;
        _aDAoverride [p]=false;
      }
    }
  }// end downstream to upstream

  // Apply time and space correction factors
  //----------------------------------------------------------------
  double time_fact = this->_pGlobalParams->GetParams()->assim_time_decay;
  double distfact  = this->_pGlobalParams->GetParams()->assim_upstream_decay/M_PER_KM; //[1/km]->[1/m]
  for(p=0; p<_nSubBasins; p++)
  {
    _aDAscale[p] =1.0+(_aDAscale[p]-1.0)*exp(-distfact*_aDAlength[p])*exp(-time_fact*_aDAtimesince[p]);
  }
}
