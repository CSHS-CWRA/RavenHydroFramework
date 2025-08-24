/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2025 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"

/*****************************************************************
   Model Streamflow Assimilation Routines
*****************************************************************/
bool IsContinuousFlowObs2(const CTimeSeriesABC *pObs,long long SBID)
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
    _aDAscale     =new double[_nSubBasins];
    _aDAscale_last=new double[_nSubBasins];
    _aDAQadjust   =new double[_nSubBasins];
    _aDADrainSum  =new double[_nSubBasins];
    _aDADownSum   =new double[_nSubBasins];
    _aDAlength    =new double[_nSubBasins];
    _aDAtimesince =new double[_nSubBasins];
    _aDAoverride  =new bool  [_nSubBasins];
    _aDAobsQ      =new double[_nSubBasins];
    _aDAobsQ2     =new double[_nSubBasins];
    _aDASinceLastBlank=new double [_nSubBasins];
    for(int p=0; p<_nSubBasins; p++) {
      _aDAscale     [p]=1.0;
      _aDAscale_last[p]=1.0;
      _aDAQadjust   [p]=0.0;
      _aDADrainSum  [p]=0.0;
      _aDADownSum   [p]=0.0;
      _aDAlength    [p]=0.0;
      _aDAtimesince [p]=0.0;
      _aDAoverride  [p]=false;
      _aDAobsQ      [p]=0.0;
      _aDAobsQ2     [p]=0.0;
      _aDASinceLastBlank[p]=0.0; //initialized to zero to ramp from IC (default)
    }
    int count=0;
    for(int p=0; p<_nSubBasins; p++) {
      if(_pSubBasins[p]->UseInFlowAssimilation())
      {
        count++;
      }
    }
    if(count==0) {
      WriteWarning("InitializeDataAssimlation: :AssimilateStreamflow command was used in .rvi file, but no observations were enabled for assimilation using :AssimilateStreamflow command in .rvt file",Options.noisy);
    }
  }

}
/////////////////////////////////////////////////////////////////
/// \brief Overrides flows with observations at gauges and propagates adjustments upstream of gauges
/// \param p [in] global subbasin index
/// \param Options [in] current model options structure
/// \param tt [in] current model time structure
//
void CModel::AssimilationOverride(const int p,const optStruct& Options,const time_struct& tt)
{
  if (!Options.assimilate_flow)    { return; }
  if (!_pSubBasins[p]->IsEnabled()){ return; }

  //Get current, up-to-date flow and calculate scaling factor
  //---------------------------------------------------------------------
  double Qobs;
  if(_aDAoverride[p])
  {
    double Qobs2,Qmod,Qmodlast;
    double alpha = _pGlobalParams->GetParams()->assimilation_fact;

    //alpha*=(1.0-exp(-0.06*_aDASinceLastBlank[p])); //shock/oscillation prevention

    Qobs    = _aDAobsQ[p];
    Qobs2   = _aDAobsQ2[p];
    Qmod    = _pSubBasins[p]->GetOutflowRate();
    Qmodlast= _pSubBasins[p]->GetLastOutflowRate();
    if(Qmod>PRETTY_SMALL) {
      _aDAscale  [p]=1.0+alpha*((Qobs-Qmod)/Qmod); //if alpha = 1, Q=Qobs in observation basin
      //_aDAQadjust[p]=alpha*(Qobs-Qmod); //Option A: end of time step flow 
      //_aDAQadjust[p]=0.5*alpha*(2.0*Qobs-Qmodlast-Qmod);//Option B: mean flow - rapidly oscillatory Q - Qmean is perfect
      if (Qobs2 == RAV_BLANK_DATA) { //Option C: aim for midpoint of two observation flows (this is the way)
        _aDAQadjust[p]=alpha*(Qobs-Qmod);
      }
      else {
        _aDAQadjust[p]=alpha*(0.5*(Qobs2+Qobs)-Qmod);
      }
    }
    else {
      _aDAscale  [p]=1.0;
      _aDAQadjust[p]=0.0;
    }
    _aDAscale_last[p]=1.0;//no need to scale previous
  }

  // actually scale flows
  //---------------------------------------------------------------------
  double mass_added=0.0;
  if (Options.assim_method==DA_RAVEN_DEFAULT){
    mass_added=_pSubBasins[p]->ScaleAllFlows(_aDAscale[p]/_aDAscale_last[p],_aDAoverride[p],Options.timestep,tt.model_time);
  }
  else if (Options.assim_method==DA_ECCC) {
    if(_aDAoverride[p]){
      mass_added=_pSubBasins[p]->AdjustAllFlows(_aDAQadjust[p],true,true,Options.timestep,tt.model_time);
    }
  }

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
  if (!Options.assimilate_flow)                                       {return;}
  if (tt.model_time<(Options.assimilation_start-Options.timestep/2.0)){return;}//assimilates all data after or on assimilation date

  int    p,pdown;
  double Qobs,Qobs2;
  double t_observationsOFF=ALMOST_INF;//Only used for debugging - keep as ALMOST_INF otherwise
  bool   ObsExists; //observation available in THIS basin

  for(p=0; p<_nSubBasins; p++) {
    _aDAscale_last[p]=_aDAscale[p];
    _aDADrainSum[p]=0.0;
  }

  // check for observations in each basin
  //----------------------------------------------------------------
  int nn=(int)((tt.model_time+TIME_CORRECTION)/Options.timestep)+1;//end of timestep index

  for(int pp=_nSubBasins-1; pp>=0; pp--)//downstream to upstream
  {
    p=GetOrderedSubBasinIndex(pp);

    pdown=GetDownstreamBasin(p);

    ObsExists=false; 

    if(_pSubBasins[p]->UseInFlowAssimilation())
    {
      for(int i=0; i<_nObservedTS; i++) //determine whether flow observation is available
      {
        if(IsContinuousFlowObs2(_pObservedTS[i],_pSubBasins[p]->GetID()))//flow observation is available and linked to this subbasin
        {
          Qobs = _pObservedTS[i]->GetSampledValue(nn); //mean timestep flow
          Qobs2 = _pObservedTS[i]->GetSampledValue(nn+1); //mean timestep flow
          ObsExists=true;
          break; //avoids duplicate observations
        }
      }
    }

    pdown=GetDownstreamBasin(p);

    // observations in this basin, determine scaling variables based upon blank/not blank
    //----------------------------------------------------------------
    if (ObsExists) {
      if((Qobs!=RAV_BLANK_DATA) && (tt.model_time<t_observationsOFF))
      {
        //_aDAscale[p] calculated live in AssimilationOverride when up-to-date modelled flow available
        _aDAlength   [p]=0.0;
        _aDAtimesince[p]=0.0;
        _aDAoverride [p]=true;
        _aDAobsQ     [p]=Qobs;
        _aDAobsQ2    [p]=Qobs2;
        _aDASinceLastBlank[p]+=1.0;
      }
      else
      { //found a blank or zero flow value
        _aDAscale    [p]=_aDAscale[p];//same adjustment as before - scaling persists
        _aDAlength   [p]=0.0;
        _aDAtimesince[p]+=Options.timestep;
        _aDAoverride [p]=false;
        _aDAobsQ     [p]=0.0;
        _aDAobsQ2    [p]=0.0;
        _aDASinceLastBlank[p]+=0.0;
      }
    }
    // no observations in this basin, get scaling from downstream
    //----------------------------------------------------------------
    else if(!ObsExists)  //observations may be downstream, propagate scaling upstream
    {
      //if ((pdown!=DOESNT_EXIST) && (!_aDAoverride[pdown])){ //alternate - allow information to pass through reservoirs
      if( (pdown!=DOESNT_EXIST) && 
          (_pSubBasins[p]->GetReservoir()==NULL) && 
          (!_aDAoverride[p]) && 
          (_pSubBasins[p]->IsEnabled()) && (_pSubBasins[pdown]->IsEnabled())) {
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

  //Calculate _aDADrainSum, sum of assimilated drainage areas upstream of a subbasin outlet
  // and _aDADownSum, drainage of nearest assimilated flow observation
  // dynamic because data can disappear mid simulation
  //-------------------------------------------------------------------
  for(int pp=0;pp<_nSubBasins; pp++)
  {
    _aDADrainSum[p]=0.0;
    _aDADownSum [p]=0.0;
  }
  for(int pp=0;pp<_nSubBasins; pp++)//upstream to downstream
  {
    p=GetOrderedSubBasinIndex(pp);
    pdown=GetDownstreamBasin(p);

    if (_aDAoverride[p]) {
      _aDADrainSum[p]=_pSubBasins[p]->GetDrainageArea();
    }
    /*
    else if (_pSubBasins[p]->GetReservoir()!=NULL){ //??
      _aDADrainSum[p]=0.0;
    }
    */
    else if (pdown!=DOESNT_EXIST){
      _aDADrainSum[pdown]+=_aDADrainSum[p];
    }
  }
  for(int pp=_nSubBasins-1;pp>=0; pp--)// downstream to upstream
  {
    p=GetOrderedSubBasinIndex(pp);
    pdown=GetDownstreamBasin(p);

    if (_aDAoverride[p]) {
      _aDADownSum[p]=_pSubBasins[p]->GetDrainageArea();
    }
    else if (pdown!=DOESNT_EXIST){
      _aDADownSum[p]=_aDADownSum[pdown];
    }
  }

  // Apply time and space correction factors
  //----------------------------------------------------------------
  double time_fact = _pGlobalParams->GetParams()->assim_time_decay;
  double distfact  = _pGlobalParams->GetParams()->assim_upstream_decay/M_PER_KM; //[1/km]->[1/m]
  for(p=0; p<_nSubBasins; p++)
  {
    _aDAscale[p] =1.0+(_aDAscale[p]-1.0)*exp(-distfact*_aDAlength[p])*exp(-time_fact*_aDAtimesince[p]);
  }
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
void CModel::AssimilationBackPropagate(const optStruct &Options,const time_struct &tt)
{
  int p,pdown;

  if (!Options.assimilate_flow)  {return;}

  // Determine flow adjustment magnitude
  //----------------------------------------------------------------
  for(int pp=_nSubBasins-1; pp>=0; pp--)//downstream to upstream
  {
    p    =GetOrderedSubBasinIndex(pp);
    pdown=GetDownstreamBasin(p);

    if(!_aDAoverride[p]) { //observations may be downstream, get scaling from downstream
      if( (pdown!=DOESNT_EXIST                 ) && 
          (!_aDAoverride[p]                    ) && 
          (_pSubBasins[p]->GetReservoir()==NULL) && 
          (_pSubBasins[p]->IsEnabled()) && (_pSubBasins[pdown]->IsEnabled())) 
      {
        _aDAQadjust[p]= _aDAQadjust  [pdown] * (_pSubBasins[p]->GetDrainageArea() / _pSubBasins[pdown]->GetDrainageArea());
      }
      else{ //Nothing downstream or reservoir present in this basin, no adjustment
        _aDAQadjust[p]= 0.0;
      }
    }
  }

  // Apply time and space correction factors
  //----------------------------------------------------------------
  double distfact  = _pGlobalParams->GetParams()->assim_upstream_decay/M_PER_KM; //[1/km]->[1/m]
  double ECCCwt;
  for(p=0; p<_nSubBasins; p++)
  {
    if(!_aDAoverride[p]) 
    {
      pdown=GetDownstreamBasin(p);
      if (pdown!=DOESNT_EXIST){
        _aDAQadjust[p] *=exp(-distfact*_aDAlength[pdown]);
      }

      ECCCwt=1.0;
      if (_aDADrainSum[p]!=0.0){
        ECCCwt = (_pSubBasins[p]->GetDrainageArea() - _aDADrainSum[p])/(_aDADownSum[p]-_aDADrainSum[p]);
      }

      if (!_aDAoverride[p]) { //no scaling for stations being overridden, only upstream
        _aDAQadjust[p] = _aDAQadjust[p]*ECCCwt;
      }
    }
  }

  // Actually update flows
  //----------------------------------------------------------------
  for(p=0; p<_nSubBasins; p++)
  {
    _pSubBasins[p]->AdjustAllFlows(_aDAQadjust[p],false,_aDAoverride[p],Options.timestep,tt.model_time);
  }
}
