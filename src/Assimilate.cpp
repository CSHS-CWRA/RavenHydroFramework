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
    for(int p=0; p<_nSubBasins; p++) {
      _aDAscale    [p]=1.0;
      _aDAlength   [p]=0.0;
      _aDAtimesince[p]=0.0;
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
      }
    }
    WriteAdvisory("Lake stage data assimilation will lead to mass balance error estimates in the WatershedStorage output file. This is a natural side effect of assimilation.",Options.noisy);
  }

}
/////////////////////////////////////////////////////////////////
/// \brief sifts through all observation time series. Overrides flows with observations at gauges and scales all flows and channel storages upstream of those gauges. 
/// \note updates cumulative mass input/output to reflect mass added during scaling
/// \param tt [in] current model time structure
/// \param Options [in] current model options structure
//
void CModel::AssimilateStreamflow(const optStruct &Options,const time_struct &tt)
{
  if (!Options.assimilate_flow){return;}
  
  //assimilates all data after assimilation date 
  if(tt.model_time<(Options.assimilation_start+Options.timestep/2.0)){return;}
  
  int    p,pdown;
  double t_observationsOFF=ALMOST_INF;//Only used for debugging - keep as ALMOST_INF otherwise
  double Qobs,Qmod;
  double alpha     =CGlobalParams::GetParams()->assimilation_fact; //for now: 0->no assimilation 1->full override
  double time_fact =CGlobalParams::GetParams()->assim_time_decay; 

  int nn           =(int)((tt.model_time+TIME_CORRECTION)/Options.timestep);//current timestep index

  if(alpha>0.0) { alpha = exp(4.5*(alpha-1.0)); } //makes range natural, i.e., alph~0.5 means result is halfway between unassimilated modeled and observed
    
  for(int pp=_nSubBasins-1; pp>=0; pp--)//downstream to upstream
  {
    p=GetOrderedSubBasinIndex(pp);

    bool ObsExists=false; //observation available in THIS basin
    if(_pSubBasins[p]->UseInFlowAssimilation()) 
    {
      for(int i=0; i<_nObservedTS; i++) //determine whether flow observation is available
      {
        if(IsContinuousFlowObs2(_pObservedTS[i],_pSubBasins[p]->GetID()))//flow observation is available and linked to this subbasin
        //&& 
        {
          Qobs = _pObservedTS[i]->GetSampledValue(nn+1);//correction for period ending storage of hydrograph
          Qmod = _pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
          if((Qobs!=RAV_BLANK_DATA) && (tt.model_time<t_observationsOFF))
          {
            if(Qmod>PRETTY_SMALL) { _aDAscale[p]=1.0+alpha*((Qobs-Qmod)/Qmod); }
            else { _aDAscale[p]=1.0; }
            _aDAlength[p]=0.0;
            _aDAtimesince[p]=0.0;
          }
          else
          { //found a blank or zero flow value -scaling extinguishes over time
            _aDAtimesince[p]+=Options.timestep;
            _aDAscale[p]=1.0+(_aDAscale[p]-1.0)*exp(-time_fact*_aDAtimesince[p]); //move the scale towards 1.0
            _aDAlength[p]=0.0;
          }
          ObsExists=true;
          break; //avoids duplicate observations
        }
      }
    }
    pdown=GetDownstreamBasin(p);
    if(ObsExists==false){ //observations may be downstream, propagate scaling upstream 
      //if (pdown!=DOESNT_EXIST){ //alternate - allow information to pass through reservoirs 
      if((pdown!=DOESNT_EXIST) && (_pSubBasins[p]->GetReservoir()==NULL)) {
        _aDAscale    [p]=_aDAscale [pdown];
        _aDAlength   [p]=_aDAlength[pdown]+_pSubBasins[pdown]->GetReachLength();
        _aDAtimesince[p]=_aDAtimesince[pdown];
      }
      else{ //Nothing downstream or reservoir present in this basin, no assimilation
        _aDAscale    [p]=1.0; 
        _aDAlength   [p]=0.0;
        _aDAtimesince[p]=0.0;
      }
    }
  }// end downstream to upstream

  //Calculate artificial mass added/removed by data assimilation process
  // actual changing of river flows occurs in ScaleAllFlows() call
  //--------------------------------------------------------------------------------
  double mass_added=0;
  double scalefact=1.0;
  double distfact=CGlobalParams::GetParams()->assim_upstream_decay/M_PER_KM; //[1/km]->[1/m] 
  for(p=0; p<_nSubBasins; p++)
  {
    scalefact =1.0+(_aDAscale[p]-1.0)*exp(-distfact*_aDAlength[p]);
    mass_added+=_pSubBasins[p]->ScaleAllFlows(scalefact,Options.timestep); 
  }
  if(mass_added>0.0){_CumulInput +=mass_added/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;}
  else              {_CumulOutput-=mass_added/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;}
}