/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2019 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"

/*****************************************************************
   Model Streamflow Assimilation Routines
*****************************************************************/
bool IsContinuousFlowObs2(CTimeSeriesABC *pObs,long SBID)
{
 // clears up  terribly ugly repeated if statements
  if(pObs==NULL)                                   { return false; }
  if(s_to_l(pObs->GetTag().c_str()) != SBID)       { return false; }
  if(pObs->GetType() != CTimeSeriesABC::TS_REGULAR){ return false; }
  return (!strcmp(pObs->GetName().c_str(),"HYDROGRAPH")); //name ="HYDROGRAPH"      
}
void CModel::InitializeDataAssimilation(const optStruct &Options)
{
  if(!Options.assimilation_on) { return; }
  _aDAscale    =new double[_nSubBasins];
  _aDAlength   =new double[_nSubBasins];
  _aDAtimesince=new double[_nSubBasins];
  for(int p=0; p<_nSubBasins; p++) {
    _aDAscale[p]=1.0;
    _aDAlength[p]=0.0;
    _aDAtimesince[p]=0.0;
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
  if (!Options.assimilation_on){return;}
  
  //assimilates all data after assimilation date 
  if(tt.model_time<(Options.assimilation_start+Options.timestep/2.0)){return;}
  
  int    p,pdown;
  double t_observationsOFF=ALMOST_INF;//Only for debugging
  double Qobs,Qmod;
  double alpha     =CGlobalParams::GetParams()->assimilation_fact; //for now: 0->no assimilation 1->full override
  double time_fact =CGlobalParams::GetParams()->assim_time_decay; 

  if(alpha>0.0) { alpha = exp(4.5*(alpha-1.0)); } //makes range natural, i.e., alph~0.5 means result is halfway between unassimilated modeled and observed
    
  for(int pp=_nSubBasins-1; pp>=0; pp--)//downstream to upstream
  {
    p=GetOrderedSubBasinIndex(pp);
    bool ObsExists=false;
    for(int i=0; i<_nObservedTS; i++) //determine whether flow observation is available
    {
      if(IsContinuousFlowObs2(_pObservedTS[i],_pSubBasins[p]->GetID()))//flow observation is available
      {
        Qobs = _pObservedTS[i]->GetAvgValue(tt.model_time+Options.timestep,Options.timestep);//correction for period ending storage of hydrograph
        Qmod = _pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
        if((Qobs!=RAV_BLANK_DATA) && (tt.model_time<t_observationsOFF)) 
        {
          if(Qmod>=0.0) { _aDAscale[p]=1.0+alpha*((Qobs-Qmod)/Qmod);}
          else          { _aDAscale[p]=1.0;                         }
          _aDAlength   [p]=0.0;
          _aDAtimesince[p]=0.0;
        }
        else 
        { //found a blank or zero flow value -scaling extinguishes over time
          _aDAtimesince[p]+=Options.timestep;
          _aDAscale    [p]=1.0+( _aDAscale[p]-1.0)*exp(-time_fact*_aDAtimesince[p]); //move the scale towards 1.0
          _aDAlength   [p]=0.0;
        }
        ObsExists=true;
      }
    }

    pdown=GetDownstreamBasin(p);
    if(ObsExists==false){ //observations may be downstream, propagate scaling upstream
      if(pdown!=DOESNT_EXIST){ 
        _aDAscale    [p]=_aDAscale [pdown];
        _aDAlength   [p]=_aDAlength[pdown]+_pSubBasins[pdown]->GetReachLength();
        _aDAtimesince[p]=_aDAtimesince[pdown];
      }
      else{ //Nothing downstream, no assimilation
        _aDAscale    [p]=1.0; 
        _aDAlength   [p]=0.0;
        _aDAtimesince[p]=0.0;
      }
    }
  }// end downstream to upstream

  //Calculate artificial mass added/removed by data assimilation process
  //--------------------------------------------------------------------------------
  double mass_added=0;
  double scalefact=1.0;
  double distfact=CGlobalParams::GetParams()->assim_upstream_decay/M_PER_KM; //[1/km]->[1/m] //TMP DEBUG
  for(p=0; p<_nSubBasins; p++)
  {
    scalefact =1.0+(_aDAscale[p]-1.0)*exp(-distfact*_aDAlength[p]);
    //if(tt.model_time>(t_observationsOFF-1.0)) { cout<<tt.model_time<<" p:"<<p<<" "<<scalefact<<" "<<_aDAlength[p]<<" "<<_aDAtimesince[p]<<endl; }
    mass_added+=_pSubBasins[p]->ScaleAllFlows(scalefact,Options.timestep); 
  }
  if(mass_added>0.0){_CumulInput +=mass_added/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;}
  else              {_CumulOutput-=mass_added/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;}
}