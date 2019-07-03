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
/////////////////////////////////////////////////////////////////
/// \brief sifts through all observation time series. Overrides flows with observations at gauges and scales all flows and channel storages upstream of those gauges. 
/// \note updates cumulative mass input/output to reflect mass added during scaling
/// \param tt [in] current model time structure
/// \param Options [in] current model options structure
//
void CModel::AssimilateStreamflow(const time_struct &tt, const optStruct &Options)
{
  //assimilates all data prior to assimilation date 
  if(tt.model_time<(Options.assimilation_start+Options.timestep/2.0))
  {
    double Qobs,Qmod;
    int p,pdown;
    double alpha =1.0; //_AssimilationFactor //for now: 0->no assimilation 1->full override
    double *scale =new double[_nSubBasins];
    double *length=new double[_nSubBasins];

    for(int p=0; p<_nSubBasins; p++){ 
      scale[p]=1.0; 
      length[0]=0.0;
    }
    for(int pp=_nSubBasins-1; pp>=0; pp--)//downstream to upstream
    {
      p=GetOrderedSubBasinIndex(pp);
      bool found=false;
      for(int i=0; i<_nObservedTS; i++) //determine whether flow observation is available
      {
        if(IsContinuousFlowObs2(_pObservedTS[i],_pSubBasins[p]->GetID()))//flow observation is available
        {
          Qobs = _pObservedTS[i]->GetAvgValue(tt.model_time-Options.timestep,Options.timestep);
          Qmod = _pSubBasins[p]->GetIntegratedOutflow(Options.timestep);
          //Qmod = _pSubBasins[p]->GetOutflowRate(); //UNCLEAR WHICH TO USE!
          if((Qmod>=0.0) && (Qobs!=RAV_BLANK_DATA)){
            scale[p]=1.0+alpha*(Qobs-Qmod)/Qmod;
            length[0]=0.0;
            found=true;
          }
        }
      }
      pdown=GetDownstreamBasin(p);

      if(found==false){
        if(pdown!=DOESNT_EXIST){scale[p]=scale[pdown];length[p]=length[pdown]+_pSubBasins[pdown]->GetReachLength();}
        else                   {scale[p]=1.0; length[p]=0.0;}
      }
    }// end downstream to upstream

    double mass_added=0;
    double scalefact=1.0;
    double distfact=0.0; //[1/m] //_AssimilationDecayFactor
    for(int p=0; p<_nSubBasins; p++)
    {
      scalefact =1.0+(scale[p]-1.0)*exp(-distfact*length[p]);
      mass_added+=_pSubBasins[p]->ScaleAllFlows(scalefact,Options.timestep); 
    }
    if(mass_added>0.0){_CumulInput +=mass_added/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;}
    else              {_CumulOutput-=mass_added/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;}

    delete[]scale;
  }
}