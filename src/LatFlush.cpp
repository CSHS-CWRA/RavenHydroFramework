/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ------------------------------------------------------------------
  LatFlush (abstract move of all water from one compartment in one HRU to another in another HRU)
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "LateralExchangeABC.h"

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/
CmvLatFlush::CmvLatFlush(int   from_sv_ind,
                        int     to_sv_ind,
                        int   from_HRU_grp,
                        int     to_HRU_grp) : CLateralExchangeProcessABC(LAT_FLUSH)
{
  _iFlushFrom=from_sv_ind;
  _iFlushTo  =to_sv_ind;
  _kk_from   =from_HRU_grp;
  _kk_to     =to_HRU_grp;

  DynamicSpecifyConnections(0); //purely lateral flow, no vertical 

  //check for valid SVs, HRU group indices
  bool badHRU=(to_HRU_grp<0) || (to_HRU_grp>_pModel->GetNumHRUGroups()-1);
  ExitGracefullyIf(badHRU,"CmvLatFlush::unrecognized 'to' HRU group specified in :LateralFlush command",BAD_DATA_WARN);
  badHRU=(from_HRU_grp<0) || (from_HRU_grp>_pModel->GetNumHRUGroups()-1);
  ExitGracefullyIf(badHRU,"CmvLatFlush::unrecognized 'from' HRU group specified in :LateralFlush command",BAD_DATA_WARN);
  ExitGracefullyIf(from_sv_ind==DOESNT_EXIST,"CmvLatFlush::unrecognized 'from' state variable specified in :LateralFlush command",BAD_DATA_WARN);
  ExitGracefullyIf(to_sv_ind  ==DOESNT_EXIST,"CmvLatFlush::unrecognized 'to' state variable specified in :LateralFlush command",BAD_DATA_WARN);
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvLatFlush::~CmvLatFlush(){}

//////////////////////////////////////////////////////////////////
/// \brief Initialization (prior to solution)
//
void CmvLatFlush::Initialize()
{
  int nConn;
  int *kFrom=new int[_pModel->GetNumHRUs()];
  int *kTo  =new int[_pModel->GetNumHRUs()];
  string fromHRUGrp=_pModel->GetHRUGroup(_kk_from)->GetName();
  string toHRUGrp  =_pModel->GetHRUGroup(  _kk_to)->GetName();
  int q=0;
  int k;
  bool fromfound=false;

  //sift through all HRUs 
  for(int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    //find 'to' HRU (only one allowed per SB)
    int kToSB=DOESNT_EXIST;
    for(int ks=0; ks<_pModel->GetSubBasin(p)->GetNumHRUs(); ks++)
    {
      k=_pModel->GetSubBasin(p)->GetHRU(ks)->GetGlobalIndex();

      if(_pModel->IsInHRUGroup(k,toHRUGrp)){
        ExitGracefullyIf(kToSB!=DOESNT_EXIST,
          "LatFlush::Initialize - only one HRU per subbasin can recieve flush output. More than one recipient HRU found in subbasin ",BAD_DATA_WARN);
        kToSB=k;
      }
    }

    //find 'from' HRUs to make connections
    for(int ks=0; ks<_pModel->GetSubBasin(p)->GetNumHRUs(); ks++)
    {
      k=_pModel->GetSubBasin(p)->GetHRU(ks)->GetGlobalIndex();

      if(_pModel->IsInHRUGroup(k,fromHRUGrp) && (kToSB!=DOESNT_EXIST)){
        kFrom[q]=k;
        kTo  [q]=kToSB;
        fromfound=true;
        q++;
      }
    }
  }
  nConn=q;
  
  DynamicSpecifyLatConnections(nConn);
  for(int q=0;q<_nLatConnections;q++){
    _kFrom   [q]    =kFrom[q];
    _kTo     [q]    =kTo[q];
    _iFromLat[q]    =_iFlushFrom;
    _iToLat  [q]    =_iFlushTo;
    /*cout <<"latflush connections: "<<_pModel->GetHydroUnit(_kFrom[q])->GetID()<<" "<<_pModel->GetHydroUnit(_kTo[q])->GetID()<<" ";
    cout<<CStateVariable::SVTypeToString(_pModel->GetStateVarType(_iFromLat[q]),-1)<<" ";
    cout<<CStateVariable::SVTypeToString(_pModel->GetStateVarType(_iToLat[q]),-1)<<endl;*/
  }
  delete [] kFrom;
  delete [] kTo;
}


//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param se_type [in] Model of soil evaporation used
/// \param *aSV [out] Array of state variable types needed by soil evaporation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by soil evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvLatFlush::GetParticipatingStateVarList(sv_type *aSV,int *aLev,int &nSV)
{
  nSV=0;
  //user specified 'from' & 'to' compartment, Levels - not known before construction
}

void  CmvLatFlush::GetParticipatingParamList(string *aP,class_type *aPC,int &nP) const
{
  nP=0;
}
//////////////////////////////////////////////////////////////////
/// \brief returns lateral exchange rates (mm/d) between from and to HRU/SV combinations
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *exchange_rates [out] Rate of loss from "from" compartment [mm/day]
//
void CmvLatFlush::GetLateralExchange( const double * const     *state_vars, //array of all SVs for all HRUs, [k][i]
                                      const CHydroUnit * const *pHRUs,    
                                      const optStruct          &Options,
                                      const time_struct        &tt,
                                            double             *exchange_rates) const
{
  double stor;
  for(int q=0; q<_nLatConnections; q++){
    stor=state_vars[_kFrom[q]][_iFromLat[q]];
    exchange_rates[q]+=max(stor,0.0)/Options.timestep;
  }
}