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
                         int    to_sv_ind,
                         int    from_HRU_grp,
                         int    to_HRU_grp,
                         bool   constrain_to_SBs,
                         bool   divert,
                         CModel *pModel)
  :CLateralExchangeProcessABC(LAT_FLUSH, pModel)
{
  _iFlushFrom=from_sv_ind;
  _iFlushTo  =to_sv_ind;
  _kk_from   =from_HRU_grp;
  _kk_to     =to_HRU_grp;
  _constrain_to_SBs=constrain_to_SBs;
  _divert    =divert;
  _aFrac     =NULL;

  DynamicSpecifyConnections(0); //purely lateral flow, no vertical

  //check for valid SVs, HRU group indices
  bool badHRU;
  badHRU=(to_HRU_grp<0) || (to_HRU_grp>_pModel->GetNumHRUGroups()-1);
  ExitGracefullyIf(badHRU,"CmvLatFlush::unrecognized 'to' HRU group specified in :LateralFlush/:LateralDivert command",BAD_DATA_WARN);

  badHRU=(from_HRU_grp<0) || (from_HRU_grp>_pModel->GetNumHRUGroups()-1);
  ExitGracefullyIf(badHRU,"CmvLatFlush::unrecognized 'from' HRU group specified in :LateralFlush/:LateralDivert command",BAD_DATA_WARN);

  ExitGracefullyIf(from_sv_ind==DOESNT_EXIST,"CmvLatFlush::unrecognized 'from' state variable specified in :LateralFlush/:LateralDivert command",BAD_DATA_WARN);
  ExitGracefullyIf(to_sv_ind  ==DOESNT_EXIST,"CmvLatFlush::unrecognized 'to' state variable specified in :LateralFlush/:LateralDivert command",BAD_DATA_WARN);
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvLatFlush::~CmvLatFlush(){
  delete [] _aFrac;
}

//////////////////////////////////////////////////////////////////
/// \brief Initialization (prior to solution)
//
void CmvLatFlush::Initialize()
{
  int nConn;
  int *kFrom(NULL), *kTo(NULL);
  kFrom=new int[MAX_LAT_CONNECTIONS];
  kTo = new int[MAX_LAT_CONNECTIONS];
  ExitGracefullyIf(kTo == NULL, "LatFLush::Initialize (1)",OUT_OF_MEMORY);

  CHRUGroup *fromHRUGrp=_pModel->GetHRUGroup(_kk_from);
  CHRUGroup *toHRUGrp  =_pModel->GetHRUGroup(_kk_from);
  _aFrac=NULL;
  _aFrac=new double [MAX_LAT_CONNECTIONS];
  ExitGracefullyIf(_aFrac == NULL, "LatFLush::Initialize",OUT_OF_MEMORY);
  int q=0;
  int k;

  if(_constrain_to_SBs)
  {
    //sift through all HRUs
    /*for (int p = 0; p<_pModel->GetNumSubBasins(); p++)
    {
      //find 'to' HRU (only one allowed per SB)
      int kToSB=DOESNT_EXIST;
      for(int ks=0; ks<_pModel->GetSubBasin(p)->GetNumHRUs(); ks++)
      {
        k=_pModel->GetSubBasin(p)->GetHRU(ks)->GetGlobalIndex();

        if(_pModel->IsInHRUGroup(k,toHRUGrp)) {
          ExitGracefullyIf(kToSB!=DOESNT_EXIST,
            "LatFlush::Initialize - only one HRU per subbasin can recieve flush output. More than one recipient HRU found in subbasin ",BAD_DATA_WARN);
          kToSB=k;
        }
      }

      //find 'from' HRUs to make connections
      for(int ks=0; ks<_pModel->GetSubBasin(p)->GetNumHRUs(); ks++)
      {
        k=_pModel->GetSubBasin(p)->GetHRU(ks)->GetGlobalIndex();

        if(_pModel->IsInHRUGroup(k,fromHRUGrp) && (kToSB!=DOESNT_EXIST)) {
          kFrom[q]=k;
          kTo[q]=kToSB;
          q++;
        }
      }
    }*/
    //Alternate approach: multiple recipient stores
    bool source_sink_issue=false;
    int k1,k2;
    int nRecipients;
    double Asum,area;
    for(int p=0;p<_pModel->GetNumSubBasins();p++)
    {
      if (_pModel->GetSubBasin(p)->IsEnabled()) {
      
        for(int ks=0; ks<_pModel->GetSubBasin(p)->GetNumHRUs(); ks++) //sources
        {
          k1=_pModel->GetSubBasin(p)->GetHRU(ks)->GetGlobalIndex();
          Asum=0.0;
          nRecipients=0;
          if (fromHRUGrp->IsInGroup(k1)) {
          
            for(int ks2=0; ks2<_pModel->GetSubBasin(p)->GetNumHRUs(); ks2++) //recipients
            {
              k2=_pModel->GetSubBasin(p)->GetHRU(ks2)->GetGlobalIndex();

              if (toHRUGrp->IsInGroup(k2)) 
              {
                if (k1!=k2){
                  kFrom[q]=k1;
                  kTo  [q]=k2;
                  area=_pModel->GetSubBasin(p)->GetHRU(ks2)->GetArea();
                  _aFrac[q]=area;
                  Asum+=area;//sum of recipient areas
                  nRecipients++;
                
                  //cout << "ADDING CONNECTION " << q << " in subbasin "<< _pModel->GetSubBasin(p)->GetName() << ": "
                  //     << _pModel->GetHydroUnit(kFrom[q])->GetHRUID()  << " To " <<_pModel->GetHydroUnit(kTo[q])->GetHRUID() <<" "<<_aFrac[q]<<endl;
                  q++;
                  ExitGracefullyIf(q>=MAX_LAT_CONNECTIONS,"Number of HRU connections in LateralFlush or LateralDivert exceeds MAX_LAT_CONNECTIONS limit.",BAD_DATA);
                }
                else {
                  source_sink_issue=true;//throw warning if source is also sink
                }
              }
            }
          }
          for (int qq = 0; qq < nRecipients; qq++) {//once sum is known for each source, calculate fraction
            //cout << "AFRAC INDEX: "<<q-qq-1<<endl;
            _aFrac[q-qq-1]=_aFrac[q-qq-1]/Asum;
          }
        }
      }
    }

    nConn=q;
    if(nConn==0) {
      string warning="CmvLatFlush::Initialize: no connections found between from and to HRU groups in any of the basins. if INTERBASIN flag not used, only transfer within basins is supported. ";
      WriteWarning(warning,true);
    }
    if (source_sink_issue) {
      string warning="CmvLatFlush::Initialize: one or more HRUs belongs to both source and recipient HRU groups in LateralFlush or LateralDivert process. This may have unintentended consequences. ";
      WriteWarning(warning,true);
    }
  }
  else //!constrain_to_SBs
  {
    //toHRUGrp
    int kToSB=DOESNT_EXIST;
    for(k=0;k<_pModel->GetNumHRUs();k++) {
      if (toHRUGrp->IsInGroup(k)) {
        ExitGracefullyIf(kToSB!=DOESNT_EXIST,
          "LatFlush::Initialize - only one HRU in the model can recieve flush output if INTERBASIN is used. More than one recipient HRU found.",BAD_DATA_WARN);
        kToSB=k;
      }
    }
      //find 'from' HRUs to make connections
    for(k=0;k<_pModel->GetNumHRUs();k++)
    {
      if (fromHRUGrp->IsInGroup(k) && (kToSB!=DOESNT_EXIST) && _pModel->GetHydroUnit(k)->IsEnabled()) {
        kFrom[q]=k;
        kTo  [q]=kToSB;
        _aFrac[q]=1.0; //only a single recipient
        q++;
      }
    }
    nConn=q;
    if(nConn==0) {
      string warning="CmvLatFlush::Initialize: no connections found between from and to HRU groups. HRU Group population of zero?";
      WriteWarning(warning,true);
    }
  }

  DynamicSpecifyLatConnections(nConn);
  for(q=0;q<_nLatConnections;q++){
    _kFrom   [q]    =kFrom[q];
    _kTo     [q]    =kTo[q];
    _iFromLat[q]    =_iFlushFrom;
    _iToLat  [q]    =_iFlushTo;
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
  if (_divert) {
    aP[0]="DIVERT_FRACT"; aPC[0]=CLASS_LANDUSE; nP++;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief returns lateral exchange rates (mm/d) between from and to HRU/SV combinations
/// \param **state_vars [in] 2D array of current state variables [nHRUs][nSVs]
/// \param **pHRUs [in] array of pointers to HRUs
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *exchange_rates [out] Rate of loss from "from" compartment [mm-km2/day]
//
void CmvLatFlush::GetLateralExchange( const double * const     *state_vars, //array of all SVs for all HRUs, [k][i]
                                      const CHydroUnit * const *pHRUs,
                                      const optStruct          &Options,
                                      const time_struct        &tt,
                                            double             *exchange_rates) const
{
  double stor,Afrom,Ato;
  double to_stor,max_to_stor,max_rate;
  double mult=1.0; //default: flush 100%
  int    p;

  for(int q=0; q<_nLatConnections; q++)
  {
    stor   =state_vars[_kFrom[q]][_iFromLat[q]];
    to_stor=state_vars[_kTo  [q]][_iToLat[q]];
    Afrom  =pHRUs[_kFrom[q]]->GetArea();
    Ato    =pHRUs[_kTo  [q]]->GetArea();
    max_to_stor=pHRUs[_kTo  [q]]->GetStateVarMax(_iToLat[q],state_vars[_kTo[q]],Options);
    max_rate   =max(max_to_stor-to_stor,0.0)/Options.timestep*Ato;

    if (_divert) {
      mult=pHRUs[_kFrom[q]]->GetSurfaceProps()->divert_fract;
      p   =pHRUs[_kFrom[q]]->GetSubBasinIndex();
      mult*=_pModel->GetSubBasin(p)->GetDivertFract();
      //cout<<" Lat Flush "<<q<<" kFrom: "<<_kFrom[q]<<" "<<_kTo[q] <<" iFrom: "<<_iFromLat[q]<<" "<< mult << " " << _aFrac[q] << " " << stor << endl;
    }

    exchange_rates[q]=max(mult*stor,0.0)/Options.timestep*Afrom*_aFrac[q]; //[mm-m2/d]
    exchange_rates[q]=min(exchange_rates[q],max_rate); //constrains so that it does not overfill receiving compartment
  }
}
