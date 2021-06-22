/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ------------------------------------------------------------------
  LatEquilibrate (redistribute water in set of adjacent storage units to equilibrate storage)
  probably more clever than it needs to be
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "LateralExchangeABC.h"

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/
CmvLatEquilibrate::CmvLatEquilibrate(int    sv_ind,
                                     int    HRU_grp,
                                     double mixing_rate, 
                                     bool   constrain_to_SBs) : CLateralExchangeProcessABC(LAT_EQUIL)
{
  _iSV    = sv_ind;
  _kk     = HRU_grp;
  _constrain_to_SBs=constrain_to_SBs;
  _Asum        = NULL;
  _halfconn    = NULL;
  _mixing_rate = mixing_rate; // [%/day]

  DynamicSpecifyConnections(0); //purely lateral flow, no vertical 

  //check for valid SVs, HRU group indices
  bool badHRU=(HRU_grp<0) || (HRU_grp>_pModel->GetNumHRUGroups()-1);
  ExitGracefullyIf(badHRU,"CmvLatEquilibrate::unrecognized HRU group specified in :LatEquilibrate command",BAD_DATA_WARN);
  
  ExitGracefullyIf(sv_ind==DOESNT_EXIST,"CmvLatEquilibrate::unrecognized state variable specified in :LatEquilibrate command",BAD_DATA_WARN);
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvLatEquilibrate::~CmvLatEquilibrate(){
  delete[] _Asum;
  delete[] _halfconn;
}

//////////////////////////////////////////////////////////////////
/// \brief Initialization (prior to solution)
//
void CmvLatEquilibrate::Initialize()
{
  int nConn;
  int *kFrom=new int[_pModel->GetNumHRUs()];
  int *kTo  =new int[_pModel->GetNumHRUs()];
  _Asum = new double[_pModel->GetNumSubBasins()];
  _halfconn = new int [_pModel->GetNumSubBasins()];
  string HRUGrp=_pModel->GetHRUGroup(_kk)->GetName();
  int q=0;
  int k;
  bool fromfound=false;
  for (int p = 0; p < _pModel->GetNumSubBasins(); p++)
  {
    _Asum[p] = 0;
    _halfconn[p] = 0;
  }
  if(_constrain_to_SBs)
  {  
    //sift through all subbasins 
    for(int p=0;p<_pModel->GetNumSubBasins();p++)
    {
      _Asum[p] = 0;
      _halfconn[p] = 0;
      int kToSB=DOESNT_EXIST;
      for(int ks=0; ks<_pModel->GetSubBasin(p)->GetNumHRUs(); ks++)
      {
        k=_pModel->GetSubBasin(p)->GetHRU(ks)->GetGlobalIndex();
        //find first HRU belonging to HRU group
        if ((_pModel->IsInHRUGroup(k, HRUGrp)) && (kToSB == DOESNT_EXIST)) {
          kToSB = k; //first found is target HRU
          _Asum[p] += _pModel->GetHydroUnit(k)->GetArea();
        }
        //find other HRUs to make connections -creation only made if 2+ HRUs in same basin
        else if (_pModel->IsInHRUGroup(k, HRUGrp) && (kToSB != DOESNT_EXIST)) {
          kFrom[q]=k;
          kTo  [q]=kToSB;
          fromfound=true;
          _Asum[p] += _pModel->GetHydroUnit(k)->GetArea();
          _halfconn[p]++;
          q++;
        }
      }
    }
    nConn=2*q; //one to aggregate, one to disaggregate (for constituent mixing)
    if(nConn==0) {
      string warning="CmvLatEquilibrate::Initialize: no connections found between from and to HRU groups in any of the basins. if INTERBASIN flag not used, only transfer within basins is supported. ";
      WriteWarning(warning,true);
    }
  }
  else //!constrain_to_SBs
  {
    int kToSB=DOESNT_EXIST;
    _Asum[0] = 0;    
    for(k=0;k<_pModel->GetNumHRUs();k++) {
      if((_pModel->IsInHRUGroup(k,HRUGrp)) && (kToSB == DOESNT_EXIST)) {
        kToSB=k; //first in group becomes aggregation target
        _Asum[0]+= _pModel->GetHydroUnit(k)->GetArea();
      }
    }

    //find 'from' HRUs to make connections
    for(k=0;k<_pModel->GetNumHRUs();k++)
    {
      if(_pModel->IsInHRUGroup(k,HRUGrp) && (kToSB!=DOESNT_EXIST)) {
        kFrom[q]=k;
        kTo  [q]=kToSB;
        fromfound=true;
        _Asum[0] += _pModel->GetHydroUnit(k)->GetArea();
        _halfconn[0]++;
        q++;
      }
    }
    nConn=2*q;//one to aggregate, one to disaggregate (for constituent mixing)
  }

  DynamicSpecifyLatConnections(nConn);
  int halfconn = _nLatConnections / 2;
  for (q = 0; q < _nLatConnections; q++) {
    if (q < halfconn) {
      _kFrom[q] = kFrom[q];
      _kTo  [q] = kTo  [q];
    }
    else{
      _kFrom[q] = kTo  [q-halfconn];
      _kTo  [q] = kFrom[q-halfconn];
    }
    //cout << " conn[" << q << " kFrom: " << _kFrom[q] << " kTo: " << _kTo[q] << endl;
    _iFromLat[q]    =_iSV;
    _iToLat  [q]    =_iSV;
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
void CmvLatEquilibrate::GetParticipatingStateVarList(sv_type *aSV,int *aLev,int &nSV)
{
  nSV=0;
  //user specified 'from' & 'to' compartment, Levels - not known before construction
}

void  CmvLatEquilibrate::GetParticipatingParamList(string *aP,class_type *aPC,int &nP) const
{
  nP=0;
}
//////////////////////////////////////////////////////////////////
/// \brief returns lateral exchange rates (mm/d) between from and to HRU/SV combinations
/// \param **state_vars [in] 2D array of current state variables [nHRUs][nSVs]
/// \param **pHRUs [in] array of pointers to HRUs
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *exchange_rates [out] Rate of loss from "from" compartment [mm-km2/day]
//
void CmvLatEquilibrate::GetLateralExchange( const double * const     *state_vars, //array of all SVs for all HRUs, [k][i]
                                      const CHydroUnit * const *pHRUs,    
                                      const optStruct          &Options,
                                      const time_struct        &tt,
                                            double             *exchange_rates) const
{
  double stor,Afrom,Ato;
  double to_stor;
  int p;
  int halfconn = _nLatConnections / 2;

  double mix = min(_mixing_rate*Options.timestep,1.0);  //_mixing rate is in units of % mixed/day * [days/timestep]-> % mixed per timestep

  int qstart = 0;
  double* sum = new double[_pModel->GetNumSubBasins()];
  for (int p = 0; p < _pModel->GetNumSubBasins(); p++) { sum[p] = 0.0; }
  int plast = -1;
  for (int q = 0; q < _nLatConnections; q++)
  {
    stor    = state_vars[_kFrom[q]][_iFromLat[q]];
    to_stor = state_vars[_kTo[q]][_iToLat[q]];
    Afrom   = pHRUs[_kFrom[q]]->GetArea();
    Ato     = pHRUs[_kTo[q]]->GetArea();
    p = _pModel->GetHydroUnit(_kFrom[q])->GetSubBasinIndex();
    if (!_constrain_to_SBs) { p = 0; }

    if (p!=plast) {
      sum[p] += mix * max(to_stor, 0.0) * Ato;
      plast = p;
      qstart = q;
    }
    if ((q >= qstart) && (q < qstart+_halfconn[p])) //step 1: mix % of water into single vessel HRU 
    {
      exchange_rates[q] = mix*max(stor, 0.0) / Options.timestep * Afrom; //[mm-m2/d] //a fraction of storage moved to mixing unit per time step
      sum[p] += exchange_rates[q]*Options.timestep; //[mm-m2] cumulative water mixed
    }
    else //step 2: redistribute mixed water back to contributing HRUs to even out water levels 
    {  
      exchange_rates[q] = sum[p] *(Ato / _Asum[p]) / Options.timestep; //[mm-m2/d]fraction empties out
    }
    
  }
  //JRC: NOTE SOME BIAS BASED UPON CHOICE OF MIXING UNIT if mix>1, WHICH WILL ALWAYS EQUILIBRATE FASTER
  delete[] sum;
}