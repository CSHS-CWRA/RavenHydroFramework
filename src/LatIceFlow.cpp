/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2026 the Raven Development Team
------------------------------------------------------------------
CmvLatIceFlow - simulates lateral flow of glacier ice
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "IceFlow.h"

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

CmvLatIceFlow::CmvLatIceFlow(CModel *pModel)
                   : CLateralExchangeProcessABC(LAT_ICE_FLOW, pModel)
{
  DynamicSpecifyConnections(0); // purely lateral flow, no vertical
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvLatIceFlow::~CmvLatIceFlow()
{
}

bool sortByElevation(const CHydroUnit *pA,const CHydroUnit *pB)
{
  return (pA->GetElevation() > pB->GetElevation());
}
//////////////////////////////////////////////////////////////////
/// \brief Initialization (prior to solution)
//
void CmvLatIceFlow::Initialize()
{
  // Calculate the number of lateral connections
  // All subbasins which include at least one glacier HRU get an array of lateral connections
  bool *found;
  found=new bool[_pModel->GetNumSubBasins()];
  for (int p=0;p<_pModel->GetNumSubBasins();p++){
    found[p]=false;
  }
  for (int p=0;p<_pModel->GetNumSubBasins();p++){
    for (int k=0;k<_pModel->GetSubBasin(p)->GetNumHRUs();k++){
      if (_pModel->GetSubBasin(p)->GetHRU(k)->GetHRUType()==HRU_GLACIER){
        found[p]=true;
      }
    }
    if (found[p]){
      _nLatConnections+=(_pModel->GetSubBasin(p)->GetNumHRUs()-1);
    }
  }

  // Allocate memory for the arrays
  _kFrom    = new int[_nLatConnections];
  _kTo      = new int[_nLatConnections];
  _iFromLat = new int[_nLatConnections];
  _iToLat   = new int[_nLatConnections];

  int iGlacierIce=_pModel->GetStateVarIndex(GLACIER_ICE);

  // Initialize the arrays using the methods from CLatConnect
  int q=0;
  for (int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    if (found[p])
    {
      int nHRUs=_pModel->GetSubBasin(p)->GetNumHRUs();
      const CHydroUnit **sortHRUs=new CHydroUnit const *[nHRUs];

      for (int k=0;k<nHRUs;k++) //sort from high elev to low elev
      {
        sortHRUs[k]=_pModel->GetSubBasin(p)->GetHRU(k);
      }

      sortPointerArray<CHydroUnit>(sortHRUs,nHRUs,sortByElevation);

      for (int k=0;k<nHRUs-1;k++)
      { 
        _kFrom   [q] = sortHRUs[k  ]->GetGlobalIndex();
        _kTo     [q] = sortHRUs[k+1]->GetGlobalIndex();

        cout<<" connection : "<<sortHRUs[k]->GetHRUID()<<" (elev= "<<sortHRUs[k]->GetElevation()<<") to "<<sortHRUs[k+1]->GetHRUID()<<" (elev= "<<sortHRUs[k+1]->GetElevation()<<")"<<endl;
        _iFromLat[q] = iGlacierIce;
        _iToLat  [q] = iGlacierIce;
        q++;
      }
      delete [] sortHRUs;
    }
  }
  delete [] found;

  //store initial glacier bottom for all HRUs
  _aInitGlacierHeight=new double [_pModel->GetNumHRUs()];
  for (int k=0;k<_pModel->GetNumHRUs();k++){
    _aInitGlacierHeight[k]=(DENSITY_WATER/DENSITY_ICE)*_pModel->GetHydroUnit(k)->GetStateVarValue(iGlacierIce);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \param se_type [in] Model of soil evaporation used
/// \param *aSV [out] Array of state variable types needed by soil evaporation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by soil evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvLatIceFlow::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
    nSV = 0;
    // user specified 'from' & 'to' compartment, Levels - not known before construction
}

void CmvLatIceFlow::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
    nP = 0;
}

//////////////////////////////////////////////////////////////////
/// \brief returns lateral exchange rates (mm/d) between from and to HRU/SV combinations
/// \param **state_vars [in] 2D array of current state variables [nHRUs][nSVs]
/// \param **pHRUs [in] array of pointers to HRUs
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *exchange_rates [out] Rate of loss from "from" compartment [mm-km2/day]
//
void CmvLatIceFlow::GetLateralExchange(const double *const *state_vars,
                                            const CHydroUnit *const *pHRUs,
                                            const optStruct &Options,
                                            const time_struct &tt,
                                            double *exchange_rates) const
{
  //because HRUs are sorted from high elev to low, 'from' and 'to' are directed downhill, not necessarily in flow direction 
  double stor_from,stor_to,Afrom,Ato;
  double elev_from,elev_to;   //elevation is actually *initial* elevation, including glacier cover, as consistent with most DEMs
  double gelev_from,gelev_to; //ground elevation, [masl]
  double len_from, len_to;    //HRU lengths in flow direction, [m]
  double Hfrom,Hto;           //Glacier ice thickness [m]

  const double A=2.4e-24;     //Pa^-3/s^-1, from Cuffey and Paterson (2010)
  const int    n=3;           //Glen's flow law exponent
  double distance,slope,H,velocity,surf_grad,width;

  int iGlacierIce=_pModel->GetStateVarIndex(GLACIER_ICE);
  for(int q=0; q<_nLatConnections; q++)
  {
    stor_from =state_vars[_kFrom[q]][iGlacierIce];
    stor_to   =state_vars[_kTo  [q]][iGlacierIce];
    Afrom     =pHRUs[_kFrom[q]]->GetArea();
    Ato       =pHRUs[_kTo  [q]]->GetArea();
    elev_from =pHRUs[_kFrom[q]]->GetElevation();
    elev_to   =pHRUs[_kTo  [q]]->GetElevation();
    len_from  =pHRUs[_kFrom[q]]->GetFlowLength();
    len_to    =pHRUs[_kTo  [q]]->GetFlowLength();
    gelev_from=elev_from-_aInitGlacierHeight[_kFrom[q]]/MM_PER_METER;
    gelev_to  =elev_to  -_aInitGlacierHeight[_kTo  [q]]/MM_PER_METER;
    Hfrom     =stor_from/MM_PER_METER*(DENSITY_WATER/DENSITY_ICE);
    Hto       =stor_to  /MM_PER_METER*(DENSITY_WATER/DENSITY_ICE);

    distance=0.5*(len_from+len_to);

    surf_grad=((Hto+gelev_to)-(Hfrom+gelev_from))/distance;

    //slope=(gelev_from-gelev_to)/distance; //bottom surface slope
    slope=surf_grad; //if interpreted as surface gradient (JRC: I think this is appropriate interpretation)

    if (surf_grad<0){H=Hfrom;} 
    else            {H=Hto;  }
    //H=0.5*(Hfrom+Hto); //alternative option

    width=0.5*(Afrom/len_from+Ato/len_to);
    
    velocity = ((2*A)/(n+2))*pow(DENSITY_ICE*GRAVITY*fabs(slope),n)*pow(H,n+1); //[m/s]
    velocity*=SEC_PER_DAY; //[m/d]

    velocity*=(-surf_grad)/fabs(surf_grad); //to get direction of flow right (opposite gradient)

    exchange_rates[q]=velocity*width*H*MM_PER_METER; //[mm-m2/d]

    if (surf_grad<0){exchange_rates[q]=max(exchange_rates[q],+stor_from*Afrom/Options.timestep);}
    else            {exchange_rates[q]=min(exchange_rates[q],-stor_to  *Ato  /Options.timestep);}
  }
}
