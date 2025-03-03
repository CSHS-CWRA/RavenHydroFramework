/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ------------------------------------------------------------------
  Soil Balance
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "SoilWaterMovers.h"
#include "Model.h"

/*****************************************************************
   SoilBalance Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the standard constructor
/// \param cr_type [in] Model of capillary rise selected
/// \param In_index [in] Soil storage unit index from which water is lost
/// \param Out_index [in] Soil storage unit index to which water rises
//
CmvSoilBalance::CmvSoilBalance(soilbal_type sb_type,
                               CModelABC    *pModel)
  :CHydroProcessABC(SOIL_BALANCE, pModel)
{
  _type =sb_type;

  if(_type==SOILBAL_SACSMA)
  {
    CHydroProcessABC::DynamicSpecifyConnections(13);

    int iPond,iSW,iUZT,iADIM,iLZT,iUZF,iLZFS,iLZFP,iGW;
    iPond=pModel->GetStateVarIndex(PONDED_WATER);
    iSW  =pModel->GetStateVarIndex(SURFACE_WATER);
    iUZT =pModel->GetStateVarIndex(SOIL,0);
    iUZF =pModel->GetStateVarIndex(SOIL,1);
    iLZT =pModel->GetStateVarIndex(SOIL,2);
    iLZFP=pModel->GetStateVarIndex(SOIL,3);
    iLZFS=pModel->GetStateVarIndex(SOIL,4);
    iADIM=pModel->GetStateVarIndex(SOIL,5);
    iGW  =pModel->GetStateVarIndex(SOIL,6);

    iFrom[ 0]=iPond;       iTo[ 0]=iSW;        //rates[0]:  PONDED_WATER->SURFACE_WATER
    iFrom[ 1]=iPond;       iTo[ 1]=iUZT;       //rates[1]:  PONDED_WATER->UZT
    iFrom[ 2]=iPond;       iTo[ 2]=iADIM;      //rates[2]:  PONDED_WATER->ADIMC
    iFrom[ 3]=iLZFP;       iTo[ 3]=iSW;        //rates[3]:  LZFP->SURFACE_WATER
    iFrom[ 4]=iLZFP;       iTo[ 4]=iGW;        //rates[4]:  LZFP->GW
    iFrom[ 5]=iLZFS;       iTo[ 5]=iSW;        //rates[5]:  LZFS->SURFACE_WATER
    iFrom[ 6]=iLZFS;       iTo[ 6]=iGW;        //rates[6]:  LZFS->GW
    iFrom[ 7]=iUZF;        iTo[ 7]=iLZT;       //rates[7]:  UZF->LZT
    iFrom[ 8]=iUZF;        iTo[ 8]=iLZFS;      //rates[8]:  UZF->LZFS
    iFrom[ 9]=iUZF;        iTo[ 9]=iLZFP;      //rates[9]:  UZF->LZFP
    iFrom[10]=iUZF;        iTo[10]=iSW;        //rates[10]: UZF->SURFACE_WATER
    iFrom[11]=iPond;       iTo[11]=iUZF;       //rates[11]: PONDED->UZF
    iFrom[12]=iADIM;       iTo[12]=iADIM;      //rates[12]: ADIMC->ADIMC [violates mass balance!]

  }
  //-----------------------------------------------------------------
  else
  {
    ExitGracefully("CmvSoilBalance::constructor: undefined soil balance type",BAD_DATA);
  }

}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of default destructor
//
CmvSoilBalance::~CmvSoilBalance(){}

//////////////////////////////////////////////////////////////////
/// \brief Verifies iFrom - iTo connectivity
//
void CmvSoilBalance::Initialize()
{

  //check for porosity!=1.0
  if(_type==SOILBAL_SACSMA)
  {
    ExitGracefullyIf(pModel->GetNumSoilLayers()<7,
      "CmvSoilBalance::Initialize: must have at least 7 soil layers to use SACSMA Soil balance accounting",BAD_DATA_WARN);
    for(int m=0;m<pModel->GetNumSoilLayers();m++) {
      //TMP DEBUG
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for capillary rise algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by capillary rise  algorithm (size of aP[] and aPC[])
//
void CmvSoilBalance::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  nP=0;
  if (_type==SOILBAL_SACSMA)
  { ///< SACSMA MODEL
    aP[nP]="MAX_SAT_AREA_FRAC";     aPC[nP]=CLASS_LANDUSE;  nP++;
    aP[nP]="IMPERMEABLE_FRAC";      aPC[nP]=CLASS_LANDUSE;  nP++;
    aP[nP]="BF_LOSS_FRACTION";      aPC[nP]=CLASS_LANDUSE;  nP++;
    aP[nP]="BASEFLOW_COEFF";        aPC[nP]=CLASS_SOIL;     nP++;
    aP[nP]="SAC_PERC_ALPHA";        aPC[nP]=CLASS_SOIL;     nP++;
    aP[nP]="SAC_PERC_EXPON";        aPC[nP]=CLASS_SOIL;     nP++;
    aP[nP]="SAC_PERC_PFREE";        aPC[nP]=CLASS_SOIL;     nP++;
  }
  else
  {
    ExitGracefully("CmvSoilBalance::GetParticipatingParamList: undefined soil balance algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \details User specifies from and to compartments, levels not known before construction
///
/// \param crise_type [in] Model of capillary rise used
/// \param *aSV [out] Reference to state variable types needed by cr algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by cr algorithm (size of aSV[] and aLev[] arrays)
//
void CmvSoilBalance::GetParticipatingStateVarList(soilbal_type  sb_type,sv_type *aSV, int *aLev, int &nSV)
{
  if(sb_type==SOILBAL_SACSMA)
  {
    nSV=9;
    for (int i=0;i<7;i++){aSV[i]=SOIL; aLev[i]=i; }
    aSV[7]=PONDED_WATER;  aLev[7]=DOESNT_EXIST;
    aSV[8]=SURFACE_WATER; aLev[8]=DOESNT_EXIST;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Distributes ponded water to various soil and/or conceptual subsurface stores
///
/// \param *storage [in] Reference to soil storage from which water rises
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss of water from soil to another soil layer [mm/day]
//
void   CmvSoilBalance::GetRatesOfChange(const double      *state_vars,
                                        const CHydroUnit  *pHRU,
                                        const optStruct   &Options,
                                        const time_struct &tt,
                                              double      *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//no Lakes & glaciers

  if (_type==SOILBAL_SACSMA)
  {
    // Sacramento soil moisture accounting emulation
    //  based upon original FORTRAN code file fland1.F by Dr. Eric Anderson, HRL

    double Adimp       =pHRU->GetSurfaceProps()->max_sat_area_frac; //ADIMP
    double Aimp        =pHRU->GetSurfaceProps()->impermeable_frac; //AIMP
    double cc_frac     =1.0-pHRU->GetSurfaceProps()->bf_loss_fraction; //bf_loss_fraction=1.0-1.0/(1.0+SIDE)
    double uzk         =pHRU->GetSoilProps(1)->baseflow_coeff; //UZK
    double lzpk        =pHRU->GetSoilProps(3)->baseflow_coeff; //LZPK
    double lzsk        =pHRU->GetSoilProps(4)->baseflow_coeff; //LZSK
    double zp          =pHRU->GetSoilProps(1)->SAC_perc_alpha; //ZPERC
    double rexp        =pHRU->GetSoilProps(1)->SAC_perc_expon; //REXP
    double pfree       =pHRU->GetSoilProps(1)->SAC_perc_pfree; //PFREE

    double uzt_stor_max =pHRU->GetSoilCapacity(0); //UZTM
    double uzf_stor_max =pHRU->GetSoilCapacity(1); //UZFM
    double lzt_stor_max =pHRU->GetSoilCapacity(2); //LZTM
    double lzfp_stor_max=pHRU->GetSoilCapacity(3); //LZFPM
    double lzfs_stor_max=pHRU->GetSoilCapacity(4); //LZFSM
    double adimc_stor_max=uzt_stor_max+lzt_stor_max;

    double Aperv = 1.0 - Adimp - Aimp;

    double ponded    =state_vars[pModel->GetStateVarIndex(PONDED_WATER)];
    double uzt_stor  =state_vars[pModel->GetStateVarIndex(SOIL,0)]/Aperv;
    double uzf_stor  =state_vars[pModel->GetStateVarIndex(SOIL,1)]/Aperv;
    double lzt_stor  =state_vars[pModel->GetStateVarIndex(SOIL,2)]/Aperv;
    double lzfp_stor =state_vars[pModel->GetStateVarIndex(SOIL,3)]/Aperv;
    double lzfs_stor =state_vars[pModel->GetStateVarIndex(SOIL,4)]/Aperv;
    double adimc_stor=state_vars[pModel->GetStateVarIndex(SOIL,5)]/Adimp;
    if (Adimp==0){adimc_stor=0;}

    if (ponded<0){ponded=0;} //TMP DEBUG - should not be required.

    //- compute impervious area runoff  ---------------------------------------
    rates[0]=ponded * Aimp/Options.timestep; //PONDED_WATER->SURFACE_WATER
    if (Aimp==1.0){return;}

    //- compute infiltration and runoff amounts  ------------------------------
    double fill=min(ponded,uzt_stor_max- uzt_stor);
    uzt_stor  +=fill;
    adimc_stor+=fill;
    rates[1]=fill*Aperv/Options.timestep;  //PONDED_WATER->UZT
    rates[2]=fill*Adimp/Options.timestep;  //PONDED_WATER->ADIMC

    // determine computational time increments for the basic time interval
    int    N = (int)(1.0 + 0.2 * (uzf_stor + (ponded-fill))); // number of increments in time interval
    N=min(N,100);
    double dt   = Options.timestep / (double)N;   // time increment in days
    double pinc = (ponded-fill)    / (double)N;   // moisture available in each increment

    double duz   = 1.0 - pow((1.0 - uzk ),dt);    // compute free water depletion fractions for the time increment
    double dlzp  = 1.0 - pow((1.0 - lzpk),dt);
    double dlzs  = 1.0 - pow((1.0 - lzsk),dt);

    for(int n=0;n<N;n++)
    {
      //- compute direct runoff (from Adimp area) ----------------------------
      double Frunoff = pow(max(adimc_stor - uzt_stor,0.0) / lzt_stor_max,2);

      rates[0]+=pinc * Adimp * Frunoff /Options.timestep;  //PONDED->SURFACE_WATER

      //- compute baseflow ---------------------------------------------------
      double bfp = dlzp* lzfp_stor;
      lzfp_stor-= bfp;
      rates[3]+=bfp*Aperv*(    cc_frac)/Options.timestep;  //LZFP->SURFACE_WATER
      rates[4]+=bfp*Aperv*(1.0-cc_frac)/Options.timestep;  //LZFP->GW

      double bfs = dlzs*lzfs_stor;
      lzfs_stor-= bfs;
      rates[5]+=bfs*Aperv*(    cc_frac)/Options.timestep;  //LZFS->SURFACE_WATER
      rates[6]+=bfs*Aperv*(1.0-cc_frac)/Options.timestep;  //LZFS->GW

      // if ufz_stor is dry then skip
      if((pinc + uzf_stor) > 0.01)
      {
        //- compute percolation -----------------------------------------------
        double defr; //the lower zone moisture deficiency ratio
        double deficit;
        double perc;

        defr   = max(1.0 - (lzt_stor+lzfp_stor+lzfs_stor)/(lzt_stor_max+lzfp_stor_max+lzfs_stor_max),0.0);
        perc =(dlzp * lzfp_stor_max  + dlzs* lzfs_stor_max) * (uzf_stor / uzf_stor_max);
        perc*=(1.0 + zp * pow(defr,rexp));
        perc = min(perc,uzf_stor); //can't remove more than is there

        deficit=max(lzt_stor_max-lzt_stor,0.0)+max(lzfp_stor_max-lzfp_stor,0.0)+max(lzfs_stor_max-lzfs_stor,0.0);//check-perc
        perc=min(perc,deficit);   //can't fill more than the target

        // split percolated water into the lower zones
        // tension water is filled first except for the pfree area.
        double percf = perc*(    pfree);     // percf is going to free water - splits to percs and percp and excess which goes to LZT
        double perct = perc*(1.0-pfree);     // perct is percolation to tension water

        percf+=max(perct-max(lzt_stor_max-lzt_stor,0.0),0.0); //anything left after filling deficit goes to free
        perct =min(perct,max(lzt_stor_max-lzt_stor,0.0));      //fills tension water deficit

        lzt_stor += perct;
        uzf_stor -= perct;
        rates[7]+=perct *Aperv/Options.timestep;  //UZF->LZT

        if(percf > 0.)
        {
          double hpl;   //relative size of the primary storage  as compared with total lower zone free water storage
          double ratlp; //relative fullness of primary storage
          double ratls; //relative fullness of supplemental storage
          double fracp; //fraction going to primary
          double percs; //amount of excess perc going to supplemental storage
          double percp; //amount of excess perc going to primary storage
          double excess;//leftover excess goes to tension storage

          hpl   = lzfp_stor_max / (lzfp_stor_max + lzfs_stor_max);
          ratlp = lzfp_stor / lzfp_stor_max;
          ratls = lzfs_stor / lzfs_stor_max;
          fracp = min(hpl * 2.0 * (1.0-ratlp) / (1.0-ratlp+1.0-ratls),1.0);

          percs = min(percf *(1- fracp),(lzfs_stor_max-lzfs_stor));
          percp = min(percf - percs    ,(lzfp_stor_max-lzfp_stor));
          excess= max(0.0,percf-percp-percs);

          lzfs_stor += percs;
          uzf_stor  -= percs;
          lzfp_stor += percp;
          uzf_stor  -= percp;
          lzt_stor  += excess;
          uzf_stor  -= excess;

          rates[8]+=percs *Aperv/Options.timestep;  //UZF->LZFS
          rates[9]+=percp *Aperv/Options.timestep;  //UZF->LZFP
          rates[7]+=excess*Aperv/Options.timestep;  //UZF->LZT
        } //end if percf>0

        //- compute interflow --------------------------------------------------
        double delta = max(duz * uzf_stor,0.0);
        uzf_stor -= delta;
        rates[10]+=delta*Aperv/Options.timestep; //UZF->SURFACE_WATER

        //- compute overflow from full UZF--------------------------------------
        double overflow=max(pinc-(uzf_stor_max-uzf_stor),0.0);
        pinc-=overflow;
        rates[0]+=overflow*Aperv                /Options.timestep; //PONDED_WATER->SURFACE_WATER
        rates[0]+=overflow*Adimp*(1.0 - Frunoff)/Options.timestep; //PONDED_WATER->SURFACE_WATER

        // check - what if Frunoff >1?
        // overflow * Frunoff = amount of surface runoff which comes from that portion of Adimp which is not currently generating direct runoff.
      }

      uzf_stor   += pinc;  //any uninfiltrated water that didnt overflow
      rates[11]  += pinc*Aperv/Options.timestep; // PONDED->UZF

      //JRC: NOT SURE THIS IS WHOLLY CONSISTENT WITH SUM_OVERFLOW CALCULATION ABOVE
      //shouldn't there be something else for water if adimc_stor is full and overflowing?? (not a huge problem - this will stay as ponded water)
      double to_adimc=min(pinc* (1.0 - Frunoff),max(adimc_stor_max-adimc_stor,0.0));
      adimc_stor += to_adimc;
      rates[2]   += to_adimc*Adimp/Options.timestep; // PONDED->ADIMC

    } // end for (int n=0;n<N;n++)

    // JRC: if this is not automatically satisfied, then this is a mass balance error, by definition - nowhere to get this water from
    // BUT, this was included in original SAC-SMA code.
    // if(adimc_stor < uzt_stor) { adimc_stor = uzt_stor; }
    // rates[12]+=max(uzt_stor-adimc_stor,0.0)*Adimp/Options.timestep; //VIOLATES MASS BALANCE
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep.
///
/// \param *storage [in] state variable array for this HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from soil to other soil layer [mm/day]
//
void   CmvSoilBalance::ApplyConstraints(const double      *storage,
                                        const CHydroUnit  *pHRU,
                                        const optStruct   &Options,
                                        const time_struct &tt,
                                              double      *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  //Handled internally to soil balance routine
}
