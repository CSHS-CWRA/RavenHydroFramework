/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2026 the Raven Development Team
----------------------------------------------------------------
CAgeTracer routines related to age tracer transport

Raven transports is in terms of days after introduction to system as precip
Spinup period needed to handle initialization of age
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "AgeTracers.h"

//////////////////////////////////////////////////////////////////
/// \brief age tracer model constructor
//
CAgeTracer::CAgeTracer(CModel* pMod,CTransportModel* pTMod,string name,const int c)
  :CConstituentModel(pMod,pTMod,name,AGE_TRACER,false,c)
{
  _t_diff=0.0;
  _model_time=0.0;
}

//////////////////////////////////////////////////////////////////
/// \brief age tracer  model destructor
//
CAgeTracer::~CAgeTracer(){}

//////////////////////////////////////////////////////////////////
/// \brief initializes age tracer transport model
//
void CAgeTracer::Initialize(const optStruct& Options)
{
  // Initialize base class members
  //--------------------------------------------------------------------
  CConstituentModel::Initialize(Options);
}
//////////////////////////////////////////////////////////////////
/// \brief initializes isotope transport model
//
void CAgeTracer::UpdateInitialConditions(const optStruct& Options)
{
  //Reasonable estimates of mean lake age = V/Q
  CSubBasin *pBasin;
  CHydroUnit *pHRU;
  int nSB=_pModel->GetNumSubBasins();
  double Qout,Cinit;
  for(int p=0;p<nSB;p++)
  {
    pBasin=_pModel->GetSubBasin(p);
    if ((pBasin->GetReservoir()!=NULL) && (_aMres[p]==0.0))//if not initialized from .rvc
    {
      Qout=pBasin->GetReservoir()->GetOutflowRate();
      Cinit=0.0;
      if (Qout>0){
        //Cinit=pBasin->GetReservoir()->GetStorage()/Qout; //age ~=mean residence time
      }

      _aMres     [p]=Cinit*pBasin->GetReservoir()->GetStorage();
      _aMres_last[p]=_aMres[p];

      _aMsed     [p]=Cinit*_pModel->GetHydroUnit(pBasin->GetReservoir()->GetHRUIndex())->GetArea()*M2_PER_KM2*pBasin->GetReservoir()->GetLakebedThickness();
      _aMsed_last[p]=_aMsed[p];

      _aMout_res      [p]=Cinit*pBasin->GetReservoir()->GetOutflowRate()*SEC_PER_DAY;
      _aMout_res_last [p]=Cinit*pBasin->GetReservoir()->GetOldOutflowRate()*SEC_PER_DAY;
    }
  }

  //convert mass (representing entry time relative to previous model start time) to new t==0 time datum:
  int i,iStor;

  double V,MassAge;
  for (int k=0;k<_pModel->GetNumHRUs();k++){
    pHRU=_pModel->GetHydroUnit(k);

    for (int ii=0; ii<_pTransModel->GetNumWaterCompartments();ii++){
      iStor=_pTransModel->GetStorWaterIndex(ii);
      V=pHRU->GetStateVarValue(iStor);

      i=_pTransModel->GetStorIndex(_constit_index,ii);

      MassAge=pHRU->GetStateVarValue(i)-V*_t_diff;
      pHRU->SetStateVarValue(i,MassAge);
    }
  }

  double Vres;
  for(int p=0;p<nSB;p++)
  {
    pBasin=_pModel->GetSubBasin(p);
    if ((pBasin->GetReservoir()!=NULL) && (_aMres[p]==0.0))//if not initialized from .rvc
    {
      V=pBasin->GetReservoir()->GetStorage();
      _aMres     [p]-=V*_t_diff;

      V=pBasin->GetReservoir()->GetOldStorage();
      _aMres_last[p]-=V*_t_diff;

      Vres=_pModel->GetHydroUnit(pBasin->GetReservoir()->GetHRUIndex())->GetArea()*M2_PER_KM2*pBasin->GetReservoir()->GetLakebedThickness();
      _aMsed     [p]-=Vres*_t_diff;
      _aMsed_last[p]-=Vres*_t_diff;

      _aMout_res      [p]-=_t_diff*pBasin->GetReservoir()->GetOutflowRate()   *SEC_PER_DAY;
      _aMout_res_last [p]-=_t_diff*pBasin->GetReservoir()->GetOldOutflowRate()*SEC_PER_DAY;
    }
  }
}
void   CAgeTracer::SetModelTime(const double &t){
  _model_time=t;
}
void   CAgeTracer::SetInitialTimeOffset(const double &tdiff){
  _t_diff=tdiff;
}

//////////////////////////////////////////////////////////////////
/// \brief returns age (rather than model time of entry, which is what is tracked)
/// \param M [in] mass [mg/m2]
/// \param V [in] volume [mm]
//
double CAgeTracer::CalculateReportingConcentration(const double &mass,const double &vol) const
{
  if(fabs(vol)<1e-6) {return 0.0;} //{return RAV_BLANK_DATA;}
  double C=CConstituentModel::CalculateReportingConcentration(mass,vol); //d
  return _model_time-C; //d
}
//////////////////////////////////////////////////////////////////
/// \brief returns outflow age from subbasin p at current point in time
/// \notes only used for reporting; calculations exclusively in terms of mass/energy
/// only called at end of timestep for output/diagnostics
/// \param p [in] subbasin index
//
double CAgeTracer::GetOutflowConcentration(const int p) const
{
  double C=CConstituentModel::GetOutflowConcentration(p); //d
  if (C==0.0){return 0.0;}
  return _model_time-C;
}
//////////////////////////////////////////////////////////////////
/// \brief Converts age units to days after entry
/// \param age [in] - age of water, [d]
/// \returns days after entry
//
double CAgeTracer::ConvertConcentration(const double &age) const
{
  return age+_model_time; //[d after entry]
}

//////////////////////////////////////////////////////////////////
/// \brief returns advection correction factor for age
/// \param pHRU       [in] pointer to HRU of interest
/// \param iFromWater [in] index of "from" water storage state variable
/// \param iToWater   [in] index of "to" water storage state variable
/// \param mass       [in] source mass prior to advection [mg]
/// \param vol        [in] source volume prior to advection [mm] or [mm-m2]
/// \param Q          [in] flow between source and recipient [mm/d] or [mm-m2]
///
double CAgeTracer::GetAdvectionCorrection(const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& mass, const double &vol, const double &Q) const
{
// does nothing for now. Only needed if newer or older waters should exit first.
  return 1.0;
}
//////////////////////////////////////////////////////////////////
/// \brief Write transport output file headers in .tb0 format
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CAgeTracer::WriteEnsimOutputFileHeaders(const optStruct& Options)
{
  this->WriteOutputFileHeaders(Options); //Ensim format not supported by age tracer model
}

//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output to ensim formatted file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams; called from CModel::WriteMinorOutput()
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
//
void CAgeTracer::WriteEnsimMinorOutput(const optStruct& Options,const time_struct& tt)
{
  this->WriteMinorOutput(Options,tt);//Ensim format not supported by age tracer model
}
