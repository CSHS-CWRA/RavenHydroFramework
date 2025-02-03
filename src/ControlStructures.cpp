/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#include "ControlStructures.h"
#include "Model.h"
//////////////////////////////////////////////////////////////////
/// \brief control structure constructor, creates empty featureless outflow control structure
/// \param name [in] structure name
/// \param SBID [in] subbasin ID of reservoir
/// \param downID [in] downstream subbasin
/// \notes assumes SBID basin has valid reservoir
//
CControlStructure::CControlStructure(string name, const long long SBID,const long long downID) {
  _name=name;
  _SBID=SBID;
  _target_SBID=downID;
  _dRefElev=0.0;
  _nRegimes=0;
  _aRegimes=NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief control structure destructor
//
CControlStructure::~CControlStructure()
{
  delete [] _aRegimes;
}
//////////////////////////////////////////////////////////////////
/// \brief returns control structure name
/// \returns name of control structure
//
string CControlStructure::GetName         () const {return _name;}


//////////////////////////////////////////////////////////////////
/// \brief returns target subbasin ID to which outflow is directed
/// \returns target subbasin ID
//
long long CControlStructure::GetTargetBasinID() const {return _target_SBID;}

//////////////////////////////////////////////////////////////////
/// \brief sets target subbasin ID to which outflow is directed
/// \param SBID [in] - valid subbasin ID (or DOESNT_EXIST, if flow is directed outside model)
//
void   CControlStructure::SetTargetBasin(const long long SBID)
{
  _target_SBID=SBID;
}

//////////////////////////////////////////////////////////////////
/// \brief sets downstream reference elevation for control structure
/// \param dRefElev [in] - downstream reference elevation [m]
//
void   CControlStructure::SetDownstreamRefElevation(const double dRefElev)
{
  _dRefElev=dRefElev;
}

//////////////////////////////////////////////////////////////////
/// \brief adds outflow regime to list of regime conditions
/// \param pRegime [in] - pointer to valid outflow regime
//
void CControlStructure::AddRegime(COutflowRegime* pRegime)
{
  if (!DynArrayAppend((void**&)(_aRegimes),(void*)(pRegime),_nRegimes)){
    ExitGracefully("CControlStructure::AddRegime: trying to add NULL regime",BAD_DATA);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief returns current regime name
/// \notes deterimines which regime is applicable AND applies regime flow constraints
/// \param tt      [in] - time structure corresponding to END of current timestep
/// \returns current regime name
//
string CControlStructure::GetCurrentRegimeName(const time_struct &tt) const {
  for (int i = 0; i < _nRegimes; i++)
  {
    if (_aRegimes[i]->AreConditionsMet(tt))
    {
      return _aRegimes[i]->GetName();
    }
  }
  return "NONE";
}

//////////////////////////////////////////////////////////////////
/// \brief returns outflow from control structure at end of timestep
/// \notes deterimines which regime is applicable AND applies regime flow constraints
/// the priority goes from first operating  regime to last - the first whose conditions are satisfied
/// \param stage   [in] - reservoir stage at end of timestep
/// \param Q_start [in] - total reservoir outflow at end of timestep
/// \param Q_end   [in] - current best estimate of total resevoir outflow at end of timestep
/// \param tt      [in] - time structure corresponding to END of current timestep
/// \returns outflow from control structure at end of timestep
//
double CControlStructure::GetOutflow(const double& stage, const double& stage_start, const double& Q_start, const time_struct& tt) const
{
  for (int i = 0; i < _nRegimes; i++)
  {
    if (_aRegimes[i]->AreConditionsMet(tt))
    {
      return _aRegimes[i]->GetOutflow(stage,stage_start,Q_start,_target_SBID,_dRefElev);
    }
  }
  return 0.0; //no regime active, no flow
}

//////////////////////////////////////////////////////////////////
/// \brief outflow regime constructor, creates empty featureless outflow regime
/// \param name [in] structure name
/// \param pMod  [in] pointer to model
//
COutflowRegime::COutflowRegime(const CModel *pMod, const string name)
{
  _pModel=pMod;
  _name=name;

  _nConditions=0;
  _pConditions=NULL;
  _nConstraints=0;
  _pConstraints=NULL;
  _pCurve=NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief outflow regime destructor
//
COutflowRegime::~COutflowRegime()
{
  delete [] _pConstraints;
  delete [] _pConditions;
}

//////////////////////////////////////////////////////////////////
/// \brief get regime name
/// \returns regime name
//
string COutflowRegime::GetName() const
{
  return _name;
}
//////////////////////////////////////////////////////////////////
/// \brief get stage discharge relation
/// \returns stage discharge relation
//
const CStageDischargeRelationABC* COutflowRegime::GetCurve() const {
  return _pCurve;
}
//////////////////////////////////////////////////////////////////
/// \brief checks all conditions and constraints, issues errors and warnings as needed
//
void    COutflowRegime::CheckConditionsAndConstraints() const
{
  string warn;
  for (int i = 0; i < _nConditions; i++) {

  }
  for (int j = 0; j < _nConstraints; j++) {
    string      var=_pConstraints[j]->variable;
    if ((var=="FLOW") && (_pConstraints[j]->compare_val1<=0)){
        WriteWarning("COutflowRegime::CheckConditionsAndConstraints: FLOW Constraints must be greater than zero",false);
    }
    if (!((var=="FLOW_DELTA") || (var=="FLOW"))){
        warn="COutflowRegime::CheckConditionsAndConstraints: unrecognized constraint "+var;
        ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief specifies active stage discharge relation for outflow regime
/// \param pCurv [in] pointer to valid stage discharge relation
//
void    COutflowRegime::SetCurve(CStageDischargeRelationABC *pCurv)
{
  ExitGracefullyIf(pCurv==NULL,"COutflowRegime::SetCurve: null curve",BAD_DATA_WARN);
  _pCurve=pCurv;
}

//////////////////////////////////////////////////////////////////
/// \brief adds regime condition to outflow regime
/// \param pCond [in] pointer to valid regime condition
//
void    COutflowRegime::AddRegimeCondition(RegimeCondition* pCond)
{
  if (!DynArrayAppend((void**&)(_pConditions),(void*)(pCond),_nConditions)){
    ExitGracefully("COutflowRegime::AddRegimeCondition: trying to add NULL condition",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief adds regime constraint to outflow regime
/// \param pCond [in] pointer to valid regime constraint
//
void    COutflowRegime::AddRegimeConstraint(RegimeConstraint* pCond)
{
  if (!DynArrayAppend((void**&)(_pConstraints),(void*)(pCond),_nConstraints)){
    ExitGracefully("COutflowRegime::AddRegimeConstraint: trying to add NULL constraint",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief evaluates general logical conditional statement
/// \param cond [in] comparison type
/// \param v [in] value to be compared to v1 (and sometimes also v2)
/// \param v1 [in] comparator value
/// \param v2 [in] secondary comparator value (only used for COMPARE_BETWEEN)
/// \returns true if comparison evaluates to true, false otherwise
//
bool EvaluateCondition(const comparison cond, const double& v, const double& v1, const double& v2)
{
  if      (cond==COMPARE_IS_EQUAL   ){return fabs(v-v1)<REAL_SMALL; }
  else if (cond==COMPARE_NOT_EQUAL  ){return fabs(v-v1)>REAL_SMALL; }
  else if (cond==COMPARE_GREATERTHAN){return (v>v1);  }
  else if (cond==COMPARE_LESSTHAN   ){return (v<v1);  }
  else if (cond==COMPARE_BETWEEN    ){return (v>=v1) && (v<=v2); }//inclusive
  return false;
}

//////////////////////////////////////////////////////////////////
/// \brief returns true if regime conditions are met
/// \param tt [in] time structure corresponding to END of timestep
/// \returns true if regime conditions are met
//
bool      COutflowRegime::AreConditionsMet(const time_struct& tt) const
{
  for (int i = 0; i < _nConditions; i++)
  {
    comparison comp=_pConditions[i]->compare;
    long long  SBID=_pConditions[i]->basinID;
    double v1=_pConditions[i]->compare_val1;
    double v2=_pConditions[i]->compare_val2;
    string var=_pConditions[i]->variable;

    if      (var == "DATE") {
      if (!EvaluateCondition(comp, tt.model_time, v1, v2)) {return false;}
    }
    else if (var == "YEAR") {
      if (!EvaluateCondition(comp, (double)(tt.year), v1, v2)) {return false;}
    }
    else if ((var == "DOY") || (var=="DAY_OF_YEAR")) {
      if (comp!=COMPARE_BETWEEN){
        if (!EvaluateCondition(comp, tt.julian_day, v1, v2)) {return false;}
      }
      else {
        if (!(
             ((tt.julian_day >= v1) && (tt.julian_day <= v2)) || //regular between
             (((v2 < v1) && ((tt.julian_day <= v2) || (tt.julian_day >= v1)))) //wrap between
           ))
        {
          return false;
        }
      }
    }
    else if (var=="STAGE"){
      double h=_pModel->GetSubBasinByID(SBID)->GetReservoir()->GetResStage();
      if (!EvaluateCondition(comp, h, v1, v2)) {return false;}
    }
    else if (var == "FLOW") {
      double Q=_pModel->GetSubBasinByID(SBID)->GetOutflowRate();
      if (!EvaluateCondition(comp, Q, v1, v2)) {return false;}
    }
    else if (var == "RIVER_DEPTH") {
      double h=_pModel->GetSubBasinByID(SBID)->GetRiverDepth();
      if (!EvaluateCondition(comp, h, v1, v2)) {return false;}
    }
    else if (var == "WATER_LEVEL") {
      double h=_pModel->GetSubBasinByID(SBID)->GetWaterLevel();
      if (!EvaluateCondition(comp, h, v1, v2)) {return false;}
    }
    else if (var == "STAGE_CHANGE") {
      double h     =_pModel->GetSubBasinByID(SBID)->GetReservoir()->GetResStage();
      double h_last=_pModel->GetSubBasinByID(SBID)->GetReservoir()->GetOldStage();
      double r=(h-h_last)/_pModel->GetOptStruct()->timestep;
      if (!EvaluateCondition(comp, r, v1, v2)) {return false;}
    }
    else {
      string warn="COutflowRegime::AreConditionsMet: unrecognized condition variable "+var+" in :Condition command";
      ExitGracefully(warn.c_str(),BAD_DATA);
    }
  }
  return true; //only reach here if all conditions met (or if zero conditions)
}

//////////////////////////////////////////////////////////////////
/// \brief returns outflow from control structure for given stage, with prioritized constraints applied
/// \param h [in] stage, in [m]
/// \param Q_start [in] control structure outflow at start of timestep [m3/s]
/// \returns outflow, in [m3/s]
//
double          COutflowRegime::GetOutflow(const double &h, const double &h_start, const double &Q_start, const long long &target_SBID, const double &drefelev) const
{
  double tstep    = _pModel->GetOptStruct()->timestep;
  double rivdepth = _pModel->GetSubBasinByID(target_SBID)->GetRiverDepth();

  double Q = _pCurve->GetDischarge(h, h_start, Q_start, rivdepth, drefelev);

  double dQdt = (Q-Q_start)/tstep; //m3/s/d

  //apply constraints here
  for (int j = _nConstraints-1; j >=0; j--) //priority from first to last
  {
    string      var=_pConstraints[j]->variable;
    comparison comp=_pConstraints[j]->compare_typ;
    double       v1=_pConstraints[j]->compare_val1;
    double       v2=_pConstraints[j]->compare_val2;

    if      (var== "FLOW") {

      if      (comp == COMPARE_GREATERTHAN) {if (Q<=v1){Q=v1;}}
      else if (comp == COMPARE_LESSTHAN   ) {if (Q>=v1){Q=v1;}}
      else if (comp == COMPARE_BETWEEN    ) {if (Q<=v1){Q=v1;} if (Q>=v2){Q=v2;}}
      dQdt=(Q-Q_start)/tstep;
    }
    else if (var == "FLOW_DELTA") {
      if      (comp == COMPARE_GREATERTHAN) {if (dQdt<=v1){dQdt=v1;}}
      else if (comp == COMPARE_LESSTHAN   ) {if (dQdt>=v1){dQdt=v1;}}
      else if (comp == COMPARE_BETWEEN    ) {if (dQdt<=v1){dQdt=v1;} if (dQdt>=v2){dQdt=v2;}}
      Q=Q_start+dQdt*tstep;
    }
    else{
      ExitGracefully("COutflowRegime::GetOutflow: invalid state in constraint (should be FLOW or FLOW_DELTA)",BAD_DATA);
    }
  }
  return Q;
}
