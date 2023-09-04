#include "RiverReach.h"

//-- River Reach Constructor
CRiverReach::CRiverReach(CModel *pModel, CGroundwaterModel *pGWModel, int reachNo, double length,
                         double wp, double river_bot, LeakanceMethod leakance_method)
{
  //-- Values passed
  _pModel          = pModel;
  _pGWModel        = pGWModel;
  _reachNo         = reachNo;
  _length          = length;
  _wettedPerimeter = wp;
  _leakance_method = leakance_method;

  //-- Create Leakance Object
  if (_leakance_method == LEAK_CONSTANT) {
    _pRiverLeakance = new CRiverLeakanceABC;
  }
  else if (_leakance_method == LEAK_MOREL_SEYTOUX) {
    _pRiverLeakance = new CRiverLeakanceMorelSeytoux;
  }
  else {
    // MODFLOW as Default
    _pRiverLeakance = new CRiverLeakanceMF;
  }
}

//-- River Reach Destructor
CRiverReach::~CRiverReach()
{
  delete _pRiverLeakance; _pRiverLeakance = NULL;
}

double CRiverReach::_calcConductance()
{
  double leakance = _pRiverLeakance->getLeakance();
  return leakance * _length * _wettedPerimeter;
}

void CRiverReach::addNode(int nodeid)
{
  //-- Add node id to vector
  _vNodes.push_back(nodeid);
}

vector<int>* CRiverReach::getNodes()
{
  return &_vNodes;
}

double CRiverReach::getDischarge(double head_cell, double head_river, double tstep)
{
  double q    = 0.0;
  double cond = _calcConductance();
  //-- Gaining river (cell water to river)
  if (head_cell > _river_bot){
    q = cond * (head_river - head_cell);
  }
  else if (head_cell <= _river_bot){
    //-- Losing river (river percolation to cell)
    q = cond * (head_river - _river_bot);
  }
  return q;
}

/////////////////////////////////////////////////////////////////
//-- River Leakance Base Class

void CRiverLeakanceABC::initialize(double leakance)
{
  _leakance = leakance;
}

double CRiverLeakanceABC::getLeakance()
{
  return _leakance;
}

void CRiverLeakanceABC::update() {}

/////////////////////////////////////////////////////////////////
//-- River Leakance MODFLOW
void CRiverLeakanceMF::initialize(double k_riverbed, double thickness)
{
  _leakance = k_riverbed / thickness;
}

double CRiverLeakanceMF::getLeakance()
{
  return _leakance;
}
