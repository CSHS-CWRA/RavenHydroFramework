#include "RavenInclude.h"
#include "GroundwaterModel.h"
#include <vector>
//#include "TimeSeries.h"

#ifndef RIVERREACH_H
#define RIVERREACH_H

//-- Declaration of classes for use in CRiverReach (Definitions below)
enum LeakanceMethod
{
  LEAK_CONSTANT     ,      ///< Leakance coefficient given as a single value
  LEAK_MODFLOW      ,      ///< MODFLOW River leakance, hydraulic conductivity of clogging layer over its thickness
  LEAK_MOREL_SEYTOUX,      ///< Morel-Seytoux et al. (2018) analytically derived, physically based coefficient
};

class CRiverLeakanceABC;

class CRiverReach
{
private:/*------------------------------------------------------*/
  CModel                *_pModel;
  CGroundwaterModel     *_pGWModel;
  LeakanceMethod         _leakance_method;
  CRiverLeakanceABC     *_pRiverLeakance;
  vector<int>            _vNodes;
  int                    _nNodes;
  int                    _reachNo;
  double                 _length;
  double                 _wettedPerimeter;
  double                 _river_bot;

  //-- Private methods
  double _calcConductance ();

public:/*-------------------------------------------------------*/
  //-- Constructor/Destructor
  CRiverReach              (CModel *pModel, CGroundwaterModel *pGWModel, int reachNo, double length,
                            double wp, double river_bot, LeakanceMethod leakance_method);
 ~CRiverReach              ();

  //-- Setters
  void         addNode          (int nodeid);
  //-- Accessors
  vector<int>* getNodes         ();
  double       getDischarge     (double head_cell, double head_river, double tstep);
};


class CRiverLeakanceABC
{
private:/*------------------------------------------------------*/
  double _leakance;

public:/*-------------------------------------------------------*/
  void   initialize(double leakance);
  double getLeakance();
  void   update();  // For leakance methods where the values vary in time
};

class CRiverLeakanceMF : public CRiverLeakanceABC
{
private:/*------------------------------------------------------*/
  double _leakance;
  //double _k_riverbed;       ///< Riverbed "Clogging Layer" Hydraulic Conductivity [m/d]
  //double _thickness;        ///< Riverbed "Clogging Layer" Thickness [m]
public:/*-------------------------------------------------------*/
  void   initialize(double k_riverbed, double thickness);
  double getLeakance();
};

class CRiverLeakanceMorelSeytoux : public CRiverLeakanceABC
{
private:/*------------------------------------------------------*/
  double _seepage_discharge;
  double _length;
  double _hk;
public:/*-------------------------------------------------------*/
};

#endif
