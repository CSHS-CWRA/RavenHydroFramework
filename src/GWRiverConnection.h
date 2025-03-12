/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------*/

#ifndef GWRIVERCONN_H
#define GWRIVERCONN_H

#include <set>
#include <map>
#include "RavenInclude.h"
#include "GroundwaterModel.h"

class CGWRiverConnection
{
protected:/*-------------------------------------------------------*/
  // Local Enumerated Types
  enum calc_mode { RIVER, DRAIN };
  enum cond_type { CONDUCTANCE, UNITCOND };

  CGroundwaterModel* _pGWModel;              ///< pointer to groundwater model
  calc_mode     _mode;
  cond_type     _condtype;
  int           _nSeg;                      ///< Number of gw model segments
  int         **_aNodes;                    ///< 2D array indices of segment nodes (3x each)
  int          *_aSubBasin;                 ///< array of subbasin index for each segment
  double      **_aWeights;                  ///< 2D array of node barycentric weights (3x start, 3x end)
  double       *_aLength;                   ///< array of Segment lengths (only one per segment!)
  double      **_aElevs;                    ///< 2D array of elevations (segment start, end)
  double      **_aCond;                     ///< 2D array of conductances (start, end)
  double       *_aGWFlux;                   ///< array of GW fluxes to river, by segment

  int           _NB;                        ///< Local storage of subbasin count
  int*          _anSegmentsBySB;            ///< Number of segments in each subbasin [NB] (SB index, not SBID)
  int**         _aSegmentsBySB;             ///< 2D array of segments, arranged as [NB][segments]

public: /*-------------------------------------------------------*/
  // Constructors/Destructor
  CGWRiverConnection(CGroundwaterModel* pGWModel);
  ~CGWRiverConnection();

  void Initialize(int nSegments, int NB);
  void InitializePBJ();

  // Getters
  int GetNumSegments();
  int* GetSubBasinSegments(int p) const;
    // Setters
  void SetMode      (const char *s);
  void SetCondType  (const char *s);
  void AddSegment   (const int segid, const int *nodes, const double *weights,
                     const double *elev, const double *cond, const double length,
                     const int p);

  // General Methods
  void   UpdateRiverLevelsBySB(const int p, const CSubBasin *pBasin) const;
  void   UpdateRiverFlux      () const;
  double CalcRiverFlowBySB    (const int p) const;
};

#endif
