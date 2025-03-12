/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
------------------------------------------------------------------
*/

#include "GWRiverConnection.h"
#include "MFUSGpp.h"

/*****************************************************************
   GWRiverConnection Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CGWRiverConnection constructor
/// \param pGWModel [in] pointer to groundwater model
//
CGWRiverConnection::CGWRiverConnection(CGroundwaterModel* pGWModel)
{
  _nSeg      = 0;
  _NB        = 0;
  _aNodes    = NULL;
  _aSubBasin = NULL;
  _aWeights  = NULL;
  _aLength   = NULL;
  _aElevs    = NULL;
  _aCond     = NULL;
  _aGWFlux   = NULL;
  _anSegmentsBySB = NULL;
  _aSegmentsBySB  = NULL;

  _pGWModel=pGWModel;

  _mode = RIVER;
  _condtype = CONDUCTANCE;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CGWRiverConnection::~CGWRiverConnection()
{
  delete [] _aNodes   ;
  delete [] _aSubBasin;
  delete [] _aWeights ;
  delete [] _aLength  ;
  if (_aElevs != NULL) {
    for (int i = 0; i < _nSeg; i++) { delete[] _aElevs[i]; delete[] _aCond[i]; } delete[] _aElevs; delete[] _aCond;
  }
  delete [] _aGWFlux  ;
  delete [] _anSegmentsBySB;
  delete [] _aSegmentsBySB;
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes arrays
/// \param nSegments [in] number of segments
/// \param NB        [NB] number of subbasins in model
//
void CGWRiverConnection::Initialize(int nSegments, int NB)
{
  _nSeg      = nSegments;
  _NB        = NB;
  _aNodes    = new int*    [_nSeg];
  _aSubBasin = new int     [_nSeg];
  _aWeights  = new double* [_nSeg];
  _aLength   = new double  [_nSeg];
  _aGWFlux   = new double  [_nSeg];
  _aElevs    = new double *[_nSeg];
  _aCond     = new double *[_nSeg];
  for (int i = 0; i < _nSeg; i++) {
    _aElevs[i] = new double[2];
    _aCond [i] = new double[2];
  }
  _anSegmentsBySB = new int  [NB];
  _aSegmentsBySB  = new int* [NB];

  //-- Set segments in basin counter to zero
  for (int i=0; i<NB; i++){
    _anSegmentsBySB[i] = 0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Copies values read in Raven into the memory for the
///  MODFLOW PBJ Package. Calls FORTRAN subroutines.
//
void CGWRiverConnection::InitializePBJ()
{
#ifdef _MODFLOW_USG_
  int iseg, i, p, iSB;
  int *pNBCounter;
  //-- Process mode, condtype. Explicit recast
  int pbjmode = 2;
  int condtype = 0;
  if (_mode == DRAIN) {pbjmode = 1;}
  if (_condtype == UNITCOND) {condtype = 1;}

  //-- Initialize aSegmentsBySB for each subbasin
  //-- & counter for filling operation
  pNBCounter = new int[_NB];
  for (i=0; i<_NB; i++) {
    _aSegmentsBySB[i] = new int[_anSegmentsBySB[i]];
    pNBCounter[i] = 0;
  }

  //-- Loop over all segments, filling aSegmentsBySB
  // We're using an array to count so we can fill all SB slots at once
  for (iseg=0; iseg<_nSeg; iseg++) {
    p = _aSubBasin[iseg];
    iSB = pNBCounter[p];
    _aSegmentsBySB[p][iSB] = iseg;
    pNBCounter[p] += 1;
  }
  delete [] pNBCounter;

  MFUSG::init_pbj_settings(&_nSeg, &pbjmode, &condtype);

  //-- Add Segments (One-by-one, since Fortran & C store 2D arrays differently)
  for (i=0; i<_nSeg; i++) {
    iseg = i+1;
    MFUSG::add_pbj_segment(&iseg, &_aNodes[i][0],   &_aNodes[i][1],   &_aNodes[i][2],
                                  &_aWeights[i][0], &_aWeights[i][1], &_aWeights[i][2],
                                  &_aWeights[i][3], &_aWeights[i][4], &_aWeights[i][5],
                                  &_aElevs[i][0], &_aElevs[i][1],
                                  &_aCond [i][0], &_aCond [i][1],
                                  &_aLength[i]);
  }
  //-- Let PBJ package finish up
  MFUSG::init_pbj_values();
#else
  ExitGracefully("CGWRiverConnection::Initialize PBJ connection: PBJ can only be used in Raven w/Modflow USG", RUNTIME_ERR);
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Updates stream segment elevations in MODFLOW
/// \details Currently passing the depth to MF and letting it calculate the elevation.
///
/// \param p      [in] subbasin global index
/// \param pBasin [in] pointer to subbasin
///
//
void CGWRiverConnection::UpdateRiverLevelsBySB(const int p, const CSubBasin *pBasin) const
{
  //-- Only happens for River Mode - Leave if Drain
  if (_mode == DRAIN) {return;}
#ifdef _MODFLOW_USG_
  double    river_level;

  //-- Set values
  river_level = pBasin->GetRiverDepth();

  //-- Pass to MF
  // Does _aSegmentsBySB need to be contiguous memory???
  MFUSG::update_pbj_by_basin(&_anSegmentsBySB[p], _aSegmentsBySB[p], &river_level);
#endif

}

//////////////////////////////////////////////////////////////////
/// \brief Updates _aGWFlux array with latest fluxes from MODFLOW
//
void CGWRiverConnection::UpdateRiverFlux() const
{
#ifdef _MODFLOW_USG_
  MFUSG::get_pbj_segment_flows(_aGWFlux);
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Updates _aGWFlux array with latest fluxes from MODFLOW
///
/// \param p [in] SubBasin global index
//
double CGWRiverConnection::CalcRiverFlowBySB(const int p) const
{
  int    iseg;
  double sumflux = 0.0;

  if (_nSeg == 0) {return 0.0; }

  for (int i=0; i < _anSegmentsBySB[p]; i++)
  {
    iseg = _aSegmentsBySB[p][i];
    sumflux += _aGWFlux[iseg];
  }

  return(sumflux);
}

//////////////////////////////////////////////////////////////////
/// \brief Sets whether the segments will act as Rivers or Drains in MODFLOW
///
/// \param s [in] name of calculation type
//
void CGWRiverConnection::SetMode(const char *s)
{
  if      (!strcmp(s, "RIVER")){ _mode = RIVER;}
  else if (!strcmp(s, "DRAIN")){ _mode = DRAIN;}
  else {
    ExitGracefully("CGWRiverConnection: Invalid Type",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets whether the conductance is per segment or per segment unit length
/// \param s [in] name of conductance type
//
void CGWRiverConnection::SetCondType(const char *s)
{
  if      (!strcmp(s, "CONDUCTANCE")){ _condtype = CONDUCTANCE;}
  else if (!strcmp(s,        "COND")){ _condtype = CONDUCTANCE;}
  else if (!strcmp(s,    "UNITCOND")){ _condtype = UNITCOND;   }
  else {
    ExitGracefully("CGWRiverConnection: Invalid Condutance Type",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Adds a river segment to the object. Called during parsing.
///
/// \param segid [in] sequential number of segment
/// \param *nodes [in] pointer to nodeIDs of segment
/// \param *weights [in] pointer to barycentric weights of segment start/end
/// \param *elev [in] pointer to elevation of segment start/end
/// \param *cond [in] pointer to conductance of segment start/end
/// \param length [in] segment length
/// \param SBID [in] BasinID associated with segment
//
void CGWRiverConnection::AddSegment(const int segid, const int *nodes, const double *weights,
                                 const double *elev, const double *cond, const double length,
                                 const int p)
{
  //-- Initialize
  int i;
  _aNodes  [segid] = new int[3];
  _aWeights[segid] = new double[6];

  //-- Fill
  _aLength  [segid] = length;
  _aSubBasin[segid] = p;
  for (i=0; i<3; i++) {
    _aNodes  [segid][i]   = nodes  [i];
    _aWeights[segid][i]   = weights[i];
    _aWeights[segid][i+3] = weights[i+3];
  }
  for (i=0; i<2; i++) {
    _aElevs  [segid][i]   = elev[i];
    _aCond   [segid][i]   = cond[i];
  }

  //-- Increment
  _anSegmentsBySB[p]++;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the number of river segments involved in the process
/// \return Number of river segments
//
int CGWRiverConnection::GetNumSegments() {
  return _nSeg;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns pointer to array of segment indexes for the given Basin
/// \return Segments present in subbasin p
///
/// \param p [in] SubBasin index
//
int* CGWRiverConnection::GetSubBasinSegments(int p) const
{
  return _aSegmentsBySB[p];
}
