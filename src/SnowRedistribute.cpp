/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2025 the Raven Development Team
------------------------------------------------------------------
LatRedistribute (redistribute snow based on slope and snow SWE)
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "LatConnect.h"

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

CmvLatRedistribute::CmvLatRedistribute(int sv_ind,
                                       double max_snow_height,
                                       redist_method method,
                                       CModel *pModel)
                   : CLateralExchangeProcessABC(LAT_REDISTRIBUTE, pModel)
{
    _max_snow_height=max_snow_height;
    _method=method;
    _iRedistributeFrom = sv_ind;
    _iRedistributeTo   = sv_ind;

    _pLatConnect=NULL;
    _nLatConnect=0;

    DynamicSpecifyConnections(0); // purely lateral flow, no vertical

    // check for valid SV index
    ExitGracefullyIf(sv_ind == DOESNT_EXIST, "CmvLatRedistribute::unrecognized state variable specified in :SnowRedistribute command", BAD_DATA_WARN);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvLatRedistribute::~CmvLatRedistribute()
{
  for(int n=0; n<_nLatConnect;n++) { delete _pLatConnect[n]; }  delete[] _pLatConnect;
}

//////////////////////////////////////////////////////////////////
/// \brief Copies connnection array - called during RVH file read
//
void CmvLatRedistribute::AddConnectionArray(CLatConnect **pConnections, const int nConnections)
{
  _nLatConnect=nConnections;

  _pLatConnect=new CLatConnect *[_nLatConnect];
  for (int n=0; n<_nLatConnect;n++){
    _pLatConnect[n]=new CLatConnect(*pConnections[n]);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Initialization (prior to solution)
//
void CmvLatRedistribute::Initialize()
{

  if (_nLatConnect==0){
    ExitGracefully("CmvLatRedistribute::Initialize(): no lateral connections specified. Missing the :LateralConnections command?",BAD_DATA_WARN);
    return;
  }
  // Calculate the number of lateral connections
  _nLatConnections = _nLatConnect;

  // Allocate memory for the arrays
  _kFrom    = new int[_nLatConnections];
  _kTo      = new int[_nLatConnections];
  _iFromLat = new int[_nLatConnections];
  _iToLat   = new int[_nLatConnections];

  // Initialize the arrays using the methods from CLatConnect
  for (int q = 0; q < _nLatConnections; q++)
  {
    CLatConnect* connection = _pLatConnect[q];

    _kFrom   [q] = _pModel->GetHRUByID(_pLatConnect[q]->GetSourceHRUID()   )->GetGlobalIndex();
    _kTo     [q] = _pModel->GetHRUByID(_pLatConnect[q]->GetRecipientHRUID())->GetGlobalIndex();

    _iFromLat[q] = _iRedistributeFrom;
    _iToLat  [q] = _iRedistributeTo;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \param se_type [in] Model of soil evaporation used
/// \param *aSV [out] Array of state variable types needed by soil evaporation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by soil evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvLatRedistribute::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
    nSV = 0;
    // user specified 'from' & 'to' compartment, Levels - not known before construction
}

void CmvLatRedistribute::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
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
void CmvLatRedistribute::GetLateralExchange(const double *const *state_vars,
                                            const CHydroUnit *const *pHRUs,
                                            const optStruct &Options,
                                            const time_struct &tt,
                                            double *exchange_rates) const
{
  const double TAN_80_DEGREES = tan(80.0 * PI / 180.0);
  const double RAD_TO_DEG = 180.0 / PI;
  const double _k_param = 0.075; // Bernhardt & Schulz (2010) SnowSlide threshold method parameter

  // Snow density conversion factor
  const double SNOW_DENSITY = 250.0;     // kg/m³
  const double WATER_DENSITY = 1000.0;   // kg/m³
  const double SWE_TO_DEPTH_FACTOR = WATER_DENSITY / SNOW_DENSITY;  // = 4.0

  // Define minimum slope for lakes (30 degrees)
  const double LAKE_MIN_SLOPE_DEG = 30.0;
  const double LAKE_MIN_SLOPE_RAD = LAKE_MIN_SLOPE_DEG * PI / 180.0;

  for (int q = 0; q < _nLatConnections; q++)
  {
    double weight = _pLatConnect[q]->GetWeight();

    // Get the areas of both HRUs, slope and SWE from source HRU
    double areaFrom  = pHRUs[_kFrom[q]]->GetArea();
    double areaTo    = pHRUs[_kTo  [q]]->GetArea();
    double slope_rad = pHRUs[_kFrom[q]]->GetSlope();

    // Check if source HRU is a lake, if so, increase slope
    //HRU_type hru_type = pHRUs[_kFrom[q]]->GetHRUType();
    //bool is_lake = (hru_type == HRU_LAKE);

    // Override slope for lakes
    //if (is_lake) {
    //  slope_rad = LAKE_MIN_SLOPE_RAD;  // Use minimum slope for lakes
    //}

    double slope_deg = slope_rad * RAD_TO_DEG;
    double snowSWE = state_vars[_kFrom[q]][_iRedistributeFrom];
    double snowMoved = 0.0;

    if (_method == CONTINUOUS_REDIST)
    {
      // Original continuous method - no conversion to snow depth
      double snowTransport = min(0.5, slope_deg / TAN_80_DEGREES) * min(1.0, snowSWE / _max_snow_height);
      snowMoved = snowTransport * snowSWE;
    }
    else if (_method == THRESHOLD_REDIST)
    {
      // Bernhardt & Schulz (2010) SnowSlide threshold method
      const double SNOWSLIDE_LIMIT_ANGLE = 60.0; // Maximum angle in degrees

      // Convert SWE to snow depth (both in mm)
      double snowDepth = snowSWE * SWE_TO_DEPTH_FACTOR;

      // Original implementation (keeping as comment for reference)
      // double threshold_snow_height = _max_snow_height * exp(-_k_param * slope_deg);

      // FSM2 snow holding depth calculation
      double slope_thres = slope_deg;
      if (slope_thres < 10.0) {
          slope_thres = 10.0; // limit to >10 degrees to avoid infinite values
      }

      // Calculate threshold in meters (FSM2 formula works directly in meters)
      double threshold_snow_height_m = 3178.4 * pow(slope_thres, -1.998);

      // Convert from meters to mm for consistency with other calculations
      double threshold_snow_height = threshold_snow_height_m * MM_PER_METER;

      // Ensure no accumulation above SNOWSLIDE_LIMIT_ANGLE
      if (slope_deg > SNOWSLIDE_LIMIT_ANGLE) {
          threshold_snow_height = 0.0;
      }

      // Calculate excess snow to be redistributed
      if (snowDepth > threshold_snow_height) {
          double excessSnowDepth = snowDepth - threshold_snow_height;
          // Convert back to SWE
          snowMoved = excessSnowDepth / SWE_TO_DEPTH_FACTOR;
      }
      //std::cout << " slope=" << slope_deg << "° "
      //<< " threshold=" << threshold_snow_height << " mm"
      //<< " snowSWE=" << snowSWE << " mm"
      //<< " snowDepth=" << snowDepth << " mm"
      //<< " comparison: " << (snowDepth > threshold_snow_height ? "snowDepth > threshold" : "snowDepth <= threshold")
      //<< " excessDepth=" << (snowDepth > threshold_snow_height ? (snowDepth - threshold_snow_height) : 0) << " mm";
    }

    // Calculate the exchange rate considering area
    exchange_rates[q] = (snowMoved * weight * areaFrom) / Options.timestep; // [mm-km2/d]

    // Debug info - uncomment when needed

    // Debug info with comprehensive output
    //std::cout << "Connection " << q << ": fromHRU=" << fromHRU << " toHRU=" << toHRU
    //<< " method=" << (_method == CONTINUOUS_REDIST ? "CONTINUOUS" : "THRESHOLD");

    //std::cout << " snowMoved=" << snowMoved << " mm"
    //    << " exchange_rate=" << exchange_rates[q] << " mm-km²/d" << std::endl;
  }
}
