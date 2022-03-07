/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------
  ChannelXSect.h
  ------------------------------------------------------------------
  defines Channel Cross section & rating curve calculations
  ----------------------------------------------------------------*/
#ifndef CHANNELX_H
#define CHANNELX_H

#include "RavenInclude.h"

/*****************************************************************
   Class CChannelXSect
------------------------------------------------------------------
   Data Abstraction for river/stream channel cross section
   with associated rating curve
******************************************************************/
class CChannelXSect
{
private:/*-------------------------------------------------------*/
  string       _name;               /// <nickname for XSect

  double       _bedslope;           /// <slope of river channel
  //double     _conductivity;       /// <conductivity of channel bottom sediments [m/d]
  bool         _is_closed_channel;  /// <true for closed channels/conduits 

  int          _nSurveyPts;         /// <# of survey points
  /// <nSurveyPts=0 if rating curve not generated from surveyed data
  double      *_aX;                 /// <survey points in local coordinates [m]
  double      *_aElev;              /// <bottom channel elevation at survey pts (size nSurveyPts)[m above arbitrary datum]
  double      *_aMann;              /// <Mannings roughness for each segment (between aX[i] and aX[i+1]) (size: nSurveyPts-1)
  double       _min_mannings;       /// <minimum roughness coefficient for all segments
  double       _min_stage_elev;     /// <minimum elevation of channel profile

  int          _nPoints;            /// <number of points on rating curves
  double      *_aQ;                 /// <Rating curve for flow rates [m3/s]
  double      *_aStage;             /// <Rating curve for stage elevation [m]
  double      *_aTopWidth;          /// <Rating curve for top width [m]
  double      *_aXArea;             /// <Rating curve for X-sectional area [m2]
  double      *_aPerim;             /// <Rating curve for wetted perimeter [m]

  static CChannelXSect **pAllChannelXSects;
  static int             NumChannelXSects;

  void Construct                        (const string name);
  void GenerateRatingCurvesFromProfile  ();
  void GenerateRatingCurvesFromPowerLaw ();

  //void Interpolate        (const double &Q,double &interp, int &i) const;
  void GetPropsFromProfile(const double &elev,
                           double &Q,   double &width,
                           double &A,   double &P);
  void GetFlowCorrections(const double &SB_slope,const double &SB_n,
                          double &slope_mult,double &Q_mult) const;

public:/*-------------------------------------------------------*/
  //Constructors:
  CChannelXSect( const string  name,          //constructor from profile points 
                 const int     NumSurveyPts,
                 const double *X,
                 const double *Elev,
                 const double *ManningsN,
                 const double  bedslope);
  CChannelXSect( const string  name,          //constructor from raw data
                 const int     array_size,
                 const double *flow,
                 const double *stage,
                 const double *width,
                 const double *area,
                 const double *perim,
                 const double  bedslope);
  CChannelXSect( const string  name,          //constructor for trapezoid
                 const double  bottom_w,
                 const double  sidewall_angle,
                 const double  bottom_elev,
                 const double  mannings_n,
                 const double  bedslope);
  CChannelXSect( const string  name,         //constructor for pipe (of course, not really a channel)
                 const double  diameter,
                 const double  bottom_elev,
                 const double  mannings_n,
                 const double  bedslope);
  ~CChannelXSect();

  //Accessors
  string              GetName       ()                const;
  double              GetBedslope   ()                const;
  double              GetMinMannings()                const;
  double              GetTopWidth   (const double &Q, const double &SB_slope,const double &SB_n) const;
  double              GetArea       (const double &Q, const double &SB_slope,const double &SB_n) const;
  double              GetStageElev  (const double &Q, const double &SB_slope,const double &SB_n) const;
  double              GetWettedPerim(const double &Q, const double &SB_slope,const double &SB_n) const;
  double              GetDepth      (const double &Q, const double &SB_slope,const double &SB_n) const;
  double              GetCelerity   (const double &Qref, const double &SB_slope,const double &SB_n) const;
  double              GetDiffusivity(const double &Q, const double &SB_slope, const double &SB_n) const;

  void                CheckReferenceFlow(const double& Qref,const double& SB_slope,const double& SB_n,const long SBID) const;

  //static accessors, destructor
  static int                 GetNumChannelXSects       ();
  static const CChannelXSect*StringToChannelXSect      (const string s);
  static void                CheckForDuplicates        (const optStruct& Options);
  static void                DestroyAllChannelXSections();
  static void                SummarizeToScreen         ();
  static void                WriteRatingCurves         (const optStruct& Options);
};
#endif
