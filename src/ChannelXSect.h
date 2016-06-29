/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
------------------------------------------------------------------
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
    string       tag;                 /// <nickname for XSect
    
    double       bedslope;            /// <slope of river channel
    //double     conductivity;        /// <conductivity of channel bottom sediments [m/d]

    int          nSurveyPts;          /// <# of survey points
                                      /// <nSurveyPts=0 if rating curve not generated from surveyed data
    double      *aX;                  /// <survey points in local coordinates [m]
    double      *aElev;               /// <bottom channel elevation at survey pts (size nSurveyPts)[m above arbitrary datum]
    double      *aMann;               /// <Mannings roughness for each segment (between aX[i] and aX[i+1]) (size: nSurveyPts-1)
    
    int          N;                   /// <number of points on rating curves
    double      *aQ;                  /// <Rating curve for flow rates [m3/s]
    double      *aStage;              /// <Rating curve for stage elevation [m]
    double      *aTopWidth;           /// <Rating curve for top width [m]
    double      *aXArea;              /// <Rating curve for X-sectional area [m2]
    double      *aPerim;              /// <Rating curve for wetted perimeter [m]

    static CChannelXSect **pAllChannelXSects;
    static int             NumChannelXSects;

    void Construct                        (const string name);
    void GenerateRatingCurvesFromProfile  ();
    void GenerateRatingCurvesFromPowerLaw (); 

    void Interpolate        (const double &Q, 
                                   double &interp, int &i) const;
    void GetPropsFromProfile(const double &elev, 
                                   double &Q,   double &width, 
                                   double &A,   double &P);

  public:/*-------------------------------------------------------*/
    //Constructors:
    CChannelXSect( const string  name, 
                   const int     NumSurveyPts, 
                   const double *X, 
                   const double *Elev,
                   const double *ManningsN,
                   const double  bedslope);
    CChannelXSect( const string  name, 
                   const int     array_size, 
                   const double *flow,
                   const double *stage, 
                   const double *width,
                   const double *area,
                   const double *perim,
                   const double  bedslope);
    CChannelXSect (const string  name,          //constructor for power law
                   const double  a1, 
                   const double  b1, 
                   const double  a2,
                   const double  b2,
                   const double  a3,
                   const double  b3);
    CChannelXSect( const string  name,          //constructor for trapezoid
                   const double  bottom_w, 
                   const double  sidewall_angle, 
                   const double  bottom_elev,
                   const double  mannings_n,
                   const double  bedslope);
    ~CChannelXSect();

    //Accessors
    string              GetTag        ()                const;
    double              GetBedslope   ()                const;
    double              GetTopWidth   (const double &Q) const;
    double              GetArea       (const double &Q) const;
    double              GetStageElev  (const double &Q) const;
    double              GetWettedPerim(const double &Q) const;
    double              GetDepth      (const double &Q) const;
    double              GetCelerity   (const double &Qref) const;
    double              GetDiffusivity(const double &Q) const;

    //static accessors, destructor
    static int                 GetNumChannelXSects       ();
    static const CChannelXSect*StringToChannelXSect      (const string s);
    static void                DestroyAllChannelXSections();
    static void                SummarizeToScreen         ();
    static void                WriteRatingCurves         (); 
};
#endif