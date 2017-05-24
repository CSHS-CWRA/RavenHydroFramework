/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include "ChannelXSect.h"

//////////////////////////////////////////////////////////////////
/// \brief Utility method to assign parameter name to cross section "nickname", add to static array of all x-sections
/// \param name [in] Nickname for cross section
//
void CChannelXSect::Construct(const string name)
{
  tag=name;
  if (!DynArrayAppend((void**&)(pAllChannelXSects),(void*)(this),NumChannelXSects)){
    ExitGracefully("CChannelXSect::Constructor: creating NULL channel profile",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Constructor implementation if channel profile is specified from survey data
/// \param name [in] Nickname for cross section
/// \param NumSurveyPts [in] number of survey points taken in data
/// \param *X [in] Array of coordinates of survey points
/// \param *Elev [in] Array of elevations of riverbed at survey points
/// \param *ManningsN [in] Mannings roughness for each segment (effective size: NumSurveyPts-1)
/// \param slope [in] Riverbed slope
//
CChannelXSect::CChannelXSect(const string  name,
                             const int     NumSurveyPts,
                             const double *X,
                             const double *Elev,
                             const double *ManningsN,
                             const double  slope)
{
  Construct(name);
  int i;
  nSurveyPts=NumSurveyPts;
  aX   =new double [nSurveyPts];
  aElev=new double [nSurveyPts];
  aMann=new double [nSurveyPts];
  for (i=0;i<nSurveyPts;i++)
  {
    aX   [i]=X        [i];//cout<<aX[i]<<" ";
    aElev[i]=Elev     [i];//cout<<aElev[i]<<" ";
    aMann[i]=ManningsN[i];//cout<<aMann[i]<<endl;

    //check for valid mannings n
    if(aMann[i]<0){
      ExitGracefully("CChannelXSect::Constructor: invalid mannings n",BAD_DATA);
    }
  }

  N=20; //Default
  aQ        =NULL;
  aStage    =NULL;
  aTopWidth =NULL;
  aXArea    =NULL;
  aPerim    =NULL;

  bedslope=slope;
  if(bedslope<0.0){
    ExitGracefully("CChannelXSect Constructor: channel profile bedslope must be greater than zero",BAD_DATA);
  }

  GenerateRatingCurvesFromProfile(); //All the work done here
}
//////////////////////////////////////////////////////////////////
/// \brief Constructor implementation if channel profile is specified from rating curves
/// \param name [in] Nickname for cross section
/// \param array_size [in] number of survey points taken in data
/// \param *flow [in] Array of flow rates corresponding to stages, in cms [size:array_size]
/// \param *stage [in] Array of stage elevations, in meters [size:array_size]
/// \param *width [in] Array of channel top widths corresponding to stages, in m [size:array_size]
/// \param *area [in] Array of channel x-sectional areas corresponding to stages, in m2 [size:array_size]
/// \param *perim [in] Array of channel wetted perimeters corresponding to stages, in m [size:array_size]
/// \param slope [in] Riverbed slope
//
CChannelXSect::CChannelXSect(const string  name,
                             const int     array_size,
                             const double *flow,
                             const double *stage,
                             const double *width,
                             const double *area,
                             const double *perim,
                             const double  slope)
{
  Construct(name);
  nSurveyPts=0;
  aX   =NULL;
  aElev=NULL;
  aMann=NULL;

  N=array_size;
  aQ       =new double [N];
  aStage   =new double [N];
  aTopWidth=new double [N];
  aXArea   =new double [N];
  aPerim   =new double [N];
  for (int i=0;i<N;i++)
  {
    aQ       [i]=flow [i];
    aStage   [i]=stage[i];
    aTopWidth[i]=width[i];
    aXArea   [i]=area [i];
    aPerim   [i]=perim[i];
  }
  bedslope=slope;
  ExitGracefullyIf(bedslope<=0.0,
                   "CChannelXSect Constructor: channel profile bedslope must be greater than zero",BAD_DATA);
}
//////////////////////////////////////////////////////////////////
/// \brief Constructor implementation if routing curves are explicitly defined using power law coefficients
///
/// \param name [in] Nickname for cross section
/// \param aS [in] power law coefficient for stage S=aS*Q^bS
/// \param bS [in] power law coefficient for stage S=aS*Q^bS
/// \param aW [in] power law coefficient for top width W=aW*Q^bW
/// \param bW [in] power law coefficient for top width W=aW*Q^bW
/// \param aP [in] power law coefficient for wetted perimeter P=aP*Q^bP
/// \param bP [in] power law coefficient for wetted perimeter P=aP*Q^bP
//
CChannelXSect::CChannelXSect(const string  name,          //constructor for power law
                             const double  aS,const double  bS,
                             const double  aW,const double  bW,
                             const double  aP,const double  bP)
{
  Construct(name);
  /// \todo [add funct] code power law profile constructor
  ExitGracefully("CChannelXSect::Constructor (Power law)",STUB);
}

//////////////////////////////////////////////////////////////////
/// \brief Constructor implementation if channel is a simple trapezoid
///
/// \param name [in] Nickname for cross section
/// \param bottom_w [in] Bottom wall length [m]
/// \param sidewall_angle [in] Trapezoid wall angle [rad?]
/// \param bottom_elev [in] Elevation of riverbed [m]
/// \param mannings_n [in] Mannings roughness
/// \param slope [in] Riverbed slope
//
CChannelXSect::CChannelXSect(const string  name,          //constructor for trapezoid
                             const double  bottom_w,
                             const double  sidewall_angle,
                             const double  bottom_elev,
                             const double  mannings_n,
                             const double  slope)
{
  Construct(name);
  bedslope=slope;
  /// \todo [add funct] code trapezoidal profile constructor
  ExitGracefully("CChannelXSect::Constructor (Trapezoid)",STUB);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CChannelXSect::~CChannelXSect()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING CHANNEL CROSS-SECTION "<<endl;}
  delete [] aX;
  delete [] aElev;
  delete [] aMann;

  delete [] aQ;
  delete [] aStage;
  delete [] aTopWidth;
  delete [] aXArea;
  delete [] aPerim;
}

/*****************************************************************
   Accessors
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns Cross-section nickname
/// \return Cross section nickname
//
string      CChannelXSect::GetTag()      const {return tag;}

//////////////////////////////////////////////////////////////////
/// \brief Returns Riverbed slope
/// \return Riverbed slope [m/m]
//
double      CChannelXSect::GetBedslope() const {return bedslope;}

/*****************************************************************
   Rating Curve Interpolation functions
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Interpolate value between discrete points on rating curve
/// \todo [optimize] CChannelXSect::Interpolate: This needs optimization
///
/// \param &Q [in] Flow rate [m3/s]
/// \param &interp [out] Interpolation weight between aQ[i] and aQ[i+1]
/// \param &i [out] index of aQ array corresponding to value just less than Q
//
void CChannelXSect::Interpolate(const double &Q, double &interp, int &i) const
{
  //This needs OPTIMIZATION!
  static int ilast;
  if (Q<0){interp= 0.0;i=0;ilast=i;return;}
  ExitGracefullyIf(aQ==NULL,
                   "CChannelXSect::Interpolate: Rating curves not yet generated",RUNTIME_ERR);

  //standard version-unoptimized
  //should start out assuming current Q close to or same as old Q
  i=0; while ((Q>aQ[i+1]) && (i<(N-2))){i++;}//Dumb Search
  //SmartIntervalSearch(Q,aQ,N,i,ilast);
  interp=(Q-aQ[i])/(aQ[i+1]-aQ[i]); //1<, >0, unless i=N-1, then >1
  ilast=i;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns Top width of channel [m]
/// \param &Q [in] Profile flowrate [m3/s]
/// \return Width of channel at input flowrate
//
double  CChannelXSect::GetTopWidth (const double &Q) const//Q in m3/s
{
  int i;
  double interp;
  Interpolate(Q,interp,i);
  return aTopWidth[i]+interp*(aTopWidth[i+1]-aTopWidth[i]);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns cross-sectional area of channel [m]
/// \param &Q [in] Profile area [m3/s]
/// \return Cross-sectional area of channel at input flowrate
//
double  CChannelXSect::GetArea     (const double &Q) const//Q in m3/s
{
  int i;
  double interp;
  Interpolate(Q,interp,i);
  return aXArea[i]+interp*(aXArea[i+1]-aXArea[i]);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns stage elevation at a point [m]
/// \param &Q [in] Profile flowrate [m3/s]
/// \return Stage elevation at specified flowrate
//
double  CChannelXSect::GetStageElev(const double &Q) const//Q in m3/s
{
  int i;
  double interp;
  Interpolate(Q,interp,i);
  return aStage[i]+interp*(aStage[i+1]-aStage[i]);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns channel depth at a point [m]
/// \details subtracts stage elevation at the specified discharge by the lowest elevation given
/// \param &Q [in] Profile flowrate [m3/s]
/// \return channel depth at specified flowrate
//
double  CChannelXSect::GetDepth(const double &Q) const//Q in m3/s
{
  double min_elev = aElev[0];
  for (int i=0;i<nSurveyPts;i++){
    min_elev=min(min_elev,aElev[i+1]);
  }
  return (GetStageElev(Q)-min_elev);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns Wetted perimeter of channel [m]
/// \param &Q [in] Profile flowrate [m3/s]
/// \return Wetted perimeter of channel at specified flowrate
//
double  CChannelXSect::GetWettedPerim(const double &Q) const//Q in m3/s
{
  int i;
  double interp;
  Interpolate(Q,interp,i);
  return aPerim[i]+interp*(aPerim[i+1]-aPerim[i]);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns celerity of channel [m/s]
/// \param &Qref [in] Reference flowrate [m3/s]
/// \return Celerity of channel corresponding to reference flowrate
//
double  CChannelXSect::GetCelerity(const double &Qref) const//Qref in m3/s, celerity in m/s
{
  //returns dQ/dA|_Qref ~ wave celerity in reach at given reference flow Qref
  //interpolated from rating curves

  ExitGracefullyIf(aQ==NULL,"CChannelXSect::GetCelerity: Rating curves not yet generated",RUNTIME_ERR);
  if (Qref<0){return 0.0;}
  int i=0; while ((Qref>aQ[i+1]) && (i<(N-2))){i++;}

  return (aQ[i+1]-aQ[i])/(aXArea[i+1]-aXArea[i]);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns diffusivity of channel [m/s]
/// \param &Q [in] Channel flowrate [m3/s]
/// \param &bedslope [in] Slope of riverbed [m/m]
/// \return Diffusivity of channel at specified flowrate and slope [m2/s]
//
double CChannelXSect::GetDiffusivity(const double &Q) const//Q in m2/s
{
  ExitGracefullyIf(Q<=0,"CChannelXSect::Invalid channel flowrate",BAD_DATA);
  ///< diffusivity from Roberson et al. 1995, Hydraulic Engineering \cite Roberson1998
  return 0.5*Q/GetDepth(Q)/bedslope;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns varous properties from profile, given "elev"
/// \details Given wetted stage 'elev', returns flowrate (Q), top width (w),
/// wetted perimeter (P), and x-sectional area (A).
/// assumes only that survey points are ordered
/// \author contributions from Susan Huang, Aug-2012
///
/// \param &elev [in] Profile elevation [m]
/// \param &Q [out] Flowrate [m3/s]
/// \param w [out] Top width [m]
/// \param &A [out] Cross sectional area [m2]
/// \param &P [out] Wetted perimeter [m]
//
void CChannelXSect::GetPropsFromProfile(const double &elev,
                                        double &Q,
                                        double &w, double &A, double &P)
{
  double zl,zu; //upper and lower bottom elevation for surveyed channel segment
  double dx;    //length of surveyed channel segment
  double Ai,Pi,wi;

  A=0;P=0;Q=0;w=0;
  for (int i=0;i<nSurveyPts-1;i++)
  {
    zl=min(aElev[i],aElev[i+1]);
    zu=max(aElev[i],aElev[i+1]);
    dx=fabs(aX[i+1]-aX[i]);
    if ((elev<=zl) || (dx==0))//dry part of reach, ignored
    {
      Ai=0;
      Pi=0;
      wi=0;
    }
    else if (elev>=zu) //fully wet part of profile
    {
      wi=dx;
      Ai=((elev-zu)*dx+0.5*dx*(zu-zl)); //trapezoidal section
      Pi=sqrt(dx*dx+(zu-zl)*(zu-zl));
      if ((i>0) && (aElev[i-1]>aElev[i]) && (aX[i-1]==aX[i])) //handles straight adjacent sides (left)
      {
        if(elev<=aElev[i-1]){Pi+=(elev      -aElev[i]);}
        else                {Pi+=(aElev[i-1]-aElev[i]);}
      }
      if ((i<(nSurveyPts-2)) && (aElev[i+2]>aElev[i+1]) && (aX[i+2]==aX[i+1])) //handles straight adjacent sides (right)
      {
        if(elev<=aElev[i+2]){Pi+=(elev      -aElev[i+1]);}
        else                {Pi+=(aElev[i+2]-aElev[i+1]);}
      }
      if (i==0)             {Pi+=(elev-aElev[i]  );}
      if (i==nSurveyPts-2)  {Pi+=(elev-aElev[i+1]);}
    }
    else  //partially wet part of profile (includes riverbank)
    {
      double ddx=(elev-zl)/(zu-zl)*dx; //width of wetted portion
      if ((zu-zl)<REAL_SMALL) {ddx=0.0;}//essentially flat reach, elev=zu=zl
      Ai=0.5*(elev-zl)*ddx; //triangular section
      Pi=sqrt(ddx*ddx+(elev-zl)*(elev-zl));
      wi=ddx;
    }
    A+=Ai;
    P+=Pi;
    w+=wi;
    //using Mannings equation to calculate Q
    if (Ai>0){
      Q+=pow(bedslope,0.5)*Ai*pow(Ai/Pi,2.0/3.0)/aMann[i];
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Generates rating curves from profile (utility method)
/// \details Generates the rating curves for Stage, Top width, area, and wetted
/// perimeter at a discrete set of points
//
void CChannelXSect::GenerateRatingCurvesFromProfile()
{
  aPerim    =NULL;
  aQ        =new double [N];
  aStage    =new double [N];
  aTopWidth =new double [N];
  aXArea    =new double [N];
  aPerim    =new double [N];
  ExitGracefullyIf(aPerim==NULL,
                   "GenerateRatingCurvesFromProfile",OUT_OF_MEMORY);
  ExitGracefullyIf(aElev==NULL,
                   "GenerateRatingCurvesFromProfile: bad profile array",BAD_DATA);

  double maxe=aElev[0];
  double mine=aElev[0];
  for (int i=0;i<nSurveyPts;i++)
  {
    upperswap(maxe,aElev[i]);
    lowerswap(mine,aElev[i]);
  }
  ExitGracefullyIf((maxe-mine)<REAL_SMALL,
                   "CChannelXSect::GenerateRatingCurvesFromProfile: profile survey points all have same elevation",BAD_DATA);

  double dz=(maxe-mine)/(N-1);

  int i=0;
  for (double z=mine;z<maxe+0.5*dz;z+=dz)
  {
    aStage[i]=z;
    GetPropsFromProfile(z,aQ[i],aTopWidth[i],aXArea[i],aPerim[i]);
    i++;
  }
}
/*****************************************************************
   Static Initialization, Accessors, Destructors
*****************************************************************/
CChannelXSect **CChannelXSect::pAllChannelXSects=NULL;
int             CChannelXSect::NumChannelXSects =0;

//////////////////////////////////////////////////////////////////
/// \brief Returns number of channel cross-sections in model
/// \return Number of channel cross-sections in model
//
int   CChannelXSect::GetNumChannelXSects       (){return NumChannelXSects;}

//////////////////////////////////////////////////////////////////
/// \brief Deletes all channel profiles in model
//
void  CChannelXSect::DestroyAllChannelXSections()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL CHANNEL PROFILES"<<endl;}
  for (int p=0; p<NumChannelXSects;p++){
    delete pAllChannelXSects[p];
  }
  delete [] pAllChannelXSects;
}

//////////////////////////////////////////////////////////////////
/// \brief Summarize profile information to screen
//
void CChannelXSect::SummarizeToScreen         ()
{
  cout<<"==================="<<endl;
  cout<<"Channel Profile Summary:"<<NumChannelXSects<<" profiles in database"<<endl;
  for (int p=0; p<NumChannelXSects;p++)
  {
    cout<<"-Channel profile \""<<pAllChannelXSects[p]->GetTag()<<"\" "<<endl;
    cout<<"           slope: " <<pAllChannelXSects[p]->bedslope <<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Write rating curves to file rating_curves.csv
//
void CChannelXSect::WriteRatingCurves()
{
  ofstream CURVES;
  CURVES.open("rating_curves.csv");
  if (CURVES.fail()){
    ExitGracefully("CChannelXSect::WriteRatingCurves: Unable to open output file rating_curves.csv for writing.",FILE_OPEN_ERR);
  }
  int i;
  for (int p=0; p<NumChannelXSects;p++)
  {
    const CChannelXSect *pP=pAllChannelXSects[p];
    CURVES<<pP->tag <<"----------------"<<endl;
    CURVES<<"Flow Rate [m3/s],";    for(i=0;i<pP->N;i++){CURVES<<pP->aQ[i]       <<",";}CURVES<<endl;
    CURVES<<"Stage Height [m],";    for(i=0;i<pP->N;i++){CURVES<<pP->aStage[i]   <<",";}CURVES<<endl;
    CURVES<<"Top Width [m],";       for(i=0;i<pP->N;i++){CURVES<<pP->aTopWidth[i]<<",";}CURVES<<endl;
    CURVES<<"X-sect area [m2],";    for(i=0;i<pP->N;i++){CURVES<<pP->aXArea[i]   <<",";}CURVES<<endl;
    CURVES<<"Wetted Perimeter [m],";for(i=0;i<pP->N;i++){CURVES<<pP->aPerim[i]   <<",";}CURVES<<endl;
  }
  CURVES.close();
}

//////////////////////////////////////////////////////////////////
/// \brief Converts string (e.g., "X2305" in Basin file) to channel profile
/// \param s [in] String to be converted to a channel profile
/// \return Pointer to channel profile with name equivalent to string, or NULL if string is invalid
//
const CChannelXSect*CChannelXSect::StringToChannelXSect(const string s)
{
  string sup=StringToUppercase(s);
  for (int p=0;p<NumChannelXSects;p++)
  {
    if (!sup.compare(StringToUppercase(pAllChannelXSects[p]->GetTag()))){return pAllChannelXSects[p];}
    //else if (s_to_i(s.c_str())==(p+1))         {return pAllChannelXSects[p];}
  }
  return NULL;
}




