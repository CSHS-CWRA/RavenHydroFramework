/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#include "ChannelXSect.h"
void TestManningsInfluence(const CChannelXSect *pChan,const double &Qref);
//////////////////////////////////////////////////////////////////
/// \brief Utility method to assign parameter name to cross section "nickname", add to static array of all x-sections
/// \param name [in] Nickname for cross section
//
void CChannelXSect::Construct(const string name)
{
  _name=name;
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
  _aMann=new double [nSurveyPts];
  _min_mannings=ALMOST_INF;
  for (i=0;i<nSurveyPts;i++)
  {
    aX    [i]=X        [i];
    aElev [i]=Elev     [i];
    _aMann[i]=ManningsN[i];
    lowerswap(_min_mannings,_aMann[i]);
    //check for valid mannings n
    if(_aMann[i]<0){
      ExitGracefully("CChannelXSect::Constructor: invalid mannings n",BAD_DATA);
    }
  }
  _min_stage_elev = aElev[0];
  for (i=1;i<nSurveyPts;i++){
    _min_stage_elev=min(_min_stage_elev,aElev[i]);
  }
  N=20; //Default
  aQ        =NULL;
  aStage    =NULL;
  aTopWidth =NULL;
  aXArea    =NULL;
  aPerim    =NULL;

  _bedslope=slope;
  ExitGracefullyIf(_bedslope<0.0,"CChannelXSect Constructor: channel profile bedslope must be greater than zero",BAD_DATA);

  GenerateRatingCurvesFromProfile(); //All the work done here

 // TestManningsInfluence(this,20.0);
 // ExitGracefully("TestManningsInfluence Unit Testing",SIMULATION_DONE);
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
  _aMann=NULL;

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
  _min_stage_elev = aStage[0];
  for (int i=1;i<N;i++){
    _min_stage_elev=min(_min_stage_elev,aStage[i]);
  }
  _bedslope=slope;
  _min_mannings=0.01;
  ExitGracefullyIf(_bedslope<=0.0,
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
                             const double  slope,
                             const double  mannings,
                             const double  aS,const double  bS,
                             const double  aW,const double  bW,
                             const double  aP,const double  bP)
{
  Construct(name);  ExitGracefully("CChannelXSect::Constructor (Power law)",STUB);

  nSurveyPts=0;
  aX   =NULL;
  aElev=NULL;
  _aMann=NULL;

  N=32;
  aQ       =new double [N];
  aStage   =new double [N];
  aTopWidth=new double [N];
  aXArea   =new double [N];
  aPerim   =new double [N];
  for (int i=0;i<N;i++)
  {
    aQ       [i]=pow(10.0,(double)(i)/(double)(N)*4-1);//from 0 to 1000 m3/s
    aQ       [0]=0.0;
    aStage   [i]=aS*pow(aQ[i],bS);
    aTopWidth[i]=aW*pow(aQ[i],bW);
    aPerim   [i]=aP*pow(aQ[i],bP);
    aXArea   [i]=aStage   [i]*aTopWidth[i]; //approximate - may wish to revise
    
  }
  _min_stage_elev = aStage[0];
  for (int i=1;i<N;i++){
    _min_stage_elev=min(_min_stage_elev,aStage[i]);
  }
  _bedslope=slope;
  _min_mannings=mannings;
  ExitGracefullyIf(_bedslope<=0.0,
                   "CChannelXSect Constructor: channel profile bedslope must be greater than zero",BAD_DATA);

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
  _bedslope=slope;
  ExitGracefullyIf(_bedslope<=0.0,
                  "CChannelXSect Constructor: channel profile bedslope must be greater than zero",BAD_DATA);
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
  delete [] _aMann;

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
string      CChannelXSect::GetName()      const {return _name;}

//////////////////////////////////////////////////////////////////
/// \brief Returns minimum mannings n along channel
/// \return minimum mannings n along channel
//
double      CChannelXSect::GetMinMannings()      const { return _min_mannings; }


//////////////////////////////////////////////////////////////////
/// \brief Returns Riverbed slope
/// \return Riverbed slope [m/m]
//
double      CChannelXSect::GetBedslope() const {return _bedslope;}

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
double  CChannelXSect::GetTopWidth (const double &Q,const double &SB_slope,const double &SB_n) const//Q in m3/s
{
  int i;
  double interp;
  double junk,Q_mult;
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);
  Interpolate(Q/Q_mult,interp,i);
  return aTopWidth[i]+interp*(aTopWidth[i+1]-aTopWidth[i]);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns cross-sectional area of channel [m]
/// \param &Q [in] Profile area [m3/s]
/// \return Cross-sectional area of channel at input flowrate
//
double  CChannelXSect::GetArea     (const double &Q,const double &SB_slope,const double &SB_n) const//Q in m3/s
{
  int i;
  double interp;
  double junk,Q_mult;
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);
  Interpolate(Q/Q_mult,interp,i);
  return aXArea[i]+interp*(aXArea[i+1]-aXArea[i]);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns stage elevation at a point [m]
/// \param &Q [in] Profile flowrate [m3/s]
/// \return Stage elevation at specified flowrate
//
double  CChannelXSect::GetStageElev(const double &Q,const double &SB_slope,const double &SB_n) const//Q in m3/s
{
  int i;
  double interp;
  double junk,Q_mult;
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);
  Interpolate(Q/Q_mult,interp,i);
  return aStage[i]+interp*(aStage[i+1]-aStage[i]);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns channel depth at a point [m]
/// \details subtracts stage elevation at the specified discharge by the lowest elevation given
/// \param &Q [in] Profile flowrate [m3/s]
/// \return channel depth at specified flowrate
//
double  CChannelXSect::GetDepth(const double &Q,const double &SB_slope,const double &SB_n) const//Q in m3/s
{
  return (GetStageElev(Q,SB_slope,SB_n)-_min_stage_elev);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns Wetted perimeter of channel [m]
/// \param &Q [in] Profile flowrate [m3/s]
/// \return Wetted perimeter of channel at specified flowrate
//
double  CChannelXSect::GetWettedPerim(const double &Q,const double &SB_slope,const double &SB_n) const//Q in m3/s
{
  int i;
  double interp;
  double junk,Q_mult;
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);
  Interpolate(Q/Q_mult,interp,i);
  return aPerim[i]+interp*(aPerim[i+1]-aPerim[i]);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns correction terms for subbasin-specific slope and manning's n
/// \param &SB_slope - subbasin local slope (or AUTO_COMPUTE if channel slope should be used) [-]
/// \param &SB_n - subbasin local Manning's n (or AUTO_COMPUTE if channel n should be used) [-]
/// \return &slope_mult - correction multiplier for slope (S_actual = mult * Q_channel) [-]
/// \return &Q_mult - correction multiplier for flow (Q_actual = mult * Q_channel) [-]
//
void CChannelXSect::GetFlowCorrections(const double &SB_slope,
                                       const double &SB_n,
                                       double &slope_mult,
                                       double &Q_mult) const
{
  slope_mult=1.0;
  Q_mult=1.0;
  if(SB_slope!=AUTO_COMPUTE){
    slope_mult=(SB_slope/_bedslope);
    Q_mult    =pow(SB_slope/_bedslope,0.5); //Mannings formula correction
  }
  if(SB_n    !=AUTO_COMPUTE){
    Q_mult    *=(_min_mannings/SB_n);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Returns celerity of channel [m/s]
/// \param &Qref [in] Reference flowrate [m3/s]
/// \return Celerity of channel corresponding to reference flowrate [m/s]
//
double  CChannelXSect::GetCelerity(const double &Qref, const double &SB_slope,const double &SB_n) const//Qref in m3/s, celerity in m/s
{
  //returns dQ/dA|_Qref ~ wave celerity in reach at given reference flow Qref
  //interpolated from rating curves

  ExitGracefullyIf(aQ==NULL,"CChannelXSect::GetCelerity: Rating curves not yet generated",RUNTIME_ERR);
  double Q_mult;
  double junk;//unused output
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);

  if (Qref/Q_mult<0){return 0.0;}
  int i=0; while ((Qref/Q_mult>aQ[i+1]) && (i<(N-2))){i++;}

  return (aQ[i+1]-aQ[i])/(aXArea[i+1]-aXArea[i]);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns diffusivity of channel [m/s]
/// \param &Q [in] Channel flowrate [m3/s]
/// \param &bedslope [in] Slope of riverbed [m/m]
/// \return Diffusivity of channel at specified flowrate and slope [m3/s]
//
double CChannelXSect::GetDiffusivity(const double &Q, const double &SB_slope, const double &SB_n) const
{
  ExitGracefullyIf(Q<=0,"CChannelXSect::Invalid channel flowrate",BAD_DATA);
  ///< diffusivity from Roberson et al. 1995, Hydraulic Engineering \cite Roberson1998
  double slope_mult=1.0;
  double Q_mult=1.0;
  GetFlowCorrections(SB_slope,SB_n,slope_mult,Q_mult);
  return 0.5*(Q)/GetTopWidth(Q,SB_slope,SB_n)/(slope_mult*_bedslope);
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
      Q+=pow(_bedslope,0.5)*Ai*pow(Ai/Pi,2.0/3.0)/_aMann[i];
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

  double dz=(maxe-mine)/((int)(N)-1.0);

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
    cout<<"-Channel profile \""<<pAllChannelXSects[p]->GetName()<<"\" "<<endl;
    cout<<"           slope: " <<pAllChannelXSects[p]->_bedslope <<endl;
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
    CURVES<<pP->_name <<"----------------"<<endl;
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
    if (!sup.compare(StringToUppercase(pAllChannelXSects[p]->GetName()))){return pAllChannelXSects[p];}
  }
  return NULL;
}


void TestManningsInfluence(const CChannelXSect *pChan, const double &Qref) {
  ofstream TEST;
  TEST.open("ManningsTest.csv");

  TEST<<pChan->GetName()<<endl;
  TEST<<"ManningsN,1/n,Diffusivity,Celerity,Qref*2 Diff, Qref*2 Cel, Qref*3 Diff, Qref*3 Cel"<<endl;
  for(double overn=5;overn<100;overn+=5.0) {
    double n=1.0/overn;
    TEST<<n<<",";
    TEST<<overn<<",";
    TEST<<pChan->GetDiffusivity(1*Qref,AUTO_COMPUTE,n)<<",";
    TEST<<pChan->GetCelerity   (1*Qref,AUTO_COMPUTE,n)<<",";
    TEST<<pChan->GetDiffusivity(2*Qref,AUTO_COMPUTE,n)<<",";
    TEST<<pChan->GetCelerity   (2*Qref,AUTO_COMPUTE,n)<<",";
    TEST<<pChan->GetDiffusivity(3*Qref,AUTO_COMPUTE,n)<<",";
    TEST<<pChan->GetCelerity   (3*Qref,AUTO_COMPUTE,n)<<",";
    TEST<<endl;
  }
  TEST.close();
}

