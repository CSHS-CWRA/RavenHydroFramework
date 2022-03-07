/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------*/
#include "ChannelXSect.h"
void TestManningsInfluence(const CChannelXSect *pChan,const double &Qref);
string FilenamePrepare(string filebase,const optStruct& Options); //Defined in StandardOutput.cpp
//////////////////////////////////////////////////////////////////
/// \brief Utility method to assign parameter name to cross section "nickname", add to static array of all x-sections
/// \param name [in] Nickname for cross section
//
void CChannelXSect::Construct(const string name)
{
  _name=name;
  _nSurveyPts=0;
  _aX   =NULL;
  _aElev=NULL;
  _aMann=NULL;
  _is_closed_channel=false;
  
  if (!DynArrayAppend((void**&)(pAllChannelXSects),(void*)(this),NumChannelXSects)){
    ExitGracefully("CChannelXSect::Constructor: creating NULL channel profile",BAD_DATA_WARN);
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
  _nSurveyPts=NumSurveyPts;
  _aX   =new double [_nSurveyPts];
  _aElev=new double [_nSurveyPts];
  _aMann=new double [_nSurveyPts];
  _min_mannings=ALMOST_INF;
  for (i=0;i<_nSurveyPts;i++)
  {
    _aX    [i]=X        [i];
    _aElev [i]=Elev     [i];
    _aMann[i]=ManningsN[i];
    lowerswap(_min_mannings,_aMann[i]);
    //check for valid mannings n
    if(_aMann[i]<0){
      ExitGracefully("CChannelXSect::Constructor: invalid mannings n",BAD_DATA);
    }
  }
  _min_stage_elev = _aElev[0];
  for (i=1;i<_nSurveyPts;i++){
    _min_stage_elev=min(_min_stage_elev,_aElev[i]);
  }
  _nPoints=20; //Default
  _aQ        =NULL;
  _aStage    =NULL;
  _aTopWidth =NULL;
  _aXArea    =NULL;
  _aPerim    =NULL;

  _bedslope=slope;
  ExitGracefullyIf(_bedslope<0.0,"CChannelXSect Constructor: channel profile bedslope must be greater than zero",BAD_DATA_WARN);

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
  
  _nPoints  =array_size;
  _aQ       =new double [_nPoints];
  _aStage   =new double [_nPoints];
  _aTopWidth=new double [_nPoints];
  _aXArea   =new double [_nPoints];
  _aPerim   =new double [_nPoints];
  for (int i=0;i<_nPoints;i++)
  {
    _aQ       [i]=flow [i];
    _aStage   [i]=stage[i];
    _aTopWidth[i]=width[i];
    _aXArea   [i]=area [i];
    _aPerim   [i]=perim[i];
  }
  _min_stage_elev = _aStage[0];
  for (int i=1;i<_nPoints;i++){
    _min_stage_elev=min(_min_stage_elev,_aStage[i]);
  }
  _bedslope=slope;
  _min_mannings=0.01;
  ExitGracefullyIf(_bedslope<=0.0,
                   "CChannelXSect Constructor: channel profile bedslope must be greater than zero",BAD_DATA_WARN);
}
//////////////////////////////////////////////////////////////////
/// \brief Constructor implementation if channel is a simple trapezoid
///
/// \param name [in] Nickname for cross section
/// \param bottom_w [in] Bottom wall length [m]
/// \param sidewall_ratio [in] Trapezoid wall ratio (0 for vertical walls, 1 for 45 deg angle, etc.)
/// \param bottom_elev [in] Elevation of riverbed [m]
/// \param mannings_n [in] Mannings roughness
/// \param slope [in] Riverbed slope
//
CChannelXSect::CChannelXSect(const string  name,          //constructor for trapezoid
                             const double  bottom_w,
                             const double  sidewall_ratio,
                             const double  bottom_elev,
                             const double  mannings_n,
                             const double  slope)
{
  Construct(name);
  _bedslope=slope;
  _min_mannings=mannings_n;
  ExitGracefullyIf(_bedslope<=0.0,
                  "CChannelXSect Constructor: channel profile bedslope must be greater than zero",BAD_DATA_WARN);
  ExitGracefullyIf(_min_mannings<=0.0,
                  "CChannelXSect Constructor: Manning's n must be greater than zero",BAD_DATA_WARN);

  _nPoints  =50;
  _aQ       =new double [_nPoints];
  _aStage   =new double [_nPoints];
  _aTopWidth=new double [_nPoints];
  _aXArea   =new double [_nPoints];
  _aPerim   =new double [_nPoints];
  double dh=0.1;//m -to 5 m
  double h;
  for (int i=0;i<_nPoints;i++)
  {
    h=i*dh;
    if (i==_nPoints-1){h=20;}//m
    _aStage   [i]=bottom_elev+h;
    _aTopWidth[i]=bottom_w+2.0*sidewall_ratio*h;
    _aXArea   [i]=(bottom_w+sidewall_ratio)*h;
    _aPerim   [i]=bottom_w+2*sqrt(h*h*(1+sidewall_ratio*sidewall_ratio));
    _aQ       [i]=sqrt(_bedslope)*_aXArea[i]*pow(_aXArea[i]/_aPerim[i],2.0/3.0)/_min_mannings;
  }
  _min_stage_elev = _aStage[0];
}
//////////////////////////////////////////////////////////////////
/// \brief Constructor implementation if channel is a circular pipe
///
/// \param name [in] Nickname for cross section
/// \param diameter [in] diameter [m]
/// \param bottom_elev [in] Elevation of riverbed [m]
/// \param mannings_n [in] Mannings roughness
/// \param slope [in] Riverbed slope
//
CChannelXSect::CChannelXSect( const string  name,         //constructor for pipe (of course, not really a channel)
                              const double  diameter,
                              const double  bottom_elev,
                              const double  mannings_n,
                              const double  bedslope) {
  Construct(name);
  _bedslope=bedslope;
  _min_mannings=mannings_n;
  ExitGracefullyIf(_bedslope<=0.0,
                  "CChannelXSect Constructor: channel profile bedslope must be greater than zero",BAD_DATA_WARN);
  ExitGracefullyIf(_min_mannings<=0.0,
                  "CChannelXSect Constructor: Manning's n must be greater than zero",BAD_DATA_WARN);

  _nPoints  =30;
  _aQ       =new double [_nPoints];
  _aStage   =new double [_nPoints];
  _aTopWidth=new double [_nPoints];
  _aXArea   =new double [_nPoints];
  _aPerim   =new double [_nPoints];
  double dh=diameter/(_nPoints-1);
  double h,theta;
  double r=diameter/2;
  for (int i=0;i<_nPoints;i++)
  {
    h=i*dh;
    theta=2.0*acos((r-h)/r);
    _aStage   [i]=bottom_elev+h;
    _aTopWidth[i]=sqrt((r-h/2.0)*h*8.0);
    _aXArea   [i]=r*r*(theta-sin(theta))/2.0;
    _aPerim   [i]=r*theta;
    _aQ       [i]=sqrt(_bedslope)*_aXArea[i]*pow(_aXArea[i]/_aPerim[i],2.0/3.0)/_min_mannings;
    if (_aPerim[i]==0){_aQ[i]=0.0;}
    /*cout<<"Conduit: "<<_aStage[i]<<" "<<_aTopWidth[i]<<" "<<_aXArea[i]<<" "<<_aPerim[i]<<" "<<_aQ[i]<<endl;
    if ((i > 0) && (_aQ[i] < _aQ[i-1])) {
      ExitGracefully("Circular Conduit constructor: decreasing flow rate",RUNTIME_ERR);
    }*/ //TO BE EXPECTED - BUT CAUSES CELERITY CALCULATION ISSUES
  }
  _min_stage_elev = bottom_elev;
  _is_closed_channel=true;
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CChannelXSect::~CChannelXSect()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING CHANNEL CROSS-SECTION "<<endl;}
  delete [] _aX;
  delete [] _aElev;
  delete [] _aMann;

  delete [] _aQ;
  delete [] _aStage;
  delete [] _aTopWidth;
  delete [] _aXArea;
  delete [] _aPerim;
}

/*****************************************************************
   Accessors
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns Cross-section nickname
/// \return Cross section nickname
//
string      CChannelXSect::GetName()             const { return _name;}

//////////////////////////////////////////////////////////////////
/// \brief Returns minimum mannings n along channel
/// \return minimum mannings n along channel
//
double      CChannelXSect::GetMinMannings()      const { return _min_mannings; }

//////////////////////////////////////////////////////////////////
/// \brief Returns Riverbed slope
/// \return Riverbed slope [m/m]
//
double      CChannelXSect::GetBedslope()         const { return _bedslope;}

/*****************************************************************
   Rating Curve Interpolation functions
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns Top width of channel [m]
/// \param &Q [in] Profile flowrate [m3/s]
/// \return Width of channel at input flowrate
//
double  CChannelXSect::GetTopWidth (const double &Q,const double &SB_slope,const double &SB_n) const//Q in m3/s
{
  double junk,Q_mult;
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);
  return InterpolateCurve(Q/Q_mult,_aQ,_aTopWidth,_nPoints,true);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns cross-sectional area of channel [m]
/// \param &Q [in] Profile area [m3/s]
/// \return Cross-sectional area of channel at input flowrate
//
double  CChannelXSect::GetArea     (const double &Q,const double &SB_slope,const double &SB_n) const//Q in m3/s
{
  double junk,Q_mult;
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);
  return InterpolateCurve(Q/Q_mult,_aQ,_aXArea,_nPoints,true);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns stage elevation at a point [m]
/// \param &Q [in] Profile flowrate [m3/s]
/// \return Stage elevation at specified flowrate
//
double  CChannelXSect::GetStageElev(const double &Q,const double &SB_slope,const double &SB_n) const//Q in m3/s
{
  double junk,Q_mult;
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);
  return InterpolateCurve(Q/Q_mult,_aQ,_aStage,_nPoints,true);
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
  double junk,Q_mult;
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);  
  return InterpolateCurve(Q/Q_mult,_aQ,_aPerim,_nPoints,true);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns correction terms for subbasin-specific slope and manning's n
/// \param &SB_slope - subbasin local slope (or AUTO_COMPUTE if channel slope should be used) [-]
/// \param &SB_n - subbasin local Manning's n (or AUTO_COMPUTE if channel n should be used) [-]
/// \return &slope_mult - correction multiplier for slope (S_actual = mult * S_channel) [-]
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
/// \brief Returns celerity of channel at reference flow [m/s]
/// \param &Q [in]  flowrate [m3/s]
/// \return Celerity of channel corresponding to flowrate [m/s]
//
double  CChannelXSect::GetCelerity(const double &Q, const double &SB_slope,const double &SB_n) const
{
  //returns dQ/dA|_Qref ~ wave celerity in reach at given flow Q
  //interpolated from rating curves

  ExitGracefullyIf(_aQ==NULL,"CChannelXSect::GetCelerity: Rating curves not yet generated",RUNTIME_ERR);

  double junk,Q_mult;
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);
  if ((_is_closed_channel) && (Q > Q_mult*_aQ[_nPoints-1])) {
    ExitGracefully("CChannelXSect::GetCelerity: reference flowrate exceeds closed channel maximum. Conduit area is too small",BAD_DATA);
  }
  if (Q<0){return 0.0;}
  int i=0; while ((Q>Q_mult*_aQ[i+1]) && (i<(_nPoints-2))){i++;}
  return Q_mult*(_aQ[i+1]-_aQ[i])/(_aXArea[i+1]-_aXArea[i]);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns diffusivity of channel [m/s]
/// \param &Q [in] Channel flowrate [m3/s]
/// \param &bedslope [in] Slope of riverbed [m/m]
/// \return Diffusivity of channel at specified flowrate and slope [m3/s]
//
double CChannelXSect::GetDiffusivity(const double &Q, const double &SB_slope, const double &SB_n) const
{
  ExitGracefullyIf(Q<=0,"CChannelXSect::GetDiffusivity: Invalid channel flowrate",BAD_DATA);
 
  ///< diffusivity from Roberson et al. 1995, Hydraulic Engineering \cite Roberson1998
  double slope_mult=1.0;
  double Q_mult    =1.0;
  GetFlowCorrections(SB_slope,SB_n,slope_mult,Q_mult);
  return 0.5*(Q)/GetTopWidth(Q,SB_slope,SB_n)/(slope_mult*_bedslope);
}

//////////////////////////////////////////////////////////////////
/// \brief checks to see if reference flow is in excess of channel maximum; throws warning if it is
/// \param &Qref [in] Reference flowrate [m3/s]
/// \param SB_slope [in] subbasin slope (or AUTO_COMPUTE, if channel slope to be used)
/// \param SB_n [in] subbasin mannings  (or AUTO_COMPUTE, if channel mannings to be used)
/// \param SBID [in] subbasin identifier
//
void  CChannelXSect::CheckReferenceFlow(const double& Qref,const double& SB_slope,const double& SB_n,const long SBID) const
{
  string warn;
  double junk,Q_mult;
  GetFlowCorrections(SB_slope,SB_n,junk,Q_mult);
  if((Qref > Q_mult*_aQ[_nPoints-1])) {
    warn="CChannelXSect::CheckReferenceFlow: reference flow exceeds channel maximum flow in subbasin "+to_string(SBID)+". Need to specify larger range of stages in channel profile.";
    WriteWarning(warn.c_str(),BAD_DATA);
  }
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
  for (int i=0;i<_nSurveyPts-1;i++)
  {
    zl=min(_aElev[i],_aElev[i+1]);
    zu=max(_aElev[i],_aElev[i+1]);
    dx=fabs(_aX[i+1]-_aX[i]);
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
      if ((i>0) && (_aElev[i-1]>_aElev[i]) && (_aX[i-1]==_aX[i])) //handles straight adjacent sides (left)
      {
        if(elev<=_aElev[i-1]){Pi+=(elev      -_aElev[i]);}
        else                {Pi+=(_aElev[i-1]-_aElev[i]);}
      }
      if ((i<(_nSurveyPts-2)) && (_aElev[i+2]>_aElev[i+1]) && (_aX[i+2]==_aX[i+1])) //handles straight adjacent sides (right)
      {
        if(elev<=_aElev[i+2]){Pi+=(elev      -_aElev[i+1]);}
        else                {Pi+=(_aElev[i+2]-_aElev[i+1]);}
      }
      if (i==0)             {Pi+=(elev-_aElev[i]  );}
      if (i==_nSurveyPts-2)  {Pi+=(elev-_aElev[i+1]);}
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
  _aPerim    =NULL;
  _aQ        =new double [_nPoints];
  _aStage    =new double [_nPoints];
  _aTopWidth =new double [_nPoints];
  _aXArea    =new double [_nPoints];
  _aPerim    =new double [_nPoints];
  ExitGracefullyIf(_aPerim==NULL,
                   "GenerateRatingCurvesFromProfile",OUT_OF_MEMORY);
  ExitGracefullyIf(_aElev==NULL,
                   "GenerateRatingCurvesFromProfile: bad profile array",BAD_DATA);

  double maxe=_aElev[0];
  double mine=_aElev[0];
  for (int i=0;i<_nSurveyPts;i++)
  {
    upperswap(maxe,_aElev[i]);
    lowerswap(mine,_aElev[i]);
  }
  ExitGracefullyIf((maxe-mine)<REAL_SMALL,
                   "CChannelXSect::GenerateRatingCurvesFromProfile: profile survey points all have same elevation",BAD_DATA);

  double dz=(maxe-mine)/((int)(_nPoints)-1.0);

  int i=0;
  for (double z=mine;z<maxe+0.5*dz;z+=dz)
  {
    _aStage[i]=z;
    GetPropsFromProfile(z,_aQ[i],_aTopWidth[i],_aXArea[i],_aPerim[i]);
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
/// \brief Check for duplicate channel names
//
void CChannelXSect::CheckForDuplicates(const optStruct &Options)
{
  for(int p=0; p<NumChannelXSects;p++)
  {
    for(int pp=0; pp<p;pp++)
    {
      if(pAllChannelXSects[p]->GetName()==pAllChannelXSects[pp]->GetName()) {
        string warn=" CChannelXSect::CheckForDuplicates: found duplicated channel name: "+pAllChannelXSects[p]->GetName();
        WriteWarning(warn.c_str(),Options.noisy);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Write rating curves to file rating_curves.csv
//
void CChannelXSect::WriteRatingCurves(const optStruct& Options)
{
  ofstream CURVES;
  string tmpFilename=FilenamePrepare("rating_curves.csv",Options);
  CURVES.open(tmpFilename.c_str());


  if (CURVES.fail()){
    ExitGracefully("CChannelXSect::WriteRatingCurves: Unable to open output file rating_curves.csv for writing.",FILE_OPEN_ERR);
  }
  int i;
  for (int p=0; p<NumChannelXSects;p++)
  {
    const CChannelXSect *pP=pAllChannelXSects[p];
    CURVES<<pP->_name <<"----------------"<<endl;
    CURVES<<"Flow Rate [m3/s],";    for(i=0;i<pP->_nPoints;i++){CURVES<<pP->_aQ[i]       <<",";}CURVES<<endl;
    CURVES<<"Stage Height [m],";    for(i=0;i<pP->_nPoints;i++){CURVES<<pP->_aStage[i]   <<",";}CURVES<<endl;
    CURVES<<"Top Width [m],";       for(i=0;i<pP->_nPoints;i++){CURVES<<pP->_aTopWidth[i]<<",";}CURVES<<endl;
    CURVES<<"X-sect area [m2],";    for(i=0;i<pP->_nPoints;i++){CURVES<<pP->_aXArea[i]   <<",";}CURVES<<endl;
    CURVES<<"Wetted Perimeter [m],";for(i=0;i<pP->_nPoints;i++){CURVES<<pP->_aPerim[i]   <<",";}CURVES<<endl;
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

