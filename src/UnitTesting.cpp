/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "ParseLib.h"
#include "Radiation.h"
#include "GlobalParams.h"
#include "UnitTesting.h"
void AddTimeTest();
void RavenUnitTesting(const optStruct &Options)
{
  //cout<<"RAVEN UNIT TESTING MODE"<<endl;

  //uncomment one for use:
  //DateTest();
  //OpticalAirMassTest();
  //ClearSkyTest();
  //ShortwaveTest();
  //ShortwaveGenerator();
  //JulianConvertTest();
  //SmartLookupUnitTest();
  //SmartIntervalTest();
  //AddTimeTest();
}
/////////////////////////////////////////////////////////////////
/// \brief Tests DateStringToTimeStruct() function
//
void DateTest()
{
  time_struct tt;
  cout<<"2012-03-27 12:43:02.01"<<endl;
  tt=DateStringToTimeStruct("2012-03-27","12:43:02.01");
  cout<<tt.date_string<<" month, day, year: "<<tt.month <<","<<tt.day_of_month<<","<<tt.year<<" julian: "<<tt.julian_day<<endl;

  cout<<"2000-12-31 12:00:00"<<endl;
  tt=DateStringToTimeStruct("2000-12-31","12:00:00");
  cout<<tt.date_string<<" month, day, year: "<<tt.month <<","<<tt.day_of_month<<","<<tt.year<<" julian: "<<tt.julian_day<<endl;

  cout<<"1997/12/31 23:00:00.0000"<<endl;
  tt=DateStringToTimeStruct("1997/12/31","23:00:00.0000");
  cout<<tt.date_string<<" month, day, year: "<<tt.month <<","<<tt.day_of_month<<","<<tt.year<<" julian: "<<tt.julian_day<<endl;

  cout<<"0000/01/01 24:00:00.0000"<<endl;
  tt=DateStringToTimeStruct("0000/01/01","23:59:00.0000");
  cout<<tt.date_string<<" month, day, year: "<<tt.month <<","<<tt.day_of_month<<","<<tt.year<<" julian: "<<tt.julian_day<<endl;
  ExitGracefully("DateTest",SIMULATION_DONE);
}
void AddTimeTest() {
  time_struct tt,tt2;
  double outday;int outyear;
  AddTime(360,1999,5,outday,outyear);
  JulianConvert(0.0,360,1999,tt);
  JulianConvert(0.0,outday,outyear,tt2);
  cout<<tt.date_string<<" plus 5: "<<tt2.date_string<<endl;

  AddTime(360,1999,365,outday,outyear);
  JulianConvert(0.0,360,1999,tt);
  JulianConvert(0.0,outday,outyear,tt2);
  cout<<tt.date_string<<" plus 365: "<<tt2.date_string<<endl;

  AddTime(360,1999,731,outday,outyear);
  JulianConvert(0.0,360,1999,tt);
  JulianConvert(0.0,outday,outyear,tt2);
  cout<<tt.date_string<<" plus 731: "<<tt2.date_string<<endl;

  AddTime(360,1999,-5,outday,outyear);
  JulianConvert(0.0,360,1999,tt);
  JulianConvert(0.0,outday,outyear,tt2);
  cout<<tt.date_string<<" minus 5: "<<tt2.date_string<<endl;

  AddTime(360,1999,-365,outday,outyear);
  JulianConvert(0.0,360,1999,tt);
  JulianConvert(0.0,outday,outyear,tt2);
  cout<<tt.date_string<<" minus 365: "<<tt2.date_string<<endl;

  AddTime(360,1999,-731,outday,outyear);
  JulianConvert(0.0,360,1999,tt);
  JulianConvert(0.0,outday,outyear,tt2);
  cout<<tt.date_string<<" minus 731: "<<tt2.date_string<<endl;

}
//////////////////////////////////////////////////////////////////
/// \brief Tests JulianConvert method
//
void JulianConvertTest()
{
  time_struct tt;
  JulianConvert(0.0,0.0,2000,tt);
  cout<<"Jan 1, 2000 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(0.0,154,2004,tt);
  cout<<"Jun 2, 2004 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(366.0,0.0,2000,tt);
  cout<<"Jan 1, 2001 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(400.0,0.5,2000,tt);
  cout<<"Feb 4, 2001 @ 12:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(0.0,0.0,1999,tt);
  cout<<"Jan 1, 1999 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(366.0,0.0,1999,tt);
  cout<<"Jan 2, 2000 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(397.0,3.75,1999,tt);
  cout<<"Feb 5, 2000 @ 18:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(1.0,0.0,2000,tt);
  cout<<"Jan 2, 2000 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  ExitGracefully("JulianConvertTest",SIMULATION_DONE);
}
/////////////////////////////////////////////////////////////////
/// \brief Tests DecDaysToHours() function
//
void DecDaysTest()
{
  cout<<"1800.0   -->"<<DecDaysToHours(1800.0)<<endl;
  cout<<"37.5     -->"<<DecDaysToHours(37.5)<<endl;
  cout<<"18.25    -->"<<DecDaysToHours(18.25)<<endl;
  cout<<"37.041667-->"<<DecDaysToHours(37.041667)<<endl;
  cout<<"37.999   -->"<<DecDaysToHours(37.999)<<endl;
  cout<<"37.6293  -->"<<DecDaysToHours(37.6293)<<endl;
  ExitGracefully("DecDaysTest",SIMULATION_DONE);
}
/////////////////////////////////////////////////////////////////
/// \brief Tests SmartLookup() function
//
void SmartLookupUnitTest(){
  int n;
  double *aVals=new double [40];
  for (int i=0;i<40;i++){
    aVals[i]=(double)(i)*3.0;
  }
  n=SmartLookup(60.1,20,aVals,40);
  cout<<"Guess: 20, lookup=60: n="<<n<<" between "<<aVals[n]<<" and "<<aVals[n+1]<<endl;
  n=SmartLookup(60.1,19,aVals,40);
  cout<<"Guess: 19, lookup=60: n="<<n<<" between "<<aVals[n]<<" and "<<aVals[n+1]<<endl;
  n=SmartLookup(94,19,aVals,40);
  cout<<"Guess: 19, lookup=94: n="<<n<<" between "<<aVals[n]<<" and "<<aVals[n+1]<<endl;
  n=SmartLookup(-12.0,19,aVals,40);
  cout<<"Guess: 19, lookup=-12: n="<<n<<" between "<<aVals[n]<<" and "<<aVals[n+1]<<endl;
  n=SmartLookup(-12.0,120,aVals,40);
  cout<<"Guess: 120, lookup= -12: n="<<n<<" between "<<aVals[n]<<" and "<<aVals[n+1]<<endl;
  n=SmartLookup(94,0,aVals,40);
  cout<<"Guess: 0, lookup=94: n="<<n<<" between "<<aVals[n]<<" and "<<aVals[n+1]<<endl;
  n=SmartLookup(94,38,aVals,40);
  cout<<"Guess: 0, lookup=38: n="<<n<<" between "<<aVals[n]<<" and "<<aVals[n+1]<<endl;
  delete [] aVals;
  ExitGracefully("SmartLookupUnitTest",SIMULATION_DONE);
}
/////////////////////////////////////////////////////////////////
/// \brief Tests FixTimestep() function
//
void FixTimestepTest()
{
  cout.precision(12);
  double xx;
  xx = 0.041667; cout << "Fix Timestep :" << xx << " " << FixTimestep(xx) << " " << 1.0 / FixTimestep(xx) << endl;
  xx = 1.0;      cout << "Fix Timestep :" << xx << " " << FixTimestep(xx) << " " << 1.0 / FixTimestep(xx) << endl;
  xx = 0.333;    cout << "Fix Timestep :" << xx << " " << FixTimestep(xx) << " " << 1.0 / FixTimestep(xx) << endl;
  xx = 0.25;     cout << "Fix Timestep :" << xx << " " << FixTimestep(xx) << " " << 1.0 / FixTimestep(xx) << endl;
  xx = 0.003472; cout << "Fix Timestep :" << xx << " " << FixTimestep(xx) << " " << 1.0 / FixTimestep(xx) << endl;
  xx = 0.4;      cout << "Fix Timestep :" << xx << " " << FixTimestep(xx) << " " << 1.0 / FixTimestep(xx) << endl;
  ExitGracefully("FixTimestepTest", SIMULATION_DONE);
}
/////////////////////////////////////////////////////////////////
/// \brief Tests SmartIntervalSearch() function
//
void SmartIntervalTest()
{
  double *arr=new double[7];
  arr[0]=-12;
  arr[1]=-5.32;
  arr[2]=-2;
  arr[3]=-0.5;
  arr[4]=1.6;
  arr[5]=1.7;
  arr[6]=32;

  for(int ig=0;ig<7;ig++){
    cout<<"odd-sized array: "<< ig<<" "<<SmartIntervalSearch(3.2,arr,7,ig)<<endl;
    cout<<"odd-sized array: "<< ig<<" "<<SmartLookup(3.2,ig,arr,7)<<endl;
  }
  for(int ig=0;ig<6;ig++){
    cout<<"even-sized array: "<< ig<<" "<<SmartIntervalSearch(-1.7,arr,6,ig)<<endl;
    cout<<"even-sized array: "<< ig<<" "<<SmartLookup(-1.7,ig,arr,6)<<endl;
  }
  for(int ig=0;ig<6;ig++){
    cout<<"out of array: "<< ig<<" "<<SmartIntervalSearch(-17,arr,6,ig)<<endl;
    cout<<"out of array: "<< ig<<" "<<SmartLookup(-17,ig,arr,6)<<endl;
  }
  ExitGracefully("SmartIntervalTest", SIMULATION_DONE);
}
/////////////////////////////////////////////////////////////////
/// \brief Tests TriCumDist() and GammaCumDist2() functions
//
void GammaTest()
{
  ofstream GAMMA; GAMMA.open("GammaTest.csv");
  double dt=0.2;
  for(double t=0;t<10.0;t+=dt){
    GAMMA<<t<<",";

    for(double mu=1;mu<4;mu+=1.0){
      GAMMA<<(TriCumDist(t+dt,10,mu)-TriCumDist(t,10,mu))<<",";
    }
    //  for (double mu=1;mu<4;mu+=1.0){GAMMA<<(GammaCumDist2((t+dt)/mu*3.0,3)-GammaCumDist2(t/mu*3.0,3))<<",";}
    GAMMA<<endl;
  }
  ExitGracefully("GammaTest",SIMULATION_DONE);
}
/////////////////////////////////////////////////////////////////
/// \brief Tests ADRCumDist() functions
//
void ADRCumDistTest()
{
  ofstream OUT;
  OUT.open("ADRCumDistTest.csv");
  for (double t=0;t<10;t+=0.125){
    OUT<<t<<",";
    for (double D=0.03; D<0.3;D+=0.03){
      //OUT<<ADRCumDist(t,5,1.0,D)<<",";
      OUT<<ADRCumDist(t+0.125,5,1.0,D)-ADRCumDist(t,5,1.0,D)<<",";
    }
    OUT<<endl;
  }
  OUT.close();
  ExitGracefully("ADRCumDistTest",SIMULATION_DONE);
}
//////////////////////////////////////////////////////////////////
/// \brief A unit test function for clear sky radiation routine
//
void ClearSkyTest()
{
  double day_angle,declin,ecc,day_length;
  double slope=0.0;
  double solar_noon=0.0;
  double rad;
  ofstream CST;
  CST.open("ClearSkyTest.csv");
  for (double day=0;day<365;day++){
    for (double lat=0;lat<90;lat+=1.0){
      day_angle=CRadiation::DayAngle(day,1999);
      declin=   CRadiation::SolarDeclination(day_angle);
      ecc   = CRadiation::EccentricityCorr(day_angle);
      day_length = CRadiation::DayLength(lat*PI/180.0,declin);
      rad   =CRadiation::CalcETRadiation(lat*PI/180.0,lat*PI/180.0,declin,ecc,slope,solar_noon,day_length,0.0,true);//[MJ/m2/d]
      CST<<day<<","<<lat<<","<<rad<<endl;
    }
  }
  CST.close();
  ExitGracefully("ClearSkyTest",SIMULATION_DONE);
}

//////////////////////////////////////////////////////////////////
/// \brief A unit test function for optical air mass calculation
//
void OpticalAirMassTest()
{
  double lat,dec;
  ofstream OAM;
  ofstream OAM2;
  OAM.open ("OpticalAirMassTest.csv");
  /*OAM<<"dec,lat,OM"<<endl;
    double OM;
    double latstep=90/100.0;
    double decstep=(22.5*2.0)/100.0;
    for (dec=-22.5;dec<(22.5+0.5*decstep);dec+=decstep)
    {
    cout<<dec<<endl;
    for (lat=0;lat<90;lat+=latstep){
    OM=OpticalAirMass(lat*PI/180.0,dec*PI/180.0,0.0,true);
    OAM<<dec<<","<<lat<<","<<OM<<endl;
    }
    }*/
  //---------------------------------------------------
  // testing hourly for each month:
  //---------------------------------------------------
  lat=40;
  double jday;

  OAM<<"date,t,tsol,dec,OM(t),OM_avg"<<endl;
  for (double t=0;t<365;t+=1.0/24.0)
  {
    time_struct tt;
    JulianConvert(t,0.0,2001,tt);
    tt.model_time=t;

    if (tt.day_of_month==21)//21st day of the month
    {
      cout<<tt.month<<endl;
      double day_angle =CRadiation::DayAngle(jday,1999);
      dec=CRadiation::SolarDeclination(day_angle);
      double day_length = CRadiation::DayLength(lat*PI/180.0,dec*PI/180.0);
      double tsol=t-floor(t)-0.5;
      double OM=CRadiation::OpticalAirMass(lat*PI/180.0,dec*PI/180.0,day_length,tsol,false);
      OAM<<tt.date_string<<","<<t<<","<<tsol<<","<<dec<<","<<OM<<",";
      OAM<<CRadiation::OpticalAirMass(lat*PI/180,dec,day_length,tsol,true)<<endl;
    }
  }
  //---------------------------------------------------
  OAM.close();
  ExitGracefully("OpticalAirMassTest",SIMULATION_DONE);
}

//////////////////////////////////////////////////////////////////
/// \brief a unit test for shortwave radiation calculation
//
void ShortwaveTest()
{
  ofstream SHORT;
  SHORT.open("ShortwaveTest_ET.csv");
  if(SHORT.fail()){ ExitGracefully("cannot open file.",BAD_DATA); }
  int year=2001; //no leap year
  double SW;
  double slope;
  double aspect;
  double dew_pt=GetDewPointTemp(10,0.5);
  double ET_rad;
  double day_angle,declin,ecc,day_length;
  double latrad;
  double SW2(0),lateq,solar_noon;
  double tstep=1.0/100.0;

  //header
  //double lat=39.7555; //Golden colorado
  //double lat=44.0521; //Eugene oregon
  double lat=70;

  latrad=lat*PI/180.0;
  SHORT<<"doy,";
  for(double a=0;a<=270; a+=90){
    for(double s=0;s<=90;s+=15.0){
      SHORT<<"slope_"<<s<< "(asp="<<a<<"),";
    }
  }
  double conv=MJ_PER_D_TO_WATT;
  SHORT<<endl;
  for (double day=0;day<365;day+=tstep){
    if(int(day*24+0.001) % 24==0){cout<<day<<endl;}
    SHORT<<day<<",";
    //These aren't impacted by slope / aspect, and are OK
    day_angle  = CRadiation::DayAngle(day,year);
    declin     = CRadiation::SolarDeclination(day_angle);
    ecc        = CRadiation::EccentricityCorr(day_angle);
    day_length = CRadiation::DayLength(latrad,declin);
    for(double a=0;a<=270; a+=90){
      aspect=a*PI/180;//conv to rads
      for(double s=0;s<=90;s+=15.0){
        slope=s*PI/180.0;//conv to rads
        lateq      = CRadiation::CalculateEquivLatitude(latrad,slope,aspect);
        solar_noon = CRadiation::CalculateSolarNoon    (latrad,slope,aspect); //relative to actual location
        
        SW=CRadiation::ClearSkySolarRadiation(day+0.001,tstep,latrad,lateq,slope,aspect,day_angle,day_length,solar_noon,dew_pt,ET_rad,(tstep==1.0));
        SW=ET_rad;
        SHORT<<SW*conv<<",";
        //SHORT<<solar_noon*12/PI<<","; //report solar noon, in hrs
        //SHORT<<solar_noon*12/PI<<","; //generate solar noon
      }
    }
    SHORT<<endl;
  }


//check against Lee, 1964: confirms that dingman E23 is correct and Lee eqn 6 is wrong (asin, not acos)
  cout<<"test1: (should be 26.21666667	28.83333333)"<<endl;
  latrad=37.76666667*PI/180; aspect=-336.0*PI/180.0; slope=31.33333333*PI/180;
  lateq=CRadiation::CalculateEquivLatitude(latrad,slope,aspect);
  solar_noon=CRadiation::CalculateSolarNoon(latrad,slope,aspect);
  cout<<-(latrad-lateq)*180.0/PI<<" "<<solar_noon*EARTH_ANG_VEL*180.0/PI<<endl;

  cout<<"test2: (should be -22.83333333	-28.93333333)"<<endl;
  latrad=37.76666667*PI/180;	aspect=-124*PI/180.0;	slope=34.33333333*PI/180.0;
  CRadiation::CalculateEquivLatitude(latrad,slope,aspect);
  solar_noon=CRadiation::CalculateSolarNoon(latrad,slope,aspect);
  cout<<-(latrad-lateq)*180.0/PI<<" "<<solar_noon*EARTH_ANG_VEL*180.0/PI<<endl;

  cout<<"test3: (should be 30.06666667	-40.83333333)"<<endl;
  latrad=37.76666667*PI/180;	aspect=-24*PI/180.0;	slope=37.5*PI/180.0;
  CRadiation::CalculateEquivLatitude(latrad,slope,aspect);
  solar_noon=CRadiation::CalculateSolarNoon(latrad,slope,aspect);
  cout<<-(latrad-lateq)*180.0/PI<<" "<<solar_noon*EARTH_ANG_VEL*180.0/PI<<endl;

  cout<<"test4: (should be -23.33333333	-21.3)"<<endl;
  latrad=37.76666667*PI/180;	aspect=-135*PI/180.0;	slope=30*PI/180.0;
  CRadiation::CalculateEquivLatitude(latrad,slope,aspect);
  solar_noon=CRadiation::CalculateSolarNoon(latrad,slope,aspect);
  cout<<-(latrad-lateq)*180.0/PI<<" "<<solar_noon*EARTH_ANG_VEL*180.0/PI<<endl;




  //---------------------------------------------------
  // testing hourly for each month:
  //---------------------------------------------------
  /*double lat=80;
    double dday,jday;
    int dmon,junk;double ET_rad;
    SHORT<<"date,t,tfake,SHORT(t),SHORT_avg"<<endl;
    for (double t=0;t<365;t+=1.0/48.0)
    {
    time_struct tt;
    JulianConvert(t,0.0,2001,tt);
    string thisdate=tt.date_string;
    if (ceil(dday)==21)//21st day of the month
    {
    cout<<dmon<<endl;
    SW=ClearSkySolarRadiation(jday,tstep,year,lat,slope,aspect,dew_pt,ET_rad,false);
    SHORT<<thisdate<<","<<t<<","<<t-floor(t)+dmon-1<<","<<SW<<",";
    SHORT<<ClearSkySolarRadiation(jday,year,lat,slope,aspect,dew_pt,ET_rad,true)<<endl;
    }
    }*/
  //---------------------------------------------------
  SHORT.close();
  ExitGracefully("ShortwaveTest",SIMULATION_DONE);
}



//////////////////////////////////////////////////////////////////
/// \brief a routine for taking input of day, slope,aspect, dewpoint, and albedo and calculating shortwave rafiation parameters
/// \details input file (ShortwaveInput.csv in working directory) in the following format:
/// input file
/// # of entries [tstep]
/// day,year,lat(dec),slope(m/m),aspect (degrees from north),dewpoint,albedo
/// output file
/// day, year, lat,slope,aspect,declin,ecc,solar_time,OAM,Ketp,Ket
//
void ShortwaveGenerator()
{
  double day,latrad,slope,aspect,lateq,declin,ecc,Ketp,Ket;
  double day_angle,Mopt,t_sol,TIR,albedo, dew_pt;
  double day_length, solar_noon;
  int year;
  ifstream IN;
  ofstream SHORT;

  IN.open   ("ShortwaveInput.csv");
  SHORT.open("ShortwaveGenerator.csv");
  SHORT<<"day[d], year, lat[dec],slope,aspect,declin[rad],ecc[rad],solar_time[d],albedo,dew_pt,Mopt[-],Ketp[MJ/m2/d],Ket[MJ/m2/d],Kcs[MJ/m2/d]"<<endl;
  int   Len,line(0);double ET_rad;
  char *s[MAXINPUTITEMS];
  CParser *p=new CParser(IN,line);
  p->Tokenize(s,Len);
  int NumLines=s_to_i(s[0]);
  double tstep=s_to_d(s[1]);
  p->Tokenize(s,Len);//header
  for (int i=0;i<NumLines;i++)
  {
    cout.precision(2);
    if ((NumLines>10) && (i%(NumLines/10)==0)){cout<<(double)(i)/NumLines*100<<"%"<<endl;}
    p->Tokenize(s,Len);
    day     =s_to_d(s[0]);
    year    =s_to_i(s[1]);
    latrad  =s_to_d(s[2])*PI/180;
    slope   =atan(s_to_d(s[3]));
    aspect  =s_to_d(s[4])*PI/180;
    lateq   = asin(cos(slope)*sin(latrad) + sin(slope)*cos(latrad)*cos(aspect));
    dew_pt  =s_to_d(s[5]);
    albedo  =s_to_d(s[6]);

    day_angle= CRadiation::DayAngle(day,year);
    declin   = CRadiation::SolarDeclination(day_angle);
    ecc      = CRadiation::EccentricityCorr(day_angle);
    day_length = CRadiation::DayLength(latrad,declin);
    t_sol    = day-floor(day)-0.5;

    double denom=cos(slope)*cos(latrad) - sin(slope)*sin(latrad)*cos(aspect);
    if (denom==0.0){denom = REAL_SMALL;}
    solar_noon = -atan(sin(slope)*sin(aspect)/denom)/EARTH_ANG_VEL;
    if (solar_noon> 0.5){solar_noon-=1.0;}
    if (solar_noon<-0.5){solar_noon+=1.0;}


    Mopt =CRadiation::OpticalAirMass (latrad,declin,                       day_length,t_sol,false);
    Ketp =CRadiation::CalcETRadiation(latrad,lateq,declin,ecc,slope,solar_noon,day_length,t_sol,false);
    Ket  =CRadiation::CalcETRadiation(latrad,lateq,declin,ecc,0.0  ,0.0       ,day_length,t_sol,false);

    TIR  =CRadiation::ClearSkySolarRadiation(day,tstep,latrad,lateq,tan(slope),aspect,day_angle,day_length,solar_noon,dew_pt,ET_rad,tstep==1.0);

    SHORT<<day<<","<<year<<","<<latrad*180.0/PI<<","<<tan(slope)<<","<<aspect*180/PI;
    SHORT<<","<<declin<<","<<ecc<<","<<t_sol<<",";
    SHORT<<albedo<<","<<dew_pt<<","<<Mopt<<","<<Ketp<<","<<Ket<<","<<TIR<<endl;
  }
  cout<<"100% - done"<<endl;
  IN.close();
  SHORT.close();
  ExitGracefully("ShortwaveGenerator",SIMULATION_DONE);
}

