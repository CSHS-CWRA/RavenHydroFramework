/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "ParseLib.h"
#include "Radiation.h"
#include "GlobalParams.h"
#include "UnitTesting.h"

void RavenUnitTesting(const optStruct &Options)
{
  //cout<<"RAVEN UNIT TESTING MODE"<<endl;
  //for (int i=0;i<23;i++){cout<<ravenASCII[i];}
  //uncomment one for use:
  //DateTest();
  //OpticalAirMassTest();
  //ClearSkyTest();
  //ShortwaveTest();
  //ShortwaveGenerator();
  //JulianConvertTest();
  //SmartIntervalTest();
  //AddTimeTest( );
  //GammaTest();
  //BarycentricWeights();
  //TestInversion();
  //ADRCumDistTest();
  /*cout<<"STRINGISLONG TEST: 1023 "<< StringIsLong("1023")<<endl;
  cout<<"STRINGISLONG TEST: -1 "<< StringIsLong("-1")<<endl;
  cout<<"STRINGISLONG TEST: hamburger "<< StringIsLong("hamburger")<<endl;
  cout<<"STRINGISLONG TEST: 10ham "<< StringIsLong("10ham")<<endl;*/
  //TestEnthalpyTempConvert();
  //TestConvectionSolution();
  //TestGammaSampling();
  //TestWetBulbTemps();
  //TestDateStrings();

}
/////////////////////////////////////////////////////////////////
/// \brief Tests DateStringToTimeStruct() function
//
void DateTest()
{
  time_struct tt;
  cout<<"2012-03-27 12:43:02.01"<<endl;
  tt=DateStringToTimeStruct("2012-03-27","12:43:02.01",StringToCalendar("PROLEPTIC_GREGORIAN"));
  cout<<tt.date_string<<" month, day, year: "<<tt.month <<","<<tt.day_of_month<<","<<tt.year<<" julian: "<<tt.julian_day<<endl;

  cout<<"2000-12-31 12:00:00"<<endl;
  tt=DateStringToTimeStruct("2000-12-31","12:00:00",StringToCalendar("PROLEPTIC_GREGORIAN"));
  cout<<tt.date_string<<" month, day, year: "<<tt.month <<","<<tt.day_of_month<<","<<tt.year<<" julian: "<<tt.julian_day<<endl;

  cout<<"1997/12/31 23:00:00.0000"<<endl;
  tt=DateStringToTimeStruct("1997/12/31","23:00:00.0000",StringToCalendar("PROLEPTIC_GREGORIAN"));
  cout<<tt.date_string<<" month, day, year: "<<tt.month <<","<<tt.day_of_month<<","<<tt.year<<" julian: "<<tt.julian_day<<endl;

  cout<<"0000/01/01 24:00:00.0000"<<endl;
  tt=DateStringToTimeStruct("0000/01/01","23:59:00.0000",StringToCalendar("PROLEPTIC_GREGORIAN"));
  cout<<tt.date_string<<" month, day, year: "<<tt.month <<","<<tt.day_of_month<<","<<tt.year<<" julian: "<<tt.julian_day<<endl;
  ExitGracefully("DateTest",SIMULATION_DONE);
}
void TestDateStrings() {
  string test;
  cout<<"April 23 should be 112, April 2 should be 91"<<endl;
  test="Apr-23";
  cout<<test<<" julian: "<<GetJulianDayFromMonthYear(test, StringToCalendar("PROLEPTIC_GREGORIAN"))<<endl;
  test="Apr-2";
  cout<<test<<" julian: "<<GetJulianDayFromMonthYear(test, StringToCalendar("PROLEPTIC_GREGORIAN"))<<endl;
  test="04-23";
  cout<<test<<" julian: "<<GetJulianDayFromMonthYear(test, StringToCalendar("PROLEPTIC_GREGORIAN"))<<endl;
  test="4-23";
  cout<<test<<" julian: "<<GetJulianDayFromMonthYear(test, StringToCalendar("PROLEPTIC_GREGORIAN"))<<endl;
  test="4-2";
  cout<<test<<" julian: "<<GetJulianDayFromMonthYear(test, StringToCalendar("PROLEPTIC_GREGORIAN"))<<endl;
  test="4-02";
  cout<<test<<" julian: "<<GetJulianDayFromMonthYear(test, StringToCalendar("PROLEPTIC_GREGORIAN"))<<endl;
  test="92";
  cout<<test<<" julian: "<<GetJulianDayFromMonthYear(test, StringToCalendar("PROLEPTIC_GREGORIAN"))<<endl;
  test="Apr-30";
  cout<<test<<" julian: "<<GetJulianDayFromMonthYear(test, StringToCalendar("PROLEPTIC_GREGORIAN"))<<endl;
  ExitGracefully("TestDateStrings",SIMULATION_DONE);
}
void AddTimeTest( ) {

  time_struct tt,tt2;
  double      outday;
  int         outyear;

  AddTime(360,1999,5,StringToCalendar("PROLEPTIC_GREGORIAN"),outday,outyear);
  JulianConvert(0.0,360,1999,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("PROLEPTIC_GREGORIAN"),tt2);
  cout<<tt.date_string<<" plus 5:      "<<tt2.date_string<<"   (expected   0th of 2000 (  leap) = 2000-01-01)"<<endl;

  AddTime(360,1999,365,StringToCalendar("PROLEPTIC_GREGORIAN"),outday,outyear);
  JulianConvert(0.0,360,1999,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("PROLEPTIC_GREGORIAN"),tt2);
  cout<<tt.date_string<<" plus 365:    "<<tt2.date_string<<"   (expected 360th of 2000 (  leap) = 2000-12-26)"<<endl;

  AddTime(360,1999,731,StringToCalendar("PROLEPTIC_GREGORIAN"),outday,outyear);
  JulianConvert(0.0,360,1999,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("PROLEPTIC_GREGORIAN"),tt2);
  cout<<tt.date_string<<" plus 731:    "<<tt2.date_string<<"   (expected 360th of 2001 (noleap) = 2001-12-27)"<<endl;

  AddTime(360,1999,-5,StringToCalendar("PROLEPTIC_GREGORIAN"),outday,outyear);
  JulianConvert(0.0,360,1999,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("PROLEPTIC_GREGORIAN"),tt2);
  cout<<tt.date_string<<" minus 5:     "<<tt2.date_string<<"   (expected 355th of 1999 (noleap) = 1999-12-22)"<<endl;

  AddTime(360,1999,-365,StringToCalendar("PROLEPTIC_GREGORIAN"),outday,outyear);
  JulianConvert(0.0,360,1999,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("PROLEPTIC_GREGORIAN"),tt2);
  cout<<tt.date_string<<" minus 365:   "<<tt2.date_string<<"   (expected 360th of 1998 (noleap) = 1998-12-27)"<<endl;

  AddTime(360,1999,-731,StringToCalendar("PROLEPTIC_GREGORIAN"),outday,outyear);
  JulianConvert(0.0,360,1999,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("PROLEPTIC_GREGORIAN"),tt2);
  cout<<tt.date_string<<" minus 731:   "<<tt2.date_string<<"   (expected 359th of 1997 (noleap) = 1997-12-26)"<<endl;

  AddTime(0,1,733774,StringToCalendar("PROLEPTIC_GREGORIAN"),outday,outyear);  // must be 2010-01-03
  JulianConvert(0.0,0,1,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("PROLEPTIC_GREGORIAN"),tt2);
  cout<<tt.date_string<<" plus 733774: "<<tt2.date_string<<"   (expected   2th of 2010          = 2010-01-03)"<<endl;

  // GREGORIAN has 10 days missing (4 October 1582 is followed by 15 October 1582)
  // and leap years where different before 1582
  AddTime(0,1,733774,StringToCalendar("GREGORIAN"),outday,outyear);  // must be 2010-01-01
  JulianConvert(0.0,0,1,StringToCalendar("GREGORIAN"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("GREGORIAN"),tt2);
  cout<<tt.date_string<<" plus 733774: "<<tt2.date_string<<"   (expected   0th of 2010          = 2010-01-01)"<<endl;

  // STANDARD is same as GREGORIAN
  AddTime(0,1,733774,StringToCalendar("STANDARD"),outday,outyear);  // must be 2010-01-01
  JulianConvert(0.0,0,1,StringToCalendar("STANDARD"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("STANDARD"),tt2);
  cout<<tt.date_string<<" plus 733774: "<<tt2.date_string<<"   (expected   0th of 2010          = 2010-01-01)"<<endl;


  AddTime(0,1,689945.0,StringToCalendar("GREGORIAN"),outday,outyear);  //must be 1890-01-01
  JulianConvert(0.0,0,1,StringToCalendar("GREGORIAN"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("GREGORIAN"),tt2);
  cout<<tt.date_string<<" plus 689945: "<<tt2.date_string<<"   (expected   0th of 1890          = 1890-01-01)"<<endl;

  AddTime(0,1,489800.0,StringToCalendar("GREGORIAN"),outday,outyear);  //1342, 1, 1, 0, 0, 0, 0, 1, 1)
  JulianConvert(0.0,0,1,StringToCalendar("GREGORIAN"),tt);
  JulianConvert(0.0,outday,outyear,StringToCalendar("GREGORIAN"),tt2);
  cout<<tt.date_string<<" plus 489800: "<<tt2.date_string<<"   (expected   0th of 1342          = 1342-01-01)"<<endl;

}
//////////////////////////////////////////////////////////////////
/// \brief Tests JulianConvert method
//
void JulianConvertTest()
{
  time_struct tt;
  JulianConvert(0.0,0.0,2000,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  cout<<"Jan 1, 2000 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(0.0,154,2004,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  cout<<"Jun 2, 2004 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(366.0,0.0,2000,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  cout<<"Jan 1, 2001 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(400.0,0.5,2000,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  cout<<"Feb 4, 2001 @ 12:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(0.0,0.0,1999,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  cout<<"Jan 1, 1999 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(366.0,0.0,1999,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  cout<<"Jan 2, 2000 @ 0:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(397.0,3.75,1999,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
  cout<<"Feb 5, 2000 @ 18:00: "<<tt.month<<" "<<tt.day_of_month<<", "<<tt.year<<" julian:"<<tt.julian_day<<" "<<tt.date_string<<endl;

  JulianConvert(1.0,0.0,2000,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
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
  }
  for(int ig=0;ig<6;ig++){
    cout<<"even-sized array: "<< ig<<" "<<SmartIntervalSearch(-1.7,arr,6,ig)<<endl;
  }
  for(int ig=0;ig<6;ig++){
    cout<<"out of array: "<< ig<<" "<<SmartIntervalSearch(-17,arr,6,ig)<<endl;
  }
  ExitGracefully("SmartIntervalTest", SIMULATION_DONE);
}
/////////////////////////////////////////////////////////////////
/// \brief Tests TriCumDist() and GammaCumDist2() functions
//
void GammaTest()
{
  ofstream GAMMA; GAMMA.open("GammaTest.csv");
  double dt=0.1;
  for(double t=0;t<20.0;t+=dt){
    GAMMA<<t<<",";

    /*for(double mu=1;mu<4;mu+=1.0){
      GAMMA<<(TriCumDist(t+dt,10,mu)-TriCumDist(t,10,mu))<<",";
    }*/
    for (double mu=1;mu<=4;mu+=1.0){
      GAMMA<<(GammaCumDist(t+dt,3,(3-1)/mu)-GammaCumDist(t,3,(3-1)/mu))<<","<<GammaCumDist(t+dt,3,(3-1)/mu)<<",";
    }
    GAMMA<<endl;
  }
  ExitGracefully("GammaTest",SIMULATION_DONE);
}
/////////////////////////////////////////////////////////////////
/// \brief Tests ADRCumDist() functions
//
void ADRCumDistTest()
{
  ofstream TEST;
  TEST.open("ADRCumDistTest.csv");
  double *v=new double[10];
  for(int j=0;j<10;j++) {
    v[j]=1.0+5.0*((double)(rand())/RAND_MAX);
  }
  for(double t=0;t<10;t+=0.125) {
    TEST<<t<<",";
    for(double D=0.03; D<0.3;D+=0.03) {
      //TEST<<ADRCumDist(t,5,1.0,D)<<",";
      TEST<<ADRCumDist(t+0.125,5,1.0,D)-ADRCumDist(t,5,1.0,D)<<",";
    }
    TEST<<endl;
  }
  TEST<<endl;
  for(int j=0;j<10;j++) {
    TEST<<v[j]<<",";
  }
  TEST<<endl;
  for(double t=0;t<10;t+=0.125) {
    TEST<<t<<",";
    for(double D=0.12; D<1.2;D+=0.12) {
      //TEST<< TimeVaryingADRCumDist(t      ,5,v,10,D,1.0)<<",";
      TEST<<TimeVaryingADRCumDist(t+0.125,5,v,10,D,1.0)-
           TimeVaryingADRCumDist(t      ,5,v,10,D,1.0)<<",";
    }
    TEST<<endl;

  }
  TEST<<endl;
  //paper figure
  for(int j=0;j<10;j++) {
    v[j]=1.0;
  }
  v[1]=4; v[2]=2; v[3]=5; v[4]=4; v[5]=3; v[6]=3; v[7]=5;
  for(double t=0;t<6.0;t+=0.025) {
    TEST<<t<<",";
    double D=1.0;
    TEST<< TimeVaryingADRCumDist(t      ,5,v,10,D,1.0)<<",";
    TEST<< TimeVaryingADRCumDist(t+0.025,5,v,10,D,1.0)-
          TimeVaryingADRCumDist(t      ,5,v,10,D,1.0)<<",";
    TEST<<endl;

  }
  TEST.close();
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
      day_angle  = CRadiation::DayAngle(day,1999,StringToCalendar("PROLEPTIC_GREGORIAN"));
      declin     = CRadiation::SolarDeclination(day_angle);
      ecc        = CRadiation::EccentricityCorr(day_angle);
      day_length = CRadiation::DayLength(lat*DEGREES_TO_RADIANS,declin);
      rad        = CRadiation::CalcETRadiation(lat*DEGREES_TO_RADIANS,lat*DEGREES_TO_RADIANS,declin,ecc,slope,solar_noon,day_length,0.0,true);//[MJ/m2/d]
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
    OM=OpticalAirMass(lat*DEGREES_TO_RADIANS,dec*DEGREES_TO_RADIANS,0.0,true);
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
    JulianConvert(t,0.0,2001,StringToCalendar("PROLEPTIC_GREGORIAN"),tt);
    tt.model_time=t;

    if (tt.day_of_month==21)//21st day of the month
    {
      cout<<tt.month<<endl;

      double day_angle  = CRadiation::DayAngle(jday,1999,StringToCalendar("PROLEPTIC_GREGORIAN"));
      dec               = CRadiation::SolarDeclination(day_angle);
      double day_length = CRadiation::DayLength(lat*DEGREES_TO_RADIANS,dec*DEGREES_TO_RADIANS);
      double tsol       = t-floor(t)-0.5;
      double OM         = CRadiation::OpticalAirMass(lat*DEGREES_TO_RADIANS,dec*DEGREES_TO_RADIANS,day_length,tsol,false);

      OAM<<tt.date_string<<","<<t<<","<<tsol<<","<<dec<<","<<OM<<",";
      OAM<<CRadiation::OpticalAirMass(lat*DEGREES_TO_RADIANS,dec,day_length,tsol,true)<<endl;
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
  double ET_rad,ET_flat;
  double day_angle,declin,ecc,day_length;
  double latrad;
  double lateq,solar_noon;
  double tstep=1.0/100.0;

  //header
  //double lat=39.7555; //Golden colorado
  //double lat=44.0521; //Eugene oregon
  double lat=70;

  latrad=lat*DEGREES_TO_RADIANS;
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
    day_angle  = CRadiation::DayAngle(day,year,StringToCalendar("PROLEPTIC_GREGORIAN"));
    declin     = CRadiation::SolarDeclination(day_angle);
    ecc        = CRadiation::EccentricityCorr(day_angle);
    day_length = CRadiation::DayLength(latrad,declin);
    for(double a=0;a<=270; a+=90){
      aspect=a*DEGREES_TO_RADIANS;//conv to rads
      for(double s=0;s<=90;s+=15.0){
        slope=s*DEGREES_TO_RADIANS;//conv to rads
        lateq      = CRadiation::CalculateEquivLatitude(latrad,slope,aspect);
        solar_noon = CRadiation::CalculateSolarNoon    (latrad,slope,aspect); //relative to actual location

        SW=CRadiation::ClearSkySolarRadiation(day+0.001,tstep,latrad,lateq,slope,aspect,day_angle,day_length,solar_noon,dew_pt,ET_rad,ET_flat,(tstep==1.0));
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
  latrad=37.76666667*DEGREES_TO_RADIANS; aspect=-336.0*DEGREES_TO_RADIANS; slope=31.33333333*DEGREES_TO_RADIANS;
  lateq=CRadiation::CalculateEquivLatitude(latrad,slope,aspect);
  solar_noon=CRadiation::CalculateSolarNoon(latrad,slope,aspect);
  cout<<-(latrad-lateq)*RADIANS_TO_DEGREES<<" "<<solar_noon*EARTH_ANG_VEL*RADIANS_TO_DEGREES<<endl;

  cout<<"test2: (should be -22.83333333	-28.93333333)"<<endl;
  latrad=37.76666667*DEGREES_TO_RADIANS;	aspect=-124*DEGREES_TO_RADIANS;	slope=34.33333333*DEGREES_TO_RADIANS;
  CRadiation::CalculateEquivLatitude(latrad,slope,aspect);
  solar_noon=CRadiation::CalculateSolarNoon(latrad,slope,aspect);
  cout<<-(latrad-lateq)*RADIANS_TO_DEGREES<<" "<<solar_noon*EARTH_ANG_VEL*RADIANS_TO_DEGREES<<endl;

  cout<<"test3: (should be 30.06666667	-40.83333333)"<<endl;
  latrad=37.76666667*DEGREES_TO_RADIANS;	aspect=-24*DEGREES_TO_RADIANS;	slope=37.5*DEGREES_TO_RADIANS;
  CRadiation::CalculateEquivLatitude(latrad,slope,aspect);
  solar_noon=CRadiation::CalculateSolarNoon(latrad,slope,aspect);
  cout<<-(latrad-lateq)*RADIANS_TO_DEGREES<<" "<<solar_noon*EARTH_ANG_VEL*RADIANS_TO_DEGREES<<endl;

  cout<<"test4: (should be -23.33333333	-21.3)"<<endl;
  latrad=37.76666667*DEGREES_TO_RADIANS;	aspect=-135*DEGREES_TO_RADIANS;	slope=30*DEGREES_TO_RADIANS;
  CRadiation::CalculateEquivLatitude(latrad,slope,aspect);
  solar_noon=CRadiation::CalculateSolarNoon(latrad,slope,aspect);
  cout<<-(latrad-lateq)*RADIANS_TO_DEGREES<<" "<<solar_noon*EARTH_ANG_VEL*RADIANS_TO_DEGREES<<endl;




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
    SW=ClearSkySolarRadiation(jday,tstep,year,lat,slope,aspect,dew_pt,ET_rad,ET_flat,false);
    SHORT<<thisdate<<","<<t<<","<<t-floor(t)+dmon-1<<","<<SW<<",";
    SHORT<<ClearSkySolarRadiation(jday,year,lat,slope,aspect,dew_pt,ET_rad,ET_flat,true)<<endl;
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
  ifstream SHIN;
  ofstream SHORT;

  SHIN.open   ("ShortwaveInput.csv");
  SHORT.open("ShortwaveGenerator.csv");
  SHORT<<"day[d], year, lat[dec],slope,aspect,declin[rad],ecc[rad],solar_time[d],albedo,dew_pt,Mopt[-],Ketp[MJ/m2/d],Ket[MJ/m2/d],Kcs[MJ/m2/d]"<<endl;
  int   Len,line(0);double ET_rad,ET_flat;
  char *s[MAXINPUTITEMS];
  CParser *p=new CParser(SHIN,line);
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
    latrad  =s_to_d(s[2])*DEGREES_TO_RADIANS;
    slope   =atan(s_to_d(s[3]));
    aspect  =s_to_d(s[4])*DEGREES_TO_RADIANS;
    lateq   = asin(cos(slope)*sin(latrad) + sin(slope)*cos(latrad)*cos(aspect));
    dew_pt  =s_to_d(s[5]);
    albedo  =s_to_d(s[6]);

    day_angle= CRadiation::DayAngle(day,year,StringToCalendar("PROLEPTIC_GREGORIAN"));
    declin   = CRadiation::SolarDeclination(day_angle);
    ecc      = CRadiation::EccentricityCorr(day_angle);
    day_length = CRadiation::DayLength(latrad,declin);
    t_sol    = day-floor(day)-0.5;

    double denom=cos(slope)*cos(latrad) - sin(slope)*sin(latrad)*cos(aspect);
    if (denom==0.0){denom = REAL_SMALL;}
    solar_noon = -atan(sin(slope)*sin(aspect)/denom)/EARTH_ANG_VEL;
    if (solar_noon> 0.5){solar_noon-=1.0;}
    if (solar_noon<-0.5){solar_noon+=1.0;}


    Mopt =CRadiation::OpticalAirMass (latrad,declin,                           day_length,t_sol,false);
    Ketp =CRadiation::CalcETRadiation(latrad,lateq,declin,ecc,slope,solar_noon,day_length,t_sol,false);
    Ket  =CRadiation::CalcETRadiation(latrad,lateq,declin,ecc,0.0  ,0.0       ,day_length,t_sol,false);

    TIR  =CRadiation::ClearSkySolarRadiation(day,tstep,latrad,lateq,tan(slope),aspect,day_angle,day_length,solar_noon,dew_pt,ET_rad,ET_flat,tstep==1.0);

    SHORT<<day<<","<<year<<","<<latrad*RADIANS_TO_DEGREES<<","<<tan(slope)<<","<<aspect*RADIANS_TO_DEGREES;
    SHORT<<","<<declin<<","<<ecc<<","<<t_sol<<",";
    SHORT<<albedo<<","<<dew_pt<<","<<Mopt<<","<<Ketp<<","<<Ket<<","<<TIR<<endl;
  }
  cout<<"100% - done"<<endl;
  SHIN.close();
  SHORT.close();
  ExitGracefully("ShortwaveGenerator",SIMULATION_DONE);
}

void BarycentricWeights() {
  double sum=0.0;
  int    N=5;
  double _aWeights[5];
  double aVals[4];
  ofstream TEST;
  TEST.open("BarycentricWeights.csv");
  for(int m=0; m<10000;m++) {
    sum=0;
    for(int i=0;i<4;i++) { aVals[i]=rand()/(double)(RAND_MAX); }
    for(int q=0; q<N-1;q++) {
      _aWeights[q]=(1.0-sum)*(1.0-pow(1.0-aVals[q],1.0/(N-q)));
      sum+=_aWeights[q];
    }
    _aWeights[N-1]=1.0-sum;
    for(int q=0; q<N;q++) {
      TEST<<_aWeights[q]<<",";
    }TEST<<endl;
  }
  TEST.close();
  ExitGracefully("BarycentricWeights",SIMULATION_DONE);
}

void TestWetBulbTemps() {
  double T,RH;
  double P=100;
  cout<<"T ,RH ,Tw"<<endl;
  T=30; RH=0.5;
  cout<<T<<" ,"<<RH<<" "<<GetWetBulbTemperature(P,T,RH)<<endl;
  T=20; RH=0.9;
  cout<<T<<" ,"<<RH<<" "<<GetWetBulbTemperature(P,T,RH)<<endl;
  T=20; RH=1.0;
  cout<<T<<" ,"<<RH<<" "<<GetWetBulbTemperature(P,T,RH)<<endl;
  T=5; RH=0.35;
  cout<<T<<" ,"<<RH<<" "<<GetWetBulbTemperature(P,T,RH)<<endl;
  ExitGracefully("UnitTesting:: TestWetBulbTemps",SIMULATION_DONE);
}
