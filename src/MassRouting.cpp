/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
------------------------------------------------------------------
	routing of mass in catchment and channel
----------------------------------------------------------------*/

#include "SubBasin.h"
#include "Model.h"
#include "Transport.h"
//////////////////////////////////////////////////////////////////
/// \brief Initialize mass routing variables
///
void CTransportModel::InitializeRoutingVars()
{
  int nSB=pModel->GetNumSubBasins();
  aMinHist=new double **[nSB];
  aMlatHist=new double **[nSB];
  aMout    =new double **[nSB];
  for (int p=0;p<nSB;p++)
  {
    aMinHist [p]=new double *[nConstituents];
    aMlatHist[p]=new double *[nConstituents];
    aMout    [p]=new double *[nConstituents];
    int nMinHist=pModel->GetSubBasin(p)->GetInflowHistorySize();
    int nMlatHist=pModel->GetSubBasin(p)->GetLatHistorySize();
    int nSegments=pModel->GetSubBasin(p)->GetNumSegments();
    for (int c=0;c<nConstituents;c++){
      aMinHist [p][c]=new double [nMinHist];
      aMlatHist[p][c]=new double [nMlatHist];
      aMout    [p][c]=new double [nSegments];
      for (int i=0; i<nMinHist;i++){aMinHist[p][c][i]=0.0;}
      for (int i=0; i<nMinHist;i++){aMinHist[p][c][i]=0.0;}
      for (int i=0; i<nSegments;i++){aMout  [p][c][i]=0.0;}
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Delete mass routing variables
///
void CTransportModel::DeleteRoutingVars()
{
  int nSB=pModel->GetNumSubBasins();
  if (aMinHist!=NULL){
    for (int p=0;p<nSB;p++)
    {
      for (int c=0;c<nConstituents;c++){
        delete [] aMinHist [p][c];
        delete [] aMlatHist[p][c];
        delete [] aMout    [p][c];
      }
      delete [] aMinHist [p];
      delete [] aMlatHist[p];
      delete [] aMout    [p];
    }
    delete [] aMinHist; aMinHist =NULL;
    delete [] aMlatHist;aMlatHist=NULL;
    delete [] aMout;    aMout    =NULL;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Set upstream mass loading to primary channel and updates mass loading history
///
/// \param p    subbasin index
/// \param aMin array of rates of upstream mass loading for the current timestep [mg/d] [size: nConstituents]
//
void   CTransportModel::SetMassInflows    (const int p, const double *aMin)
{
  int nMinHist=pModel->GetSubBasin(p)->GetInflowHistorySize();
  for (int c=0;c<nConstituents;c++){
    for (int n=nMinHist-1;n>0;n--){
      aMinHist[p][c][n]=aMinHist[p][c][n-1];
    }
    aMinHist[p][c][0]=aMin[c];
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets lateral mass loading to primary channel and updates mass loading history
///
/// \param p    subbasin index
/// \param aMin array of lateral mass loading rates for the current timestep [mg/d] [size: nConstituents]
//
void   CTransportModel::SetLateralInfluxes(const int p, const double *aMlat)
{
  int nMlatHist=pModel->GetSubBasin(p)->GetLatHistorySize();
  for (int c=0;c<nConstituents;c++){
    for (int n=nMlatHist-1;n>0;n--){
      aMlatHist[p][c][n]=aMlatHist[p][c][n-1];
    }
    aMlatHist[p][c][0]=aMlat[c];
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Creates aMoutnew [m^3/s], an array of point measurements for outflow at downstream end of each river segment
/// \details Array represents measurements at the end of the current timestep. if (catchment_routing==distributed), 
/// lateral flow Qlat (e.g., runoff, interflow, or baseflow) is routed as part of channel routing routine otherwise, 
/// lateral flow is all routed to most downstream outflow segment . Assumes uniform time step
///
/// \param p    subbasin index
/// \param aMout_new[][] Array of mass outflows at downstream end of each segment at end of current timestep [mg/d] [size: nsegs x nConstituents]
/// \param &Options [in] Global model options information
//
void   CTransportModel::RouteMass(const int p, double **aMout_new,const optStruct &Options) const
{
  int c;//,seg;
  const double * aUnitHydro =pModel->GetSubBasin(p)->GetUnitHydrograph();
  const double * aRouteHydro=pModel->GetSubBasin(p)->GetRoutingHydrograph();
  int nMlatHist   =pModel->GetSubBasin(p)->GetLatHistorySize();
  int nMinHist    =pModel->GetSubBasin(p)->GetInflowHistorySize();
  int nSegments   =pModel->GetSubBasin(p)->GetNumSegments();
  double seg_fraction=1.0/(double)(nSegments);

  double *Mlat_new=new double [nConstituents];
  //==============================================================
  // route from catchment
  //==============================================================
  for (int c=0;c<nConstituents;c++)
  {
	  Mlat_new[c]=0.0;
    for (int n=0;n<nMlatHist;n++){
      Mlat_new[c]+=aUnitHydro[n]*aMlatHist[p][c][n];
    }
  }
  //==============================================================
  // route from channel
  //==============================================================
  if ((Options.routing==ROUTE_PLUG_FLOW) 
        || (Options.routing==ROUTE_DIFFUSIVE_WAVE))
  {	//Simple convolution
    for (c=0;c<nConstituents;c++){
		  aMout_new[nSegments-1][c]=0.0;
      for (int n=0;n<nMinHist;n++){
        aMout_new[nSegments-1][c]+=aRouteHydro[n]*aMinHist[p][c][n];
      } 
    }
  }
  //==============================================================
  else if (Options.routing==ROUTE_NONE)
  {//In channel routing instantaneous
		for (c=0;c<nConstituents;c++){
			aMout_new[nSegments-1][c]=aMinHist[p][c][0];//spits out the mass that just entered
    }
  }
  //==============================================================
  else{
    ExitGracefully("Unrecognized or unsupported constiuent routing method",STUB);
  }

	if (Options.distrib_lat_inflow==false)
  {//all fluxes from catchment are routed directly to basin outlet
    for (c=0;c<nConstituents;c++){
		  aMout_new[nSegments-1][c]+=Mlat_new[c];
    }
  }
	else{//only last segments worth
    for (c=0;c<nConstituents;c++){
		  aMout_new[nSegments-1][c]+=Mlat_new[c]*seg_fraction;
    }
	}

  return;
}
//////////////////////////////////////////////////////////////////
/// \brief Sets mass outflow from primary channel and updates flow history
/// \details Also recalculates channel storage (uglier than desired)
/// \remark Called *after* RouteMass routine is called in Solver.cpp, to reset
/// the mass outflow rates for this basin
///
/// \param **aMoutnew [in] Array of new mass outflows [mg/d] [size:  nsegments x nConstituents]
/// \param &Options [in] Global model options information
/// \param initialize Flag to indicate if flows are to only be initialized
//
void   CTransportModel::UpdateMassOutflows(const int p,  double **aMoutnew,double dummy_var,const optStruct &Options,bool initialize)
{
  double tstep=Options.timestep;

  CSubBasin *pBasin=pModel->GetSubBasin(p);

  //basin_struct *M=aBasinStructs[p][m];
  
  for (int c=0;c<nConstituents;c++)
  {
    //Update flows
    //------------------------------------------------------
	  //MoutLast[p][c]=aMout[p][c][pBasin->GetNumSegments()-1];
    //M.MoutLast=M.aMout[pBasin->GetNumSegments()-1];
	  for (int seg=0;seg<pBasin->GetNumSegments();seg++){
		  aMout[p][c][seg]=aMoutnew[seg][c];
	  }

	  if (initialize){return;}//entering initial conditions

	  //Update channel storage
	  //------------------------------------------------------
    double dt=tstep*SEC_PER_DAY;
	  double dM=0.0;
    double Mlat_new(0.0);
    for (int n=0;n<pBasin->GetLatHistorySize();n++){
		  Mlat_new+=pBasin->GetUnitHydrograph()[n]*aMlatHist[p][c][n];
	  }	
  
	  //volume change from linearly varying upstream inflow over time step 
	  dM+=0.5*(aMinHist[p][c][0]+aMinHist[p][c][1])*dt; 

	  //volume change from linearly varying downstream outflow over time step 
	  //dM-=0.5*(aMout[p][c][pBasin->GetNumSegments()-1]+MoutLast[p][c])*dt; 
	  
	  //volume change from lateral inflows
	  //dM+=0.5*(Mlat_new+MlatLast[p][c])*dt; 
  
	  //_channel_storage[p][c]+=dM;//[mg]
  
	  //Update rivulet storage
	  //------------------------------------------------------
    dM=0.0;

    //inflow from ground surface (integrated over time step)
    dM+=aMlatHist[p][c][0]*dt;

    //outflow to channel (integrated over time step)
    //dM-=0.5*(Mlat_new[p][c]+MlatLast[p][c])*dt;

    //_rivulet_storage[p][c]+=dM;//[mg]

    //MlatLast[p][c]=Mlat_new[p][c];
  }
}

double CTransportModel::GetOutflowConcentration(const int p, const int c) const
{
  double mgd=aMout[p][c][pModel->GetSubBasin(p)->GetNumSegments()-1];
  double flow=pModel->GetSubBasin(p)->GetOutflowRate();
  if (flow<=0){return 0.0;}
  double C=mgd/flow/SEC_PER_DAY/LITER_PER_M3; //[mg/L]
  return C;
}
double CTransportModel::GetIntegratedMassOutflow(const int p, const int c) const
{
  double mgd=aMout[p][c][pModel->GetSubBasin(p)->GetNumSegments()-1];

  return mgd; //[mg/d]
}