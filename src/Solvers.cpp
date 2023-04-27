/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/

#include "RavenInclude.h"
#include "Model.h"
#include "GWRiverConnection.h"

///////////////////////////////////////////////////////////////////
/// \brief Solves system of energy and mass balance ODEs/PDEs for one timestep
/// \remark This is the heart of Raven
///
/// \param *pModel [in & out] Model
/// \param &Options [in] Global model options information
/// \param &tt [in] Time step over which energy and mass balance equations are to be solved
//
void MassEnergyBalance( CModel            *pModel,
                        const optStruct   &Options,
                        const time_struct &tt)
{
  int i,j,k,p,pp,pTo,q,qs,c;                   //counters
  int NS,NB,nHRUs,nConnections=0,nProcesses;   //array sizes (local copies)
  int nConstituents;                           //
  int iSW, iAtm, iAET, iGW, iRO;               //Surface water, atmospheric precip, used PET, runoff indices

  int                iFrom          [MAX_CONNECTIONS]; //arrays used to pass values through GetRatesOfChange routines
  int                iTo            [MAX_CONNECTIONS];
  double             rates_of_change[MAX_CONNECTIONS];

  double             tstep;       //[d] timestep
  double             t;           //[d] model time
  time_struct        tt_end;      //time at end of timestep
  
  CHydroUnit        *pHRU;        //pointer to current HRU
  CSubBasin         *pBasin;      //pointer to current SubBasin
  CGroundwaterModel *pGWModel;    //pointer to GW model
  CGWRiverConnection*pGW2River;   //pointer to GW model river connection

  static double    **aPhi=NULL;   //[mm;C;mg/m2;MJ/m2] state variable arrays at initial, intermediate times;
  static double    **aPhinew;     //[mm;C;mg/m2;MJ/m2] state variable arrays at end of timestep; value after convergence
  static double    **aPhiPrevIter;

  static double     *aQinnew;     //[m3/s] inflow rate to subbasin reach p at t+dt [size=_nSubBasins]
  static double     *aQoutnew;    //[m3/s] final outflow from reach segment seg at time t+dt [size=MAX_RIVER_SEGS]
  static double     *aRouted;     //[m3]

  static double     *aMinnew;     //[mg/d] or [MJ/d] mass/energy loading of constituents to subbasin reach p at t+dt [size=_nSubBasins]
  static double     *aMoutnew;    //[mg/d] or [MJ/d] final mass/energy output from reach segment seg at time t+dt [size= MAX_RIVER_SEGS]
  static double     *aRoutedMass; //[mg/d] or [MJ/d] amount of mass/energy [size= _nSubBasins]

  static double    **rate_guess;  //need to set first array to nProcesses

  static int        *kFrom;
  static int        *kTo; 
  static double     *exchange_rates=NULL; 
  
  //local shorthand for often-used variables
  NS           =pModel->GetNumStateVars();
  NB           =pModel->GetNumSubBasins();
  nHRUs        =pModel->GetNumHRUs();
  nProcesses   =pModel->GetNumProcesses();
  nConstituents=pModel->GetTransportModel()->GetNumConstituents();
  tstep        =Options.timestep;
  t            =tt.model_time;

  JulianConvert(t+tstep,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt_end);

  //Reserve static memory ===========================================
  //(only gets called once in course of simulation)
  if (aPhi==NULL)
  {

    aPhi        =new double *[nHRUs];
    aPhinew     =new double *[nHRUs];
    aPhiPrevIter=new double *[nHRUs];

    for (k=0;k<nHRUs;k++)
    {
      aPhi[k]        =new double [NS];
      aPhinew[k]     =new double [NS];
      aPhiPrevIter[k]=new double [NS];
    }

    aQoutnew    =NULL;
    aQinnew     =new double [NB];
    aRouted     =new double [NB];
    aQoutnew    =new double [MAX_RIVER_SEGS];
    ExitGracefullyIf(aQoutnew==NULL,"MassEnergyBalance",OUT_OF_MEMORY);

    aMinnew     =NULL;
    aMoutnew    =NULL;
    aRoutedMass =NULL;
    if(nConstituents>0)
    {
      aMinnew     =new double [NB];
      aRoutedMass =new double [NB];
      aMoutnew    =new double[MAX_RIVER_SEGS];
      ExitGracefullyIf(aMoutnew==NULL,"MassEnergyBalance(2)",OUT_OF_MEMORY);
      for(p=0;p<NB;p++) {
        aMinnew    [p] =0;
        aRoutedMass[p] =0;
      }
      for(i=0;i<MAX_RIVER_SEGS;i++) {
        aMoutnew   [i]=0.0;
      }
    }

    if(Options.sol_method==ITERATED_HEUN)
    {
      rate_guess = new double *[nProcesses];    //need to set first array to numProcesses
      for (j=0;j<nProcesses;j++){
        rate_guess[j]=new double [NS*NS];       //maximum number of connections possible
      }
    }
    //For lateral flow processes
    kFrom         =new int   [MAX_LAT_CONNECTIONS];
    kTo           =new int   [MAX_LAT_CONNECTIONS];
    exchange_rates=new double[MAX_LAT_CONNECTIONS];
  }//end static memory if

  if(Options.modeltype == MODELTYPE_COUPLED)
  {
    // Get pointer to GW model
    pGWModel = pModel->GetGroundwaterModel();
    pGW2River = pGWModel->GetRiverConnection();
  }
  // Initialize variables============================================
  for (i=0;i<MAX_CONNECTIONS;i++)
  {
    iFrom          [i]=DOESNT_EXIST;
    iTo            [i]=DOESNT_EXIST;
    rates_of_change[i]=0.0;
  }
  for (k=0;k<nHRUs;k++)
  {
    pHRU=pModel->GetHydroUnit(k);
    for (i=0;i<NS ;i++)
    {
      aPhi        [k][i]=pHRU->GetStateVarValue(i);
      aPhinew     [k][i]=aPhi[k][i];
      aPhiPrevIter[k][i]=aPhi[k][i];
    }
  }

  iSW  =pModel->GetStateVarIndex(SURFACE_WATER);
  iAtm =pModel->GetStateVarIndex(ATMOS_PRECIP);

  // Used PET and runoff reboots to zero every timestep==============
  iAET=pModel->GetStateVarIndex(AET);
  iRO =pModel->GetStateVarIndex(RUNOFF);
  if(iAET!=DOESNT_EXIST) {
    for(k=0;k<nHRUs;k++)
    {
      aPhi        [k][iAET]=0.0;
      aPhinew     [k][iAET]=0.0;
      aPhiPrevIter[k][iAET]=0.0;
    }
  }
  if(iRO!=DOESNT_EXIST) {
    for(k=0;k<nHRUs;k++)
    {
      aPhi        [k][iRO ]=0.0;
      aPhinew     [k][iRO ]=0.0;
      aPhiPrevIter[k][iRO ]=0.0;
    }
  }

  // GW reboots each time step=======================================
  iGW=pModel->GetStateVarIndex(GROUNDWATER);
  if(iGW!=DOESNT_EXIST) {
    for(k=0;k<nHRUs;k++)
    {
      aPhi        [k][iGW]=0.0;
      aPhinew     [k][iGW]=0.0;
      aPhiPrevIter[k][iGW]=0.0;
    }
  }

  //=================================================================
  //==Standard (in series) approach==================================
  // -order is critical!
  if (Options.sol_method==ORDERED_SERIES)
  {
    for (k=0;k<nHRUs;k++)
    {
      pHRU=pModel->GetHydroUnit(k);
      pModel->ApplyLocalParamOverrrides(k, false); 

      if(pHRU->IsEnabled()) 
      {
        qs=0;
        for(j=0;j<nProcesses;j++)
        {
          nConnections=0;
          if(pModel->ApplyProcess(j,aPhinew[k],pHRU,Options,tt,iFrom,iTo,nConnections,rates_of_change)) //note aPhinew is newest state variable vector
          {
#ifdef _STRICTCHECK_
            if(nConnections>MAX_CONNECTIONS) {
              cout<<nConnections<<endl;
              ExitGracefully("MassEnergyBalance:: Maximum number of connections exceeded. Please contact author.",RUNTIME_ERR); }
#endif
            for(q=0;q<nConnections;q++)//each process may have multiple connections
            {
              sv_type typ=pModel->GetStateVarType(iFrom[q]);
              if(iTo[q]!=iFrom[q]) {
                aPhinew[k][iFrom[q]]-=rates_of_change[q]*tstep;//mass/energy balance maintained
                aPhinew[k][iTo  [q]]+=rates_of_change[q]*tstep;//change is an exchange of energy or mass, which must be preserved
              }
              else if (CStateVariable::IsWaterStorage(typ) && (typ!=CONVOLUTION)){
                rates_of_change[q]=0.0;
                aPhinew[k][iTo  [q]]+=0.0; //likely from redirect - water moves back to itself
              }
              else {
                aPhinew[k][iTo  [q]]+=rates_of_change[q]*tstep;//for state vars that are not storage compartments
              }
              pModel->IncrementBalance(qs,k,rates_of_change[q]*tstep);   //this is only this easy for Euler/Ordered!
              qs++;
            }//end for q=0 to nConnections
          }// end if (pModel->ApplyProcess
          else
          {
            for(q=0;q<nConnections;q++)
            {
              pModel->IncrementBalance(qs,k,0.0);
              qs++;
            }
          }
        }//end for j=0 to nProcesses
      }
      pModel->ApplyLocalParamOverrrides(k, true); 
    }//end for k=0 to nHRUs

  }//end if Options.sol_method==ORDERED_SERIES

  //=================================================================
  //==Simple Euler Method ===========================================
  // -order of processes doesn't matter
  else if (Options.sol_method==EULER)
  {
    for (k=0;k<nHRUs;k++)
    {
      pHRU=pModel->GetHydroUnit(k);

      //model all hydrologic processes occuring at HRU scale
      //-----------------------------------------------------------------
      qs=0;
      for (j=0;j<nProcesses;j++)
      {
        nConnections=0;

        if (pModel->ApplyProcess(j,aPhi   [k],pHRU,Options,tt,iFrom,iTo,nConnections,rates_of_change))//note aPhi is info from start of timestep
        {
#ifdef _STRICTCHECK_
          if(nConnections>MAX_CONNECTIONS) {
            cout<<nConnections<<endl;
            ExitGracefully("MassEnergyBalance:: Maximum number of connections exceeded. Please contact author.",RUNTIME_ERR);}
#endif
          for (q=0;q<nConnections;q++)//each process may have multiple connections
          {
            sv_type typ=pModel->GetStateVarType(iFrom[q]);
            if (iTo[q]!=iFrom[q]){
              aPhinew[k][iFrom[q]]-=rates_of_change[q]*tstep;//mass/energy balance maintained
              aPhinew[k][iTo  [q]]+=rates_of_change[q]*tstep;//change is an exchange of energy or mass, which must be preserved
            }
            else if (CStateVariable::IsWaterStorage(typ) && (typ!=CONVOLUTION)){
              rates_of_change[q]=0.0;
              aPhinew[k][iTo  [q]]+=0.0; //likely from redirect - water moves back to itself
            }
            else{
              aPhinew[k][iTo  [q]]+=rates_of_change[q]*tstep;//for state vars that are not storage compartments
            }
            pModel->IncrementBalance(qs,k,rates_of_change[q]*tstep);//this is only this easy for Euler/Ordered!
            qs++;
          }//end for q=0 to nConnections
        }
        else
        {
          for(q=0;q<nConnections;q++)
          {
            pModel->IncrementBalance(qs,k,0.0);
            qs++;
          }
        }
      }//end for j=0 to nProcesses
    }//end for k=0 to nHRUs
  }//end if Options.sol_method==EULER

  //===================================================================
  //==Iterated Heun Method ============================================
  // -order of processes doesn't matter, converges to specified criteria
  else if(Options.sol_method==ITERATED_HEUN)
  {
    qs= 0;
    int    iter = 0;              //iteration counter
    bool   converg = false;
    double converg_check = 0.0;
    double rate1[MAX_CONNECTIONS];
    double rate2[MAX_CONNECTIONS];

    //Go through all HRUs
    for (k=0;k<nHRUs;k++)
    {
      pHRU=pModel->GetHydroUnit(k);

      do  //Iterate
      {
        iter++;           // iteration counter

        for(i=0;i<NS;i++)  //loop through all state variables
        {
          aPhiPrevIter[k][i]=aPhinew[k][i];       //iteration k-1 value
          aPhinew     [k][i]=aPhi   [k][i];
        }

        //model all other hydrologic processes occuring at HRU scale
        //-----------------------------------------------------------------
        for (j=0;j<nProcesses;j++)
        {
          // ROC 1 - uses initial state var values
          // ROC 2 - uses previous iteration values
          if (pModel->ApplyProcess(j,aPhi[k]        ,pHRU,Options,tt     ,iFrom,iTo,nConnections,rate1))
          {
            pModel->ApplyProcess(j,aPhiPrevIter[k],pHRU,Options,tt_end ,iFrom,iTo,nConnections,rate2);

            if(nConnections>MAX_CONNECTIONS) {
              cout<<nConnections<<endl;
              ExitGracefully("MassEnergyBalance:: Maximum number of connections exceeded. Please contact author.",RUNTIME_ERR);
            }

            for (q=0;q<nConnections;q++)//each process may have multiple connections
            {
              sv_type typ=pModel->GetStateVarType(iFrom[q]);
              rate_guess[j][q] = 0.5*(rate1[q] + rate2[q]);

              if(iFrom[q]==iAtm){               //check if water is coming from precipitation
                rate_guess[j][q] = rate1[q];    //sets the rate of change to be the original (prevents over filling of SV's)
              }

              if (iTo[q]!=iFrom[q]){
                aPhinew[k][iFrom[q]]  -= rate_guess[j][q]*tstep;//mass/energy balance maintained
                aPhinew[k][iTo  [q]]  += rate_guess[j][q]*tstep;//change is an exchange of energy or mass, which must be preserved
              }
              else if (CStateVariable::IsWaterStorage(typ) && (typ!=CONVOLUTION)){
                rates_of_change[q]=0.0;
                aPhinew[k][iTo  [q]]+=0.0; //likely from redirect - water moves back to itself
              }
              else{   //correction for state vars that are not storage compartments
                aPhinew[k][iTo  [q]]  += rate_guess[j][q]*tstep;
              }
            }//end for q=0 to nConnections
          }
        }//end for j=0 to nProcesses

        //Calculate convegence criterion
        for(i=0;i<NS;i++)
        {
          //converg_check += (2*(fabs(aPhinew[k][i] - aPhiPrevIter[i])/(aPhinew[k][i] + aPhiPrevIter[i]))); //possible converg check #1
          converg_check += (fabs(aPhinew[k][i] - aPhiPrevIter[k][i])/NS);    //possible converg check #2
          //converg_check += (fabs(aPhinew[k][i]-aPhiPrevIter[i])/aPhinew[k][i]); //possible converg check #3
          //converg_check = pow((converg_check + pow((aPhinew[k][i]-aPhiPrecIter[i]),2)),0.5);  //possible converg check #4
        }
        converg = false;

        if((converg_check <= Options.convergence_crit) ||
           (iter          == Options.max_iterations))   //convergence check
        {
          qs  =0;
          iter=0;
          converg = true;

          for(j=0;j<nProcesses;j++)
          {
            nConnections=pModel->GetNumConnections(j);
            for(q=0;q<nConnections;q++)
            {
              pModel->IncrementBalance(qs,k,rate_guess[j][q]*tstep);
              qs++;
            }
          }
        }//end of (converg_check <=...)

        converg_check = 0.0;

      } while(converg != true);  //end do loop
    }//end of for k=0 to nHRUs
  }//end iterated Heun

  //==Other Methods Below ============================================
  else
  {
    ExitGracefully("MassEnergyBalance",STUB);
  }

  //-----------------------------------------------------------------
  //      LATERAL EXCHANGE PROCESSES
  //-----------------------------------------------------------------
  int    qss=0;//index of lateral process connection
  int    nLatConnections;
  double Afrom,Ato;
   
  for (q=0;q<MAX_LAT_CONNECTIONS;q++)
  {
    kFrom[q]=DOESNT_EXIST;
    kTo  [q]=DOESNT_EXIST;
    exchange_rates[i]=0.0;
  }
  for (j=0;j<nProcesses;j++)
  {
    if (pModel->ApplyLateralProcess(j,aPhinew,Options,tt,kFrom,kTo,iFrom,iTo,nLatConnections,exchange_rates))
    {       
#ifdef _STRICTCHECK_
      if(nLatConnections>MAX_LAT_CONNECTIONS) {
        cout<<nLatConnections<<endl;
        ExitGracefully("MassEnergyBalance:: Maximum number of lateral connections exceeded. Please contact author.",RUNTIME_ERR);
      }
#endif 
      for (q=0;q<nLatConnections;q++)
      {
        Afrom=pModel->GetHydroUnit(kFrom[q])->GetArea();
        Ato  =pModel->GetHydroUnit(kTo[q]  )->GetArea();
        aPhinew[kFrom[q]][iFrom[q]]-=exchange_rates[q]/Afrom*tstep;  
        aPhinew[  kTo[q]][  iTo[q]]+=exchange_rates[q]/Ato  *tstep;
       
        pModel->IncrementLatBalance(qss,exchange_rates[q]*tstep); 

        qss++;
      }
    }
    else{
      for(q=0;q<nLatConnections;q++)
      {
        pModel->IncrementLatBalance(qss,0.0);
        qss++;
      }
    }
  }

  //-----------------------------------------------------------------
  //      GROUNDWATER SOLVER
  //-----------------------------------------------------------------
  // Following solution to SW system at end of timestep, solve GW system for lateral flow
  if (Options.modeltype == MODELTYPE_COUPLED)
	{
    // Update River water levels
    if (pGW2River->GetNumSegments() > 0){
      for (p = 0; p < NB; p++)
      {
        pBasin = pModel->GetSubBasin(p);
        pGW2River->UpdateRiverLevelsBySB(p, pBasin);
      }
    }
    // Add in Raven Flux for each HRU (non-GWSW Process contribution)
    for (k=0;k<nHRUs;k++) {
      pHRU=pModel->GetHydroUnit(k);
      pGWModel->FluxToGWEquation(pHRU, aPhinew[k][iGW]);
    }
    
    pGWModel->Solve(t);               // Run MODFLOW-USG
    pGWModel->PostSolve(tstep);       // Post-solve MFUSG Routines, Budget update, etc
    pGW2River->UpdateRiverFlux();     // Update River Fluxes (for routing)
    pGWModel->ClearMatrix();
  } // End of Groundwater processes


  //-----------------------------------------------------------------
  //      ROUTING
  //-----------------------------------------------------------------
  double res_ht,res_outflow;
  double down_Q,irr_Q,div_Q, Qwithdrawn, SWvol;
  int    pDivert;
  res_constraint res_const;
  const int MAX_CONTROL_STRUCTURES=10;
  double *res_Qstruct=new double [MAX_CONTROL_STRUCTURES];
  //determine total outflow from HRUs into respective basins (aRouted[p])
  for (p=0;p<NB;p++)
  {
    aRouted[p]=0.0;
    aQinnew[p]=0.0;
  }
  for (k=0;k<nHRUs;k++)
  {
    pHRU=pModel->GetHydroUnit(k);
    if(pHRU->IsEnabled())
    {
      p   =pHRU->GetSubBasinIndex();
      SWvol=(aPhinew[k][iSW]/MM_PER_METER)*(pHRU->GetArea()*M2_PER_KM2);//[m3] 
      if(pHRU->IsLinkedToReservoir()) { 
        pModel->GetSubBasin(p)->GetReservoir()->SetPrecip(SWvol);//[SW is treated as precip on reservoir]
      }
      else{
        //surface water moved instantaneously from HRU to basin reach/channel storage
        aRouted[p]+=SWvol;
      }
      aPhinew[k][iRO]=aPhinew[k][iSW]; //track net runoff [mm]
      aPhinew[k][iSW]=0.0;             //zero out surface water storage
    }
  }
  // Identify magnitude of flow diversions, calculate inflows  
  for(p=0;p<NB;p++)
  {
    //Regular Diversions
    pBasin=pModel->GetSubBasin(p);
    for(int i=0; i<pBasin->GetNumDiversions();i++) {
      div_Q=pBasin->GetDiversionFlow(i,pBasin->GetOutflowRate(),Options,tt,pDivert); //diversions based upon flows at start of timestep
      if(pDivert!=DOESNT_EXIST)
      {
        aQinnew[pDivert]+=div_Q;
      }
      else {
        // \todo[funct] - must handle this in mass balance - should go to cumulative output
	      //Qdivloss+=div_Q;
      }
    }
    //Diversions from reservoir control structures
    CReservoir *pRes=pBasin->GetReservoir();
    if (pRes != NULL) {
      for (int i = 0; i < pRes->GetNumControlStructures(); i++) {
        if (pRes->GetControlFlowTarget(i) != pBasin->GetDownstreamID()) {
          pDivert=pModel->GetSubBasinIndex(pRes->GetControlFlowTarget(i)); //p, not SBID
          if (pDivert!=DOESNT_EXIST){
            aQinnew[pDivert]+=pRes->GetControlOutflow(i);
          }
          else {
            //\todo[funct] - must handle this in mass balance - should go to cumulative output
	          //Qdivloss+=pRes->GetControlOutflow(i);
          }
        }
      }
    }
    //User-specified inflows to upstream end 
    aQinnew[p]+=pBasin->GetSpecifiedInflow(t+tstep);
  }
  // Route water over timestep
  // ----------------------------------------------------------------------------------------
  // calculations performed in order from upstream (pp=0) to downstream (pp=nSubBasins-1)
  for (pp=0;pp<NB;pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp); //p refers to actual index of basin, pp is ordered list index upstream to down

    pBasin=pModel->GetSubBasin(p);
    if(pBasin->IsEnabled())
    {
      pBasin->UpdateSubBasin(tt,Options);            // also used to assimilate lake levels and update routing hydrograph for timestep
      
      pBasin->UpdateInflow(aQinnew[p]);              // from upstream, diversions, and specified flows

      down_Q=pBasin->GetDownstreamInflow(t);         // treated as additional runoff (period starting)

      if (Options.modeltype == MODELTYPE_COUPLED)
	    {
        aRouted[p]+= pGW2River->CalcRiverFlowBySB(p)*tstep;      // [m3]
      }

      pBasin->UpdateLateralInflow(aRouted[p]/(tstep*SEC_PER_DAY)+down_Q);//[m3/d]->[m3/s]

      pBasin->RouteWater    (aQoutnew,res_ht,res_outflow,res_const,res_Qstruct,Options,tt);      //Where everything happens!

      Qwithdrawn=0;
      irr_Q=pBasin->ApplyIrrigationDemand(t+tstep,aQoutnew[pBasin->GetNumSegments()-1]);
      Qwithdrawn+=irr_Q;

      for(int i=0; i<pBasin->GetNumDiversions();i++) {
        div_Q=pBasin->GetDiversionFlow(i,pBasin->GetOutflowRate(),Options,tt,pDivert); //diversions based upon flows at start of timestep
        Qwithdrawn+=div_Q;
      }

      pBasin->UpdateOutflows(aQoutnew,irr_Q,res_ht,res_outflow,res_const,res_Qstruct,Options,tt,false);//actually updates flow values here

      pModel->AssimilationOverride(p,Options,tt); //modifies flows using assimilation, if needed

      pTo   =pModel->GetDownstreamBasin(p);
      if(pTo!=DOESNT_EXIST)//update downstream inflows
      {
        aQinnew[pTo]+=pBasin->GetOutflowRate()-Qwithdrawn;
      }
      else {
        //still need to remove Qwithdrawn from somewhere if downstream outflow doesn't exist!
        //Qdivloss+=Qwithdrawn;
      }

      if(pBasin->GetReservoir()!=NULL) {//update AET for reservoir-linked HRUs
        k=pBasin->GetReservoir()->GetHRUIndex();
        if ((k!=DOESNT_EXIST) && (iAET!=DOESNT_EXIST)){
          aPhinew[k][iAET]=pBasin->GetReservoir()->GetAET();//[mm/d]
        }
      }
    }
  }//end for pp...
  delete [] res_Qstruct;

  //-----------------------------------------------------------------
  //      CONSTITUENT (MASS OR ENERGY) ROUTING
  //-----------------------------------------------------------------
   //determine total mass/energy loading from HRUs into respective basins (aRoutedMass[p])
  double ResMass    =0;
  double ResSedMass =0;
  double MassOutflow=0;
  double Ploading;
  int    iSWmass,m;
  CConstituentModel *pConstitModel;
  for(c=0;c<nConstituents;c++)
  {
    pConstitModel=pModel->GetTransportModel()->GetConstituentModel(c);
    for(p=0;p<NB;p++) {
      aRoutedMass[p]=0.0;
      aMinnew    [p]=0.0;
    }
    for(k=0;k<nHRUs;k++)
    {
      pHRU=pModel->GetHydroUnit(k);
      if(pHRU->IsEnabled())
      {
        p       =pHRU->GetSubBasinIndex();
        m       =pModel->GetTransportModel()->GetLayerIndex(c,iSW);
        iSWmass =pModel->GetStateVarIndex(CONSTITUENT,m);
        
        Ploading=(aPhinew[k][iSWmass])*(pHRU->GetArea()*M2_PER_KM2)/tstep;//[mg/d] or [MJ/d] [iSWmass is ALL precip mass/energy on HRU]
        if(pHRU->IsLinkedToReservoir()) {
          pConstitModel->SetReservoirPrecipLoad(p,Ploading);
        }
        else{
          aRoutedMass[p]+=Ploading;
        }
        aPhinew[k][iSWmass]=0.0; //empty out from landscape storage
      }
    }
    
    //Route mass over timestep
    //calculations performed in order from upstream (pp=0) to downstream (pp=nSubBasins-1)
    for(pp=0;pp<NB;pp++)
    {
      p     =pModel->GetOrderedSubBasinIndex(pp); //p refers to actual index of basin, pp is ordered list
      pBasin=pModel->GetSubBasin(p);

      if(pBasin->IsEnabled())
      {
        pConstitModel->ApplySpecifiedMassInflows(p,t+tstep,aMinnew[p]); //overrides or supplements mass loadings

        pConstitModel->SetMassInflows    (p,aMinnew[p]);
        pConstitModel->SetLateralInfluxes(p,aRoutedMass[p]);
        pConstitModel->RouteMass         (p,aMoutnew,ResMass,ResSedMass,Options,tt);  //Where everything happens!
        pConstitModel->UpdateMassOutflows(p,aMoutnew,ResMass,ResSedMass,MassOutflow,Options,tt,false); //actually updates mass flow values here

        pTo   =pModel->GetDownstreamBasin(p);
        if(pTo!=DOESNT_EXIST)
        {
          aMinnew[pTo]+=MassOutflow;
        }
      }
    }//end for pp...
  }//end (c=0;c<nConstituents;c++)

  //update state variable values=====================================
  for (k=0;k<nHRUs;k++){
    pHRU=pModel->GetHydroUnit(k);
    if(pHRU->IsEnabled())
    {
      for(i=0;i<NS;i++){
        pHRU->SetStateVarValue(i,aPhinew[k][i]);
      }
    }
  }

  //delete static arrays (only called once)=========================
  if(t>=Options.duration-Options.timestep)
  {
    if(DESTRUCTOR_DEBUG) { cout<<"DELETING STATIC ARRAYS IN MASSENERGYBALANCE"<<endl; }
    for(k=0;k<nHRUs;k++) { delete[] aPhi[k];         } delete[] aPhi;         aPhi=NULL;
    for(k=0;k<nHRUs;k++) { delete[] aPhinew[k];      } delete[] aPhinew;      aPhinew=NULL;
    for(k=0;k<nHRUs;k++) { delete[] aPhiPrevIter[k]; } delete[] aPhiPrevIter; aPhiPrevIter=NULL;
    if(Options.sol_method == ITERATED_HEUN)
    {
      for(j=0;j<nProcesses;j++) { delete[] rate_guess[j]; }  delete[] rate_guess; rate_guess=NULL;
    }
    delete[] aQinnew;      aQinnew     = NULL;
    delete[] aQoutnew;     aQoutnew    = NULL;
    delete[] aRouted;      aRouted     = NULL;
    //delete transport static arrays.
    if(nConstituents>0)
    {
      delete[] aMinnew;
      delete[] aRoutedMass;
      delete[] aMoutnew;
    }
  }
}

