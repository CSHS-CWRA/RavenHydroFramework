/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------*/

#include "RavenInclude.h"
#include "Model.h"

//Defined in GWSolvers.cpp
void SolveGWSystem    (      CModel      *pModel, 
                       const optStruct   &Options,
					             const time_struct &tt,             //time
                             double      *GW_head_prev,
                             double      *GW_head_new);         

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
  int         i,j,k,p,pp,pTo,q,qs,c;                 //counters
  int         NS,NB,nHRUs,nConnections=0,nProcesses; //array sizes (local copies)
  int         nConstituents;                         //
  int         iSW, iAtm, iAET, iGW;                  //Surface water, atmospheric precip, used PET indices

  int         iFrom[MAX_CONNECTIONS];                //arrays used to pass values through GetRatesOfChange routines
  int         iTo  [MAX_CONNECTIONS];
  double      rates_of_change[MAX_CONNECTIONS];

  double      tstep;              //[day] timestep
  double      t;                  //model time
  time_struct tt_end;             //time at end of timestep
  
  double      stor_diff;
  double      poro;               //aquifer porosity
  CHydroUnit *pHRU;               //pointer to current HRU
  CSubBasin  *pBasin;             //pointer to current SubBasin

  static double    **aPhi=NULL;   //[mm;C] state variable arrays at initial, intermediate times;
  static double    **aPhinew;     //[mm;C] state variable arrays at end of timestep; value after convergence
  static double    **aPhiPrevIter;

  static double     *aQinnew;     //[m3/s] inflow rate to subbasin reach p at t+dt [size=_nSubBasins]
  static double     *aQoutnew;    //[m3/s] final outflow from reach segment seg at time t+dt [size=MAX_RIVER_SEGS]
  static double     *aRouted;     //[m3]

  static double    **aMinnew;     //[mg/d] mass loading of constituents to subbasin reach p at t+dt [size=_nSubBasins x nConstituents ]
  static double    **aMoutnew;    //[mg/d] final mass output from reach segment seg at time t+dt [size= MAX_RIVER_SEGS  x nConstituents ]
  static double    **aRoutedMass; //[mg/d] amount of mass [size= _nSubBasins xnConstituents]
  static double     *aResMass;    //[mg]   amount of mass in reservoir [size= nConstitutents]
  static double     *aMassOutflow;//[mg/d] rate of mass outflow from subbasin at time t+dt [size=nConstituents]

  static double    **rate_guess;  //need to set first array to nProcesses

  double             *GW_stor     =NULL; //[mm] groundwater storage at start of timestep
  double             *GW_stor_new =NULL; //[mm] groundwater storage after processes have moved water
  double             *GW_stor_old =NULL;
  double             *GW_heads    =NULL; //[m] groundwater heads at start of timestep
  double             *GW_heads_new=NULL; //[m] groundwater heads after processes have moved water
  double             *head_change =NULL;

  
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

    aMoutnew    =NULL;
    aRoutedMass =NULL;
    aResMass    =NULL;
    aMassOutflow=NULL;
    if (nConstituents>0)
    {
      aMinnew     =new double *[NB];
      aRoutedMass =new double *[NB];
      for (p=0;p<NB;p++){
        aRoutedMass[p]=NULL;
        aMinnew    [p] =new double [nConstituents];
        aRoutedMass[p] =new double [nConstituents];
        ExitGracefullyIf(aRoutedMass[p]==NULL,"MassEnergyBalance(2)",OUT_OF_MEMORY);
      }

      aMoutnew    =new double *[MAX_RIVER_SEGS];
      ExitGracefullyIf(aMoutnew==NULL,"MassEnergyBalance(3)",OUT_OF_MEMORY);
      for (i=0;i<MAX_RIVER_SEGS;i++){
        aMoutnew[i]=NULL;
        aMoutnew[i]=new double [nConstituents];
        ExitGracefullyIf(aMoutnew[i]==NULL,"MassEnergyBalance(4)",OUT_OF_MEMORY);
      }
      aResMass=new double [nConstituents];
      aMassOutflow=new double [nConstituents];
    }

    if(Options.sol_method==ITERATED_HEUN)
    {
      rate_guess = new double *[nProcesses];      //need to set first array to numProcesses
      for (j=0;j<nProcesses;j++){
        rate_guess[j]=new double [NS*NS]; //maximum number of connections possible
      }
    }
  }

  if((Options.modeltype == MODELTYPE_COUPLED) || (Options.modeltype == MODELTYPE_GROUNDWATER))
  {
    //GW storage pre-process
    int nNodes=(int)(pModel->GetGroundwaterModel()->GetGWGeom()->GetGWGeoProperty("NUM_NODES"));
    GW_stor         = new double [nNodes];
    GW_stor_new     = new double [nNodes];
    GW_stor_old     = new double [nNodes];
    GW_heads        = new double [nNodes];
    GW_heads_new    = new double [nNodes];
    head_change     = new double [nNodes];
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

  // Used PET reboots to zero every timestep=========================
  iAET=pModel->GetStateVarIndex(AET);
  if(iAET!=DOESNT_EXIST) {
    for(k=0;k<nHRUs;k++)
    {
      aPhi        [k][iAET]=0.0;
      aPhinew     [k][iAET]=0.0;
      aPhiPrevIter[k][iAET]=0.0;
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

      if(pHRU->IsEnabled()) {
        //model all hydrologic processes occuring at HRU scale
        //-----------------------------------------------------------------
        qs=0;
        for(j=0;j<nProcesses;j++)
        {
          nConnections=0;

          if(pModel->ApplyProcess(j,aPhinew[k],pHRU,Options,tt,iFrom,iTo,nConnections,rates_of_change)) //note aPhinew is newest version
          {
            if(nConnections>MAX_CONNECTIONS) {
              cout<<nConnections<<endl;
              ExitGracefully("MassEnergyBalance:: Maximum number of connections exceeded. Please contact author.",RUNTIME_ERR);
            }

            for(q=0;q<nConnections;q++)//each process may have multiple connections
            {
              if(iTo[q]!=iFrom[q]) {
                aPhinew[k][iFrom[q]]-=rates_of_change[q]*tstep;//mass/energy balance maintained
                aPhinew[k][iTo[q]]+=rates_of_change[q]*tstep;//change is an exchange of energy or mass, which must be preserved
              }
              else {
                aPhinew[k][iTo[q]]+=rates_of_change[q]*tstep;//for state vars that are not storage compartments
              }
              pModel->IncrementBalance(qs,k,rates_of_change[q]*tstep);//this is only this easy for Euler/Ordered!
              qs++;
            }//end for q=0 to nConnections
          }// end if (pModel->ApplyProcess
          else
          {
            qs+=nConnections;
          }
        }//end for j=0 to nProcesses
      }
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
          if(nConnections>MAX_CONNECTIONS) {
            cout<<nConnections<<endl;
            ExitGracefully("MassEnergyBalance:: Maximum number of connections exceeded. Please contact author.",RUNTIME_ERR);
          }

          for (q=0;q<nConnections;q++)//each process may have multiple connections
          {
            if (iTo[q]!=iFrom[q]){
              aPhinew[k][iFrom[q]]-=rates_of_change[q]*tstep;//mass/energy balance maintained
              aPhinew[k][iTo  [q]]+=rates_of_change[q]*tstep;//change is an exchange of energy or mass, which must be preserved
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
          qs+=nConnections;
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
    double rate1[MAX_STATE_VARS];
    double rate2[MAX_STATE_VARS];

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
              rate_guess[j][q] = 0.5*(rate1[q] + rate2[q]);

              if(iFrom[q]==iAtm){               //check if water is coming from precipitation
                rate_guess[j][q] = rate1[q];    //sets the rate of change to be the original (prevents over filling of SV's)
              }

              if (iTo[q]!=iFrom[q]){
                aPhinew[k][iFrom[q]]  -= rate_guess[j][q]*tstep;//mass/energy balance maintained
                aPhinew[k][iTo  [q]]  += rate_guess[j][q]*tstep;//change is an exchange of energy or mass, which must be preserved
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
  int    *kFrom,*kTo;
  
  double *exchange_rates=NULL;
  kFrom         =new int   [MAX_LAT_CONNECTIONS];
  kTo           =new int   [MAX_LAT_CONNECTIONS];
  exchange_rates=new double[MAX_LAT_CONNECTIONS];
  ExitGracefullyIf(exchange_rates==NULL,"MassEnergyBalance(5)",OUT_OF_MEMORY);
   
  for (int q=0;q<MAX_LAT_CONNECTIONS;q++)
  {
    kFrom[q]=DOESNT_EXIST;
    kTo  [q]=DOESNT_EXIST;
    exchange_rates[i]=0.0;
  }
  for (j=0;j<nProcesses;j++)
  {
    if (pModel->ApplyLateralProcess(j,aPhinew,Options,tt,kFrom,kTo,iFrom,iTo,nLatConnections,exchange_rates))
    {       
      if(nLatConnections>MAX_LAT_CONNECTIONS) {
        cout<<nLatConnections<<endl;
        ExitGracefully("MassEnergyBalance:: Maximum number of lateral connections exceeded. Please contact author.",RUNTIME_ERR);
      }

      for (int q=0;q<nLatConnections;q++)
      {
#ifdef _STRICTCHECK_
        ExitGracefullyIf(kFrom[q]==DOESNT_EXIST,"MassEnergyBalance: kFrom[q] doesnt exist",RUNTIME_ERR);
        ExitGracefullyIf(kTo  [q]==DOESNT_EXIST,"MassEnergyBalance: kFrom[q] doesnt exist",RUNTIME_ERR);
#endif 
        Afrom=pModel->GetHydroUnit(kFrom[q])->GetArea();
        Ato  =pModel->GetHydroUnit(kTo[q]  )->GetArea();
        aPhinew[kFrom[q]][iFrom[q]]-=exchange_rates[q]/Afrom*tstep;  
        aPhinew[  kTo[q]][  iTo[q]]+=exchange_rates[q]/Ato  *tstep;
       
        pModel->IncrementLatBalance(qss,exchange_rates[q]*tstep); 

        qss++;
      }
    }
    else{
      qss+=nLatConnections;
    }
  }
  delete [] kFrom;
  delete [] kTo;
  delete [] exchange_rates;

  //-----------------------------------------------------------------
  //      GROUNDWATER SOLVER
  //-----------------------------------------------------------------
  //GW MIGRATE - FULL REVIEW NEEDED
  //Following solution to SW system at end of timestep, solve GW system for lateral flow
  if (Options.modeltype == MODELTYPE_GROUNDWATER || Options.modeltype == MODELTYPE_COUPLED)
	{
    CGroundwaterModel *pGWmod=pModel->GetGroundwaterModel();
    COverlapExchangeClass *pOE;
    string param;
    double bot;
    double stor_sum = 0;
    double stor_sum2 = 0;
    double weight;

    for(k = 0; k < nHRUs; k++)
    {
      head_change[k] = 0;
    }

    //convert storage to Heads
    for (k=0;k<nHRUs;k++)
    {  
      pHRU            = pModel->GetHydroUnit(k);                           
      poro            = pHRU->GetAquiferProps(0)->porosity;
      bot             = pGWmod->GetGWGeom()->GetGWGeoProperty("BOT_ELEV",k,0);
      iGW             = pModel->GetStateVarIndex(GROUNDWATER,0);                 //*****FIX - needs to work with multiple layers
      
      GW_stor[k]      = aPhinew[k][iGW];                                         //post process application storage value [mm]
      GW_stor_old[k]  = aPhi   [k][iGW];
      stor_diff       = GW_stor[k] - GW_stor_old[k];
      //cout<<"stor: "<<GW_stor[k]<<"   stor_old: "<<GW_stor_old[k]<<endl;
      if(Options.overlap_type == AREA_WEIGHTED)
      {
        pOE         = pModel->GetGroundwaterModel()->GetOEs(0);
        GW_heads[k] = pGWmod->GetGWGeom()->GetGWGeoProperty("GW_HEAD",k,0);       //*****FIX - needs to work with multiple layers

        //get OE table
        for(i = 0; i < pOE->GetOEProperty("NGWCELLS"); i++)
        {
          weight = pOE->GetOEStruct()->overlap[k][i];
          head_change[i] += (weight * stor_diff);
          //cout<<"w: "<<weight<<"   sd: "<<stor_diff<<endl;
        }
      }
      else
      {
        GW_heads[k] = pGWmod->GetGWGeom()->GetGWGeoProperty("GW_HEAD",k,0);
        GW_heads_new[k] = ((GW_stor[k] / poro) / MM_PER_METER) + bot;         //convert to head [m]
      }
      stor_sum += GW_stor[k];
    }
    if(Options.overlap_type == AREA_WEIGHTED)
    {
      for(i = 0; i<nHRUs; i++)
      {
        GW_heads_new[i] = GW_heads[i] + ((head_change[i] / poro) / MM_PER_METER);
        //cout<<"i: "<<i<<"   ho: "<<GW_heads[i]<<"   hc: "<<((head_change[i] / poro) / MM_PER_METER)<<"   hn: "<<GW_heads_new[i]<<endl;
      }
    }

    //Solve for lateral flow
    SolveGWSystem(pModel,Options, tt, GW_heads, GW_heads_new);

    //convert heads back to storage & update statevars with storag and head values
    for (k=0;k<nHRUs;k++)
    {
      param = "GW_HEAD";
      pHRU = pModel->GetHydroUnit(k);                          
      poro = pHRU->GetAquiferProps(0)->porosity;
      bot  = pGWmod->GetGWGeom()->GetGWGeoProperty("BOT_ELEV",k,0);

      GW_stor_new[k]     = (((GW_heads_new[k] - bot) * MM_PER_METER) * poro);      //convert solved heads into change in storage [mm]
      //cout<<"stor: "<<GW_stor_new[k]<<"   head: "<<GW_heads_new[k]<<endl;                      
      iGW  = pModel->GetStateVarIndex(GROUNDWATER,0);                 //*****FIX - needs to work with multiple layers
      //pHRU->SetStateVarValue(iGW, GW_stor_new[k]);
      aPhinew[k][iGW] = GW_stor_new[k];                               //sets aPhinew to redistributed storage value for later updating after routing [mm]
      pGWmod->GetGWGeom()->SetGWGeoProperty(param,GW_heads_new[k],k,0);         //*****FIX - needs to work with multiple layers

      stor_sum2 += GW_stor_new[k];
      /*for(j = 0; j < inc_bal_cnt[k]; j++)
      {
        pModel->IncrementBalance(inc_bal[k][j],k,(GW_stor_new[k] - pHRU->GetStateVarValue(iGW))*tstep);
      }*/
    }
    //for (k=0;k<nHRUs;k++)
    {
      //cout<<setprecision(12)<<"s_preGW: "<<GW_stor[k]<<"   s_postGW: "<<GW_stor_new[k]<<endl;
    }
    //cout<<setprecision(12)<<"s_preGW: "<<stor_sum<<"   s_postGW: "<<stor_sum2<<"  diff: "<<stor_sum2-stor_sum<<endl;
    //cout<<endl;
  }


  //-----------------------------------------------------------------
  //      ROUTING
  //-----------------------------------------------------------------
  double res_ht,res_outflow;
  double down_Q,irr_Q,div_Q, Qwithdrawn;
  int    pDivert;
  res_constraint res_const;

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

      //surface water moved instantaneously from HRU to basin reach/channel storage
      aRouted[p]+=(aPhinew[k][iSW]/MM_PER_METER)*(pHRU->GetArea()*M2_PER_KM2);//aRouted=[m3]

      aPhinew[k][iSW]=0.0;//zero out surface water storage
    }
  }
  // Identify magnitude of flow diversions, calculate inflows  
  for(p=0;p<NB;p++)
  {
    pBasin=pModel->GetSubBasin(p);
    for(int i=0; i<pBasin->GetNumDiversions();i++) {
      div_Q=pBasin->GetDiversionFlow(i,pBasin->GetOutflowRate(),Options,tt,pDivert); //diversions based upon flows at start of timestep
      if(pDivert!=DOESNT_EXIST)
      {
        aQinnew[pDivert]+=div_Q;
      }
      else {
        // \todo[funct] - must handle this in mass balance - should go to cumulative output
		    //Qdiverted+=div_Q;
      }
    }
    aQinnew[p]+=pBasin->GetSpecifiedInflow(t+tstep);
  }
  // Route water over timestep
  // calculations performed in order from upstream (pp=0) to downstream (pp=nSubBasins-1)
  for (pp=0;pp<NB;pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp); //p refers to actual index of basin, pp is ordered list index upstream to down

    pBasin=pModel->GetSubBasin(p);
    if(pBasin->IsEnabled())
    {
      pBasin->UpdateFlowRules(tt,Options);
      
      pBasin->SetInflow(aQinnew[p]);                 // from upstream, diversions, and specified flows

      down_Q=pBasin->GetDownstreamInflow(t);         // treated as additional runoff (period starting)

      pBasin->SetLateralInflow(aRouted[p]/(tstep*SEC_PER_DAY)+down_Q);//[m3/d]->[m3/s]

      pBasin->RouteWater    (aQoutnew,res_ht,res_outflow,res_const,Options,tt);      //Where everything happens!

      Qwithdrawn=0;
      irr_Q=pBasin->ApplyIrrigationDemand(t+tstep,aQoutnew[pBasin->GetNumSegments()-1]);
      Qwithdrawn+=irr_Q;

      for(int i=0; i<pBasin->GetNumDiversions();i++) {
        div_Q=pBasin->GetDiversionFlow(i,pBasin->GetOutflowRate(),Options,tt,pDivert); //diversions based upon flows at start of timestep
        Qwithdrawn+=div_Q;
      }

      pBasin->UpdateOutflows(aQoutnew,irr_Q,res_ht,res_outflow,res_const,Options,tt,false);//actually updates flow values here

      pTo   =pModel->GetDownstreamBasin(p);
      if(pTo!=DOESNT_EXIST)//update downstream inflows
      {
        aQinnew[pTo]+=pBasin->GetOutflowRate()-Qwithdrawn;
      }

    }
  }//end for pp...

  //-----------------------------------------------------------------
  //      MASS CONSTITUENT ROUTING
  //-----------------------------------------------------------------

  //determine total mass loading from HRUs into respective basins (aRoutedMass[p])
  if (nConstituents>0)
  {
    for (p=0;p<NB;p++){
      for (c=0;c<nConstituents;c++){
        aRoutedMass[p][c]=0.0;
        aMinnew    [p][c]=0.0;
      }
    }
    for (k=0;k<nHRUs;k++)
    {
      pHRU=pModel->GetHydroUnit(k);
      if(pHRU->IsEnabled())
      {
        p   =pHRU->GetSubBasinIndex();

        for(c=0;c<nConstituents;c++)
        {
          int iSWmass =pModel->GetStateVarIndex(CONSTITUENT,pModel->GetTransportModel()->GetLayerIndex(c,iSW));
          aRoutedMass[p][c]+=(aPhinew[k][iSWmass])*(pHRU->GetArea()*M2_PER_KM2)/tstep;//aRoutedMass=[mg/d]
          aPhinew[k][iSWmass]=0.0; //empty out
        }
      }
    }

    //Route mass over timestep
    //calculations performed in order from upstream (pp=0) to downstream (pp=nSubBasins-1)
    for (pp=0;pp<NB;pp++)
    {
      p=pModel->GetOrderedSubBasinIndex(pp); //p refers to actual index of basin, pp is ordered list

      pBasin=pModel->GetSubBasin(p);

      if(pBasin->IsEnabled())
      {
        pModel->GetTransportModel()->SetMassInflows    (p,aMinnew[p]);
        pModel->GetTransportModel()->SetLateralInfluxes(p,aRoutedMass[p]); 
        pModel->GetTransportModel()->RouteMass         (p,aMoutnew,aResMass,Options,tt);  //Where everything happens!
        pModel->GetTransportModel()->UpdateMassOutflows(p,aMoutnew,aResMass,aMassOutflow, Options,tt,false); //actually updates mass flow values here

        pTo   =pModel->GetDownstreamBasin(p);
        if(pTo!=DOESNT_EXIST)
        {
          for(int c=0;c<nConstituents;c++){
            aMinnew[pTo][c]+=aMassOutflow[c];
          }
        }
      }
    }//end for pp...
  }//end if (nConstituents>0)

  //-----------------------------------------------------------------
  //          AGGREGATION ACROSS HRU GROUPS
  // \todo[funct]  Should replace with special lateral exchange process
  //-----------------------------------------------------------------
  double agg_phi;
  double agg_area,area;
  double avg_phi;
  int    nHRUgroups=pModel->GetNumHRUGroups();
  for (i=0;i<NS;i++)
  {
    for (int kk=0;kk<nHRUgroups;kk++)
    {
      if (pModel->GetHRUGroup(kk)->IsAggregatorGroup(i))
      {
        CHRUGroup *pAggGroup=pModel->GetHRUGroup(kk);
        //calculate average state variable value
        agg_phi=0.0;
        agg_area=0.0;
        for (int k2=0;k2<pAggGroup->GetNumHRUs();k2++)
        {
          pHRU=pAggGroup->GetHRU(k2);
          k   =pHRU->GetGlobalIndex();
          area=pHRU->GetArea();
          agg_phi +=area*aPhinew[k][i];
          agg_area+=area;
        }
        avg_phi=agg_phi/agg_area;

        //assign average state variable value to all HRUs in aggregation group
        for (int k2=0;k2<pAggGroup->GetNumHRUs();k2++)
        {
          k=pAggGroup->GetHRU(k2)->GetGlobalIndex();
          aPhinew[k][i]=avg_phi;
          //No need to increment cumulative mass balance?
        }
      }
    }
  }

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
  if (t>=Options.duration-Options.timestep)
  {
    if (DESTRUCTOR_DEBUG){cout<<"DELETING STATIC ARRAYS IN MASSENERGYBALANCE"<<endl;}
    for (k=0;k<nHRUs;k++)     {delete [] aPhi   [k];}     delete [] aPhi;       aPhi=NULL;
    for (k=0;k<nHRUs;k++)     {delete [] aPhinew[k];}     delete [] aPhinew;    aPhinew=NULL;
    for (k=0;k<nHRUs;k++)     {delete [] aPhiPrevIter[k];}delete [] aPhiPrevIter; aPhiPrevIter=NULL;
    if(Options.sol_method == ITERATED_HEUN)
    {
      for (j=0;j<nProcesses;j++){delete [] rate_guess[j];}  delete [] rate_guess; rate_guess=NULL;
    }
    delete [] aQinnew;      aQinnew     =NULL;
    delete [] aQoutnew;     aQoutnew    =NULL;
    delete [] aRouted;      aRouted     =NULL;
	delete [] GW_stor;      GW_stor     = NULL; //GWMIGRATE - should name aGW_stor,...
    delete [] GW_stor_new;  GW_stor_new = NULL;
    delete [] GW_heads;     GW_heads    = NULL;
    delete [] GW_heads_new; GW_heads_new= NULL;
    //delete transport static arrays.
    if (nConstituents>0){
      for (p=0;p<NB;p++){
        delete [] aMinnew[p];
        delete [] aRoutedMass[p];
      }

      for (i=0;i<MAX_RIVER_SEGS;i++){
        delete [] aMoutnew[i];
      }
    }
    delete [] aMinnew;      aMinnew     =NULL;
    delete [] aRoutedMass;  aRoutedMass =NULL;
    delete [] aMoutnew;     aMoutnew    =NULL;
    delete [] aResMass;
    delete [] aMassOutflow;
  }
}

