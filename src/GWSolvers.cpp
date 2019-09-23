/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2015 the Raven Development Team
----------------------------------------------------------------*/

#include "RavenInclude.h"
#include "Model.h"
#include "GroundwaterClass.h"
#include "SparseMatrix.h"

#ifdef _ARMADILLO_
#include<armadillo>
using namespace arma;
#endif


double CalcFhup(double h_int,double top,double bot);
void CreateJacobian(const optStruct   &Options, CGWGeometryClass  *pGWGeomet, CGWStressPeriodClass *pGWSP, double *Heads_new,double *Heads_prev, double *HeadsPrevT, double timestep,int *&ija,CSparseMatrix *&sNR,double *&matB);

///////////////////////////////////////////////////////////////////
/// \brief Solves system of mass balance ODEs/PDEs for one timestep for GW System
///
/// \param *pModel [in & out] Model
/// \param &Options [in] Global model options information
/// \param &tt [in] Time step over which energy and mass balance equations are to be solved
//
void SolveGWSystem     ( CModel            *pModel, 
                          const optStruct  &Options,
                        const time_struct  &tt,
                                   double  *GW_head_prev,
                                   double  *GW_head_new)
{

#ifdef _ARMADILLO_
  int i;
  int outer_iter_cnt= 0;
  int max_outer_iter= 15;              //****FIX - get from input files
  //int max_inner_iter= 50;            //****FIX - get from input files

  CGWGeometryClass  *pGWGeomet=pModel->GetGroundwaterModel()->GetGWGeos(0);
  CGWStressPeriodClass  *pGWSP=pModel->GetGroundwaterModel()->GetGWSPs(0); //GWMIGRATE - only works with single stress period?

  double outer_tol_check;
  double outer_tol_check_numer;
  double outer_tol_check_denom;
  double relax_fact = 1;//0.85;
  double outer_head_delta, outer_head_crit, outer_head_initial;

  double tstep  = Options.timestep;      //timestep
  double t;
  t             = tt.model_time;

  static double *aHeadPrevT;                //heads at end of previous timestep
  static double *aHeadTemp;                 //temporary head storage post inner solver step
  static double *aHeadPrevIter;             //heads from previous iteration
  static double *aHeadCurrIter;             //heads for current iteration
  
  int    *SparseIndex;
  int     NonZero;
  double *MatB;

  CSparseMatrix *pSpMat = NULL;             //empty sparse matrix - double *spA, int *indices, double *b

  int nGWnodes;

  nGWnodes = int(pGWGeomet->GetGWGeoProperty("num_nodes"));

  aHeadPrevT    = new double [nGWnodes];
  aHeadTemp     = new double [nGWnodes];
  aHeadPrevIter = new double [nGWnodes];
  aHeadCurrIter = new double [nGWnodes];
 
  aHeadPrevT    = GW_head_new;
  aHeadPrevIter = GW_head_new;
  aHeadCurrIter = GW_head_new;

  
  
  //Nonlinear solvers surrounding calls to linear solvers
  if(Options.gw_solver_outer == GWSOL_NEWTONRAPHSON)
  {
    arma::vec b_vec;
    arma::vec x_vec;
    double err;
    double inner_tol = 1.0e-7;
    double outer_tol = 1.0e-10;                      //****FIX - read in from input files
    outer_head_crit   = 0.001;
    
    NonZero     = nGWnodes + (int) (pGWGeomet->GetGWGeoStruct()->num_connections);
    MatB        = new double [nGWnodes];
    SparseIndex = new int [NonZero + 1];
    b_vec       = vec(nGWnodes);
    x_vec       = vec(nGWnodes);
  
    //outer solver (nonlinear) loop until convergence or max iteration count
    do
    {
      double head_sum1 = 0;
      double head_sum2 = 0;

      outer_head_initial = 0;
      outer_head_delta   = 0;
      for(i = 0; i < nGWnodes; i++){outer_head_initial += aHeadCurrIter[i];}
      outer_iter_cnt++;                               //increase outer iteration counter

      for(i = 0; i < 16; i++)
      {
        head_sum1 += aHeadCurrIter[i];
      }
      CreateJacobian(Options, pGWGeomet, pGWSP, aHeadCurrIter,aHeadPrevIter, aHeadPrevT,tstep,*&SparseIndex,*&pSpMat,*&MatB);
      
      for(i = 0; i < nGWnodes; i++)
      {
        b_vec[i] = MatB[i];                           //assign rhs and heads to vectors
        x_vec[i] = aHeadCurrIter[i];
      }
      
      if(Options.gw_solver_inner == GWSOL_BICGSTAB)
      {
        //inner solver (linear) loop until convergence or max iteration count
        pSpMat->BCGstab(b_vec,*&x_vec,nGWnodes,err,inner_tol,pSpMat);

      }
      else if(Options.gw_solver_inner == GWSOL_PCGU)
      {
        //To be added....

        ExitGracefully("GWSolversFile: PCGU Linear Solver not implemented yet",STUB);
      }
      else
      {
        ExitGracefully("GWSolversFile: Invalid Linear Solver Option",STUB);
      }      

      for(i = 0; i < nGWnodes; i++)
      {
        aHeadTemp[i] = (x_vec[i] * relax_fact) + (aHeadPrevIter[i] * (1 - relax_fact));//x_vec[i];                 //update current heads to current iteration values
        head_sum2   += x_vec[i];//(x_vec[i] * relax_fact) + (aHeadPrevIter[i] * (1 - relax_fact)); //aHeadCurrIter[i];
        //cout<<"prevh: "<<aHeadPrevIter[i]<<" currh: "<<aHeadCurrIter[i]<<" temph: "<<aHeadTemp[i]<<endl;
      }

      aHeadPrevIter = aHeadCurrIter;          //update prev heads to current iteration values
      aHeadCurrIter = aHeadTemp;              //update current heads to end of iteration values

      for(i = 0; i < nGWnodes; i++)
      {
        //cout<<"prevh2: "<<aHeadPrevIter[i]<<" currh2: "<<aHeadCurrIter[i]<<endl;
      }
      //cout<<setprecision(12)<<"s_preB: "<<((head_sum1 + (366*16)) * 1000 * 0.363)<<"   s_postB: "<<((head_sum2 + (366*16)) * 1000 * 0.363)<<"   diff: "<<((head_sum2 + (366*16)) * 1000 * 0.363)-((head_sum1 + (366*16)) * 1000 * 0.363)<<endl;
      //cout<<setprecision(12)<<"h_preB: "<<head_sum1<<"   h_postB: "<<head_sum2<<"   diff: "<<head_sum2-head_sum1<<endl;
      //cout<<endl;

      //check for outer iteration convergence
      outer_tol_check = 0;
      outer_tol_check_numer = 0;
      outer_tol_check_denom = 0;

      for(i = 0; i < nGWnodes; i++)
      {
        //outer_tol_check += aHeadCurrIter[i]-aHeadPrevIter[i];           
        outer_tol_check_numer += pow((aHeadCurrIter[i]-aHeadPrevIter[i]),2);                     //norm(delta h)/norm(current h)
        outer_tol_check_denom += pow((aHeadCurrIter[i]),2);
        outer_head_delta += aHeadCurrIter[i];
        //cout<<"   prev: "<<aHeadPrevIter[i]<<"   curr: "<<aHeadCurrIter[i]<<"   diff: "<<(aHeadCurrIter[i]-aHeadPrevIter[i])<<endl;
      }
      outer_head_delta = outer_head_delta - outer_head_initial;       //check for head change tolerance

      outer_tol_check_numer = sqrt(outer_tol_check_numer);
      outer_tol_check_denom = sqrt(outer_tol_check_denom);
      outer_tol_check = outer_tol_check_numer / outer_tol_check_denom;    //check for convergence
      //cout<<setprecision(12)<<"outer_iter: "<<outer_iter_cnt<<"   err: "<<outer_tol_check<<"   head del: "<<outer_head_delta<<endl;
      //cout<<endl;
    }while(outer_tol_check > outer_tol && outer_iter_cnt < max_outer_iter && abs(outer_head_delta) > outer_head_crit);   //end of outer loop
  }
  else if(Options.gw_solver_outer == GWSOL_PICARD)
  {
    //To be added...

    ExitGracefully("GWSolversFile: Picard Nonlinear Solver not implemented yet",STUB);
  }
  else
  {
    ExitGracefully("GWSolversFile: Invalid Nonlinear Solver Option",STUB);
  }

  if(outer_iter_cnt == max_outer_iter)
  {
    cout<<"Outer Solver failed to converge!"<<endl;
  }

  for(i = 0; i < nGWnodes; i++)
  {
    GW_head_new[i] = aHeadCurrIter[i];
  }
  /*for(i =0; i<16; i++)
  {
    cout<<"head["<<i<<"]: "<<GW_head_new[i]<<endl;
  }*/
  //delete static arrays (only called once)=========================
  /*if (t>=Options.duration-Options.timestep)
  {
    if (DESTRUCTOR_DEBUG){cout<<"DELETING STATIC ARRAYS IN GWSOLVER"<<endl;}
    //delete and release memory
    delete [] aHeadPrevT;     aHeadPrevT    = NULL;
    delete [] aHeadPrevIter;  aHeadPrevIter = NULL;
    delete [] aHeadCurrIter;  aHeadCurrIter = NULL;
    //delete [] MatB;           MatB          = NULL;
    //delete [] SparseIndex;    SparseIndex   = NULL;
  }*/
#else
  ExitGracefully("Groundwater models cannot be run without Armadillo library",RUNTIME_ERR);
#endif
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Create Jacobian Matrix for Newton-Raphson Method
/// 			
///
/// <remarks>	asnowdon, 3/18/2015. </remarks>
///
/// <param name="pGWGeomet">		GW Geometry. </param>
/// <param name="Heads">		[in,out] Perturbed heads</param>
/// <param name="Heads_prev">		[in,out] Previous interation heads. </param>
/// <param name="timestep">		timestep size. </param>
/// <param name="ija">	[in,out] sparse matrix indices. </param>
/// <param name="sNR">	[in,out] sparse matrix. </param>
/// <param name="matB">	[in,out] sparse matrix vector b. </param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateJacobian(const optStruct   &Options, CGWGeometryClass  *pGWGeomet, CGWStressPeriodClass *pGWSP, double *Heads_new,double *Heads_prev, double *HeadsPrevT, double timestep,int *&ija,CSparseMatrix *&sNR,double *&matB)
{
#ifdef _ARMADILLO_
  int i, j, k, m;
  int GWFSize = 0;
  int nNonZero = 0;
  int n_conn;
  
  double *cell_volume;
  double *spec_stor;
  double *hydCond;
  double *Hcof;            //sum of all terms that are coefficients of hn for cell n
  double top_n, top_m, bot_n, bot_m;
  double *RHS_prev;             //the right hand side of the balance  equation
  double *RHS_new;
  double len_nm;          //length between nodes n and m of adjoining cells [m]
  double C_nm;              //conductance coefficient 
  double width;
  double sat_thick_avg;
  double *R;
  double *R_diag_new;
  double *R_diag_prev;
  double *sA;
  double R_offdiag_new;
  double fh_diag;
  double fh_offdiag;
  double fhup;
  double fh_diag_prev;
  double fh_diag_new;
  double fh_offdiag_prev;
  double fh_offdiag_new;
  double K_avg;
  double R_offdiag_temp;
  //double R_offdiag_R;
  
  bool n_ups;               //switch to determine if n or m is upstream

  //Get number of nodes and horizontal connections to determine number of nonzero elements for sparse matrix
  GWFSize = int (pGWGeomet->GetGWGeoStruct()->num_nodes);
  nNonZero = GWFSize + (int) (pGWGeomet->GetGWGeoStruct()->num_connections);

  sA              = new double [nNonZero + 1];
  ija             = new int    [nNonZero + 1];
  matB            = new double [GWFSize];
  cell_volume     = new double [GWFSize];
  RHS_prev        = new double [GWFSize];
  RHS_new         = new double [GWFSize];
  Hcof            = new double [GWFSize];
  R               = new double [GWFSize];
  R_diag_prev     = new double [GWFSize];
  R_diag_new      = new double [GWFSize];
  spec_stor       = new double [GWFSize];
  hydCond         = new double [GWFSize];
  

  for(i = 0; i < GWFSize; i++)
  {
    //Calculate cell volumes
    bot_n          = pGWGeomet->GetGWGeoStruct()->bot_elev[i];
    top_n          = pGWGeomet->GetGWGeoStruct()->top_elev[i];
    cell_volume[i] = pGWGeomet->GetGWGeoStruct()->cell_area[i] * (top_n - bot_n);
    //get conductivities and storage coefficient
    hydCond[i]     = pGWSP->GetGWProperty("HYD_COND",i);        // K in m/d   ****FIX - take from soil profile or groundwater input????
    spec_stor[i]   = pGWSP->GetGWProperty("SS",i);             // Ss in 1/m
  }

  //define size of sparse index and set index k of sparse array after diagonals
  ija[0] = GWFSize+1;
  k = GWFSize;

  //loop through connections between cells and store coeffecients in sparse format
  for(i = 0; i < GWFSize; i++)
  {
    RHS_prev[i]  = ((-spec_stor[i] * cell_volume[i]  * HeadsPrevT[i])/ timestep);
    RHS_new[i]   = ((-spec_stor[i] * cell_volume[i]  * HeadsPrevT[i])/ timestep);
    Hcof[i]      = ((-spec_stor[i] * cell_volume[i] ) / timestep);

    /*if((Heads_new[i] - Heads_prev[i]) == 0)
    {
      dHead = 0;
    }
    else
    {
      dHead = (RHS_new[i] - RHS_prev[i])/(Heads_new[i] - Heads_prev[i]);
      //cout<<"dh: "<<dHead<<endl;
    }*/

    R[i]         = RHS_new[i];//(-Hcof[i] * Heads_new[i]) + RHS_new[i];// - (dHead * Heads_new[i]);         //NOT SURE ABOUT dRHS/dh value
    R_diag_new[i]= Hcof[i];// - dHead;
    //cout<<" hn: "<<Heads_new[i]<<"   hp: "<<Heads_prev[i]<<"   ht: "<<HeadsPrevT[i]<<endl;
    n_conn = (int) pGWGeomet->GetGWGeoStruct()->num_cell_connections[i] - 1;
    //cout<<"rp: "<<RHS_prev[i]<<"   rn: "<<RHS_new[i]<<"   hcof: "<<Hcof[i]<<"   r: "<<R[i]<<"   rdn: "<<R_diag_new[i]<<endl;

    for(j = 0; j < n_conn; j++)
    {
      bool up_test = false;
      k++;
      m = int (pGWGeomet->GetGWGeoStruct()->cell_connections[i][j]) - 1;      //index of gw cell array (not cell number)

      bot_n = pGWGeomet->GetGWGeoStruct()->bot_elev[i];     //top or current cell
      top_n = pGWGeomet->GetGWGeoStruct()->top_elev[i];     //bottom of current cell
      top_m = pGWGeomet->GetGWGeoStruct()->top_elev[m];    //top of connecting cell
      bot_m = pGWGeomet->GetGWGeoStruct()->bot_elev[m];    //bottom of connecting cell

      len_nm    = pGWGeomet->GetGWGeoStruct()->node_interface_length[i][j];                            //****FIX - need to access second node length to be most robust
      width     = pGWGeomet->GetGWGeoStruct()->face_width[i][j];

      //determine upstream direction
      if(Heads_new[i] > Heads_new[m]){n_ups = true;} else if(Heads_new[i] < Heads_new[m]){n_ups = false;} else{up_test = true;}

      if(n_ups == true)
      {
        fh_offdiag   = 0;
        fh_diag_prev = CalcFhup(Heads_prev[i], top_n, bot_n);
        fh_diag_new  = CalcFhup(Heads_new[i], top_n, bot_n);
        if((Heads_new[i] - Heads_prev[i]) == 0){fh_diag = 0;} else{fh_diag = (fh_diag_new - fh_diag_prev) / (Heads_new[i] - Heads_prev[i]);}
        fhup = CalcFhup(Heads_new[i], top_n, bot_n);

        //sat_thick_avg = Heads_new[i] - bot_n;
      }
      else
      {
        fh_diag         = 0;
        fh_offdiag_prev = CalcFhup(Heads_prev[m], top_m, bot_m);
        fh_offdiag_new  = CalcFhup(Heads_new[m], top_m, bot_m);
        if((Heads_new[m] - Heads_prev[m]) == 0){fh_offdiag = 0;} else{fh_offdiag = (fh_offdiag_new - fh_offdiag_prev) / (Heads_new[m] - Heads_prev[m]);}
        fhup = CalcFhup(Heads_new[m], top_m, bot_m);

        //sat_thick_avg = Heads_new[m] - bot_m;
      }

      if(up_test == true)
      {
        fh_diag = 1;
        fh_offdiag = 1;
        fhup = CalcFhup(Heads_new[i], top_n, bot_n);
        //fhup_m = 0;

        //sat_thick_avg = ((top_n - bot_n) + (top_m - bot_m)) / 2;//((Heads_new[i] - bot_n) + (Heads_new[m] - bot_m)) / 2;
      }

      //Calculate inter-cell conductances
      //sat_thick_avg = min((Heads_new[i] - bot_n),(Heads_new[m] - bot_m));      //nested cells sat
      //sat_thick_avg = ((Heads_new[i] - bot_n) + (Heads_new[m] - bot_m)) / 2;//((top_n - bot_n) + (top_m - bot_m)) / 2;
      sat_thick_avg = ((top_n - bot_n) + (top_m - bot_m)) / 2;
      K_avg         = 2*len_nm * (hydCond[i] * hydCond[m]) / (len_nm * hydCond[i] + len_nm * hydCond[m]);     //harmonic mean of conductivities
      C_nm          = (sat_thick_avg * width * K_avg) / (len_nm + len_nm);

      //cout<<"r: "<<R[i]<<"   c: "<<C_nm<<"   fhd: "<<fh_diag<<"   hprevm: "<<Heads_prev[m]<<"   hprevn: "<<Heads_prev[i]<<"   fhoff: "<<fh_offdiag<<endl;
      //find residuals
      R[i] += ((((C_nm * fh_diag * Heads_new[m]) - (C_nm * fh_diag * Heads_new[i])) * Heads_new[i]) + (((C_nm * fh_offdiag * Heads_new[m]) - (C_nm * fh_offdiag * Heads_new[i])) * Heads_new[m]));
      //R[i] += -C_nm * fhup_n * Heads_new[m] + C_nm * fhup_n * Heads_new[i];
      //R[i] += ((C_nm * fh_diag * Heads_new[m] * Heads_new[i]) - (C_nm * fh_diag * Heads_new[i] * Heads_new[i])) + ((C_nm * fh_offdiag * Heads_new[m] * Heads_new[m]) - (C_nm * fh_offdiag * Heads_new[i] * Heads_new[m]));
      //R[i] += ((C_nm * fh_diag * Heads_new[m] - C_nm * fh_diag * Heads_new[i]) * Heads_new[i]) + ((C_nm * fh_offdiag * Heads_new[m] - C_nm * fh_offdiag * Heads_new[i]) * Heads_new[m]);
      //R[i] += (-C_nm * Heads_new[m] + C_nm * Heads_new[i]) + ((-C_nm * Heads_new[i]) * Heads_new[i]) + ((C_nm * Heads_new[m]) * Heads_new[m]);         //no smoothing function
      //R1 += -C_nm;         //no smoothing function
      //R2 += C_nm * Heads_new[m];

      //cout<<"r: "<<R[i]<<endl;
      //calculate diagonal terms
      R_diag_new[i]  += (-C_nm * fhup) + ((C_nm * fh_diag * Heads_new[m]) - (C_nm * fh_diag * Heads_new[i]));
      //R_diag_new[i]  += -C_nm * Heads_new[i];         //no smoothing function
      //R_diag_new[i]  += -C_nm;         //no smoothing function

      //calculate off-diagonal terms
      R_offdiag_new  = (C_nm * fhup) + ((C_nm * fh_offdiag * Heads_new[m]) - (C_nm * fh_offdiag * Heads_new[i]));
      //R_offdiag_new  = C_nm * Heads_new[m];         //no smoothing function
      //R_offdiag_new  = C_nm;         //no smoothing function
      R_offdiag_temp = R_offdiag_new * Heads_new[m];
      // R_offdiag_R   += R_offdiag_temp;

      //fill off diagonals of matrix in sparse format
      sA[k] = R_offdiag_new;
      ija[k] = int (pGWGeomet->GetGWGeoStruct()->cell_connections[i][j]) - 1;
    } //end of j if loop
    //R[i] = (R1 * Heads_new[i]) + R2;
    //R_diag_new[i] = R_diag_new[i] * Heads_new[i];

    //fill diagonals of matrix in sparse format
    sA[i] = R_diag_new[i];
    ija[i+1] = k + 1;

    //fill initial matrix b
    matB[i] =  R[i];
    //matB[i] =  R[i] + (R_diag_new[i] * Heads_new[i]) + R_offdiag_R;
  }  //end of i if loop

  //make new sparse matrix class and store sparse matrix in it
  sNR = new CSparseMatrix(sA,ija,matB,GWFSize);

  //free memory
  /*sA              = new double [NonZero + 1];
  ija             = new int    [NonZero + 1];
  matB            = new double [GWFSize];
  cell_volume     = new double [GWFSize];
  RHS_prev        = new double [GWFSize];
  RHS_new         = new double [GWFSize];
  Hcof            = new double [GWFSize];
  R               = new double [GWFSize];
  R_diag_prev     = new double [GWFSize];
  R_diag_new      = new double [GWFSize];
  spec_stor       = new double [GWFSize];
  hydCond         = new double [GWFSize];*/
#endif
  return;
}

double CalcFhup(double h_int,double top,double bot)
{
  double discont = 1e-06;
  double A_omega;
  double Sf;
  double fhup;

  A_omega = 1/(1-discont);
  Sf = (h_int - bot) / (top - bot);

  if(Sf < 0)                                {fhup = 0;}
  else if(Sf < discont && Sf >= 0)          {fhup = (A_omega/(2*discont))*pow(Sf,2);}
  else if(Sf < (1-discont) && Sf >= discont){fhup = A_omega*Sf + 0.5*(1-A_omega);}
  else if(Sf < 1 && Sf >= (1-discont))      {fhup = 1 - (A_omega/(2*discont)) * pow((1-Sf),2);}
  else                                      {fhup = 1;}

  return fhup;
}