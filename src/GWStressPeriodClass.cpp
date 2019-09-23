/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
----------------------------------------------------------------*/
#include "Properties.h"
#include "GroundwaterClass.h"

/*****************************************************************
   Constructor / Destructor
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the groundwater stress period class constructor
/// \param name [in] String nickname for groundwater stress period profile 
//
CGWStressPeriodClass::CGWStressPeriodClass(const string name)
{
  tag=name;
   if (!DynArrayAppend((void**&)(pAllGWClasses),(void*)(this),NumGWClasses)){
     ExitGracefully("CGroundwater::Constructor: creating NULL groundwater stress period class",BAD_DATA);};
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CGWStressPeriodClass::~CGWStressPeriodClass(){
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING GROUNDWATER STRESS PERIOD CLASS "<<endl;}
}

///////////////////////////////////////////////////////////////////
/// \brief Returns reference to groundwater stress period properties
/// \return Groundwater stress period properties associated with groundwater class
//
const gw_struct *CGWStressPeriodClass::GetGWStruct() const{return &W;}

///////////////////////////////////////////////////////////////////
/// \brief Return nick name identifier of groundwater stress period class
/// \return nick name identifier of groundwater stress period class
//
string                CGWStressPeriodClass::GetTag					 () const{return tag;}

/*****************************************************************
   Static Initialization, Accessors, Destructors
*****************************************************************/
CGWStressPeriodClass **CGWStressPeriodClass::pAllGWClasses=NULL; 
int                    CGWStressPeriodClass::NumGWClasses=0;     


//////////////////////////////////////////////////////////////////
/// \brief Return number of groundwater stress period classes
/// \return Number of groundwater classes
//
int CGWStressPeriodClass::GetNumGWSPClasses(){
  return NumGWClasses;
}


//////////////////////////////////////////////////////////////////
/// \brief Summarize groundwater stress period class information to screen
//
void CGWStressPeriodClass::SummarizeToScreen()
{
  cout<<"Groundwater Stress Period Summary: "<<endl;
  cout<<"     " <<NumGWClasses<<" groundwater properties class in database"<<endl;
  for (int c=0; c<NumGWClasses;c++){
    cout<<"    Groundwater Properties class \""<<pAllGWClasses[c]->GetTag()<<"\": "<<endl;
    cout<<"    Stress period \"" << pAllGWClasses[c]->GetGWStruct()->stressperiod << "\" properties:" << endl;
    const sp SPout = pAllGWClasses[c]->GetGWStruct()->sp_props;
    cout<<"                   Number of ET Records: "<<SPout.numrec_ET<<endl;
    cout<<"                       Number of Drains: "<<SPout.numrec_D<<endl;
    cout<<"              Number of Recharge Values: "<<SPout.numrec_R<<endl;
    cout<<"Number of Hydraulic Conductivity Values: "<<SPout.numrec_K<<endl;
    cout<<"               Number of Storage Values: "<<SPout.numrec_S<<endl;
    cout<<endl;
  }
}


//////////////////////////////////////////////////////////////////
/// \brief Write groundwater stress period class properties to file
/// \param &OUT [out] Output stream to which information is written
/// \todo [add funct] repair such that only parameters that are used by model are written to parameters.csv.
/// one way to do this would be to store the template groundwater structure with class, and use NOT_NEEDED status of variables as indicator.
//
void CGWStressPeriodClass::WriteParamsToFile(ofstream &OUT)
{
  CGWStressPeriodClass *s;
  const gw_struct *w;
  OUT<<endl<<"---Groundwater Stress Period Class Parameters---------------------"<<endl;
  OUT<<"CLASS,";
  OUT<<"STRESS_PERIOD,";
  OUT<<"NUM_REC_ET,NUM_REC_D,NUM_REC_R,NUM_REC_K,NUM_REC_S,CELL_ET,CELL_DRAIN,CELL_RECHARGE,PET_MAX,EXT_DEPTH,DRAIN_ELEV,CONDUCT,RECHARGE,HYD_COND,SS";
  OUT<<endl;
  for (int c=0; c<NumGWClasses;c++)
  {
    s=pAllGWClasses[c];
    w=s->GetGWStruct();
    OUT<<s->GetTag()<<",";
    OUT<<w->stressperiod <<",";
    OUT<<w->sp_props.numrec_ET <<","<<w->sp_props.numrec_D  <<","<<w->sp_props.numrec_R  <<","<<w->sp_props.numrec_K <<","<<w->sp_props.numrec_S <<","<<w->sp_props.numcell_ET  <<","<<w->sp_props.numcell_Drains   <<","<<w->sp_props.numcell_Rech  <<","<<w->sp_props.pet_max  <<","<<w->sp_props.ext_depth   <<","<<w->sp_props.drain_elev   <<","<<w->sp_props.conductance   <<","<<w->sp_props.recharge   <<","<<w->sp_props.hyd_cond   <<","<<w->sp_props.ss   <<",";
    OUT<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all groundwater stress period classes
//
void CGWStressPeriodClass::DestroyAllGWClasses()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL GROUNDWATER CLASSES"<<endl;}
  for (int c=0; c<NumGWClasses;c++){
    delete pAllGWClasses[c];
  }
  delete [] pAllGWClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the groundwater class corresponding to passed string
/// \details Converts string to gwclass
///  if string is invalid, returns NULL
/// \param w [in] GW class identifier (tag or index)
/// \return Reference to groundwater class corresponding to identifier w
//
CGWStressPeriodClass *CGWStressPeriodClass::StringToGWClass(const string s)
{
  string sup=StringToUppercase(s);
  for (int c=0;c<NumGWClasses;c++)
  {
    if (!sup.compare(StringToUppercase(pAllGWClasses[c]->GetTag()))){return pAllGWClasses[c];}
    else if (s_to_i(s.c_str())==(c+1))           {return pAllGWClasses[c];}
  }
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns the groundwater stress period class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] GW class index 
/// \return Reference to groundwater stress period class corresponding to index c
//
const CGWStressPeriodClass *CGWStressPeriodClass::GetGWSPClass(int c)
{
  if ((c<0) || (c>=CGWStressPeriodClass::NumGWClasses)){return NULL;}
  return pAllGWClasses[c];
}

//////////////////////////////////////////////////////////////////
/// \brief Allocate struct array sizes
//
void CGWStressPeriodClass::AllocateArrays(string &param_name, int nodes)
{
  AllocateArrays(W, param_name, nodes);
}

//////////////////////////////////////////////////////////////////
/// \brief Allocate struct array sizes
//
void CGWStressPeriodClass::AllocateArrays(gw_struct &W, string param_name, int nodes)
{
/*  W.sp_props.numcell_ET     = new double[nodes];
  W.sp_props.numcell_Drains = new double[nodes];
  W.sp_props.numcell_Rech   = new double[nodes];
  W.sp_props.numcell_K      = new double[nodes];
  W.sp_props.numcell_S      = new double[nodes];

  W.sp_props.pet_max        = new double[nodes];
  W.sp_props.ext_depth      = new double[nodes];
  W.sp_props.drain_elev     = new double[nodes];
  W.sp_props.conductance    = new double[nodes];            //****FIX - ALL OF THIS NEEDS TO BE DELETED!
  W.sp_props.recharge       = new double[nodes];
  W.sp_props.hyd_cond       = new double[nodes];
  W.sp_props.ss             = new double[nodes];*/

  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("NUMREC_ET"         ))
  {
    W.sp_props.numcell_ET     = new double[nodes];
    W.sp_props.pet_max        = new double[nodes];           //****FIX - NEED TO RELEASE MEMORY
    W.sp_props.ext_depth      = new double[nodes];
  }
  else if (!name.compare("NUMREC_D"     ))
  {
    W.sp_props.numcell_Drains = new double[nodes];
    W.sp_props.drain_elev     = new double[nodes];
    W.sp_props.conductance    = new double[nodes];
  }
  else if (!name.compare("NUMREC_R"       ))
  {
    W.sp_props.numcell_Rech   = new double[nodes];
    W.sp_props.recharge       = new double[nodes];
  }
  else if (!name.compare("NUMREC_K"       ))
  {
    W.sp_props.numcell_K      = new double[nodes];
    W.sp_props.hyd_cond       = new double[nodes];
  }
  else if (!name.compare("NUMREC_S"       ))
  {
    W.sp_props.numcell_S      = new double[nodes];
    W.sp_props.ss             = new double[nodes];
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets default groundwater stress period properties
/// \param &W [out] Groundwater stress period properties class
/// \param is_template [in] True if the default value being set is for the template class
//
void CGWStressPeriodClass::InitializeGWProperties(gw_struct &W, bool is_template)//static
{
  W.stressperiod = 1;      

  W.sp_props.numcell_ET      = NULL;
  W.sp_props.numcell_Drains  = NULL;
  W.sp_props.numcell_Rech    = NULL;
  W.sp_props.numcell_K       = NULL;
  W.sp_props.numcell_S       = NULL;
  W.sp_props.numrec_ET       = 10;
  W.sp_props.numrec_D        = 10;
  W.sp_props.numrec_R        = 10;
  W.sp_props.numrec_K        = 10;
  W.sp_props.numrec_S        = 10;
  W.sp_props.pet_max         = NULL;
  W.sp_props.ext_depth       = NULL;
  W.sp_props.drain_elev      = NULL;
  W.sp_props.conductance     = NULL;
  W.sp_props.recharge        = NULL;
  W.sp_props.hyd_cond        = NULL;
  W.sp_props.ss              = NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the groundwater stress period property corresponding to param_name
/// \param param_name [in] Parameter identifier
//
void  CGWStressPeriodClass::SetGWProperty(string &param_name, 
                                          const double value, int i)
{
  SetGWProperty(W,param_name,value,i);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the gw stress period property corresponding to param_name
/// \note This is declared as a static member because groundwater class 
/// is not instantiated prior to read of .rvg file
/// \param &W [out] Groundwater properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  CGWStressPeriodClass::SetGWProperty(gw_struct &W, 
                                          string     param_name, 
                                          const double value,
                                          int i)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("NUMCELL_ET"         )){W.sp_props.numcell_ET[i]=value;}
  else if (!name.compare("NUMCELL_DRAINS"     )){W.sp_props.numcell_Drains[i]=value;}
  else if (!name.compare("NUMCELL_RECH"       )){W.sp_props.numcell_Rech[i]=value;}
  else if (!name.compare("NUMCELL_K"          )){W.sp_props.numcell_K[i]=value;}
  else if (!name.compare("NUMCELL_S"          )){W.sp_props.numcell_S[i]=value;}
  else if (!name.compare("NUMREC_ET"          )){W.sp_props.numrec_ET=value;}
  else if (!name.compare("NUMREC_D"           )){W.sp_props.numrec_D=value;}
  else if (!name.compare("NUMREC_R"           )){W.sp_props.numrec_R=value;}
  else if (!name.compare("NUMREC_K"           )){W.sp_props.numrec_K=value;}
  else if (!name.compare("NUMREC_S"           )){W.sp_props.numrec_S=value;}
  else if (!name.compare("STRESSPERIOD"       )){W.stressperiod=value;}
  else if (!name.compare("PET_MAX"            )){W.sp_props.pet_max[i]=value;}
  else if (!name.compare("EXT_DEPTH"          )){W.sp_props.ext_depth[i]=value;}
  else if (!name.compare("DRAIN_ELEV"         )){W.sp_props.drain_elev[i]=value;}
  else if (!name.compare("CONDUCTANCE"        )){W.sp_props.conductance[i]=value;}
  else if (!name.compare("RECHARGE"           )){W.sp_props.recharge[i]=value;}
  else if (!name.compare("HYD_COND"           )){W.sp_props.hyd_cond[i]=value;}
  else if (!name.compare("SS"                 )){W.sp_props.ss[i]=value;}
  else{
    WriteWarning("CGWStressPeriodClass::SetGWProperty: Unrecognized/invalid groundwater stress period parameter name ("+name+") in .rvg file",false);
  }
}
///////////////////////////////////////////////////////////////////////////
/// \brief Returns gw stress period property value corresponding to param_name
/// \param param_name [in] Parameter name
/// \return GW properties corresponding to parameter name
//
double CGWStressPeriodClass::GetGWProperty(string param_name, int i) const
{
  return GetGWProperty(W,param_name, i);
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns groundwater stress period property value corresponding to param_name from structure provided
/// \param gw_struct [in] Groundwater structure
/// \param param_name [in] Parameter name
/// \return Groundwater stress period properties corresponding to parameter name
//
double CGWStressPeriodClass::GetGWProperty(const gw_struct &W, string param_name, int i)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("NUMCELL_ET"           )){return W.sp_props.numcell_ET[i];}
  else if (!name.compare("NUMCELL_D"            )){return W.sp_props.numcell_Drains[i];}
  else if (!name.compare("NUMCELL_RECH"         )){return W.sp_props.numcell_Rech[i];}
  else if (!name.compare("NUMCELL_K"            )){return W.sp_props.numcell_K[i];}
  else if (!name.compare("NUMCELL_S"            )){return W.sp_props.numcell_S[i];}
  else if (!name.compare("NUMREC_ET"            )){return W.sp_props.numrec_ET;}
  else if (!name.compare("NUMREC_D"             )){return W.sp_props.numrec_D;}
  else if (!name.compare("NUMREC_R"             )){return W.sp_props.numrec_R;}
  else if (!name.compare("NUMREC_K"             )){return W.sp_props.numrec_K;}
  else if (!name.compare("NUMREC_S"             )){return W.sp_props.numrec_S;}
  else if (!name.compare("STRESSPERIOD"         )){return W.stressperiod;}
  else if (!name.compare("PET_MAX"              )){return W.sp_props.pet_max[i];}
  else if (!name.compare("EXT_DEPTH"            )){return W.sp_props.ext_depth[i];}
  else if (!name.compare("DRAIN_ELEV"           )){return W.sp_props.drain_elev[i];}
  else if (!name.compare("CONDUCTANCE"          )){return W.sp_props.conductance[i];}
  else if (!name.compare("RECHARGE"             )){return W.sp_props.recharge[i];}
  else if (!name.compare("HYD_COND"             )){return W.sp_props.hyd_cond[i];}
  else if (!name.compare("SS"                   )){return W.sp_props.ss[i];}
  else{
    string msg="CGWStressPeriodClass::GetGWProperty: Unrecognized/invalid groundwater stress period parameter name in .rvg file: "+name;
    ExitGracefully(msg.c_str(),BAD_DATA);
    return 0.0;
  }
}