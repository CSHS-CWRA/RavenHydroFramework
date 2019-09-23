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
COverlapExchangeClass::COverlapExchangeClass()
{
   if (!DynArrayAppend((void**&)(pAllOEClasses),(void*)(this),NumOEClasses)){
     ExitGracefully("COverlapExchange::Constructor: creating NULL Overlap/Exchange class",BAD_DATA);};
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
COverlapExchangeClass::~COverlapExchangeClass(){
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING OVERLAP/EXCHANGE CLASS "<<endl;}
}

///////////////////////////////////////////////////////////////////
/// \brief Returns reference to groundwater stress period properties
/// \return Groundwater stress period properties associated with groundwater class
//
const oe_struct *COverlapExchangeClass::GetOEStruct() const{return &E;}

///////////////////////////////////////////////////////////////////
/// \brief Return nick name identifier of groundwater stress period class
/// \return nick name identifier of groundwater stress period class
//
string                COverlapExchangeClass::GetTag					 () const{return tag;}

/*****************************************************************
   Static Initialization, Accessors, Destructors
*****************************************************************/
COverlapExchangeClass **COverlapExchangeClass::pAllOEClasses=NULL; 
int                     COverlapExchangeClass::NumOEClasses=0;     


//////////////////////////////////////////////////////////////////
/// \brief Return number of groundwater stress period classes
/// \return Number of groundwater classes
//
int COverlapExchangeClass::GetNumOEClasses(){
  return NumOEClasses;
}


//////////////////////////////////////////////////////////////////
/// \brief Summarize groundwater stress period class information to screen
//
void COverlapExchangeClass::SummarizeToScreen()
{
  string nHRU = "NHRU";
  string NGWC = "NGWCELLS";
  double over_per;
  COverlapExchangeClass *e;
  const oe_struct *OEout;

  cout<<"Overlap/Exchange Summary: "<<endl;
  cout<<"     " <<NumOEClasses<<" Overlap/Exchange class in database"<<endl;
  for (int c=0; c<NumOEClasses;c++){
    cout<<"       Overlap Properties" << endl;
    e=pAllOEClasses[c];
    OEout=e->GetOEStruct();
        
    for(int i = 0; i < e->GetOEProperty(nHRU);i++)
    {
      int conns = 0;
      for( int j = 0; j < e->GetOEProperty(NGWC); j++)
      {
        over_per = OEout->overlap[i][j];
        if(over_per > 0){conns++;}
        //cout<<"p: "<<over_per<<"   c:"<<conns<<endl;
      }
        cout<<"                HRU #: "<<i+1<<endl;
        cout<<"   Connected GW Cells: "<<conns<<endl;
      
    }
    cout<<endl;
    cout<<"    Exchange Definition for cell \"";
    if(OEout->cell == -9999){cout<<"ALL";}else{cout<< OEout->cell;}
    cout<< "\" properties:" << endl;
    cout<<"              ET Relationship Type: "<<OEout->oe_props.ET_type<<endl;
    cout<<"          Drains Relationship Type: "<<OEout->oe_props.D_type<<endl;
    cout<<"  Saturated Area Relationship Type: "<<OEout->oe_props.SA_type<<endl;
    cout<<endl;
  }
}


//////////////////////////////////////////////////////////////////
/// \brief Write groundwater stress period class properties to file
/// \param &OUT [out] Output stream to which information is written
/// \todo [add funct] repair such that only parameters that are used by model are written to parameters.csv.
/// one way to do this would be to store the template groundwater structure with class, and use NOT_NEEDED status of variables as indicator.
//
void COverlapExchangeClass::WriteParamsToFile(ofstream &OUT)
{
  COverlapExchangeClass *s;
  const oe_struct *e;

  OUT<<endl<<"---Overlap/ Exchange Class Parameters---------------------"<<endl;
  OUT<<"NHRU,NGWCELLS,OVERLAP,CELL,ET RELATIONSHIP TYPE,DRAIN RELATIONSHIP TYPE, SAT. AREA RELATIONSHIP TYPE";
  OUT<<endl;
  for (int c=0; c<NumOEClasses;c++)
  {
    s=pAllOEClasses[c];
    e=s->GetOEStruct();
    OUT<<e->nHRU <<","<<e->nGWCells  <<","<<"ADD OVERLAP OUTPUT"<<",";
    if(e->cell == -9999){OUT<<"ALL";}else{OUT<<e->cell;}
    OUT<<","<<e->oe_props.ET_type<<","<<e->oe_props.D_type<<","<<e->oe_props.SA_type<<",";
    OUT<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all overlap/exchange classes
//
void COverlapExchangeClass::DestroyAllOEClasses()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL OVERLAP/EXCHANGE CLASSES"<<endl;}
  for (int c=0; c<NumOEClasses;c++){
    delete pAllOEClasses[c];
  }
  delete [] pAllOEClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the overlap/exchange class corresponding to passed string
/// \details Converts string to oeclass
///  if string is invalid, returns NULL
/// \param e [in] OE class identifier (tag or index)
/// \return Reference to overlap/exchange class corresponding to identifier e
//
COverlapExchangeClass *COverlapExchangeClass::StringToOEClass(const string s)
{
  string sup=StringToUppercase(s);
  for (int c=0;c<NumOEClasses;c++)
  {
    if (!sup.compare(StringToUppercase(pAllOEClasses[c]->GetTag()))){return pAllOEClasses[c];}
    else if (s_to_i(s.c_str())==(c+1))           {return pAllOEClasses[c];}
  }
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns the overlap/exchange class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] OE class index 
/// \return Reference to overlap/exchange class corresponding to index c
//
const COverlapExchangeClass *COverlapExchangeClass::GetOEClass(int c)
{
  if ((c<0) || (c>=COverlapExchangeClass::NumOEClasses)){return NULL;}
  return pAllOEClasses[c];
}

//////////////////////////////////////////////////////////////////
/// \brief Allocate struct array sizes
//
void COverlapExchangeClass::AllocateOEArrays(int nodes)
{
  AllocateOEArrays(E,nodes);
}

//////////////////////////////////////////////////////////////////
/// \brief Allocate struct array sizes
//
void COverlapExchangeClass::AllocateOEArrays(oe_struct &E, int nodes)
{
  E.oe_props.ET_a       = new double[nodes];
  E.oe_props.ET_b       = new double[nodes];
  E.oe_props.D_a        = new double[nodes];
  E.oe_props.D_b        = new double[nodes];
  E.oe_props.SA_a       = new double[nodes];
  E.oe_props.SA_b       = new double[nodes];

  E.overlap             = new double*[nodes];
  
  for(int i = 0; i < nodes; i++)
  {
    E.overlap[i]        = new double[nodes];
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets default overlap/exchange properties
/// \param &E [out] Overlap/Exchange properties class
/// \param is_template [in] True if the default value being set is for the template class
//
void COverlapExchangeClass::InitializeOEProperties(oe_struct &E, bool is_template)//static
{
  E.nHRU                = 1;
  E.nGWCells            = 1;
  E.overlap             = NULL;

  E.cell                = -9999;
  E.oe_props.ET_type    = "POWER_LAW";
  E.oe_props.D_type     = "POWER_LAW";
  E.oe_props.SA_type    = "POWER_LAW";

  E.oe_props.ET_a       = NULL;
  E.oe_props.ET_b       = NULL;
  E.oe_props.D_a        = NULL;
  E.oe_props.D_b        = NULL;
  E.oe_props.SA_a       = NULL;
  E.oe_props.SA_b       = NULL;


}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the overlap/exchange property corresponding to param_name
/// \param param_name [in] Parameter identifier
//
void  COverlapExchangeClass::SetOEProperty(string &param_name, 
                                          const double value, int i, int j)
{
  SetOEProperty(E,param_name,value,i,j);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the overlap/exchange property corresponding to param_name
/// \param param_name [in] Parameter identifier
//
void  COverlapExchangeClass::SetOEProperty(string &param_name, 
                                          const double value)
{
  SetOEProperty(E,param_name,value);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the overlap/exchange property corresponding to param_name
/// \param param_name [in] Parameter identifier
//
void  COverlapExchangeClass::SetOEProperty(string &param_name, 
                                          const double value,
                                          int i)
{
  SetOEProperty(E,param_name,value, i);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the overlap/exchange property corresponding to param_name
/// \param param_name [in] Parameter identifier
//
void  COverlapExchangeClass::SetOEProperty(string &param_name, 
                                          const int value)
{
  SetOEProperty(E,param_name,value);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the overlap/exchange property corresponding to param_name
/// \param param_name [in] Parameter identifier
//
void  COverlapExchangeClass::SetOEProperty(string &param_name, 
                                          string value)
{
  SetOEProperty(E,param_name,value);
}

//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the overlap/exchange property corresponding to param_name
/// \note This is declared as a static member because overlap/exchange class 
/// is not instantiated prior to read of .rvv/.rve files
/// \param &W [out] Overlap/Exchange properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  COverlapExchangeClass::SetOEProperty(oe_struct &E, 
                                          string     param_name, 
                                          const double value)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("NHRU"               )){E.nHRU=(int)value;}
  else if (!name.compare("CELL"               )){E.cell=(int)value;}
  else if (!name.compare("NGWCELLS"           )){E.nGWCells=(int)value;}

  else{
    WriteWarning("COverlapExchangeClass::SetOEProperty: Unrecognized/invalid overlap/exchange parameter name ("+name+") in .rvv/.rve file",false);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the overlap/exchange property corresponding to param_name
/// \note This is declared as a static member because overlap/exchange class 
/// is not instantiated prior to read of .rvv/.rve files
/// \param &W [out] Overlap/Exchange properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  COverlapExchangeClass::SetOEProperty(oe_struct &E, 
                                          string     param_name, 
                                          const double value, int i)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("ET_A"               )){E.oe_props.ET_a[i]=value;}
  else if (!name.compare("ET_B"               )){E.oe_props.ET_b[i]=value;}
  else if (!name.compare("D_A"                )){E.oe_props.D_a[i]=value;}
  else if (!name.compare("D_B"                )){E.oe_props.D_b[i]=value;}
  else if (!name.compare("SA_A"               )){E.oe_props.SA_a[i]=value;}
  else if (!name.compare("SA_B"               )){E.oe_props.SA_b[i]=value;}


  else{
    WriteWarning("COverlapExchangeClass::SetOEProperty: Unrecognized/invalid overlap/exchange parameter name ("+name+") in .rvv/.rve file",false);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the overlap/exchange property corresponding to param_name
/// \note This is declared as a static member because overlap/exchange class 
/// is not instantiated prior to read of .rvv/.rve files
/// \param &W [out] Overlap/Exchange properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  COverlapExchangeClass::SetOEProperty(oe_struct &E, 
                                          string     param_name, 
                                          const string value)
{
  string name, val;
  name = StringToUppercase(param_name);
  val = StringToUppercase(value);

  if      (!name.compare("ET_TYPE"            )){E.oe_props.ET_type=value;}
  else if (!name.compare("D_TYPE"             )){E.oe_props.D_type=value;}
  else if (!name.compare("SA_TYPE"            )){E.oe_props.SA_type=value;}

  else{
    WriteWarning("COverlapExchangeClass::SetOEProperty: Unrecognized/invalid overlap/exchange parameter name ("+name+") in .rvv/.rve file",false);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the overlap/exchange property corresponding to param_name
/// \note This is declared as a static member because overlap/exchange class 
/// is not instantiated prior to read of .rvv/.rve files
/// \param &W [out] Overlap/Exchange properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  COverlapExchangeClass::SetOEProperty(oe_struct &E, 
                                          string     param_name, 
                                          const double value,
                                          int i, int j)
{
  string name;
  name = StringToUppercase(param_name);

  if (!name.compare("OVERLAP"            )){E.overlap[i][j]=value;}

  else{
    WriteWarning("COverlapExchangeClass::SetOEProperty: Unrecognized/invalid overlap/exchange parameter name ("+name+") in .rvv/.rve file",false);
  }
}
///////////////////////////////////////////////////////////////////////////
/// \brief Returns overlap/exchange property value corresponding to param_name
/// \param param_name [in] Parameter name
/// \return Overlap/Exchange properties corresponding to parameter name
//
double COverlapExchangeClass::GetOEProperty(string param_name) const
{
  return GetOEProperty(E,param_name);
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns overlap/exchange property value corresponding to param_name
/// \param param_name [in] Parameter name
/// \return Overlap/Exchange properties corresponding to parameter name
//
string COverlapExchangeClass::GetOEProperty(string param_name, string param_name1, int i)
{
  return GetOEProperty(E,param_name, param_name1, 0);
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns overlap/exchange property value corresponding to param_name
/// \param param_name [in] Parameter name
/// \return Overlap/Exchange properties corresponding to parameter name
//
double COverlapExchangeClass::GetOEProperty(string param_name, int i, int j)
{
  return GetOEProperty(E,param_name, i, j);
}
///////////////////////////////////////////////////////////////////////////
/// \brief Returns overlap/exchange property value corresponding to param_name
/// \param param_name [in] Parameter name
/// \return Overlap/Exchange properties corresponding to parameter name
//
double COverlapExchangeClass::GetOEProperty(string param_name, int i)
{
  return GetOEProperty(E,param_name, i);
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns overlap/exchange property value corresponding to param_name from structure provided
/// \param oe_struct [in] Overlap/Exchange structure
/// \param param_name [in] Parameter name
/// \return Overlap/Exchange properties corresponding to parameter name
//
double COverlapExchangeClass::GetOEProperty(const oe_struct &E, string param_name)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("NHRU"               )){return E.nHRU;}
  else if (!name.compare("NGWCELLS"           )){return E.nGWCells;}
  else if (!name.compare("CELL"               )){return E.cell;}
  
  else{
    string msg="COverlapExchangeClass::GetOEProperty: Unrecognized/invalid overlap/exchange parameter name in .rvv/.rve file: "+name;
    ExitGracefully(msg.c_str(),BAD_DATA);
    return 0.0;
  }
}
///////////////////////////////////////////////////////////////////////////
/// \brief Returns overlap/exchange property value corresponding to param_name from structure provided
/// \param oe_struct [in] Overlap/Exchange structure
/// \param param_name [in] Parameter name
/// \return Overlap/Exchange properties corresponding to parameter name
//
double COverlapExchangeClass::GetOEProperty(const oe_struct &E, string param_name, int i)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("ET_A"               )){return E.oe_props.ET_a[i];}
  else if (!name.compare("ET_B"               )){return E.oe_props.ET_b[i];}
  else if (!name.compare("D_A"                )){return E.oe_props.D_a[i];}
  else if (!name.compare("D_B"                )){return E.oe_props.D_b[i];}
  else if (!name.compare("SA_A"               )){return E.oe_props.SA_a[i];}
  else if (!name.compare("SA_B"               )){return E.oe_props.SA_b[i];}
  
  else{
    string msg="COverlapExchangeClass::GetOEProperty: Unrecognized/invalid overlap/exchange parameter name in .rvv/.rve file: "+name;
    ExitGracefully(msg.c_str(),BAD_DATA);
    return 0.0;
  }
}
///////////////////////////////////////////////////////////////////////////
/// \brief Returns overlap/exchange property value corresponding to param_name from structure provided
/// \param oe_struct [in] Overlap/Exchange structure
/// \param param_name [in] Parameter name
/// \return Overlap/Exchange properties corresponding to parameter name
//
string COverlapExchangeClass::GetOEProperty(const oe_struct &E, string param_name, string param_name1, int i)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("ET_TYPE"            )){return E.oe_props.ET_type;}
  else if (!name.compare("D_TYPE"             )){return E.oe_props.D_type;}
  else if (!name.compare("SA_TYPE"            )){return E.oe_props.SA_type;}
  
  else{
    string msg="COverlapExchangeClass::GetOEProperty: Unrecognized/invalid overlap/exchange parameter name in .rvv/.rve file: "+name;
    ExitGracefully(msg.c_str(),BAD_DATA);
    return "";
  }
}
///////////////////////////////////////////////////////////////////////////
/// \brief Returns overlap/exchange property value corresponding to param_name from structure provided
/// \param oe_struct [in] Groundwater structure
/// \param param_name [in] Parameter name
/// \return Overlap/Exchange properties corresponding to parameter name
//
double COverlapExchangeClass::GetOEProperty(const oe_struct &E, string param_name, int i, int j)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("OVERLAP"               )){return E.overlap[i][j];}

  else{
    string msg="COverlapExchangeClass::GetOEProperty: Unrecognized/invalid overlap/exchange parameter name in .rvv/.rve file: "+name;
    ExitGracefully(msg.c_str(),BAD_DATA);
    return 0.0;
  }
}