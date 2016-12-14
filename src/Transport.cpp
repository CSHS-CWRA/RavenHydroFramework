/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
------------------------------------------------------------------
	Master Transport/Tracer class
  coordinates information about constituent storage
----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "Model.h"
#include "Transport.h"

string FilenamePrepare(string filebase, const optStruct &Options); //Defined in StandardOutput.cpp

//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Transport constructor
/// \param pModel [in] Model object
//
CTransportModel::CTransportModel(CModel *pMod)
{
  pModel=pMod;

  nAdvConnections=0;
  iFromWater=NULL;
  iToWater=NULL;
  js_indices=NULL;

  nWaterCompartments=0;
  iWaterStorage=NULL;

  nConstituents=0;
  pConstituents=NULL;
  pConstitParams=NULL;

  pSources=NULL;
  nSources=0;

  aIndexMapping=NULL;
  aSourceIndices=NULL;

  pTransModel=this;

  aMinHist =NULL;
  aMlatHist=NULL;
  aMout    =NULL;

}
//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Transport destructor
//
CTransportModel::~CTransportModel()
{
  delete [] iFromWater;     iFromWater=NULL;
  delete [] iToWater;       iToWater=NULL;
  delete [] js_indices;     js_indices=NULL;
  delete [] iWaterStorage;  iWaterStorage=NULL;
  delete [] pConstituents;  pConstituents=NULL;
  delete [] pConstitParams; pConstitParams=NULL;
  delete [] aIndexMapping;  aIndexMapping=NULL;
  for (int i=0;i<nSources;i++     ){delete pSources[i];         } delete [] pSources;
  if (aSourceIndices!=NULL){
  for (int c=0;c<nConstituents;c++){delete [] aSourceIndices[c];} delete [] aSourceIndices;
  }
  DeleteRoutingVars();
}

//Static declaration
CTransportModel* CTransportModel::pTransModel=NULL;

//////////////////////////////////////////////////////////////////
/// \brief converts model layer index (i.e., m in CONSTITUENT[m]) to 
/// local constituent index c and index of corresponding local constiuent storage index j
//
void CTransportModel::m_to_cj(const int m, int &c, int &j) const
{
  ExitGracefullyIf((m<0) || (m>nConstituents*nWaterCompartments),
    "CTransportModel::LayerIndToConstitData: invalid layer index",BAD_DATA);

  j=m % nWaterCompartments;
  c=(m-j) / nWaterCompartments;
}

//////////////////////////////////////////////////////////////////
/// \brief returns layer index m corresponding to constituent c in water storage unit i_stor
//
int CTransportModel::GetLayerIndex(const int c, const int i_stor) const 
{
  int j=aIndexMapping[i_stor];
  if (j==DOESNT_EXIST){return DOESNT_EXIST;}
  if (c==DOESNT_EXIST){return DOESNT_EXIST;}
  /*ExitGracefullyIf(j==DOESNT_EXIST,
    "CTransportModel::GetLayerIndex: constituent storage unit not found. Invalid index passed",RUNTIME_ERR);*/
  //cout<<" layer index = "<<c*nWaterCompartments+j<<" istor: "<<i_stor<<" j: "<<j<<endl;
  return c*nWaterCompartments+j;
}

//////////////////////////////////////////////////////////////////
/// \brief returns layer index m corresponding to a specific constituent storage compartment (static version)
/// \input name in format !Nitrogen|SOIL[2]
/// \input comp_m layer index of storage compartment (e.g., 2 in above example)
//
int CTransportModel::GetLayerIndexFromName(const string name,const int comp_m) //static
{
  return pTransModel->GetLayerIndexFromName2(name,comp_m);
}

//////////////////////////////////////////////////////////////////
/// \brief returns layer index m corresponding to a specific constituent storage compartment (non-static version)
/// \input name in format !Nitrogen|SOIL ([2] already trimmed at this point)
/// \input comp_m layer index of storage compartment (e.g., 2 in above example)
//
int CTransportModel::GetLayerIndexFromName2(const string name,const int comp_m) const //local
{
  string constituent_name;
  string compartment_name;
  string tmp=name;
  int k;

  if (name.substr(0,1)!="!"){return DOESNT_EXIST;}//bad format, no leading "!"

  tmp=name.substr(1,name.length()-1);             //trim leading '!'
  k=(int)(tmp.find_first_of("|"));                //find char index of "|"
  if (k==DOESNT_EXIST){return DOESNT_EXIST;}      //bad format, no "|"
  
  constituent_name=tmp.substr(0,k);
  compartment_name=tmp.substr(k+1,tmp.length()-k);
  //cout<<":GetLayerIndexFromName2:  "<<name<<" -> "<<tmp<<" -> |"<<constituent_name<<"| + |"<<compartment_name<<"|"<<endl;

  int     cc,layer_index,i_stor;
  sv_type typ;
  cc    =pTransModel->GetConstituentIndex(constituent_name);
  typ   =CStateVariable::StringToSVType(compartment_name,layer_index,true);
  i_stor=pModel->GetStateVarIndex(typ,comp_m);

  return GetLayerIndex(cc,i_stor);
}

//////////////////////////////////////////////////////////////////
/// \brief returns full name of constitutent e.g.,  "Nitrogen in Soil Water[2]"
/// static version allowing calls without passing transport model
//
string CTransportModel::GetConstituentLongName(const int m)
{
  return pTransModel->GetConstituentName(m);
}

//////////////////////////////////////////////////////////////////
/// \brief returns full name of constitutent e.g., "Nitrogen in Soil Water[2]"
//
string CTransportModel::GetConstituentName(const int m) const
{
  int c,j;
  m_to_cj(m,c,j);
  sv_type typ=pModel->GetStateVarType(iWaterStorage[j]);
  int     ind=pModel->GetStateVarLayer(iWaterStorage[j]);

  return pConstituents[c]->name+" in "+CStateVariable::GetStateVarLongName(typ,ind);
}
//////////////////////////////////////////////////////////////////
/// \brief returns name of constitutent type e.g., "Nitrogen"
/// \remark static routine 
/// \param c index of constituent storage 
//
string CTransportModel::GetConstituentTypeName2(const int c) 
{
  return pTransModel->GetConstituent(c)->name;
}
//////////////////////////////////////////////////////////////////
/// \brief returns global sv index of water storage unit corresponding to CONSTITUENT[m]
//
int CTransportModel::GetWaterStorIndexFromLayer(const int m) const
{
  int c,j;
  m_to_cj(m,c,j);
  return iWaterStorage[j];
}
//////////////////////////////////////////////////////////////////
/// \brief returns full name of constitutent e.g., "!Nitrogen|SOIL[2]"
//
string CTransportModel::GetConstituentShortName(const int m) const
{
  int c,j;
  m_to_cj(m,c,j);
  sv_type typ=pModel->GetStateVarType(iWaterStorage[j]);
  int     ind=pModel->GetStateVarLayer(iWaterStorage[j]);

  return "!"+pConstituents[c]->name+"|"+CStateVariable::SVTypeToString(typ,ind);
}

//////////////////////////////////////////////////////////////////
/// \brief returns name of constitutent type e.g., "Nitrogen"
/// \remark static routine 
/// \param m layer index of constituent storage 
//
string CTransportModel::GetConstituentTypeName(const int m) 
{
  int c,j;
  pTransModel->m_to_cj(m,c,j);
  return pTransModel->GetConstituent(c)->name;
}

//////////////////////////////////////////////////////////////////
/// \brief returns number of water compartments in model
//
int    CTransportModel::GetNumWaterCompartments() const{return nWaterCompartments;}

//////////////////////////////////////////////////////////////////
/// \brief returns number of water transport connections in model
//
int    CTransportModel::GetNumAdvConnections() const {return nAdvConnections;}

//////////////////////////////////////////////////////////////////
/// \brief returns number of constituents transported in model
//
int    CTransportModel::GetNumConstituents() const{return nConstituents;}

//////////////////////////////////////////////////////////////////
/// \brief returns constituent c in model 
//
const constituent *CTransportModel::GetConstituent(const int c) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((c<0)||(c>=nConstituents),
    "CTransportModel::GetConstituent: invalid index",BAD_DATA);
#endif
  return pConstituents[c];
}
//////////////////////////////////////////////////////////////////
/// \brief returns constituent c in model 
//
const transport_params *CTransportModel::GetConstituentParams(const int c) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((c<0) || (c >= nConstituents),
    "CTransportModel::GetConstituent: invalid index", BAD_DATA);
#endif
  return pConstitParams[c];
}

//////////////////////////////////////////////////////////////////
/// \brief returns index c of constituents transported in model
//
int    CTransportModel::GetConstituentIndex(const string name) const
{
  for (int c=0; c<nConstituents;c++){
    if (StringToUppercase(pConstituents[c]->name)==StringToUppercase(name)){return c;}
  }
  return DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of "from" constituent mass compartment
/// \param c [in] constituent index
/// \param q [in] local index of connection
//
int    CTransportModel::GetFromIndex(const int c,const int q) const
{
  int j=aIndexMapping[iFromWater[q]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*nWaterCompartments+j);
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of "to" constituent mass compartment
/// \param c [in] constituent index
/// \param q [in] local index of connection
//
int    CTransportModel::GetToIndex  (const int c,const int q) const
{
  int j=aIndexMapping[iToWater[q]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*nWaterCompartments+j);
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index of  constituent mass compartment
/// \param c [in] constituent index
/// \param ii [in] local index of water storage unit
//
int    CTransportModel::GetStorIndex  (const int c,const int ii) const
{
  int j=aIndexMapping[iWaterStorage[ii]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*nWaterCompartments+j);
}

//////////////////////////////////////////////////////////////////
/// \brief returns state variable index i of "from" water compartment
/// \param q [in] local index of connection
//
int    CTransportModel::GetFromWaterIndex(const int q) const
{
  return iFromWater[q]; 
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of "to" water compartment
/// \param q [in] index of connection
//
int    CTransportModel::GetToWaterIndex  (const int q) const
{
  return iToWater[q]; 
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of water compartment ii
/// \param ii [in] index of water storage (from 0 to nWaterCompartments-1)
//
int    CTransportModel::GetStorWaterIndex  (const int ii) const
{
  return iWaterStorage[ii]; 
}

//////////////////////////////////////////////////////////////////
/// \brief returns master process index js of advection connection
/// \param q [in] local index of connection
//
int CTransportModel::GetJsIndex(const int q) const
{
  return js_indices[q];
}
//////////////////////////////////////////////////////////////////
/// \brief returns effective retardation factor for constituent c 
/// being transported from storage compartment iFromWater to storage compartment iToWater
/// \param c [in] constituent index
/// \param iFromWater [in] index of "from" water storage state variable
/// \param iToWater [in] index of "to" water storage state variable
//
double CTransportModel::GetRetardationFactor(const int c,const CHydroUnit	*pHRU, const  int iFromWater,const int iToWater) const
{
  sv_type fromType,toType;
  //int soil_ind;
  fromType=pModel->GetStateVarType(iFromWater);
  toType  =pModel->GetStateVarType(iToWater);
  int    m=pModel->GetStateVarLayer(iFromWater);
#ifndef _STRICTCHECK_
    ExitGracefullyIf(m==-1,"GetRetardationFactor:invalid iFromWater",RUNTIME_ERR);
#endif
  if (fromType==SOIL) 
  {
    if (toType!=ATMOSPHERE)
    {
      return pHRU->GetSoilProps(m)->retardation[c];
      //return pHRU->GetSoilTransportProps(m)->retardation[c];
    }
    else if (toType==ATMOSPHERE)
    {
      return 1.0;
      //return pHRU->GetVegProps()->uptake_fact[c];
    }
    return 1.0;
  }
  else 
  {
    return 1.0;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief returns effective retardation factor for constituent c 
/// being transported from storage compartment iFromWater to storage compartment iToWater
/// \param c [in] constituent index
/// \param iFromWater [in] index of "from" water storage state variable
/// \param iToWater [in] index of "to" water storage state variable
//
double CTransportModel::GetDecayCoefficient(const int c,const CHydroUnit	*pHRU,const int iStorWater) const
{
  sv_type storType;
  double decay_coeff = GetConstituentParams(c)->decay_coeff;

  storType=pModel->GetStateVarType(iStorWater);

  //add special decay_coefficients from other processes
  if (storType == SOIL)
  {
    int m = pModel->GetStateVarLayer(iStorWater);
#ifndef _STRICTCHECK_
    ExitGracefullyIf(m==-1,"GetDecayCoefficient:invalid iFromWater",RUNTIME_ERR);
#endif
    //decay_coeff += pHRU->GetSoilProps(m)->mineralization_rate[c];
    //decay_coeff += pHRU->GetSoilProps(m)->loss_rate[c]; //e.g., denitrification

  }
  return decay_coeff;
}
//////////////////////////////////////////////////////////////////
/// \brief adds new transportable constituent to model
/// \note adds corresponding state variables to model
/// \param name [in] name of constituent
/// \param is_tracer [in] boolean - true if this is a tracer
//
void   CTransportModel::AddConstituent(string name, bool is_tracer)
{
  constituent *pConstit=new constituent;

  //Initialize constituent members
  pConstit->name=name;
  pConstit->is_tracer=is_tracer;
  if (is_tracer){
    pConstit->can_evaporate=true;
  }
  pConstit->initial_mass=0;
  pConstit->cumul_input =0;
  pConstit->cumul_output=0;

  //add to master list of constituents
  if (!DynArrayAppend((void**&)(pConstituents),(void*)(pConstit),nConstituents)){
   ExitGracefully(" CTransportModel::AddConstituent: adding NULL constituent",BAD_DATA);} 

  int c=(nConstituents-1);//constit. index of current constituent

  //Add corresponding constituent storage state variables to model
  sv_type *aSV =new sv_type [nWaterCompartments];
  int     *aLev=new int     [nWaterCompartments];
  for (int j=0;j<nWaterCompartments; j++)
  { 
    int m=j+c*nWaterCompartments;
    aSV [j]=CONSTITUENT;
    aLev[j]=m;
  }
  pModel->AddStateVariables(aSV,aLev,nWaterCompartments);

  //Add corresponding constituent source & sink state variables to model
  aSV [0]=CONSTITUENT_SRC;
  aLev[0]=c;
  pModel->AddStateVariables(aSV,aLev,1);

  aSV [0]=CONSTITUENT_SW;
  aLev[0]=c;
  pModel->AddStateVariables(aSV,aLev,1);

  //Parameter initialization
  transport_params *pP = new transport_params;
  InitializeConstitParams(pP);
  int junk = nConstituents - 1;

  //add to master list of constituent parameters
  if (!DynArrayAppend((void**&)(pConstitParams), (void*)(pP), junk)){
    ExitGracefully(" CTransportModel::AddConstituent: adding NULL constituent parameter set", BAD_DATA);
  }

  delete [] aSV;
  delete [] aLev;

}
//////////////////////////////////////////////////////////////////
/// \brief Initialization of transport parameter structure
/// \param pP [in] valid pointer to transport_params structure pP
//
void CTransportModel::InitializeConstitParams(transport_params *pP)
{
  pP->decay_coeff = 0.0;
  for (int m = 0; m < MAX_SOIL_CLASSES; m++){
    pP->retardation[m] = 1.0;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Preparation of all transport variables
/// \note  called after all waterbearing state variables & processes have been generated, but before constiuents created
/// \param Options [in] Global model options structure
//
void CTransportModel::Prepare(const optStruct &Options)
{
  CHydroProcessABC *pProc;
  int js=0;
  
  //determine number of connections
  nAdvConnections=0;
  for (int j=0; j<pModel->GetNumProcesses(); j++)
  {
    pProc=pModel->GetProcess(j); 
    for (int q=0; q<pProc->GetNumConnections(); q++)
    {
      int iF=pProc->GetFromIndices()[q];
      int iT=pProc->GetToIndices()  [q];
      sv_type typF=pModel->GetStateVarType(iF);
      sv_type typT=pModel->GetStateVarType(iT);
      if ((iF!=iT) && (CStateVariable::IsWaterStorage(typF)) && (CStateVariable::IsWaterStorage(typT)))
      { //This is a process that may advect a constituent
        nAdvConnections++;
      }
    }
  }
  
  //determine iFromWater, iToWater indices and process index for each and every connection
  js=0;
  int qq=0;
  iFromWater=new int [nAdvConnections];
  iToWater  =new int [nAdvConnections];
  js_indices=new int [nAdvConnections];
  for (int j=0; j<pModel->GetNumProcesses(); j++)
  {
    pProc=pModel->GetProcess(j);
    
    for (int q=0; q<pProc->GetNumConnections(); q++)
    {
      int iF=pProc->GetFromIndices()[q];
      int iT=pProc->GetToIndices()  [q];
      sv_type typF=pModel->GetStateVarType(iF);
      sv_type typT=pModel->GetStateVarType(iT);
      if ((iF!=iT) && (CStateVariable::IsWaterStorage(typF)) && (CStateVariable::IsWaterStorage(typT)))
      { //This is a process that may advect a constituent
        iFromWater[qq]=iF;
        iToWater  [qq]=iT;
        js_indices[qq]=js;
        qq++;
      }
      else{
        //cout<<"NOT A WATER CONNECTION:"<<endl;
        //cout<<pProc->GetProcessType()<<"  from "<<CStateVariable::SVTypeToString(pModel->GetStateVarType(iF),0)<<" to "<<CStateVariable::SVTypeToString(pModel->GetStateVarType(iT),0)<<endl;
      }
      js++;
    }
  }
  
  //identify number of water storage compartments
  nWaterCompartments=0;
  for (int i=0;i<pModel->GetNumStateVars(); i++)
  {
    if (CStateVariable::IsWaterStorage(pModel->GetStateVarType(i)))
    {
      nWaterCompartments++;
    }
  }
  // identify all state variables which are water storage compartments
  iWaterStorage=new int [nWaterCompartments]; 
  aIndexMapping=new int [pModel->GetNumStateVars()];//current # of state variables does not include transport constituents
  int j=0; //j sifts through storage compartments j=0..nWaterCompartments
  for (int i=0;i<pModel->GetNumStateVars(); i++)
  {
    aIndexMapping[i]=DOESNT_EXIST;
    if (CStateVariable::IsWaterStorage(pModel->GetStateVarType(i)))
    {
      iWaterStorage[j]=i;
      aIndexMapping[i]=j;
      j++;
    }
  }
  
  /// \todo [QA/QC]: check for two constituents with same name?
  if (nConstituents==0){return;}/// all of the above work necessary even with no transport?

  //Synopsis
	if (!Options.silent){
	  cout<<"===TRANSPORT MODEL SUMMARY=============================="<<endl;
	  cout<<"   number of compartments: "<<nWaterCompartments<<endl;
	  cout<<"   number of constituents: "<<nConstituents<<endl;
	  cout<<"   number of connections:  "<<nAdvConnections<<endl;
	  cout<<"   number of source terms: "<<nSources<<endl;
	  // TMP DEBUG below===================================================
	  if (false){
	    for (int i=0;i<nAdvConnections;i++)
	    {
	      cout<<i<<": moves ";
	      cout<< " from "<<CStateVariable::SVTypeToString(pModel->GetStateVarType(iFromWater[i]),pModel->GetStateVarLayer(iFromWater[i]));
	      cout<< " to "  <<CStateVariable::SVTypeToString(pModel->GetStateVarType(iToWater  [i]),pModel->GetStateVarLayer(iToWater  [i]));
	      cout<<endl;
	    }
	     for (int i=0;i<pModel->GetNumStateVars(); i++)
	      {
	        if (aIndexMapping[i]!=DOESNT_EXIST){
	          cout<<"constituent "<< i;
	          cout<<" corresponds to "<<CStateVariable::SVTypeToString(pModel->GetStateVarType(aIndexMapping[i]),pModel->GetStateVarLayer(aIndexMapping[i]));
	          cout<<endl;
	        }
	        else
	        {
	          cout<<"SV "<<i<<" has no corresponding mapping"<<endl;
	        }
	     }
	     int c=0;
	     for (int j=0;j<nWaterCompartments;j++)
	     {
	       int m=c*nWaterCompartments+j;
	       cout<<"water compartment "<< CStateVariable::SVTypeToString(pModel->GetStateVarType(iWaterStorage[j]),pModel->GetStateVarLayer(iWaterStorage[j]));
	       cout<<" corresponds to "<<CStateVariable::SVTypeToString(CONSTITUENT,m);
	       cout<<endl;
	    }
	  }
	  // ====================================================================
	  cout<<"========================================================"<<endl;
	}
}

//////////////////////////////////////////////////////////////////
/// \brief adds dirichlet source
/// \param const_name [in] constituent name
/// \param i_stor [in] global index of water storage state variable
/// \param kk [in] HRU group index (or -1 if this applies to all HRUs)
/// \param Cs [in] Dirichlet source concentration [mg/m2]
//
void   CTransportModel::AddDirichletCompartment(const string const_name, const int i_stor, const int kk, const double Cs)
{
  static constit_source *pLast;
  constit_source *pSource=new constit_source();
  
  pSource->constit_index=DOESNT_EXIST;
  for (int c=0; c<nConstituents;c++){
    if (StringToUppercase(const_name)==StringToUppercase(pConstituents[c]->name)){
      pSource->constit_index =c;
    }
  }
  pSource->dirichlet    =true;
  pSource->concentration=Cs;
  pSource->flux         =0.0;
  pSource->i_stor       =i_stor;
  pSource->kk           =kk;
  pSource->pTS          =NULL;

  ExitGracefullyIf(pSource->constit_index==DOESNT_EXIST,"AddDirichletCompartment: invalid constituent name",BAD_DATA_WARN);
  ExitGracefullyIf(pSource->i_stor       ==DOESNT_EXIST,"AddDirichletCompartment: invalid storage compartment index",BAD_DATA_WARN);

  if (!DynArrayAppend((void**&)(pSources),(void*)(pSource),nSources)){
   ExitGracefully("CTransportModel::AddDirichletCompartment: adding NULL source",BAD_DATA);} 

  pLast=pSource;//so source is not deleted upon leaving this routine
}
//////////////////////////////////////////////////////////////////
/// \brief adds dirichlet source time series
/// \param const_name [in] constituent name
/// \param i_stor [in] global index of water storage state variable
/// \param kk [in] HRU group index (or -1 if this applies to all HRUs)
/// \param pTS [in] Time series of Dirichlet source concentration [mg/m2]
//
void   CTransportModel::AddDirichletTimeSeries(const string const_name, const int i_stor, const int kk, const CTimeSeries *pTS)
{
  static constit_source *pLast;
  constit_source *pSource=new constit_source();
  pSource->constit_index=DOESNT_EXIST;
  for (int c=0; c<nConstituents;c++){
    if (StringToUppercase(const_name)==StringToUppercase(pConstituents[c]->name)){
      pSource->constit_index =c;
    }
  }
  pSource->dirichlet    =true;
  pSource->concentration=DOESNT_EXIST;
  pSource->flux         =0.0;
  pSource->i_stor       =i_stor;
  pSource->kk           =kk;
  pSource->pTS          =pTS;

  ExitGracefullyIf(pSource->constit_index==DOESNT_EXIST,("AddDirichletTimeSeries: invalid constituent name "+const_name).c_str(),BAD_DATA_WARN);
  ExitGracefullyIf(pSource->i_stor       ==DOESNT_EXIST,"AddDirichletTimeSeries: invalid storage compartment index",BAD_DATA_WARN);

  if (!DynArrayAppend((void**&)(pSources),(void*)(pSource),nSources)){
   ExitGracefully("CTransportModel::AddDirichletCompartment: adding NULL source",BAD_DATA);} 

  pLast=pSource; //so source is not deleted upon leaving this routine
}

//////////////////////////////////////////////////////////////////
/// \brief adds influx (Neumann) source
/// \param const_name [in] constituent name
/// \param i_stor [in] global index of water storage state variable
/// \param kk [in] HRU group index (or -1 if this applies to all HRUs)
/// \param flux [in] specified fixed mass influx rate [mg/m2/d]
//
void   CTransportModel::AddInfluxSource(const string const_name, const int i_stor, const int kk, const double flux)
{
  static constit_source *pLast;
  constit_source *pSource=new constit_source();
  
  pSource->constit_index=DOESNT_EXIST;
  for (int c=0; c<nConstituents;c++){
    if (StringToUppercase(const_name)==StringToUppercase(pConstituents[c]->name)){
      pSource->constit_index =c;
    }
  }
  pSource->dirichlet    =false;
  pSource->concentration=0.0;
  pSource->flux         =flux;
  pSource->i_stor       =i_stor;
  pSource->kk           =kk;
  pSource->pTS          =NULL;

  ExitGracefullyIf(pSource->constit_index==DOESNT_EXIST,"AddInfluxSource: invalid constituent name",BAD_DATA_WARN);
  ExitGracefullyIf(pSource->i_stor       ==DOESNT_EXIST,"AddInfluxSource: invalid storage compartment index",BAD_DATA_WARN);

  if (!DynArrayAppend((void**&)(pSources),(void*)(pSource),nSources)){
   ExitGracefully("CTransportModel::AddInfluxSource: adding NULL source",BAD_DATA);} 

  pLast=pSource;//so source is not deleted upon leaving this routine
}
//////////////////////////////////////////////////////////////////
/// \brief adds influx (Neumann) source source time series
/// \param const_name [in] constituent name
/// \param i_stor [in] global index of water storage state variable
/// \param kk [in] HRU group index (or -1 if this applies to all HRUs)
/// \param pTS [in] Time series of Dirichlet source flux rate [mg/m2/d]
//
void   CTransportModel::AddInfluxTimeSeries(const string const_name, const int i_stor, const int kk, const CTimeSeries *pTS)
{
  static constit_source *pLast;
  constit_source *pSource=new constit_source();
  pSource->constit_index=DOESNT_EXIST;
  for (int c=0; c<nConstituents;c++){
    if (StringToUppercase(const_name)==StringToUppercase(pConstituents[c]->name)){
      pSource->constit_index =c;
    }
  }
  pSource->dirichlet    =false;
  pSource->concentration=0.0;
  pSource->flux         =DOESNT_EXIST;
  pSource->i_stor       =i_stor;
  pSource->kk           =kk;
  pSource->pTS          =pTS;

  ExitGracefullyIf(pSource->constit_index==DOESNT_EXIST,("AddDirichletTimeSeries: invalid constituent name "+const_name).c_str(),BAD_DATA_WARN);
  ExitGracefullyIf(pSource->i_stor       ==DOESNT_EXIST,"AddDirichletTimeSeries: invalid storage compartment index",BAD_DATA_WARN);

  if (!DynArrayAppend((void**&)(pSources),(void*)(pSource),nSources)){
   ExitGracefully("CTransportModel::AddDirichletCompartment: adding NULL source",BAD_DATA);} 

  pLast=pSource; //so source is not deleted upon leaving this routine
}

//////////////////////////////////////////////////////////////////
/// \brief Initialization of all transport variables
/// \note determines initial conditions for all constituents, initializes routing variables
/// \note  called after all constituents have been added by CModel::Initialize
//
void CTransportModel::Initialize()
{
  double area=pModel->GetWatershedArea();
  for (int c=0;c<nConstituents;c++)
  {
    pConstituents[c]->cumul_input=0;
    pConstituents[c]->cumul_output=0;
    for (int ii=0;ii<nWaterCompartments;ii++)
    {
      int m=c*nWaterCompartments+ii;
      int i=pModel->GetStateVarIndex(CONSTITUENT,m);

      // zero out initial mass in all storage units with zero volume (to avoid absurdly high concentrations)
      double watstor;
      int i_stor=GetStorWaterIndex(ii);
      for (int k = 0; k < pModel->GetNumHRUs(); k++){
        watstor = pModel->GetHydroUnit(k)->GetStateVarValue(i_stor);
        if (watstor<1e-9){pModel->GetHydroUnit (k)->SetStateVarValue(i,0.0);}
      }
      //update initial model mass
      pConstituents[c]->initial_mass+=pModel->GetAvgStateVar(i)*(area*M2_PER_KM2);
    }
  }
  


  /// \todo[funct]: calculate initial mass from Dirichlet cells 


  //populate array of source indices
  // \todo [funct] will have to revise to support different sources in different HRUs (e.g., aSourceIndices[c][i_stor][k])
  aSourceIndices=NULL;
  aSourceIndices=new int *[nConstituents];
  ExitGracefullyIf(aSourceIndices==NULL,"CTransport::Initialize",OUT_OF_MEMORY);
  for (int c=0;c<nConstituents;c++)
  {
    aSourceIndices[c]=NULL;
    aSourceIndices[c]=new int [pModel->GetNumStateVars()];
    ExitGracefullyIf(aSourceIndices[c]==NULL,"CTransport::Initialize",OUT_OF_MEMORY);
    for (int i_stor=0;i_stor<pModel->GetNumStateVars();i_stor ++){
      aSourceIndices[c][i_stor]=DOESNT_EXIST;
      for (int i=0;i<nSources;i++){
        if ((pSources[i]->i_stor==i_stor) && (pSources[i]->constit_index==c))
        {
          if (aSourceIndices[c][i_stor] != DOESNT_EXIST){
            WriteWarning("CTransportModel::Intiialize: cannot currently have more than one constitutent source per constituent/storage combination",false);
          }
          aSourceIndices[c][i_stor]=i; //each 
        }
      }
    }
  }

  InitializeRoutingVars();
}
//////////////////////////////////////////////////////////////////
/// \brief Increment cumulative input of mass to watershed.
/// \note called within solver to track mass balance
/// \param Options [in] Global model options structure
/// \param tt [in] current model time structure
//
void CTransportModel::IncrementCumulInput (const optStruct &Options, const time_struct &tt)
{
  for (int c=0;c<nConstituents;c++)
  {
    pConstituents[c]->cumul_input+=0;// \todo [funct]: increment cumulative input [mg]
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Increment cumulative output of mass from watershed.
/// \note called within solver to track mass balance
/// \param Options [in] Global model options structure
//
void CTransportModel::IncrementCumulOutput(const optStruct &Options)
{
  for (int c=0;c<nConstituents;c++)
  {
    for (int p=0;p<pModel->GetNumSubBasins();p++)
    {
      if (pModel->GetSubBasin(p)->GetDownstreamID()==DOESNT_EXIST)//outlet does not drain into another subbasin
      { 
        pConstituents[c]->cumul_output+=GetIntegratedMassOutflow(p,c)*Options.timestep;
      }
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Test for whether a dirichlet condition applies to a certain compartment, place, and time
/// \note called within solver to track mass balance
/// \returns true if dirichlet source applies
/// \returns Cs, source concentration
/// \param i_stor [in] storage index of water compartment
/// \param c [in] constituent index
/// \param k [in] global HRU index
/// \param tt [in] current time structure
/// \param Cs [out] Dirichlet source concentration
//
bool  CTransportModel::IsDirichlet(const int i_stor, const int c, const int k, const time_struct tt, double &Cs) const
{
  Cs=0.0;
  
  int i_source=aSourceIndices[c][i_stor];
  if (i_source!=DOESNT_EXIST) {
    if (!pSources[i_source]->dirichlet){return false;}
    if (pSources[i_source]->kk==DOESNT_EXIST)
    {
      Cs = pSources[i_source]->concentration; 
      if (Cs != DOESNT_EXIST){return true;}
      else{//time series 
        Cs = pSources[i_source]->pTS->GetValue(tt.model_time );
        return true;
      }
    }
    else { //Check if we are in HRU Group
      if (pModel->GetHRUGroup(pSources[i_source]->kk)->IsInGroup(k))
      {
        Cs=pSources[i_source]->concentration; 
        if (Cs != DOESNT_EXIST){return true;}
        else{//time series 
          Cs = pSources[i_source]->pTS->GetValue(tt.model_time );
          return true;
        }
      }
    }
  }
  return false;
}
//////////////////////////////////////////////////////////////////
/// \brief returns specified mass flux for given constitutent and water storage unit at time tt
/// \returns source flux in mg/m2/d
/// \param i_stor [in] storage index of water compartment
/// \param c [in] constituent index
/// \param k [in] global HRU index
/// \param tt [in] current time structure
//
double  CTransportModel::GetSpecifiedMassFlux(const int i_stor, const int c, const int k, const time_struct tt) const
{
  double flux;
  bool retrieve=false;
  int i_source=aSourceIndices[c][i_stor];
  if (i_source == DOESNT_EXIST) {return 0.0;}
  if (pSources[i_source]->dirichlet){return 0.0;}

  if (pSources[i_source]->kk==DOESNT_EXIST)//no HRU gorup
  {
    retrieve=true;
  }
  else { //Check if we are in HRU Group
    retrieve=pModel->GetHRUGroup(pSources[i_source]->kk)->IsInGroup(k);
  }
  if (retrieve)
  {
    flux=pSources[i_source]->concentration; 
    if (flux == DOESNT_EXIST){//get from time series
      flux=pSources[i_source]->pTS->GetValue(tt.model_time );
    }
    return flux;
  }
  else {return 0.0;}

}
//////////////////////////////////////////////////////////////////
/// \brief Write transport output file headers
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CTransportModel::WriteOutputFileHeaders(const optStruct &Options) const
{
  string filename;
  ofstream OUT;

  int iCumPrecip=pModel->GetStateVarIndex(ATMOS_PRECIP);

  string kg,mgL,kgd;
 
  for (int c=0;c<nConstituents;c++)
  {
    
    //units names
    kg="[kg]"; kgd="[kg/d]"; mgL="[mg/l]";
    //kg="[mg/m2]"; kgd="[mg/m2/d]"; mgL="[mg/m2]";//TMP DEBUG OUTPUT OVERRIDE
    if (pConstituents[c]->is_tracer){
      kg="[-]"; kgd="[-]"; mgL="[-]";
    }
    
    //Concentrations file
    //--------------------------------------------------------------------
    filename=pConstituents[c]->name+"Concentrations.csv";
    filename=FilenamePrepare(filename,Options);

    pConstituents[c]->OUTPUT.open(filename.c_str());
    if (pConstituents[c]->OUTPUT.fail()){
      ExitGracefully(("CTransportModel::WriteOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
    }

    pConstituents[c]->OUTPUT<<"time[d],date,hour,influx"<<kgd<<",Channel Storage"<<kg<<",Rivulet Storage"<<kg;  
    for (int i=0;i<pModel->GetNumStateVars();i++)
    {
      if ((CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))) && (i!=iCumPrecip)){
         pConstituents[c]->OUTPUT<<","<<CStateVariable::GetStateVarLongName(pModel->GetStateVarType(i),pModel->GetStateVarLayer(i))<<" "<<mgL;
      }
    }
    pConstituents[c]->OUTPUT<<", Total Mass "<<kg<<", Cum. Loading "<<kg<<", Cum. Mass Lost "<<kg<<", MB Error "<<kg<<endl; 

    //Pollutograph file
    //--------------------------------------------------------------------
    filename=pConstituents[c]->name+"Pollutographs.csv";
    filename=FilenamePrepare(filename,Options);

    pConstituents[c]->POLLUT.open(filename.c_str());
    if (pConstituents[c]->POLLUT.fail()){
      ExitGracefully(("CTransportModel::WriteOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    pConstituents[c]->POLLUT<<"time[d],date,hour";
    const CSubBasin *pBasin;
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      pBasin=pModel->GetSubBasin(p);
      if (pBasin->IsGauged()){
        string name;
        if (pBasin->GetName()==""){pConstituents[c]->POLLUT<<",ID="<<pBasin->GetID()  <<" "<<mgL;}
        else                      {pConstituents[c]->POLLUT<<","   <<pBasin->GetName()<<" "<<mgL;}
      }
    }
    pConstituents[c]->POLLUT<<endl;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Write transport output file headers in .tb0 format
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CTransportModel::WriteEnsimOutputFileHeaders(const optStruct &Options) const
{

  string filename;
  ofstream OUT;

  int iCumPrecip=pModel->GetStateVarIndex(ATMOS_PRECIP);

  string kg,mgL,kgd;
 
  int i;
  time_struct tt,tt2;
  JulianConvert(0.0,Options.julian_start_day,Options.julian_start_year,tt);
  JulianConvert(Options.timestep, Options.julian_start_day, Options.julian_start_year, tt2);//end of the timestep

  for (int c=0;c<nConstituents;c++)
  {
    //units names
    kg="kg"; kgd="kg/d"; mgL="mg/l";
    //kg="mg/m2"; kgd="mg/m2/d"; mgL="mg/m2";//TMP DEBUG OUTPUT OVERRIDE
    if (pConstituents[c]->is_tracer){
      kg="none"; kgd="none"; mgL="none";
    }

    //Concentrations file
    //--------------------------------------------------------------------
    filename=pConstituents[c]->name+"Concentrations.tb0";
    filename=FilenamePrepare(filename,Options);

    pConstituents[c]->OUTPUT.open(filename.c_str());
    if (pConstituents[c]->OUTPUT.fail()){
      ExitGracefully(("CTransportModel::WriteEnsimOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    pConstituents[c]->OUTPUT<<"#########################################################################"<<endl;
    pConstituents[c]->OUTPUT<<":FileType tb0 ASCII EnSim 1.0"<<endl;
    pConstituents[c]->OUTPUT<<"#"<<endl;
    pConstituents[c]->OUTPUT<<":Application   Raven"<<endl;
    pConstituents[c]->OUTPUT<<":Version       "<<Options.version<<endl;
    pConstituents[c]->OUTPUT<<":CreationDate  "<<GetCurrentTime()<<endl;
    pConstituents[c]->OUTPUT<<"#"<<endl;
    pConstituents[c]->OUTPUT<<"#------------------------------------------------------------------------"<<endl;
    pConstituents[c]->OUTPUT<<"#"<<endl;
    pConstituents[c]->OUTPUT<<":RunName       "<<Options.run_name<<endl;
    pConstituents[c]->OUTPUT<<":Format        Instantaneous" << endl;
    pConstituents[c]->OUTPUT<<"#"<<endl;

    if (Options.suppressICs){
    pConstituents[c]->OUTPUT<< ":StartTime " << tt2.date_string << " " << DecDaysToHours(tt2.julian_day) << endl;
    }
    else{
    pConstituents[c]->OUTPUT<< ":StartTime " << tt.date_string << " " << DecDaysToHours(tt.julian_day) << endl;
    }

    if (Options.timestep!=1.0){pConstituents[c]->OUTPUT<<":DeltaT " <<DecDaysToHours(Options.timestep)<<endl;}
    else                      {pConstituents[c]->OUTPUT<<":DeltaT 24:00:00.00"  <<endl;}
    pConstituents[c]->OUTPUT<<"#"<<endl;

    pConstituents[c]->OUTPUT<<":ColumnMetaData"<<endl;
    pConstituents[c]->OUTPUT<<"  :ColumnName influx \"Channel storage\" \"Rivulet storage\"";
    for (i=0;i<pModel->GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))) && (i!=iCumPrecip)){
          pConstituents[c]->OUTPUT<<" \""<<CStateVariable::GetStateVarLongName(pModel->GetStateVarType(i),pModel->GetStateVarLayer(i))<<"\"";}}
    pConstituents[c]->OUTPUT<<" \"Total Mass\" \"Cum. Loading\" \"Cum. Mass lost\" \"MB error\""<<endl;

    pConstituents[c]->OUTPUT<<"  :ColumnUnits "<<kgd<<" "<<kg <<" "<< kg; 
    for (i=0;i<pModel->GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))) && (i!=iCumPrecip)){
          pConstituents[c]->OUTPUT<<" mgL";}}
    pConstituents[c]->OUTPUT<<" "<<kg<<" "<<kg<<" "<<kg<<" "<<kg<<endl;

    pConstituents[c]->OUTPUT<<"  :ColumnType float float float";
    for (i=0;i<pModel->GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))) && (i!=iCumPrecip)){
          pConstituents[c]->OUTPUT<<" float";}}
    pConstituents[c]->OUTPUT<<" float float float float"<<endl;

    pConstituents[c]->OUTPUT<<"  :ColumnFormat -1 0 0";
    for (i=0;i<pModel->GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))) && (i!=iCumPrecip)){
          pConstituents[c]->OUTPUT<<" 0";}}
    pConstituents[c]->OUTPUT<<" 0 0 0 0"<<endl;

    pConstituents[c]->OUTPUT<<":EndColumnMetaData"<<endl;
    pConstituents[c]->OUTPUT<<":EndHeader"<<endl;

    //Pollutograph file
    //--------------------------------------------------------------------
    filename=pConstituents[c]->name+"Pollutographs.tb0";
    filename=FilenamePrepare(filename,Options);

    pConstituents[c]->POLLUT.open(filename.c_str());
    if (pConstituents[c]->OUTPUT.fail()){
      ExitGracefully(("CTransportModel::WriteEnsimOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
	  pConstituents[c]->POLLUT<<"#########################################################################"<<endl;
    pConstituents[c]->POLLUT<<":FileType tb0 ASCII EnSim 1.0"<<endl;
    pConstituents[c]->POLLUT<<"#"<<endl;
    pConstituents[c]->POLLUT<<":Application   Raven"<<endl;
    pConstituents[c]->POLLUT<<":Version       "<<Options.version<<endl;
	  pConstituents[c]->POLLUT<<":CreationDate  "<<GetCurrentTime()<<endl;
	  pConstituents[c]->POLLUT<<"#"<<endl;
	  pConstituents[c]->POLLUT<<"#------------------------------------------------------------------------"<<endl;
	  pConstituents[c]->POLLUT<<"#"<<endl;
    pConstituents[c]->POLLUT<<":RunName       "<<Options.run_name<<endl;
    pConstituents[c]->POLLUT<<":Format        Instantaneous" << endl;
	  pConstituents[c]->POLLUT<<"#"<<endl;

    if (Options.suppressICs){
    pConstituents[c]->POLLUT<< ":StartTime " << tt2.date_string << " " << DecDaysToHours(tt2.julian_day) << endl;
    }
    else{
    pConstituents[c]->POLLUT<< ":StartTime " << tt.date_string << " " << DecDaysToHours(tt.julian_day) << endl;
    }

    if (Options.timestep!=1.0){pConstituents[c]->POLLUT<<":DeltaT " <<DecDaysToHours(Options.timestep)<<endl;}
    else                      {pConstituents[c]->POLLUT<<":DeltaT 24:00:00.00"  <<endl;}
    pConstituents[c]->POLLUT<<"#"<<endl;

    const CSubBasin *pBasin;

    pConstituents[c]->POLLUT<<":ColumnMetaData"<<endl;
    pConstituents[c]->POLLUT<<"  :ColumnName";
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      pBasin=pModel->GetSubBasin(p);
      if (pBasin->IsGauged()){
        if (pBasin->GetName()==""){pConstituents[c]->POLLUT<<" ID="<<pBasin->GetID()  ;}
        else                      {pConstituents[c]->POLLUT<<" "   <<pBasin->GetName();}
      }
    }pConstituents[c]->POLLUT<<endl;

    pConstituents[c]->POLLUT<<"  :ColumnUnits"; 
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      if (pModel->GetSubBasin(p)->IsGauged()){pConstituents[c]->POLLUT<<" "<<mgL;}
    }
    pConstituents[c]->POLLUT<<endl;

    pConstituents[c]->POLLUT<<"  :ColumnType"; 
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      if (pModel->GetSubBasin(p)->IsGauged()){pConstituents[c]->POLLUT<<" float";}
    }
    pConstituents[c]->POLLUT<<endl;

    pConstituents[c]->POLLUT << "  :ColumnFormat";
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      if (pModel->GetSubBasin(p)->IsGauged()){pConstituents[c]->POLLUT<<" 0";}
    }
    pConstituents[c]->POLLUT<< endl;

    pConstituents[c]->POLLUT<<":EndColumnMetaData"<<endl;
    pConstituents[c]->POLLUT<<":EndHeader"<<endl;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output to file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams; called from CModel::WriteMinorOutput()
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
//
void CTransportModel::WriteMinorOutput(const optStruct &Options, const time_struct &tt) const
{
  double currentMass,CumInflux,CumOutflux,initMass; //[kg]
  double M; //[mg/m2]
  double V; //[mm] 
  double concentration; //[mg/L]
  int    iCumPrecip;

  string thisdate=tt.date_string;
  string thishour=DecDaysToHours(tt.julian_day);

  double area=pModel->GetWatershedArea(); //[km2]

  iCumPrecip=pModel->GetStateVarIndex(ATMOS_PRECIP);

  double convert;//
  convert=1.0/MG_PER_KG; //[mg->kg]
  //convert=1.0/(area*M2_PER_KM2); //[mg->mg/m2]//TMP DEBUG OUTPUT OVERRIDE

  if ((Options.suppressICs) && (tt.model_time==0.0)) { return; }

  for (int c=0;c<nConstituents;c++)
  {
    // Concentrations.csv
    //----------------------------------------------------------------
    double influx      =0;//GetAverageInflux(c)*(area*M2_PER_KM2);//[mg/d] // \todo [funct]: create GetAverageInflux() routine
    double channel_stor=0;//GetTotalChannelConstituentStorage(c);//[mg]// \todo [funct]: create GetTotalChannelConstituentStorage() routine
    double rivulet_stor=pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT_SW,c))*(area*M2_PER_KM2);//GetTotalRivuletConstituentStorage();//[mg]// \todo [funct]: create GetTotalRivuletConstituentStorage() routine

    pConstituents[c]->OUTPUT<<tt.model_time <<","<<thisdate<<","<<thishour;

    if (tt.model_time!=0.0){pConstituents[c]->OUTPUT<<","<<influx*convert;}
    else                   {pConstituents[c]->OUTPUT<<",---";}
    pConstituents[c]->OUTPUT<<","<<channel_stor*convert<<","<<rivulet_stor*convert;

    currentMass=0.0; 
    double atmos_prec=0;
    for (int j=0;j<nWaterCompartments;j++)
    { 
      //Get constituent concentration
      int m=j+c*nWaterCompartments;
      M=pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT,m)); //mg/m2
      V=pModel->GetAvgStateVar(iWaterStorage[j]); //mm
      if (fabs(V)<=1e-6){concentration=0.0;}
      else              {concentration=(M/V)*(MM_PER_METER/LITER_PER_M3);}//[mg/mm/m2]->[mg/L]
         
      //concentration=M;//TMP DEBUG OUTPUT OVERRIDE [mg/m2]

      if (iWaterStorage[j]!=iCumPrecip)
      {
        pConstituents[c]->OUTPUT<<","<<concentration;   //print column entry

        currentMass+=M*(area*M2_PER_KM2); //mg          //increment total mass in system
      }
      else{
        atmos_prec+=M*(area*M2_PER_KM2); //mg
      }
    }
    currentMass+=channel_stor+rivulet_stor;

    CumInflux =-pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT_SRC,c))*(area*M2_PER_KM2);//mg
    CumInflux +=-atmos_prec;//mg
    CumOutflux=pConstituents[c]->cumul_output;
    initMass  =pConstituents[c]->initial_mass;

    pConstituents[c]->OUTPUT<<","<<currentMass*convert; 
    pConstituents[c]->OUTPUT<<","<<CumInflux*convert;
    pConstituents[c]->OUTPUT<<","<<CumOutflux*convert;
    pConstituents[c]->OUTPUT<<","<<((currentMass-initMass)+(CumOutflux-CumInflux))*convert; 
    pConstituents[c]->OUTPUT<<endl;

    // Pollutographs.csv
    //----------------------------------------------------------------
    pConstituents[c]->POLLUT<<tt.model_time<<","<<thisdate<<","<<thishour;  
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      if (pModel->GetSubBasin(p)->IsGauged())
      {
        pConstituents[c]->POLLUT<<","<<GetOutflowConcentration(p,c);
      }
    }
    pConstituents[c]->POLLUT<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output to ensim formatted file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams; called from CModel::WriteMinorOutput()
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
/// \todo [reorg] merge with WriteMinorOutput - too much complex code repetition here when only difference is (1) delimeter and (2) timestep info included in the .csv file
//
void CTransportModel::WriteEnsimMinorOutput(const optStruct &Options, const time_struct &tt) const
{
  double currentMass,CumInflux,CumOutflux,initMass; //[kg]
  double M; //[mg/m2]
  double V; //[mm] 
  double concentration; //[mg/L]
  int    iCumPrecip;

  string thisdate=tt.date_string;
  string thishour=DecDaysToHours(tt.julian_day);

  double area=pModel->GetWatershedArea(); //[km2]

  iCumPrecip=pModel->GetStateVarIndex(ATMOS_PRECIP);

  double convert;//
  convert=1.0/MG_PER_KG; //[mg->kg]
  //convert=1.0/(area*M2_PER_KM2); //[mg->mg/m2]//TMP DEBUG OUTPUT OVERRIDE

  for (int c=0;c<nConstituents;c++)
  {
    // Concentrations.tb0
    //----------------------------------------------------------------
    double influx      =0;//GetAverageInflux(c)*(area*M2_PER_KM2);//[mg/d] 
    double channel_stor=0;//GetTotalChannelConstituentStorage();//[mg]
    double rivulet_stor=pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT_SW,c))*(area*M2_PER_KM2);//GetTotalRivuletConstituentStorage();//[mg]
  

    if (tt.model_time!=0){pConstituents[c]->OUTPUT<<" "<<influx*convert;}
    else                 {pConstituents[c]->OUTPUT<<" 0.0";}
    pConstituents[c]->OUTPUT<<" "<<channel_stor*convert<<" "<<rivulet_stor*convert;

    currentMass=0.0; 
    double atmos_prec=0;
    for (int j=0;j<nWaterCompartments;j++)
    { 
      //Get constituent concentration
      int m=j+c*nWaterCompartments;
      M=pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT,m)); //mg/m2
      V=pModel->GetAvgStateVar(iWaterStorage[j]); //mm
      if (fabs(V)<=1e-6){concentration=0.0;}
      else              {concentration=(M/V)*(MM_PER_METER/LITER_PER_M3);}//[mg/mm/m2]->[mg/L]
         
      //concentration=M;//TMP DEBUG OUTPUT OVERRIDE [mg/m2]

      if (iWaterStorage[j]!=iCumPrecip)
      {
        pConstituents[c]->OUTPUT<<" "<<concentration;   //print column entry

        currentMass+=M*(area*M2_PER_KM2); //mg          //increment total mass in system
      }
      else{
        atmos_prec+=M*(area*M2_PER_KM2); //mg
      }
    }
    currentMass+=channel_stor+rivulet_stor;
    

    CumInflux =-pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT_SRC,c))*(area*M2_PER_KM2);//mg
    CumInflux +=-atmos_prec;                       //mg
    CumOutflux=pConstituents[c]->cumul_output;     //mg
    initMass  =pConstituents[c]->initial_mass;     //mg

    pConstituents[c]->OUTPUT<<" "<<currentMass*convert; //kg
    pConstituents[c]->OUTPUT<<" "<<CumInflux*convert;   //kg
    pConstituents[c]->OUTPUT<<" "<<CumOutflux*convert;  //kg
    pConstituents[c]->OUTPUT<<" "<<((currentMass-initMass)+(CumOutflux-CumInflux))*convert; //kg
    pConstituents[c]->OUTPUT<<endl;

    // Pollutographs.tb0
    //----------------------------------------------------------------
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      if (pModel->GetSubBasin(p)->IsGauged())
      {
        pConstituents[c]->POLLUT<<" "<<GetOutflowConcentration(p,c);
      }
    }
    pConstituents[c]->POLLUT<<endl;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Close transport output files
//
void CTransportModel::CloseOutputFiles() const
{
  for (int c=0;c<nConstituents;c++)
  {
    pConstituents[c]->OUTPUT.close();
    pConstituents[c]->POLLUT.close();
  }
}