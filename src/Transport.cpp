/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------
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

  _nAdvConnections=0;
  _iFromWater=NULL;
  _iToWater=NULL;
  _js_indices=NULL;

  _nLatConnections=0;
  _iLatFromHRU=NULL;
  _iLatToHRU=NULL;
  _iLatFromWater=NULL;
  _iLatToWater=NULL;
  _latqss_indices=NULL;

  _nWaterCompartments=0;
  _iWaterStorage=NULL;

  _nConstituents=0;
  _pConstituents=NULL;
  _pConstitParams=NULL;

  pSources=NULL;
  nSources=0;

  _aIndexMapping=NULL;
  _aSourceIndices=NULL;

  _pTransModel=this;

  _aMinHist =NULL;
  _aMlatHist=NULL;
  _aMout    =NULL;

}
//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Transport destructor
//
CTransportModel::~CTransportModel()
{
  delete [] _iFromWater;     _iFromWater=NULL;
  delete [] _iToWater;       _iToWater=NULL;
  delete [] _js_indices;     _js_indices=NULL;
  delete [] _iLatFromWater;  _iLatFromHRU=NULL;
  delete [] _iLatToWater;    _iLatToWater=NULL;
  delete [] _iLatFromHRU;    _iLatFromHRU=NULL;
  delete [] _iLatToHRU;      _iLatToHRU=NULL;
  delete [] _iWaterStorage;  _iWaterStorage=NULL;
  delete [] _pConstituents;  _pConstituents=NULL;
  delete [] _pConstitParams; _pConstitParams=NULL;
  delete [] _aIndexMapping;  _aIndexMapping=NULL;
  for (int i=0;i<nSources;i++     ){delete pSources[i];         } delete [] pSources;
  if (_aSourceIndices!=NULL){
    for (int c=0;c<_nConstituents;c++){delete [] _aSourceIndices[c];} delete [] _aSourceIndices;
  }
  DeleteRoutingVars();
}

//Static declaration
CTransportModel* CTransportModel::_pTransModel=NULL;

//////////////////////////////////////////////////////////////////
/// \brief converts model layer index (i.e., m in CONSTITUENT[m]) to
/// local constituent index c and index of corresponding local constiuent storage index j
//
void CTransportModel::m_to_cj(const int m, int &c, int &j) const
{
  ExitGracefullyIf((m<0) || (m>_nConstituents*_nWaterCompartments),
                   "CTransportModel::LayerIndToConstitData: invalid layer index",BAD_DATA);

  j=m % _nWaterCompartments;
  c=(m-j) / _nWaterCompartments;
}

//////////////////////////////////////////////////////////////////
/// \brief returns layer index m corresponding to constituent c in water storage unit i_stor
//
int CTransportModel::GetLayerIndex(const int c, const int i_stor) const
{
  int j=_aIndexMapping[i_stor];
  if (j==DOESNT_EXIST){return DOESNT_EXIST;}
  if (c==DOESNT_EXIST){return DOESNT_EXIST;}
  /*ExitGracefullyIf(j==DOESNT_EXIST,
    "CTransportModel::GetLayerIndex: constituent storage unit not found. Invalid index passed",RUNTIME_ERR);*/
  //cout<<" layer index = "<<c*_nWaterCompartments+j<<" istor: "<<i_stor<<" j: "<<j<<endl;
  return c*_nWaterCompartments+j;
}

//////////////////////////////////////////////////////////////////
/// \brief returns layer index m corresponding to a specific constituent storage compartment (static version)
/// \input name in format !Nitrogen|SOIL[2]
/// \input comp_m layer index of storage compartment (e.g., 2 in above example)
//
int CTransportModel::GetLayerIndexFromName(const string name,const int comp_m) //static
{
  return _pTransModel->GetLayerIndexFromName2(name,comp_m);
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
  cc    =_pTransModel->GetConstituentIndex(constituent_name);
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
  return _pTransModel->GetConstituentName(m);
}

//////////////////////////////////////////////////////////////////
/// \brief returns full name of constitutent e.g., "Nitrogen in Soil Water[2]"
//
string CTransportModel::GetConstituentName(const int m) const
{
  int c,j;
  m_to_cj(m,c,j);
  sv_type typ=pModel->GetStateVarType(_iWaterStorage[j]);
  int     ind=pModel->GetStateVarLayer(_iWaterStorage[j]);

  return _pConstituents[c]->name+" in "+CStateVariable::GetStateVarLongName(typ,ind);
}
//////////////////////////////////////////////////////////////////
/// \brief returns name of constitutent type e.g., "Nitrogen"
/// \remark static routine
/// \param c index of constituent storage
//
string CTransportModel::GetConstituentTypeName2(const int c)
{
  return _pTransModel->GetConstituent(c)->name;
}
//////////////////////////////////////////////////////////////////
/// \brief returns global sv index of water storage unit corresponding to CONSTITUENT[m]
//
int CTransportModel::GetWaterStorIndexFromLayer(const int m) const
{
  int c,j;
  m_to_cj(m,c,j);
  return _iWaterStorage[j];
}
//////////////////////////////////////////////////////////////////
/// \brief returns full name of constitutent e.g., "!Nitrogen|SOIL[2]"
//
string CTransportModel::GetConstituentShortName(const int m) const
{
  int c,j;
  m_to_cj(m,c,j);
  sv_type typ=pModel->GetStateVarType(_iWaterStorage[j]);
  int     ind=pModel->GetStateVarLayer(_iWaterStorage[j]);

  return "!"+_pConstituents[c]->name+"|"+CStateVariable::SVTypeToString(typ,ind);
}

//////////////////////////////////////////////////////////////////
/// \brief returns name of constitutent type e.g., "Nitrogen"
/// \remark static routine
/// \param m layer index of constituent storage
//
string CTransportModel::GetConstituentTypeName(const int m)
{
  int c,j;
  _pTransModel->m_to_cj(m,c,j);
  return _pTransModel->GetConstituent(c)->name;
}

//////////////////////////////////////////////////////////////////
/// \brief returns number of water compartments in model
//
int    CTransportModel::GetNumWaterCompartments() const{return _nWaterCompartments;}

//////////////////////////////////////////////////////////////////
/// \brief returns number of water transport connections in model
//
int    CTransportModel::GetNumAdvConnections() const {return _nAdvConnections;}

//////////////////////////////////////////////////////////////////
/// \brief returns number of lateral water transport connections in model
//
int    CTransportModel::GetNumLatAdvConnections() const { return _nLatConnections; }

//////////////////////////////////////////////////////////////////
/// \brief returns number of constituents transported in model
//
int    CTransportModel::GetNumConstituents() const{return _nConstituents;}

//////////////////////////////////////////////////////////////////
/// \brief returns constituent c in model
//
const constituent *CTransportModel::GetConstituent(const int c) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((c<0)||(c>=_nConstituents),"CTransportModel::GetConstituent: invalid index",BAD_DATA);
#endif
  return _pConstituents[c];
}
//////////////////////////////////////////////////////////////////
/// \brief returns constituent c in model
//
const transport_params *CTransportModel::GetConstituentParams(const int c) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((c<0) || (c >= _nConstituents),"CTransportModel::GetConstituent: invalid index", BAD_DATA);
#endif
  return _pConstitParams[c];
}

//////////////////////////////////////////////////////////////////
/// \brief returns index c of constituents transported in model
//
int    CTransportModel::GetConstituentIndex(const string name) const
{
  for (int c=0; c<_nConstituents;c++){
    if (StringToUppercase(_pConstituents[c]->name)==StringToUppercase(name)){return c;}
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
  int j=_aIndexMapping[_iFromWater[q]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*_nWaterCompartments+j);
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of "to" constituent mass compartment
/// \param c [in] constituent index
/// \param q [in] local index of connection
//
int    CTransportModel::GetToIndex  (const int c,const int q) const
{
  int j=_aIndexMapping[_iToWater[q]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*_nWaterCompartments+j);
}
//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of lateral "from" constituent mass compartment
/// \param c [in] constituent index
/// \param qq [in] local index of lateral connection
//
int    CTransportModel::GetLatFromIndex(const int c,const int qq) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf(_iLatFromWater==NULL,"CTransportModel::GetLatFromIndex: NULL array",RUNTIME_ERR);
#endif 
  int j=_aIndexMapping[_iLatFromWater[qq]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*_nWaterCompartments+j);
}
//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of lateral "to" constituent mass compartment
/// \param c [in] constituent index
/// \param qq [in] local index of lateral connection
//
int    CTransportModel::GetLatToIndex(const int c,const int qq) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf(_iLatToWater==NULL,"CTransportModel::GetLatFromIndex: NULL array",RUNTIME_ERR);
#endif 
  int j=_aIndexMapping[_iLatToWater[qq]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*_nWaterCompartments+j);

}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index of  constituent mass compartment
/// \param c [in] constituent index
/// \param ii [in] local index of water storage unit
//
int    CTransportModel::GetStorIndex  (const int c,const int ii) const
{
  int j=_aIndexMapping[_iWaterStorage[ii]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*_nWaterCompartments+j);
}

//////////////////////////////////////////////////////////////////
/// \brief returns state variable index i of "from" water compartment
/// \param q [in] local index of connection
//
int    CTransportModel::GetFromWaterIndex(const int q) const
{
  return _iFromWater[q];
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of "to" water compartment
/// \param q [in] index of connection
//
int    CTransportModel::GetToWaterIndex  (const int q) const
{
  return _iToWater[q];
}
//////////////////////////////////////////////////////////////////
/// \brief returns state variable index i of lateral "from" water compartment
/// \param q [in] local index of connection
//
int    CTransportModel::GetLatFromWaterIndex(const int qq) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf(_iLatFromWater==NULL,"CTransportModel::GetLatFromIndex: NULL array",RUNTIME_ERR);
#endif 
  return _iLatFromWater[qq];
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of "to" water compartment
/// \param q [in] index of connection
//
int    CTransportModel::GetLatToWaterIndex(const int qq) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf(_iLatToWater==NULL,"CTransportModel::GetLatFromIndex: NULL array",RUNTIME_ERR);
#endif 
  return _iLatToWater[qq];
}
//////////////////////////////////////////////////////////////////
/// \brief returns state variable index i of lateral "from" water compartment
/// \param q [in] local index of connection
//
int    CTransportModel::GetLatFromHRU(const int qq) const
{
  return _iLatFromHRU[qq];
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of "to" water compartment
/// \param q [in] index of connection
//
int    CTransportModel::GetLatToHRU(const int qq) const
{
  return _iLatToHRU[qq];
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of water compartment ii
/// \param ii [in] index of water storage (from 0 to _nWaterCompartments-1)
//
int    CTransportModel::GetStorWaterIndex  (const int ii) const
{
  return _iWaterStorage[ii];
}

//////////////////////////////////////////////////////////////////
/// \brief returns master process index js of advection connection
/// \param q [in] local index of connection
//
int CTransportModel::GetJsIndex(const int q) const
{
  return _js_indices[q];
}
//////////////////////////////////////////////////////////////////
/// \brief returns master process index js of advection connection
/// \param q [in] local index of connection
//
int CTransportModel::GetLatqsIndex(const int qq) const
{
  return _latqss_indices[qq];
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if constituent c is passive 
//
bool  CTransportModel::ConstituentIsPassive(const int c) const {
  return _pConstituents[c]->is_passive;
}
//////////////////////////////////////////////////////////////////
/// \brief returns effective retardation factor for constituent c
/// being transported from storage compartment _iFromWater to storage compartment _iToWater
/// \param c [in] constituent index
/// \param _iFromWater [in] index of "from" water storage state variable
/// \param _iToWater [in] index of "to" water storage state variable
//
double CTransportModel::GetRetardationFactor(const int c,const CHydroUnit *pHRU, const  int _iFromWater,const int _iToWater) const
{
  sv_type fromType,toType;
  fromType=pModel->GetStateVarType(_iFromWater);
  toType  =pModel->GetStateVarType(_iToWater);

  if (ConstituentIsPassive(c)){return ALMOST_INF;}
  if (fromType==SOIL)
  {
    int    m=pModel->GetStateVarLayer(_iFromWater);
#ifdef _STRICTCHECK_
    ExitGracefullyIf(m==DOESNT_EXIST,"GetRetardationFactor:invalid _iFromWater",RUNTIME_ERR);
#endif
    if (toType!=ATMOSPHERE){return pHRU->GetSoilProps(m)->retardation[c];}
    else                   {return pHRU->GetVegetationProps()->uptake_moderator[c];}
  }
  return 1.0;
}
//////////////////////////////////////////////////////////////////
/// \brief returns decay coefficient for constituent c
/// \param c [in] constituent index
/// \param iStorWater [in] index of water storage state variable
//
double CTransportModel::GetDecayCoefficient(const int c, const CHydroUnit *pHRU,const int iStorWater) const
{
  sv_type storType;
  double decay_coeff = GetConstituentParams(c)->decay_coeff;


  storType=pModel->GetStateVarType(iStorWater);

  //add special decay_coefficients from other processes
  if (storType == SOIL)
  {
    int m = pModel->GetStateVarLayer(iStorWater);
#ifdef _STRICTCHECK_
    ExitGracefullyIf(m==DOESNT_EXIST,"GetDecayCoefficient:invalid _iFromWater",RUNTIME_ERR);
    ExitGracefullyIf((c<0) || (c>MAX_CONSTITUENTS),"GetDecayCoefficient:invalid constit",RUNTIME_ERR);
#endif
    decay_coeff += pHRU->GetSoilProps(m)->mineraliz_rate[c];
    decay_coeff += pHRU->GetSoilProps(m)->loss_rate     [c];//e.g., denitrification
  }
  return decay_coeff;
}

//////////////////////////////////////////////////////////////////
/// \brief returns transformation coefficient for constituent c
/// \param c [in] constituent index
/// \param iStorWater [in] index of water storage state variable
//
double CTransportModel::GetTransformCoefficient(const int c, const int c2, const CHydroUnit *pHRU,const int iStorWater) const
{
  sv_type storType=pModel->GetStateVarType(iStorWater);

  if (storType == SOIL)
  {
    int m = pModel->GetStateVarLayer(iStorWater);
#ifdef _STRICTCHECK_
    ExitGracefullyIf(m==DOESNT_EXIST,"GetTransformCoefficient:invalid iStorWater",RUNTIME_ERR);
    ExitGracefullyIf((c<0) || (c>MAX_CONSTITUENTS),"GetTransformCoefficient:invalid constit",RUNTIME_ERR);
    ExitGracefullyIf((c2<0) || (c2>MAX_CONSTITUENTS),"GetTransformCoefficient:invalid constit2",RUNTIME_ERR);
#endif
    return pHRU->GetSoilProps(m)->transf_coeff[c][c2];
  }
  return 0.0;
}
//////////////////////////////////////////////////////////////////
/// \brief returns stoichiometric transformation coefficient for constituent c
/// \param c [in] constituent index
/// \param iStorWater [in] index of water storage state variable
//
double CTransportModel::GetStoichioCoefficient(const int c, const int c2, const CHydroUnit *pHRU,const int iStorWater) const
{
  sv_type storType=pModel->GetStateVarType(iStorWater);

  if (storType == SOIL)
  {
    int m = pModel->GetStateVarLayer(iStorWater);
#ifdef _STRICTCHECK_
    ExitGracefullyIf(m==DOESNT_EXIST,"GetStoichioCoefficient:invalid iStorWater",RUNTIME_ERR);
    ExitGracefullyIf((c<0) || (c>MAX_CONSTITUENTS),"GetStoichioCoefficient:invalid constit",RUNTIME_ERR);
    ExitGracefullyIf((c2<0) || (c2>MAX_CONSTITUENTS),"GetStoichioCoefficient:invalid constit2",RUNTIME_ERR);
#endif
    return pHRU->GetSoilProps(m)->stoichio_coeff[c][c2];
  }
  return 0.0;
}
//////////////////////////////////////////////////////////////////
/// \brief adds new transportable constituent to model
/// \note adds corresponding state variables to model
/// \param name [in] name of constituent
/// \param is_tracer [in] boolean - true if this is a tracer
//
void   CTransportModel::AddConstituent(string name, bool is_tracer,bool is_passive)
{
  constituent *pConstit=new constituent;

  //Initialize constituent members
  pConstit->name=name;
  pConstit->is_tracer=is_tracer;
  pConstit->can_evaporate=false; //default behaviour - evaporation impossible for most contaminants
  if (is_tracer){
    pConstit->can_evaporate=true;
  }
  pConstit->is_passive=false;
  pConstit->initial_mass=0;
  pConstit->cumul_input =0;
  pConstit->cumul_output=0;

  //add to master list of constituents
  if (!DynArrayAppend((void**&)(_pConstituents),(void*)(pConstit),_nConstituents)){
    ExitGracefully(" CTransportModel::AddConstituent: adding NULL constituent",BAD_DATA);}

  int c=(_nConstituents-1);//constit. index of current constituent

  //Add corresponding constituent storage state variables to model
  sv_type *aSV =new sv_type [_nWaterCompartments];
  int     *aLev=new int     [_nWaterCompartments];
  for (int j=0;j<_nWaterCompartments; j++)
  {
    int m=j+c*_nWaterCompartments;
    aSV [j]=CONSTITUENT;
    aLev[j]=m;
  }
  pModel->AddStateVariables(aSV,aLev,_nWaterCompartments);

  //Add corresponding constituent source & sink state variables to model
  aSV [0]=CONSTITUENT_SRC;
  aLev[0]=c;
  pModel->AddStateVariables(aSV,aLev,1);

  aSV [0]=CONSTITUENT_SINK;
  aLev[0]=c;
  pModel->AddStateVariables(aSV,aLev,1);

  aSV [0]=CONSTITUENT_SW;
  aLev[0]=c;
  pModel->AddStateVariables(aSV,aLev,1);

  //Parameter initialization
  transport_params *pP = new transport_params;
  InitializeConstitParams(pP);
  int junk = _nConstituents - 1;

  //add to master list of constituent parameters
  if (!DynArrayAppend((void**&)(_pConstitParams), (void*)(pP), junk)){
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
  //----------------------------------------------------------------------------
  _nAdvConnections=0;
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
        _nAdvConnections++;
      }
    }
  }

  //determine _iFromWater, _iToWater indices and process index for each and every connection
  //----------------------------------------------------------------------------
  js=0;
  int qq=0;
  _iFromWater=new int [_nAdvConnections];
  _iToWater  =new int [_nAdvConnections];
  _js_indices=new int [_nAdvConnections];
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
        _iFromWater[qq]=iF;
        _iToWater  [qq]=iT;
        _js_indices[qq]=js;
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
  //----------------------------------------------------------------------------
  _nWaterCompartments=0;
  for (int i=0;i<pModel->GetNumStateVars(); i++)
  {
    if (CStateVariable::IsWaterStorage(pModel->GetStateVarType(i)))
    {
      _nWaterCompartments++;
    }
  }
  // identify all state variables which are water storage compartments
  _iWaterStorage=new int [_nWaterCompartments];
  _aIndexMapping=new int [pModel->GetNumStateVars()];//current # of state variables does not include transport constituents
  int j=0; //j sifts through storage compartments j=0.._nWaterCompartments
  for (int i=0;i<pModel->GetNumStateVars(); i++)
  {
    _aIndexMapping[i]=DOESNT_EXIST;
    if (CStateVariable::IsWaterStorage(pModel->GetStateVarType(i)))
    {
      _iWaterStorage[j]=i;
      _aIndexMapping[i]=j;
      j++;
    }
  }

   /// \todo [QA/QC]: check for two constituents with same name?
  if (_nConstituents==0){return;}/// all of the above work necessary even with no transport?

  //Synopsis
  //----------------------------------------------------------------------------
  if (!Options.silent){
    cout<<"===TRANSPORT MODEL SUMMARY=============================="<<endl;
    cout<<"   number of compartments: "<<_nWaterCompartments<<endl;
    cout<<"   number of constituents: "<<_nConstituents<<endl;
    cout<<"   number of connections:  "<<_nAdvConnections<<endl;
    cout<<"   number of lat. connect.:"<<_nLatConnections<<endl;
    cout<<"   number of source terms: "<<nSources<<endl;
    // TMP DEBUG below===================================================
    if (false){
      for (int i=0;i<_nAdvConnections;i++)
      {
        cout<<i<<": moves ";
        cout<< " from "<<CStateVariable::SVTypeToString(pModel->GetStateVarType(_iFromWater[i]),pModel->GetStateVarLayer(_iFromWater[i]));
        cout<< " to "  <<CStateVariable::SVTypeToString(pModel->GetStateVarType(_iToWater  [i]),pModel->GetStateVarLayer(_iToWater  [i]));
        cout<<endl;
      }
      for (int i=0;i<pModel->GetNumStateVars(); i++)
      {
        if (_aIndexMapping[i]!=DOESNT_EXIST){
          cout<<"constituent "<< i;
          cout<<" corresponds to "<<CStateVariable::SVTypeToString(pModel->GetStateVarType(_aIndexMapping[i]),pModel->GetStateVarLayer(_aIndexMapping[i]));
          cout<<endl;
        }
        else
        {
          cout<<"SV "<<i<<" has no corresponding mapping"<<endl;
        }
      }
      int c=0;
      for (int j=0;j<_nWaterCompartments;j++)
      {
        int m=c*_nWaterCompartments+j;
        cout<<"water compartment "<< CStateVariable::SVTypeToString(pModel->GetStateVarType(_iWaterStorage[j]),pModel->GetStateVarLayer(_iWaterStorage[j]));
        cout<<" corresponds to "<<CStateVariable::SVTypeToString(CONSTITUENT,m);
        cout<<endl;
      }
    }
    // ====================================================================
    cout<<"========================================================"<<endl;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief generates member arrays of information about lateral connections between HRUs
/// \note must be called AFTER initialization of lateral water flow processes
//
void   CTransportModel::CalculateLateralConnections()
{
  CHydroProcessABC           *pProc;  
  CLateralExchangeProcessABC *pLatProc;

  if (_nConstituents==0){return;}

  //determine number of lateral connections
  _nLatConnections=0;
  for(int j=0; j<pModel->GetNumProcesses(); j++)
  {
    pProc=pModel->GetProcess(j);
    pLatProc=(CLateralExchangeProcessABC*)pProc;//re-cast
    if(pLatProc->GetNumLatConnections()>0){
      //cout<<"CTransport: CalculateLateralConnections: process "<< j<<" #conns: " <<pLatProc->GetNumLatConnections()<<endl;
    }
    for(int q=0; q<pLatProc->GetNumLatConnections(); q++)
    {
      int iF=pLatProc->GetLateralFromIndices()[q];
      int iT=pLatProc->GetLateralToIndices()[q];
      sv_type typF=pModel->GetStateVarType(iF);
      sv_type typT=pModel->GetStateVarType(iT);
      if((iF!=iT) && (CStateVariable::IsWaterStorage(typF)) && (CStateVariable::IsWaterStorage(typT)))
      { //This is a process that may advect a constituent
        _nLatConnections++;
      }
    }
  }

  //determine _iLatFromWater, _iLatToWater indices and process index for each and every lateral connection
  int qss=0;
  int qq=0;
  _iLatFromWater =new int[_nLatConnections];
  _iLatToWater   =new int[_nLatConnections];
  _iLatFromHRU   =new int[_nLatConnections];
  _iLatToHRU     =new int[_nLatConnections];
  _latqss_indices=new int[_nLatConnections];
  for(int j=0; j<pModel->GetNumProcesses(); j++)
  {
    pProc=pModel->GetProcess(j);
    pLatProc=(CLateralExchangeProcessABC*)pProc;//re-cast

    for(int q=0; q<pProc->GetNumLatConnections(); q++)
    {
      int iF=pLatProc->GetLateralFromIndices()[q];
      int iT=pLatProc->GetLateralToIndices()[q];
      int kF=pLatProc->GetFromHRUIndices()[q];
      int kT=pLatProc->GetToHRUIndices()[q];
      sv_type typF=pModel->GetStateVarType(iF);
      sv_type typT=pModel->GetStateVarType(iT);
      if((iF!=iT) && (CStateVariable::IsWaterStorage(typF)) && (CStateVariable::IsWaterStorage(typT)))
      { //This is a process that may advect a constituent
        _iLatFromWater [qq]=iF;
        _iLatToWater   [qq]=iT;
        _iLatToHRU     [qq]=kT;
        _iLatFromHRU   [qq]=kF;
        _latqss_indices[qq]=qss;
        //cout<<"CalculateLateralConnections: ADDED CONNECTION (iF,iT): ("<<iF<<","<<iT<<") (kF,kT): ("<<kF<<","<<kT<<") qss:"<<qss<<endl;
        qq++;
      }
      ExitGracefullyIf(qq>_nLatConnections,"CTransport::CalculateLateralConnections",RUNTIME_ERR);
      qss++;
    }
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
  for (int c=0; c<_nConstituents;c++){
    if (StringToUppercase(const_name)==StringToUppercase(_pConstituents[c]->name)){
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
  for (int c=0; c<_nConstituents;c++){
    if (StringToUppercase(const_name)==StringToUppercase(_pConstituents[c]->name)){
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
/// \param kk [in] HRU group index (or DOESNT_EXIST if this applies to all HRUs)
/// \param flux [in] specified fixed mass influx rate [mg/m2/d]
//
void   CTransportModel::AddInfluxSource(const string const_name, const int i_stor, const int kk, const double flux)
{
  static constit_source *pLast;
  constit_source *pSource=new constit_source();

  pSource->constit_index=DOESNT_EXIST;
  for (int c=0; c<_nConstituents;c++){
    if (StringToUppercase(const_name)==StringToUppercase(_pConstituents[c]->name)){
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
/// \param kk [in] HRU group index (or DOESNT_EXIST if this applies to all HRUs)
/// \param pTS [in] Time series of Dirichlet source flux rate [mg/m2/d]
//
void   CTransportModel::AddInfluxTimeSeries(const string const_name, const int i_stor, const int kk, const CTimeSeries *pTS)
{
  static constit_source *pLast;
  constit_source *pSource=new constit_source();
  pSource->constit_index=DOESNT_EXIST;
  for (int c=0; c<_nConstituents;c++){
    if (StringToUppercase(const_name)==StringToUppercase(_pConstituents[c]->name)){
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
/// \brief Set transport parameter value for specified constituent
//
void   CTransportModel::SetGlobalParameter(const string const_name,const string param_name,const double &value, bool noisy)
{
  int constit_ind=DOESNT_EXIST;
  for (int c=0; c<_nConstituents;c++){
    if (StringToUppercase(const_name)==StringToUppercase(_pConstituents[c]->name)){
      constit_ind=c;
    }
  }
  if(constit_ind==DOESNT_EXIST){
    WriteWarning("CTransportModel::SetGlobalParameter: Unrecognized constituent name",noisy);
  }

  if(StringToUppercase(param_name)=="DECAY_COEFF"){
    _pConstitParams[constit_ind]->decay_coeff=value;
    ExitGracefullyIf(value<0.0," CTransportModel::SetGlobalParameter: decay coefficient cannot be negative",BAD_DATA_WARN);
  }
  else{
    WriteWarning("CTransportModel::SetGlobalParameter: Unrecognized parameter name",noisy);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Initialization of all transport variables
/// \note determines initial conditions for all constituents, initializes routing variables
/// \note  called after all constituents have been added by CModel::Initialize
//
void CTransportModel::Initialize()
{
  double area=pModel->GetWatershedArea();
  for (int c=0;c<_nConstituents;c++)
  {
    _pConstituents[c]->cumul_input=0;
    _pConstituents[c]->cumul_output=0;
    for (int ii=0;ii<_nWaterCompartments;ii++)
    {
      int m=c*_nWaterCompartments+ii;
      int i=pModel->GetStateVarIndex(CONSTITUENT,m);

      // zero out initial mass in all storage units with zero volume (to avoid absurdly high concentrations)
      double watstor;
      int i_stor=GetStorWaterIndex(ii);
      for (int k = 0; k < pModel->GetNumHRUs(); k++){
        watstor = pModel->GetHydroUnit(k)->GetStateVarValue(i_stor);
        if (watstor<1e-9){pModel->GetHydroUnit (k)->SetStateVarValue(i,0.0);}
      }
      //update initial model mass (initialized from .rvc file)
      _pConstituents[c]->initial_mass+=pModel->GetAvgStateVar(i)*(area*M2_PER_KM2); //mg
    }
  }



  /// \todo [funct]: calculate initial mass from Dirichlet cells


  //populate array of source indices
  // \todo [funct] will have to revise to support different sources in different HRUs (e.g., _aSourceIndices[c][i_stor][k])
  _aSourceIndices=NULL;
  _aSourceIndices=new int *[_nConstituents];
  ExitGracefullyIf(_aSourceIndices==NULL,"CTransport::Initialize",OUT_OF_MEMORY);
  for (int c=0;c<_nConstituents;c++)
  {
    _aSourceIndices[c]=NULL;
    _aSourceIndices[c]=new int [pModel->GetNumStateVars()];
    ExitGracefullyIf(_aSourceIndices[c]==NULL,"CTransport::Initialize",OUT_OF_MEMORY);
    for (int i_stor=0;i_stor<pModel->GetNumStateVars();i_stor ++){
      _aSourceIndices[c][i_stor]=DOESNT_EXIST;
      for (int i=0;i<nSources;i++){
        if ((pSources[i]->i_stor==i_stor) && (pSources[i]->constit_index==c))
        {
          if (_aSourceIndices[c][i_stor] != DOESNT_EXIST){
            WriteWarning("CTransportModel::Intiialize: cannot currently have more than one constitutent source per constituent/storage combination",false);
          }
          _aSourceIndices[c][i_stor]=i; //each
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
  //double influx;
  for (int c=0;c<_nConstituents;c++)
  {
    _pConstituents[c]->cumul_input+=0;// \todo [funct]: increment cumulative input [mg]

    //Most of this is handled using the CONSTITUENT_SRC storage compartment
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Increment cumulative output of mass from watershed exiting via stream channel
/// \note called within solver to track mass balance
/// \param Options [in] Global model options structure
//
void CTransportModel::IncrementCumulOutput(const optStruct &Options)
{
  for (int c=0;c<_nConstituents;c++)
  {
    for (int p=0;p<pModel->GetNumSubBasins();p++)
    {
      if (pModel->GetSubBasin(p)->GetDownstreamID()==DOESNT_EXIST)//outlet does not drain into another subbasin
      {
        _pConstituents[c]->cumul_output+=GetIntegratedMassOutflow(p,c,Options.timestep);
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
bool  CTransportModel::IsDirichlet(const int i_stor, const int c, const int k, const time_struct &tt, double &Cs) const
{
  Cs=0.0;

  int i_source=_aSourceIndices[c][i_stor];
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
double  CTransportModel::GetSpecifiedMassFlux(const int i_stor, const int c, const int k, const time_struct &tt) const
{
  double flux;
  bool retrieve=false;
  int i_source=_aSourceIndices[c][i_stor];
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
    flux=pSources[i_source]->concentration; //'concentration' stores mass flux [mg/m2/d] if this is a flux-source
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

  for (int c=0;c<_nConstituents;c++)
  {

    //units names
    kg="[kg]"; kgd="[kg/d]"; mgL="[mg/l]";
    if (_pConstituents[c]->is_tracer){
      kg="[-]"; kgd="[-]"; mgL="[-]";
    }
    if(Options.write_constitmass){
      kg="[mg/m2]"; kgd="[mg/m2/d]"; //mgL="[mg/m2]";
    }

    //Concentrations file
    //--------------------------------------------------------------------
    filename=_pConstituents[c]->name+"Concentrations.csv";
    if(Options.write_constitmass){ filename=_pConstituents[c]->name+"Mass.csv"; }
    filename=FilenamePrepare(filename,Options);

    _pConstituents[c]->OUTPUT.open(filename.c_str());
    if (_pConstituents[c]->OUTPUT.fail()){
      ExitGracefully(("CTransportModel::WriteOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
    }

    _pConstituents[c]->OUTPUT<<"time[d],date,hour,influx"<<kgd<<",Channel Storage"<<kg<<",Rivulet Storage"<<kg;
    for (int i=0;i<pModel->GetNumStateVars();i++)
    {
      if ((CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))) && (i!=iCumPrecip)){
        _pConstituents[c]->OUTPUT<<","<<CStateVariable::GetStateVarLongName(pModel->GetStateVarType(i),pModel->GetStateVarLayer(i))<<" "<<mgL;
      }
    }
    _pConstituents[c]->OUTPUT<<", Total Mass "<<kg<<", Cum. Loading "<<kg<<", Cum. Mass Lost "<<kg<<", MB Error "<<kg<<endl;

    //Pollutograph file
    //--------------------------------------------------------------------
    filename=_pConstituents[c]->name+"Pollutographs.csv";
    filename=FilenamePrepare(filename,Options);

    _pConstituents[c]->POLLUT.open(filename.c_str());
    if (_pConstituents[c]->POLLUT.fail()){
      ExitGracefully(("CTransportModel::WriteOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    _pConstituents[c]->POLLUT<<"time[d],date,hour";
    const CSubBasin *pBasin;
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      pBasin=pModel->GetSubBasin(p);
      if (pBasin->IsGauged() && (pBasin->IsEnabled())){
        string name;
        if (pBasin->GetName()==""){_pConstituents[c]->POLLUT<<",ID="<<pBasin->GetID()  <<" "<<mgL;}
        else                      {_pConstituents[c]->POLLUT<<","   <<pBasin->GetName()<<" "<<mgL;}
      }
    }
    _pConstituents[c]->POLLUT<<endl;
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
  JulianConvert(0.0,              Options.julian_start_day, Options.julian_start_year, Options.calendar, tt);
  JulianConvert(Options.timestep, Options.julian_start_day, Options.julian_start_year, Options.calendar, tt2);//end of the timestep

  for (int c=0;c<_nConstituents;c++)
  {
    //units names
    kg="kg"; kgd="kg/d"; mgL="mg/l";
    //kg="mg/m2"; kgd="mg/m2/d"; mgL="mg/m2";//TMP DEBUG OUTPUT OVERRIDE
    if (_pConstituents[c]->is_tracer){
      kg="none"; kgd="none"; mgL="none";
    }

    //Concentrations file
    //--------------------------------------------------------------------
    filename=_pConstituents[c]->name+"Concentrations.tb0";
    filename=FilenamePrepare(filename,Options);

    _pConstituents[c]->OUTPUT.open(filename.c_str());
    if (_pConstituents[c]->OUTPUT.fail()){
      ExitGracefully(("CTransportModel::WriteEnsimOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    _pConstituents[c]->OUTPUT<<"#########################################################################"<<endl;
    _pConstituents[c]->OUTPUT<<":FileType tb0 ASCII EnSim 1.0"<<endl;
    _pConstituents[c]->OUTPUT<<"#"<<endl;
    _pConstituents[c]->OUTPUT<<":Application   Raven"<<endl;
    if(!Options.benchmarking){
      _pConstituents[c]->OUTPUT<<":Version       "<<Options.version<<endl;
      _pConstituents[c]->OUTPUT<<":CreationDate  "<<GetCurrentTime()<<endl;
    }
    _pConstituents[c]->OUTPUT<<"#"<<endl;
    _pConstituents[c]->OUTPUT<<"#------------------------------------------------------------------------"<<endl;
    _pConstituents[c]->OUTPUT<<"#"<<endl;
    _pConstituents[c]->OUTPUT<<":RunName       "<<Options.run_name<<endl;
    _pConstituents[c]->OUTPUT<<":Format        Instantaneous" << endl;
    _pConstituents[c]->OUTPUT<<"#"<<endl;

    if (Options.suppressICs){
      _pConstituents[c]->OUTPUT<< ":StartTime " << tt2.date_string << " " << DecDaysToHours(tt2.julian_day) << endl;
    }
    else{
      _pConstituents[c]->OUTPUT<< ":StartTime " << tt.date_string << " " << DecDaysToHours(tt.julian_day) << endl;
    }

    if (Options.timestep!=1.0){_pConstituents[c]->OUTPUT<<":DeltaT " <<DecDaysToHours(Options.timestep)<<endl;}
    else                      {_pConstituents[c]->OUTPUT<<":DeltaT 24:00:00.00"  <<endl;}
    _pConstituents[c]->OUTPUT<<"#"<<endl;

    _pConstituents[c]->OUTPUT<<":ColumnMetaData"<<endl;
    _pConstituents[c]->OUTPUT<<"  :ColumnName influx \"Channel storage\" \"Rivulet storage\"";
    for (i=0;i<pModel->GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))) && (i!=iCumPrecip)){
        _pConstituents[c]->OUTPUT<<" \""<<CStateVariable::GetStateVarLongName(pModel->GetStateVarType(i),pModel->GetStateVarLayer(i))<<"\"";}}
    _pConstituents[c]->OUTPUT<<" \"Total Mass\" \"Cum. Loading\" \"Cum. Mass lost\" \"MB error\""<<endl;

    _pConstituents[c]->OUTPUT<<"  :ColumnUnits "<<kgd<<" "<<kg <<" "<< kg;
    for (i=0;i<pModel->GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))) && (i!=iCumPrecip)){
        _pConstituents[c]->OUTPUT<<" mgL";}}
    _pConstituents[c]->OUTPUT<<" "<<kg<<" "<<kg<<" "<<kg<<" "<<kg<<endl;

    _pConstituents[c]->OUTPUT<<"  :ColumnType float float float";
    for (i=0;i<pModel->GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))) && (i!=iCumPrecip)){
        _pConstituents[c]->OUTPUT<<" float";}}
    _pConstituents[c]->OUTPUT<<" float float float float"<<endl;

    _pConstituents[c]->OUTPUT<<"  :ColumnFormat -1 0 0";
    for (i=0;i<pModel->GetNumStateVars();i++){
      if ((CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))) && (i!=iCumPrecip)){
        _pConstituents[c]->OUTPUT<<" 0";}}
    _pConstituents[c]->OUTPUT<<" 0 0 0 0"<<endl;

    _pConstituents[c]->OUTPUT<<":EndColumnMetaData"<<endl;
    _pConstituents[c]->OUTPUT<<":EndHeader"<<endl;

    //Pollutograph file
    //--------------------------------------------------------------------
    filename=_pConstituents[c]->name+"Pollutographs.tb0";
    filename=FilenamePrepare(filename,Options);

    _pConstituents[c]->POLLUT.open(filename.c_str());
    if (_pConstituents[c]->OUTPUT.fail()){
      ExitGracefully(("CTransportModel::WriteEnsimOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
    }
    _pConstituents[c]->POLLUT<<"#########################################################################"<<endl;
    _pConstituents[c]->POLLUT<<":FileType tb0 ASCII EnSim 1.0"<<endl;
    _pConstituents[c]->POLLUT<<"#"<<endl;
    _pConstituents[c]->POLLUT<<":Application   Raven"<<endl;
    if(!Options.benchmarking){
      _pConstituents[c]->POLLUT<<":Version       "<<Options.version<<endl;
      _pConstituents[c]->POLLUT<<":CreationDate  "<<GetCurrentTime()<<endl;
    }
    _pConstituents[c]->POLLUT<<"#"<<endl;
    _pConstituents[c]->POLLUT<<"#------------------------------------------------------------------------"<<endl;
    _pConstituents[c]->POLLUT<<"#"<<endl;
    _pConstituents[c]->POLLUT<<":RunName       "<<Options.run_name<<endl;
    _pConstituents[c]->POLLUT<<":Format        Instantaneous" << endl;
    _pConstituents[c]->POLLUT<<"#"<<endl;

    if (Options.suppressICs){
      _pConstituents[c]->POLLUT<< ":StartTime " << tt2.date_string << " " << DecDaysToHours(tt2.julian_day) << endl;
    }
    else{
      _pConstituents[c]->POLLUT<< ":StartTime " << tt.date_string << " " << DecDaysToHours(tt.julian_day) << endl;
    }

    if (Options.timestep!=1.0){_pConstituents[c]->POLLUT<<":DeltaT " <<DecDaysToHours(Options.timestep)<<endl;}
    else                      {_pConstituents[c]->POLLUT<<":DeltaT 24:00:00.00"  <<endl;}
    _pConstituents[c]->POLLUT<<"#"<<endl;

    const CSubBasin *pBasin;

    _pConstituents[c]->POLLUT<<":ColumnMetaData"<<endl;
    _pConstituents[c]->POLLUT<<"  :ColumnName";
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      pBasin=pModel->GetSubBasin(p);
      if (pBasin->IsGauged() && (pBasin->IsEnabled())){
        if (pBasin->GetName()==""){_pConstituents[c]->POLLUT<<" ID="<<pBasin->GetID()  ;}
        else                      {_pConstituents[c]->POLLUT<<" "   <<pBasin->GetName();}
      }
    }
    _pConstituents[c]->POLLUT<<endl;

    _pConstituents[c]->POLLUT<<"  :ColumnUnits";
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      pBasin=pModel->GetSubBasin(p);
      if (pBasin->IsGauged() && (pBasin->IsEnabled())){_pConstituents[c]->POLLUT<<" "<<mgL;}
    }
    _pConstituents[c]->POLLUT<<endl;

    _pConstituents[c]->POLLUT<<"  :ColumnType";
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      pBasin=pModel->GetSubBasin(p);
      if(pBasin->IsGauged() && (pBasin->IsEnabled())) { _pConstituents[c]->POLLUT<<" float";}
    }
    _pConstituents[c]->POLLUT<<endl;

    _pConstituents[c]->POLLUT << "  :ColumnFormat";
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      pBasin=pModel->GetSubBasin(p);
      if(pBasin->IsGauged() && (pBasin->IsEnabled())){_pConstituents[c]->POLLUT<<" 0";}
    }
    _pConstituents[c]->POLLUT<< endl;

    _pConstituents[c]->POLLUT<<":EndColumnMetaData"<<endl;
    _pConstituents[c]->POLLUT<<":EndHeader"<<endl;
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
  if(Options.write_constitmass){ convert=1.0/(area*M2_PER_KM2); } //[mg->mg/m2]

  if ((Options.suppressICs) && (tt.model_time==0.0)) { return; }

  for (int c=0;c<_nConstituents;c++)
  {
    // Concentrations.csv
    //----------------------------------------------------------------
    double influx      =0;//GetAverageInflux(c)*(area*M2_PER_KM2);//[mg/d] // \todo [funct]: create GetAverageInflux() routine
    double channel_stor=0;//GetTotalChannelConstituentStorage(c);//[mg]// \todo [funct]: create GetTotalChannelConstituentStorage() routine
    double rivulet_stor=pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT_SW,c))*(area*M2_PER_KM2);//GetTotalRivuletConstituentStorage();//[mg]// \todo [funct]: create GetTotalRivuletConstituentStorage() routine

    _pConstituents[c]->OUTPUT<<tt.model_time <<","<<thisdate<<","<<thishour;

    if (tt.model_time!=0.0){_pConstituents[c]->OUTPUT<<","<<influx*convert;}
    else                   {_pConstituents[c]->OUTPUT<<",---";}
    _pConstituents[c]->OUTPUT<<","<<channel_stor*convert<<","<<rivulet_stor*convert;

    currentMass=0.0;
    double atmos_prec=0;
    for (int j=0;j<_nWaterCompartments;j++)
    {
      //Get constituent concentration
      int m=j+c*_nWaterCompartments;
      M=pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT,m)); //mass- mg/m2
      V=pModel->GetAvgStateVar(_iWaterStorage[j]); //mm
      if (fabs(V)<=1e-6){concentration=0.0;}
      else              {concentration=(M/V)*(MM_PER_METER/LITER_PER_M3);}//[mg/mm/m2]->[mg/L]
      if(Options.write_constitmass){ concentration=M; }//[mg/m2]

      if (_iWaterStorage[j]!=iCumPrecip)
      {
        _pConstituents[c]->OUTPUT<<","<<concentration;   //print column entry

        currentMass+=M*(area*M2_PER_KM2); //mg          //increment total mass in system
      }
      else{
        atmos_prec+=M*(area*M2_PER_KM2); //mg
      }
    }
    currentMass+=channel_stor+rivulet_stor;

    double sink  = pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT_SINK,c))*(area*M2_PER_KM2);//mg
    double source=-pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT_SRC ,c))*(area*M2_PER_KM2);//mg
    CumInflux =source;
    CumInflux +=-atmos_prec;//mg
    CumOutflux=_pConstituents[c]->cumul_output;//outflow from system (mg)
    CumOutflux+=sink;

    initMass  =_pConstituents[c]->initial_mass;

    _pConstituents[c]->OUTPUT<<","<<currentMass*convert;
    _pConstituents[c]->OUTPUT<<","<<CumInflux*convert;
    _pConstituents[c]->OUTPUT<<","<<CumOutflux*convert;
    _pConstituents[c]->OUTPUT<<","<<((currentMass-initMass)+(CumOutflux-CumInflux))*convert;
    _pConstituents[c]->OUTPUT<<endl;

    // Pollutographs.csv
    //----------------------------------------------------------------
    _pConstituents[c]->POLLUT<<tt.model_time<<","<<thisdate<<","<<thishour;
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      CSubBasin *pBasin=pModel->GetSubBasin(p);
      if(pBasin->IsGauged() && (pBasin->IsEnabled()))
      {
        _pConstituents[c]->POLLUT<<","<<GetOutflowConcentration(p,c);
      }
    }
    _pConstituents[c]->POLLUT<<endl;
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

  for (int c=0;c<_nConstituents;c++)
  {
    // Concentrations.tb0
    //----------------------------------------------------------------
    double influx      =0;//GetAverageInflux(c)*(area*M2_PER_KM2);//[mg/d]
    double channel_stor=0;//GetTotalChannelConstituentStorage();//[mg]
    double rivulet_stor=pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT_SW,c))*(area*M2_PER_KM2);//GetTotalRivuletConstituentStorage();//[mg]


    if (tt.model_time!=0){_pConstituents[c]->OUTPUT<<" "<<influx*convert;}
    else                 {_pConstituents[c]->OUTPUT<<" 0.0";}
    _pConstituents[c]->OUTPUT<<" "<<channel_stor*convert<<" "<<rivulet_stor*convert;

    currentMass=0.0;
    double atmos_prec=0;
    for (int j=0;j<_nWaterCompartments;j++)
    {
      //Get constituent concentration
      int m=j+c*_nWaterCompartments;
      M=pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT,m)); //mg/m2
      V=pModel->GetAvgStateVar(_iWaterStorage[j]); //mm
      if (fabs(V)<=1e-6){concentration=0.0;}
      else              {concentration=(M/V)*(MM_PER_METER/LITER_PER_M3);}//[mg/mm/m2]->[mg/L]

      //concentration=M;//TMP DEBUG OUTPUT OVERRIDE [mg/m2]

      if (_iWaterStorage[j]!=iCumPrecip)
      {
        _pConstituents[c]->OUTPUT<<" "<<concentration;   //print column entry

        currentMass+=M*(area*M2_PER_KM2); //mg          //increment total mass in system
      }
      else{
        atmos_prec+=M*(area*M2_PER_KM2); //mg
      }
    }
    currentMass+=channel_stor+rivulet_stor;


    CumInflux =-pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT_SRC,c))*(area*M2_PER_KM2);//mg
    CumInflux +=-atmos_prec;                       //mg
    CumOutflux=_pConstituents[c]->cumul_output;     //mg
    CumOutflux+=pModel->GetAvgStateVar(pModel->GetStateVarIndex(CONSTITUENT_SINK,c))*(area*M2_PER_KM2);//mg
    initMass  =_pConstituents[c]->initial_mass;     //mg

    _pConstituents[c]->OUTPUT<<" "<<currentMass*convert; //kg
    _pConstituents[c]->OUTPUT<<" "<<CumInflux*convert;   //kg
    _pConstituents[c]->OUTPUT<<" "<<CumOutflux*convert;  //kg
    _pConstituents[c]->OUTPUT<<" "<<((currentMass-initMass)+(CumOutflux-CumInflux))*convert; //kg
    _pConstituents[c]->OUTPUT<<endl;

    // Pollutographs.tb0
    //----------------------------------------------------------------
    for (int p=0;p<pModel->GetNumSubBasins();p++){
      CSubBasin *pBasin=pModel->GetSubBasin(p);
      if (pBasin->IsGauged() && pBasin->IsEnabled())
      {
        _pConstituents[c]->POLLUT<<" "<<GetOutflowConcentration(p,c);
      }
    }
    _pConstituents[c]->POLLUT<<endl;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Close transport output files
//
void CTransportModel::CloseOutputFiles() const
{
  for (int c=0;c<_nConstituents;c++)
  {
    _pConstituents[c]->OUTPUT.close();
    _pConstituents[c]->POLLUT.close();
  }
}
