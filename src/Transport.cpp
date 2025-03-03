/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2022 the Raven Development Team
----------------------------------------------------------------
Master Transport/Tracer class
coordinates information about constituent storage
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Model.h"
#include "Transport.h"
#include "HeatConduction.h"
#include "EnergyTransport.h"
#include "IsotopeTransport.h"

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
  _pConstitModels=NULL;
  _pEnthalpyModel=NULL;

  _pGeochemParams=NULL;
  _nGeochemParams=0;

  _aIndexMapping=NULL; _nIndexMapping=0;

  _aProcessNames=NULL;
  _nProcessNames=0;

  // _pTransModel=this;
}
//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Transport destructor
//
CTransportModel::~CTransportModel()
{
  delete[] _iFromWater;     _iFromWater=NULL;
  delete[] _iToWater;       _iToWater=NULL;
  delete[] _js_indices;     _js_indices=NULL;
  delete[] _iLatFromWater;  _iLatFromHRU=NULL;
  delete[] _iLatToWater;    _iLatToWater=NULL;
  delete[] _iLatFromHRU;    _iLatFromHRU=NULL;
  delete[] _iLatToHRU;      _iLatToHRU=NULL;
  delete[] _iWaterStorage;  _iWaterStorage=NULL;
  delete[] _pConstitModels; _pConstitModels=NULL;
  delete[] _aIndexMapping;  _aIndexMapping=NULL;
  if(_pGeochemParams!=NULL) {
    for(int i=0;i<_nGeochemParams;i++) { delete _pGeochemParams[i]; } delete[] _pGeochemParams;
  }
  delete[] _aProcessNames;  _aProcessNames=NULL;
}

//Static declaration
// CTransportModel* CTransportModel::_pTransModel=NULL;

//////////////////////////////////////////////////////////////////
/// \brief converts model layer index (i.e., m in CONSTITUENT[m]) to
/// local constituent index c and index of corresponding local constituent storage index j
//
void CTransportModel::m_to_cj(const int m,int &c,int &j) const
{
  ExitGracefullyIf((m < 0) || (m > (this->_nConstituents * this->_nWaterCompartments)),
    "CTransportModel::LayerIndToConstitData: invalid layer index",BAD_DATA);

  j=m % _nWaterCompartments;
  c=(m-j) / _nWaterCompartments;
}

//////////////////////////////////////////////////////////////////
/// \brief returns layer index m corresponding to constituent c in water storage unit i_stor
//
int CTransportModel::GetLayerIndex(const int c,const int i_stor) const
{

  ExitGracefullyIf(i_stor>=_nIndexMapping,"CTransportModel::GetLayerIndex: invalid storage index",RUNTIME_ERR);
  int ii=_aIndexMapping[i_stor];
  if(ii==DOESNT_EXIST) { return DOESNT_EXIST; }
  if(c ==DOESNT_EXIST) { return DOESNT_EXIST; }
  /*ExitGracefullyIf(j==DOESNT_EXIST,
  "CTransportModel::GetLayerIndex: constituent storage unit not found. Invalid index passed",RUNTIME_ERR);*/
  //cout<<" layer index = "<<c*_nWaterCompartments+j<<" istor: "<<i_stor<<" j: "<<j<<endl;
  return c*_nWaterCompartments+ii;
}

//////////////////////////////////////////////////////////////////
/// \brief returns layer index m corresponding to a specific constituent storage compartment (static version)
/// \input name in format !Nitrogen|SOIL[2]
/// \input comp_m layer index of storage compartment (e.g., 2 in above example)
//
/*
int CTransportModel::GetLayerIndexFromName(const string name,const int comp_m) //static
{
  return _pTransModel->GetLayerIndexFromName2(name,comp_m);
}
*/

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

  if(name.substr(0,1)!="!") { return DOESNT_EXIST; }//bad format, no leading "!"

  tmp=name.substr(1,name.length()-1);             //trim leading '!'
  k=(int)(tmp.find_first_of("|"));                //find char index of "|"
  if(k==DOESNT_EXIST) { return DOESNT_EXIST; }      //bad format, no "|"

  constituent_name=tmp.substr(0,k);
  compartment_name=tmp.substr(k+1,tmp.length()-k);
  //cout<<":GetLayerIndexFromName2:  "<<name<<" -> "<<tmp<<" -> |"<<constituent_name<<"| + |"<<compartment_name<<"|"<<endl;

  int     cc,layer_index,i_stor;
  sv_type typ;
  cc    = this->GetConstituentIndex(constituent_name);
  typ   = pModel->GetStateVarInfo()->StringToSVType(compartment_name, layer_index, true);
  i_stor = pModel->GetStateVarIndex(typ, comp_m);

  return GetLayerIndex(cc,i_stor);
}

//////////////////////////////////////////////////////////////////
/// \brief returns name of constitutent type e.g., "Nitrogen"
/// \remark static routine
/// \param m layer index of constituent storage
//
string CTransportModel::GetConstituentTypeName(const int m)
{
  int c, j;
  this->m_to_cj(m, c, j);
  return (this->GetConstituentModel(c)->GetName());
}
//////////////////////////////////////////////////////////////////
/// \brief returns name of constitutent type e.g., "Nitrogen"
/// \remark static routine
/// \param c index of constituent storage
//
string CTransportModel::GetConstituentTypeName2(const int c) const
{
  return this->GetConstituentModel(c)->GetName();
}
//////////////////////////////////////////////////////////////////
/// \brief returns full name of constitutent e.g.,  "Nitrogen in Soil Water[2]"
/// static version allowing calls without passing transport model
//
string CTransportModel::GetConstituentLongName(const int m) const
{
  return this->GetConstituentLongName_loc(m);
}
//////////////////////////////////////////////////////////////////
/// \brief returns full name of constitutent e.g., "Nitrogen in Soil Water[2]"
//
string CTransportModel::GetConstituentLongName_loc(const int m) const
{
  int c,j;
  m_to_cj(m,c,j);
  sv_type typ = pModel->GetStateVarType (_iWaterStorage[j]);
  int     ind = pModel->GetStateVarLayer(_iWaterStorage[j]);

  if(_pConstitModels[c]->GetType()==ENTHALPY) {
    return _pConstitModels[c]->GetName()+" of "+ CStateVariable::GetStateVarLongName(typ, ind, this);
  }
  else {
    return _pConstitModels[c]->GetName() + " in " + CStateVariable::GetStateVarLongName(typ, ind, this);
  }
}
//////////////////////////////////////////////////////////////////g_z
/// \brief returns full name of constitutent e.g., "!Nitrogen|SOIL[2]"
//
string CTransportModel::GetConstituentShortName(const int m) const
{
  int c,j;
  m_to_cj(m,c,j);
  sv_type typ = pModel->GetStateVarType(_iWaterStorage[j]);
  int     ind = pModel->GetStateVarLayer(_iWaterStorage[j]);

  return "!" + _pConstitModels[c]->GetName() + "|" + pModel->GetStateVarInfo()->SVTypeToString(typ, ind);
}

//////////////////////////////////////////////////////////////////
/// \brief returns global sv index of water storage unit corresponding to CONSTITUENT[m], or DOESNT_EXIST if layer index is invalid
//
int CTransportModel::GetWaterStorIndexFromLayer(const int m) const
{
  int c,j;
  if (m>_nConstituents*_nWaterCompartments){return DOESNT_EXIST;}
  m_to_cj(m,c,j);
  return _iWaterStorage[j];
}

//////////////////////////////////////////////////////////////////
/// \brief returns number of water compartments in model
//
int    CTransportModel::GetNumWaterCompartments() const { return _nWaterCompartments; }

//////////////////////////////////////////////////////////////////
/// \brief returns number of water transport connections in model
//
int    CTransportModel::GetNumAdvConnections() const    { return _nAdvConnections; }

//////////////////////////////////////////////////////////////////
/// \brief returns number of lateral water transport connections in model
//
int    CTransportModel::GetNumLatAdvConnections() const { return _nLatConnections; }

//////////////////////////////////////////////////////////////////
/// \brief returns number of constituents transported in model
//
int    CTransportModel::GetNumConstituents() const      { return _nConstituents; }


//////////////////////////////////////////////////////////////////
/// \brief returns constituent model c in model
//
CConstituentModel *CTransportModel::GetConstituentModel(const int c) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((c<0) || (c >= _nConstituents),"CTransportModel::GetConstituentModel: invalid index",BAD_DATA);
#endif
  return _pConstitModels[c];
}

//////////////////////////////////////////////////////////////////
/// \brief returns constituent model c in model (const version)
//
CConstituentModel *CTransportModel::GetConstituentModel2(const int c) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((c<0) || (c >= _nConstituents),"CTransportModel::GetConstituentModel2: invalid index",BAD_DATA);
#endif
  return _pConstitModels[c];
}

//////////////////////////////////////////////////////////////////
/// \brief returns index c of constituents transported in model
//
int    CTransportModel::GetConstituentIndex(const string name) const
{
  for(int c=0; c<_nConstituents;c++) {
    if(StringToUppercase(_pConstitModels[c]->GetName())==StringToUppercase(name)) { return c; }
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
  int ii=_aIndexMapping[_iFromWater[q]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*_nWaterCompartments+ii);
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index i of "to" constituent mass compartment
/// \param c [in] constituent index
/// \param q [in] local index of connection
//
int    CTransportModel::GetToIndex(const int c,const int q) const
{
  int ii=_aIndexMapping[_iToWater[q]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*_nWaterCompartments+ii);
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
  int ii=_aIndexMapping[_iLatFromWater[qq]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*_nWaterCompartments+ii);
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
  int ii=_aIndexMapping[_iLatToWater[qq]];
  return pModel->GetStateVarIndex(CONSTITUENT,c*_nWaterCompartments+ii);

}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index of  constituent mass compartment
/// \param c [in] constituent index
/// \param ii [in] local index of water storage unit
//
int    CTransportModel::GetStorIndex(const int c,const int ii) const
{
  int j=_aIndexMapping[_iWaterStorage[ii]]; //JRC: shouldnt these be inverted so this is just j==ii?
  return pModel->GetStateVarIndex(CONSTITUENT,c*_nWaterCompartments+j);
}

//////////////////////////////////////////////////////////////////
/// \brief returns global state variable index of  constituent mass compartment
/// \param i [in] global index of water storage
/// \returns ii, local index of water storage unit
//
int    CTransportModel::GetWaterStorIndexFromSVIndex(const int i) const
{
  //return _aIndexMapping[i ]; //JRC - check - this should be the same, but faster??

  for(int ii=0;ii<_nWaterCompartments;ii++) {
    if(_iWaterStorage[ii]==i) { return ii; }
  }
  return DOESNT_EXIST;
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
int    CTransportModel::GetToWaterIndex(const int q) const
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
int    CTransportModel::GetStorWaterIndex(const int ii) const
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
/// \brief returns pointer to Enthalpy model, or NULL if none exists
//
const CEnthalpyModel *CTransportModel::GetEnthalpyModel() const
{
  return _pEnthalpyModel;
}
//////////////////////////////////////////////////////////////////
/// \brief returns concentration (or temperature, or isotopic composition) in HRU k corresponding to transport state variable with i=sv_index
/// \param k - global HRU index of interest
/// \param sv_index - state variable index (i)
//
double CTransportModel::GetConcentration(const int k,const int sv_index) const
{
  int         m=pModel->GetStateVarLayer(sv_index);
  int    i_stor=GetWaterStorIndexFromLayer(m);
  double    vol=pModel->GetHydroUnit(k)->GetStateVarValue(i_stor);
  double   mass=pModel->GetHydroUnit(k)->GetStateVarValue(sv_index);
  int       c,j;
  m_to_cj(m,c,j);

  return _pConstitModels[c]->CalculateReportingConcentration(mass,vol);
}

//////////////////////////////////////////////////////////////////
/// \brief adds new transportable constituent to model
/// \note adds corresponding state variables to model
/// \param name [in] name of constituent
/// \param type [in] constit_type - type of constituent (mass/energy/tracer)
//
double CTransportModel::GetAdvectionCorrection(const int c,const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& C) const
{
  return _pConstitModels[c]->GetAdvectionCorrection(pHRU,iFromWater,iToWater,C);
}

//////////////////////////////////////////////////////////////////
/// \brief Return index of process, or DOESNT_EXIST if name not found in process list
//
int   CTransportModel::GetProcessIndex(const string name) const
{
  for(int i=0;i<_nProcessNames;i++) {
    if (StringToUppercase(_aProcessNames[i])==StringToUppercase(name)){return i;}
  }
  return DOESNT_EXIST;
}
//////////////////////////////////////////////////////////////////
/// \brief returns geochemical parameter value
/// \notes needs optimization! - perhaps hand in iguess?
/// note - if procind==DOESNT_EXIST, first value of parameter regardless of process name is used (means that (e.g.) global decay rate should be specified before local decay rates in .rvp file)
/// This could be a problem when we have dozens of geochemistry parameters
/// other solution - make 2D array of parameters _pGeochemParams[gtyp][i] with size _nGeochemParams[gtyp]
//
double CTransportModel::GetGeochemParam(const gparam_type gtyp, const int c,const int ii,const int procind,const CHydroUnit* pHRU) const
{
  //cout<<" GEOCHEM PAR 1"<<_nGeochemParams<<endl;
  for(int i=0;i<_nGeochemParams;i++) {
    if (_pGeochemParams[i]->type!=gtyp){ continue;}
    if ((_pGeochemParams[i]->ii_water_stor!=DOESNT_EXIST) && (_pGeochemParams[i]->ii_water_stor!=ii)) { continue;}
    if ((procind!=DOESNT_EXIST) && (_pGeochemParams[i]->process_ind!=procind)) { continue; }
    if (_pGeochemParams[i]->constit_ind!=c){ continue;}
    //TODO - handle soil type
    return _pGeochemParams[i]->value;
  }
  return NOT_SPECIFIED;
}
double CTransportModel::GetGeochemParam(const gparam_type gtyp,const int c,const int c2,const int ii,const int procind,const CHydroUnit* pHRU) const
{
  //cout<<" GEOCHEM PAR "<<_nGeochemParams<<endl;
  for(int i=0;i<_nGeochemParams;i++) {
    //geochem_param *p=_pGeochemParams[i];
    //cout<<" GEOCHEM PAR "<<p->type<<" "<<p->constit_ind<<" "<<p->constit_ind2<<" "<<p->ii_water_stor<<" "<<p->value<<endl;
    if(_pGeochemParams[i]->type!=gtyp) { continue; }
    if((_pGeochemParams[i]->ii_water_stor!=DOESNT_EXIST) && (_pGeochemParams[i]->ii_water_stor!=ii)) { continue; }
    if((procind!=DOESNT_EXIST) && (_pGeochemParams[i]->process_ind!=procind)) { continue; }
    if(_pGeochemParams[i]->constit_ind!=c) { continue; }
    if(_pGeochemParams[i]->constit_ind2!=c2) { continue; }
    //TODO - handle soil type
    return _pGeochemParams[i]->value;
  }
  return NOT_SPECIFIED;
}
//////////////////////////////////////////////////////////////////
/// \brief adds new transportable constituent to model
/// \note adds corresponding state variables to model
/// \param name [in] name of constituent
/// \param type [in] constit_type - type of constituent (mass/energy/tracer)
//
void   CTransportModel::AddConstituent(string name,constit_type typ,bool is_passive)
{
  CConstituentModel *pConstitModel;
  if (typ==ENTHALPY){
    _pEnthalpyModel = new CEnthalpyModel(pModel, this, name, _nConstituents);
    pConstitModel=_pEnthalpyModel;
  }
  else if (typ == ISOTOPE) {
    pConstitModel = new CIsotopeModel(pModel, this, name, _nConstituents);
  }
  else {
    pConstitModel = new CConstituentModel(pModel, this, name, typ, is_passive, _nConstituents);
  }
  //add to master list of constituents
  if(!DynArrayAppend((void**&)(_pConstitModels),(void*)(pConstitModel),_nConstituents)) {
    ExitGracefully(" CTransportModel::AddConstituent: adding NULL constituent model",BAD_DATA);
  }

  int c=(_nConstituents-1);//constit. index of current constituent

  //Add corresponding constituent storage state variables to model
  sv_type *aSV =new sv_type[_nWaterCompartments];
  int     *aLev=new int    [_nWaterCompartments];
  for(int ii=0;ii<_nWaterCompartments; ii++)
  {
    int m=ii+c*_nWaterCompartments;
    aSV [ii]=CONSTITUENT;
    aLev[ii]=m;
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

  delete[] aSV;
  delete[] aLev;

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
  for(int j=0; j<pModel->GetNumProcesses(); j++)
  {
    pProc=pModel->GetProcess(j);
    for(int q=0; q<pProc->GetNumConnections(); q++)
    {
      int iF=pProc->GetFromIndices()[q];
      int iT=pProc->GetToIndices()[q];
      sv_type typF=pModel->GetStateVarType(iF);
      sv_type typT=pModel->GetStateVarType(iT);
      if((iF!=iT) && (CStateVariable::IsWaterStorage(typF,false)) && (CStateVariable::IsWaterStorage(typT,false)))
      { //This is a process that may advect a constituent
        _nAdvConnections++;
      }
    }
  }

  //determine _iFromWater, _iToWater indices and process index for each and every connection
  //----------------------------------------------------------------------------
  js=0;
  int qq=0;
  _iFromWater=new int[_nAdvConnections];
  _iToWater  =new int[_nAdvConnections];
  _js_indices=new int[_nAdvConnections];
  for(int j=0; j<pModel->GetNumProcesses(); j++)
  {
    pProc=pModel->GetProcess(j);

    for(int q=0; q<pProc->GetNumConnections(); q++)
    {
      int iF=pProc->GetFromIndices()[q];
      int iT=pProc->GetToIndices()[q];
      sv_type typF=pModel->GetStateVarType(iF);
      sv_type typT=pModel->GetStateVarType(iT);
      if((iF!=iT) && (CStateVariable::IsWaterStorage(typF,false)) && (CStateVariable::IsWaterStorage(typT,false)))
      { //This is a process that may advect a constituent
        _iFromWater[qq]=iF;
        _iToWater  [qq]=iT;
        _js_indices[qq]=js;
        qq++;
      }
      else {
        //cout<<"NOT A WATER CONNECTION:"<<endl;
        //cout<<pProc->GetProcessType()<<"  from "<<CStateVariable::SVTypeToString(pModel->GetStateVarType(iF),0)<<" to "<<CStateVariable::SVTypeToString(pModel->GetStateVarType(iT),0)<<endl;
      }
      js++;
    }
  }

  //identify number of water storage compartments
  //----------------------------------------------------------------------------
  _nWaterCompartments=0;
  for(int i=0;i<pModel->GetNumStateVars(); i++)
  {
    if(CStateVariable::IsWaterStorage(pModel->GetStateVarType(i),false))
    {
      _nWaterCompartments++;
    }
  }
  // identify all state variables which are water storage compartments
  _iWaterStorage=new int[_nWaterCompartments];
  _aIndexMapping=new int[pModel->GetNumStateVars()];//current # of state variables does not include transport constituents
  int ii=0; //ii sifts through storage compartments ii=0.._nWaterCompartments
  for(int i=0;i<pModel->GetNumStateVars(); i++)
  {
    _aIndexMapping[i]=DOESNT_EXIST;
    if(CStateVariable::IsWaterStorage(pModel->GetStateVarType(i),false))
    {
      _iWaterStorage[ii]=i;
      _aIndexMapping[i ]=ii;
      ii++;
    }
  }
  _nIndexMapping=pModel->GetNumStateVars();


  if(_nConstituents==0) { return; }/// all of the above work necessary even with no transport. Not sure why.

  /// QA/QC
  //----------------------------------------------------------------------------
  // check for two constituents with same name
  for (int c = 0; c < _nConstituents; c++) {
    for (int cc = 0; cc < c; cc++) {
      if (_pConstitModels[c]->GetName() == _pConstitModels[cc]->GetName()) {
         ExitGracefully("CTransport::Prepare: a model cannot have two transported constituents with the same name",BAD_DATA_WARN);
      }
    }
  }

  //Synopsis
  //----------------------------------------------------------------------------
  if(!Options.silent) {
    cout<<"===TRANSPORT MODEL SUMMARY=============================="<<endl;
    cout<<"   number of compartments: "<<_nWaterCompartments<<endl;
    cout<<"   number of constituents: "<<_nConstituents<<endl;
    cout<<"   number of connections:  "<<_nAdvConnections<<endl;
    cout<<"   number of lat. connect.:"<<_nLatConnections<<endl;
    //cout<<"   number of source terms: "<<_nSources<<endl;
    cout<<"========================================================"<<endl;
    // TMP DEBUG below===================================================
    /*if (false){
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
    }*/
    // ====================================================================
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

  if(_nConstituents==0) { return; }

  //determine number of lateral connections
  _nLatConnections=0;
  for(int j=0; j<pModel->GetNumProcesses(); j++)
  {
    pProc=pModel->GetProcess(j);
    pLatProc=(CLateralExchangeProcessABC*)pProc;//re-cast
    if(pLatProc->GetNumLatConnections()>0) {
      //cout<<"CTransport: CalculateLateralConnections: process "<< j<<" #conns: " <<pLatProc->GetNumLatConnections()<<endl;
    }
    for(int q=0; q<pLatProc->GetNumLatConnections(); q++)
    {
      int iF=pLatProc->GetLateralFromIndices()[q];
      int iT=pLatProc->GetLateralToIndices()[q];
      sv_type typF=pModel->GetStateVarType(iF);
      sv_type typT=pModel->GetStateVarType(iT);
      if((iF!=iT) && (CStateVariable::IsWaterStorage(typF,false)) && (CStateVariable::IsWaterStorage(typT,false)))
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
      if((iF!=iT) && (CStateVariable::IsWaterStorage(typF,false)) && (CStateVariable::IsWaterStorage(typT,false)))
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
/// \brief Add geochem parameter
/// \param typ [in] geochemical parameter type
/// \param constit_name [in] name of constitutent
/// \param constit_name2 [in] name of 2nd constituent (or "" if not relevant)
/// \param procname [in] process name
/// \param watcompart [in] SV index i of water compartment state variable (or DOESNT_EXIST if parameter not linked to compartment)
/// \param pclass [in] parameter class name string (or "" if parameter not linked to class)
/// \param value [in] parameter value
//
void   CTransportModel::AddGeochemParam(const gparam_type typ,const string constit_name,const string constit_name2,
                                        const string procname,const int watcompart,const string pclass,const double& value)
{
  AddProcessName(procname);

  geochem_param* par=new geochem_param(typ);
  par->value=value;
  par->constit_ind=GetConstituentIndex(constit_name);
  if (constit_name2!=""){ par->constit_ind2=GetConstituentIndex(constit_name2); }
  par->class_name =pclass;
  par->process_ind=GetProcessIndex(procname);
  if(watcompart!=DOESNT_EXIST) {
    par->ii_water_stor= GetWaterStorIndexFromSVIndex(watcompart);
  }
  ExitGracefullyIf(par->constit_ind==DOESNT_EXIST,"AddGeochemParam: invalid constituent name in :GeochemParameter command",BAD_DATA_WARN);

  DynArrayAppend((void**&)(_pGeochemParams),(void*)(par),_nGeochemParams);
}

//////////////////////////////////////////////////////////////////
/// \brief add process name to process list
/// \param name [in] name of process
//
void   CTransportModel::AddProcessName(string name)
{
  if (GetProcessIndex(name)!=DOESNT_EXIST){return;}//name on list already

  string *tmp=new string[_nProcessNames+1];
  for(int i=0;i<_nProcessNames;i++) {
    tmp[i]=_aProcessNames[i];
  }
  tmp[_nProcessNames]=name;
  delete [] _aProcessNames;
  _aProcessNames=tmp;
  _nProcessNames++;
}
//////////////////////////////////////////////////////////////////
/// \brief Initialization of all transport variables
/// \note determines initial conditions for all constituents, initializes routing variables
/// \note  called after all constituents have been added by CModel::Initialize
//
void CTransportModel::Initialize(const optStruct &Options)
{
  for(int c=0;c<_nConstituents;c++){
    _pConstitModels[c]->Initialize(Options);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Initialization of all transport parameter variables
/// \note  called after all constituents have been added by CModel::Initialize but BEFORE .rvp file read
//
void CTransportModel::InitializeParams(const optStruct& Options)
{

}

//////////////////////////////////////////////////////////////////
/// \brief Increment cumulative mass/energy added to watershed
/// \param Options [in] Global model options structure
/// \param tt [in] time structure
//
void   CTransportModel::IncrementCumulInput(const optStruct &Options,const time_struct &tt)
{
  for(int c=0;c<_nConstituents; c++) {
    _pConstitModels[c]->IncrementCumulInput(Options,tt);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Increment cumulative mass/energy lost from watershed
/// \param Options [in] Global model options structure
//
void   CTransportModel::IncrementCumulOutput(const optStruct &Options)
{
  for(int c=0;c<_nConstituents; c++) {
    _pConstitModels[c]->IncrementCumulOutput(Options);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Write transport output file headers
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CTransportModel::WriteOutputFileHeaders(const optStruct &Options) const
{
  if (Options.output_format==OUTPUT_STANDARD){
    for(int c=0;c<_nConstituents;c++) {
      _pConstitModels[c]->WriteOutputFileHeaders(Options);
    }
  }
  else if (Options.output_format == OUTPUT_ENSIM) {
    for(int c=0;c<_nConstituents;c++) {
       _pConstitModels[c]->WriteEnsimOutputFileHeaders(Options);
    }
  }
  else if (Options.output_format == OUTPUT_NETCDF) {
    for(int c=0;c<_nConstituents;c++) {
       _pConstitModels[c]->WriteNetCDFOutputFileHeaders(Options);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output to file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams; called from CModel::WriteMinorOutput()
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
//
void CTransportModel::WriteMinorOutput(const optStruct &Options,const time_struct &tt) const
{
  if (Options.output_format==OUTPUT_STANDARD){
    for(int c=0;c<_nConstituents;c++) {
      _pConstitModels[c]->WriteMinorOutput(Options,tt);
    }
  }
  else if (Options.output_format == OUTPUT_ENSIM) {
    for(int c=0;c<_nConstituents;c++) {
       _pConstitModels[c]->WriteEnsimMinorOutput(Options,tt);
    }
  }
  else if (Options.output_format == OUTPUT_NETCDF) {
    for(int c=0;c<_nConstituents;c++) {
       _pConstitModels[c]->WriteNetCDFMinorOutput(Options,tt);
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Close transport output files
//
void CTransportModel::CloseOutputFiles() const
{
  for(int c=0;c<_nConstituents;c++){
    _pConstitModels[c]->CloseOutputFiles();
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Write major output
//
void  CTransportModel::WriteMajorOutput(ofstream &RVC) const
{
  for(int c=0;c<_nConstituents;c++) {
    _pConstitModels[c]->WriteMajorOutput(RVC);
  }
}
