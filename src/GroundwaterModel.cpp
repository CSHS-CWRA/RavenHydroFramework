#include "GroundwaterClass.h"
string FilenamePrepare(string filebase,const optStruct &Options); //Defined in StandardOutput.cpp

CGroundwaterModel::CGroundwaterModel(CModel *pModel) {
  _pModel=pModel;

  _nAquiferStacks=0; _pAQStacks=NULL;
  _pGWGeo=NULL;
  _nGWSP=0;          _pGWSP=NULL;
  _nOEC=0;           _pOEC=NULL;


}
CGroundwaterModel::~CGroundwaterModel() {
  int g;
  for(g=0;g<_nAquiferStacks; g++) { delete _pAQStacks[g]; } delete[] _pAQStacks;    _pAQStacks=NULL;
  delete _pGWGeo;       _pGWGeo=NULL;
  for(g=0;g<_nGWSP; g++) { delete _pGWSP[g]; }    delete[] _pGWSP;        _pGWSP =NULL;
  for(g=0;g<_nOEC; g++) { delete _pOEC[g]; }     delete[] _pOEC;         _pOEC  =NULL;
  //CGWGeometryClass::DestroyAllGWGeoClasses();  ////  NEED TO FIX!!!!  Throws an assertion error GWMIGRATE
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific GW Geometry denoted by index k
///
/// \param c [in] GW Geo index
/// \return pointer to GW Geometry corresponding to passed index k
//
CGWGeometryClass *CGroundwaterModel::GetGWGeom() const
{
  return _pGWGeo;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns num of GW Stress Period Classes
///
/// 
/// \return int number of GW Stress Period Classes
//
int CGroundwaterModel::GetNumGWSPs() const { return _nGWSP; }

//////////////////////////////////////////////////////////////////
/// \brief Returns specific GW Stress Period Class denoted by index c
///
/// \param c [in] GW Stress Period index
/// \return pointer to GW Stress Period corresponding to passed index c
//
CGWStressPeriodClass  *CGroundwaterModel::GetGWSPs(int c) const
{
#if _STRICTCHECK_
  ExitGracefullyIf((c<0) || (c>=_nGWSP),"CModel GetGWSPClass::improper index",BAD_DATA);
#endif
  return _pGWSP[c];
}
//////////////////////////////////////////////////////////////////
/// \brief Returns num of GW Overlap/Exchange Classes
///
/// 
/// \return int number of GW Overlap/Exchange Classes
//
int CGroundwaterModel::GetNumOEs() const { return _nOEC; }

//////////////////////////////////////////////////////////////////
/// \brief Returns specific GW Overlap/Exchange Class denoted by index c
///
/// \param c [in] GW Overlap/Exchange index
/// \return pointer to GW Overlap/Exchange corresponding to passed index c
//
COverlapExchangeClass  *CGroundwaterModel::GetOEs(int c) const
{
#if _STRICTCHECK_
  ExitGracefullyIf((c<0) || (c>=_nOEC),"CModel GetOEClass::improper index",BAD_DATA);
#endif
  return _pOEC[c];
}
//////////////////////////////////////////////////////////////////
/// \brief Adds Aquifer Stack to model
/// 
/// \param *pAqStacks [in] (valid) pointer to Aquifer Stack to be added to model
//
void CGroundwaterModel::AddAquiferStack(CAquiferStack *pAquiferStack)
{
  if(!DynArrayAppend((void**&)(_pAQStacks),(void*)(pAquiferStack),_nAquiferStacks)) {
    ExitGracefully("CModel::AddAquiferStack: adding NULL Aquifer Stack",BAD_DATA);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Adds GW Stress Period to model
/// 
/// \param *pGWSP [in] (valid) pointer to stress period to be added to model
//
void CGroundwaterModel::AddGWSP(CGWStressPeriodClass *pGWStressPeriod)
{
  if(!DynArrayAppend((void**&)(_pGWSP),(void*)(pGWStressPeriod),_nGWSP)) {
    ExitGracefully("CModel::AddGWSP: adding NULL Stress Period",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Adds GW Geometry to model
/// 
/// \param *pGWGeo [in] (valid) pointer to GW Geometry to be added
//
void CGroundwaterModel::SetGWGeometry(CGWGeometryClass *pGWGeometry)
{
  _pGWGeo=pGWGeometry;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds Overlap/Exchange Class to model
/// 
/// \param *pOEc [in] (valid) pointer to Overlap/Exchange to be added
//
void CGroundwaterModel::AddOE(COverlapExchangeClass *pOEin)
{
  if(!DynArrayAppend((void**&)(_pOEC),(void*)(pOEin),_nOEC)) {
    ExitGracefully("CModel::AddOE: adding NULL Overlap/Exchange",BAD_DATA);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets Number of Aquifer Stacks layers to model
/// 
/// \param int nLayers [in] in aquifer stacks
//
void CGroundwaterModel::SetNumAquiferStack(int nAQLayers)
{
  //_nAquiferLayers = nAQLayers; //GWMIGRATE - issues here.
}
///////////////////////////////////////////////////////////////////
/// \brief Converts initial conditions of groundwater storage to hydraulic head
///
//
void CGroundwaterModel::ConvertInitialConditions(CModel *&pModel)
{
  int iGW;
  int gws = 0;
  int gwh = 0;
  int nHRUs     =pModel->GetNumHRUs();
  double poro,bot;
  double head,stor,gw_stor,gw_head;
  CHydroUnit *pHRU;
  const  soil_struct *S;
  string param;

  for(int k=0;k<nHRUs;k++)
  {
    pHRU  = pModel->GetHydroUnit(k);

    iGW   = pModel->GetStateVarIndex(GROUNDWATER,0);          //*****FIX - needs to work with multiple layers
    gw_stor = pHRU->GetStateVarValue(iGW);
    gw_head = _pGWGeo->GetGWGeoProperty("GW_HEAD",k,0);    //*****FIX - needs to work with multiple layers

    if(gw_stor != 0) { gws++; }
    if(gw_head != 0) { gwh++; }
  }

  if(gws > 0 && gwh == 0)     //if initial conditions have storage but no heads
  {
    for(int k=0;k<nHRUs;k++)
    {
      pHRU = pModel->GetHydroUnit(k);
      S    = pHRU->GetAquiferProps(0);
      poro = S->porosity;
      bot  = _pGWGeo->GetGWGeoProperty("BOT_ELEV",k,0);

      iGW  = pModel->GetStateVarIndex(GROUNDWATER,0);         //*****FIX - needs to work with multiple layers
      stor = pHRU->GetStateVarValue(iGW);
      head = ((stor / poro) / MM_PER_METER) + bot;
      param = "GW_HEAD";
      _pGWGeo->SetGWGeoProperty(param,head,k,0);        //*****FIX - needs to work with multiple layers
    }
  }
  else if(gws == 0 && gwh > 0)    //if initial conditions have heads but no storage
  {
    for(int k=0;k<nHRUs;k++)
    {
      pHRU = pModel->GetHydroUnit(k);
      S    = pHRU->GetAquiferProps(0);
      poro = S->porosity;
      bot  = _pGWGeo->GetGWGeoProperty("BOT_ELEV",k,0);

      head = _pGWGeo->GetGWGeoProperty("GW_HEAD",k,0);      //*****FIX - needs to work with multiple layers
      stor = (((head - bot) * MM_PER_METER) * poro);
      iGW  = pModel->GetStateVarIndex(GROUNDWATER,0);         //*****FIX - needs to work with multiple layers
      pHRU->SetStateVarValue(iGW,stor);
    }
  }
}

void CGroundwaterModel::WriteOutputFileHeaders(const optStruct &Options) 
{
  string tmpFilename;

  //GroundwaterHeads.csv
  //--------------------------------------------------------------
  if(Options.write_gwhead)
  {
    ofstream _GWHEAD;
    tmpFilename=FilenamePrepare("GroundwaterHeads.csv",Options);
    _GWHEAD.open(tmpFilename.c_str());

    _GWHEAD<<"time[d],date,hour";
    for(int i = 0;i<_pGWGeo->GetGWGeoStruct()->num_nodes;i++)
    {
      _GWHEAD<<", Cell "<<i+1 <<" [m]";
    }
    _GWHEAD<<endl;
  }

  //GroundwaterFlows.csv
  //--------------------------------------------------------------
  if(Options.write_gwflow)
  {
    ofstream _GWFLOW;
    tmpFilename=FilenamePrepare("GroundwaterFlows.csv",Options);
    _GWFLOW.open(tmpFilename.c_str());

    _GWFLOW<<"time[d],date,hour";
    for(int i = 0;i<_pGWGeo->GetGWGeoStruct()->num_nodes;i++)
    {
      for(int j = 0; j<_pGWGeo->GetGWGeoStruct()->num_cell_connections[i]-1;j++)
      {
        _GWFLOW<<","<<i+1<<" -> "<<_pGWGeo->GetGWGeoStruct()->cell_connections[i][j]<<" [mm/d]";
      }
    }
    _GWFLOW<<endl;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Writes minor output to file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time *at the end of* the pertinent time step
//
void CGroundwaterModel::WriteMinorOutput(const optStruct &Options,const time_struct &tt)
{
  string tmpFilename;
  string thisdate=tt.date_string;
  string thishour=DecDaysToHours(tt.julian_day);

  //Write groundwater heads to GroundwaterHeads.csv
  //----------------------------------------------------------------
  if(Options.write_gwhead)
  {
    _GWHEAD<<tt.model_time<<","<<thisdate<<","<<thishour;

    for(int k = 0; k<_pModel->GetNumHRUs();k++)
    {
      _GWHEAD<<","<<_pGWGeo->GetGWGeoProperty("GW_HEAD",k,0);
    }
    _GWHEAD<<endl;
  }

  //Write groundwater flows to GroundwaterFlows.csv
  //----------------------------------------------------------------
  if(Options.write_gwflow)
  {
    double h_i,h_j,c_ij;
    _GWFLOW.open(tmpFilename.c_str(),ios::app);

    _GWFLOW<<tt.model_time<<","<<thisdate<<","<<thishour<<",";     
    
    //need to send from GW SOLVER
    for(int i = 0;i<_pGWGeo->GetGWGeoStruct()->num_nodes;i++)
    {
      h_i=_pGWGeo->GetGWGeoProperty("GW_HEAD",i,0);
      for(int j = 0; j<_pGWGeo->GetGWGeoStruct()->num_cell_connections[i]-1;j++)
      {
        h_j=0.0;//_pGWGeo->GetGWGeoProperty("GW_HEAD",_pGWGeo->GetGWGeoStruct()->cell_connections[i][j],0); //JRC : maybe?
        c_ij=0.0; //GetConductance ???
        _GWFLOW<<","<<c_ij*(h_i-h_j);
      }
    }
  }
}
void CGroundwaterModel::CloseOutputFiles()
{
  if(_GWHEAD.is_open()) { _GWHEAD.close(); }
  if(_GWFLOW.is_open()) { _GWFLOW.close(); }
}