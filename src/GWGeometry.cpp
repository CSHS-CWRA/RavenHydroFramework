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
/// \brief Implementation of the groundwater geometry class constructor
/// \param name [in] String nickname for groundwater geometry profile 
//
CGWGeometryClass::CGWGeometryClass(const string name)
{
  tag=name;
   if (!DynArrayAppend((void**&)(pAllGWGeoClasses),(void*)(this),NumGWGeoClasses)){
     ExitGracefully("CGWGeometry::Constructor: creating NULL groundwater geometry class",BAD_DATA);};
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CGWGeometryClass::~CGWGeometryClass(){
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING GROUNDWATER GEOMETRY CLASS "<<endl;}
}

///////////////////////////////////////////////////////////////////
/// \brief Returns reference to groundwater geometry properties
/// \return Groundwater geometry properties associated with groundwater class
//
const gw_dis_struct *CGWGeometryClass::GetGWGeoStruct() const{return &W;}

///////////////////////////////////////////////////////////////////
/// \brief Return nick name identifier of groundwater geometry class
/// \return nick name identifier of groundwater geometry class
//
string                CGWGeometryClass::GetTag					 () const{return tag;}

/*****************************************************************
   Static Initialization, Accessors, Destructors
*****************************************************************/
CGWGeometryClass **CGWGeometryClass::pAllGWGeoClasses=NULL; 
int                CGWGeometryClass::NumGWGeoClasses=0;     


//////////////////////////////////////////////////////////////////
/// \brief Return number of groundwater geometry classes
/// \return Number of groundwater classes
//
int CGWGeometryClass::GetNumGWGeoClasses(){
  return NumGWGeoClasses;
}


//////////////////////////////////////////////////////////////////
/// \brief Summarize groundwater geometry class information to screen
//
void CGWGeometryClass::SummarizeToScreen()
{
  string base;
  //cout<<"==================="<<endl;
  cout<<"Groundwater Geometry Summary: "<<endl;
  cout<<"      " << NumGWGeoClasses <<" groundwater geometry in database"<<endl;
  for (int c=0; c<NumGWGeoClasses;c++){
    cout<<"   Groundwater Geometry class \""<<pAllGWGeoClasses[c]->GetTag()<<"\": "<<endl;
    cout<<"             Number of Layers: "<<pAllGWGeoClasses[c]->GetGWGeoStruct()->num_layers<<endl;
    cout<<"     Number of GW Cells/Layer: "<<pAllGWGeoClasses[c]->GetGWGeoStruct()->node_layer<<endl;
    cout<<"     Number of GW Connections: "<<pAllGWGeoClasses[c]->GetGWGeoStruct()->num_connections<<endl;

    if((pAllGWGeoClasses[c]->GetGWGeoStruct()->basement)==0){base = "Confined";} else{base = "Unconfined";}
    cout<<"               GW Basement is: "<< base <<endl;
    cout<<endl;
  }
  return;
}


//////////////////////////////////////////////////////////////////
/// \brief Write groundwater geometry class properties to file
/// \param &OUT [out] Output stream to which information is written
/// \todo [add funct] repair such that only parameters that are used by model are written to parameters.csv.
/// one way to do this would be to store the template groundwater structure with class, and use NOT_NEEDED status of variables as indicator.
//
void CGWGeometryClass::WriteParamsToFile(ofstream &OUT)
{
  CGWGeometryClass *s;
  const gw_dis_struct *w;
  OUT<<endl<<"---Groundwater Geometry Class Parameters---------------------"<<endl;
  OUT<<"CLASS,";
  OUT<<"NODES,LAYERS,CONNECTIONS,VERT_DIS,MATRIX,ORIGIN,";
  OUT<<"BASEMENT,NODE/LAYER,TOP_ELEV,BOT_ELEV,CELL_AREA,";
  OUT<<"NUM_CELL_CONNECTIONS,CELL_CONNECTIONS,VERT_CELL_CONN,";
  OUT<<"NODE_INTERFACE_LENGTH,FACE_WIDTH,INTERFACE_AREA";
  OUT<<endl;
  for (int c=0; c<NumGWGeoClasses;c++)
  {
    s=pAllGWGeoClasses[c];
    w=s->GetGWGeoStruct();
    OUT<<s->GetTag()<<",";
	OUT << w->num_nodes << "," << w->num_layers << "," << w->num_connections << ",";
	OUT << w->vert_dis << "," << w->matrix_format << "," << w->grid_origin << ",";
	OUT << w->basement << "," << w->node_layer << "," << w->top_elev << "," << w->bot_elev << ",";
	OUT << w->cell_area << "," << w->num_cell_connections << "," << w->cell_connections << "," << w->vert_cell_connections << ",";
    OUT<<w->node_interface_length <<","<<w->face_width          <<","<<w->interface_area      <<",";
    OUT<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all groundwater geometry classes
//
void CGWGeometryClass::DestroyAllGWGeoClasses()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL GROUNDWATER GEOMETRY CLASSES"<<endl;}
  for (int c=0; c<NumGWGeoClasses;c++){
    delete pAllGWGeoClasses[c];
  }
  delete [] pAllGWGeoClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Allocate struct array sizes
//
void CGWGeometryClass::AllocateArrays(int nodes)
{
  AllocateArrays(W,nodes);
}

//////////////////////////////////////////////////////////////////
/// \brief Allocate struct array sizes
//
void CGWGeometryClass::AllocateArrays(gw_dis_struct &W, int nNodes)
{
  W.top_elev             = new double[nNodes];
  W.bot_elev             = new double[nNodes];
  W.min_topog            = new double[nNodes];
  W.cell_area            = new double[nNodes];
  W.num_cell_connections = new double[nNodes];
  W.groundwater_head     = new double[nNodes];

  W.cell_connections     = new double*[nNodes];
  W.vert_cell_connections= new double*[nNodes];
  W.node_interface_length= new double*[nNodes];
  W.face_width           = new double*[nNodes];            //****FIX - ALL OF THIS NEEDS TO BE DELETED!
  W.interface_area       = new double*[nNodes];

  for(int i = 0; i < nNodes; i++)
  {
    W.cell_connections[i]     = new double[10];
    W.vert_cell_connections[i]= new double[10];
    W.node_interface_length[i]= new double[10];
    W.face_width[i]           = new double[10];
    W.interface_area[i]       = new double[10];
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the groundwater class corresponding to passed string
/// \details Converts string to gwclass
///  if string is invalid, returns NULL
/// \param w [in] GW class identifier (tag or index)
/// \return Reference to groundwater class corresponding to identifier w
//
CGWGeometryClass *CGWGeometryClass::StringToGWGeoClass(const string s)
{
  string sup=StringToUppercase(s);
  for (int c=0;c<NumGWGeoClasses;c++)
  {
    if (!sup.compare(StringToUppercase(pAllGWGeoClasses[c]->GetTag()))){return pAllGWGeoClasses[c];}
    else if (s_to_i(s.c_str())==(c+1))           {return pAllGWGeoClasses[c];}
  }
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns the groundwater geometry class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] GW class index 
/// \return Reference to groundwater geometry class corresponding to index c
//
const CGWGeometryClass *CGWGeometryClass::GetGWGeoClass(int c)
{
  if ((c<0) || (c>=CGWGeometryClass::NumGWGeoClasses)){return NULL;}
  return pAllGWGeoClasses[c];
}

//////////////////////////////////////////////////////////////////
/// \brief Sets default groundwater geometry properties
/// \param &W [out] Groundwater geometry properties class
/// \param is_template [in] True if the default value being set is for the template class
//
void CGWGeometryClass::InitializeGWGeoProperties(gw_dis_struct &W, bool is_template)//static
{
  //W.ID                    ="";
  W.num_nodes             = 1;   
  W.num_layers            = 1;    
  W.num_connections       = 1;    
  W.vert_dis              =-1;          
  W.matrix_format         = 0;  
  W.grid_origin           = 1;   

  W.basement              = 0;   
  W.node_layer            = 1;
  W.top_elev              = NULL;                
  W.bot_elev              = NULL;
  W.min_topog             = NULL;
  W.cell_area             = NULL;
  
  W.num_cell_connections  = NULL;
  W.cell_connections      = NULL;
  W.vert_cell_connections = NULL;

  W.node_interface_length = NULL;
  W.face_width            = NULL;
  W.interface_area        = NULL; 

  W.groundwater_head      = NULL;

}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the groundwater geometry property corresponding to param_name
/// \param param_name [in] Parameter identifier
//
void  CGWGeometryClass::SetGWGeoProperty(string &param_name, 
                                          const double value)
{
  SetGWGeoProperty(W,param_name,value);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the groundwater geometry property corresponding to param_name
/// \param param_name [in] Parameter identifier
//
void  CGWGeometryClass::SetGWGeoProperty(string &param_name, 
                                          const double value,int i, int j)
{
  SetGWGeoProperty(W,param_name,value, i, j);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the gw geometry property corresponding to param_name
/// \note This is declared as a static member because groundwater class 
/// is not instantiated prior to read of .rvd file
/// \param &W [out] Groundwater properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  CGWGeometryClass::SetGWGeoProperty(gw_dis_struct &W, 
                                          string     param_name, 
                                          const double value)
{
  string name;
  name = StringToUppercase(param_name);

  //if      (!name.compare("ID"                     )){W.ID=value;}
  if      (!name.compare("NUM_NODES"              )){W.num_nodes=value;}
  else if (!name.compare("NUM_LAYERS"             )){W.num_layers=value;}
  else if (!name.compare("NUM_CONNECTIONS"        )){W.num_connections=value;}
  else if (!name.compare("VERT_DIS"               )){W.vert_dis=value;}
  else if (!name.compare("MATRIX_FORMAT"          )){W.matrix_format=value;}
  else if (!name.compare("GRID_ORIGIN"            )){W.grid_origin=value;}
  else if (!name.compare("BASEMENT"               )){W.basement=value;}
  else if (!name.compare("NODES_LAYER"            )){W.node_layer=value;}
  else{
    WriteWarning("CGWGeometryClass::SetGWGeoProperty: Unrecognized/invalid groundwater geometry parameter name ("+name+") in .rvd file",false);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the gw geometry property corresponding to param_name
/// \note This is declared as a static member because groundwater class 
/// is not instantiated prior to read of .rvd file
/// \param &W [out] Groundwater properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  CGWGeometryClass::SetGWGeoProperty(gw_dis_struct &W, 
                                          string     param_name, 
                                          const double value,
                                          int i, int j)
{
  string name;
  name = StringToUppercase(param_name);
  
  if      (!name.compare("TOP_ELEV"               )){W.top_elev[i]=value;}
  else if (!name.compare("BOT_ELEV"               )){W.bot_elev[i]=value;}
  else if (!name.compare("MIN_TOPOG"              )){W.min_topog[i]=value;}
  else if (!name.compare("CELL_AREA"              )){W.cell_area[i]=value;}
  else if (!name.compare("NUM_CELL_CONNECTIONS"   )){W.num_cell_connections[i]=value;}
  else if (!name.compare("CELL_CONNECTIONS"       )){W.cell_connections[i][j]=value;}
  else if (!name.compare("VERT_CELL_CONNECTIONS"  )){W.vert_cell_connections[i][j]=value;}
  else if (!name.compare("NODE_INTERFACE_LENGTH"  )){W.node_interface_length[i][j]=value;}
  else if (!name.compare("FACE_WIDTH"             )){W.face_width[i][j]=value;}
  else if (!name.compare("INTERFACE_AREA"         )){W.interface_area[i][j]=value;}
  else if (!name.compare("GW_HEAD"                )){W.groundwater_head[i]=value;}
  else{
    WriteWarning("CGWGeometryClass::SetGWGeoProperty: Unrecognized/invalid groundwater geometry parameter name ("+name+") in .rvd file",false);
  }
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns gw geometry property value corresponding to param_name
/// \param param_name [in] Parameter name
/// \return GW properties corresponding to parameter name
//
double CGWGeometryClass::GetGWGeoProperty(string param_name) const
{
  return GetGWGeoProperty(W,param_name);
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns gw geometry property value corresponding to param_name
/// \param param_name [in] Parameter name
/// \return GW properties corresponding to parameter name
//
double CGWGeometryClass::GetGWGeoProperty(string param_name,int i,int j) const
{
  return GetGWGeoProperty(W,param_name,i,j);
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns groundwater geometry property value corresponding to param_name from structure provided
/// \param gw_struct [in] Groundwater discretization structure
/// \param param_name [in] Parameter name
/// \return Groundwater discretization corresponding to parameter name
//
double CGWGeometryClass::GetGWGeoProperty(const gw_dis_struct &W, string param_name)
{
  string name;
  name = StringToUppercase(param_name);

  //if      (!name.compare("ID"                     )){return W.ID;}
  if      (!name.compare("NUM_NODES"              )){return W.num_nodes;}
  else if (!name.compare("NUM_LAYERS"             )){return W.num_layers;}
  else if (!name.compare("NUM_CONNECTIONS"        )){return W.num_connections;}
  else if (!name.compare("VERT_DIS"               )){return W.vert_dis;}
  else if (!name.compare("MATRIX_FORMAT"          )){return W.matrix_format;}
  else if (!name.compare("GRID_ORIGIN"            )){return W.grid_origin;}
  else if (!name.compare("BASEMENT"               )){return W.basement;}
  else if (!name.compare("NODES_LAYER"            )){return W.node_layer;}
  else{
    string msg="CGWGeometryClass::GetGWGeoProperty: Unrecognized/invalid groundwater geometry parameter name in .rvd file: "+name;
    ExitGracefully(msg.c_str(),BAD_DATA);
    return 0.0;
  }
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns groundwater geometry property value corresponding to param_name from structure provided
/// \param gw_struct [in] Groundwater structure
/// \param param_name [in] Parameter name
/// \return Groundwater stress period properties corresponding to parameter name
//
double CGWGeometryClass::GetGWGeoProperty(const gw_dis_struct &W, string param_name, int i, int j)
{
  //\todo[optimize] : passing strings here is really expensive; should use enumerated type
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("TOP_ELEV"               )){return W.top_elev[i];}
  else if (!name.compare("BOT_ELEV"               )){return W.bot_elev[i];}
  else if (!name.compare("MIN_TOPOG"              )){return W.min_topog[i];}
  else if (!name.compare("CELL_AREA"              )){return W.cell_area[i];}
  else if (!name.compare("NUM_CELL_CONNECTIONS"   )){return W.num_cell_connections[i];}
  else if (!name.compare("CELL_CONNECTIONS"       )){return W.cell_connections[i][j];}
  else if (!name.compare("VERT_CELL_CONNECTIONS"  )){return W.vert_cell_connections[i][j];}
  else if (!name.compare("NODE_INTERFACE_LENGTH"  )){return W.node_interface_length[i][j];}
  else if (!name.compare("FACE_WIDTH"             )){return W.face_width[i][j];}
  else if (!name.compare("INTERFACE_AREA"         )){return W.interface_area[i][j];}
  else if (!name.compare("GW_HEAD"                )){return W.groundwater_head[i];}
  else{
    string msg="CGWGeometryClass::GetGWGeoProperty: Unrecognized/invalid groundwater geometry parameter name in .rvd file: "+name;
    ExitGracefully(msg.c_str(),BAD_DATA);
    return 0.0;
  }
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns aquifer cell thickness corresponding to indices
/// \indice cell index [in] cell index
/// \return cell thickness
//
double CGWGeometryClass::GetAQThickness(int index)
{
  double thickness;
  thickness = (pAllGWGeoClasses[0]->GetGWGeoStruct()->top_elev[index] - pAllGWGeoClasses[0]->GetGWGeoStruct()->bot_elev[index]) * MM_PER_METER;

  return thickness;
}

