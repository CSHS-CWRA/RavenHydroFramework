/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
------------------------------------------------------------------
  Class CGroundwater
----------------------------------------------------------------*/
#ifndef GROUNDWATER_H
#define GROUNDWATER_H

#include "RavenInclude.h"
#include "SparseMatrix.h"
#include "Model.h"

class CModel;
class CAquiferStack;

///////////////////////////////////////////////////////////////////
/// \brief Structure that, when instantiated, contains defining information that reflects a specific gw stress period type.
/// \details GW stress period properties calculated. All properties are static
/// and should not change over time.
// GWMIGRATE - should be moved to header of stress period class 
// GWMIGRATE - why are some of these referring to a stress period - some of these are cell properties?
//
struct sp {
  //ET parameters
  double numrec_ET;         ///<        number of cells being given ET parameters
  double *numcell_ET;       ///<         cell number of current parameter
  double *pet_max;          ///< [mm/d]  PET MAX for cell
  double *ext_depth;        ///< [m]     extinction depth for cell

                            //Drain parameters
  double numrec_D;          ///<        number of cells being given drain parameters
  double *numcell_Drains;   ///<         cell number of current parameter
  double *drain_elev;       ///< [m]     elevation of drain
  double *conductance;      ///< [mm/d]  conductance of drain

                            //Recharge parameters
  double numrec_R;          ///<        number of cells being given recharge parameters
  double *numcell_Rech;     ///<         cell number of current parameter
  double *recharge;         ///< [mm/d]  recharge rate for cell

                            //Additional GW params
  double numrec_K;          ///<         number of cells being given hydraulic conductivity values
  double *numcell_K;        ///<        cell number of current parameter
  double *hyd_cond;         ///< [mm/d] hydraulic conductivity value for cell

  double numrec_S;          ///<         number of cells being given storage values
  double *numcell_S;        ///<        cell number of current parameter
  double *ss;               ///< [1/m]   specific storage value for cell
};
struct gw_struct
{
  double stressperiod;        ///< stress period number
                              //Required Parameters
  sp sp_props;
};

///////////////////////////////////////////////////////////////////
/// \brief Structure that, when instantiated, contains defining information that reflects a gw geometry.
/// \details GW discretization properties calculated. All properties are static
/// and should not change over time.
//
struct gw_dis_struct 
{
  double num_nodes;                      ///<       number of nodes/cells in unstructured grid
  double num_layers;                     ///<       number of layers in model
  double num_connections;                ///<       total number of connections between cells
 
  double vert_dis;                       ///<       type of vertical discretization (0=no sub-dis, 1=sub-dis, -1=no sub-dis or horiz sub-dis)
  double matrix_format;                  ///<       type of matrix supplied (0=full, 1=upper tri)
  double grid_origin;                    ///<       location of grid cell 1 (1=bottom left, 0=top right) //JRC:why is this a double?
  double basement;                       ///<       indicates if bottom of model is confined or unconfined (0=confined, 1=unconfined)
  double node_layer;                     ///<       number of nodes per layer

  double *top_elev;                      ///< [m]   top elevations of cells
  double *bot_elev;                      ///< [m]   bottom elevations of cells
  double *min_topog;                     ///< [m]   minimum topographic elevation in a cell
  double *cell_area;                     ///< [m^2] area of cells

  //double *hyd_cond; 
  //double *specific_stor;               //JRC: I would store these here...
  //double *specific_yield;
  //

  double *num_cell_connections;          ///<       number of connections/node + 1
  double **cell_connections;             ///<       list of cells connected //JRC: are these indices of adjacent cells?
  double **vert_cell_connections;        ///<       indicates if cell_connections are vertical or horizontal (1=vertical, 0=horizontal)

  double **node_interface_length;        ///< [m]   length between node and interface between connecting nodes
  double **face_width;                   ///< [m]   width of interface between 2 connected cells
  double **interface_area;               ///< [m^2] area of interface between cells 

  double *groundwater_head;              ///< [m]   groundwater heads stored in metres [GWMIGRATE - HOW IS THIS STATIC!!?]
};

///////////////////////////////////////////////////////////////////
/// \brief Structure that, when instantiated, contains defining information that reflects a specific overlap/exchange type.
/// \details Overlap/Exchange properties calculated. All properties are static
/// and should not change over time.
//
struct oe {

  string ET_type;         ///< ET relationship type (ie. POWER LAW)
  string D_type;          ///< Drain relationship type (ie. POWER LAW)
  string SA_type;         ///< Saturated Area relationship type (ie. POWER LAW)

  int RelationshipCount;  ///< Number of  UFR relationships in file

  double *ET_a;            ///< Coefficient A of power law relationship for ET
  double *ET_b;            ///< Power B of power law relationship for ET
  double *D_a;             ///< Coefficient A of power law relationship for Drains
  double *D_b;             ///< Power B of power law relationship for Drains
  double *SA_a;            ///< Coefficient A of power law relationship for Saturated Area
  double *SA_b;            ///< Power B of power law relationship for Saturated Area



};
struct oe_struct
{
  int nHRU;            ///<  number of HRUs
  int nGWCells;        ///<  number of GW cells in connection with surface soil layers

  double **overlap;     ///< array of percentage breakdown between HRUs and GW cells

  int cell;            ///< cell that this type is applied to (ie. ALL, 1, 2, ...)
                       //Required Parameters
  oe oe_props;
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for groundwater stress period classification
//
class CGWStressPeriodClass
{
  protected:/*----------------------------------------------------*/

    string                   tag;                 ///< nickname for groundwater
    gw_struct                W;                   ///< corresponding properties for gw

    static CGWStressPeriodClass **pAllGWClasses;     ///< array of pointers to all gw classes that have been created
    static int                    NumGWClasses;      ///< number of gw classes

  public:/*-------------------------------------------------------*/
    //Constructors:
    CGWStressPeriodClass(const string name);
   ~CGWStressPeriodClass();

    //Accessors
    string                          GetTag() const;
    const gw_struct                *GetGWStruct() const;
    double                          GetGWProperty(string param_name, int i) const;
    void                            SetGWProperty(string &param_name, const double value,int i); 
       
    static int                          GetNumGWSPClasses();
    static const CGWStressPeriodClass  *GetGWSPClass(int c);
    static       CGWStressPeriodClass  *StringToGWClass(const string s);
    static void                         DestroyAllGWClasses();
    static void                         SetGWProperty         (gw_struct &W, string param_name, const double value,int i);
    static double                       GetGWProperty         (const gw_struct &W, string param_name, int i);
    static void                         InitializeGWProperties(gw_struct &W, bool is_template);

    static void                     SummarizeToScreen();
    static void                     WriteParamsToFile(ofstream &PARAMS);

    void                            AllocateArrays(string &param_name, int nodes);
    void                            AllocateArrays(gw_struct &W, string param_name, int nodes);
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for groundwater stress period classification
//
class CGWGeometryClass
{
  protected:/*----------------------------------------------------*/

    string                   tag;                ///< nickname for groundwater
    gw_dis_struct            W;                  ///< corresponding discretization for gw

    static CGWGeometryClass **pAllGWGeoClasses;     ///< array of pointers to all gw classes that have been created
    static int                NumGWGeoClasses;      ///< number of gw classes

  public:/*-------------------------------------------------------*/
    //Constructors:
    CGWGeometryClass(const string name);
   ~CGWGeometryClass();

    //Accessors
    string                          GetTag() const;
    const gw_dis_struct            *GetGWGeoStruct() const;
    double                          GetGWGeoProperty(string param_name) const;
    double                          GetGWGeoProperty(string param_name,int i,int j) const;
    void                            SetGWGeoProperty(string &param_name, const double value);
    void                            SetGWGeoProperty(string &param_name, const double value, int i, int j);
       
    static int                      GetNumGWGeoClasses();
    static const CGWGeometryClass  *GetGWGeoClass(int c);
    static       CGWGeometryClass  *StringToGWGeoClass(const string s);
    static void                     DestroyAllGWGeoClasses();

    static void                     SetGWGeoProperty         (gw_dis_struct &W, string param_name, const double value);
    static void                     SetGWGeoProperty         (gw_dis_struct &W, string param_name, const double value, int i, int j);
    static double                   GetGWGeoProperty         (const gw_dis_struct &W, string param_name);
    static double                   GetGWGeoProperty         (const gw_dis_struct &W, string param_name, int i, int j);
    static void                     InitializeGWGeoProperties(gw_dis_struct &W, bool is_template);

    static void                     SummarizeToScreen();
    static void                     WriteParamsToFile(ofstream &PARAMS);

    void                            AllocateArrays(int nodes);
    void                            AllocateArrays(gw_dis_struct &W, int nodes);

    double                          GetAQThickness(int index);
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for groundwater stress period classification
//
class COverlapExchangeClass
{
  protected:/*----------------------------------------------------*/

    string                   tag;                 ///< nickname for groundwater
    oe_struct                E;                   ///< corresponding properties for gw

    static COverlapExchangeClass **pAllOEClasses;     ///< array of pointers to all gw classes that have been created
    static int                    NumOEClasses;      ///< number of gw classes

  public:/*-------------------------------------------------------*/
    //Constructors:
    COverlapExchangeClass();
   ~COverlapExchangeClass();

    //Accessors
    string                          GetTag() const;
    const oe_struct                *GetOEStruct() const;
    double                          GetOEProperty(string param_name) const;
    string                          GetOEProperty(string param_name, string param_name1, int i);
    double                          GetOEProperty(string param_name, int i, int j);
    double                          GetOEProperty(string param_name, int i);
    void                            SetOEProperty(string param_name, const double value,int i, int j); 
    void                            SetOEProperty(string param_name, const double value);
    void                            SetOEProperty(string param_name, const double value, int i);
    void                            SetOEProperty(string param_name, const int value);
    void                            SetOEProperty(string param_name, string value);
       
    static int                          GetNumOEClasses();
    static const COverlapExchangeClass  *GetOEClass(int c);
    static       COverlapExchangeClass  *StringToOEClass(const string s);
    static void                         DestroyAllOEClasses();
    static void                         SetOEProperty         (oe_struct &E, string param_name, const double value);
    static void                         SetOEProperty         (oe_struct &E, string param_name, const string value);
    static void                         SetOEProperty         (oe_struct &E, string param_name, const double value,int i);
    static void                         SetOEProperty         (oe_struct &E, string param_name, const double value,int i,int j);
    static double                       GetOEProperty         (const oe_struct &E, string param_name);
    string                              GetOEProperty         (const oe_struct &E, string param_name, string param_name1, int i);
    static double                       GetOEProperty         (const oe_struct &E, string param_name, int i,int j);
    static double                       GetOEProperty         (const oe_struct &E, string param_name, int i);
    static void                         InitializeOEProperties(oe_struct &E, bool is_template);

    static void                     SummarizeToScreen();
    static void                     WriteParamsToFile(ofstream &PARAMS);

    void                            AllocateOEArrays(int nodes);
    void                            AllocateOEArrays(oe_struct &E, int nodes);
};
///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for groundwater model
//
class CGroundwaterModel
{
private:/*----------------------------------------------------*/

  CModel               *_pModel;

  //int                 _nNodes;  ///< JRC should be stored here in addition to geometry class
  //double   **_aOverlapWeights;  ///< 2D matrix of weights relating which HRUs are connected to which GW cells 
                                  ///< =A_ki/A_i where A_ki is area of overlap between HRU k and cell i and A_i is area of cell i

  //double             *_aHeads;  ///< Current groundwater heads [size: _nNodes]

  int           _nAquiferStacks;  ///< number of aquifer layers  //JRC : should be removed 

  CAquiferStack    **_pAQStacks;  ///< array of pointers to aquifer stacks
  CGWGeometryClass     *_pGWGeo;  ///< pointer to gw geometry 
  int                    _nGWSP;  ///< number of stress periods
  CGWStressPeriodClass **_pGWSP;  ///< array of pointers to stress period information
  int                     _nOEC;  ///< number of overlap/exchange classes
  COverlapExchangeClass **_pOEC;  ///< array of pointers to overlap/exchange information

  ofstream               _GWHEAD; ///< output file stream for GroundwaterHeads.csv
  ofstream               _GWFLOW; ///< output file stream for GroundwaterFlows.csv

public:/*-------------------------------------------------------*/

       //Constructors:
  CGroundwaterModel(CModel *pModel);
  ~CGroundwaterModel();

  CGWGeometryClass      *GetGWGeom() const;  
  int                    GetNumGWSPs() const;
  CGWStressPeriodClass  *GetGWSPs(int c) const;
  int                    GetNumOEs() const;
  COverlapExchangeClass *GetOEs(int c) const;

  //int GetNodes() const;
  //void SetOverlapWeights (const double **aOverlapWts, const int nCells, const int nHRUs);

  void    ConvertInitialConditions  (CModel *&pModel);

  void    SetNumAquiferStack        (int nAQLayers);
  void    AddAquiferStack           (CAquiferStack             *pAquiferStack);
  void    AddGWSP                   (CGWStressPeriodClass      *pGWStressPeriod);
  void    SetGWGeometry             (CGWGeometryClass          *pGWGeometry);
  void    AddOE                     (COverlapExchangeClass     *pOEin);

  void    WriteOutputFileHeaders    (const optStruct &Options);
  void    WriteMinorOutput(const optStruct &Options,const time_struct &tt);
  void    CloseOutputFiles ();
};
#endif