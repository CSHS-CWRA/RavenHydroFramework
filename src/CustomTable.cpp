/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/

#include "CustomOutput.h"

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CCustomTable constructor
/// \param filename          [in] specified filename
/// \param pp                [in] Subbasin group index
/// \param *pMod             [in] Pointer to model
//
CCustomTable::CCustomTable(string filename, int pp, const CModel* pMod)
{
  _filename=filename;
  _sb_grp_ind=pp;
  _aForcings=NULL;
  _aSVtypes=NULL;
  _aLayers=NULL;
  _nCols=0;
  _pModel=pMod;
}
CCustomTable::~CCustomTable()
{
  delete [] _aForcings;
  delete [] _aSVtypes;
  delete [] _aLayers;
}
//////////////////////////////////////////////////////////////////
/// \brief adds column to table
//
void  CCustomTable::ExpandArrays()
{
  sv_type       *tmp=new sv_type      [_nCols+1];
  forcing_type *tmp2=new forcing_type [_nCols+1];
  int          *tmp3=new int          [_nCols+1];

  for (int i = 0; i < _nCols; i++) {
    tmp [i]=_aSVtypes[i];
    tmp2[i]=_aForcings[i];
    tmp3[i]=_aLayers[i];
  }
  tmp [_nCols]=UNRECOGNIZED_SVTYPE;
  tmp2[_nCols]=F_UNRECOGNIZED;
  tmp3[_nCols]=DOESNT_EXIST;
  if (_nCols>0){delete [] _aSVtypes;delete [] _aForcings;delete [] _aLayers; }
  _aSVtypes=tmp;
  _aForcings=tmp2;
  _aLayers=tmp3;
  _nCols++;
}
//////////////////////////////////////////////////////////////////
/// \brief adds state variable column to table
//
void  CCustomTable::AddStateVariable(const sv_type& sv, const int lay) {
  ExpandArrays();
  _aSVtypes[_nCols-1]=sv;
  _aLayers [_nCols-1]=lay;
}
//////////////////////////////////////////////////////////////////
/// \brief adds forcing column to table
//
void   CCustomTable::AddForcing(const forcing_type& ff)
{
  ExpandArrays();
  _aForcings[_nCols-1]=ff;
}
//////////////////////////////////////////////////////////////////
/// \brief Open a stream to the file and write header info
//
void  CCustomTable::WriteFileHeader(const optStruct& Options)
{
  _CUSTTAB.open(_filename.c_str());
  if (_CUSTTAB.fail()){
    WriteWarning("CCustomTable::WriteFileHeader: Unable to create file "+_filename,Options.noisy);
  }
  _CUSTTAB<<",Datetime";
  for (int i = 0; i < _nCols; i++) {
    if (_aSVtypes[i] != UNRECOGNIZED_SVTYPE) {
      _CUSTTAB<<","<<_pModel->GetStateVarInfo()->SVTypeToString(_aSVtypes[i],_aLayers[i]);
    }
    else {
      _CUSTTAB<<","<<ForcingToString(_aForcings[i]);
    }
  }
  _CUSTTAB<<endl;
}
//////////////////////////////////////////////////////////////////
/// \brief Write custom table to file by querying the model
/// \param &tt [in] Current model time
/// \param &Options [in] Global model options information
//
void  CCustomTable::WriteCustomTable(const time_struct& tt, const optStruct& Options)
{
  int pp=_sb_grp_ind;
  _CUSTTAB<<tt.model_time<<","<<tt.date_string; //TODO - fix Date String
  for (int i = 0; i < _nCols; i++) {
    if (_aSVtypes[i] != UNRECOGNIZED_SVTYPE) {
      _CUSTTAB<<","<<_pModel->GetSubBasinGroup(pp)->GetAvgStateVar(_pModel->GetStateVarIndex(_aSVtypes[i],_aLayers[i]));
    }
    else {

      _CUSTTAB<<","<<_pModel->GetSubBasinGroup(pp)->GetAvgForcing(_aForcings[i]);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Closes output stream after all information written to file
//
void  CCustomTable::CloseFile(const optStruct& Options)
{
  _CUSTTAB.close();
}
