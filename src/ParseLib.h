/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Teams
  ----------------------------------------------------------------*/
#ifndef PARSELIB_H
#define PARSELIB_H
#define _CRT_SECURE_NO_DEPRECATE 1
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <string>
#include <stdlib.h>
#include <string.h>

using namespace std;

//type definitions-------------------
typedef const double *               Unchangeable1DArray; ///< unchangeable but movable 1d array
typedef       double * const         Writeable1DArray;    ///< unmoveable but changeble 1d array
typedef const double * const         Ironclad1DArray;     ///< unmodifiable 1d array

typedef const double * const *       Unchangeable2DArray; ///< unchangeable but moveable 2d array
typedef       double *       * const Writeable2DArray;    ///< unmoveable but changeable 2d array
typedef const double * const * const Ironclad2DArray;     ///< unmodifiable 2d array

//Strings/Parser
//************************************************************************
const int    MAXINPUTITEMS  =    500;  ///< maximum delimited input items per line
const int    MAXCHARINLINE  =   6000;  ///< maximum characters in line
const bool   parserdebug=false;   ///< turn to true for debugging of parser

///////////////////////////////////////////////////////
/// \brief Types of parsing errors
//
enum parse_error
{
  PARSE_BAD,        ///< Corrupt file
  PARSE_NOT_ENOUGH, ///< Not enough parameters
  PARSE_TOO_MANY,   ///< Too many parameters
  PARSE_GOOD,       ///< No error
  PARSE_EOF         ///< End of file error
};

///////////////////////////////////////////////////////////////////
/// \brief Class for parsing data from file
//
class CParser
{
private:

  ifstream *INPUT;      ///current input file
  int      l;                     ///current line in input file
  string   filename;///current input filename

  bool     comma_only;//true if spaces & tabs ignored in tokenization

public:

  CParser(ifstream &FILE, const int init_line_num);
  CParser(ifstream &FILE, string filename, const int init_line_num);
  ~CParser(){}

  void   SetLineCounter(int i);
  int    GetLineNumber ();
  string GetFilename();
  void   ImproperFormat(char **s);
  void   IgnoreSpaces  (bool ignore){comma_only=ignore;}

  bool   Tokenize(char **tokens, int &numwords);

  void         SkipLine         ();

  /* [double]                   */
  parse_error    Parse_dbl                              (double &v1);
  /* [double] [double]          */
  parse_error    Parse_dbl                              (double &v1, double &v2);
  /* [double] [double] [double] */
  parse_error    Parse_dbl                              (double &v1, double &v2, double &v3);
  /* [double] [double] [double] [double] */
  parse_error    Parse_dbl                              (double &v1, double &v2, double &v3, double &v4);
  /* [int   ] [double] [double] */
  parse_error    Parse_intdbldbl        (   int &v1, double &v2, double &v3);
  /* [int   ]                   */
  parse_error    Parse_int                              (   int &v1);
  /* [int   ] [int  ]                  */
  parse_error    Parse_int                              (   int &v1,    int &v2);

  /* [double]
     ...
     [double]  (fixed (known) array size)
     & (if optfollow=true)*/
  parse_error    ParseArray_dbl   (Writeable1DArray v,  int numv, int &optfollow);

  /* [double] [double]
     ...      ...
     [double] [double] (fixed (known) array size)
     & (if optfollow=true)*/
  parse_error    ParseArray_dbl   (Writeable1DArray v1,
                                   Writeable1DArray v2, int numv, int &optfollow);
  /* [double] [double] [double]
     ...      ...
     [double] [double] [double] (fixed (known) array size)
     & (if optfollow=true) */
  parse_error    ParseArray_dbl   (Writeable1DArray v1,
                                   Writeable1DArray v2,
                                   Writeable1DArray v3, int numv, int &optfollow);
  /* [double]
     ...
     [double]  (dynamic (unknown) array size, max entries maxv)
     & (if optfollow=true) */
  parse_error    ParseArray_dbl_dyn(Writeable1DArray v, int &numv, const int maxv, int &optfollow);

  /* [double]  [double]
     ...
     [double]  [double] (dynamic (unknown) array sizes, max entries maxv)
     & (if optfollow=true) */
  parse_error    ParseArray_dbldbl_dyn(Writeable1DArray v1, Writeable1DArray v2, int &numv, const int maxv, int &optfollow);

  //special routine for borehole format
  parse_error  Parse2DArray_dbl (Writeable1DArray v1,
                                 Writeable1DArray v2,
                                 Writeable2DArray v3, int numv, int numcol, int &optfollow);
  /* [double] [double] ... [double]
     ...       ...   ...    ...
     [double] [double] ... [double]              (fixed (known) 2D array size)
     & (if optfollow=true) */
  parse_error  Parse2DArray_dbl  (Writeable2DArray v3, int numv, int numcol, int &optfollow);

  /* [double] [double] ... [double] ... [double]
     ...       ...   ...    ...   ...   ...
     [double] [double] ... [double]              (fixed (known) array size)
     & (if optfollow=true) */
  parse_error    ParseBigArray_dbl (Writeable1DArray v,  int numv);
};

#endif
