/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------*/

#include "ParseLib.h"
#include "RavenInclude.h"
inline int      s_to_i (char *s1)            {return (int)atof(s1);   }
inline double   s_to_d (char *s1)            {return atof(s1);        }
inline bool     s_to_b (char *s1)            {return ((int)atof(s1)!=0);   }

/*----------------------------------------------------------------
  Constructor
  -----------------------------------------------------------------------*/
CParser::CParser(ifstream &FILE, const int i)
{
  filename="";
  INPUT =&FILE;
  l=i;
  comma_only=false;
}
//-----------------------------------------------------------------------
CParser::CParser(ifstream &FILE, string _filename, const int i)
{
  filename=_filename;
  INPUT =&FILE;
  l=i;
  comma_only=false;
}
/*----------------------------------------------------------------
  Basic Member Functions
  -----------------------------------------------------------------------*/
void   CParser::SetLineCounter(int i)    {l=i;}
//-----------------------------------------------------------------------
int    CParser::GetLineNumber ()         {return l;}
//-----------------------------------------------------------------------
string CParser::GetFilename() {return filename;}
/*----------------------------------------------------------------
  Tokenize
  ----------------------------------------------------------------
  tokenizes a sentence delimited by  spaces, tabs & return characters

  parameters:
  out is the array of strings in the line
  numwords is the number of strings in the line
  returns true if file has ended
  -------------------------------------------------------------------------*/
bool CParser::Tokenize(char **out, int &numwords){

  static char wholeline     [MAXCHARINLINE];
  static char *tempwordarray[MAXINPUTITEMS];
  char *p;
  //char *junk=NULL;
  int ct(0),w;
  char delimiters[6];
  delimiters[0]=' '; //space
  delimiters[1]='\t';//tab
  delimiters[2]=','; //comma
  delimiters[3]='\r';//carriage return (ignored in visual c++)
  delimiters[4]='\n';//newline -necc?
  delimiters[5]='\0';//last delimiter
  if (comma_only){
    delimiters[0]=delimiters[1]=',';
  }
  (*wholeline)=0;
  while ((*wholeline)==0){//while loop handles blank lines
    if (INPUT->eof()){return true;}
    INPUT->getline(wholeline,MAXCHARINLINE);            //get entire line as 1 string
    if (INPUT->fail()){
      //cout<<"failed: "<<filename<<" line "<<l<<"|"<<wholeline<<"|"<<INPUT->ios::eofbit<<endl;
      //ExitGracefully("Too many characters in line or (maybe) using Mac-style carriage return line endings.",BAD_DATA);
    }
    l++;
    if ((parserdebug) && ((*wholeline)!=0)){cout <<wholeline<<endl;}
  }                 //if line is null, break, return true
  p=strtok(wholeline, delimiters);
  while (p){                                         //sift through words, place in temparray, count line length
    tempwordarray[ct]=p;
    //p=strtok_s(NULL, delimiters,&junk);
    p=strtok(NULL, delimiters);
    ct++;
    if ((p!=NULL) && (p[0]=='#')){break;} //ignore all content after '#'
  }
  for (w=0; w<ct; w++){                              //copy temp array of words into out[]
    if (w>MAXINPUTITEMS){numwords=ct;
      string warn="Tokenizeline:: exceeded maximum number of items in single line in file "+filename;
      ExitGracefully(warn.c_str(),BAD_DATA);
      return true;
    }
    out[w]=tempwordarray[w];
    //cout<<out[w]<<"|";
  }
  //cout <<endl;
  numwords=ct;
  return false;
}
/*----------------------------------------------------------------*/
void   CParser::ImproperFormat(char **s)
{
  cout <<"line "<< l << " in file "<<filename<< " is wrong length"<<endl; //FOR NOW
  cout <<"line "<< l << ": "<<s[0]<<endl;
}
/*----------------------------------------------------------------*/
void CParser::SkipLine()
{
  int      Len;
  char    *s[MAXINPUTITEMS];
  if (Tokenize(s,Len)){}
}
/*----------------------------------------------------------------
  Parse_dbl
  ----------------------------------------------------------------
  Parses a single line from an input file expected to have 1,2, or 3 DOUBLE
  input parameters
  ----------------------------------------------------------------*/
parse_error CParser::Parse_dbl(double &v1)
{
  int      Len;
  char    *s[MAXINPUTITEMS];

  if (Tokenize(s,Len))                     {return PARSE_EOF;   }

  if      (Len==1) {v1=s_to_d(s[0]);                              return PARSE_GOOD;}
  else             {ImproperFormat(s);            return PARSE_BAD;       }
}
//------------------------------------------------------------------------------
parse_error CParser::Parse_dbl(double &v1, double &v2)
{
  int      Len;
  char    *s[MAXINPUTITEMS];

  if (Tokenize(s,Len))                     {return PARSE_EOF;   }

  if      (Len==2) {v1=s_to_d(s[0]);
    v2=s_to_d(s[1]);                                return PARSE_GOOD;}
  else             {ImproperFormat(s);            return PARSE_BAD;       }
}
//------------------------------------------------------------------------------
parse_error CParser::Parse_dbl(double &v1, double &v2, double &v3)
{
  int      Len;
  char    *s[MAXINPUTITEMS];

  if (Tokenize(s,Len))                     {return PARSE_EOF;   }

  if      (Len==3) {v1=s_to_d(s[0]);
    v2=s_to_d(s[1]);
    v3=s_to_d(s[2]);                                return PARSE_GOOD;}
  else             {ImproperFormat(s);            return PARSE_BAD;       }

}
//------------------------------------------------------------------------------
parse_error CParser::Parse_dbl(double &v1, double &v2, double &v3, double &v4)
{
  int      Len;
  char    *s[MAXINPUTITEMS];

  if (Tokenize(s,Len))                     {return PARSE_EOF;   }

  if      (Len==4) {v1=s_to_d(s[0]);
    v2=s_to_d(s[1]);
    v3=s_to_d(s[2]);
    v4=s_to_d(s[3]);        return PARSE_GOOD;}
  else             {ImproperFormat(s);            return PARSE_BAD;       }

}
//------------------------------------------------------------------------------
parse_error CParser::Parse_intdbldbl(int &v1, double &v2, double &v3)
{
  int      Len;
  char    *s[MAXINPUTITEMS];

  if (Tokenize(s,Len))                     {return PARSE_EOF;   }

  if      (Len==3) {v1=s_to_i(s[0]);
    v2=s_to_d(s[1]);
    v3=s_to_d(s[2]);                                return PARSE_GOOD;}
  else             {ImproperFormat(s);            return PARSE_BAD;       }

}
//------------------------------------------------------------------------------
parse_error CParser::Parse_int(int &v1)
{
  int      Len;
  char    *s[MAXINPUTITEMS];

  if (Tokenize(s,Len))                     {return PARSE_EOF;   }

  if   (Len==1){v1=s_to_i(s[0]);                          return PARSE_GOOD;}
  else         {ImproperFormat(s);                return PARSE_BAD;       }
}
//------------------------------------------------------------------------------
parse_error CParser::Parse_int(int &v1, int &v2)
{
  int      Len;
  char    *s[MAXINPUTITEMS];

  if (Tokenize(s,Len))                     {return PARSE_EOF;   }

  if      (Len==2) {v1=s_to_i(s[0]);
    v2=s_to_i(s[1]);                                return PARSE_GOOD;}
  else             {ImproperFormat(s);            return PARSE_BAD;       }
}
/*-------------------------------------------------------------------------
  ParseArray_dbl
  -------------------------------------------------------------------------
  Parses a set of numv lines from an input file, each expected to have 1,2, or 3 DOUBLE
  input parameters
  Ended with "&" [optint]
  Ex.:
  0.1 2.3 4.5
  0.1 2.3 4.5
  0.1 2.3 4.5
  0.1 2.3 4.5
  & 5
  -------------------------------------------------------------------------*/
parse_error      CParser::ParseArray_dbl   (Writeable1DArray v,  int numv, int &optfollow)
{
  bool             done(false);
  int                      count(0),Len;
  char    *s[MAXINPUTITEMS];

  do
  {
    if (Tokenize(s,Len)){return PARSE_EOF;  }

    if ((Len==2) && (!strcmp(s[0],"&")))
    {
      optfollow=s_to_i(s[1]);
      done=true;
    }
    else if  ((Len==1) && (strcmp(s[0],"&")))
    {
      if (count>=numv){return PARSE_TOO_MANY;}
      v[count]=s_to_d(s[0]);
      count++;
    }
    else if ((Len==1) && (!strcmp(s[0],"&")))
    {
      done=true;
    }
    else if  (Len!=0)
    {
      ImproperFormat(s); return PARSE_BAD;
    }
  } while (!done);

  if (count!=numv){return PARSE_NOT_ENOUGH;}

  return PARSE_GOOD;
}
//------------------------------------------------------------------------------
parse_error      CParser::ParseArray_dbl   (Writeable1DArray v1,
                                            Writeable1DArray v2, int numv, int &optfollow)
{
  bool             done(false);
  int                      count(0),Len;
  char    *s[MAXINPUTITEMS];

  do
  {
    if (Tokenize(s,Len)){return PARSE_EOF;  }

    if ((Len==2) && (!strcmp(s[0],"&")))
    {
      optfollow=s_to_i(s[1]);
      done=true;
    }
    else if  (Len==2)
    {
      if (count>=numv){return PARSE_TOO_MANY;}
      v1[count]=s_to_d(s[0]);
      v2[count]=s_to_d(s[1]);
      count++;
    }
    else if ((Len==1) && (!strcmp(s[0],"&")))
    {
      done=true;
    }

    else if  (Len!=0)
    {
      ImproperFormat(s); return PARSE_BAD;
    }
  } while (!done);

  if (count!=numv){return PARSE_NOT_ENOUGH;}

  return PARSE_GOOD;
}
//------------------------------------------------------------------------------
parse_error      CParser::ParseArray_dbl (Writeable1DArray v1,
                                          Writeable1DArray v2,
                                          Writeable1DArray v3, int numv, int &optfollow){
  bool             done(false);
  int                      count(0),Len;
  char    *s[MAXINPUTITEMS];

  do
  {
    if (Tokenize(s,Len)){return PARSE_EOF;  }

    if  (Len==3)
    {
      if (count>=numv){return PARSE_TOO_MANY;}
      v1[count]=s_to_d(s[0]);
      v2[count]=s_to_d(s[1]);
      v3[count]=s_to_d(s[2]);
      count++;
    }
    else if ((Len==1) && (!strcmp(s[0],"&")))
    {
      done=true;
    }
    else if ((Len==2) && (!strcmp(s[0],"&")))
    {
      optfollow=s_to_i(s[1]);
      done=true;
    }
    else if  (Len!=0)
    {
      ImproperFormat(s); return PARSE_BAD;
    }
  } while (!done);

  if (count!=numv){return PARSE_NOT_ENOUGH;}

  return PARSE_GOOD;
}

/*-------------------------------------------------------------------------
  ParseBigArray_dbl
  -------------------------------------------------------------------------
  Parses an unknown number lines (max length=20) from an input file until numv items have been read
  list ends with "&
  Ended with "&" [optint]
  Ex.:
  0.1 2.3 4.5
  0.1 2.3 4.5 6.7
  0.1 2.3 4.5 6.7 7.8
  0.1
  &
  -------------------------------------------------------------------------*/
parse_error CParser::ParseBigArray_dbl(Writeable1DArray v, int numv)
{
  bool             done(false);
  int                      count(0),Len,i;
  char    *s[MAXINPUTITEMS];

  do
  {
    if (Tokenize(s,Len)){return PARSE_EOF;  }

    if  ((Len<=20) && (Len>=1) && (strcmp(s[0],"&")))
    {
      for (i=0;i<Len;i++)
      {
        if (count>=numv){return PARSE_TOO_MANY;}
        v[count]=s_to_d(s[i]);
        count++;
      }
    }
    else if ((Len==1) && (!strcmp(s[0],"&")))
    {
      if (count!=numv){return PARSE_NOT_ENOUGH;}
      done=true;
    }
    else if  (Len!=0)
    {
      ImproperFormat(s); return PARSE_BAD;
    }

  } while (!done);

  return PARSE_GOOD;
}
/*-------------------------------------------------------------------------
  Parse2DArray_dbl
  -------------------------------------------------------------------------*
  Parses an known number lines (length=2+numcol) from an input file until numv items have been read
  list ends with "&
  Ended with "&" [optint]
  Ex.:
  0.1 2.3 4.5 4.5 4.5 4.5
  0.1 2.3 4.5 4.5 4.5 4.5
  0.1 2.3 4.5 4.5 4.5 4.5
  & 6
  (for numv=3, numcol=4)
  -------------------------------------------------------------------------*/
parse_error  CParser::Parse2DArray_dbl  (Writeable1DArray v1,
                                         Writeable1DArray v2,
                                         Writeable2DArray v3, int numv, int numcol, int &optfollow){
  bool             done(false);
  int                      count(0),Len;
  char    *s[MAXINPUTITEMS];
  do
  {
    if (Tokenize(s,Len))  {return PARSE_EOF;};
    if (IsComment(s[0],Len)){}
    else if  (Len==(2+numcol)){
      if (count>=numv){return PARSE_TOO_MANY;}
      v1[count]=s_to_d(s[0]);
      v2[count]=s_to_d(s[1]);
      for (int j=0;j<numcol;j++){
        v3[count][j]=s_to_d(s[2+j]);
      }
      count++;
    }
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
      if (count!=numv){return PARSE_NOT_ENOUGH;}
      done=true;
    }
    else if ((Len==2) && (!strcmp(s[0],"&"))) {
      if (count!=numv){return PARSE_NOT_ENOUGH;}
      optfollow=s_to_i(s[1]);
      done=true;
    }
    else if  (Len!=0){
      return PARSE_BAD;
    }
  } while (!done);

  return PARSE_GOOD;
}
/*-------------------------------------------------------------------------
  Parse2DArray_dbl
  -------------------------------------------------------------------------
  Parses an known number lines (length=numcol) from an input file until numv items have been read
  list ends with "&
  Ended with "&" [optint]
  Ex.:
  4.5 4.5 4.5 4.5
  4.5 4.5 4.5 4.5
  4.5 4.5 4.5 4.5
  & 6
  (for numv=3, numcol=4)
  -------------------------------------------------------------------------*/
parse_error  CParser::Parse2DArray_dbl  (Writeable2DArray v3, int numv, int numcol, int &optfollow){
  bool             done(false);
  int                      count(0),Len;
  char    *s[MAXINPUTITEMS];
  do
  {
    if (Tokenize(s,Len))  {return PARSE_EOF;};
    if  (Len==(numcol)){
      if (count>=numv){return PARSE_TOO_MANY;}
      for (int j=0;j<numcol;j++){
        v3[count][j]=s_to_d(s[j]);
      }
      count++;
    }
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
      if (count!=numv){return PARSE_NOT_ENOUGH;}
      done=true;
    }
    else if ((Len==2) && (!strcmp(s[0],"&"))) {
      if (count!=numv){return PARSE_NOT_ENOUGH;}
      optfollow=s_to_i(s[1]);
      done=true;
    }
    else if  (Len!=0){
      return PARSE_BAD;
    }
  } while (!done);

  return PARSE_GOOD;
}
parse_error     CParser::ParseArray_dbl_dyn(Writeable1DArray v, int &numv, const int maxv, int &optfollow)
{
  bool             done(false);
  int                      count(0),Len;
  char    *s[MAXINPUTITEMS];

  do
  {
    if (Tokenize(s,Len)){return PARSE_EOF;};
    if  ((Len==1) && (strcmp(s[0],"&")))
    {
      if (count>=maxv)    {return PARSE_TOO_MANY;}
      v[count]=s_to_d(s[0]);
      count++;
    }
    else if ((Len<=2) && (!strcmp(s[0],"&")))
    {
      if (Len==2){
        optfollow= s_to_i(s[1]);
      }
      done=true;
    }
    else if (Len!=0)
    {
      ImproperFormat(s);
      return PARSE_BAD;
    }
  } while (!done);

  numv=count;
  return PARSE_GOOD;
}
parse_error     CParser::ParseArray_dbldbl_dyn(Writeable1DArray v1,
                                               Writeable1DArray v2,
                                               int &numv,
                                               const int maxv,
                                               int &optfollow)
{
  bool             done(false);
  int                      count(0),Len;
  char    *s[MAXINPUTITEMS];

//      cout <<"ParseArray_dbldbl_dyn"<<endl;
  do
  {
    if (Tokenize(s,Len)){return PARSE_EOF;};
    if  ((Len==2) && (strcmp(s[0],"&")))
    {
//                      cout <<"len=2!"<<endl;
      if (count>=maxv)    {return PARSE_TOO_MANY;}
      v1[count]=s_to_d(s[0]);
      v2[count]=s_to_d(s[1]);
      count++;
    }
    else if ((Len<=2) && (!strcmp(s[0],"&")))
    {
      if (Len==2){
        optfollow= s_to_i(s[1]);
      }
      done=true;
    }
    else if (Len!=0)
    {
      ImproperFormat(s);
      return PARSE_BAD;
    }
  } while (!done);

  numv=count;
  return PARSE_GOOD;
}
