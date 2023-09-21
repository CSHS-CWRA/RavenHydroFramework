//SparseMatrix.cpp
#include "SparseMatrix.h"
#include "RavenInclude.h"

#ifdef _ARMADILLO_
#include<armadillo>
using namespace arma;


//-----------------------------------------------------------------------
CSparseMatrix::CSparseMatrix(){
	sA      =NULL;
	ija     =NULL;
	nNonZero=0;
	size    =0;
	buffersize=0;
	dirichlet=NULL;
	p =NULL;
	pp=NULL;
	r =NULL;
	rr=NULL;
	z =NULL;
	zz=NULL;

	// Intialize armadillo library vectors;     R Simms Mar 21 2011
#ifdef _ARMADILLO_
	pArma = new arma::vec(1);
	ppArma = new arma::vec(1);
	rArma = new arma::vec(1);
	rrArma = new arma::vec(1);
	zArma = new arma::vec(1);
	zzArma = new arma::vec(1);
#endif
}
//-----------------------------------------------------------------------
CSparseMatrix::CSparseMatrix(double *spA, int *indices, double *bvec){
	sA      =spA;
	ija     =indices;
  b       =bvec;
	//nNonZero=0;
	//buffersize=0;
	dirichlet=NULL;
	p =NULL;
	pp=NULL;
	r =NULL;
	rr=NULL;
	z =NULL;
	zz=NULL;

	// Intialize armadillo library vectors;     R Simms Mar 21 2011
#ifdef _ARMADILLO_
	pArma = new arma::vec(1);
	ppArma = new arma::vec(1);
	rArma = new arma::vec(1);
	rrArma = new arma::vec(1);
	zArma = new arma::vec(1);
	zzArma = new arma::vec(1);
#endif
}
//-----------------------------------------------------------------------
CSparseMatrix::CSparseMatrix(double *spA, int *indices, double *bvec, int sizen){
	sA      =spA;
	ija     =indices;
  b       =bvec;
	//nNonZero=0;
  size = sizen;
	//buffersize=0;

  dirichlet=new bool   [size];
	p        =new double [size];
	pp       =new double [size];
	r        =new double [size];
	rr       =new double [size];
	z        =new double [size];
	zz       =new double [size];

	// Intialize armadillo library vectors;
#ifdef _ARMADILLO_
	pArma = new arma::vec(size);
	ppArma = new arma::vec(size);
	rArma = new arma::vec(size);
	rrArma = new arma::vec(size);
	zArma = new arma::vec(size);
	zzArma = new arma::vec(size);
#endif
  for (int j=0; j<size;j++){
		p[j]=pp[j]=r[j]=rr[j]=z[j]=zz[j]=0.0;
		dirichlet[j]=false;
	}
}
//-----------------------------------------------------------------------
CSparseMatrix::CSparseMatrix(const CSparseMatrix &A){

	int j,k;
	this->size    =A.size;
	this->nNonZero=A.nNonZero;
	this->buffersize=this->nNonZero+1;
	//Allocate memory
	this->sA =new double [this->nNonZero+1];
	this->ija=new int    [this->nNonZero+1];

	for (k=0; k<=nNonZero;k++){
		this->sA [k]=A.sA [k];
		this->ija[k]=A.ija[k];
	}
	p        =new double [size];
	pp       =new double [size];
	r        =new double [size];
	rr       =new double [size];
	z        =new double [size];
	zz       =new double [size];
	dirichlet=new bool   [size];

	// Intialize armadillo library vectors;     R Simms Mar 21 2011
#ifdef _ARMADILLO_
	pArma = new arma::vec(size);
	ppArma = new arma::vec(size);
	rArma = new arma::vec(size);
	rrArma = new arma::vec(size);
	zArma = new arma::vec(size);
	zzArma = new arma::vec(size);
#endif
	for (j=0; j<size;j++){
		p[j]=pp[j]=r[j]=rr[j]=z[j]=zz[j]=0.0;
		dirichlet[j]=false;
	}
}
//-----------------------------------------------------------------------
CSparseMatrix::~CSparseMatrix(){
	//if (globaldebug){cout<<"DELETING SPARSE MATRIX"<<endl;}
  cout<<"DELETING SPARSE MATRIX"<<endl;
	delete [] sA;
	delete [] ija;
	nNonZero=0;
	size    =0;
	buffersize=0;
	delete [] p;
	delete [] pp;
	delete [] r;
	delete [] rr;
	delete [] z;
	delete [] zz;
	delete [] dirichlet;

	delete pArma;
	delete ppArma;
	delete rArma;
	delete rrArma;
	delete zArma;
	delete zzArma;
}
/************************************************************************
	Translates A[N][N] to sparse matrix S
------------------------------------------------------------------------*/
CSparseMatrix::CSparseMatrix(Ironclad2DArray A, const int N, const double thresh){
	int i,j,k;
	size    =N;
	nNonZero=N;
	//count non-zero elements
	for (i=0; i<size; i++){
		for (j=0;j<size; j++){
			if ((fabs(A[i][j])>=thresh) && (i!=j)){nNonZero++;}
		}
	}
	//Allocate memory
	sA =new double [nNonZero+1];
	ija=new int    [nNonZero+1];
	if (ija==NULL){MatExitGracefully("CSparseMatrix::Constructor::Out of memory",OUT_OF_MEMORY_MAT);}
	buffersize=nNonZero+1;

	ija[0]=size+1;//index to row 0 of off-diagonal element
	for (j=0;j<size; j++){sA[j]=A[j][j];}
	k=size;
	for (i=0; i<size; i++){
		for (j=0;j<size; j++){
			if ((fabs(A[i][j])>thresh) && (i!=j)){
				if (k+1>nNonZero){MatExitGracefully("sprsTranslate: bad index calc",RUNTIME_ERR_MAT);}
				k++;
				//++k;
				sA [k]=A[i][j];
				ija[k]=j;
			}
		}
		if (i>nNonZero){MatExitGracefully("sprsTranslate: bad index calc",RUNTIME_ERR_MAT);}
		ija[i+1]=k+1;
	}
	p       =new double [size];
	pp      =new double [size];
	r       =new double [size];
	rr      =new double [size];
	z       =new double [size];
	zz      =new double [size];
	dirichlet=new bool   [size];

	// Intialize armadillo library vectors;     R Simms Mar 21 2011
#ifdef _ARMADILLO_
	pArma = new arma::vec(size);
	ppArma = new arma::vec(size);
	rArma = new arma::vec(size);
	rrArma = new arma::vec(size);
	zArma = new arma::vec(size);
	zzArma = new arma::vec(size);
#endif
	for (j=0; j<size;j++){
		p[j]=pp[j]=r[j]=rr[j]=z[j]=zz[j]=0.0;
		dirichlet[j]=false;
	}
}
//-----------------------------------------------------------------------
int CSparseMatrix::GetNumEntries() const{return nNonZero;}
/************************************************************************
 DynamicInitialize:
	Initializes dynamic sparse array with zero diagonal values
------------------------------------------------------------------------*/
void CSparseMatrix::DynamicInitialize(const int N){
	if (size!=0) {
		MatExitGracefully("CSparseMatrix::DynamicInitialize::Already initialized",RUNTIME_ERR_MAT);}
	size    =N;
	nNonZero=N;
	//Allocate memory
	sA =new double [size+1];
	ija=new int    [size+1];
	if (ija==NULL){
		MatExitGracefully("CSparseMatrix::DynamicInitialize::Out of memory",OUT_OF_MEMORY_MAT);}
	buffersize=size+1;
	ija[0]=size+1;
	for (int i=0;i<size; i++){
		sA[i]   =0.0;
		ija[i+1]=size+1;
	}
	sA[size]=0.0; //Arbitrary value (doesn't matter)
	p       =new double [size];
	pp      =new double [size];
	r       =new double [size];
	rr      =new double [size];
	z       =new double [size];
	zz      =new double [size];
	dirichlet=new bool   [size];


	// Intialize armadillo library vectors;     R Simms Mar 21 2011
#ifdef _ARMADILLO_
	pArma = new arma::vec(size);
	ppArma = new arma::vec(size);
	rArma = new arma::vec(size);
	rrArma = new arma::vec(size);
	zArma = new arma::vec(size);
	zzArma = new arma::vec(size);
#endif
	for (int j=0; j<size;j++){
		p[j]=pp[j]=r[j]=rr[j]=z[j]=zz[j]=0.0;
		dirichlet[j]=false;
	}
}
/************************************************************************
 DynamicInitialize:
	Initializes dynamic sparse array with diagonal values only
------------------------------------------------------------------------*/
void CSparseMatrix::DynamicInitialize(Ironclad1DArray Adiag,
																      const int       N){
	if (size!=0) {MatExitGracefully("CSparseMatrix::DynamicInitialize::Already initialized",RUNTIME_ERR_MAT);}
	size    =N;
	nNonZero=N;
	//Allocate memory
	sA =new double [size+1];
	ija=new int    [size+1];
	if (ija==NULL){MatExitGracefully("CSparseMatrix::DynamicInitialize::Out of memory",OUT_OF_MEMORY_MAT);}
	buffersize=size+1;

	ija[0]=size+1;
	for (int i=0;i<size; i++){
		sA [i]  =Adiag[i];
		ija[i+1]=size+1;
	}
	sA[size]=0.00001; //Arbitrary value (doesn't matter)
	p =new double [size];
	pp=new double [size];
	r =new double [size];
	rr=new double [size];
	z =new double [size];
	zz=new double [size];
  dirichlet=new bool [size];

  // Intialize armadillo library vectors;     R Simms Mar 21 2011
#ifdef _ARMADILLLO_
	pArma = new arma::vec(size);
	ppArma = new arma::vec(size);
	rArma = new arma::vec(size);
	rrArma = new arma::vec(size);
	zArma = new arma::vec(size);
	zzArma = new arma::vec(size);
#endif

	for (int j=0; j<size;j++){
		p[j]=pp[j]=r[j]=rr[j]=z[j]=zz[j]=0.0;
		dirichlet[j]=false;
	}
}

/************************************************************************
 DynamicAdd:
	Dynamically adds a value to row=i, col=j, if that entry does not exist it gets created
	could get costly with large matrices;
	but is helped with buffering of the matrices-rebuffering is called 4 times for each bandwidth increase

	Created by R.Simms, Mar 15, 2011
------------------------------------------------------------------------*/
void CSparseMatrix::DynamicAdd(const double &a, const int i, const int j){

	int buffer=(size/4);

	int    *tmpija;
	double *tmpsA;
	int     knew, i1,k;

	if (size==0)                  {MatExitGracefully("CSparseMatrix::DynamicAdd::not initialized correctly",RUNTIME_ERR_MAT);}
	if ((sA==NULL) || (ija==NULL)){MatExitGracefully("CSparseMatrix::DynamicAdd: not initialized correctly",RUNTIME_ERR_MAT);}
	if ((i>=size)  || (j>=size))  {MatExitGracefully("CSparseMatrix::DynamicAdd: bad indices specified",RUNTIME_ERR_MAT);}
	if ((i<0)      || (j<0))      {MatExitGracefully("CSparseMatrix::DynamicAdd: bad indices specified(2)",RUNTIME_ERR_MAT);}

	if (i==j){sA[i]+=a; return;}

	if (a==0.0){return;}                     //do not dynamically add if a zero entry (should this be here?}

	//find k index of inserted coefficient--------------------------------
	knew=ija[i+1];                           //default- guess last coefficient in row i
	for (k=ija[i]; k<ija[i+1]; k++){           //searches through row i
		if      (ija[k]==j){sA[k]+=a;return;}     //if already exists, add a to the value located at this location
		else if (ija[k]> j){knew=k; break;}      //otherwise, adopt knew
	}

	//copy all terms in array and insert----------------------------------
	nNonZero++;                              //item added to total array

	if (knew>nNonZero){cout << "bad knew"<<endl;}

	if (buffersize<=nNonZero){
		tmpija=new int    [nNonZero+1+buffer];
		tmpsA =new double [nNonZero+1+buffer];
		if (tmpsA==NULL){
			MatExitGracefully("CSparseMatrix::DynamicAdd::Out of memory",OUT_OF_MEMORY_MAT);}
		buffersize=nNonZero+1+buffer;

		for (k=0; k<knew; k++){               //before inserted coefficient
			tmpija[k]=ija[k];
			tmpsA [k]= sA[k];
		}
		for (k=nNonZero; k>=knew+1; k--){     //after inserted coefficient
			tmpija[k]=ija[k-1];
			tmpsA [k]= sA[k-1];
		}
		tmpija[knew]=j;                       //inserted coefficient
		tmpsA [knew]=a;

		delete [] ija;                        //delete old values
		delete [] sA;

		ija=tmpija;                           //redirect array pointers
		sA =tmpsA;
	}
	else{  //no need to update buffer

		for (k=nNonZero; k>=knew+1; k--){     //after inserted coefficient //still pretty slow!
			ija[k]=ija[k-1];
			sA [k]= sA[k-1];
		}
		ija[knew]=j;                          //inserted coefficient
		sA [knew]=a;
	}

	for (i1=i; i1<size; i1++){ //update indexing information
		ija[i1+1]++;             //item added to row i (ija[i+1] is index of lowest term in row i+1)
	}

	return;
}

/************************************************************************
 DynamicInsert:
	Dynamically inserts new coefficient row=i, col=j, into sparse matrix
	could get costly with large matrices;
	but is helped with buffering of the matrices-rebuffering is called 4 times for each bandwidth increase

	created by J.Craig
	renamed for accuracy by R.Simms, Mar 15 2011
------------------------------------------------------------------------*/
void CSparseMatrix::DynamicInsert(const double &a, const int i, const int j){

	int buffer=(size/4);

	int    *tmpija;
	double *tmpsA;
	int     knew, i1,k;

	if (size==0)                  {MatExitGracefully("CSparseMatrix::DynamicAdd::not initialized correctly",RUNTIME_ERR_MAT);}
	if ((sA==NULL) || (ija==NULL)){MatExitGracefully("CSparseMatrix::DynamicAdd: not initialized correctly",RUNTIME_ERR_MAT);}
	if ((i>=size)  || (j>=size))  {MatExitGracefully("CSparseMatrix::DynamicAdd: bad indices specified",RUNTIME_ERR_MAT);}
	if ((i<0)      || (j<0))      {MatExitGracefully("CSparseMatrix::DynamicAdd: bad indices specified(2)",RUNTIME_ERR_MAT);}

	if (i==j){sA[i]=a; return;}

	//find k index of inserted coefficient--------------------------------
	knew=ija[i+1];                           //default- guess last coefficient in row i
	for (k=ija[i]; k<ija[i+1]; k++){           //searches through row i
		if      (ija[k]==j){sA[k]=a;return;}     //if already exists, merely replace coefficient (same as SetAij)
		else if (ija[k]> j){knew=k; break;}      //otherwise, adopt knew
	}
	if (a==0.0){return;}                     //do not dynamically add if a zero entry (should this be here?}

	//copy all terms in array and insert----------------------------------
	nNonZero++;                              //item added to total array

	if (knew>nNonZero){cout << "bad knew"<<endl;}

	if (buffersize<=nNonZero){
		tmpija=new int    [nNonZero+1+buffer];
		tmpsA =new double [nNonZero+1+buffer];
		if (tmpsA==NULL){
			MatExitGracefully("CSparseMatrix::DynamicAdd::Out of memory",OUT_OF_MEMORY_MAT);}
		buffersize=nNonZero+1+buffer;

		for (k=0; k<knew; k++){               //before inserted coefficient
			tmpija[k]=ija[k];
			tmpsA [k]= sA[k];
		}
		for (k=nNonZero; k>=knew+1; k--){     //after inserted coefficient
			tmpija[k]=ija[k-1];
			tmpsA [k]= sA[k-1];
		}
		tmpija[knew]=j;                       //inserted coefficient
		tmpsA [knew]=a;

		delete [] ija;                        //delete old values
		delete [] sA;

		ija=tmpija;                           //redirect array pointers
		sA =tmpsA;
	}
	else{  //no need to update buffer

		for (k=nNonZero; k>=knew+1; k--){     //after inserted coefficient //still pretty slow!
			ija[k]=ija[k-1];
			sA [k]= sA[k-1];
		}
		ija[knew]=j;                          //inserted coefficient
		sA [knew]=a;
	}

	for (i1=i; i1<size; i1++){ //update indexing information
		ija[i1+1]++;             //item added to row i (ija[i+1] is index of lowest term in row i+1)
	}

	return;
}
/************************************************************************
 MakeDirichlet:
	marks all rows and colums specified by rows[] as dirichlet
	affects BCG, MatVectMultiply, SolveEst, Print and CalculateAdjustedNorm
------------------------------------------------------------------------*/
void CSparseMatrix::MakeDirichlet(const int *rows, const int ndirichlet){
	for (int k=0; k<size; k++){
		dirichlet[k]=false;
	}
	for (int j=0; j<ndirichlet; j++){
		dirichlet[rows[j]]=true;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Sets rows to zero except for the diagonal which is set to 1
///
/// 			There is a strong potential for this method to be unstable because
/// 			the diagonal isnt normalized.</summary>
///
/// <remarks>	rsimms, 3/16/2011. </remarks>
///
/// <param name="rows">		 	The index numbers of the dirichlet nodes. </param>
/// <param name="ndirichlet">	The nnumber of dirichlet nodes. </param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CSparseMatrix::MakeDirichletWorkAround(const int *rows, const int ndirichlet){
	for(int i = 0; i < ndirichlet; i++){
		int row = rows[i];
		// set the off diagonls
		for(int j =0; j < size; j++){
			SetAij(0.0,row,j);
		}
		SetAij(1,row,row); //set the diagonals
	}
}
/************************************************************************
 MakeEssential1:
	revises matrix to meet essential conditions X_{j1}+a*X_{k1}=b
	affects matrix only- if condition is removed, matrix must be remade!!!!
	current implementation is pretty slow- CAN OPTIMIZE!!
	revised from Akin, pg. 403
------------------------------------------------------------------------*/
void CSparseMatrix::MakeEssential1(const int j1, const int k1, const double a, const double b){
	int i,j;
	if ((j1<0)      || (k1<0))      {
		MatExitGracefully("CSparseMatrix::MakeEssential1: bad indices specified",RUNTIME_ERR_MAT);}
	if ((j1==k1) && (a!=0.0)) {return;} //cannot force constraint on itself

	double Sjj=this->GetAij(j1,j1);
	double *Skrow=new double [this->size];
	double *Sjrow=new double [this->size];
	double *Sij=new double [this->size];
	double *Sik=new double [this->size];
	for (i=0; i<size; i++){
		Skrow[i]=this->GetAij(k1,i);
		Sjrow[i]=this->GetAij(j1,i);
    Sij[i]=this->GetAij(i,j1);
		Sik[i]=this->GetAij(i,k1);
	}

	for (i=0; i<size; i++){ //for each row in matrix

		if ((i==k1) && (a!=0.0)){ //special row (constraint for X_k1)
			for (j=0; j<size; j++){ //for each column
				double Sk1=Skrow[j];//this->GetAij(k1,j);
				double Sj1=Sjrow[j];//this->GetAij(j1,j);

				if      (j==j1){this->DynamicAdd(a,k1,j1);}
				else if (j==k1){this->DynamicAdd(Skrow[k1]+2.0*a*(Sjrow[k1]-Sjj)+a*a,k1,k1);}//
				else           {this->DynamicAdd(Skrow[j]-a*Sjrow[j],k1,j);}
			}
		}

		else if (i==j1){ //special row (constraint for X_j1)
			for (j=0; j<size; j++){ //for each column
				this->SetAij(0.0,i,j); //does nothing for most sparse matrix rows
				if (j==j1){this->SetAij(1.0,j1,j1);}
				if (a!=0.0){if (j==k1){this->DynamicAdd(a,j1,k1);}}
			}
		}

		else 	{		//All other rows
			if (a!=0.0){this->DynamicAdd(Sik[i]-a*Sij[i],i,k1);}
			this->DynamicAdd(0.0,i,j1);
		}
	}
	delete [] Skrow;
	delete [] Sjrow;
	delete [] Sik;
	delete [] Sij;
}
/************************************************************************
 MakeEssential1:
	revises matrix to meet essential conditions X_{j1}+a*X_{k1}=b
	affects matrix only- if condition is removed, matrix must be remade!!!!
	current implementation is pretty slow- CAN OPTIMIZE!!
	revised from Akin, pg. 403
------------------------------------------------------------------------*/
void CSparseMatrix::MakeEssential2(const CSparseMatrix *A, const int j1, const int k1, const double a, const double b){
	int i,j;
	if ((j1<0)      || (k1<0))      {MatExitGracefully("CSparseMatrix::MakeEssential2: bad indices specified",RUNTIME_ERR_MAT);}
	if ((j1==k1) && (a!=0.0)) {return;} //cannot force constraint on itself


	//cout <<"CSparseMatrix::MakeEssential2: "<<a<<endl;
	for (i=0; i<size; i++){ //for each row in matrix

		if ((i==k1) && (a!=0.0)){ //special row (constraint for X_k1)
			for (j=0; j<size; j++){ //for each column

				if      (j==j1){this->DynamicAdd(a,k1,j1);}
				else if (j==k1){this->DynamicAdd(A->GetAij(k1,k1)-2.0*a*(A->GetAij(j1,k1)-A->GetAij(j1,j1))+a*a,k1,k1);}//
				else           {this->DynamicAdd(A->GetAij(k1,j)-a*A->GetAij(j1,j),k1,j);}
			}
		}

		else if (i==j1){ //special row (constraint for X_j1)
			for (j=0; j<size; j++){ //for each column
				this->SetAij(0.0,i,j); //does nothing for most sparse matrix rows
				if (j==j1){this->SetAij(1.0,j1,j1);}
				if (a!=0.0){if (j==k1){this->DynamicAdd(a,j1,k1);}}
			}
		}

		else 	{		//All other rows
			if (a!=0.0){this->DynamicAdd(A->GetAij(i,k1)-a*A->GetAij(i,j1),i,k1);}
			this->DynamicAdd(0.0,i,j1);
		}
	}

}
/************************************************************************
 IJ_TO_K:
	Translates i (row) and j (column) index of equivalent matrix A[size][size]
	to single index [k] of sparse matrix sA[nNonZero]
	The k returned skips over k=size
------------------------------------------------------------------------*/
int CSparseMatrix::IJ_to_K   (const int i, const int j) const{
	if (i==j){return i;}
	else{
		for (int k=ija[i]; k<ija[i+1]; k++){ //searches through row i
			if (ija[k]==j){return k-1;}
		}
		return -1; //must always check
	}
}
/************************************************************************
 K_to_IJ:
	Translates k index of sparse matrix sA[nNonZero] to i,j of equivalent matrix A[size][size]
------------------------------------------------------------------------*/
bool   CSparseMatrix::K_to_IJ    (const int k, int &i, int &j) const{
	i=-2;
	int ktmp;
	if ((k<0) || (k>nNonZero)){i=j=0; return false;}
	if (k<size){i=j=k;}
	else{
		ktmp=k+1;
    j=ija[ktmp];
		for (int itmp=0; itmp<=size; itmp++){ //searches for correct row i
			if (ija[itmp]> ktmp){i=itmp-1;return true;}
      if (ija[itmp]==ktmp){i=itmp;  return true;}
			i=itmp;//takes care of k==nNonZero
		}
	}

	return true;
}

int CSparseMatrix::GetNumOffDiag (const int row) const{return ija[row+1]-ija[row];}
//-----------------------------------------------------------------------
int CSparseMatrix::GetNumRows    () const{return size;}
//-----------------------------------------------------------------------
double CSparseMatrix::GetAvgDiagCoeff  () const{
	double avg(0.0);
  for (int i=0; i<size;i++){
		avg+=sA[i];
	}
	return avg/size;
}
//-----------------------------------------------------------------------
double CSparseMatrix::GetMaxCoeff      () const{
	double max=0.0;
	for (int k=0; k<nNonZero+1; k++){
		if (k!=size){upperswap(max,fabs(sA[k]));}
	}
	return max;
}
/************************************************************************
 GetAij:
	Gets entry of equivalent matrix A, A[i][j]
	if A[i][j] not represented, returns 0.0;
	i=row, j=column
------------------------------------------------------------------------*/
double  CSparseMatrix::GetAij(const int i, const int j) const{
	if (i==j){return sA[i];}
	else{
		for (int k=ija[i]; k<ija[i+1]; k++){ //searches through row i
			//if (k==size){ExitGracefully("CSparseMatrix::GetAij: badly initialized matrix",RUNTIME_ERR_MAT);}
			if (ija[k]==j){return sA[k];}
		}
		return 0.0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Gets entry in sparse matrix sA(k) </summary>
///
/// <remarks>	rsimms, 6/21/2011. </remarks>
///
/// <param name="k">	The sparse matrix index. </param>
///
/// <returns>	The value at the specified index </returns>
////////////////////////////////////////////////////////////////////////////////////////////////////
double  CSparseMatrix::GetAk(const int k) const{
	if (k > nNonZero){return 0.0;};
	return sA[k];
}

/************************************************************************
 SetAij:
	Overwrites entry of equivalent matrix A, A[i][j] with a
	if A[i][j] not represented, returns false;
------------------------------------------------------------------------*/
bool CSparseMatrix::SetAij(const double &a, const int i, const int j){
	if (i==j){sA[i]=a;return true;}
	else {
		for (int k=ija[i]; k<ija[i+1]; k++){ //searches through row i
			if (ija[k]==j){
				sA[k]=a;
				return true;
			}
		}
		if (a==0.0){return true;}
		else       {return false;}
	}
}
/************************************************************************
 AddToAij:
	Updates entry of equivalent matrix A, A[i][j] by adding a
	if A[i][j] not represented, returns false;
------------------------------------------------------------------------*/
bool CSparseMatrix::AddToAij(const double &a, const int i, const int j){
	if (i==j){sA[i]+=a;return true;}
	else{
		for (int k=ija[i]; k<ija[i+1]; k++){ //searches through row i
			if (ija[k]==j){
				sA[k]+=a;
				return true;}
		}
	  if (a==0.0){return true;}
		else       {return false;}
	}
}
/************************************************************************
 MatVectMult:
	Sparse Matrix vector multiply

	Turned on the dirichlet condition for the !transpose case (R. Simms, Mar 16 2011)
------------------------------------------------------------------------*/
void CSparseMatrix::MatVectMult(Ironclad1DArray  x,
															  Writeable1DArray b,
																const int        N,
																const bool       transpose,
																const bool       ignore) const {
	if (size<=N){
		MatExitGracefully("CSparseMatrix::MatVectMult: Cannot multiply matrix and vector of different sizes",BAD_DATA_MAT);}

	int i,k;
	if (!transpose){
		for (i=0;i<size;i++){
			if (!(dirichlet[i])){
				b[i]=sA[i]*x[i];
				for (k=ija[i];k<ija[i+1];k++){
					//if (!(ignore && dirichlet[ija[k]])){
						b[i]+=sA[k]*x[ija[k]];
					//}
				}
			}
			else{
				b[i]=x[i];
			}
		}
	}
	else{
		int j;
		for (i=0;i<size;i++){
			if (!dirichlet[i]){
				b[i]=sA[i]*x[i];
			}
			else{
				b[i]=x[i];
			}
		}
		for (i=0;i<size;i++){
			//if (!dirichlet[i]){ //not sure if activity handled correctly here
				for (k=ija[i];k<ija[i+1];k++){
					j=ija[k];
					//if (!dirichlet[j]){
						b[j]+=sA[k]*x[i];
					//}
				}
			//}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Mat vect multiply. Wrapper using arma::vec </summary>
///
/// <remarks>	rsimms, 2/13/2011.
/// 			rsimms, 2/13/2011: Rewrote so that this code used armadillo library</remarks>
///
/// <param name="x">	The x in Ax=b. </param>
/// <param name="b">	[out] The b in Ax=b. </param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void   CSparseMatrix::MatVectMult (	vec  x, vec *b ) const {

	//Use mat vec multiply directly
	if (size<=(int)x.n_rows ){
		MatExitGracefully("CSparseMatrix::MatVectMult: Cannot multiply matrix and vector of different sizes",BAD_DATA_MAT);}

	(*b).set_size(x.n_rows);
	(*b).fill(0);

	int i,k;
	for (i=0;i<size;i++){
		if (!(dirichlet[i])){
			(*b)[i]=sA[i]*x[i];
			for (k=ija[i];k<ija[i+1];k++){
				//if (!(ignore && dirichlet[ija[k]])){
					(*b)[i]+=sA[k]*x[ija[k]];
				//}
			}
		}
		else{
			(*b)[i]=x[i];
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Mat vect multiply. Wrapper using arma::vec </summary>
///
/// <remarks> asnowdon,  </remarks>
///
/// <param name="x">	The x in Ax=b. </param>
/// <param name="b">	[out] The b in Ax=b. </param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void   CSparseMatrix::MatVectMult (	vec  x, vec *b, CSparseMatrix *SpNR ) const {

	//Use mat vec multiply directly
	if ((SpNR->size)<=(int)x.n_rows ){
		MatExitGracefully("CSparseMatrix::MatVectMult: Cannot multiply matrix and vector of different sizes",BAD_DATA_MAT);}

	(*b).set_size(x.n_rows);
	(*b).fill(0);

	int i,k;
	for (i=0;i<(SpNR->size);i++){
		if (!(SpNR->dirichlet[i])){
			(*b)[i]=(SpNR->sA[i])*x[i];
      //cout<<"i: "<<i<<"   sa: "<<SpNR->sA[i]<<"   x:"<<x[i]<<"   b:"<<(*b)[i]<<endl;
			for (k=(SpNR->ija[i]);k<(SpNR->ija[i+1]);k++){
				//if (!(ignore && dirichlet[ija[k]])){
					(*b)[i]+=(SpNR->sA[k])*x[(SpNR->ija[k])];
          //cout<<"i: "<<i<<"   ijai:: "<<SpNR->ija[i]<<"   ijak:: "<<SpNR->ija[k]<<"   sa: "<<SpNR->sA[k]<<"   x:"<<x[SpNR->ija[k]]<<"   b:"<<(*b)[i]<<endl;
				//}
			}
		}
		else{
			(*b)[i]=x[i];
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Add the components of the input matrix to the current matrix</summary>
///
/// <remarks>	rsimms, 6/21/2011.
/// 			OPTIMIZE: could probably find a better routine than using the k's and converting to IJs</remarks>
///
/// <param name="A"> The matrix to be added </param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CSparseMatrix::MatMatAdd(const CSparseMatrix &A){
	if(A.GetNumRows()!=GetNumRows()){
		MatExitGracefully("CSparseMatrix::MatVectMult: Cannot add matrices of different sizes",BAD_DATA_MAT);
	}

	// go through all the entries in the matrix that was passed in and add them
	for(int k =0; k < A.GetNumEntries(); k++){
		int i,j;
		A.K_to_IJ(k,i,j);

		DynamicAdd(A.GetAk(k),i,j);
	}
}


/************************************************************************
	ScalarMult:
		multiplies sparse matrix by a scalar value w
------------------------------------------------------------------------*/
void CSparseMatrix::ScalarMult (const double w){
	for (int k=0; k<nNonZero+1; k++){
		if (k!=size){
			sA[k]*=w;
		}
	}
}
/************************************************************************
	Print
------------------------------------------------------------------------*/
void CSparseMatrix::Print(const bool full, int sample) const{
	if (sample>size){sample=size;}
	for (int i=0;i<sample;i++){//for each row
		PrintRow(full,i);
	}
}
//-----------------------------------------------------------------------
void CSparseMatrix::PrintRow(const bool full, const int i) const{
	double tmpA;
	if (!dirichlet[i]){
		for (int j=0;j<size;j++){//for each column
			tmpA=GetAij(i,j);
			if ((full) || (tmpA!=0.0)){
				if (i==j){printf (" [%6.5f]   ", tmpA);}
				else     {printf ("  %6.5f    ", tmpA);}
			}

		}
	}
	else{
		for (int j=0;j<size;j++){
			if (i==j){tmpA=1.0;}else{tmpA=0.0;}
			if ((full) || (tmpA!=0.0)){
				if (i==j){printf (" [%6.5f]   ", tmpA);}
				else     {printf ("  %6.5f    ", tmpA);}
			}
		}
	}
	cout <<endl;
}
/************************************************************************
	SolveEst:
		result,{x}, is {b} normalized by the diagonal of A
		S is sparse matrix
		same as asolve from Press et al
		can be changed to speed up convergence (using off diagonal terms?)
------------------------------------------------------------------------*/
void CSparseMatrix::SolveEst(Ironclad1DArray b, Writeable1DArray x, const int size, const bool transpose) const{

		for (int i=0; i<size;i++){
			if (!dirichlet[i]){x[i]=(sA[i]!=0.0 ? b[i]/sA[i] : b[i]);}
			else              {x[i]=b[i];                            }
		}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Solve est:
/// 			same as James's old version
///
/// 			Jacobi/Diagonal Preconditioner
/// 			</summary>
///
/// <remarks>	rsimms, 3/18/2011. </remarks>
///
/// <param name="b">			The b. </param>
/// <param name="x">			[out] If non-null, the x resultant of b normaliazed against A</param>
/// <param name="size">			The size of the vectors. </param>
/// <param name="transpose">	The transpose yes/no. THE TRANSPOSE IS NOT IMPLEMENTED</param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CSparseMatrix::SolveEst(vec b, vec * x, const int size, const bool transpose) const{
		for (int i=0; i<size;i++){
			if (!dirichlet[i]){(*x)[i]=(sA[i]!=0.0 ? b[i]/sA[i] : b[i]);}
			else              {(*x)[i]=b[i];                            }
		}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Solve est:
/// 			same as James's old version
///
/// 			Jacobi/Diagonal Preconditioner
/// 			</summary>
///
/// <remarks>	asnowdon, 3/18/2015. </remarks>
///
/// <param name="b">			The b. </param>
/// <param name="x">			[out] If non-null, the x resultant of b normaliazed against A</param>
/// <param name="size">			The size of the vectors. </param>
/// <param name="transpose">	The transpose yes/no. THE TRANSPOSE IS NOT IMPLEMENTED</param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CSparseMatrix::SolveEst(vec b, vec * x, const int size, const bool transpose, CSparseMatrix *SpNR) const{
		for (int i=0; i<(SpNR->size);i++){
			if (!(SpNR->dirichlet[i])){(*x)[i]=((SpNR->sA[i])!=0.0 ? b[i]/(SpNR->sA[i]) : b[i]);
      //cout<<"i: "<<i<<"   size: "<<SpNR->size<<"   x: "<<(*x)[i]<<"   sA: "<<SpNR->sA[i]<<"   b: "<<b[i]<<endl;
      }
			else              {(*x)[i]=b[i];                            }
		}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Converts this matrix to an ILU matrix
/// 			DESTROYS THE ORIGINAL MATRIX
///
/// 			ILU(0) Preconditioner
/// 			Implementation as specified in
/// 			Iterative Methods for Sparse Linear Systems
/// 			Yousef Saad, 2000
/// 			page 275</summary>
///
/// <remarks>	rsimms, 3/28/2012. </remarks>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CSparseMatrix::ConvertToILU() {
	int N = GetNumRows();

	for (int i = 1; i < N; i++){
		for(int Kiter =ija[i]; Kiter<ija[i+1];Kiter++){
			int k = ija[Kiter];
			if(k >= i){
				break;
			}
			if(GetAij(k,k) == 0){
				// the solution to guassian elimination will be singular, I'm not going to implement pivoting in this
				//return a base solution
				MatExitGracefully("BCGstab: zero in diagonal",BAD_DATA_MAT);
				return;
			}

			double aik = sA[Kiter]/GetAij(k,k);

			if(aik!=0){
				sA[Kiter] = aik;

				for(int Jiter = ija[i]; Jiter < ija[i+1]; Jiter++){
					int j = ija[Jiter];
					if(j > k){
						sA[Jiter] = sA[Jiter] - aik * GetAij(k,j);
					} //else do nothing
				}
				SetAij(GetAij(i,i) - aik*GetAij(k,i),i,i); //need to solve for the diagonal separately

			} // else {} do nothing
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Solves Ax=b where A is represented in LU form
/// 			THIS MATRIX MUST ALREADY BE IN LU FORM FOR THIS METHOD TO WORK
/// 			Works with complete and incomplete LU factorizations</summary>
///
///				Solves: Ly=b for y then...
///						Ux=y for x
///
/// <remarks>	rsimms, 3/30/2012. </remarks>
///
/// <param name="b">	The rhs. </param>
///
/// <returns>	The solution vector. </returns>
////////////////////////////////////////////////////////////////////////////////////////////////////
vec CSparseMatrix::SolveLUSystem(vec b) const{
	int N = (int) (b.n_rows);
	vec x = b;

	//START OF FORWARD SOLVE
	x(0) = x(0) / GetAij(0,0);
	for(int i =1; i< N;i++){// row iterator
		for(int Kiter =ija[i]; Kiter<ija[i+1];Kiter++){
			int k = ija[Kiter];
			if(k >= i){
				break;
			}
			x(i) = x(i) - sA[Kiter] * x(k);
		}
	}

	//START OF BACK SOLVE
	for(int i =N-1; i>=0;i--){// row iterator
		for(int Kiter =ija[i]; Kiter<ija[i+1];Kiter++){
			int k = ija[Kiter];
			if(k > i){
				x(i) = x(i) - sA[Kiter] * x(k);
			}
		}
		x(i) = x(i) / sA[i];
	}

	return x;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Solves Ax=b where A is represented in LU form
/// 			THIS MATRIX MUST ALREADY BE IN LU FORM FOR THIS METHOD TO WORK
/// 			Works with complete and incomplete LU factorizations</summary>
///
///				Solves: Ly=b for y then...
///						Ux=y for x
///
/// <remarks>	asnowdon, 3/30/2015. </remarks>
///
/// <param name="b">	The rhs. </param>
///
/// <returns>	The solution vector. </returns>
////////////////////////////////////////////////////////////////////////////////////////////////////
vec CSparseMatrix::SolveLUSystem(vec b, CSparseMatrix *SpNR) const{
	int N = (int)( b.n_rows);
	vec x = b;

	//START OF FORWARD SOLVE
	x(0) = x(0) / GetAij(0,0);
	for(int i =1; i< N;i++){// row iterator
		for(int Kiter =SpNR->ija[i]; Kiter<SpNR->ija[i+1];Kiter++){
			int k = SpNR->ija[Kiter];
			if(k >= i){
				break;
			}
			x(i) = x(i) - SpNR->sA[Kiter] * x(k);
		}
	}

	//START OF BACK SOLVE
	for(int i =N-1; i>=0;i--){// row iterator
		for(int Kiter =SpNR->ija[i]; Kiter<SpNR->ija[i+1];Kiter++){
			int k = SpNR->ija[Kiter];
			if(k > i){
				x(i) = x(i) - SpNR->sA[Kiter] * x(k);
			}
		}
		x(i) = x(i) / SpNR->sA[i];
	}

	return x;
}

/************************************************************************
 CalculateAdjustedNorm:
	Calculates norm of a vector (dirichlet terms removed)
------------------------------------------------------------------------*/
double CSparseMatrix::CalculateAdjustedNorm(Ironclad1DArray b, int size, const normtype type) const{
	int i,isamax;
    double norm;
	if      (type==VECTOR_MAGNITUDE_NORM){

		norm=0.0;for (i=0;i<size;i++){norm+=b[i]*b[i];}//if (!dirichlet[i]){norm+=b[i]*b[i];}}
		return sqrt(norm);
	}
	else if (type==LARGEST_COMPONENT_NORM){

		isamax=0;for (i=1;i<size;i++){if (fabs(b[i])>fabs(b[isamax])){isamax=i;}}//if (!dirichlet[i]){if (fabs(b[i])>fabs(b[isamax])){isamax=i;}}}
		return fabs(b[isamax]);
	}
	else{
		return 0.0;
	}
}
/************************************************************************
 BCG:
	Sparse Preconditioned Biconjugate Gradient Solver
	S is sparse matrix representing A[size][size] matrix
	x[size] is output vector
	b[size] is RHS
	BCGtype is testing criteria
		if RELATIVE_RESIDUAL,       |Ax-b|/|b|         <BCG_tol
		if RELATIVE_TRANS_RESIDUAL, |At(Ax-b)|/|(At*b)|<BCG_tol
		if NORM_ERROR,              |xerr|/|x|         <BCG_tol
		if MAX_ERROR,               xerr_max/|x|_max   <BCG_tol

	MUST OPTIMIZE!!!!!
------------------------------------------------------------------------*/
void CSparseMatrix::BCG(						Ironclad1DArray     b,
												Writeable1DArray    x,
												const int           size,
												const BCGtestparam  BCGtype,
												double             &err,
//												double              normalize,
												const double        BCG_tol) const{

    double ak,akden;
    double bk,bkden(1.0),bknum;
    double bnorm,dxnorm,xnorm,zm1norm,znorm;
	int    j;
	bkden=1.0;
	const double eps=1e-14;

	if (size!=this->size){MatExitGracefully("BCG: Matrix and vector not same size",BAD_DATA_MAT);}

	//double *bcopy=new double [size];
	/*for (j=0; j<size;j++){
		bcopy[j]=b[j];*/
		/*for (int i=0;i<size;i++){
			if (dirichlet[i]){
				bcopy[j]-=b[i]*GetAij(i,j); //do not want to do this!!
			}
		}*/
	/*}*/

	int iter=0;
	MatVectMult(x,r,size,false,true);

	for (j=0; j<size;j++){
		if (!dirichlet[j]){
			r [j]=b[j]-r[j];
			rr[j]=r[j];
		}
	}
	if (CalculateAdjustedNorm(r,size,VECTOR_MAGNITUDE_NORM)==0.0){
		cout <<"converged"<<endl;
		return;
		//ExitGracefully("BCG: residual norm is zero (already converged)",BAD_DATA_MAT);
	}

	//MatVectMult(r,rr,size,false,true); //minimum residual variant if uncommented (blows up dirichlet conditions)

	if      (BCGtype==RELATIVE_RESIDUAL){

		bnorm=CalculateAdjustedNorm(b,size,VECTOR_MAGNITUDE_NORM);
		if (bnorm==0.0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}

		SolveEst(r,z,size,false);
	}
	else if (BCGtype==RELATIVE_TRANS_RESIDUAL){

		SolveEst(b,z,size,false);
		bnorm=CalculateAdjustedNorm(z,size,VECTOR_MAGNITUDE_NORM);
		if (bnorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}

		SolveEst(r,z,size,false);
	}
	else if (BCGtype==NORM_ERROR){

		SolveEst(b,z,size,false);
		bnorm=CalculateAdjustedNorm(z,size,VECTOR_MAGNITUDE_NORM);
		if (bnorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}

		SolveEst(r,z,size,false);

		znorm=CalculateAdjustedNorm(z,size,VECTOR_MAGNITUDE_NORM);
		if (znorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}
	}
	else if (BCGtype==MAX_ERROR){

		SolveEst(b,z,size,false);
		bnorm=CalculateAdjustedNorm(z,size,LARGEST_COMPONENT_NORM);
		if (bnorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}

		SolveEst(r,z,size,false);
		znorm=CalculateAdjustedNorm(z,size,LARGEST_COMPONENT_NORM);
		if (znorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}
	}
	else {MatExitGracefully("BCG: Illegal value for tolerance",RUNTIME_ERR_MAT);}

	//cout <<fixed <<setprecision(6);

	while (iter<=BCG_max_iter){
		iter++;
		SolveEst(rr,zz,size,true); //true indicates transpose matrix At

		bknum=0.0;
		for (j=0;j<size;j++){bknum+=z[j]*rr[j];}
		if (bknum==0.0){cout <<z[0]<<" "<<rr[0]<<endl;MatExitGracefully("BCG: beta numerator is zero",BAD_DATA_MAT);}

		//Calculate coeff bk and direction vectors p and pp--------------------------
		if (iter==1){
			for (j=0;j<size;j++){
				if (!dirichlet[j]){
					p [j]=z [j];
					pp[j]=zz[j];
				}
			}
		}
		else{
			bk=bknum/bkden;
			for (j=0;j<size;j++){
				if (!dirichlet[j]){
					p [j]=bk*p [j]+z [j];
					pp[j]=bk*pp[j]+zz[j];
				}
				else{
					p[j]=0.0;
					pp[j]=0.0;
				}
			}
		}

		//Calculate coeff ak, new iterate x and new residuals r and rr---------------
		bkden=bknum;

		MatVectMult(p,z,size,false,true);

		akden=0.0;
		for (j=0;j<size;j++){if (!dirichlet[j]){akden+=z[j]*pp[j];}} //{S*p}*{S/rr}
		if (akden==0){MatExitGracefully("BCG: alpha denominator is zero",BAD_DATA_MAT);}

		ak=bknum/akden;

		MatVectMult(pp,zz,size,true,true);
		for (j=0;j<size;j++){
			if (!dirichlet[j]){
				x [j]+=ak*p [j];
				r	[j]-=ak*z [j];
				rr[j]-=ak*zz[j];
			}
			else{
				x[j]=b[j];
				r[j]=0.0;
				rr[j]=0.0;
			}
		}

		//Solve A*z=r and check stopping criteria------------------------------------
		SolveEst(r,z,size,false);

		if       (BCGtype==RELATIVE_RESIDUAL){

			err=CalculateAdjustedNorm(r,size,VECTOR_MAGNITUDE_NORM)/bnorm;
		}
		else if  (BCGtype==RELATIVE_TRANS_RESIDUAL){

			err=CalculateAdjustedNorm(z,size,VECTOR_MAGNITUDE_NORM)/bnorm;
		}
		else if ((BCGtype==NORM_ERROR) ||
			       (BCGtype==MAX_ERROR)){

			normtype ntype;
			if (BCGtype==NORM_ERROR){ntype=VECTOR_MAGNITUDE_NORM;}
			else                    {ntype=LARGEST_COMPONENT_NORM;}

			zm1norm=znorm;
			znorm  =CalculateAdjustedNorm(z,size,ntype);
			if (fabs(zm1norm-znorm)>eps*znorm){
				dxnorm=fabs(ak)*CalculateAdjustedNorm(p,size,ntype);
				err=znorm/fabs(zm1norm-znorm)*dxnorm;
			}
			else{
				err=znorm/bnorm;
				continue;
			}
			xnorm=CalculateAdjustedNorm(x,size,ntype);
			if (err<=0.5*xnorm){err/=xnorm;}
			else{
				err=znorm/bnorm;
				continue;
			}
		}
		//cout <<"BCG iter: " <<iter<<" "<<err<<endl;
		if (err<=BCG_tol){break;}

	}//end while
	//cout <<"BCG_iter:"<<iter<<endl;

	//cout <<"BCG Condition number"<<endl;
	if (iter>=BCG_max_iter){MatExitGracefully("BCG: Unable to solve system of equations",RUNTIME_ERR_MAT);}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	BCGstabalized implementation
///
/// 			as taken from van der Vortst (1992)
///
/// 			only wise to use relative residual error at the moment</summary>
///
/// <remarks>	rsimms, 3/18/2011. </remarks>
///
/// <param name="b">		The rhs. </param>
/// <param name="x">		[in,out] The Ax=b. Fed in as the initial guess</param>
/// <param name="size">		The size of b and x. </param>
/// <param name="err">		[out] The error. </param>
/// <param name="BCG_tol">	The bcg tolerance. </param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CSparseMatrix::BCGstab(
	arma::vec				    b,
	arma::vec					&x,
	const int           size,
	double				&err,
	const double        BCG_tol,
	preconditioner preconditionerType,
	CSparseMatrix*		ILUmatrix
	) const{
	if (size!=this->size){MatExitGracefully("BCG: Matrix and vector not same size",BAD_DATA_MAT);}
	//check the initial conditions to make sure x is an initial condition
	if(x.n_rows != b.n_rows){ x.set_size(b.n_rows); x.fill(0.0);}

	//begin the algorthim
	double bnorm = norm(b,2);
  if(bnorm == 0){bnorm = 1;}
	int n =(int)( b.n_rows );

  double head_sum1 = 0;
  double head_sum2 = 0;
  for(int i = 0; i<16; i++)
  {
    head_sum1 += x[i];
  }

	vec productAx;
	productAx = x;

	MatVectMult(x,&productAx);//,ILUmatrix); // equivelent to Ax0 product
  //cout<<"Ax: "<<productAx<<endl;
	vec res = b - productAx;
  //cout<<"res: "<<res<<endl;
	vec rr = res;

	/*!!! step 1,2,3 complete: x0 = x, r0 = b-Ax0, and r^ = r0   !!!*/

	// for the i-1's
	double rho_1 = 1;
	double alpha = 1;
	double omega = 1;

	// for the i-2's
	double rho_2 = 1;

	/*!!! step 4 complete: alpha, omega, rho1 and rho2 set   !!!*/

	if (norm(res,2)==0.0){
		//cout <<"converged at res normalization"<<endl;
		return;
		//ExitGracefully("BCG: residual norm is zero (already converged)",BAD_DATA_MAT);
	}

	vec p(n);
	p.fill(0.0);
	vec v(n);
	v.fill(0.0);

	/*!!! step 5 complete: o and v vectors initialized to 0  !!!*/

	int iter=0;
	while (iter<=BiCGSTAB_max_iter){ //!!! begining of step 6
		iter++;

		rho_1 = dot(rr,res);

		if (iter == 1){
			p = res;
		} else
		{
			double beta = (rho_1 / rho_2) * (alpha/omega);
			p = res + beta * (p - omega *v);
		}

		/*!!! step 6,7 complete: rho, beta and p calculated !!!*/

		// Apply preconditioner to p
		vec phat(n);
		switch (preconditionerType){
		case ILU_ZERO:
      phat = ILUmatrix->SolveLUSystem(p,ILUmatrix);
			//phat = ILUmatrix->SolveLUSystem(p);
			break;
		case NONE_MAT:
			phat = p;
			break;
		case DIAGONAL:
      SolveEst(p,&phat,n,false,ILUmatrix);
			//SolveEst(p,&phat,n,false);
			break;
		default: //NONE
			phat = p;
			break;
		}

		MatVectMult(phat,&v);//,ILUmatrix); //change to phat for preconditioning
		/*!!! step 8 complete						 !!!*/

		alpha = rho_1 / (dot(rr,v));
		/*!!! step 9 complete						 !!!*/

		vec s = res - alpha*v;
		/*!!! step 10 complete						 !!!*/

		// Apply preconditioner to s
		vec shat(n);
		switch (preconditionerType){
		case ILU_ZERO:
      shat = ILUmatrix->SolveLUSystem(s,ILUmatrix);
			//shat = ILUmatrix->SolveLUSystem(s);
			break;
		case NONE_MAT:
			shat = s;
			break;
		case DIAGONAL:
      SolveEst(s,&shat,n,false,ILUmatrix);
			//SolveEst(s,&shat,n,false);
			break;
		default: //NONE
			shat = s;
			break;
		}

		vec t(n);
		MatVectMult(shat,&t,ILUmatrix); // change to shat for preconditioning
		/*!!! step 11 complete						 !!!*/

		omega = dot(t,s) / dot(t,t);
		/*!!! step 12 complete						 !!!*/

		x = x + alpha*phat+omega*shat; // these would be the preconditioned p and s if we were using a preconditioner...
		/*!!! step 13 complete						 !!!*/
    //cout<<"iter: "<<iter<<"   x: "<<x<<endl;
		res = s - omega * t;
		/*!!! step 14 complete						 !!!*/

		// update the vectors and coeffecients
		rho_2 = rho_1;

		err=norm(res,2)/bnorm;
		/*!!! step 15 complete						 !!!*/

		//cout <<"BCG iter: " <<iter<<" err: "<<err<<endl;
		if (err<=BCG_tol){break;}
	}//end while
	//cout <<"BCG_iter:"<<iter<<endl;

  for(int i = 0; i<16; i++)
  {
    head_sum2 += x[i];
  }
  //cout<<setprecision(12)<<"s_B1: "<<((head_sum1 + (366*16)) * 1000 * 0.363)<<"   s_B2: "<<((head_sum2 + (366*16)) * 1000 * 0.363)<<"   diff: "<<((head_sum2 + (366*16)) * 1000 * 0.363)-((head_sum1 + (366*16)) * 1000 * 0.363)<<endl;
  //cout<<setprecision(12)<<"h_B1: "<<head_sum1<<"   h_B2: "<<head_sum2<<"   diff: "<<head_sum2-head_sum1<<endl;


	//cout <<"BCG Condition number"<<endl;
	if (iter>=BiCGSTAB_max_iter){MatExitGracefully("BCG: Unable to solve system of equations",RUNTIME_ERR_MAT);}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	BCGstabalized implementation
/// 			defaults to no preconditioner
///
/// 			as taken from van der Vortst (1992)
///
/// 			only wise to use relative residual error at the moment</summary>
///
/// <remarks>	rsimms, 3/18/2011. </remarks>
///
/// <param name="b">		The rhs. </param>
/// <param name="x">		[in,out] The Ax=b. Fed in as the initial guess</param>
/// <param name="size">		The size of b and x. </param>
/// <param name="err">		[out] The error. </param>
/// <param name="BCG_tol">	The bcg tolerance. </param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CSparseMatrix::BCGstab(
	arma::vec				    b,
	arma::vec					&x,
	const int           size,
	double				&err,
	const double        BCG_tol
	) const{
		CSparseMatrix* nullMatrix = new CSparseMatrix();
		BCGstab(b,x,size,err,BCG_tol,NONE_MAT,nullMatrix);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	BCGstabalized implementation
/// 			defaults to no preconditioner
///
/// 			as taken from van der Vortst (1992)
///
/// 			only wise to use relative residual error at the moment</summary>
///
/// <remarks>	asnowdon, 3/18/2015. </remarks>
///
/// <param name="b">		The rhs. </param>
/// <param name="x">		[in,out] The Ax=b. Fed in as the initial guess</param>
/// <param name="size">		The size of b and x. </param>
/// <param name="err">		[out] The error. </param>
/// <param name="BCG_tol">	The bcg tolerance. </param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CSparseMatrix::BCGstab(
	arma::vec				   b,
	arma::vec					&x,
	const int          size,
	double				    &err,
	const double       BCG_tol,
  CSparseMatrix     *SparseNewtonMat
	) const{

		BCGstab(b,x,size,err,BCG_tol,NONE_MAT,SparseNewtonMat);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Bicongugate gradient method (use instead of the wrapper)</summary>
///
/// <remarks>	rsimms, 2/9/2011</remarks>
///
/// <param name="b">		The rhs. </param>
/// <param name="x">		[out] The Ax=b. </param>
/// <param name="size">		The size of b and x. </param>
/// <param name="BCGtype">	The bcg type. </param>
/// <param name="err">		[out] The error. </param>
/// <param name="BCG_tol">	The bcg tolerance. </param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CSparseMatrix::BCG(
												arma::vec				    b,
												arma::vec					&x,
												const int           size,
												const BCGtestparam  BCGtype,
												double				&err,
											//	double              normalize,
												const double        BCG_tol
												) const{



    double ak,akden;
    double bk,bkden(1.0),bknum;
    double bnorm,dxnorm,xnorm,zm1norm,znorm;
	int    j;
	bkden=1.0;
	const double eps=1e-14;

	if (size!=this->size){MatExitGracefully("BCG: Matrix and vector not same size",BAD_DATA_MAT);}

	int iter=0;
	MatVectMult(x,rArma);

	for (j=0; j<size;j++){
		if (!dirichlet[j]){
			(*rArma) [j]=b[j]-(*rArma)[j];
			(*rrArma)[j]=(*rArma)[j];
		}
	}
	if (norm(*rArma,2)==0.0){
		cout <<"converged"<<endl;
		return;
		//ExitGracefully("BCG: residual norm is zero (already converged)",BAD_DATA_MAT);
	}

	//MatVectMult(r,rr,size,false,true); //minimum residual variant if uncommented (blows up dirichlet conditions)

	if      (BCGtype==RELATIVE_RESIDUAL){

		bnorm=norm(b,2);
		if (bnorm==0.0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}

		SolveEst(*rArma,zArma,size,false);
	}
	else if (BCGtype==RELATIVE_TRANS_RESIDUAL){

		SolveEst(b,zArma,size,false);
		bnorm=norm(*zArma,2);
		if (bnorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}

		SolveEst(*rArma,zArma,size,false);
	}
	else if (BCGtype==NORM_ERROR){

		SolveEst(b,zArma,size,false);
		bnorm=norm(*zArma,2);
		if (bnorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}

		SolveEst(*rArma,zArma,size,false);

		znorm=norm(*zArma,2);
		if (znorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}
	}
	else if (BCGtype==MAX_ERROR){

		SolveEst(b,zArma,size,false);
		bnorm=norm(*zArma,"inf");
		if (bnorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}

		SolveEst(*rArma,zArma,size,false);
		znorm=norm(*zArma,"inf");
		if (znorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA_MAT);}
	}
	else {MatExitGracefully("BCG: Illegal value for tolerance",RUNTIME_ERR_MAT);}

	//cout <<fixed <<setprecision(6);

	while (iter<=BCG_max_iter){
		iter++;
		SolveEst(*rrArma,zzArma,size,true); //true indicates transpose matrix At

		bknum=0.0;
		for (j=0;j<size;j++){bknum+=(*zArma)[j]*(*rrArma)[j];}

		if (bknum==0.0){cout <<z[0]<<" "<<rr[0]<<endl;MatExitGracefully("BCG: beta numerator is zero",BAD_DATA_MAT);}

		//Calculate coeff bk and direction vectors p and pp--------------------------
		if (iter==1){
			*pArma = *zArma;
			*ppArma = *zzArma;
		}
		else{
			bk=bknum/bkden;

			*pArma = bk*(*pArma) + (*zArma);
			*ppArma = bk*(*ppArma) + (*zzArma);

			for (j=0;j<size;j++){
				if (dirichlet[j]){
					(*pArma)[j]=0.0;
					(*ppArma)[j]=0.0;
				}
			}
		}

		//Calculate coeff ak, new iterate x and new residuals r and rr---------------
		bkden=bknum;

		MatVectMult(*pArma,zArma);

		akden=0.0;
		for (j=0;j<size;j++){if (!dirichlet[j]){akden+=(*zArma)[j]*(*ppArma)[j];}} //{S*p}*{S/rr}
		if (akden==0){MatExitGracefully("BCG: alpha denominator is zero",BAD_DATA_MAT);}

		ak=bknum/akden;

		MatVectMult(pp,zz,size,true,true);
		for (j=0;j<size;j++){
			if (!dirichlet[j]){
				x [j]+=ak*(*pArma) [j];
				(*rArma)	[j]-=ak*(*zArma) [j];
				(*rrArma)   [j]-=ak*(*zzArma)[j];
			}
			else{
				x[j]=b[j];
				(*rArma)[j]=0.0;
				(*rrArma)[j]=0.0;
			}
		}

		//Solve A*z=r and check stopping criteria------------------------------------
		SolveEst(*rArma,zArma,size,false);

		if       (BCGtype==RELATIVE_RESIDUAL){

			err=norm(*rArma,2)/bnorm;
		}
		else if  (BCGtype==RELATIVE_TRANS_RESIDUAL){

			err=norm(*zArma,2)/bnorm;
		}
		else if ((BCGtype==NORM_ERROR) ||
			       (BCGtype==MAX_ERROR)){

			normtype ntype;
			if (BCGtype==NORM_ERROR){ntype=VECTOR_MAGNITUDE_NORM;
				zm1norm=znorm;
				znorm  =norm(*zArma,2);
				if (fabs(zm1norm-znorm)>eps*znorm){
					dxnorm=fabs(ak)*norm(*pArma,2);
					err=znorm/fabs(zm1norm-znorm)*dxnorm;
				}
				else{
					err=znorm/bnorm;
					continue;
				}
				xnorm=norm(x,2);
				if (err<=0.5*xnorm){err/=xnorm;}
				else{
					err=znorm/bnorm;
					continue;
				}
			}
			else                    {ntype=LARGEST_COMPONENT_NORM;

				zm1norm=znorm;
				znorm  =norm(*zArma,"inf");
				if (fabs(zm1norm-znorm)>eps*znorm){
					dxnorm=fabs(ak)*norm(*pArma,"inf");
					err=znorm/fabs(zm1norm-znorm)*dxnorm;
				}
				else{
					err=znorm/bnorm;
					continue;
				}
				xnorm=norm(x,"inf");
				if (err<=0.5*xnorm){err/=xnorm;}
				else{
					err=znorm/bnorm;
					continue;
				}
			}
		}
		//cout <<"BCG iter: " <<iter<<" "<<err<<endl;
		if (err<=BCG_tol){break;}

	}//end while

	if (iter>=BCG_max_iter){MatExitGracefully("BCG: Unable to solve system of equations",RUNTIME_ERR_MAT);}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Triangle band solve. Only for one off diagonal. Solves A\b=x </summary>
///
/// <remarks>	rsimms, 8/8/2011.   NO SAFETY TO CHECK FOR OFF DIAGONALS</remarks>
///
/// <param name="b">	The right hand side. </param>
/// <param name="x">	[out] The vector to be solved for. </param>
////////////////////////////////////////////////////////////////////////////////////////////////////
void CSparseMatrix::TriBandSolve(
	arma::vec const			    b,
	arma::vec					&x
	)
{
	if (b.n_rows!=this->size){MatExitGracefully("TriBandSolve: Matrix and vector not same size",BAD_DATA_MAT);}

	//copy over the b vector which will be used to build the solution
	x = b;

	// create a new matrix that is reduced in size, use reduced notation
	mat A(size,3);
	//A.fill(0);
	A(0,1) = this ->GetAij(0,0);
	A(0,2) = this ->GetAij(0,1);

	// left = off diagonal, mid = diagonal, right = off diagonal positive
	for(int k = 1; k< size-1;k++){
		A(k,0) = GetAij(k,k-1);
		A(k,1) =GetAij(k,k);
		A(k,2)=GetAij(k,k+1);
	}
	A(size-1,0)=GetAij(size-1,size-2);
	A(size-1,1)=GetAij(size-1,size-1);

	// Solve this beast, forward first
	for(int k = 1; k < size; k++){
		double factor = A(k,0)/A(k-1,1);
		A(k,1) = A(k,1) - factor * A(k-1,2);
		x(k) = x(k) - factor*x(k-1);
	}

	// back solve
	for(int k = size-2; k >= 0;k--){
		double factor = A(k,2)/A(k+1,1);
		x(k) = x(k) - factor * x(k+1);
	}

	//normalize to ones on the diagonal
	for(int k = 0; k < size; k++){
		x(k) = x(k)/A(k,1);
	}

	//done like bacon
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// <summary>	Gets the dense version of this sparse matrix
/// 			Useful for debuging </summary>
///
/// <remarks>	rsimms, 1/27/2012. </remarks>
///
/// <returns>	A dense version of the original matrix. </returns>
////////////////////////////////////////////////////////////////////////////////////////////////////
arma::mat CSparseMatrix::ConvertToDense(){
	mat denseMatrix(GetNumRows(),GetNumRows());

	for(int j = 0; j < GetNumRows(); j++){
		for(int i = 0; i < GetNumRows(); i++){
			denseMatrix(i,j) = GetAij(i,j);
		}
	}

	return denseMatrix;
}

void   MatExitGracefully(std::string statement, badcode code)
{
//TMP DEBUG
}
void   MatExitGracefullyIf(bool condition, std::string statement, badcode code)
{
//TMP DEBUG
}
#endif
