//SparseMatrix.h
#ifndef SPARSEMATRIX
#define SPARSEMATRIX
#include "RavenInclude.h"

#ifdef _ARMADILLO_
#include<armadillo>
using namespace arma;



//#include "GeneralInclude.h"

const int    BCG_max_iter=60;   //Maximum iterations for sparse BCG solver
const int    BiCGSTAB_max_iter=60;  //Maximum iterations for sparse BiCGSTAB solver

enum BCGtestparam{RELATIVE_RESIDUAL,RELATIVE_TRANS_RESIDUAL,NORM_ERROR,MAX_ERROR};

enum normtype  {
	VECTOR_MAGNITUDE_NORM,
	LARGEST_COMPONENT_NORM
};
enum preconditioner {
	NONE_MAT,
	DIAGONAL,
	ILU_ZERO
};
enum  badcode {//error codes
	NO_ERROR,
	COMPLETED,
	RUNTIME_ERR_MAT,
	TOO_MANY,
	NOT_ENOUGH,
	BAD_DATA_MAT,
	BAD_SOL_READ,
	SINGMAT,
	OUT_OF_MEMORY_MAT,
	NO_INPUT_FILE,
	VIRTUAL_ERROR,
	MESH_GEN_ERR,
	STUB_MAT,
	OTHER
};

  //type definitions-------------------
typedef const double *               Unchangeable1DArray; //unchangeable but movable
typedef       double * const         Writeable1DArray;    //unmoveable but changeble
typedef const double * const         Ironclad1DArray;     //unmodifiable

typedef const double * const *       Unchangeable2DArray;
typedef       double *       * const Writeable2DArray;
typedef const double * const * const Ironclad2DArray;


/*!************************************************
 \class CSparseMatrix
 \brief Sparse Matrix Data Abstraction

 Adapted from Press et al. Numerical Recipes for c++, pg.81.

 Dynamic Creation developed by James R. Craig
 BiCGstab, modification to use armadillo for vector math added by Richard B. Simms
 ***************************************************/
class CSparseMatrix
{
 private: /*------------------------------------------------*/
	double *sA;         ///Storage
	int    *ija;
  double *b;
	int     size;				///actual size of original square matrix
	int     nNonZero;		///number of non-zero elements
  int     buffersize; ///buffer used to speed up memory allocation
	bool   *dirichlet;  ///dirichlet entries (used for dirichlet conditions)

	double *p,*pp,*r,*rr,*z,*zz;

	//residual vector
		arma::vec * rArma;
		arma::vec * zArma;
		arma::vec * rrArma;
		arma::vec * zzArma;
		arma::vec * pArma;
		arma::vec * ppArma;

	double CalculateAdjustedNorm(Ironclad1DArray b, int size, const normtype type) const;
	void   SolveEst   (Ironclad1DArray b, Writeable1DArray x, const int size, const bool transpose) const;
	void   SolveEst(arma::vec b, arma::vec * x, const int size, const bool transpose) const;
  void   SolveEst(arma::vec b, arma::vec * x, const int size, const bool transpose, CSparseMatrix *SpNR) const;

 public: /*------------------------------------------------*/

  CSparseMatrix();
  CSparseMatrix(double *spA, int *indices, double *bvec);
  CSparseMatrix(double *spA, int *indices, double *bvec, int sizen);
	CSparseMatrix(const CSparseMatrix &A);
	CSparseMatrix(Ironclad2DArray A, const int size, const double thresh);
 ~CSparseMatrix();

	int    GetNumEntries() const;

	bool   K_to_IJ    (const int k, int &i, int &j) const;
	int    IJ_to_K    (const int i, const int j) const;
	double GetAij     (const int i, const int j) const;
	double GetAk	  (const int k) const;

	int    GetNumOffDiag    (const int row) const;
	int    GetNumRows       () const;
	double GetAvgDiagCoeff  () const;
	double GetMaxCoeff      () const;

	bool   AddToAij         (const double &a, const int i, const int j);
	bool   SetAij           (const double &a, const int i, const int j);

	void   MakeDirichlet    (const int *rows, const int ndirichlet);
	void   MakeDirichletWorkAround(const int *rows, const int ndirichlet);
	void   MakeEssential1   (const int j, const int k, const double ac, const double bc  );
	void   MakeEssential2   (const CSparseMatrix *A, const int j1, const int k1, const double a, const double b);

	void   DynamicInitialize(const int N);
  void   DynamicInitialize(Ironclad1DArray Adiag,  const int N);
  void   DynamicAdd       (const double &a, const int i, const int j);
  void   DynamicInsert    (const double &a, const int i, const int j);

	void   ScalarMult       (const double     w);
	void   MatVectMult      (Ironclad1DArray  x,
		                       Writeable1DArray b,
													 const int        size,
													 const bool       transpose,
													 const bool       ignore) const;

	void   MatVectMult      (	arma::vec  x,
								arma::vec *b ) const;
  void   MatVectMult (	vec  x, vec *b, CSparseMatrix *SpNR ) const;
	void   MatMatAdd		(   const CSparseMatrix &A);

	void   Print            (const bool full, int sample) const;
	void   PrintRow         (const bool full, const int i) const;

	void   BCG              (						Ironclad1DArray     b,
													Writeable1DArray    x,
													const int           size,
													const BCGtestparam  BCGtype,
													double				&err,
													const double        tolerance) const;

	void   BCG              (						arma::vec					b,
													arma::vec					&x,
													const int           size,
													const BCGtestparam  BCGtype,
													double				&err,
													const double        tolerance) const;

	void   BCGstab         (						arma::vec					b,
													arma::vec					&x,
													const int           size,
													double				&err,
													const double        tolerance,
													preconditioner		precondition,
													CSparseMatrix*		ILUmatrix=0) const;

	void   BCGstab         (						arma::vec					b,
													arma::vec					&x,
													const int           size,
													double				&err,
													const double        tolerance) const;

  void   BCGstab         (arma::vec				    b,
	                        arma::vec					 &x,
	                        const int           size,
	                        double				     &err,
	                        const double        BCG_tol,
                          CSparseMatrix      *SparseNewtonMat) const;

	void TriBandSolve			(					arma::vec const			    b,
													arma::vec					&x );

	arma::mat ConvertToDense	();
	void ConvertToILU();
	arma::vec SolveLUSystem			(				arma::vec b) const;
  arma::vec SolveLUSystem     (arma::vec b, CSparseMatrix *SpNR) const;

};

/*!**************************************************
 *  Class CSparseMatrix
 *  Sparse Matrix Data Abstraction
 */


  //library exit strategy functions (functions taken from MatrixSolvers.cpp)------

	void   MatExitGracefully(std::string statement, badcode code);

	void   MatExitGracefullyIf(bool condition, std::string statement, badcode code);

#else
//Dummy class
class CSparseMatrix
{

};
#endif

#endif
