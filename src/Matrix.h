/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2022 the Raven Development Team
----------------------------------------------------------------*/

#ifndef MATRIX
#define MATRIX
#include "RavenInclude.h"

typedef const double*        Unchangeable1DArray; //unchangeable but movable
typedef       double* const  Writeable1DArray;    //unmoveable but changeble
typedef const double* const  Ironclad1DArray;     //unmodifiable

typedef const double* const* Unchangeable2DArray;
typedef       double** const Writeable2DArray;
typedef const double* const* const Ironclad2DArray;

void   MatVectMult    (Ironclad2DArray  A,
											 Ironclad1DArray  x,
											 const int N,
											 bool transpose,
	                     Writeable1DArray b);
void   MatVectMult    (Ironclad2DArray  A,
											 Ironclad1DArray  x,
											 const int N,
											 const int M,
	                     Writeable1DArray b);
void   MatMult        (Ironclad2DArray  A,
											 Ironclad2DArray  B,
											 const int				N,
											 const int				M,
											 const int				P,
											 Writeable2DArray C);
void   ScalarMatMult  (Ironclad2DArray  A,
	                     const double &muliplier,
											 const int				N,
									     const int				M,
											 Writeable2DArray C);
void   MatAdd         (Ironclad2DArray  A,
											 Ironclad2DArray  B,
											 const int				N,
											 const int				M,
											 Writeable2DArray C);
void  CopyMatrix      (Ironclad2DArray  A,
											 const int				N,
											 const int				M,
										 	 Writeable2DArray C);
void TransposeMat     (Ironclad2DArray  A,
										   const int        N,
								       const int        M,
											 Writeable2DArray AT);
void AllocateMatrix(const int N, const int M, double ** &A);
void DeleteMatrix  (const int N, const int M, double ** A);
bool						  SVD(Ironclad2DArray A,
	                    Ironclad1DArray b,
											Writeable1DArray x,
											const int        size,
											const double     SVTolerance);
#endif
