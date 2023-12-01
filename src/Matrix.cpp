/*----------------------------------------------------------------
	Raven Library Source Code
	Copyright (c) 2008-2022 the Raven Development Team
	----------------------------------------------------------------*/
#include "Matrix.h"
/************************************************************************
 MatMult:
	Multiplies NxN square matrix and Nx1 vector A*x. Returns b
-----------------------------------------------------------------------*/
void   MatVectMult(Ironclad2DArray A,Ironclad1DArray x,Writeable1DArray b,const int N,bool transpose)
{
	int i,j;
	if(!transpose) {
		for(i=0;i<N;i++) {
			b[i]=0.0;  for(j=0; j<N;j++) { b[i]+=A[i][j]*x[j]; }
		}
	}
	else {
		for(i=0;i<N;i++) {
			b[i]=0.0;  for(j=0; j<N;j++) { b[i]+=A[j][i]*x[j]; }
		}
	}
}
/************************************************************************
 MatVectMult:
	Multiplies non-square NxM matrix and Mx1 vector A*x. Returns b (Nx1)
-----------------------------------------------------------------------*/
void   MatVectMult(Ironclad2DArray A,Ironclad1DArray x,Writeable1DArray b,const int N,const int M)
{
	for(int i=0;i<N;i++) {
		b[i]=0.0;  for(int j=0; j<M;j++) { b[i]+=A[i][j]*x[j]; }
	}
}
/************************************************************************
 MatMult:
	Multiplies NxM matrix A times MxP matrix B. Returns C (NxP)
-----------------------------------------------------------------------*/
void   MatMult(Ironclad2DArray A,Ironclad2DArray B,const int N,const int M,const int P,Writeable2DArray C)
{
	for(int n=0; n<N; n++) {
		for(int p=0; p<P; p++) {
			C[n][p]=0.0;
			for(int m=0; m<M; m++) {
				C[n][p]+=A[n][m]*B[m][p];
			}
		}
	}
}

/************************************************************************
 ScalarMatMult:
	Multiplies NxM matrix A times scalara multiplier. Returns C (NxM)
	matrix C may be same as matrix A
-----------------------------------------------------------------------*/
void   ScalarMatMult(Ironclad2DArray  A,const double& muliplier,const int	N,const int M,Writeable2DArray C)
{
	for(int i=0;i<N;i++) {
		for(int j=0; j<M;j++) {
			C[i][j]=A[i][j]*muliplier;
		}
	}
}
/************************************************************************
 MatAdd:
	Adds NxM matrix A to NxM matrix B. Returns C (NxM)
	matrix C may be same as matrix A or B
-----------------------------------------------------------------------*/
void   MatAdd(Ironclad2DArray  A,Ironclad2DArray  B,const int	N,const int M,Writeable2DArray C)
{
	for(int i=0;i<N;i++) {
		for(int j=0; j<M;j++) {
			C[i][j]=A[i][j]+B[i][j];
		}
	}
}
/************************************************************************
 CopyMatrix:
	Copies NxM matrix A to NxM matrix C
-----------------------------------------------------------------------*/
void   CopyMatrix(Ironclad2DArray  A,const int	N,const int M,Writeable2DArray C)
{
	for(int i=0;i<N;i++) {
		for(int j=0; j<M;j++) {
			C[i][j]=A[i][j];
		}
	}
}
/************************************************************************
 TransposeMat:
	Transposes NxM matrix A to MxN matrix AT
-----------------------------------------------------------------------*/
void TransposeMat(Ironclad2DArray  A,const int N,const int M,Writeable2DArray AT)
{
	for(int i=0;i<N;i++) {
		for(int j=0; j<M;j++) {
			AT[j][i]=A[i][j];
		}
	}
}
/************************************************************************
 AllocateMatrix:
	allocates memory to NxM matrix A , initializes to all zeroes
-----------------------------------------------------------------------*/
void AllocateMatrix(const int N,const int M,double**&A) {
	A=NULL;
	A=new double *[N];
	ExitGracefullyIf(A==NULL,"Matrix allocation (1)",OUT_OF_MEMORY);
	for(int i=0;i<N;i++) {
		A[i]=NULL;
		A[i]=new double [M];
		ExitGracefullyIf(A[i]==NULL,"Matrix allocation",OUT_OF_MEMORY);
		for(int j=0; j<M;j++) {
			A[i][j]=0.0;
		}
	}
}
/************************************************************************
 AllocateMatrix:
	deletes 2D NxM matrix A
-----------------------------------------------------------------------*/
void DeleteMatrix(const int N,const int M,double** A) {
	for(int i=0;i<N;i++) {
		delete [] A[i]; A[i]=NULL;
	}
	delete [] A; A=NULL;
}
inline double SIGN(const double& a,const double& b) {
	return b>=0.0 ? (a>=0.0 ? a : -a): (a>=0.0 ? -a : a);
}
//------------------------------------------------------------------------------
inline double    pythag(const double a,const double b) {
	//not subject to overflow errs c=sqrt(a^2+b^2)
	double absa=fabs(a);
	double absb=fabs(b);
	if(absa>absb) { return absa*sqrt(1.0+(absb/absa)*(absb/absa)); }
	else { return (absb==0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb))); }
}
/************************************************************************
 Singular Value Decomposition:
Returns solution, x and rank
	this operation destroys the A matrix
	from Numerical Recipes (Press et al)
-----------------------------------------------------------------------*/
bool SVD(Ironclad2DArray AA,
				 Ironclad1DArray b,
				 Writeable1DArray x,
				 const int size,
				 const double SVTolerance)
{
	bool flag;
	int i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,xx,y,z;

	int m=size;
	int n=size;
	double* tmp=new double[n];
	double* w  =new double[n];
	double** v=NULL;
	AllocateMatrix(n,n,v);
	double **A=NULL;
	AllocateMatrix(n,n,A);
	CopyMatrix(AA,n,n,A);
	/*	cout <<endl;
		for (j=0;j<size;j++){
			for (i=0;i<size;i++){
				cout << A[i][j] <<" ";
			}
			cout << "| "<<b[j]<<" | "<<endl;
		}
		cout <<endl;*/

	g=scale=anorm=0.0;
	for(i=0; i<n; i++) {
		l=i+2;
		tmp[i]=scale*g;
		g=s=scale=0.0;
		if(i<m) {
			for(k=i;k<m;k++) { scale+=fabs(A[k][i]); }
			if(scale!=0.0) {
				for(k=i; k<m; k++) {
					A[k][i]/=scale;
					s+=A[k][i]*A[k][i];
				}
				f=A[i][i];
				g=-SIGN(sqrt(s),f);
				h=f*g-s;
				A[i][i]=f-g;
				for(j=l-1;j<n;j++) {
					s=0.0;
					for(k=i;k<m;k++) { s+=A[k][i]*A[k][j]; }
					f=s/h;
					for(k=i;k<m;k++) { A[k][j]+=f*A[k][i]; }
				}
				for(k=i;k<m;k++) { A[k][i]*=scale; }
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if((i+1<=m) && (i!=n)) {
			for(k=l-1;k<n;k++) { scale+=fabs(A[i][k]); }
			if(scale!=0.0) {
				for(k=l-1;k<n;k++) {
					A[i][k]/=scale;
					s+=A[i][k]*A[i][k];
				}
				f=A[i][l-1];
				g=-SIGN(sqrt(s),f);
				h=f*g-s;
				A[i][l-1]=f-g;
				for(k=l-1;k<n;k++) { tmp[k]=A[i][k]/h; }
				for(j=l-1;j<m;j++) {
					s=0.0;
					for(k=l-1;k<n;k++) { s      +=A[j][k]*A[i][k]; }
					for(k=l-1;k<n;k++) { A[j][k]+=s*tmp[k]; }
				}
				for(k=l-1;k<n;k++) { A[i][k]*=scale; }
			}
		}
		upperswap(anorm,(fabs(w[i])+fabs(tmp[i])));
	}
	for(i=n-1; i>=0; i--) {
		if(i<n-1) {
			if(g!=0.0) {
				for(j=l;j<n;j++) { v[j][i]=(A[i][j]/A[i][l])/g; }
				for(j=l;j<n;j++) {
					s=0.0;
					for(k=l;k<n;k++) { s 		 +=A[i][k]*v[k][j]; }
					for(k=l;k<n;k++) { v[k][j]+=s*v[k][i]; }
				}
			}
			for(j=l;j<n;j++) { v[i][j]=v[j][i]=0.0; }
		}
		v[i][i]=1.0;
		g=tmp[i];
		l=i;
	}
	for(i=min(n,m)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for(j=l;j<n;j++) { A[i][j]=0.0; }
		if(g!=0.0) {
			g=1.0/g;
			for(j=l;j<n;j++) {
				s=0.0;
				for(k=l;k<m;k++) { s+=A[k][i]*A[k][j]; }
				f=(s/A[i][i])*g;
				for(k=i;k<m;k++) { A[k][j]+=f*A[k][i]; }
			}
			for(j=i;j<m;j++) { A[j][i]*=g; }
		}
		else {
			for(j=i;j<m;j++) { A[j][i]=0.0; }
		}
		++A[i][i];
		//A[i][i]++;
	}
	for(k=n-1;k>=0;k--) {
		for(its=0; its<30;its++) {
			//cout <<"SVD iter "<< its <<endl;
			flag=true;
			for(l=k;l>=0;l--) {
				nm=l-1;
				if(fabs(tmp[l])+anorm==anorm) {
					flag=false;
					break;
				}
				if(fabs(w[nm])+anorm==anorm) { break; }
			}//end for l
			if(flag) {
				c=0.0;
				s=1.0;
				for(i=l-1;i<k+1;i++) {
					f=s*tmp[i];
					tmp[i]*=c;
					if(fabs(f)+anorm==anorm) { break; }
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s=-f*h;
					for(j=0;j<m;j++) {
						y=A[j][nm];
						z=A[j][i];
						A[j][nm]=y*c+z*s;
						A[j][i]=z*c-y*s;
					}
				}
			}//end if flag
			z=w[k];
			if(l==k) {
				if(z<0.0) {
					w[k]=-z;
					for(j=0;j<n;j++) { v[j][k]=-v[j][k]; }
				}
				break;
			}
			if(its==29) { ExitGracefully("Singular Value Decomposition did not converge",RUNTIME_ERR); }
			xx=w[l];
			nm=k-1;
			y=w[nm];
			g=tmp[nm];
			h=tmp[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((xx-z)*(xx+z)+h*((y/(f+SIGN(g,f)))-h))/xx;
			c=s=1.0;
			for(j=l; j<=nm; j++) {
				i=j+1;
				g=tmp[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				tmp[j]=z;
				c=f/z;
				s=h/z;
				f=xx*c+g*s;
				g=g*c-xx*s;
				h=y*s;
				y*=c;
				for(jj=0;jj<n;jj++) {
					xx=v[jj][j];
					z=v[jj][i];
					v[jj][j]=xx*c+z*s;
					v[jj][i]=z*c-xx*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if(z!=0.0) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				xx=c*y-s*g;
				for(jj=0;jj<m;jj++) {
					y=A[jj][j];
					z=A[jj][i];
					A[jj][j]=y*c+z*s;
					A[jj][i]=z*c-y*s;
				}
			}//end for j=l to nm
			tmp[l]=0.0;
			tmp[k]=f;
			w[k]=xx;
		}//end for its
	}
	//end sub svdcmp (Press et al)

	double maxw(-ALMOST_INF);
	double minw(ALMOST_INF);

	for(j=0; j<n; j++) { upperswap(maxw,w[j]); }                //get condition number info
	for(j=0; j<n; j++) {
		if(fabs(w[j]/maxw)<SVTolerance) { w[j]=0.0; }             //zero out bad terms
	}
	for(j=0; j<n; j++) { if(w[j]!=0) { lowerswap(minw,w[j]); } } //reevaluate min

	//back substitute
	for(j=0; j<n; j++) {
		s=0.0;
		if(w[j]!=0.0) {
			for(i=0; i<m; i++) { s+=A[i][j]*b[i]; }
			s/=w[j];
		}
		tmp[j]=s;
	}

	for(j=0; j<n; j++) {
		s=0.0;
		for(jj=0;jj<n;jj++) { s+=v[j][jj]*tmp[jj]; }
		x[j]=s;
	}

	/*	for (j=0;j<size;j++){
			cout << "| "<<x[j]<<" | "<< w[j]<<endl;
		}*/

		//calculate error

		//cout << "n: "<<n<< " m: "<<m<<endl;
	delete[] tmp;
	delete[] w;
	DeleteMatrix(n,n,v);
	DeleteMatrix(n,n,A);

	//cout <<"SVD CONDITION #: " <<maxw/minw<<endl;
	return true; //should be dependent upon condition number (max wj/minwj) should be low
}
