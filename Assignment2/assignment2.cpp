#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<iomanip>
#include<iostream>
#include"../lib/matrix.h"
using namespace std;
const double eps=1e-12;
const int N=10,L=100;
Matrix<double> A(10,10),P(10,10),Q(10,10),T(10,10),AAA(10,10);
void Init(){
	for (int i=1; i<=N; ++i)
		for (int j=1; j<=N; ++j)
			if (i==j)
				A[i][j]=1.52*cos(i+1.2*j);
			else
				A[i][j]=sin(0.5*i+0.2*j);
}
inline int min(int a,int b){ return a<b?a:b; }
inline double sqr(double x){ return x*x; }
inline double sign(double x){
	if (fabs(x)<eps) return 0;
	return x/fabs(x);
}

namespace QuasiTriangulation{
	bool couldSkip(int lineIndex){
		for (int i=lineIndex+2; i<=N; ++i)
			if (fabs(A[i][lineIndex])>eps)
				return 0;
		return 1;
	}
	void main(){
		for (int lineIndex=1; lineIndex<=N-2; ++lineIndex){
			if (couldSkip(lineIndex)) continue;
			double d=0;
			for (int i=lineIndex+1; i<=N; ++i)
				d+=sqr(A[i][lineIndex]);
			double c=-sign(A[lineIndex+1][lineIndex])*sqrt(d);
			if (fabs(A[lineIndex+1][lineIndex])<eps) c=sqrt(d);
			double u[N+1],h=sqr(c)-c*A[lineIndex+1][lineIndex];
			for (int i=1; i<=N; ++i)
				if (i<=lineIndex)
					u[i]=0;
				else
					u[i]=A[i][lineIndex]-c*(i==lineIndex+1);
			double p[N+1],q[N+1],w[N+1],t=0;
			T=A.Transposition();
			for (int i=1; i<=N; ++i){
				p[i]=q[i]=w[i]=0;
				for (int j=1; j<=N; ++j){
					p[i]+=T[i][j]*u[j];
					q[i]+=A[i][j]*u[j];
				}
				p[i]/=h; q[i]/=h;
				t+=p[i]*u[i];
			}
			t/=h;
			for (int i=1; i<=N; ++i)
				w[i]=q[i]-t*u[i];
			for (int i=1; i<=N; ++i)
				for (int j=1; j<=N; ++j){
					P[i][j]=w[i]*u[j];
					Q[i][j]=u[i]*p[j];
				}
			A=A-P-Q;
		}
		for (int i=3; i<=N; ++i)
			for (int j=1; j<i-1; ++j)
				A[i][j]=0;
		A.print(); AAA=A;
	}
}
#define Re(x) lambda[x][0]
#define Im(x) lambda[x][1]
double lambda[N+1][2];int cnt=0;
namespace QRMethod{
	void calculate(double delta,double a,double b){
		if (delta>eps){
			++cnt;
			Re(cnt)=(-b+sqrt(delta))/(2*a);
			Im(cnt)=0;
			++cnt;
			Re(cnt)=(-b-sqrt(delta))/(2*a);
			Im(cnt)=0;
			return;
		}
		++cnt;
		Re(cnt)=-b/(2*a);
		Im(cnt)=sqrt(-delta)/(2*a);
		++cnt;
		Re(cnt)=-b/(2*a);
		Im(cnt)=-sqrt(-delta)/(2*a);
	}
	Matrix<double> M(N,N),Q(N,N),R(N,N),I(N,N),C(N,N),B(N,N);
	double u[N+1],v[N+1],p[N+1],q[N+1],w[N+1];
	void main1(){
		for (int k=1,m=N; ; ++k){
			if (!m) break;
			if (m==1){
				++cnt;
				Re(cnt)=A[m][m];
				Im(cnt)=0;
				break;
			}
			if (fabs(A[m][m-1])<=eps){
				++cnt;
				Re(cnt)=A[m][m];
				Im(cnt)=0;
				--m;
				continue;
			}
			double a=1,b=-(A[m-1][m-1]+A[m][m]),c=A[m-1][m-1]*A[m][m]-A[m][m-1]*A[m-1][m];
			double delta=sqr(b)-4*a*c;
			if (fabs(A[m-1][m-2])<=eps){
				calculate(delta,a,b);
				m-=2;
				continue;
			}
			if (k==L) break;
			double s=A[m-1][m-1]+A[m][m],t=c;
			for (int i=1; i<=N; ++i)
				I[i][i]= i<=m?1:0;
			M=A*A-A*s+I*t;
			for (int i=1; i<=m; ++i)
				for (int j=1; j<=m; ++j)
					B[i][j]=M[i][j];
			C=A;
			for (int r=1; r<=m-1; ++r){
				int flag=1;
				for (int i=r+1; i<=m; ++i)
					if (fabs(B[i][r])>eps){
						flag=0;
						break;
					}
				if (flag) continue;
				double sigma=0;
				for (int i=r; i<=m; ++i)
					sigma+=sqr(B[i][r]);
				double d=sqrt(sigma),c,h;
				if (fabs(B[r][r])<=eps) c=d; else c=-sign(B[r][r])*d;
				h=sqr(c)-c*B[r][r];
				for (int i=1; i<=m; ++i)
					u[i]=i<r?0:B[i][r]-c*(i==r);
				for (int i=1; i<=m; ++i){
					v[i]=0;
					for (int j=1; j<=m; ++j)
						v[i]+=u[j]*B[j][i];
					v[i]/=h;
				}
				for (int i=1; i<=m; ++i)
					for (int j=1; j<=m; ++j)
						B[i][j]-=u[i]*v[j];
				for (int i=1; i<=m; ++i){
					p[i]=q[i]=0;
					for (int j=1; j<=m; ++j){
						p[i]+=u[j]*C[j][i];
						q[i]+=u[j]*C[i][j];
					}
					p[i]/=h; q[i]/=h;
				}
				double t=0;
				for (int i=1; i<=m; ++i)
					t+=p[i]*u[i];
				t/=h;
				for (int i=1; i<=m; ++i)
					w[i]=q[i]-t*u[i];
				for (int i=1; i<=m; ++i)
					for (int j=1; j<=m; ++j)
						C[i][j]-=w[i]*u[j]+u[i]*p[j];
			}
			A=C;
		}
		for (int i=cnt; i; --i)
			cout<<setprecision(12)<<"Characteristic value No. "<<cnt-i+1<<":  Re  "<<Re(i)<<"  Im  "<<Im(i)<<endl;
	}
	
	bool couldSkip(int lineIndex){
		for (int i=lineIndex+1; i<=N; ++i)
			if (fabs(A[i][lineIndex])>eps)
				return 0;
		return 1;
	}
	void main2(){
		A=AAA;Q=Matrix<double> (N,N);
		for (int i=1; i<=N; ++i)
			Q[i][i]=1;
		for (int lineIndex=1; lineIndex<=N-1; ++lineIndex){
			if (couldSkip(lineIndex)) continue;
			double d=0;
			for (int i=lineIndex; i<=N; ++i)
				d+=sqr(A[i][lineIndex]);
			d=sqrt(d);
			double c;
			if (A[lineIndex][lineIndex]<-eps) c=d; else c=-d;
			double u[N+1],h=sqr(c)-c*A[lineIndex][lineIndex];
			for (int i=1; i<=N; ++i)
				if (i<=lineIndex-1)
					u[i]=0;
				else
					u[i]=A[i][lineIndex]-c*(i==lineIndex);
			double p[N+1],q[N+1],w[N+1];
			for (int i=1; i<=N; ++i){
				w[i]=0;
				for (int j=1; j<=N; ++j)
					w[i]+=Q[i][j]*u[j];
			}
			for (int i=1; i<=N; ++i)
				for (int j=1; j<=N; ++j)
					Q[i][j]-=w[i]*u[j]/h;
			for (int i=1; i<=N; ++i){
				p[i]=0;
				for (int j=1; j<=N; ++j)
					p[i]+=A[j][i]*u[j];
				p[i]/=h;
			}
			for (int i=1; i<=N; ++i)
				for (int j=1; j<=N; ++j)
					A[i][j]-=u[i]*p[j];
		}
		cout<<"----------------------------------------R---------------------------------------------"<<endl;
		A.print();
		cout<<"----------------------------------------Q---------------------------------------------"<<endl;
		Q.print();
		cout<<"---------------------------------------R*Q---------------------------------------------"<<endl;
		T=A*Q; T.print();
		cout<<"---------------------------------------Q*R---------------------------------------------"<<endl;
		T=Q*A; T.print();
	}
}
namespace GaussMethod{
	Matrix<double> AA(N,N+1);
	double x[N+1];
	void MakeEquation(double val){
		AA=A; AA.M=N+1;
		for (int i=1; i<=N; ++i){
			AA[i][i]-=val;
			AA[i][N+1]=0;
		}
	}
	void Gauss(){
		for (int i=1; i<=N-1; ++i){
			if (fabs(AA[i][i])<eps) break;
			int L=i;
			for (int j=i+1; j<=N; ++j)
				if (fabs(AA[j][i])>fabs(AA[L][i])) L=j;
			AA.swapline(i,L);
			for (int j=i+1; j<=N; ++j){
				double tmp=AA[j][i]/AA[i][i];
				for (int k=i; k<=N+1; ++k){
					AA[j][k]-=AA[i][k]*tmp;
				}
			}
		}
		for (int i=1; i<=N-1; ++i) x[i]=0;
		x[N]=1; AA[N][N+1]/=AA[N][N];
		for (int i=N-1; i>=1; --i){
			for (int j=i+1; j<=N; ++j)
				x[i]-=x[j]*AA[i][j];
			x[i]=(x[i]+AA[i][N+1])/AA[i][i];
		}
		double sigma=0;
		for (int i=1; i<=N; ++i)
			sigma+=sqr(x[i]);
		sigma=sqrt(sigma);
		for (int i=1; i<=N; ++i)
			x[i]/=sigma;
	}
	int main(){
		Init();
		for (int i=cnt; i>=1; --i){
			if (fabs(Im(i))>eps) continue;
			MakeEquation(Re(i));
			Gauss();
			cout<<"Characteristic vector No. "<<cnt-i+1<<":"<<endl;
			cout<<"[";
			for (int j=1; j<=N-1; ++j)
				cout<<setprecision(12)<<x[j]<<", ";
			cout<<setprecision(12)<<x[N]<<"]"<<endl;
		}
		return 0;
	}
}
int main(){
	Init();
	QuasiTriangulation::main();
	QRMethod::main1();
	QRMethod::main2();
	GaussMethod::main();//Use Gauss'Method to solve for character vector
	return 0;
}
