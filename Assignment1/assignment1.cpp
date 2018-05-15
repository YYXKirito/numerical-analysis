#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<iomanip>
#include<iostream>
#include<algorithm>
using namespace std;
#define NORM_2
//#define NORM_INF
const int N=501,N_MAX=501+10,inf=~0u>>1,r=2,s=2;
const double eps=1e-12,b=0.16,c=-0.064;

double lambda[N_MAX],C[r+s+1+10][N_MAX],Doo[r+s+1+10][N_MAX],lambda1,lambda501,lambda_s,det,cond;

inline int sgn(double x){
	if (x>eps) return 1;
	if (x<-eps) return -1;
	return 0;
}
inline double sqr(double x){ return x*x; }
inline double Matrix(int i,int j,double p=0){
	if (i>j+2 || i<j-2) return 0;
	if (i==j) return C[i-j+s+1][j]-p;
	return C[i-j+s+1][j];
}

void Init(){
	memset(C,0,sizeof(C));
	for (int i=1; i<=N; ++i)
		C[s+1][i]=(1.64-0.024*i)*sin(0.2*i)-0.64*exp(0.1/i);
	for (int i=2; i<=N; ++i)
		C[s+2][i-1]=b;
	for (int i=2; i<=N; ++i)
		C[s+3][i-2]=c;
	for (int i=1; i<=N-1; ++i)
		C[s][i+1]=b;
	for (int i=1; i<=N-1; ++i)
		C[s-1][i+2]=c;
}

#ifdef NORM_2
double PowerMethod(double p=0){
	double u[N_MAX],y[N_MAX],eta,beta=0;
	for (int i=1; i<=N; ++i)
		u[i]=1;
	for (; ; ){
		eta=0;
		for (int i=1; i<=N; ++i)
			eta+=sqr(u[i]);
		eta=sqrt(eta);
		for (int i=1; i<=N; ++i)
			y[i]=u[i]/eta;
		for (int i=1; i<=N; ++i){
			u[i]=0;
			for (int j=1; j<=N; ++j){
				if (j>i+2) break;
				u[i]+=Matrix(i,j,p)*y[j];
			}
		}
		double cur=0;
		for (int i=1; i<=N; ++i)
			cur+=y[i]*u[i];
		if (fabs(cur-beta)/fabs(cur)<=eps) return cur;
		beta=cur;
	}
}
#endif

#ifdef NORM_INF
double PowerMethod(double p=0){
	double u[N_MAX],y[N_MAX],beta=0;
	for (int i=1; i<=N; ++i)
		u[i]=1;
	for (; ; ){
		double hr=0;
		for (int i=1; i<=N; ++i)
			if (fabs(u[i])-fabs(hr)>eps)
				hr=u[i];
		for (int i=1; i<=N; ++i)
			y[i]=u[i]/fabs(hr);
		for (int i=1; i<=N; ++i){
			u[i]=0;
			for (int j=1; j<=N; ++j){
				if (j>i+2)) break;
				u[i]+=Matrix(i,j,p)*y[j];
			}
		}
		double cur=0;
		for (int i=1; i<=N; ++i)
			if (fabs(u[i])-fabs(cur)>eps)
				cur=u[i];
		cur=cur*sgn(hr);
		if (fabs(cur-beta)/fabs(cur)<=eps) return cur;
		beta=cur;
	}
}
#endif
inline int convert(int i,int j){
	return ((i-j+s+1>0)&&(i-j+s+1<6))?i-j+s+1:0;
}
void DoolittleInit(double p=0){
	memcpy(Doo,C,sizeof(C));
	for (int i=1; i<=N; ++i)
		Doo[convert(i,i)][i]-=p;
	for (int i=1; i<=N; ++i){
		for (int j=i; j<=min(N,i+2); ++j){
			double tmp=Doo[convert(i,j)][j];
			for (int k=1; k<=i-1; ++k)
				tmp-=Doo[convert(i,k)][k]*Doo[convert(k,j)][j];
			Doo[convert(i,j)][j]=tmp;
		}
		for (int j=i+1; j<=N; ++j){
			double tmp=Doo[convert(j,i)][i];
			for (int k=max(1,i-2); k<=i-1; ++k)
				tmp-=Doo[convert(j,k)][k]*Doo[convert(k,i)][i];
			Doo[convert(j,i)][i]=tmp/Doo[convert(i,i)][i];
		}
	}
}
void Doolittle(double u[],double y[],double p=0){
	double Y[N_MAX];
	memcpy(Y,y,sizeof(Y));
	for (int i=2; i<=N; ++i){
		double tmp=y[i];
		for (int j=max(1,i-2); j<=i-1; ++j)
			tmp-=Doo[convert(i,j)][j]*Y[j];
		Y[i]=tmp;
	}
	u[N]=Y[N]/Doo[convert(N,N)][N];
	for (int i=N-1; i>=1; --i){
		double tmp=Y[i];
		for (int j=i+1; j<=min(N,i+2); ++j)
			tmp-=(Doo[convert(i,j)][j])*u[j];
		u[i]=tmp/Doo[convert(i,i)][i];
	}
}

double InvPowerMethod(double p=0){
	double u[N_MAX],y[N_MAX],beta=inf;
	for (int i=1; i<=N; ++i)
		u[i]=1;
	for (; ; ){
		double eta=0;
		for (int i=1; i<=N; ++i)
			eta+=sqr(u[i]);
		eta=sqrt(eta);
		for (int i=1; i<=N; ++i)
			y[i]=u[i]/eta;
		Doolittle(u,y,p);
		double cur=0;
		for (int i=1; i<=N; ++i)
			cur+=y[i]*u[i];
		if (fabs(1/cur-1/beta)/fabs(1/cur)<=eps) return 1/cur;
		beta=cur;
	}
}
double Det(){
	double ret=1;
	for (int i=1; i<=N; ++i)
		ret*=Doo[convert(i,i)][i];
	return ret;
}
int main(){
	Init();
	DoolittleInit();
	det=Det();
	lambda1=PowerMethod();
	lambda_s=InvPowerMethod();
	lambda501=PowerMethod(lambda1)+lambda1;
	cond=fabs(lambda1/lambda_s);
	if (lambda1>lambda501) swap(lambda1,lambda501);
	printf("lambda_1=%.12E lambda_501=%.12E lambda_s=%.12E\n", lambda1, lambda501, lambda_s);
	for (int i=1; i<=39; ++i){
		double miu=lambda1+i*(lambda501-lambda1)/40;
		DoolittleInit(miu);
		printf("lambda_i%d=%.12E ", i, InvPowerMethod(miu)+miu);
	}
	printf("\ncond(A)=%.12E detA=%.12E\n", cond, det);
	
	return 0;
}
