#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<iomanip>
#include<iostream>
#include<algorithm>
using namespace std;
//#define NORM_2
#define NORM_INF
const int N=501,N_MAX=501+10,inf=~0u>>1;
const double eps=1e-15,b=0.16,c=-0.064;

double lambda[N_MAX];

int sgn(double x){
	if (x>eps) return 1;
	if (x<-eps) return -1;
	return 0;
}
double sqr(double x){ return x*x; }
double Matrix(int i,int j){
	if (i==j) return (1.64-0.024*i)*sin(0.2*i)-0.64*exp(0.1/i);
	else if ((i==j+1) || (i==j-1)) return b;
	else if ((i==j+2) || (i==j-2)) return c;
	return 0;
}
#ifdef NORM_2
double PowerMethod(){
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
			for (int j=1; j<=N; ++j)
				u[i]+=Matrix(i,j)*y[j];
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
double PowerMethod(){
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
			for (int j=1; j<=N; ++j)
				u[i]+=Matrix(i,j)*y[j];
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
int main(){
	printf("%12E\n", PowerMethod());
	cout<<setprecision(12)<<PowerMethod()<<endl;
	return 0;
}
