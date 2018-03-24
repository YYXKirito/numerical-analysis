#include"../lib/matrix.h"
#include<cmath>
#include<cstdio>
#include<cstring>
using namespace std;
const double eps=1e-5;
Matrix<double> A;int N;double s[N_MAX];
void MakeEquation(){}
void Doolittle(){
	for (int i=1; i<=N; ++i){
		int L=i;
		for (int j=i; j<=N; ++j){
			s[j]=A[j][i];
			for (int k=1; k<=i-1; ++k)
				s[j]-=A[j][k]*A[k][i];
			if (fabs(s[j])>fabs(s[L])+eps) L=j;
		}
		A.swapline(L,i,1,N+1); swap(s[i],s[L]);
		A[i][i]=s[i];
		for (int j=i+1; j<=N; ++j){
			double tmp=A[i][j];
			for (int k=1; k<=i-1; ++k)
				tmp-=A[i][k]*A[k][j];
			A[i][j]=tmp;
		}
		for (int j=i+1; j<=N; ++j)
			A[j][i]=s[j]/A[i][i];
	}
	for (int i=2; i<=N; ++i){
		double tmp=A[i][N+1];
		for (int j=1; j<=i-1; ++j)
			tmp-=A[i][j]*A[j][N+1];
		A[i][N+1]=tmp;
	}
	A[N][N+1]/=A[N][N];
	for (int i=N-1; i>=1; --i){
		double tmp=A[i][N+1];
		for (int j=i+1; j<=N; ++j)
			tmp-=A[i][j]*A[j][N+1];
		A[i][N+1]=tmp/A[i][i];
	}
}
int main(){
	//input
    MakeEquation();
	Doolittle();
	//output
	return 0;
}
