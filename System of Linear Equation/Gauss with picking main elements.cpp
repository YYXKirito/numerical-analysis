#include"../lib/matrix.h"
#include<cmath>
#include<cstdio>
#include<cstring>
#include<iostream>
using namespace std;
const int N_SIZE=100+10;
const double eps=1e-5;
Matrix<double> A;int N;
double sqr(double x){ return x*x; }
void MakeEquation(){}
void Gauss(){
	for (int i=1; i<=N-1; ++i){
		if (abs(A[i][i])<eps) break;
		int L=i;
		for (int j=i+1; j<=N; ++j)
			if (fabs(A[j][i])>fabs(A[L][i])) L=j;
		A.swapline(i,L);
		for (int j=i+1; j<=N; ++j){
			double tmp=A[j][i]/A[i][i];
			for (int k=i; k<=N+1; ++k){
				A[j][k]-=A[i][k]*tmp;
			}
		}
	}
	A[N][N+1]/=A[N][N];
	for (int i=N-1; i>=1; --i){
		for (int j=i+1; j<=N; ++j)
			A[i][N+1]-=A[j][N+1]*A[i][j];
		A[i][N+1]/=A[i][i];
	}
}
int main(){
	//input
	MakeEquation();
	Gauss();
	//output
	return 0;
}
