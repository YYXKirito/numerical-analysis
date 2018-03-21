#include"../lib/matrix.h"
#include<cmath>
#include<cstdio>
#include<cstring>
#include<iostream>
using namespace std;
const int N_SIZE=100+10;
const double eps=1e-5;
Matrix<double> A;int N;double p[N_SIZE][N_SIZE];
double sqr(double x){ return x*x; }
void MakeEquation(){
	A=Matrix<double>(N,N+1);
	for (int i=1; i<=N; ++i){
		for (int j=1; j<=N; ++j){
			A[i][N+1]+=sqr(p[N+1][j])-sqr(p[i][j]);
			A[i][j]=p[N+1][j]-p[i][j];
		}
		A[i][N+1]/=2;
	}
}
void Gauss(){
	for (int i=1; i<=N-1; ++i){
		if (abs(A[i][i])<eps) break;
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
