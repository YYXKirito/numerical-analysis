#include"../lib/matrix.h"
#include<cmath>
#include<cstdio>
#include<cstring>
using namespace std;
const int N_MAX=100+10;
Matrix<double> A;int N;
void MakeEquation(){}
void Crout(){
	for (int i=1; i<=N; ++i){
		for (int j=i; j<=N; ++j){
			double tmp=A[j][i];
			for (int k=1; k<=i-1; ++k)
				tmp-=A[j][k]*A[k][i];
			A[j][i]=tmp;
		}
		for (int j=i+1; j<=N; ++j){
			double tmp=A[i][j];
			for (int k=1; k<=i-1; ++k)
				tmp-=A[i][k]*A[k][j];
			A[i][j]=tmp/A[i][i];
		}
	}
	A[1][N+1]/=A[1][1];
	for (int i=2; i<=N; ++i){
		double tmp=A[i][N+1];
		for (int j=1; j<=i-1; ++j)
			tmp-=A[i][j]*A[j][N+1];
		A[i][N+1]=tmp/A[i][i];
	}
	for (int i=N-1; i>=1; --i){
		double tmp=A[i][N+1];
		for (int j=i+1; j<=N; ++j)
			tmp-=A[i][j]*A[j][N+1];
		A[i][N+1]=tmp;
	}
}
int main(){
	//input
    MakeEquation();
	Crout();
	//output
	return 0;
}
