#include<cstdio>
#include<cstring>
#include<iostream>
using namespace std;
template<typename T>class Matrix{
	static const int MAX_SIZE=1000+10;
	public:
	int N,M;T bit[MAX_SIZE][MAX_SIZE];
	Matrix(){ N=M=0;memset(bit,0,sizeof(bit)); }
	Matrix(int _N,int _M){ N=_N;M=_M;memset(bit,0,sizeof(bit)); }
	
	T* operator[](const int &p){ return bit[p]; }
	Matrix operator +(Matrix p)const{
		Matrix result(N,M);
		for (int i=1; i<=N; ++i)
			for (int j=1; j<=M; ++j)
				result[i][j]=bit[i][j]+p[i][j];
		return result;
	}
	Matrix operator -(Matrix p)const{
		Matrix result(N,M);
		for (int i=1; i<=N; ++i)
			for (int j=1; j<=M; ++j)
				result[i][j]=bit[i][j]-p[i][j];
		return result;
	}
	Matrix operator *(Matrix p)const{
		Matrix result(N,p.M);
		for (int i=1; i<=N; ++i)
			for (int j=1; j<=p.M; ++j)
				for (int k=1; k<=M; ++k)
					result[i][j]+=bit[i][k]*p[k][j];
		return result;
	}
	
	Matrix operator +=(Matrix p){
		*this=*this+p;
		return *this;
	}
	Matrix operator -=(Matrix p){
		*this=*this-p;
		return *this;
	}
	Matrix operator *=(Matrix p){
		*this=*this*p;
		return *this;
	}
	
	void swapline(int i,int j,int st=1,int ed=MAX_SIZE){
		T t;
		for (int k=st; k<=ed; ++k){
			t=bit[i][k]; bit[i][k]=bit[j][k]; bit[j][k]=t;
		}
	}
	void swapcol(int i,int j,int st=1,int ed=MAX_SIZE){
		T t;
		for (int k=st; k<=ed; ++k){
			t=bit[k][i]; bit[k][i]=bit[k][j]; bit[k][j]=t;
		}
	}
	void print(){
		printf("N=%d	M=%d\n", N, M);
		for (int i=1; i<=N; ++i){
			for (int j=1; j<=M; ++j)
				cout<<bit[i][j]<<"\t";
			cout<<endl;
		}
	}
};
