<center><font size=6>数值分析第一次大作业实验报告</font></center>
<center><font size=5>16091038 叶宇轩</font></center>
## 题目描述
&emsp;&emsp;<font size=4>给定一个501阶矩阵$A$，  
$$
\left[
\begin{matrix}
a_1 & b      & c      &        &   \\
b   & \ddots & \ddots & \ddots &   \\
c   & \ddots & \ddots & \ddots & c \\
    & \ddots & \ddots & \ddots & b \\
    &        & c      & b      & a_{501}
\end{matrix}
\right]
$$
其中，$a_i=(1.64-0.024i)\sin(0.2i)-0.64e^{0.1\over i}\ \ (i=1,2,\cdots,501)，b=0.16，c=-0.064。$记矩阵$A$的特征值为$\lambda_i\ \ (i=1,2,\cdots,501)$，并且有
<center>$\lambda_1\leq\lambda_2\leq...\leq\lambda_{501}，\ |\lambda_s|\ =\ min \{ |\lambda_i| \ :\ 1\leq i\leq 501 \}$ </center>  
1.求$\lambda_1，\lambda_{501}，\lambda_s$  
2.求$A$与$\mu_k=\lambda_1+k\frac{\lambda_{501} - \lambda_1}{40}$最接近的特征值  
3.求$A$的谱范数条件数 $cond(A)_2$和行列式$detA$ 
</font>  
___
## 要求
<font size=4>
&emsp;&emsp;1. 算法涉及的精度取$\epsilon=10^{-12}$  
&emsp;&emsp;2. $A$的零元素不储存  
&emsp;&emsp;3. 显示至少12位有效数字
</font>
___
## 算法设计
### 1.需要使用的算法
#### 1.1 幂法
<font size=4>
&emsp;&emsp;幂法主要用于计算矩阵按模最大的特征向量。  
&emsp;&emsp;设$n\times n$实矩阵$A$有$n$个线性无关的特征向量$x_1,x_2,\dots,x_n$，$x_i$对应的特征值$\lambda_i$满足
$$|\lambda_1|>|\lambda_2|\geq|\lambda_3|\geq\cdots\geq|\lambda_n|，Ax_i=\lambda_ix_i\ \ (i=1,2,\cdots,n)$$
&emsp;&emsp;任取$n$维非零向量$u_0$，有
$$u_0=\alpha_1x_1+\alpha_2x_2+\cdots+\alpha_nx_n$$
&emsp;&emsp;从$u_0$出发，令$u_k=A u_{k-1}\ \ (k=1,2,\cdots)$，可以推出
$$u_k=A u_{k-1}=\cdots=A^ku_0=$$  
$$\alpha_1 A^k x_1 + \alpha_2 A^k x_2 + \cdots + \alpha_n A^k x_n=$$  
$$\alpha_1 \lambda_1^k x_1 + \alpha_2 \lambda_2^k x_2 + \cdots + \alpha_n \lambda_n^k x_n=$$  
$$\lambda_1^k[\alpha_1 x_1 + \alpha_2 (\frac{\lambda_2}{\lambda_1})^k x_2 + \cdots + \alpha_n (\frac{\lambda_n}{\lambda_1})^k x_n]$$  
&emsp;&emsp;设$\alpha_1\neq 0$，则k充分大时，有
$$u_k \approx \lambda_1^k \alpha_1 x1$$  
&emsp;&emsp;实际计算中为了避免$u_k$的模过大，每次迭代后对$u_k$单位化，因此迭代公式为　　
$$
\begin{cases}
y_{k-1} = \frac{u_{k-1}}{\parallel u_{k-1} \parallel}\\
u_k = A y_{k-1}
\end{cases}
\ \ (k=1,2,\cdots)
$$
&emsp;&emsp;当使用二范数$\parallel \cdot \parallel_2$时，令$\beta_k = y_{k-1}^T u_k$，可以得到
$$\lim_{k \rightarrow \infty} \beta_k=\lambda_1$$
&emsp;&emsp;当使用无穷范数$\parallel \cdot \parallel_{\infty}$时，令$\beta_k=\frac{e_r^T u_k}{e_r^T y_{k-1}}$，设$u_{k-1}$模最大的分量为第$r$个分量，$e_r$是$n$维基本单位向量，则有
$$\lim_{k \rightarrow \infty} \beta_k=\frac{e_r^T A x_1}{e_r^T x_1}=\lambda_1$$
&emsp;&emsp;当$\frac{|\beta_k - \beta_{k-1}|}{|\beta_k|} \leq \epsilon$时，认为当前$\beta_k$与$\lambda_1$足够接近，算法停止。  

#### 1.2 反幂法
<font size=4>
&emsp;&emsp;反幂法主要用于计算矩阵按模最大的特征向量。  
&emsp;&emsp;设$n\times n$实矩阵$A$有$n$个线性无关的特征向量$x_1,x_2,\dots,x_n$，$x_i$对应的特征值$\lambda_i$满足
$$|\lambda_1| \geq |\lambda_2| \geq \cdots \geq |\lambda_{n-1}| > |\lambda_n|，Ax_i=\lambda_ix_i\ \ (i=1,2,\cdots,n)$$
&emsp;&emsp;由于$A$非奇异，故$\lambda_i \neq 0\ \ (i=1,2,\cdots,n)$。由$Ax_i=\lambda_ix_i$得
$$A^{-1}x_i=\frac{1}{\lambda_i}x_i$$
&emsp;&emsp;此时$\frac{1}{\lambda_n}$是$A^{-1}$的按模最大的特征值，$x_n$是对应的特征向量，对$A^{-1}$使用幂法即可求出$\frac{1}{\lambda_n}$，为了减少计算$A^{-1}$的计算量，将迭代公式变形为
$$
\begin{cases}
y_{k-1} = \frac{u_{k-1}}{\parallel u_{k-1} \parallel}\\
A u_k = y_{k-1}
\end{cases}
\ \ (k=1,2,\cdots)
$$
&emsp;&emsp;每次迭代需要对线性方程组$A u_k = y_{k-1}$进行求解，为了节省计算量，我们求解方程组时不使用高斯消去法，而是选择使用下面的Doolittle分解法。  
#### 1.3_Doolittle_分解法
&emsp;&emsp;_Doolittle_分解法主要用于求解求线性方程组$Ax=b$的解。通过将方程系数矩阵$A$分解为$A=LU$，其中$L$是下三角矩阵，$U$是上三角矩阵。此时方程组可化简成为两个容易求解的三角形方程组
$$Ly=b，Ux=y$$
先由$Ly=b$解出向量$y$，再由$Ux=y$解出向量$x$即为原方程组解向量。  
&emsp;&emsp;_Doolittle_分解计算公式：
$$
\begin{cases}
for\ i:=1 \rightarrow n \\
u_{ij}=a_{ij}- \sum_{k=1}^{i-1}l_{ik}u_{kj} \ \ (j=i,i+1,\cdots,n)\\
l_{ji}=(a_{ji} - \sum_{k=1}^{i-1} l_{jk}u_{ki})/u_{ii} \ \ (j=i+1,i+2,\cdots,n)
\end{cases}
$$
&emsp;&emsp;解三角方程组$Ly=b，Ux=y$的计算公式：
$$
\begin{cases}
y_1=b_1 \\
y_i=b_i-\sum_{j=1}^{i-1} l_{ij}y_j \ \ (i=2,3,\cdots,n) \\
x_n=y_n/u_{nn} \\
x_i=(y_i-\sum_{j=i+1}^n u_{ij}x_j)/u_{ii} \ \ (i=n-1,n-2,\cdots,1)
\end{cases}
$$  
#### 1.4带原点平移的幂法和反幂法
&emsp;&emsp;若$\lambda$是$A$的特征值，则$\lambda-p$是矩阵$A-pI$的特征值，其中$I$是单位矩阵；反之也成立，即
$$Ax=\lambda x，（A-pI）x=(\lambda - p)x$$
可互为因果。其中$p$称为原点位移。  
***
### 2.总体算法
#### 2.1求$\lambda_1,\lambda_{501},\lambda_s$
&emsp;&emsp;1.先使用反幂法求出按模最小的特征值$\lambda_s$  
&emsp;&emsp;2.使用幂法求出$A$的按模最大的特征值，记为$\lambda$  
&emsp;&emsp;3.对于$\lambda$，有
$$\lambda_1-\lambda \leq \lambda_2-\lambda \leq \cdots \leq \lambda_{501}-\lambda \leq 0$$
&emsp;&emsp;&emsp;此时按模最大的特征值为$\lambda_1-\lambda$，使用位移为$\lambda$的幂法求出这个值并求出$\lambda_1$  
&emsp;&emsp;4.若$\lambda\neq\lambda_1$,则$\lambda_{501}=\lambda$；若$\lambda = \lambda_1$，有
$$0 \leq \lambda_1+\lambda \leq \lambda_2+\lambda \leq \cdots \leq \lambda_{501}+\lambda$$
&emsp;&emsp;&emsp;此时按模最大的特征值为$\lambda_{501}+\lambda$，使用位移为$-\lambda$的幂法求出这个值并求出$\lambda_{501}$  
#### 2.2 求与$\mu_k=\lambda_1+k\frac{\lambda_{501} - \lambda_1}{40}$最接近的特征值
&emsp;&emsp;利用带原点平移的反幂法，设位移为$\mu_k$，则与其最接近的特征值$\lambda_{i_k}$满足
$$| \lambda_{i_k} | = \min_{1 \leq i \leq 501} {|\lambda_i-\mu_k|}$$
&emsp;&emsp;使用反幂法求出$\lambda_{i_k}-\mu_k$，并计算得到$\lambda_{i_k}$  
#### 2.3求$A$的谱范数条件数 $cond(A)_2$和行列式$detA$
##### 2.3.1$\ cond(A)_2$
&emsp;&emsp;对于非奇异的实对称矩阵$A$，有
$$cond(A)_2 = |\frac{\lambda_1}{\lambda_n}|$$
&emsp;&emsp;其中$\lambda_1$和$\lambda_n$分别是$A$的按模最大最小的特征值  
##### 2.3.2$\ detA$
&emsp;&emsp;由_Doolittle_分解法
$$A=LU \Longrightarrow det A = det L \ det U$$
&emsp;&emsp;而$L$的对角线上都为$1$，所以有
$$det A = \prod_{i=1}^{501}u_{ii}$$  
___
## Source Code
 ```c++
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
/*如果define NORM_2则使用2范数，如果define NORM_INF则使用无穷范数*/
const int N=501,N_MAX=501+10,inf=~0u>>1,r=2,s=2;
const double eps=1e-15,b=0.16,c=-0.064;

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
}//计算矩阵元素的值

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
}//将501阶的带状矩阵转化为r+s+1行501列的矩阵

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
}//位移为p，使用2范数的幂法（默认位移为0）
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
}//位移为p，使用无穷范数的幂法（默认位移为0）
#endif
inline int convert(int i,int j){
	return ((i-j+s+1>0)&&(i-j+s+1<6))?i-j+s+1:0;
}//坐标转换
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
}//Doolittle分解
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
}//位移为p的反幂法（默认位移为0）
double Det(){
	double ret=1;
	for (int i=1; i<=N; ++i)
		ret*=Doo[convert(i,i)][i];
	return ret;
}//计算行列式
int main(){
	Init();
	DoolittleInit();
	det=Det();
	lambda1=PowerMethod();
	lambda_s=InvPowerMethod();
	lambda501=PowerMethod(lambda1)+lambda1;
	cond=fabs(lambda1/lambda_s);
	if (lambda1>lambda501) swap(lambda1,lambda501);
    //用%.12E控制输出格式为科学计数法，且保留小数点后12位（此时有效数字至少为13位）
	printf("lambda_1=%.12E lambda_501=%.12E lambda_s=%.12E\n", lambda1, lambda501, lambda_s);
	for (int i=1; i<=39; ++i){
		double miu=lambda1+i*(lambda501-lambda1)/40;
		DoolittleInit(miu);
		printf("lambda_i%d=%.12E ", i, InvPowerMethod(miu)+miu);
	}
	printf("\ncond(A)=%.12E detA=%.12E\n", cond, det);
	
	return 0;
 ```    
 ___
 ## 结果
 