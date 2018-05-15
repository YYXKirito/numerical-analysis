#include<stdio.h>
#include<math.h>
#include<stdlib.h>
int main(){
	double A[10][10],Q[10][10],z[10][10],x[10][10],d = 0,c = 0,temp,t,h = 0,u[10] = {0},p[10] = {0},q[10] = {0},w[10] = {0};
	int i,j,k,n,r,o = 1,O = 1;
		FILE *fp;//建立一个文件操作指针
 	fp=fopen("gzw.txt","w");//以追加的方式建立或打开1.txt，默认位置在你程序的目录下面
	for(i = 0; i < 10; i++){
		for(j = 0; j < 10; j++){
			if(j==i){
				A[i][j] = 1.52*cos(2.2*(j+1));
			}
			else{
				A[i][j] = sin(0.5*(i+1) + 0.2*(j+1));
			}
		}
	} 
	for(i = 0; i < 10; i++){
		for(j =0; j < 10; j++){
			fprintf(fp,"%12.8f    ",A[i][j]);
		}
		fprintf(fp,"\n");				
	}
	fprintf(fp,"\n");	
	for(r = 0; r < 8; r++){
		O = 1;
		o = 1;
		for(i = r+2; i < 10; i++){
			if(A[i][r] != 0){
				o = 0;
			}
			O = O * o;
		}
		if(O == 0){
			d = 0;
			for(i = r+1; i < 10; i++){
				d = d + A[i][r]*A[i][r];
			}				
			d = sqrt(d);
			if(A[r+1][r] > 0){
				c = -d;
			}
			else{
				c = d;
			}
			h = c*c - c*A[r+1][r];
			for(i = 0; i <= r; i++){
				u[i] = 0;
			}
			u[r+1] = A[r+1][r] - c;
			for(i = r+2; i < 10; i++){
				u[i] = A[i][r];
			}
			for(i = 0; i < 10; i++){
				p[i] = 0;
				for(j =0; j < 10; j++){
					p[i] = p[i] + A[j][i]*u[j];
				}
				p[i] = p[i]/h;
			}
			for(i = 0; i < 10; i++){
				q[i] = 0;
				for(j =0; j < 10; j++){
					q[i] = q[i] + A[i][j]*u[j];
				}
				q[i] = q[i]/h;
			}
			t = 0;
			for(i = 0; i < 10; i++){
				t = t + p[i]*u[i];
			}
			t = t/h; 
			for(i = 0; i < 10; i++){
				w[i] = q[i] - t*u[i];
			}
			for(i = 0; i < 10; i++){
				for(j =0; j < 10; j++){
					A[i][j] = A[i][j] - w[i]*u[j] - u[i]*p[j];
				}				
			}
		}
	}
	for(i = 0; i < 10; i++){
			for(j =0; j < 10; j++){
				fprintf(fp,"%12.8f    ",A[i][j]);
			}
			fprintf(fp,"\n");				
		}
		fprintf(fp,"\n");	//打印拟上三角化矩阵 
	for(i = 0; i < 10; i++){
		for(j =0; j < 10; j++){
			Q[i][j] = 0;
		}				
	}
	for(i = 0; i < 10; i++){
		Q[i][i] = 1;			
	}
	for(r = 0; r < 9; r++){
		O = 1;
		o = 1;
		for(i = r+1; i < 10; i++){
			if(A[i][r] != 0){
				o = 0;
			}
			O = O * o;
		}
		if(O == 0){
			d = 0;
			for(i = r; i < 10; i++){
				d = d + A[i][r]*A[i][r];
			}				
			d = sqrt(d);
			if(A[r][r] > 0){
				c = -d;
			}
			else{
				c = d;
			}
			h = c*c - c*A[r][r];
			fprintf(fp,"c = %12.8f  h =  %12.8f \n", c, h);
			for(i = 0; i <= r-1; i++){
				u[i] = 0;
			}
			u[r] = A[r][r] - c;
			for(i = r+1; i < 10; i++){
				u[i] = A[i][r];
			}
			for (int i=0; i<=8; ++i)
				fprintf(fp,"%12.8f ", u[i]);
			fprintf(fp,"%12.8f\n ", u[9]);
			for(i = 0; i < 10; i++){
				w[i] = 0;
				for(j =0; j < 10; j++){
					w[i] = w[i] + Q[i][j]*u[j];
				}
			}
			for(i = 0; i < 10; i++){
				for(j =0; j < 10; j++){
					Q[i][j] = Q[i][j] - w[i]*u[j]/h;
				}				
			}
			for(i = 0; i < 10; i++){
				p[i] = 0;
				for(j =0; j < 10; j++){
					p[i] = p[i] + A[j][i]*u[j];
				}
				p[i] = p[i]/h;
			} 
			for(i = 0; i < 10; i++){
				for(j =0; j < 10; j++){
					A[i][j] = A[i][j] - u[i]*p[j];
				}				
			}
		}
	}
	for(i = 0; i < 10; i++){
			for(j =0; j < 10; j++){
				fprintf(fp,"%14.12e     ",A[i][j]);
			}
			fprintf(fp,"\n");				
		}
		fprintf(fp,"\n");	
	for(i = 0; i < 10; i++){
			for(j =0; j < 10; j++){
				fprintf(fp,"%14.12e     ",Q[i][j]);
			}
			fprintf(fp,"\n");				
		}
		fprintf(fp,"\n");
	for(i = 0; i < 10; i++){
			for(j =0; j < 10; j++){
				temp = 0;
				for(k =0; k < 10; k++){
					temp = temp + A[i][k]*Q[k][j];
 				}
 				z[i][j] = temp;
 				fprintf(fp,"%14.12e  ",z[i][j]);
			}
			fprintf(fp,"\n");				
		}
		fprintf(fp,"\n");	
	for(i = 0; i < 10; i++){
			for(j =0; j < 10; j++){
				temp = 0;
				for(k =0; k < 10; k++){
					temp = temp + Q[i][k]*A[k][j];
 				}
 				x[i][j] = temp;
 				fprintf(fp,"%14.12e  ",x[i][j]);
			}
			fprintf(fp,"\n");				
		}
		fprintf(fp,"\n");	
	fclose(fp);
	return 0;
}
