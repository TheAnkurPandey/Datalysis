#include<graphics.h>
#include<conio.h>
#include <stdlib.h>
#include<math.h>
void findCorrelation(float xx[],float yy[],int n){
	int i;
	float X[n],Y[n],meanx,meany,c,dr0=0,dr1=0,dr,nr=0,sumx=0,sumy=0;
	
	for(i=0;i<n;i++){
		sumx+=xx[i];
		sumy+=yy[i];
	}
	meanx=sumx/n;
	meany=sumy/n;
	for(i=0;i<n;i++){
		X[i]=xx[i]-meanx;
		Y[i]=yy[i]-meany;
		nr+=X[i]*Y[i];
	}
	for(i=0;i<n;i++){
		dr0+=X[i]*X[i];
		dr1+=Y[i]*Y[i];
	}
	dr=sqrt(dr0*dr1);
	c=nr/dr;
	printf("correlation coefficient = %f\n",c);
	printf("line of regression of x on y:\t");
	printf("y = %f + %fx\n",meany-((nr*meanx)/dr0),nr/dr0);
	printf("line of regression of y on x:\t");
	printf("x = %f + %fy\n",meanx-((nr*meany)/dr1),nr/dr1);           
}
float LIPF(float xx[],float yy[],int n,float x) {
	int i,j,k;
	float y=0,nr[n],dr[n];
	for(i=0;i<n;i++){
		nr[i]=1;
	}
	for(i=0;i<n;i++){
		dr[i]=1;
	}
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
			nr[i]=nr[i]*(x-xx[j]);
				
		for(k=0;k<n;k++){
			if(i!=k){
				dr[i]=dr[i]*(xx[i]-xx[k]);	
			}
		}
		y+= (nr[i]*yy[i])/(dr[i]*(x-xx[i]));
	}
	return y;
}
void equation(float xx[],float yy[],int n){
		int i,j,k;
		float dd[n][n];
		for(i=0;i<n;i++){
			dd[0][i]=yy[i];
		}
		for(i=1;i<n;i++){
    		k=i;
    		for(j=0;j<n-i;j++){
     			dd[i][j]=(dd[i-1][j+1]-dd[i-1][j])/(xx[k]-xx[j]);
     			k++;
    		}
  		}
		printf("%.2f  ",yy[0]);
		for(i=0;i<n-1;i++){
			printf(" +  ");
			for(j=0;j<=i;j++){
			printf("(x-%.2f)",xx[j]);
			}
			printf("%.2f ",dd[i][0]);
			printf("\t");		
		}
}
float f(float r,float fd[4][4]){
	float result;
	result= (r*r*fd[3][0])	/2 + r*(fd[2][0]-fd[3][0]) + (fd[1][0]-(fd[2][0])/2+ (fd[3][0])/3);
	return result;	
}
float d(float r,float fd[4][4]){
	float result;
	result= (r*fd[3][0]) + (fd[2][0]-fd[3][0]);
	return result;
}
float NRM(float fd[4][4]){
	float x1,x0,f1,d1;
	int c=0;
	printf("enter the value of x0");
	scanf("%f",&x0);
	do
	{
		f1=f(x0,fd);
		d1=d(x0,fd);
		x1 = x0 - ( f1 / d1 ) ;

		printf("x1 = %f\n",x1);
		x0=x1;
		c++;
	}while( (fabs( f(x1,fd) )>0.001) && (c<20) );
	printf("\nsolution is %f",x1);
	return x1;
}
void findOptimalPoint(float xx[],float yy[],int n){
	int i,j,k;
	float o[n-3],h=xx[1]-xx[0],fd[4][4],t,minv,min;
	for(i=0;i<=n-4;i++){
		for(j=0;j<4;j++){
			fd[0][j]=yy[i+j];
		}
		for(k=1;k<4;k++){
    		for(j=0;j<4-i;j++){
     			fd[k][j]=(fd[k-1][j+1]-fd[k-1][j]);
     		}
  		}
  		/*for(k=0;k<4;k++){
  			for(j=0;j<4-k;j++){
  			printf("%f\t",fd[k][j])		;
			}
			printf("\n");	
		}*/
		t=NRM(fd);
			//printf("\n r = %f",t);
		o[i] =xx[0] + h*t;
		
		//printf("\n h = %f",h);
		//printf("\n x = %f",o[i]);	
	}
	//printf("\n o[0] = %f",o[0]);
	minv=LIPF(xx,yy,n,o[0]);
	min=o[0];
	for(i=0;i<n-3;i++){
		if(minv<LIPF(xx,yy,n,o[i])){
			minv=LIPF(xx,yy,n,o[i]);
			min=o[i];
		}
	}
	printf("\noptimal solution is %f at %f",minv,min);
}
void equate(float xx[],float yy[],int n){
	int a[5],i,j,k;
	float c[n];
		float dd[n][n];
		for(i=0;i<n;i++){
			dd[0][i]=yy[i];
		}
		for(i=1;i<n;i++){
    		k=i;
    		for(j=0;j<n-i;j++){
     			dd[i][j]=(dd[i-1][j+1]-dd[i-1][j])/(xx[k]-xx[j]);
     			k++;
    		}
  		}
	for(i=0;i<n;i++){
		a[i]=1;
	}
	c[0]= yy[0] -(a[0]*xx[0]*dd[1][0]) + (a[1]*xx[0]*xx[1]*dd[2][0]) - (a[2]*xx[0]*xx[1]*xx[2]*dd[3][0]) + (a[3]*xx[0]*xx[1]*xx[2]*xx[3]*dd[4][0]);
	c[1]= dd[1][0] - (a[1]*(xx[0]+xx[1])*dd[2][0]) + (a[2]*(xx[0]*xx[1]+xx[0]*xx[2]+xx[1]*xx[2])*dd[3][0]) - (a[3]*(xx[0]*xx[1]*xx[2]+xx[0]*xx[1]*xx[3]+xx[0]*xx[2]*xx[3]*xx[3]+xx[1]*xx[2]*xx[3])*dd[4][0]);
	c[2]= dd[2][0] - (a[2]*(xx[0]+xx[1]+xx[2])*dd[3][0]) + (a[3]*(xx[0]*xx[1]+xx[0]*xx[2]+xx[0]*xx[3]+xx[1]*xx[2]+xx[1]*xx[3]+xx[2]*xx[3])*dd[4][0]);
	c[3]= dd[3][0] - (a[3]*(xx[0]+xx[1]+xx[2]+xx[3])*dd[4][0]);
	printf("%fx^3 + %fx^2+ %fx + %f",c[3],c[2],c[1],c[0]);
}

 
main()
{
   int gd = DETECT,gm,left=100,top=100,right=200,bottom=200,x= 300,y=150,radius=3,i,n=4;
   float yy[]={0.6221,0.6155,0.6138,0.6170},xx[]={0.60,0.65,0.7,0.75},miny=yy[2],maxy=yy[0],p[15],q[15];
 	initgraph(&gd, &gm, "C:\\TC\\BGI");
   line(left - 10, top + 150, left - 10, top - 70);
   line(left - 10, top + 150, left + 410, top + 150);
   
   for(i=0;i<10;i++){
   	circle(left-10 , top+150-22*i, radius);		
    }
   for(i=0;i<10;i++){
   	circle(left-10 + 42*i, top+150, radius);		
   	}
   
   for(i=0;i<4;i++){
   	circle(((xx[i]-xx[0])*336)/(xx[n-1]-xx[0])+132 , 228-( ((yy[i]-miny)*198)/(maxy-miny) ),radius  );
   }
   for(i=0;i<3;i++){
   	line(((xx[i]-xx[0])*336)/(xx[n-1]-xx[0])+132,228-(((yy[i]-miny)*198)/(maxy-miny)),((xx[i+1]-xx[0])*336)/(xx[n-1]-xx[0])+132,228-(((yy[i+1]-miny)*198)/(maxy-miny)));	
   }
   printf("x\t\t y\n");
   setcolor(RED);
   for(i=1;i<=10;i++){
   	p[i-1]=xx[0]+((xx[n-1]-xx[0])/(11))*i;
   	printf("%f\t",p[i-1]);
   	q[i-1]=LIPF(xx,yy,n,p[i-1]);
 	printf("%f\t",q[i-1]);
   	circle(((p[i-1]-xx[0])*336)/(xx[n-1]-xx[0])+132 , 228-( ((q[i-1]-miny)*198)/(maxy-miny) ),radius  );
   	printf("\n");
   }
   
   for(i=9;i>=1;i--){
   	line(((p[i-1]-xx[0])*336)/(xx[n-1]-xx[0])+132 , 228-( ((q[i-1]-miny)*198)/(maxy-miny) ),((p[i]-xx[0])*336)/(xx[n-1]-xx[0])+132 , 228-( ((q[i]-miny)*198)/(maxy-miny) ));	
   }
   setcolor(WHITE);
   outtextxy(left , top + 165, "x- axis");
   settextstyle(DEFAULT_FONT, VERT_DIR, 1);
   outtextxy(left -25, top + 15, "y- axis");
 
  getch();
   closegraph();
	if(n<=4)
		equate(xx,yy,4);
	else	
		equation(xx,yy,6);
	printf("\nat x=3.5 is y = %f\n ",LIPF(xx,yy,4,3.5));
	findCorrelation(xx,yy,4);
	//if(n>=4)
	//	findOptimalPoint(xx,yy,4); 	
   return 0;
}

