// Filename: tools.c
// Function: service subroutines
//
#include "stdafx.h"

#include <math.h>
#include <stdlib.h>
#include "tools.h"


// x^n
double xen(double x,  int n)
{
	int i=0, j=1, k=1;
	double retval=1.0;
	double xn[10];

	if(n<0) {k=-1; n =-n;}
	xn[0]=x;
	while(j<n) {j*=2; i++;xn[i]=xn[i-1]*xn[i-1];}
	
	for (i; i>=0; i--)
	{
		if (n/j == 1) retval *= xn[i];
		n %= j;
		j /= 2;
	}
	if (k<0) return (1.0/retval);
	else return (retval);
}

// x^n with preset x[] so faster 
double xpn(double x[],  int n)
{
	return (n<0) ? (1.0/x[-n]) : x[n];
}

void init_xpn_A(double x, double xs[], int rank)
{
	int i,m;

	m=rank>=0 ? rank : -rank;
		
	xs[0] = 1;
	for(i=1;i<=m; i++)
		xs[i]=xs[i-1] * x;
}

//initiate xpn(xs[], n), p[] is the power matrix I[] and J[]
void init_xpn_B(double x, double xs[], int p[], int rank)
{
	int i,m,n;

	m=p[0]>=0 ? p[0] : -p[0];
	for (i=1;i<rank;i++) //search for highest rank
	{
		n=(p[i]>=0) ? p[i] : -p[i];
		m=(m>=n) ? m : n;
	}
	
	xs[0] = 1;
	for(i=1;i<=m; i++)
		xs[i]=xs[i-1] * x;
}

double polynom (double para, int rank, double coef[])
{
 int i;
 double p = 0.0;

 for (i=rank; i>=0; i--)
    p  =  p * para + coef[i] ;

 return (p);
 }

// Sgm(a[i,j]*x^i*y^j)
double polynom2(double x,double y, double a[], int m, int n)
{
	int i, j;
	double sgm=0.0;

	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
		sgm +=a[i*n+j]*xen(x,i)*xen(y,j);

	return sgm;
}

// Sgm(a[i,j]*x^i*y^j) efficient version
double polynom2e(double x,double y, double a[], int m, int n)
{
	int i, j;
	double sgm=0.0, xs[MAX_N], ys[MAX_N];

	//initiate x and y power series
	init_xpn_A(x, xs, m);
	init_xpn_A(y, ys, n);


	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
		sgm +=a[i*n+j]*xpn(xs,i)*xpn(ys,j);

	return sgm;
}

// Sgm(n*x^I*y^J)
double polyXY2(double x,double y,double n[],int I[], int J[], int rank)
{
	int i;
	double sgm=0.0;

	for (i=0;i<rank;i++)
		sgm +=n[i]*xen(x,I[i])*xen(y,J[i]);

	return sgm;
}


// Sgm(n*x^I*y^J) efficient
double polyXY2e(double x,double y,double n[],int I[], int J[], int rank)
{
	int i;
	double sgm=0.0, xs[MAX_N], ys[MAX_N];

	//initiate x and y power series
	init_xpn_B(x, xs, I, rank);
	init_xpn_B(y, ys, J, rank);

	for (i=0;i<rank;i++)
		sgm +=n[i]*xpn(xs,I[i])*xpn(ys,J[i]);

	return sgm;
}



/*
  Aitkin accelaration method
  obj=f(var)   <===>   F(obj,var)=obj-f(var)=0
  var[k+1]= var[k]+F{obj,var[k]} = var[k]+obj-f{var[k]}
*/
double aitkin1(double (*fp)(double),double obj,double initial_var,
			      double precision)
{
     int i=1;
     double x1 , _x1 , _x2, err, err_x1;

	 err = fabs(initial_var) < 0.0001 ? 0.0001 : initial_var;

     _x1=initial_var+obj-(*fp)(initial_var);
     if(fabs((initial_var-_x1)/err)<=precision)
	     return( (initial_var+_x1)/2.0 );

	err_x1 = fabs(_x1) < 0.0001 ? 0.0001 : _x1;

     _x2=initial_var+obj-(*fp)(_x1);
     if(fabs((_x1-_x2)/err_x1)<=precision)
	     return( (_x1+_x2)/2.0 );
     x1=(initial_var*_x2-_x1*_x1)/(initial_var-2.0*_x1+_x2);
     while(fabs((initial_var-x1)/err)>precision)
	  {
	     if(i>1000)  {
/*		     cout<<"\r\n";
		     cout<<"Iration times have exceeded 1000 , Failed ! ";
*/		     return ERRCALC;
			  }
	     else         {
		     initial_var=x1;
		     _x1=initial_var+obj-(*fp)(initial_var);
		     _x2=initial_var+obj-(*fp)(_x1);
		     x1=(initial_var*_x2-_x1*_x1)/(initial_var-2.0*_x1+_x2);
		     i+=1;
			  }
	  }
     return( (initial_var+x1)/2.0 );
}


/*
  Aitkin accelaration method
  obj=f(const,var)   <===>   F(obj,const,var)=obj-f(const,var)=0
  var[k+1]= var[k]+F{obj,const,var[k]} = var[k]+obj-f{const,var[k]}
*/
double aitkin2(double (*fp)(double constant,double initial_var),double obj,double constant,double initial_var,
	       double precision)
{
     int i=1;
     double x1 , _x1 , _x2, err, err_x1;

	 err = fabs(initial_var) < 0.0001 ? 0.0001 : initial_var;

     _x1=initial_var+obj - (*fp)(constant,initial_var);
     if(fabs((initial_var-_x1)/err)<=precision)
	     return( (initial_var+_x1)/2.0 );

	err_x1 = fabs(_x1) < 0.0001 ? 0.0001 : _x1;
     _x2=initial_var+obj - (*fp)(constant,_x1);
     if(fabs((_x1-_x2)/err_x1)<=precision)
	     return( (_x1+_x2)/2.0 );
     x1=(initial_var*_x2-_x1*_x1)/(initial_var-2.0*_x1+_x2);
     while(fabs((initial_var-x1)/err)>precision)
	  {
	     if(i>1000)  {
//		     cout<<"\r\n";
//		     cout<<"Iration times have exceeded 1000 , Failed ! ";
		     return ERRCALC;
			  }
	     else         {
		     initial_var=x1;
		     _x1=initial_var+obj - (*fp)(constant,initial_var);
		     _x2=initial_var+obj - (*fp)(constant,_x1);
		     x1=(initial_var*_x2-_x1*_x1)/(initial_var-2.0*_x1+_x2);
		     i++;
			  }
	  }
     return( (initial_var+x1)/2.0 );
}

double aitkin2r(double (*fp)(double initial_var,double constant),double obj,double initial_var,double constant,
	       double precision)
{
     int i=1;
     double x1 , _x1 , _x2, err, err_x1;

	 err = fabs(initial_var) < 0.0001 ? 0.0001 : initial_var;


     _x1=initial_var+obj - (*fp)(initial_var,constant);
     if(fabs((initial_var-_x1)/err)<=precision)
	     return( (initial_var+_x1)/2.0 );

	 err_x1 = fabs(_x1) < 0.0001 ? 0.0001 : _x1;

     _x2=initial_var+obj - (*fp)(_x1,constant);
     if(fabs((_x1-_x2)/err_x1)<=precision)
	     return( (_x1+_x2)/2.0 );
     x1=(initial_var*_x2-_x1*_x1)/(initial_var-2.0*_x1+_x2);

     while(fabs((initial_var-x1)/err)>precision)
	  {
	     if(i>1000)  {
/*			 if (INTERACT)
			 {
		     cout<<"\r\n";
		     cout<<"Iration times have exceeded 1000 , Failed ! ";
			 }
*/		     return ERRCALC;
			  }
	     else         {
		     initial_var=x1;
		     _x1=initial_var+obj - (*fp)(initial_var,constant);
		     _x2=initial_var+obj - (*fp)(_x1,constant);
		     x1=(initial_var*_x2-_x1*_x1)/(initial_var-2.0*_x1+_x2);
		     i++;
			  }
	  }
     return( (initial_var+x1)/2.0 );
}

/* parabolic algorithm for seek root for function */
double parabola1(double (*fp)(double),double obj,double initial_var,
			      double precision)
{
     int i=1;
     double x0, x1, x2, y0, y1, y2, err;

     x0=initial_var;
     y0=(*fp)(x0);
     x1=x0+x0*(obj-y0)/y0;
     y1=(*fp)(x1);
     do
	 {
		 x2=(obj-y1)/(y0-y1)*x0 +(obj-y0)/(y1-y0)*x1;
		 y2 = (*fp)(x2);
		 
		 x0 = x1; y0=y1;
		 x1=x2; y1=y2;

		 err = fabs(x2) < 0.0001 ? 0.0001 : x2;

	 } while (fabs((x1-x2)/err)>precision && i++<1000);

     if (i>=1000) return ERRCALC;
//   cout<<"Iration times have exceeded 1000 , Failed ! ";

     return x2;
}

double parabola2(double (*fp)(double constant,double initial_var),double obj,double constant,double initial_var,
	       double precision)
{
     int i=1;
     double x0, x1, x2, y0, y1, y2, err;

     x0=initial_var;
     y0=(*fp)(constant, x0);
     x1=x0+x0*(obj-y0)/y0;
     y1=(*fp)(constant, x1);
     do
	 {
		 x2=(obj-y1)/(y0-y1)*x0 +(obj-y0)/(y1-y0)*x1;
		 y2 = (*fp)(constant, x2);
		 
		 x0 = x1; y0=y1;
		 x1=x2; y1=y2;

		 err = fabs(x2) < 0.0001 ? 0.0001 : x2;
	 } while (fabs((x1-x2)/err)>precision && i++<1000);

     if (i>=1000) return ERRCALC;
//   cout<<"Iration times have exceeded 1000 , Failed ! ";

     return x2;
}

double parabola2r(double (*fp)(double initial_var,double constant),double obj,double initial_var,double constant,
	       double precision)
{
     int i=1;
     double x0, x1, x2, y0, y1, y2, err;

     x0=initial_var;
     y0=(*fp)(x0, constant);
     x1=x0+x0*(obj-y0)/y0;
     y1=(*fp)(x1, constant);
     do
	 {
		 x2=(obj-y1)/(y0-y1)*x0 +(obj-y0)/(y1-y0)*x1;
		 y2 = (*fp)(x2, constant);
		 
		 x0 = x1; y0=y1;
		 x1=x2; y1=y2;

		 err = fabs(x2) < 0.0001 ? 0.0001 : x2;
	 } while (fabs((x1-x2)/err)>precision && i++<1000);

     if (i>=1000) return ERRCALC;
//   cout<<"Iration times have exceeded 1000 , Failed ! ";

     return x2;
}

/* half divid search root */
double divid1(double (*fp)(double),double obj,double r1, double r2,
			      double precision)
{
	int i=0;
    double x1, x2, x, y, err;

	err = fabs(obj) < 0.0001 ? 0.0001 : obj;

	y = (*fp)(r1);
	if (fabs((y-obj)/err) < precision) return y;
	x = (*fp)(r2);
	if (fabs((x-obj)/err) < precision) return x;

	if (y > x)
	{
		x1=r1;
		x2=r2;
	}
	else
	{
		x1 = r2;
		x2 = r1;
	}

	do 
	{
		x = (x1+x2)/2;
		y = (*fp)(x);
		if (y>obj) x1=x;
		else x2 = x;
	} while (fabs((y-obj)/err)>precision && i++<30);

	if (i>=30) return ERRCALC;

    return x;
}

double divid2(double (*fp)(double constant,double initial_var),double obj,double constant,double r1, double r2,
	       double precision)
{
	int i=0;
    double x1, x2, x, y, err;
	err = fabs(obj) < 0.0001 ? 0.0001 : obj;

	y = (*fp)(constant, r1);
	if (fabs((y-obj)/err) < precision) return y;
	x = (*fp)(constant, r2);
	if (fabs((x-obj)/err) < precision) return x;

	if (y > x)
	{
		x1=r1;
		x2=r2;
	}
	else
	{
		x1 = r2;
		x2 = r1;
	}

	do 
	{
		x = (x1+x2)/2;
		y = (*fp)(constant, x);
		if (y>obj) x1=x;
		else x2 = x;
	} while (fabs((y-obj)/err)>precision && i++<40);

	if (i>=40) return ERRCALC;

    return x;
}

double divid2r(double (*fp)(double initial_var,double constant),double obj, double r1, double r2,double constant,
	       double precision)
{
	int i=0;
    double x1, x2, x, y, err;
	
	err = fabs(obj) < 0.0001 ? 0.0001 : obj;

	y = (*fp)(r1, constant);
	if (fabs((y-obj)/err) < precision) return y;
	x = (*fp)(r2, constant);
	if (fabs((x-obj)/err) < precision) return x;

	if (y > x)
	{
		x1=r1;
		x2=r2;
	}
	else
	{
		x1 = r2;
		x2 = r1;
	}

	do 
	{
		x = (x1+x2)/2;
		y = (*fp)(x, constant);
		if (y>obj) x1=x;
		else x2 = x;
	} while (fabs((y-obj)/err)>precision && i++<40);

	if (i>=40) return ERRCALC;

    return x;
}



// solve multi-dim function (dissolve method)
int ciggj(double a[],int n, double b[])
  { int *js,i,j,k,is,u,v;
    double d,t;
    js=(int *)malloc(n*sizeof(int));
    for (k=0; k<=n-1; k++)
      { d=0.0;
        for (i=k; i<=n-1; i++)
        for (j=k; j<=n-1; j++)
          { t=fabs(a[i*n+j]);
            if (t>d) {d=t; js[k]=j; is=i;}
          }
        if (d+1.0==1.0)
          { free(js);  return(0);}
        if (is!=k)
          { for (j=k; j<=n-1; j++)
              { u=k*n+j; v=is*n+j;
                t=a[u]; a[u]=a[v]; a[v]=t;
              }
            t=b[k]; b[k]=b[is]; b[is]=t;
          }
        if (js[k]!=k)
          for (i=0; i<=n-1; i++)
            { u=i*n+k; v=i*n+js[k];
              t=a[u]; a[u]=a[v]; a[v]=t;
            }
        t=a[k*n+k];
        for (j=k+1; j<=n-1; j++)
          { u=k*n+j;
            if (a[u]!=0.0) a[u]=a[u]/t;
          }
        b[k]=b[k]/t;
        for (j=k+1; j<=n-1; j++)
          { u=k*n+j;
            if (a[u]!=0.0)
              { for (i=0; i<=n-1; i++)
                  { v=i*n+k;
                    if ((i!=k)&&(a[v]!=0.0))
                      { is=i*n+j;
                        a[is]=a[is]-a[v]*a[u];
                      }
                  }
              }
          }
        for (i=0; i<=n-1; i++)
          { u=i*n+k;
            if ((i!=k)&&(a[u]!=0.0))
              b[i]=b[i]-a[u]*b[k];
          }
      }
    for (k=n-1; k>=0; k--)
      if (k!=js[k])
        { t=b[k]; b[k]=b[js[k]]; b[js[k]]=t;}
    free(js);
    return(1);
  }

// interpolate from 3 point
double interpolate3(double x, double y, double c[3][3])
{
	double b[3] = {1, 1, 1};
	double ret;

	if (ciggj(&c[0][0], 3, b) != 1) return ERRCALC;

	ret = (1.0 - b[0]*x - b[1]*y)/b[2];

	return ret;
}

double interpolate4(double x, double y, double vx[], double vy[], int m, int n, double *vr)
{
	int xx=0, yy=0;
	double kx, ky;
	double v[3][3];
	
	while (x >= vx[xx] && xx<m) {xx++;};  // get rectangular area
	while (y >= vy[yy] && yy<n) {yy++;};

	if (*(vr+n*xx+yy) >0 && *(vr+n*(xx-1)+(yy-1)) >0)
	{
		v[0][0] = vx[xx];   v[0][1] = vy[yy];   v[0][2] = *(vr+n*(xx)+yy);
		v[1][0] = vx[xx-1]; v[1][1] = vy[yy-1]; v[1][2] = *(vr+n*(xx-1)+(yy-1));

		ky = (y-vy[yy-1])/(vy[yy]-vy[yy-1]);
		kx = (x-vx[xx-1])/(vx[xx]-vx[xx-1]);

		if ((ky >= kx) && *(vr+n*(xx-1)+yy)>0) 
		{
			v[2][0] = vx[xx-1]; v[2][1] = vy[yy];   v[2][2] = *(vr+n*(xx-1)+yy);
		}
		else
		{
			v[2][0] = vx[xx]; v[2][1] = vy[yy-1];   v[2][2] = *(vr+n*xx+(yy-1));
		}
	}
	else
	{
		v[0][0] = vx[xx-1]; v[0][1] = vy[yy-1]; v[0][2] = *(vr+n*(xx-1)+(yy-1));
		v[1][0] = vx[xx];   v[1][1] = vy[yy];   v[1][2] = *(vr+n*xx+yy);

		ky = (y-vy[yy-1])/(vy[yy]-vy[yy-1]);
		kx = (vx[xx]-x)/(vx[xx]-vx[xx-1]);
		
		if ((ky >= kx) && *(vr+n*xx+yy)>0) 
		{
			v[2][0] = vx[xx]; v[2][1] = vy[yy];   v[2][2] = *(vr+n*xx+yy);
		}
		else
		{
			v[2][0] = vx[xx-1]; v[2][1] = vy[yy-1];   v[2][2] = *(vr+n*(xx-1)+(yy-1));
		}
	}

		return interpolate3(x, y, v);
}

double interpolate2(double x, int n, double c[])
{
	int i=0;
	double ret;

	while (x<c[i])	i++;

	ret = (x-c[i-1])*(c[i]-c[i-1])/(c[n+i]-c[n+i-1]);
	
	return ret;
}

/* comparison
double min(double a, double b)
{
	return a<b?a:b;
}

double max(double a, double b)
{
	return a>b?a:b;
} */
