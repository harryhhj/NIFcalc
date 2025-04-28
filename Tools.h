/*
file name :  tools.h
*/

/* error limits */
#define ERRORHI4 1E-10
#define ERRORHI3 1E-8
#define ERRORHI2 1E-6
#define ERRORHI 1E-5
#define ERRORMD 0.0005
#define ERRORLO 0.005
#define ERRCALC  -9
#define MAX_N	64

int report_err(int err);

// x^n
double xen(double x,  int n);
// xpn() - x^n with preset x[] so faster for re-use, need to initiate the x[]
double xpn(double x[],  int n);
void init_xpn_A(double x, double xs[], int rank);
void init_xpn_B(double x, double xs[], int p[], int rank);

double polynom(double para, int rank, double coef[]);

// Sgm(n*x^I*y^J)
double polyXY2(double x,double y,double n[],int I[], int J[], int rank);
double polyXY2e(double x,double y,double n[],int I[], int J[], int rank); //efficient version
// Sgm(n[i,j]*x^i*y^j)
double polynom2(double x,double y, double a[], int m, int n);
double polynom2e(double x,double y, double a[], int m, int n);



/*
  Aitkin accelaration method
  obj=f(var)   <===>   F(obj,var)=obj-f(var)=0
  var[k+1]= var[k]+F{obj,var[k]} = var[k]+obj-f{var[k]}
*/
double aitkin1(double (*fp)(double),double object,double initial_var,
			      double precision);

/*
  Aitkin accelaration method
  obj=f(const,var)   <===>   F(obj,const,var)=obj-f(const,var)=0
  var[k+1]= var[k]+F{obj,const,var[k]} = var[k]+obj-f{const,var[k]}
*/
double aitkin2(double (*fp)(double, double),double object,double constant,double initial_var,
	       double precision);
/* obj=f(var, const) */
double aitkin2r(double (*fp)(double, double),double object,double initial_var,double constant,
	       double precision);

/* parabolic algorithm for seek root for function */
double parabola1(double (*fp)(double),double obj,double initial_var,
			      double precision);
double parabola2(double (*fp)(double constant,double initial_var),double obj,double constant,double initial_var,
	       double precision);
double parabola2r(double (*fp)(double initial_var,double constant),double obj,double initial_var,double constant,
	       double precision);

/* half divid search root */
double divid1(double (*fp)(double),double obj,double r1, double r2,
			      double precision);
double divid2(double (*fp)(double constant,double initial_var),double obj,double constant,double r1, double r2,
	       double precision);
double divid2r(double (*fp)(double initial_var,double constant),double obj, double r1, double r2,double constant,
	       double precision);


// solve multi-dim function (dissolve method)
int ciggj(double a[],int n, double b[]);

/* interpolatation */
double interpolate3(double x, double y, double c[3][3]);
double interpolate4(double x, double y, double vt[], double vp[], int m, int n, double *vr);
double interpolate2(double x, int n, double c[]);

/* comparison 
double min(double a, double b);
double max(double a, double b);*/
