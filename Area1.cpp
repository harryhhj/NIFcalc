/**********************************************************************/
/*	    Filename: area1.cpp                                           */
/*		Function: Define Functions for Region 1 in NIF                */
/**********************************************************************/
#include "stdafx.h"

#include <math.h>
#include "area1.h"
#include "area4.h"
#include "areas.h"
#include "tools.h"

#define TX1 1386.0       /* K */
#define PX1 16.53        /* MPa */
#define VX1 0.001    /* m^3kg^-1  */
#define W2X1 1000.0   /* m^2s^-2 */


static double GE1n[34]  =
    {  0.14632971213167,      -0.84548187169114,     -0.37563603672040e01,
       0.33855169168385e01,   -0.95791963387872,      0.15772038513228,
      -0.16616417199501e-01,   0.81214629983568e-03,  0.28319080123804e-03,
      -0.60706301565874e-03,  -0.18990068218419e-01, -0.32529748770505e-01,
      -0.21841717175414e-01,  -0.52838357969930e-04, -0.47184321073267e-03,
      -0.30001780793026e-03,   0.47661393906987e-04, -0.44141845330846e-05,
      -0.72694996297594e-15,  -0.31679644845054e-04, -0.28270797985312e-05,
      -0.85205128120103e-09,  -0.22425281908000e-05, -0.65171222895601e-06,
      -0.14341729937924e-12,  -0.40516996860117e-06, -0.12734301741641e-08,
      -0.17424871230634e-09,  -0.68762131295531e-18,  0.14478307828521e-19,
       0.26335781662795e-22,  -0.11947622640071e-22,  0.18228094581404e-23,
      -0.93537087292458e-25
     };
static int GE1I[34]  =
     { 0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,
       4,4,4,5,8,8,21,23,29,30,31,32
     };
static int GE1J[34]  =
     { -2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1,3,17,-4,
       0,6,-5,-2,10,-8,-11,-6,-29,-31,-38,-39,-40,-41
     };

// region 1. backward equation
  // p,h-> t
static double Bh1n[20]  =
  {   -0.23872489924521e03,    0.40421188637945e03,   0.11349746881718e03,
      -0.58457616048039e01,   -0.15285482413140e-03, -0.10866707695377e-05,
      -0.13391744872602e02,    0.43211039183559e02,  -0.54010067170506e02,
       0.30535892203916e02,   -0.65964749423638e01,   0.93965400878363e-02,
       0.11573647505340e-06,  -0.25858641282073e-04, -0.40644363084799e-08,
       0.66456186191635e-07,   0.80670734103027e-10, -0.93477771213947e-12,
       0.58265442020601e-14,  -0.15020185953503e-16
 };
static int Bh1I[20] =
    { 0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,3,3,4,5,6 };
static int Bh1J[20] =
    {0,1,2,6,22,32,0,1,2,3,4,10,32,10,32,10,32,32,32,32 };

  // p,s-> t
static double Bs1n[20]  =
  {   0.17478268058307e03,0.34806930892873e02,0.65292584978455e01,
      0.33039981775489,-0.19281382923196e-06,-0.24909197244573e-22,
      -0.26107636489332,0.22592965981586,-0.64256463395226e-01,
      0.78876289270526e-02,0.35672110607366e-09,0.17332496994895e-23,
      0.56608900654837e-03,-0.32635483139717e-03,0.44778286690632e-04,
      -0.51322156908507e-09,-0.42522657042207e-25,0.26400441360689e-12,
      0.78124600459723e-28,-0.30732199903668e-30
 };
static int Bs1I[20] =
    { 0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,4 };
static int Bs1J[20] =
    { 0,1,2,3,11,31,0,1,2,3,12,31,0,1,2,9,31,10,32,32 };

/* definitions. */
double pt1v(double p, double t)
{
	double ret;
	double rdp = p / PX1;
	double rdt = TX1 / t;
	double Gmt1_p= - Gmt_p(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);

	ret = rdp * Gmt1_p;
	ret *= (R*t/p); 

	return (ret*VX1);
}
    
double pt1u(double p, double t)
{
	double ret;
	double rdp = p / PX1;
	double rdt = TX1 / t;
	double Gmt1_t= Gmt_t(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);
	double Gmt1_p= - Gmt_p(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);

	ret = rdt *Gmt1_t - rdp*Gmt1_p;
    ret *= (R*t);
	
	return ret;
}
    
double pt1s(double p, double t)
{
    double ret;
	double rdp = p / PX1;
	double rdt = TX1 / t;
	double Gmt1_t= Gmt_t(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);
	double Gmt1= Gmt(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);

	ret = (R)* (rdt *Gmt1_t - Gmt1);

	return ret;
}

double pt1h(double p, double t)
{
    double ret;
	double rdp = p / PX1;
	double rdt = TX1 / t;
	double Gmt1_t= Gmt_t(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);

	ret = rdt *Gmt1_t;
	ret *= (R*t);

	return ret;
}

double pt1cp(double p, double t)
{
	double ret;
	double rdp = p / PX1;
	double rdt = TX1 / t;
	double Gmt1_tt= Gmt_tt(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);


	ret = -rdt*rdt *Gmt1_tt;
	ret *= (R); 

	return ret;
}    

double pt1cv(double p, double t)
{
	double ret;
	double rdp = p / PX1;
	double rdt = TX1 / t;
	double Gmt1_tt= Gmt_tt(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);
	double Gmt1_p= - Gmt_p(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);
	double Gmt1_pp= Gmt_pp(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);
	double Gmt1_pt= - Gmt_pt(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);


	ret = -rdt*rdt * Gmt1_tt;
	ret += xen((Gmt1_p-rdt*Gmt1_pt),2) / Gmt1_pp;
	ret *= (R); 

	return ret;
}

double pt1w(double p, double t)
{
	double ret;
	double rdp = p / PX1;
	double rdt = TX1 / t;
	double Gmt1_p= - Gmt_p(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);
	double Gmt1_tt= Gmt_tt(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);
	double Gmt1_pp= Gmt_pp(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);
	double Gmt1_pt= - Gmt_pt(7.1-rdp, rdt-1.222, GE1n, GE1I, GE1J, 34);


	ret = Gmt1_p*Gmt1_p;
	ret /= (xen((Gmt1_p-rdt*Gmt1_pt),2) / (rdt*rdt*Gmt1_tt)
		     -Gmt1_pp);
	ret *= (R*t); 

	return sqrt(ret*W2X1);
}


// backward equations
double ph1t(double p,double h)
{
   double ret;
   double Pi = p/1.0;
   double Et = h/2500.0;

   ret = polyXY2e(Pi, Et+1.0, Bh1n, Bh1I, Bh1J, 20);

   return (ret*1.0);
}

double ps1t(double p,double s)
{
   double ret;
   double Pi = p/1.0;
   double Sg = s/1.0;

   ret = polyXY2e(Pi, Sg+2.0, Bs1n, Bs1I, Bs1J, 20);

   return (ret*1.0);
}

double ts1p(double t, double s)
{
	double a[3][3] = {{370, 1.26, 2.7167},          /* ==>ln(p) */
						{500, 2.46, 4.3917},
						{600, 3.54, 2.28}};
	double p0, p;

	p0 = exp(interpolate3(t, s, a));
	p = aitkin2r(ps1t, t, p0, s,  ERRORHI);
//	p = divid2r(ps1t, t, PB100, ps(t), s,  ERRORHI3);

	return p;
}

double th1p(double t, double h)
{  /* very good */
	double p0, p, err=1;
	int sum = 0;

	p0 = h - pt1u(50, t);
	
	while (err > ERRORHI && sum++ <1000)
	{
		p = (h-pt1u(p0, t))/(pt1v(p0,t)*1000);
		err = fabs(p-p0)/p0;
		p0 = p;
	}

	return p;
}

/*	non-linear question will cause dead cycle in half-division
double th1p(double t, double h)
{
	double p;

	p = divid2r(ph1t, t, 100, 0.0005, h, ERRORHI);
	return p;
}*/


double ts1h(double t, double s)
{
	double p, h;

	p = ts1p(t, s);
	h = pt1h(p, t);

	return h;
}

double th1s(double t, double h)
{
	double p, s;

	p = th1p(t, h);
	s = pt1s(p, t);

	return s;
}

/* obsolete implementation -1999 
double hs1t(double h, double s)
{
	double t, t0;

	t0=s4ts(s);
	t = aitkin2r(th1s, s, t0, h,  ERRORHI);

	return t;
} */

/* Sep.2001 Gaithersburg, Maryland, USA */
double hs1t(double h, double s) /* implemented Aug23, 2005 */
{
	double p;

	p = hs1p(h,s);
	return ph1t(p,h);
}

/* obsolete implementation
double hs1p(double h, double s)
{
	double p, t;

	t = hs1t(h, s);
	p = th1p(t, h);

	return p;
} */

/* Sep.2001 Gaithersburg, Maryland, USA */
double hs1p(double h, double s) /* implemented Aug23, 2005 */
{
	double px = 100;
	double hx = 3400;
	double sx = 7.6;
	int I[19] = {0,0,0,0,0, 0,0,0,1,1, 1,1,2,2,2, 3,4,4,5};
	int J[19] = {0,1,2,4,5, 6,8,14,0,1, 4,6,0,1,10, 4,1,4,0};
	double n[19] = {-0.691997014660582, -0.183612548787560e2,
		-0.928332409297335e1, 0.659639569909906e2, -0.162060388912024e2,
		0.450620017338667e3, 0.854680678224170e3, 0.607523214001162e4,
		0.326487682621856e2, -0.269408844582931e2, -0.319947848334300e3,
		-0.928354307043320e3, 0.303634537455249e2, -0.650540422444146e2,
		-0.430991316516130e4, -0.747512324096068e3, 0.730000345529245e3,
		0.114284032569021e4, -0.436407041874559e3};

	double rdp, rdh, rds;

	rdh = h/hx + 0.05;
	rds = s/sx + 0.05;
	rdp = polyXY2e(rdh, rds, n, I, J, 19);

	return (rdp * px);
}

/* p,s */
double ps1h(double p, double s)
{
	double t, h;

	t = ps1t(p,s);
	h = pt1h(p, t);

	return h;
}

/* p,v */
double pv1t(double p, double v)
{
	double t, t0;

	/* v->t */
	t0 = 700.38163 - 405.79078*exp((0.000977 - v)/0.000252514); 
	/* p->t */
	t0 += -30.77384 + 1.19365*exp(log(p)/1.1843);

	t = parabola2(pt1v, v, p, t0,  ERRORHI);

	return t;
}

double tv1p(double t, double v)
{
	double p;

	p = divid2r(pt1v, v, PB100, PB0, t, ERRORHI);

	return p;
}

double ph1v(double p, double h)
{
	double t;

	t=ph1t(p,h);
	return (pt1v(p,t));
}

double vh1p(double v, double h)
{
	double p;

	p=divid2r(ph1v, v, PB0, PB100, h, ERRORHI2);

	return p;
}

double ps1v(double p, double s)
{
	double t;

	t=ps1t(p,s);
	return (pt1v(p,t));
}

double vs1p(double v, double s)
{
	double p;

	p=divid2r(ps1v, v, PB0, PB100, s, ERRORHI2);

	return p;
}
