/**********************************************************************/
/*	     Filename: area4.cpp                                          */
/*		 Function: NIF area 4 calculations                            */
/**********************************************************************/
#include "stdafx.h"

#include <math.h>
#include "area1.h"
#include "area2.h"
#include "area3.h"
#include "area4.h"
#include "areas.h"
#include "tools.h"
#define T3Aux 643.2

static double Sat4n[10]  =
  { 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
    0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
   -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849e0,
    0.65017534844798e3
  };

// Saturation line
double t4p(double t)
{
  double A, B, C, s, rdt, rdp;

 rdt = t + Sat4n[8]/(t - Sat4n[9]);
 s = rdt * rdt;
 A = s +Sat4n[0]*rdt +Sat4n[1];
 B = Sat4n[2]*s + Sat4n[3]*rdt + Sat4n[4];
 C = Sat4n[5]*s + Sat4n[6]*rdt + Sat4n[7];

 rdp = 2*C/(sqrt(B*B-4*A*C)-B);
 rdp *= rdp;
 rdp *= rdp;

 return (rdp);
}

double p4t(double p)
{
  double D, E, F, G, s, rdp, rdt;

  rdp = pow(p, 1.0/4.0);
 s = rdp * rdp;
 E = s +Sat4n[2]*rdp +Sat4n[5];
 F = Sat4n[0]*s + Sat4n[3]*rdp + Sat4n[6];
 G = Sat4n[1]*s + Sat4n[4]*rdp + Sat4n[7];

 D = 2*G/(-sqrt(F*F-4*E*G)-F);
 rdt = (Sat4n[9]+ D- sqrt((Sat4n[9]+D)*(Sat4n[9]+D) -
	 4*(Sat4n[8]+Sat4n[9]*D))) /2;

 return (rdt);
}


double t4vs1(double t)
{
	double p;
	double a[4]={-10.83596,0.05187347,-8.27838E-5,4.40493E-8};

	p=t4p(t);
	if (t<=TB2) return pt1v(p,t);
	else return pt3v(p+ERRORHI,t); //Jul 2005, Greece
}

double t4vs2(double t)
{
	double p;

	p=t4p(t);
	if (t<=TB2) return pt2v(p,t);
	else return pt3v(p-ERRORHI,t); //Jul 2005, Greece

}

double t4us1(double t)
{
	double p,u;

	p=t4p(t);
	if (t<=TB2) u=pt1u(p,t);
	else u=rt3u(1/t4vs1(t), t);
	
	return u;
}


double t4us2(double t)
{
	double p,u;

	p=t4p(t);
	if (t<=TB2) u=pt2u(p,t);
	else u=rt3u(1/t4vs2(t), t);
	
	return u;
}


double t4hs1(double t)
{
	double p,h;

	p=t4p(t);
	if (t<=TB2) h=pt1h(p,t);
	else h=rt3h(1/t4vs1(t), t);
	
	return h;
}

double t4hs2(double t)
{
	double p,h;

	p=t4p(t);
	if (t<=TB2) h=pt2h(p,t);
	else h=rt3h(1/t4vs2(t), t);
	
	return h;
}

double t4ss1(double t)
{
	double p,s;

	p=t4p(t);
	if (t<=TB2) s=pt1s(p,t);
	else s=rt3s(1/t4vs1(t), t);
	
	return s;
}

double t4ss2(double t)
{
	double p,s;

	p=t4p(t);
	if (t<=TB2) s=pt2s(p,t);
	else s=rt3s(1/t4vs2(t), t);
	
	return s;
}

double t4cps1(double t)
{
	double p,c;

	p=t4p(t);
	if (t<=TB2) c=pt1cp(p,t);
	else c=rt3cp(1/t4vs1(t), t);
	
	return c;
}

double t4cps2(double t)
{
	double p,c;

	p=t4p(t);
	if (t<=TB2) c=pt2cp(p,t);
	else c=rt3cp(1/t4vs2(t), t);
	
	return c;
}

double t4cvs1(double t)
{
	double p,c;

	p=t4p(t);
	if (t<=TB2) c=pt1cv(p,t);
	else c=rt3cv(1/t4vs1(t), t);
	
	return c;
}

double t4cvs2(double t)
{
	double p,c;

	p=t4p(t);
	if (t<=TB2) c=pt2cv(p,t);
	else c=rt3cv(1/t4vs2(t), t);
	
	return c;
}

double t4ws1(double t)
{
	double p,w;

	p=t4p(t);
	if (t<=TB2) w=pt1w(p,t);
	else w=rt3w(1/t4vs1(t), t);
	
	return w;
}

double t4ws2(double t)
{
	double p,w;

	p=t4p(t);
	if (t<=TB2) w=pt2w(p,t);
	else w=rt3w(1/t4vs2(t), t);
	
	return w;
}

double p4vs1(double p)
{
	return t4vs1(p4t(p));
}

double p4vs2(double p)
{
	return t4vs2(p4t(p));
}

double p4hs1(double p)
{
	return t4hs1(p4t(p));
}

double p4hs2(double p)
{
	return t4hs2(p4t(p));
}

double p4ss1(double p)
{
	return t4ss1(p4t(p));
}

double p4ss2(double p)
{
	return t4ss2(p4t(p));
}

double p4cps1(double p)
{
	return t4cps1(p4t(p));
}

double p4cps2(double p)
{
	return t4cps1(p4t(p));
}

double p4ws1(double p)
{
	return t4ws1(p4t(p));
}

double p4ws2(double p)
{
	return t4ws2(p4t(p));
}

double tx4v(double t, double x)
{
	if (x<ERRORHI3) return t4vs1(t);
	if (x>(1-ERRORHI3)) return t4vs2(t);

	return ((1-x)*t4vs1(t) + x*t4vs2(t));
}

double px4v(double p, double x)
{
	return tx4v(p4t(p),x);
}

double tx4h(double t, double x)
{
	if (x<ERRORHI3) return t4hs1(t);
	if (x>(1-ERRORHI3)) return t4hs2(t);

	return ((1-x)*t4hs1(t) + x*t4hs2(t));
}

double tx4u(double t, double x)
{
	if (x<ERRORHI3) return t4us1(t);
	if (x>(1-ERRORHI3)) return t4us2(t);

	return ((1-x)*t4us1(t) + x*t4us2(t));
}

double tx4cp(double t, double x)
{
	if (x<ERRORHI3) return t4cps1(t);
	if (x>(1-ERRORHI3)) return t4cps2(t);

	return ((1-x)*t4cps1(t) + x*t4cps2(t));
}

double tx4cv(double t, double x)
{
	if (x<ERRORHI3) return t4cvs1(t);
	if (x>(1-ERRORHI3)) return t4cvs2(t);

	return ((1-x)*t4cvs1(t) + x*t4cvs2(t));
}

double tx4w(double t, double x)
{
	if (x<ERRORHI3) return t4ws1(t);
	if (x>(1-ERRORHI3)) return t4ws2(t);

	return ((1-x)*t4ws1(t) + x*t4ws2(t));
}

double px4h(double p, double x)
{
	return tx4h(p4t(p),x);
}

double tx4s(double t, double x)
{
	if (x<ERRORHI3) return t4ss1(t);
	if (x>(1-ERRORHI3)) return t4ss2(t);
	return ((1-x)*t4ss1(t) + x*t4ss2(t));
}

double px4s(double p, double x)
{
	return tx4s(p4t(p),x);
}

double tv4x(double t, double v)
{
	double vs1, vs2;

	vs1=t4vs1(t);
	vs2=t4vs2(t);
	return ((v-vs1)/(vs2-vs1));
}

double th4x(double t, double h)
{
	double hs1, hs2;

	hs1=t4hs1(t);
	hs2=t4hs2(t);
	return ((h-hs1)/(hs2-hs1));
}

double ts4x(double t, double s)
{
	double ss1, ss2;

	ss1=t4ss1(t);
	ss2=t4ss2(t);
	return ((s-ss1)/(ss2-ss1));
}

double vx4t(double v, double x)
{
	double t;
/*	double a[5]={-0.74908,-1.05451,-0.04059,-0.06521,-0.03346};

	if (v>2.627) pb2=0.0626;
	else pb2=exp(polynom(log(v), 4, a));
*/
	double a[4]={630.99766,11510.41807,-2.36954E6,1.09084E8};

	if (fabs(x)<ERRORHI || fabs(x-1)<ERRORHI) 
		return v4ts(v);

	if (v>0.009) t=divid2r(tx4v, v, TC, TMIN, x, ERRORHI);
	else t=parabola2r(tx4v, v, polynom(v,3,a), x, ERRORHI);

	return t;
}

double hx4t(double h, double x)  //multi-value
{
	double t;

	if (fabs(x)<ERRORHI || fabs(x-1)<ERRORHI) 
		return h4ts(h);
	if (h>2500)	t=parabola2r(tx4h, h, TMIN, x, ERRORHI);
	else t=parabola2r(tx4h, h, 513.15, x, ERRORHI);
	if (t<T0) return ERRCALC;

	return t;
}

double sx4t(double s, double x)  //multi-value
{
	if (fabs(x)<ERRORHI || fabs(x-1)<ERRORHI) 
		return s4ts(s);

//	t=parabola2r(tx4s, s, 500, x, ERRORHI);
	return divid2r(tx4s, s, TC, TMIN, x, ERRORHI);
}

double tv4h(double t, double v)
{
	double x;

	x=tv4x(t,v);

	return tx4h(t,x);
}

double vh4t(double v, double h)
{
	double t;
	double a[4]={630.99766,11510.41807,-2.36954E6,1.09084E8};

	if (v>=0.003) t=divid2r(tv4h, h, TC, TMIN, v, ERRORHI);
	else t=parabola2r(tv4h, h, polynom(v,3,a), v, ERRORHI);

	return t;
}

double tv4s(double t, double v)
{
	double x;

	x=tv4x(t,v);

	return tx4s(t,x);
}

double vs4t(double v, double s)    // not confirmed
{
	double t;
	double a[4]={625.135,15341.99447,-3.16612E6,1.61894E8};

	if (v>=0.006) t=divid2r(tv4s, s, TC, TMIN, v, ERRORHI);
	else
	{
	//	t=aitkin2r(tv4s, s, polynom(v,3,a), v, ERRORHI);
		t=divid2r(tv4s, s, polynom(v,3,a), TMIN, v, ERRORHI);
	}

	return t;
}

double vh4x(double v, double h)
{
	double t, x;

	t=vh4t(v,h);
	x=th4x(t,h);

	if (x>1||x<0) return ERRCALC;
	return x;
}

double vs4x(double v, double s)
{
	double t;

	t=vs4t(v,s);

	return ts4x(t,s);
}

double vs4h(double v, double s)
{
	double t, x;

	t=vs4t(v,s);
	x=vs4x(v,s);

	if (x>1||x<0) return ERRCALC;
	return tx4h(t,x);
}

double vh4s(double v, double h)
{
	double t, x;

	t=vh4t(v,h);
	x=vh4x(v,h);

	if (x>1||x<0) return ERRCALC;
	return tx4s(t,x);
}


double hs4v(double h, double s)
{
	double v;

	v=divid2r(vs4h, h, 0.001, 206.14, s, ERRORHI);

	return v;
}

double ph4x(double p, double h)
{
	return th4x(p4t(p), h);
}

/* NIF97 obsolete implementation
double hs4t(double h, double s)
{
	double v;

	v = hs4v(h, s);
	
	return vs4t(v,s);
}*/

double hs4t(double h, double s) /* SatT, Kyoto, Sep 2004 */
{
	double tx = 550;
	double hx = 2800.0;
	double sx = 9.2;
	int I[36] = {0,0,0,1,1, 1,1,2,2,2, 3,3,3,3,4, 4,5,5,5,5,
		6,6,6,8,10, 10,12,14,14,16, 16,18,18,18,20, 28};
	int J[36] = {0,3,12,0,1, 2,5,0,5,8, 0,2,3,4,0, 1,1,2,4,16,
		6,8,22,1,20, 36,24,1,28,12, 32,14,22,36,24, 36};
	double n[36] = {0.179882673606601, -0.267507455199603,
		0.116276722612600e1, 0.147545428713616, -0.512871635973248,
		0.421333567697984, 0.563749522189870, 0.429274443819153,
		-0.335704552142140e1, 0.108890916499278e2, -0.248483390456012,
		0.304153221906390, -0.494819763939905, 0.107551674933261e1,
		0.733888415457688e-1, 0.140170545411085e-1, -0.106110975998808,
		0.168324361811875e-1, 0.125028363714877e1, 0.101316840309509e4,
		-0.151791558000712e1, 0.524277865990866e2, 0.230495545563912e5,
		0.249459806365456e-1, 0.210796467412137e7, 0.366836848613065e9,
		-0.144814105365163e9, -0.179276373003590e-2, 0.489955602100459e10,
		0.471262212070518e3, -0.829294390198652e11, -0.171545662263191e4,
		0.355777682973575e7, 0.586062760258436e12, -0.129887635078195e8,
		0.317247449371057e11};

	double rdt, rdh, rds;

	if (s>=5.210887825) /* Ssat vavor @623.15k - the valid range of new formula */
//	if (s>0)
	{
		rdh = h/hx - 0.119;
		rds = s/sx - 1.07;
		rdt = polyXY2e(rdh, rds, n, I, J, 36);

		return (rdt * tx);
	}
	else		/* the old implementation */
	{
		return divid2r(ts4h, h, 273.15, s4ts(s),s, ERRORHI);
	}
}

/* obsolete implementation *
double hs4x(double h, double s)
{
	double v;

	v = hs4v(h, s);
	
	return vs4x(v,s);
} */

/* Sep. 2004 Kyoto, Japan */
double hs4x(double h, double s)  /* implemented Aug22, 2005 */
{
	double t, x;

	t = hs4t(h, s);
	x = th4x(t, h);

	return x;
}


double v4ts(double v)
{
	double n[7]={-305.8777871, 1519.336161, -1003.41105, 350.9587066, -68.47912521, 7.060373845, -0.300587882};

	if (v<0.002222126971) return divid1(t4vs1, v, 643.5, TMIN, ERRORHI);
	if (v<0.00494619656943761) return polynom(v*1000,6,n);
	return divid1(t4vs2, v, 643.5, TMIN, ERRORHI);
}

double h4ts(double h)
{
	double n[6]={-36862.6161, 81527.00263, -70492.86556, 30263.94426, -6437.885626, 541.2889354};

	if (h<=1892.652303) return divid1(t4hs1, h, T3Aux, TMIN, ERRORHI);
	if (h<=2333.502066) return polynom(h/1000,5,n);
	return divid1(t4hs2, h, T3Aux, TMIN, ERRORHI);
}

double s4ts(double s)
{
	double n[5]={-50224.26746, 45162.52446, -15051.12619, 2231.851405, -124.2550044};

	if (s<=4.114169019) return divid1(t4ss1, s, T3Aux, TMIN, ERRORHI);
	if (s<=4.799621967) return polynom(s, 4, n);
	return divid1(t4ss2, s, T3Aux, TMIN, ERRORHI);
}

double ts4h(double t, double s) /* Aug18-2005 */
{
	return tx4h(t, ts4x(t,s));
}