/**********************************************************************/
/*	     Define Functions for Area 5 in NIF                       */
/**********************************************************************/
#include "stdafx.h"

#include <math.h>
#include "areas.h"
#include "area5.h"
#include "tools.h"

#define PX5 1.0       /* MPa */
#define TX5 1000.0    /* K */
#define VX5 0.001    /* m^3kg^-1  */
#define W2X5 1000.0   /* m^2s^-2 */

static double GE5no[6]  =
  { -0.13179983674201e2, 0.68540841634434e1 ,-0.24805148933466e-1,
    0.36901534980333e0, -0.31161318213925e1, -0.32961626538917e0,
  };

static int GE5Jo[6]  =
  { 0, 1, -3, -2, -1, 2
  };

/* version 97
static double GE5n[5]  =
  { -0.12563183589592e-3, 0.21774678714571e-2, -0.45942820899910e-2,
    -0.39724828359569e-5, 0.12919228289784e-6
  };

static int GE5I[5]  =
  { 1, 1, 1, 2, 3
  };

static int GE5J[5]  =
  { 0, 1, 3, 9, 3
  };

  */


static double GE5n[6]  =
  { 0.15736404855259e-2, 0.90153761673944e-3, -0.50270077677648e-2,
	0.22440037409485e-5, -0.41163275453471e-5, 0.37919454822955e-7
  };

static int GE5I[6]  =
  { 1, 1, 1, 2, 2, 3
  };

static int GE5J[6]  =
  { 1, 2, 3, 3, 9, 7
  };

/* coeff for t=f(h) */
double n_t5_h[4]={-1106.99307, 0.66613, -4.06334E-5, 1.68873E-9};


/* specific volume */
double pt5v(double p, double t)
{
	double sgm;
	double rdp = p / PX5;
	double rdt = TX5 / t;
	double Gmo5_p= Gmo_p(rdp, rdt, GE5no, GE5Jo, 6);
	double Gmt5_p= Gmt_p(rdp, rdt, GE5n, GE5I, GE5J, 6);

	sgm = rdp *(Gmo5_p+Gmt5_p);
	sgm *= (R*t/p); 

	return (sgm*VX5);
}

double pt5u(double p, double t)
{
	double sgm;
	double rdp = p / PX5;
	double rdt = TX5 / t;
	double Gmo5_t= Gmo_t(rdp, rdt, GE5no, GE5Jo, 6);
	double Gmt5_t= Gmt_t(rdp, rdt, GE5n, GE5I, GE5J, 6);
	double Gmo5_p= Gmo_p(rdp, rdt, GE5no, GE5Jo, 6);
	double Gmt5_p= Gmt_p(rdp, rdt, GE5n, GE5I, GE5J, 6);

	sgm = rdt *(Gmo5_t+Gmt5_t) - rdp*(Gmo5_p+Gmt5_p);
    sgm *= (R*t);
	
	return sgm;
}

double pt5s(double p, double t)
{
	double sgm;
	double rdp = p / PX5;
	double rdt = TX5 / t;
	double Gmo5_t= Gmo_t(rdp, rdt, GE5no, GE5Jo, 6);
	double Gmt5_t= Gmt_t(rdp, rdt, GE5n, GE5I, GE5J, 6);
	double Gmo5= Gmo(rdp, rdt, GE5no, GE5Jo, 6);
	double Gmt5= Gmt(rdp, rdt, GE5n, GE5I, GE5J, 6);

	sgm = (R)* (rdt *(Gmo5_t+Gmt5_t) - (Gmo5+Gmt5));

	return sgm;
}

double pt5h(double p, double t)
{
	double sgm;
	double rdp = p / PX5;
	double rdt = TX5 / t;
	double Gmo5_t= Gmo_t(rdp, rdt, GE5no, GE5Jo, 6);
	double Gmt5_t= Gmt_t(rdp, rdt, GE5n, GE5I, GE5J, 6);

	sgm = rdt *(Gmo5_t+Gmt5_t);
	sgm *= (R*t);

	return sgm;
}

double pt5cp(double p, double t)
{
	double sgm;
	double rdp = p / PX5;
	double rdt = TX5 / t;
	double Gmo5_tt= Gmo_tt(rdp, rdt, GE5no, GE5Jo, 6);
	double Gmt5_tt= Gmt_tt(rdp, rdt, GE5n, GE5I, GE5J, 6);
	
	sgm = -rdt*rdt *(Gmo5_tt+Gmt5_tt);
	sgm *= (R); 

	return sgm;
}    

double pt5cv(double p, double t)
{
	double sgm;
	double rdp = p / PX5;
	double rdt = TX5 / t;
	double Gmo5_tt= Gmo_tt(rdp, rdt, GE5no, GE5Jo, 6);
	double Gmt5_tt= Gmt_tt(rdp, rdt, GE5n, GE5I, GE5J, 6);
	double Gmt5_p= Gmt_p(rdp, rdt, GE5n, GE5I, GE5J, 6);
	double Gmt5_pp= Gmt_pp(rdp, rdt, GE5n, GE5I, GE5J, 6);
	double Gmt5_pt= Gmt_pt(rdp, rdt, GE5n, GE5I, GE5J, 6);

	sgm = -rdt*rdt *(Gmo5_tt+Gmt5_tt);
	sgm += -xen(1+rdp*(Gmt5_p-rdt*Gmt5_pt),2) /
		(1- rdp*rdp*Gmt5_pp);
	sgm *= (R); 

	return sgm;
}

double pt5w(double p, double t)
{
	double sgm;
	double rdp = p / PX5;
	double rdt = TX5 / t;
	double Gmt5_p= Gmt_p(rdp, rdt, GE5n, GE5I, GE5J, 6);
	double Gmo5_tt= Gmo_tt(rdp, rdt, GE5no, GE5Jo, 6);
	double Gmt5_tt= Gmt_tt(rdp, rdt, GE5n, GE5I, GE5J, 6);
	double Gmt5_pp= Gmt_pp(rdp, rdt, GE5n, GE5I, GE5J, 6);
	double Gmt5_pt= Gmt_pt(rdp, rdt, GE5n, GE5I, GE5J, 6);

	sgm = 1+rdp*Gmt5_p*(2+rdp*Gmt5_p);
	sgm /= ((1- rdp*rdp*Gmt5_pp) +
		xen(1+rdp*(Gmt5_p-rdt*Gmt5_pt),2) /
		(rdt*rdt*(Gmo5_tt+Gmt5_tt)));
	sgm *= (R*t); 

	return sqrt(sgm*W2X5);
}


/* backward functions */
double ph5t(double p,double h)
{
   double t0;
 
   /* in region 5, h mainly depend on t */
   t0=polynom(h, 3, n_t5_h);

   return aitkin2(pt5h, h, p, t0,  ERRORHI);
}

double ph5s(double p,double h) /* implemented Sep.1 2005 */
{
	return pt5s(p, ph5t(p,h));
}

double ps5t(double p,double s)
{
	double t;

	t = parabola2(pt5s, s, p, 1675, ERRORHI);

	return t;
}

double ps5h(double p,double s)
{
	double t, h;

	t = ps5t(p, s);
	h = pt5h(p, t);

	return h;
}


double ts5p(double t, double s)
{
	double pb1, p;

	pb1= 5.12811*exp(-(s-7.73853)/0.47525) -0.04874*exp(-(s-7.73853)/1.15405);
	p = divid2r(pt5s, s, PB50, pb1, t, ERRORHI2);

	return p;
}
	
double ts5h(double t, double s)
{
	double p, h;

	p = ts5p(t, s);
	h = pt5h(p, t);

	return h;
}

double th5p(double t, double h)
{
	double p;

	p = divid2r(pt5h, h, PB50, PB0, t, ERRORHI3);

	return p;
}


double hs5t(double h, double s)
{
	return ph5t(hs5p(h,s), h);
}

double hs5p(double h, double s)  /* implemented Sep.1 2005 */
{
	/* iterate s_5_ph() which can be used as backward */
	return divid2r(ph5s, s, PMIN, PB50, h, ERRORHI);

}

double pv5t(double p, double v)
{
	double t;

	// t0 = p*v/0.00046;
	// t = parabola2(pt5v, v, p, t0, ERRORHI);
	t = divid2(pt5v, v, p, TB4, TB3, ERRORHI);

	return t;
}

double tv5p(double t, double v)
{
	double p;

	p = divid2r(pt5v, v, PB50, PB0, t, ERRORHI);

	return p;
}

double ph5v(double p, double h)
{
	double t;

	t=ph5t(p,h);
	return (pt5v(p,t));
}

double vh5p(double v, double h)
{
	double p;

	p=divid2r(ph5v, v, PB0, PB50, h, ERRORHI2);

	return p;
}

double ps5v(double p, double s)
{
	double t;

	t=ps5t(p,s);
	return (pt5v(p,t));
}

double vs5p(double v, double s)
{
	double p, pb1;

//	pb1=5.00198*exp(-(s-7.7459)/0.47048)+0.01; 
	pb1= 5.12811*exp(-(s-7.73853)/0.47525) -0.04874*exp(-(s-7.73853)/1.15405);
	p=divid2r(ps5v, v, PB50, pb1, s, ERRORHI);

	return p;
}