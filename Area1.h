/**********************************************************************/
/*	    Filename: area1.h                                           */
/*		Function: Define Functions for Region 1 in NIF                */
/**********************************************************************/


double pt1v(double p,double t);
double pt1u(double p,double t);
double pt1s(double p,double t);
double pt1h(double p,double t);
double pt1cp(double p,double t);
double pt1cv(double p,double t);
double pt1w(double p,double t);

// backward equations
double ph1t(double p,double h);
double ps1t(double p,double s);

double ts1p(double t, double s);
double th1p(double t, double s);
double th1s(double t, double h);
double ts1h(double t, double s);
double hs1t(double h, double s);
double hs1p(double h, double s);  /* Gaithersburg Maryland USA Sep 2001 */

double pv1t(double p,double v);
double tv1p(double t, double v);
double ph1v(double p, double h);
double vh1p(double v, double h);
double ps1v(double p, double s);
double vs1p(double v, double s);

/* derived */
double ps1h(double p, double s); /* h_ps(p, t_ps(p,s)) */
