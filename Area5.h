/**********************************************************************/
/*	     Define Functions for Area 5 in NIF                       */
/**********************************************************************/

double pt5v(double p, double t);
double pt5u(double p, double t);
double pt5s(double p, double t);
double pt5h(double p, double t);
double pt5cp(double p, double t);
double pt5cv(double p, double t);
double pt5w(double p, double t);
/* Derived functions */
double ph5t(double p,double h); /* can be used as BACKWARD function */
double ph5s(double p,double h); /* implemented Sep.1 2005 - BACKWARD */
double ps5t(double p,double s);
double ps5h(double p,double s);

double ts5p(double t, double s);
double ts5h(double t, double s);
double th5p(double t, double h);

double hs5t(double h, double s);
double hs5p(double h, double s); /* implemented Sep.1 2005 */

double pv5t(double p, double v);
double tv5p(double t, double v);

double ph5v(double p, double h);
double vh5p(double v, double h);
double ps5v(double p, double s);
double vs5p(double v, double s);
