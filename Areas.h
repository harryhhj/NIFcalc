/*  Filename: areas.h
 *  Function: generic subroutines for areas
*/

#define T0 273.15      /* unit: K. celcius 0 degree */
#define R 0.461526     /* unit: KJKg^-1K^-1 */
#define TC 647.096      /* unit: K */
#define PC 22.064      /* unit: MPa */
#define VC 0.003106    /* unit: m^3/kg */
#define SC 4.41202148223476
#define HC 2087.55
#define ROC 322.0     /* unit: Kgm^-3 */
#define Tt 273.16     /* unit: K. Triple point where u'=s'=0, Tt, and Pt */
#define Pt 611.657e-6	  /* unit: MPa Triple point where u'=s'=0, Tt, and Pt */
#define Ht 0.611783   /* Hmin, Triple point where u'=s'=0, Tt, and Pt */ 
#define M 18.015257   /* J mol^-1 */

#define TMIN	T0
#define PMIN	6.11212677e-4	/* Ps(t=273.15) */
#define HMIN	-4.1587826e-2   /* (p= PMIN, t=TMIN) */
#define SMIN	-8.582287e-3	/* (p= PMAX, t=TMIN) */
#define PMAX	100
#define TMAX	2273.15
#define SMAX    13.904956091	/* (PMIN, TB4)  */

#define TB2  623.15
#define PSB2 16.52916 /* ps(TB2) */
#define TB3  1073.15
#define TB4  2273.15
#define TB23 863.15
#define PB0		PMIN
//#define PB10  10.0
#define PB50  50.0
#define PB100  100.0      
#define PB4  4.0        /* p=4MPa*/
#define SBC  5.85   /* boundary between 2b & 2c is iso-entropy line*/
#define H1_B34 1.670858218e3 /* H'sat(t=623.15k) */
#define H2_B34 2.563592004e3 /* H"sat(t=623.15k) */
#define S1_B34 3.778323559 /* H'sat(t=623.15k) */
#define S2_B34 5.210951570 /* H"sat(t=623.15k) */
#define S1_B14 -1.545495919e-4  /* S_sat(t=273.15) */

// Gibbs free Energy partial equations
double Gmo(double x, double y, double no[], int Jo[], int rank);
double Gmo_p(double x, double y, double no[], int Jo[], int rank);
double Gmo_pp(double x, double y, double no[], int Jo[], int rank);
double Gmo_t(double x, double y, double no[], int Jo[], int rank);
double Gmo_tt(double x, double y, double no[], int Jo[], int rank);
double Gmo_pt(double x, double y, double no[], int Jo[], int rank);

double Gmt(double x, double y, double n[], int I[], int J[], int rank);
double Gmt_p(double x, double y, double n[], int I[], int J[], int rank);
double Gmt_pp(double x, double y, double n[], int I[], int J[], int rank);
double Gmt_t(double x, double y, double n[], int I[], int J[], int rank);
double Gmt_tt(double x, double y, double n[], int I[], int J[], int rank);
double Gmt_pt(double x, double y, double n[], int I[], int J[], int rank);

/* saturation line for region purpose */
double t_sat_p(double p);
double p_sat_t(double t);

// boundary between region 2 and 3 -- B23-equation
double p_B23_t(double t);
double t_B23_p(double p);

/* describe the border of the NIF */

/* #1 = PB1 (TB1 -> TB4) */
/* #2 = TB1 (PB1 -> PB3) */
/* #3 = PB3 (TB1 -> TB3) */
/* #4 = TB3 (PB1 -> PB3) */
/* #5 = PB2 (TB3 -> TB4) */
/* #6 = TB4 (PB1 -> PB2) */
/* #7 = TB2 (Ps  -> PB3) */
/* #8 = Between area 2 & 3 */
/* if it is necessary, function maybe divided into several sector inside */

double sbd4p(double s);
double sbd3t(double s);
double hbd3v(double h);
double vbd3h(double v);
double sbd3v(double s);
double vbd8s(double v);

/* boundary of area 4 - 2004 formula implementation */
double h_satB1_s(double s);  /* boundary of area 1,4 h(s) Kyoto, Japan, Sep.2004 */
double h_satB3a_s(double s); /* boundary of area 3,4 h(s) below Critical Kyoto, Japan, Sep.2004 */

double h_satB2c3b_s(double s);
double h_satB2ab_s(double s);
double h_B4_s(double s); /*includes all new h(s) boundary formula in Kyoto 2004 release */
double p_B34_h(double h); /* Implemented Aug22,2005 */
double p_B34_s(double s); /* Implemented Aug22,2005 */

double h_B13_s(double s); /* Implemented Aug29,2005 */
double t_B23_hs(double h, double s); /* Implemented Aug29,2005 */
