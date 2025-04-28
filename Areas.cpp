//  Filename: areas.cpp
//  Function: generic subroutines for areas
#include "stdafx.h"

#include <math.h>
#include "tools.h"
#include "area1.h"
#include "area2.h"
#include "area4.h"
#include "areas.h"
#include "region.h"
 
// constants for boundary between region 2 and 3 -- B23-equation
double B23n[5] = {0.34805185628969e3, -0.11671859879975e1, 0.10192970039326e-2,
                  0.57254459862746e3, 0.13918839778870e2 };

// Gibbs free Energy partial equations
double Gmo(double x, double y, double no[], int Jo[], int rank)
{
	double sgm;
	int i;

	sgm = log(x);

	for (i=0; i<rank; i++)
		sgm += no[i] * xen(y, Jo[i]) ;

	return sgm;
}

double Gmo_p(double x, double y, double no[], int Jo[], int rank)
{
	double sgm;
	
	sgm = 1/x;

	return sgm;
}

double Gmo_pp(double x, double y, double no[], int Jo[], int rank)
{
	double sgm;
	
	sgm = -1/(x*x);

	return sgm;
}

double Gmo_t(double x, double y, double no[], int Jo[], int rank)
{
	double sgm;
	int i;

	sgm = 0;

	for (i=0; i<rank; i++)
		sgm += no[i] * Jo[i]*xen(y, Jo[i]-1) ;

	return sgm;
}

double Gmo_tt(double x, double y, double no[], int Jo[], int rank)
{
	double sgm;
	int i;

	sgm = 0;

	for (i=0; i<rank; i++)
		sgm += no[i] * Jo[i]*(Jo[i]-1)*xen(y, Jo[i]-2) ;

	return sgm;
}

double Gmo_pt(double x, double y, double no[], int Jo[], int rank)
{
	double sgm =0;

	return sgm;
}

double Gmt(double x, double y, double n[], int I[], int J[], int rank)
{
	double sgm;
	int i;

	sgm = 0;

	for (i=0; i<rank; i++)
		sgm += n[i] * xen(x, I[i])*xen(y, J[i]);

	return sgm;
}

double Gmt_p(double x, double y, double n[], int I[], int J[], int rank)
{
	double sgm;
	int i;

	sgm = 0;

	for (i=0; i<rank; i++)
		sgm += n[i]*I[i] * 
		     xen(x, I[i]-1)*xen(y, J[i]);

	return sgm;
}

double Gmt_pp(double x, double y, double n[], int I[], int J[], int rank)
{
	double sgm;
	int i;

	sgm = 0;

	for (i=0; i<rank; i++)
		sgm += n[i]*I[i]*(I[i]-1) * 
		     xen(x, I[i]-2)*xen(y, J[i]);

	return sgm;
}

double Gmt_t(double x, double y, double n[], int I[], int J[], int rank)
{
	double sgm;
	int i;

	sgm = 0;

	for (i=0; i<rank; i++)
		sgm += n[i]*J[i] * 
		     xen(x, I[i])*xen(y, J[i]-1);

	return sgm;
}

double Gmt_tt(double x, double y, double n[], int I[], int J[], int rank)
{
	double sgm;
	int i;

	sgm = 0;

	for (i=0; i<rank; i++)
		sgm += n[i]*J[i]*(J[i]-1) *
		   xen(x, I[i])*xen(y, J[i]-2);

	return sgm;
}

double Gmt_pt(double x, double y, double n[], int I[], int J[], int rank)
{
	double sgm;
	int i;

	sgm = 0;

	for (i=0; i<rank; i++)
		sgm += n[i]*I[i]*J[i] *
		   xen(x, I[i]-1)*xen(y, J[i]-1);

	return sgm;
}

/* saturation line for region purpose */
double t_sat_p(double p)
{
	return p4t(p);
}

double p_sat_t(double t)
{
	return t4p(t);
}


// boundary between region 2 and 3 -- B23-equation
double p_B23_t(double t)
{
	double sgm;
	double St = t/1.0;

	sgm = B23n[0] + B23n[1]*St + B23n[2]*St*St;

	return (1.0*sgm);
}

double t_B23_p(double p)
{
	double sgm;
	double Pi = p/1.0;

	sgm = B23n[3] + sqrt((Pi -B23n[4])/B23n[2]);

	return (1.0*sgm);
}

// area identification

double sbd4p(double s)
{
	double a1[4]={3.4217,2.53803,-0.48615,0.01655};

	return exp(polynom(s,3,a1));
}

double sbd3t(double s)
{
	double a1[5]={2.65502,44.84736,35.68935,-8.188,0.76934};

	return (polynom(s,4,a1)+T0);
}

double hbd3v(double h) /* border3 in area3 */
{
	double a1[5]={-6.9706,1.44872E-4,-5.41643E-8,8.48005E-11,-1.3693E-14};

	return exp(polynom(h,4,a1));
}

double vbd3h(double v)
{
	double a1[5]={-1.90375E6,-1.29205E6,-326915.64349,-36650.98209,-1538.20245};
	return polynom(log(v), 4, a1);
}

double sbd3v(double s)
{
	return 9.13439E-4+2.80822E-5*exp(s/1.25289);
}
/*
double vbd8s(double v)
{
	double a1[4]={5.2125,-157.82255,40906.1903,-2.5557E6};
	double a2[2][5]={{8.80066E-3,8.75794E-3,8.63202E-3,8.44859E-3,
		8.22477E-3},{5.20996,5.23473,5.24971,5.25774,5.26054}};
	
	if (v<0.003286)
		return 5.04242+0.05607*exp(-(v-0.00259)/2.81042E-4);
	if (v<0.008225)
		return polynom(v,3,a1);
	return interpolate2(v, 5, &a2[0][0]);
}
*/

/* boundary of area4 - NIF97 obsolete implementation
double h_B4_s(double s)
{
	double h;

	if (s<4.412) h = t4hs1(sx4t(s,0));
	else h = t4hs2(sx4t(s,1));

	return h;
}
*/

/* boundary of area 4 h(s) Kyoto, Japan, Sep.2004 */
double h_satB1_s(double s)
{
	double hx = 1700;
	double sx = 3.8;

	int I[27] = {0,0,1,1,2, 2,3,3,4,4, 4,5,5,7,8,
		12,12,14,14,16, 20,20,22,24,28, 32,32};
	int J[27] = {14,36,3,16,0, 5,4,36,4,16, 24,18,24,1,4,
		2,4,1,22,10, 12,28,8,3,0, 6,8};
	double n[27] = {0.332171191705237, 0.611217706323496e-3,
		-0.882092478906822e1, -0.455628192543250, -0.263483840850452e-4,
		-0.223949661148062e2, -0.428398660164013e1, -0.616679338856916,
		-0.146823031104040e2, 0.284523138727299e3, -0.113398503195444e3,
		0.115671380760859e4, 0.395551267359325e3, -0.154891257229285e1,
		0.194486637751291e2, -0.357915139457043e1, -0.335369414148819e1,
		-0.664426796332460, 0.323321885383934e5, 0.331766744667084e4,
		-0.223501257931087e5, 0.573953875852936e7, 0.173226193407919e3,
		-0.363968822121321e-1, 0.834596332878346e-6, 0.503611916682674e1,
		0.655444787064505e2};

	double rds1, rds2, rdh;

	rds1 = s/sx -1.09;
	rds2 = s/sx +0.366e-4;
	rdh = polyXY2e(rds1, rds2, n, I, J, 27);

	return (rdh * hx);
}

double h_satB3a_s(double s)
{
	double hx = 1700;
	double sx = 3.8;

	int I[19] = {0,0,0,0,2, 3,4,4,5,5, 6,7,7,7,10,
		10,10,32,32};
	int J[19] = {1,4,10,16,1, 36,3,16,20,36, 4,2,28,32,14,
		32,36,0,6};
	double n[19] = {0.822673364673336, 0.181977213534479,
		-0.112000260313624e-1, -0.746778287048033e-3, -0.179046263257381,
		0.424220110836657e-1, -0.341355823438768, -0.209881740853565e1,
		-0.822477343323596e1, -0.499684082076008e1, 0.191413958471069,
		0.581062241093136e-1, -0.165505498701029e4, 0.158870443421201e4,
		-0.850623535172818e2, -0.317714386511207e5, -0.945890406632871e5,
		-0.139273847088690e-5, 0.631052532240980};

	double rds1, rds2, rdh;

	rds1 = s/sx -1.09;
	rds2 = s/sx +0.366e-4;
	rdh = polyXY2e(rds1, rds2, n, I, J, 19);

	return (rdh * hx);
}

double h_satB2c3b_s(double s)
{
	double hx = 2800;
	double sx = 5.9;

	int I[16] = {0,0,0,1,1, 5,6,7,8,8, 12,16,22,22,24, 36};
	int J[16] = {0,3,4,0,12, 36,12,16,2,20, 32,36,2,32,7, 20};
	double n[16] = {0.104351280732769e1, -0.227807912708513e1,
		0.180535256723202e1, 0.420440834792042, -0.105721244834660e6,
		0.436911607493884e25, -0.328032702839753e12, -0.678686760804270e16,
		0.743957464645363e4, -0.356896445355761e20, 0.167590585186801e32,
		-0.355028625419105e38, 0.396611982166538e12, -0.414716268484468e41,
		0.359080103867382e19, -0.116994334851995e41};

	double rds1, rds2, rdh;

	rds1 = s/sx -1.02;
	rds2 = s/sx -0.726;
	rdh = polyXY2e(rds1, rds2, n, I, J, 16);

	return (xen(rdh,4) * hx);
}

double h_satB2ab_s(double s)
{
	double hx = 2800;
	double sx1 = 5.21;
	double sx2 = 9.2;

	int I[30] = {1,1,2,2,4, 4,7,8,8,10, 12,12,18,20,24,
		28,28,28,28,28, 32,32,32,32,32, 36,36,36,36,36};
	int J[30] = {8,24,4,32,1, 2,7,5,12,1, 0,7,10,12,32,
		8,12,20,22,24, 2,7,12,14,24, 10,12,20,22,28};
	double n[30] = {-0.524581170928788e3, -0.926947218142218e7,
		-0.237385107491666e3, 0.210770155812776e11, -0.239494562010986e2,
		0.221802480294197e3, -0.510472533393438e7, 0.124981396109147e7,
		0.200008436996201e10, -0.815158509791035e3, -0.157612685637523e3,
		-0.114200422332791e11, 0.662364680776872e16, -0.227622818296144e19,
		-0.171048081348406e32, 0.660788766938091e16, 0.166320055886021e23,
		-0.218003784381501e30, -0.787276140295618e30, 0.151062329700346e32,
		0.795732170300541e7, 0.131957647355347e16, -0.325097068299140e24,
		-0.418600611419248e26, 0.297478906557467e35, -0.953588761745473e20,
		0.166957699620939e25, -0.175407764869978e33, 0.347581490626396e35,
		-0.710971318427851e39};

	double rds1, rds2, rdh;

	rds1 = sx1/s -0.513;
	rds2 = s/sx2 - 0.524;
	rdh = polyXY2e(rds1, rds2, n, I, J, 30);

	return (exp(rdh) * hx);
}

double h_B4_s(double s)
{
	/*if (s<-1.545495919e-4)   S'(273.15k) 
		return S_OUTRANGE;*/
	if (s<3.778281340)  /* S'(623.15k) */
		return h_satB1_s(s);
	if (s<4.41202148223476) return h_satB3a_s(s);
	if (s<5.85) return h_satB2c3b_s(s);
	if (s<=9.155759395) return h_satB2ab_s(s);
	else return h_satB2ab_s(s); /* in case */
}

/* boundary of B34 p(h) - Kyoto, Japan, Sep.2004 */
double p_B34_h(double h) /* Implemented Aug22,2005 */
{
	double px = 22;
	double hx = 2600;

	int I[14] = {0,1,1,1,1, 5,7,8,14,20, 22,24,28,36};
	int J[14] = {0,1,3,4,36, 3,0,24,16,16, 3,18,8,24};
	double n[14] = {0.600073641753024, -0.936203654849857e1,
		0.246590798594147e2, -0.107014222858224e3, -0.915821315805768e14,
		-0.862332011700662e4, -0.235837344740032e2, 0.252304969384128e18,
		-0.389718771997719e19, -0.333775713645296e23, 0.356499469636328e11,
		-0.148547544720641e27, 0.330611514838798e19, 0.813641294467829e38};

	double rdh1, rdh2, rdp;

	rdh1 = h/hx -1.02;
	rdh2 = h/hx -0.608;
	rdp = polyXY2e(rdh1, rdh2, n, I, J, 14);

	return (rdp * px);
}

/* boundary of B34 p(s) Kyoto, Japan, Sep.2004 */
double p_B34_s(double s) /* Implemented Aug22,2005 */
{
	double px = 22;
	double sx = 5.2;

	int I[10] = {0,1,1,4,12, 12,16,24,28,32};
	int J[10] = {0,1,32,7,4, 14,36,10,0,18};
	double n[10] = {0.639767553612785, -0.129727445396014e2,
		-0.224595125848403e16, 0.177466741801846e7, 0.717079349571538e10,
		-0.378829107169011e18, -0.955586736431328e35, 0.187269814676188e24,
		0.119254746466473e12, 0.110649277244882e37};

	double rds1, rds2, rdp;

	rds1 = s/sx -1.03;
	rds2 = s/sx -0.699;
	rdp = polyXY2e(rds1, rds2, n, I, J, 10);

	return (rdp * px);
}

/* boundary of B13 h_B13(s) - Kyoto, Japan, Sep.2004  */
double h_B13_s(double s) /* Implemented Aug29,2005 */
{
	double hx = 1700;
	double sx = 3.8;

	int I[6] = {0,1,1,3,5, 6};
	int J[6] = {0,-2,2,-12,-4, -3};
	double n[6] = {0.913965547600543, -0.430944856041991e-4,
		0.603235694765419e2, 0.117518273082168e-17, 0.220000904781292,
		-0.690815545851641e2};

	double rds1, rds2, rdh;

	rds1 = s/sx -0.884;
	rds2 = s/sx -0.864;
	rdh = polyXY2e(rds1, rds2, n, I, J, 6);

	return (rdh * hx);
}

double t_B23_hs(double h, double s) /* Implemented Aug29,2005 */
{
	double tx = 900;
	double hx = 3000;
	double sx = 5.3;

	int I[25] = {-12,-10,-8,-4,-3, -2,-2,-2,-2,0, 1,1,1,3,3,
		5,6,6,8,8,	8,12,12,14,14};
	int J[25] = {10,8,3,4,3, -6,2,3,4,0, -3,-2,10,-2,-1,
		-5,-6,-3,-8,-2, -1,-12,-1,-12,1};
	double n[25] = {0.629096260829810e-3, -0.823453502583165e-3,
		0.515446951519474e-7, -0.117565945784945e1, 0.348519684726192e1,
		-0.507837382408313e-11, -0.284637670005479e1, -0.236092263939673e1,
		0.601492324973779e1, 0.148039650824546e1, 0.360075182221907e-3,
		-0.126700045009952e-1, -0.122184332521413e7, 0.149276502463272,
		0.698733471798484, -0.252207040114321e-1, 0.147151930985213e-1,
		-0.108618917681849e1, -0.936875039816322e-3, 0.819877897570217e2,
		-0.182041861521835e3, 0.261907376402688e-5, -0.291626417025961e5,
		0.140660774926165e-4, 0.783237062349385e7};

	double rdh, rds, rdt;

	rdh = h/hx -0.727;
	rds = s/sx -0.864;
	rdt = polyXY2e(rdh, rds, n, I, J, 25);

	return (rdt * tx);
}
