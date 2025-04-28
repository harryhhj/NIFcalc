/**********************************************************************/
/*	      Filename: area2.cpp                                         */
/*		  Function: Define Functions for Region 2 in NIF              */
/**********************************************************************/
#include "stdafx.h"

#include <math.h>
#include "areas.h"
#include "area2.h"
#include "area4.h"
#include "tools.h"

#define TX2 540        /* K */
#define PX2 1.0        /* MPa */
#define HX2 2000.0     /* kJkg^-1 */
#define VX2 0.001    /* m^3kg^-1  */
#define W2X2 1000.0   /* m^2s^-2 */


// coefficient of Gibs free energy ideal gas 
static double GE2no[9]  =
  {
   -0.96927686500217e01,  0.10086655968018e02,-0.56087911283020e-02,
    0.71452738081455e-01,-0.40710498223928e00, 0.14240819171444e01,
   -0.43839511319450e01, -0.28408632460772e00, 0.21268463753307e-01
  };
static int GE2Jo[9]  =
  { 0,1,-5,-4,-3,-2,-1,2,3};

// coefficient of Gibs free energy residual 
static double GE2n[43]  =
  {
   -0.17731742473213e-02,-0.17834862292358e-01,-0.45996013696365e-01,
   -0.57581259083432e-01,-0.50325278727930e-01,-0.33032641670203e-04,
   -0.18948987516315e-03,-0.39392777243355e-02,-0.43797295650573e-01,
   -0.26674547914087e-04,0.20481737692309e-07,0.43870667284435e-06,
   -0.32277677238570e-04,-0.15033924542148e-02,-0.40668253562649e-01,
   -0.78847309559367e-09,0.12790717852285e-07,0.48225372718507e-06,
   0.22922076337661e-05,-0.16714766451061e-10,-0.21171472321355e-02,
   -0.23895741934104e02,-0.59059564324270e-17,-0.12621808899101e-05,
   -0.38946842435739e-01, 0.11256211360459e-10,-0.82311340897998e01,
   0.19809712802088e-07, 0.10406965210174e-18,-0.10234747095929e-12,
   -0.10018179379511e-08,-0.80882908646985e-10,0.10693031879409,
   -0.33662250574171,0.89185845355421e-24,0.30629316876232e-12,
   -0.42002467698208e-05,-0.59056029685639e-25,0.37826947613457e-05,
   -0.12768608934681e-14,0.73087610595061e-28,0.55414715350778e-16,
   -0.94369707241210e-06
  };
static int  GE2I[43]  =
  { 1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,5,6,6,6,7,7,7,8,8,9,10,10,10,
    16,16,18,20,20,20,21,22,23,24,24,24
  };
static int GE2J[43]  =
  { 0,1,2,3,6,1,2,4,7,36,0,1,3,6,35,1,2,3,7,3,16,35,0,11,25,8,36,13,4,10,
  14,29,50,57,20,35,48,21,53,39,26,40,58};

// Supplementary Equation for Metastable-Vapor Region
static double sGE2no[9]  =
  {
   -0.96937268393049e01,  0.10087275970006e02,-0.56087911283020e-02,
    0.71452738081455e-01,-0.40710498223928e00, 0.14240819171444e01,
   -0.43839511319450e01, -0.28408632460772e00, 0.21268463753307e-01
  };
static double sGE2n[13]  =
  { -0.73362260186506e-02,-0.88223831943146e-01,-0.72334555213245e-01,
    -0.40813178534455e-02, 0.20097803380207e-02,-0.53045921898642e-01,
    -0.76190409086970e-02,-0.63498037657313e-02,-0.86043093028588e-01,
     0.75321581522770e-02,-0.79238375446139e-02,-0.22888160778447e-03,
    -0.26456501482810e-02
  };
static int sGE2I[13]  =
  { 1,1,1,1,2,2,2,3,3,4,4,5,5};
static int sGE2J[13]  =
  { 0,2,5,11,1,7,16,4,16,7,10,9,10};

// Backward equations
  // 2b-2c subregion boundary
static double B2bcn[5]  =
  { 0.90584278514723e03,-0.67955786399241e00,0.12809002730136e-03,
    0.26526571908428e04, 0.45257578905948e01
  };

  // 2a  p,h -> t 
static double B2an[34] =
 { 0.10898952318288e04,0.84951654495535e03, -0.10781748091826e03,
   0.33153654801263e02,-0.74232016790248e01, 0.11765048724356e02,
   0.18445749355790e01,-0.41792700549624e01, 0.62478196935812e01,
  -0.17344563108114e02,-0.20058176862096e03, 0.27196065473796e03,
  -0.45511318285818e03, 0.30919688604755e04, 0.25226640357872e06,
  -0.61707422868339e-02,-0.31078046629583e00, 0.11670873077107e02,
   0.12812798404046e09,-0.98554909623276e09, 0.28224546973002e10,
  -0.35948971410703e10, 0.17227349913197e10,-0.13551334240775e05, 
   0.12848734664650e08, 0.13865724283226e01, 0.23598832556514e06,
  -0.13105236545054e08, 0.73999835474766e04,-0.55196697030060e06,
   0.37154085996233e07, 0.19127729239660e05,-0.41535164835634e06,
  -0.62459855192507e02
 };
static int B2aI[34]  =
 { 0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,4,4,4,5,5,5,6,6,7
 };
static int B2aJ[34]  =
 { 0,1,2,3,7,20,0,1,2,3,7,9,11,18,44,0,2,7,36,38,40,42,44,24,44,
   12,32,44,32,36,42,34,44,28
 };
  // 2b  p,h -> t 
static double B2bn[38]   =
 { 0.14895041079516e04,  0.74307798314034e03, -0.97708318797837e02,
   0.24742464705674e01, -0.63281320016026e00,  0.11385952129658e01,
  -0.47811863648625e00,  0.85208123431544e-02,  0.93747147377932e00,
   0.33593118604916e01,  0.33809355601454e01,  0.16844539671904e00,
   0.73875745236695e00, -0.47128737436186e00,  0.15020273139707e00,
  -0.21764114219750e-02,-0.21810755324761e-01,-0.10829784403677e00,
  -0.46333324635812e-01, 0.71280351959551e-04, 0.11032831789999e-03,
   0.18955248387902e-03, 0.30891541160537e-02, 0.13555504554949e-02,
   0.28640237477456e-06,-0.10779857357512e-04,-0.76462712454814e-04,
   0.14052392818316e-04,-0.31083814331434e-04,-0.10302738212103e-05,
   0.28217281635040e-06, 0.12704902271945e-05, 0.73803353468292e-07,
  -0.11030139238909e-07,-0.81456365207833e-13,-0.25180545682962e-10,
  -0.17565233969407e-17, 0.86934156344163e-14
  };
static int B2bI[38]  =
  { 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4,4,5,5,5,6,7,7,9,9
  };
static int B2bJ[38] =
  { 0,1,2,12,18,24,28,40,0,2,6,12,18,24,28,40,2,8,18,40,1,2,12,24,2,12,
    18,24,28,40,18,24,40,28,2,28,1,40
  };
  // 2c  p,h -> t 
static double Bh2cn[23]   =
{ -0.32368398555242e13,  0.73263350902181e13,  0.35825089945447e12,
  -0.58340131851590e12, -0.10783068217470e11,  0.20825544563171e11,
   0.61074783564516e06,  0.85977722535580e06, -0.25745723604170e05,
   0.31081088422714e05, 0.12082315865936e04,  0.48219755109255e03,
   0.37966001272486e01, -0.10842984880077e02, -0.45364172676660e-01,
   0.14559115658698e-12, 0.11261597407230e-11,-0.17804982240686e-10,
   0.12324579690832e-06,-0.11606921130984e-05, 0.27846367088554e-04,
  -0.59270038474176e-03, 0.12918582991878e-02
};
static int Bh2cI[23]  =
{ -7,-7,-6,-6,-5,-5,-2,-2,-1,-1,0,0,1,1,2,6,6,6,6,6,6,6,6};
static int Bh2cJ[23]  =
{0,4,0,2,0,2,0,1,0,2,0,1,4,8,4,0,1,4,10,12,16,20,22};
  // 2a  p,s -> t 
static double Bs2an[46] =
{ -0.39235983861984e06,0.51526573827270e06,0.40482443161048e05,
  -0.32193790923902e03,0.96961424218694e02,-0.22867846371773e02,
  -0.44942914124357e06,-0.50118336020166e04,0.35684463560015,
  0.44235335848190e05,-0.13673388811708e05,0.42163260207864e06,
  0.22516925837475e05,0.47442144865646e03,-0.14931130797647e03,
  -0.19781126320452e06,-0.23554399470760e05,-0.19070616302076e05,
  0.55375669883164e05,0.38293691437363e04,-0.60391860580567e03,
  0.19363102620331e04,0.42660643698610e04,-0.59780636672718e04,
  -0.70401463926862e03,0.33836784107553e03,0.20862786635187e02,
  0.33834172656196e-01,-0.43124428414893e-04,0.16653791356412e03,
  -0.13986292055898e03,-0.78849547999872,0.72132411753872e-01,
  -0.59754839398283e-02,-0.12141358953904e-04,0.23227096733871e-06,
  -0.10538463566194e02,0.20718925496502e01,-0.72193155260427e-01,
  0.20749887081120e-06,-0.18340657911379e-01,0.290362723486965e-06,
  0.21037527893619e00,0.25681239729999e-03,-0.12799002933781e-01,
  -0.82198102652018e-05
};
/* real parameter should be 0.25*Bs2aI[i]  */
static int  Bs2aI[46] =
 { -6, -6, -6, -6, -6, -6, -5, -5, -5, -4,
   -4, -4, -4, -4, -4, -3, -3, -2, -2, -2,
   -2, -1, -1, -1, -1, 1, 1, 1, 1, 2,
   2, 2, 2, 2, 2, 2, 3, 3, 3, 3,
   4, 4, 5, 5, 6, 6
};
static int  Bs2aJ[46] =
{ -24, -23, -19, -13, -11, -10, -19, -15, -6, -26,
  -21, -17,  -16,  -9,  -8, -15, -14, -26, -13, -9,
   -7, -27, -25, -11,  -6,   1,   4,   8,  11,  0,
    1,   5,   6,  10,  14,  16,   0,   4,   9, 17,
    7,  18,   3,  15,   5,  18
};
  // 2b  p,s -> t 
static double Bs2bn[44] =
{ 0.31687665083497e6, 0.20864175881858e2, -0.39859399803599e6,
  -0.21816058518877e2, 0.22369785194242e6, -0.27841703445817e04,
  0.99207436071480e1, -0.75197512299157e5, 0.29708605951158e4,
  -0.34406878548526e1, 0.38815564249115e0, 0.17511295085750e5,
  -0.14237112854449e4, 0.10943803364167e01,0.89971619308495e0,
  -0.33759740098958e4, 0.47162885818355e3, -0.19188241993679e1,
  0.41078580492196e0, -0.33465378172097e0, 0.13870034777505e4,
  -0.40663326195838e3, 0.41727347159610e2, 0.21932549434532e1,
  -0.10320050009077e1, 0.35882943516703e0, 0.52511453726066e-2,
  0.12838916450705e2, -0.28642437219381e1, 0.56912683664855e0,
  -0.99962954584931e-1, -0.32632037778459e-2, 0.23320922576723e-3,
  -0.15334809857450e0, 0.29072288239902e-1, 0.37534702741167e-3,
  0.17296691702411e-2, -0.38556050844504e-3, -0.35017712292608e-4,
  -0.14566393631492e-4, 0.56420857267269e-5, 0.41286150074605e-7,
  -0.20684671118824e-7, 0.16409393674725e-8
};
static int  Bs2bI[44] =
 { -6, -6, -5, -5, -4, -4, -4, -3, -3, -3,
   -3, -2, -2, -2, -2, -1, -1, -1, -1, -1,
    0,  0,  0,  0,  0,  0,  0,  1,  1,  1,
    1,  1,  1,  2,  2,  2,  3,  3,  3,  4,
    4,  5,  5,  5
};
static int  Bs2bJ[44] =
{ 0, 11, 0, 11, 0, 1, 11, 0, 1, 11,
  12, 0, 1, 6, 10, 0, 1, 5, 8, 9,
   0, 1, 2, 4, 5, 6, 9, 0, 1, 2,
   3, 7, 8, 0, 1, 5, 0, 1, 3, 0,
   1, 0, 1, 2
};
  // 2c  p,s -> t 
static double Bs2cn[30] =
{ 0.90968501005365e3, 0.24045667088420e4, -0.59162326387130e3,
  0.54145404128074e3, -0.27098308411192e3, 0.97976525097926e3,
  -0.46966772959435e3, 0.14399274604723e2, -0.19104204230429e2,
  0.53299167111971e1, -0.21252975375934e2, -0.31147334413760e0,
  0.60334840894623e0, -0.42764839702509e-1, 0.58185597255259e-2,
  -0.14597008284753e-1, 0.56631175631027e-2, -0.76155864584577e-4,
  0.22440342919332e-3, -0.12561095013413e-4, 0.63323132660934e-6,
  -0.20541989675375e-5, 0.36405370390082e-7, -0.29759897789215e-8,
  0.10136618529763e-7, 0.59925719692351e-11, -0.20677870105164e-10,
  -0.20874278181886e-10, 0.10162166825089e-9, -0.16429828281347e-9
};
static int  Bs2cI[30] =
 { -2, -2, -1, 0, 0, 0, 0, 1, 1, 1,
   1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
   5,  5,  5,  6, 6,  7,  7,  7,  7,  7
};
static int  Bs2cJ[30] =
{ 0, 1, 0, 0, 1, 2, 3, 0, 1, 3,
  4, 0, 1, 2, 0, 1, 5, 0, 1, 4,
   0, 1, 2, 0, 1, 0, 1, 3, 4, 5
};


// function definition 
double bTT=1.0;/*K*/
double bPP=1.0; /*Mpa*/
double bHH=2000.0;/*Kj/Kg*/
double bSSa=2.0;/*Kj/KgK*/
double bSSb=0.7853;/*Kj/KgK*/
double bSSc=2.9251;/*Kj/KgK*/

// region 2  forward equations
 // basic equations
double pt2fv(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2_p= Gmo_p(rdp, rdt, GE2no, GE2Jo, 9);
	double Gmt2_p= Gmt_p(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);

	ret = rdp *(Gmo2_p+Gmt2_p);
	ret *= (R*t/p); 

	return (ret*VX2);
}
    
double pt2fu(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2_t= Gmo_t(rdp, rdt, GE2no, GE2Jo, 9);
	double Gmt2_t= Gmt_t(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);
	double Gmo2_p= Gmo_p(rdp, rdt, GE2no, GE2Jo, 9);
	double Gmt2_p= Gmt_p(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);

	ret = rdt *(Gmo2_t+Gmt2_t) -
		 rdp*(Gmo2_p+Gmt2_p);
    ret *= (R*t);
	
	return ret;
}
    
double pt2fs(double p, double t)
{
    double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2_t= Gmo_t(rdp, rdt, GE2no, GE2Jo, 9);
	double Gmt2_t= Gmt_t(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);
	double Gmo2= Gmo(rdp, rdt, GE2no, GE2Jo, 9);
	double Gmt2= Gmt(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);

	ret = (R)* (rdt *(Gmo2_t+Gmt2_t) -
		 (Gmo2+Gmt2));

	return ret;
}

double pt2fh(double p, double t)
{
    double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2_t= Gmo_t(rdp, rdt, GE2no, GE2Jo, 9);
	double Gmt2_t= Gmt_t(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);

	ret = rdt *(Gmo2_t+Gmt2_t);
	ret *= (R*t);

	return ret;
}

double pt2fcp(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2_tt= Gmo_tt(rdp, rdt, GE2no, GE2Jo, 9);
	double Gmt2_tt= Gmt_tt(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);


	ret = -rdt*rdt *(Gmo2_tt+Gmt2_tt);
	ret *= (R); 

	return ret;
}    

double pt2fcv(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2_tt= Gmo_tt(rdp, rdt, GE2no, GE2Jo, 9);
	double Gmt2_tt= Gmt_tt(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);
	double Gmt2_p= Gmt_p(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);
	double Gmt2_pp= Gmt_pp(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);
	double Gmt2_pt= Gmt_pt(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);


	ret = -rdt*rdt *(Gmo2_tt+Gmt2_tt);
	ret += -xen(1+rdp*(Gmt2_p-rdt*Gmt2_pt),2) /
		(1- rdp*rdp*Gmt2_pp);
	ret *= (R); 

	return ret;
}

double pt2fw(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmt2_p= Gmt_p(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);
	double Gmo2_tt= Gmo_tt(rdp, rdt, GE2no, GE2Jo, 9);
	double Gmt2_tt= Gmt_tt(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);
	double Gmt2_pp= Gmt_pp(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);
	double Gmt2_pt= Gmt_pt(rdp, rdt-0.5, GE2n, GE2I, GE2J, 43);


	ret = 1+rdp*Gmt2_p*(2+rdp*Gmt2_p);
	ret /= ((1- rdp*rdp*Gmt2_pp) +
		xen(1+rdp*(Gmt2_p-rdt*Gmt2_pt),2) /
		(rdt*rdt*(Gmo2_tt+Gmt2_tt)));
	ret *= (R*t); 

	return sqrt(ret*W2X2);
}

// Supplementary Equation for Metastable-Vapor Region
double pt2sv(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2s_p= Gmo_p(rdp, rdt, sGE2no, GE2Jo, 9);
	double Gmt2s_p= Gmt_p(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);

	ret = rdp *(Gmo2s_p+Gmt2s_p);
	ret *= (R*t/p); 

	return (ret*VX2);
	return (ret*VX2);
}
    
double pt2su(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2s_t= Gmo_t(rdp, rdt, sGE2no, GE2Jo, 9);
	double Gmt2s_t= Gmt_t(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);
	double Gmo2s_p= Gmo_p(rdp, rdt, sGE2no, GE2Jo, 9);
	double Gmt2s_p= Gmt_p(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);

	ret = rdt *(Gmo2s_t+Gmt2s_t) -
		 rdp*(Gmo2s_p+Gmt2s_p);
    ret *= (R*t);
	
	return ret;
}
    
double pt2ss(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2s_t= Gmo_t(rdp, rdt, sGE2no, GE2Jo, 9);
	double Gmt2s_t= Gmt_t(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);
	double Gmo2s= Gmo(rdp, rdt, sGE2no, GE2Jo, 9);
	double Gmt2s= Gmt(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);

	ret = (R)* (rdt *(Gmo2s_t+Gmt2s_t) -
		 (Gmo2s+Gmt2s));

	return ret;
}

double pt2sh(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2s_t= Gmo_t(rdp, rdt, sGE2no, GE2Jo, 9);
	double Gmt2s_t= Gmt_t(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);

	ret = rdt *(Gmo2s_t+Gmt2s_t);
	ret *= (R*t);

	return ret;
}

double pt2scp(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2s_tt= Gmo_tt(rdp, rdt, sGE2no, GE2Jo, 9);
	double Gmt2s_tt= Gmt_tt(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);


	ret = -rdt*rdt *(Gmo2s_tt+Gmt2s_tt);
	ret *= (R); 

	return ret;
} 

double pt2scv(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmo2s_tt= Gmo_tt(rdp, rdt, sGE2no, GE2Jo, 9);
	double Gmt2s_tt= Gmt_tt(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);
	double Gmt2s_p= Gmt_p(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);
	double Gmt2s_pp= Gmt_pp(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);
	double Gmt2s_pt= Gmt_pt(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);


	ret = -rdt*rdt *(Gmo2s_tt+Gmt2s_tt);
	ret += -xen(1+rdp*(Gmt2s_p-rdt*Gmt2s_pt),2) /
		(1- rdp*rdp*Gmt2s_pp);
	ret *= (R); 

	return ret;
}

double pt2sw(double p, double t)
{
	double ret;
	double rdp = p / PX2;
	double rdt = TX2 / t;
	double Gmt2s_p= Gmt_p(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);
	double Gmo2s_tt= Gmo_tt(rdp, rdt, sGE2no, GE2Jo, 9);
	double Gmt2s_tt= Gmt_tt(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);
	double Gmt2s_pp= Gmt_pp(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);
	double Gmt2s_pt= Gmt_pt(rdp, rdt-0.5, sGE2n, sGE2I, sGE2J, 13);


	ret = 1+rdp*Gmt2s_p*(2+rdp*Gmt2s_p);
	ret /= ((1- rdp*rdp*Gmt2s_pp) +
		xen(1+rdp*(Gmt2s_p-rdt*Gmt2s_pt),2) /
		(rdt*rdt*(Gmo2s_tt+Gmt2s_tt)));
	ret *= (R*t); 

	return sqrt(ret*W2X2);
}

// generic area 2 equations
double pt2v(double p, double t)
{
	double ret;

	switch (metastable(p, t))
	{
	case 1: ret = pt2sv(p,t);
		break;
	case 0: ret = pt2fv(p,t);
		break;
	default: ret = -1;
	}

	return ret;
}

double pt2u(double p, double t)
{
	double ret;

	switch (metastable(p, t))
	{
	case 1: ret = pt2su(p,t);
		break;
	case 0: ret = pt2fu(p,t);
		break;
	default: ret = -1;
	}

	return ret;
}

double pt2s(double p, double t)
{
	double ret;

	switch (metastable(p, t))
	{
	case 1: ret = pt2ss(p,t);
		break;
	case 0: ret = pt2fs(p,t);
		break;
	default: ret = -1;
	}

	return ret;
}

double pt2h(double p, double t)
{
	double ret;

	switch (metastable(p, t))
	{
	case 1: ret = pt2sh(p,t);
		break;
	case 0: ret = pt2fh(p,t);
		break;
	default: ret = -1;
	}

	return ret;
}

double pt2cp(double p, double t)
{
	double ret;

	switch (metastable(p, t))
	{
	case 1: ret = pt2scp(p,t);
		break;
	case 0: ret = pt2fcp(p,t);
		break;
	default: ret = -1;
	}

	return ret;
}

double pt2cv(double p, double t)
{
	double ret;

	switch (metastable(p, t))
	{
	case 1: ret = pt2scv(p,t);
		break;
	case 0: ret = pt2fcv(p,t);
		break;
	default: ret = -1;
	}

	return ret;
}

double pt2w(double p, double t)
{
	double ret;

	switch (metastable(p, t))
	{
	case 1: ret = pt2sw(p,t);
		break;
	case 0: ret = pt2fw(p,t);
		break;
	default: ret = -1;
	}

	return ret;
}


// Backward equations
 // 2b-2c subregion boundary
  // 2a  p,h -> t 
double t_2a_ph(double p,double h)
{
   double ret;
   double Pi= p / PX2;
   double Et= h / HX2;

   ret = polyXY2e(Pi, Et-2.1, B2an, B2aI, B2aJ, 34);

  return (1.0* ret);
}

  // 2b  p,h -> t 
double t_2b_ph(double p,double h)
{
   double ret;
   double Pi= p / PX2;
   double Et= h / HX2;

   ret = polyXY2e(Pi-2.0, Et-2.6, B2bn, B2bI, B2bJ, 38);

  return (1.0* ret);
}

  // 2c  p,h -> t 
double t_2c_ph(double p,double h)
{
   double ret;
   double Pi= p / PX2;
   double Et= h / HX2;

   ret = polyXY2e(Pi+25.0, Et-1.8, Bh2cn, Bh2cI, Bh2cJ, 23);

  return (1.0* ret);
}

  // 2a  p,s -> t 
double t_2a_ps(double p,double s)
{
   double ret;
   double Pi=pow(p / PX2, 0.25);
   double Sg= s / 2.0;

   ret = polyXY2e(Pi, Sg-2.0, Bs2an, Bs2aI, Bs2aJ, 46);

  return (1.0*ret);
}

  // 2b  p,s -> t 
double t_2b_ps(double p,double s)
{
   double ret;
   double Pi= p / PX2;
   double Sg= s / 0.7853;

   ret = polyXY2e(Pi, 10.0-Sg, Bs2bn, Bs2bI, Bs2bJ, 44);
   
  return (1.0*ret);
}

  // 2c  p,s -> t 
double t_2c_ps(double p,double s)
{
	double ret;
	double Pi= p / PX2;
	double Sg= s / 2.9251;

	ret = polyXY2e(Pi, 2.0-Sg, Bs2cn, Bs2cI, Bs2cJ, 30);
	
	return (1.0* ret);
}

//backward equation. generic 
  // B2bc-equation, is approximately s=5.81~5.85
double p_B2bc_h(double h)
{
	double ret;
	double rdh = h/1.0; /*Mpa*/

	ret = B2bcn[0]+B2bcn[1]*rdh+B2bcn[2]*rdh*rdh;

	return (PX2*ret);
}

double h_B2bc_p(double p)
{
	double ret;
	double rdp = p/PX2; /*Mpa*/

	ret = B2bcn[3]+sqrt((rdp-B2bcn[4])/B2bcn[2]);
   
	return (1.0*ret);
}

double ph2t(double p,double h)
{
	double ret;

	if (p<=4.0) ret=t_2a_ph(p,h);
	else
	{
		if (h>h_B2bc_p(p)) ret=t_2b_ph(p,h);
		else ret=t_2c_ph(p,h);
	}

	return (ret);
}

double ps2t(double p,double s)
{
   double ret;

   if (p<=4.0) ret=t_2a_ps(p,s);
   else
   {
	   if (s>5.85) ret=t_2b_ps(p,s);
	   else ret=t_2c_ps(p,s);
   }

   return (ret);
}

double ps2h(double p, double s)
{
	double t, h;

	t = ps2t(p, s);
	h = pt2h(p, t);

	return h;
}

double ts2p(double t, double s)
{
	double pu, pd;

	if (s>=9.1558)
	{
		pu=sbd4p(s);  /* border 4 TB3 */
		pd=PB0;
		return divid2r(ps2t, t, pu, pd, s, ERRORHI);
	}
	else
	{
		if (s>=6.0405)
		{
			pu=sbd4p(s);  /* border 4 TB3 */
			pd=t4p(s4ts(s));
			return divid2r(ps2t, t, pu, pd, s, ERRORHI);
		}
		else
		{
			pu=PB100;  /* PB100 */
			pd=t4p(s4ts(s));
			return divid2r(ps2t, t, pu, pd, s, ERRORHI);
		}
	}
}

double ts2h(double t, double s)
{
	double p, h;

	p = ts2p(t, s);
	h = pt2h(p, t);

	return h;
}

double th2p(double t, double h)
{
	double p, pu;

	if (t<TB2) pu=t4p(t);
	if (t<TB23) pu=p_B23_t(t);
	else pu=PB100;

	p = divid2r(ph2t, t, pu, PB0, h, ERRORHI);
	return p;
}

double th2s(double t, double h)
{
	double p,s ;

	p = th2p(t, h);
	s = pt2s(p, t);

	return s;
}

/* obsolete - 1999 
double hs2t(double h, double s)
{
	double t;

	t = divid2r(ts2h, h, TB3, TMIN, s, ERRORHI);

	return t;
} */

/* Sep.2001, Gaithersburg, Maryland, USA */
double hs2t(double h, double s) /* implemented Aug23, 2005 */
{
	double p;

	p = hs2p(h,s);

	return ph2t(p,h);
}

/* 97 implementation. obsolete 
double hs2p(double h, double s)
{	double pu, pd;

	if (s>=9.1558)
	{
		pu=sbd4p(s);  // border 4 TB3
		pd=PB0;
		return divid2r(ps2h, h, pu, pd, s, ERRORHI);
	}
	else
	{
		if (s>=6.0405)
		{
			pu=sbd4p(s);  // border 4 TB3 
			pd=t4p(s4ts(s));
			return divid2r(ps2h, h, pu, pd, s, ERRORHI);
		}
		else
		{
			pu=PB100;  // PB100 
			pd=t4p(s4ts(s));
			return divid2r(ps2h, h, pu, pd, s, ERRORHI);
		}
	}
}
obsolete NIF97 implementation */

/* New 01 formula hs2p() */
double h_B2ab_s(double s)
{
	double hx = 1;
	double sx = 1;
	double n_B2ab_hs[4]={-0.349898083432139e4, 0.257560716905876e4,
		-0.421073558227969e3, 0.276349063799944e2};

	double rds, rdh;

	rds = s/sx;
	rdh	= polynom (rds, 3, n_B2ab_hs);

	return rdh * hx;
}

double p_2a_hs(double h, double s)
{
	double px = 4.0;
	double hx = 4200.0;
	double sx = 12.0;
	int I[29] = {0,0,0,0,0, 0,1,1,1,1, 1,1,1,1,1, 1,2,2,2,3,
		3,3,3,3,4, 5,5,6,7};
	int J[29] = {1,3,6,16,20, 22,0,1,2,3, 5,6,10,16,20,
		22,3,16,20,0, 2,3,6,16,16, 3,16,3,1};
	double n[29] = {-1.82575361923032E-02, -1.25229548799536E-01,
		5.92290437320145E-01, 6.04769706185122E+00, 2.38624965444474E+02,
		-2.98639090222922E+02, 5.12250813040750E-02, -4.37266515606486E-01,
		4.13336902999504E-01, -5.16468254574773E+00, -5.57014838445711E+00,
		1.28555037824478E+01, 1.14144108953290E+01, -1.19504225652714E+02,
		-2.84777985961560E+03, 4.31757846408006E+03, 1.12894040802650E+00,
		1.97409186206319E+03, 1.51612444706087E+03, 1.41324451421235E-02,
		5.85501282219601E-01, -2.97258075863012E-00, 5.94567314847319E+00,
		-6.23656565798905E+03, 9.65986235133332E+03, 6.81500934948134E+00,
		-6.33207286824489E+03, -5.58919224465760E+00, 4.00645798472063E-02};

	double rdp, rdh, rds;

	rdh = h/hx - 0.5;
	rds = s/sx - 1.2;
	rdp = polyXY2e(rdh, rds, n, I, J, 29);

	return (xen(rdp, 4) * px);
}

double p_2b_hs(double h, double s)
{
	double px = 100.0;
	double hx = 4100.0;
	double sx = 7.9;
	int I[33] = {0,0,0,0,0, 1,1,1,1,1, 1,2,2,2,3, 3,3,3,4,4,
		5,5,6,6,6, 7,7,8,8,8, 8,12,14};
	int J[33] = {0,1,2,4,8, 0,1,2,3,5, 12,1,6,18,0, 1,7,12,1,16,
		1,12,1,8,18, 1,16,1,3,14, 18,10,16};
	double n[33] = {0.801496989929495e-1, -0.543862807146111,
		0.337455597421283, 0.890555451157450e1, 0.313840736431485e3,
		0.797367065977789, -0.121616973556240e1, 0.872803386937477e1,
		-0.169769781757602e2, -0.186552827328416e3, 0.951159274344237e5,
		-0.189168510120494e2, -0.433407037194840e4, 0.543212633012715e9,
		0.144793408386013, 0.128024559637516e3, -0.672309534071268e5,
		0.336972380095287e8, -0.586634196762720e3, -0.221403224769889e11,
		0.171606668708389e4, -0.570817595806302e9, -0.312109693178482e4,
		-0.207841384633010e7, 0.305605946157786e13, 0.322157004314333e4,
		0.326810259797295e12, -0.144104158934487e4, 0.410694867802691e3,
		0.109077066873024e12, -0.247964654258893e14, 0.188801906865134e10,
		-0.123651009018773e15};

	double rdp, rdh, rds;

	rdh = h/hx - 0.6;
	rds = s/sx - 1.01;
	rdp = polyXY2e(rdh, rds, n, I, J, 33);

	return (xen(rdp, 4) * px);
}

double p_2c_hs(double h, double s)
{
	double px = 100.0;
	double hx = 3500.0;
	double sx = 5.9;
	int I[31] = {0,0,0,0,0, 0,1,1,1,1, 1,2,2,2,2, 2,3,3,3,3,
		3,4,5,5,5, 5,6,6,10,12, 16};
	int J[31] = {0,1,2,3,4, 8,0,2,5,8, 14,2,3,7,10,
		18,0,5,8,16, 18,18,1,4,6, 14,8,18,7,7, 10};
	double n[31] = {0.112225607199012, -0.339005953606712e1,
		-0.320503911730094e2, -0.197597305104900e3, -0.407693861553446e3,
		0.132943775222331e5, 0.170846839774007e1, 0.373694198142245e2,
		0.358144365815434e4, 0.423014446424664e6, -0.751071025760063e9,
		0.523446127607898e2, -0.228351290812417e3, -0.960652417056937e6,
		-0.807059292526074e8, 0.162698017225669e13, 0.772465073604171,
		0.463929973837746e5, -0.137317885134128e8, 0.170470392630512e13,
		-0.251104628187308e14, 0.317748830835520e14, 0.538685623675312e2,
		-0.553089094625169e5, -0.102861522421405e7, 0.204249418756234e13,
		0.273918446626977e9, -0.263963146312685e16, -0.107890854108088e10,
		-0.296492620980124e11, -0.111754907323424e16};

	double rdp, rdh, rds;

	rdh = h/hx - 0.7;
	rds = s/sx - 1.1;
	rdp = polyXY2e(rdh, rds, n, I, J, 31);

	return (xen(rdp, 4) * px);
}

double hs2p(double h, double s) /* Aug-17-2005 */
{

	if (s<5.85) return p_2c_hs(h, s);
	else
	{
		if (h>h_B2ab_s(s)) return p_2b_hs(h, s);
		else return p_2a_hs(h, s);
	}
}


double pv2t(double p, double v)
{
	double t;

	t = divid2(pt2v, v, p, TB3, TMIN, ERRORHI);

	return t;
}

double tv2p(double t, double v)
{
	double p;

	p = divid2r(pt2v, v, PB100, PB0, t, ERRORHI);

	return p;
}

double ph2v(double p, double h)
{
	double t;

	t=ph2t(p,h);
	return (pt2v(p,t));
}

double vh2p(double v, double h)
{
	double p;

	p=divid2r(ph2v, v, PB0, PB100, h, ERRORHI2);

	return p;
}

double ps2v(double p, double s)
{
	double t;

	t=ps2t(p,s);
	return (pt2v(p,t));
}

double ts2v(double t, double s)
{
	double p;

	p=ts2p(t,s);

	return (pt2v(p,t));
}

double vs2p(double v, double s)
{
	double p, pu, pd;

	if (s>=9.1558)
	{
		pu=sbd4p(s);  /* border 4 TB3 */
		pd=PB0;
		return divid2r(ps2v, v, pu, pd, s, ERRORHI);
	}
	else
	{
		if (s>=6.0405)
		{
			pu=sbd4p(s);  /* border 4 TB3 */
			pd=t4p(s4ts(s));
			return divid2r(ps2v, v, pu, pd, s, ERRORHI);
		}
		else
		{
			pu=PB100;  /* PB100 */
			pd=t4p(s4ts(s));
			p=divid2r(ps2v, v, pu, pd, s, ERRORHI);
			return p;
		}
	}

}


 // area 2 metastable-vapor region 
int metastable(double p, double t)
{
	return 0;
}

/*int metastable(double p, double t)
{
	if (p>10.0) return 0;
	double ts = p4t(p);
	if (t > ts) return 0;

	double hv = pt2fh(p, ts);
	double hw = pt1h(p, ts);
	double hx = pt2fh(p, t);
	double sgm = (hx-hw)/(hv-hw);

	if ( sgm > 1.0) return 0;
	if ( sgm < 0.95) return -1;
	else return 1;
}*/