/* 
   Filename: supple.cpp
   Supplementary Equations for other Porperties
*/
#include "stdafx.h"

#include <math.h>
#include "supple.h"
#include "tools.h"
#include "region.h"
#include "generic.h"
#include "area3.h"
#include "areas.h"


/* Dynamic Viscosity -Revised Release 1985 issued in 2003- obsolete
DllExport double WINAPI dyn_visc(double p, double t, double r, int area)
{
	double Csi0, Csi1, ret;
	double Ytx, Rox, Tx, Tao, Dlt; 
	int areaflag=0;
	double n0[4]={0.1e1, 0.978197, 0.579829, -0.202354};
	int I[19]={0,0,0,0,1,1,1,1,2,2,2,3,3,3,3,4,4,5,6};
	int J[19]={0,1,4,5,0,1,2,3,0,1,2,0,1,2,3,0,3,1,3};
	double n[19]={0.51320047,0.3205656,-0.7782567,0.1885447,0.2151778,
		0.7317883,0.1241033e1,0.1476783e1,-0.2818107,-0.1070786e1,
		-0.1263184e1,0.1778064,0.4605040,0.2340379,-0.4924179,
		-0.4176610e-1,0.1600435,-0.1578386e-1,-0.3629481e-2};

	if (t>= 273.15 && t<=423.15 && p<= 500) areaflag = 1;
	if (t>= 423.15 && t<=873.15 && p<= 350) areaflag = 1;
	if (t>= 873.15 && t<=1173.15 && p<= 300) areaflag = 1;

	if (areaflag == 0) return AREAX;

	if (r<0)
	{
		r = 1/pt_v(p,t,area);
		if (r<=0) return 1/r;
	}

	Rox = 317.763; // kg m^-3
	Tx = 647.226; // K 
	Tao = Tx / t;
	Csi0  = 1.0 / (sqrt(Tao) * polynom(Tao, 3, n0));

	Dlt = r / Rox;
	Csi1 = exp(Dlt * polyXY2e((Dlt-1.0), (Tao-1.0), n, I, J, 19));

	Ytx = 55.071e-6; // Pa s
	ret = Csi0 * Csi1 * Ytx;

	return ret;
}*/

// Dynamic Viscosity - Berlin, Germany September 2008
DllExport double WINAPI dyn_visc(double p, double t, double r, int area)
{
	double tm, T_, Tx = 647.096; //K
	double Px = 22.064; //MPa
	double Ro_, Rox = 322.0; //kg/m^-3
	double Mu_, Mux=1e-6; //Pa.s
	double Mu0, Mu1; //Mu2 is omitted for industrial application
	int i, areaflag=0;

	if (r<0)
	{
		r = 1/pt_v(p,t,area);
		if (r<=0) return 1/r;
	}
	Ro_=r/Rox;
	T_ = t / Tx;

	tm=tm_p(p);
	if (p>0 && p< Pt && t>= 273.16 && t<=1173.15 ) areaflag = 1;
	if (p>=Pt && p<= 300 && t>= tm && t<=1173.15 ) areaflag = 1;
	if (p>300 && p<= 350 && t>= tm && t<=873.15 ) areaflag = 1;
	if (p>350 && p<= 500 && t>= tm && t<=433.15 ) areaflag = 1;
	if (p>500 && p<= 1000 && t>= tm && t<=373.15 ) areaflag = 1;
	if (areaflag == 0) return AREAX;

	//Mu0
	double H[4]={1.67752, 2.20462, 0.6366564, -0.241605};
	Mu0 = 100*sqrt(T_) / polynom(1/T_, 3, H);

	//Mu1
	double A[6][7]={{5.20094e-1, 2.22531e-1, -2.81378e-1, 1.61913e-1, -3.25372e-2, 0, 0},
	{8.50895e-2, 9.99115e-1, -9.06851e-1, 2.57399e-1, 0, 0, 0},
	{-1.08374, 1.88797, -7.72479e-1, 0, 0, 0, 0},
	{-2.89555e-1, 1.26613, -4.89837e-1, 0, 6.98452e-2, 0, -4.35673e-3},
	{0, 0, -2.57040e-1, 0, 0, 8.72102e-3, 0},
	{0, 1.20573e-1, 0, 0, 0, 0, -5.93264e-4}};
	double v[6];

	for (i=0;i<=5;i++)
		v[i]=polynom((Ro_-1),6,A[i]);
	Mu1=exp(Ro_ * polynom((1/T_-1),5,v));

	//Mu2 is omitted for industrial application
	Mu_ =  Mu0 * Mu1;
	return Mu_ * Mux;
}

/* Thermal Conductivity -obsolete version
DllExport double WINAPI thm_cond(double p, double t, double r, int area)
{
	double Rmd0, Rmd1, Rmd2, part1, part2, k, ret;
	double Rox, Tx, Tao, Dlt, Px, Ytx = 55.071e-6;; 
	int areaflag=0;
	double n0[4]={0.1e1, 0.6978267e1, 0.2599096e1, -0.998254};
	double n[5][6]={
		{0.13293046e1,-0.40452437,0.24409490,0.18660751e-1,-0.12961068,0.44809953e-1},
		{0.17018363e1,-0.22156845e1,0.16511057e1,-0.76736002,0.37283344,-0.11203160},
		{0.52246158e1,-0.10124111e2,0.49874687e1,-0.27297694,-0.43083393,0.13333849},
		{0.87127675e1,-0.95000611e1,0.43786606e1,-0.91783782,0,0},
		{-0.18525999e1,0.93404690,0,0,0,0}};

	if (t>= 273.15 && t<=398.15 && p<= 400) areaflag = 1;
	if (t>= 398.15 && t<=523.15 && p<= 200) areaflag = 1;
	if (t>= 523.15 && t<=673.15 && p<= 150) areaflag = 1;
	if (t>= 673.15 && t<=1073.15 && p<= 100) areaflag = 1;

	if (areaflag == 0) return AREAX;

	if (r<0)
	{
		r = 1/pt_v(p,t,area);
		if (r<=0) return 1/r;
	}

	Px = 22.115;   // MPa 
	Tx = 647.226; // K 
	Tao = Tx / t;
	Rmd0  = 1.0 / (sqrt(Tao) * polynom(Tao, 3, n0));

	Rox = 317.763; // kg m^-3 
	Dlt = r / Rox;
	Rmd1 = exp(Dlt * polynom2e((Tao-1.0), (Dlt-1.0), &n[0][0], 5, 6));

	area = region_pt(p, t);

	if (area!=AREA3)
	{
		ret = Rmd0*Rmd1;
		return (ret*0.4945);
	}
	
	// get partial value
//	Rox = ROC;     // this coefficient refer to area 3 
//	Tx = TC;
	Px = 22.115e3;

//	Dlt = r / Rox;
//	Tao = Tx / t;
	k = Rox*R*Tx/Px;
	
	part1 = k*Dlt*Dlt*(Fi_d(Dlt, Tao) - Tao*Fi_dt(Dlt, Tao));
	
	part2 = Tao/(k*(2*Fi_d(Dlt, Tao) + Dlt*Fi_dd(Dlt, Tao)));

	// other 
	Rmd2 = 0.0013848/(dyn_visc(p, t, r, area)/Ytx);

//	Rox = 317.763;
//	Tx = 647.226;

//	Dlt = r / Rox;
//	Tao = Tx / t;

	Rmd2 /= (Tao*Tao*Dlt*Dlt);
	Rmd2 *= (part1*part1);
	Rmd2 *= pow(part2, 0.4678);
	Rmd2 *= sqrt(Dlt);
	Rmd2 *= exp((-18.66*(1/Tao-1)*(1/Tao-1)) - xen((Dlt-1), 4));

	ret = Rmd0*Rmd1 + Rmd2;

	return (ret*0.4945);
}

*/

/* Pressure along the Melting and Sublimation Curves of Ordinary Water Substance (September Czech Republic 2011) */
/* required by Thermal Conductivity */
double pm_Ih_t(double t) //273.16k - 251.165k
{
	double Tx=Tt, Px=Pt, p_=1, t_; //273.16K, 611.657e-6MPa
	double a[3]={0.119539337e7, 0.808183159e5, 0.333826860e4};
	double b[3]={0.300000e1, 0.257500e2, 0.103750e3};
	int i;
	
	t_ = t/Tx;
	for (i=0;i<3;i++)
		p_ += a[i]*(1-pow(t_,b[i]));

	return p_*Px;
}

double pm_III_t(double t) //273.16k - 251.165k
{
	double Tx=251.165, Px=208.566, p_, t_; //K, MPa

	t_ = t/Tx;
	p_= 1-0.299948*(1-xen(t_,60));

	return p_*Px;
}

double pm_V_t(double t)
{
	double Tx=256.164, Px=350.1, p_, t_; //K, MPa

	t_ = t/Tx;
	p_= 1-1.18721*(1-xen(t_,8));

	return p_*Px;
}
double pm_VI_t(double t)
{
	double Tx=273.31, Px=632.4, p_, t_; //K, MPa

	t_ = t/Tx;
	p_= 1-1.07476*(1-pow(t_,4.6));

	return p_*Px;
}

double pm_VII_t(double t)
{
	double Tx=355, Px=2216, p_, t_; //K, MPa

	t_ = t/Tx;
	p_= 1.73683*(1-1/t_)-0.0544606*(1-xen(t_,5))+0.806106e-7*(1-xen(t_,22));

	return exp(p_)*Px;
}

double pm_t(double t)
{
	return 0;
}
double tm_p(double p)
{
	double Tx_Ih=Tt, Px_Ih=Pt, p_=1, t_; //273.16K, 611.657e-6MPa
	double Tx_III=251.165, Px_III=208.566; //K, MPa
	double Tx_V=256.164, Px_V=350.1; //K, MPa
	double Tx_VI=273.31, Px_VI=632.4; //K, MPa
	double Tx_VII=355, Px_VII=2216; //K, MPa

	if (p<Px_Ih) return AREAX;
	if (p<Px_III)  //in Ih
	{
		p_=p/Px_Ih;
		t_=aitkin1(pm_III_t,p_,0.95, ERRORHI);
		//t_=divid1(pm_III_t,p_,Tx_III/Tx_Ih,1, ERRORHI);
		return t_*Tx_Ih;
	}
	if (p<Px_V) //in III
	{
		p_=p/Px_III;
		t_=pow((1-(1-p_)/0.299948),1/60);
		return t_*Tx_III;
	}
	if (p<Px_VI) //in V
	{
		p_=p/Px_V;
		t_=pow((1-(1-p_)/1.18721),1/8);
		return t_*Tx_V;
	}
	if (p<Px_VII)
	{
		p_=p/Px_VI;
		t_=pow((1-(1-p_)/1.07476),1/4.6);
		return t_*Tx_VI;
	}
	if(p<20617)	//VII: 355k-715k; 2216MP-20617MP
	{
		p_=p/Px_VII;
		t_=aitkin1(pm_VII_t,p_,1.5, ERRORHI);
		return t_*Tx_VII; 
	}

	return AREAX;
}

/* Thermal Conductivity - Czech Republic September 2011 */
DllExport double WINAPI thm_cond(double p, double t, double r, int area)
{
	double Rmd0, Rmd1, Rmd2, cp, cv, k, y,  Z, dp;
	double tm, T_, Tx = 647.064; //K
	double Px = 22.064; //MPa
	double Ro_, Rox = 322.0; //kg/m^-3
	double Rmd_, Rmdx = 1e-3; //W.K^-1.m^-1
	double Mu_, Mux=1e-6; //Pa.s
	double Cp_, Rx = 0.46151805; //kJ.kg^-1.K^-1, Ytx = 55.071e-6;; 
	int i,j, areaflag=0;
	double n0[5]={2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4};
	double L[5][6]={
		{1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634,0.00609859258},
		{2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019, -0.00719201245},
		{2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.0205938816},
		{-1.21051378, 1.60812989, -0.621178141, 0.0716373224, 0, 0},
		{-2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842}};
	double v[5];

	if (r<0)
	{
		r = 1/pt_v(p,t,area);
		if (r<=0) return 1/r;
	}
	Ro_=r/Rox;

	tm=tm_p(p);
	if (p>0 && p< Pt && t>= 273.16 && t<=1173.15 ) areaflag = 1;
	if (p>=Pt && p<= 100 && t>= tm && t<=1173.15 ) areaflag = 1;
	if (p>100 && p<= 250 && t>= tm && t<=874 ) areaflag = 1;
	if (p>250 && p<= 687 && t>= tm && t<=573 ) areaflag = 1;
	if (p>687 && p<= 785 && t>= tm && t<=403 ) areaflag = 1;
	if (p>785 && p<= 1000 && t>= tm && t<=348 ) areaflag = 1;

	if (areaflag == 0) return AREAX;

	T_ = t/Tx;
	Rmd0 = sqrt(T_)/polynom(1/T_,4,n0);

	for (i=0;i<=4;i++)
		v[i]=polynom((Ro_-1),5,L[i]);
	Rmd1=exp(Ro_ * polynom((1/T_-1),4,v));

	//Rmd2 Critical enhancement - only implementation for industrial use
	//A[j][i] - transposed for easier implementation 
	double A[5][6]={{6.53786807199516, -5.61149954923348, 3.39624167361325, -2.27492629730878, 10.2631854662709, 1.97815050331519},
	{6.52717759281799, -6.30816983387575, 8.08379285492595, -9.82240510197603, 12.1358413791395, -5.54349664571295},
	{5.35500529896124, -3.96415689925446, 8.91990208918795, -12.033872950579, 9.19494865194302, -2.16866274479712},
	{1.55225959906681, 0.464621290821181, 8.93237374861479, -11.0321960061126, 6.1678099993336, -0.965458722086812},
	{1.11999926419994, 0.595748562571649, 9.8895256507892, -10.325505114704, 4.66861294457414, -0.503243546373828}};

	double Zeta0, Zeta, Chi_, Ksi;
	//range of validity of industrial equation
	areaflag = 0;
	if (p<= 100 && t>= 273.15 && t<=1073.15 ) areaflag = 1;
	if (p<= 50 && t>= 1073.15 && t<=1173.15 ) areaflag = 1;
	if (areaflag == 0) return AREAX;	

	if (Ro_ <=0.310559006) j = 0;
	if (Ro_ >0.310559006 && Ro_ <=0.776397516) j = 1; 
	if (Ro_ >0.776397516 && Ro_ <=1.242236025) j = 2; 
	if (Ro_ >1.242236025 && Ro_ <=1.863354037) j = 3; 
	if (Ro_ >1.863354037) j = 4; 

	Zeta0 = 1/polynom(Ro_,5,A[j]);
	dp = 0.0001*Px;
	Zeta = ((1/pt_v(p+dp,t,area)-1/pt_v(p-dp,t,area))/(2*dp))*(Px/Rox);
	Chi_ = Ro_*(Zeta - Zeta0 *1.5/T_); //TR_ = 1.5

	Ksi= 0.13 * pow(Chi_/0.06, 0.63/1.239); //Ksi0=0.13nm, Gama0 = 0.06, nu=0.630, gama=1.239

	
	y=2.5*Ksi; //qD_^-1 = 0.4nm
	cp = pt_cp(p, t, area);
	cv = pt_cv(p, t, area);
	k= cv/ cp; //actually k^-1, k = cp/cv
	if (y<1.2e-7) Z = 0;
	else Z=(2/(3.1415926*y))*(((1-k)*atan(y)+k*y) - (1-exp(-1/(1/y+y*y/(3*Ro_*Ro_)))));

	Mu_ = dyn_visc(p, t, r, area)/Mux;

	Cp_ = cp/Rx;
	Rmd2 = 177.8514 * Ro_ * Cp_ * T_ * Z / Mu_; //Lambda=177.8514

	Rmd_ = Rmd0*Rmd1 + Rmd2;

	return (Rmd_ * Rmdx);
}

/* Surface Tension */
DllExport double WINAPI surf_tens(double t, int area)
{
	double Sgm, Sta, ret;
	int areaflag = 0;
	
	if (t>=Tt && t<=TC) areaflag = 1;
	if (areaflag == 0) return AREAX;

	Sta = 1.0 - t/TC;
	Sgm = 235.8*pow(Sta,1.256) * (1.0 - 0.625*Sta);
	ret = Sgm * 1;

	return ret;
}

/* Static Dielectric Constant */
DllExport double WINAPI sta_dielec(double p, double t, double r, int area)
{
	double k = 1.380658e-23;  /* J K^-1 */
	double na = 6.0221367e23; /* mol^-1 */
	double Alf = 1.636e-40; /* C^2 J^-1 m^2 */
	double Eps0 = 8.854187817e-12; /* C^2 J^-1 m^-1 */
	double Mu = 6.138e-30; /* C m */
	int I[11] = {1,1,1,2,3,3,4,5,6,7,10};
	int J[11] = {1,4,10,6,6,10,8,8,20,2,40}; /* x^0.25 */
	double n[12] = {0.978224486826,-0.957771379375,0.237511794148,0.714692244396,
		-0.298217036956,-0.108863472196,0.949327488264e-1,-0.980469816509e-2,
		0.165167634970e-4,0.937359795772e-4,-0.123179218720e-9,0.196096504426e-2};
	double  Dlt, Tao, g, A, B, ret;
	int areaflag;

	if (t>= 273.15 && t<=323.15 && p<= 600) areaflag = 1;

	if (areaflag == 0) return AREAX;

	if (r<0)
	{
		r = 1/pt_v(p,t,area);
		if (r<=0) return 1/r;
	}

	Tao = pow((TC/t), 0.25);
	Dlt = r / ROC;
	g = 1+ polyXY2e(Dlt, Tao, n, I, J, 11)+ n[11]*Dlt*pow((t/228-1),-1.2);
	
	A = 1e3*na*Mu*Mu*r*g/(M*Eps0*k*t);
	B = 1e3*na*Alf*r/(3*M*Eps0);

	ret = (1+A+5*B+pow((9+2*A+18*B+A*A+10*A*B+9*B*B), 0.5))/(4-4*B);

	return ret;
}

/* Refractive Index */
DllExport double WINAPI Ref_Index(double p, double t, double r, double rmd, int area)
{
	double Sta, ret, Rmd, Dlt, A;
	double RmdUV = 0.229202; 
	double RmdIR = 5.432937;
	double a[8] = {0.244257733, 0.974634476e-2,-0.373234996e-2,0.268678472e-3,
		0.158920570e-2,0.245934259e-2,0.900704920,-0.166626219e-1};
	int areaflag = 0;
	
	if (t>=261.15 && t<=773.15) areaflag = 1;
	if (rmd>=0.2 && rmd<=1.1) areaflag = 1;
	if (areaflag == 0) return AREAX;

	if (r<0)
	{
		r = 1/pt_v(p,t,area);
		if (r<=0) return 1/r;
	}

	if (r<=0 || r>1060) return AREAX;
	Dlt = r/1000;
	Rmd = rmd / 0.589;
	Sta = t/T0;

	Rmd *= Rmd;
	RmdUV *= RmdUV;
	RmdIR *= RmdIR;
	A = Dlt*(a[0]+a[1]*Dlt+a[2]*Sta+a[3]*Rmd*Sta+a[4]/Rmd+a[5]/(Rmd-RmdUV)
		+a[6]/(Rmd-RmdIR)+a[7]*Dlt*Dlt);
	ret = pow((2*A+1)/(1-A), 0.5);
	
	return ret;
}
