/* Filename: region.c */
#include "stdafx.h"

#include <math.h>
#include "region.h"
#include "tools.h"
#include "area5.h"
#include "area4.h"
#include "area3.h"
#include "area2.h"
#include "area1.h"
#include "areas.h"

DllExport int WINAPI  region_pt(double p, double t)
{
	if (t<TMIN) return (T_OUTRANGE);

	if (t<=TB2)
	{
		if (p<=PB0) return P_OUTRANGE;
		if (p<=t4p(t)) return AREA2;
		if (p<=PB100) return AREA1;
		return P_OUTRANGE;
	}

	if (t<TB23)
	{
		if (p<=PB0) return P_OUTRANGE;
		if (p<=p_B23_t(t)) return AREA2;
		if (p<=PB100) return AREA3;
		return P_OUTRANGE;
	}

	if (t<=TB3)
	{
		if (p<=PB0) return P_OUTRANGE;
		if (p<=PB100) return AREA2;
		return P_OUTRANGE;
	}

	if (t<=TB4)
	{
		if (p<=PB0) return P_OUTRANGE;
		if (p<=PB50) return AREA5;
		return AREAX;
	}

	return T_OUTRANGE;
}

DllExport int WINAPI  region_ph(double p, double h)
{
	if (p<=PB0) return P_OUTRANGE;

	if (p<=t4p(TB2))
	{
		if (h<pt1h(p, TMIN)) return H_OUTRANGE;
		if (h<=p4hs1(p)) return AREA1;
		if (h<p4hs2(p)) return AREA4;
		if (h<=pt2h(p, TB3)) return AREA2;
		return AREAX;
	}
	
	if (p<PC)
	{
		if (h<pt1h(p, TMIN)) return H_OUTRANGE;
	    if (h<=pt1h(p, TB2)) return AREA1;
		if (h<=p4hs1(p)) return AREA3;
		if (h<p4hs2(p)) return AREA4;
		if (h<pt2h(p, t_B23_p(p))) return AREA3;
		if (h<=pt2h(p, TB3)) return AREA2;
		return AREAX;
	}

	if (p<=PB50)
	{
		if (h<pt1h(p, TMIN)) return H_OUTRANGE;
	    if (h<=pt1h(p, TB2)) return AREA1;
		if (h<pt2h(p, t_B23_p(p))) return AREA3;
		if (h<=pt2h(p, TB3)) return AREA2;
		if (h<=pt5h(p, TB4)) return AREA5;
		return H_OUTRANGE;
	}

	if (p<=PB100)
	{
		if (h<pt1h(p, TMIN)) return H_OUTRANGE;
	    if (h<=pt1h(p, TB2)) return AREA1;
		if (h<pt2h(p, t_B23_p(p))) return AREA3;
		if (h<=pt2h(p, TB3)) return AREA2;
		return AREAX;
	}

	return P_OUTRANGE;
}

DllExport int WINAPI  region_ps(double p, double s)
{
	if (p<=PB0) return P_OUTRANGE;

	if (p<=t4p(TB2))
	{
		if (s<pt1s(p, TMIN)) return S_OUTRANGE;
		if (s<=p4ss1(p)) return AREA1;
		if (s<p4ss2(p)) return AREA4;
		if (s<=pt2s(p, TB3)) return AREA2;
		return AREAX;
	}
	
	if (p<PC)
	{
		if (s<pt1s(p, TMIN)) return S_OUTRANGE;
	    if (s<=pt1s(p, TB2)) return AREA1;
		if (s<=p4ss1(p)) return AREA3;
		if (s<p4ss2(p)) return AREA4;
		if (s<pt2s(p, t_B23_p(p))) return AREA3;
		if (s<=pt2s(p, TB3)) return AREA2;
		return AREAX;
	}

	if (p<=PB50)
	{
		if (s<pt1s(p, TMIN)) return S_OUTRANGE;
	    if (s<=pt1s(p, TB2)) return AREA1;
		if (s<pt2s(p, t_B23_p(p))) return AREA3;
		if (s<=pt2s(p, TB3)) return AREA2;
		if (s<=pt5s(p, TB4)) return AREA5;
		return S_OUTRANGE;
	}

	if (p<=PB100)
	{
		if (s<pt1s(p, TMIN)) return S_OUTRANGE;
	    if (s<=pt1s(p, TB2)) return AREA1;
		if (s<pt2s(p, t_B23_p(p))) return AREA3;
		if (s<=pt2s(p, TB3)) return AREA2;
		return AREAX;
	}

	return P_OUTRANGE;
}


DllExport int WINAPI  region_th(double t, double h)
{
	if (t<TMIN) return T_OUTRANGE;

	if (t<=TB2)
	{
		if (h<0 || h>=pt2h(PB0, t)) return H_OUTRANGE;
		if (h<=pt1h(PB100,t)&&h>=tx4h(t,0)) return AREA1;
		if (h>=t4hs2(t)) return AREA2;
		return AREA4;
	}

	if (t<TC)
	{
		if (h<pt3h(PB100, t) || h>=pt2h(PB0, t)) return H_OUTRANGE;
		if (h>=pt2h(p_B23_t(t),t)) return AREA2;
		if (h<=t4hs1(t)||h>=t4hs2(t)) return AREA3;
		return AREA4;
	}

	if (t<TB23)
	{
		if (h<pt3h(PB100, t) || h>=pt2h(PB0, t)) return H_OUTRANGE;
		if (h>=pt2h(p_B23_t(t),t)) return AREA2;
		return AREA3;
	}
	
	if (t<=TB3)
	{
		if (h<pt2h(PB100, t) || h>=pt2h(PB0, t)) return H_OUTRANGE;
		return AREA2;
	}

	if (t<=TB4)
	{
		if (h<pt5h(PB50, t)) return AREAX;
		if (h>=pt5h(PB0, t)) return H_OUTRANGE;
		return AREA5;
	}

	return T_OUTRANGE;
}

DllExport int WINAPI  region_ts(double t, double s)
{
	if (t<TMIN) return T_OUTRANGE;

	if (t<=TB2)
	{
		if (s<pt1s(PB100, t) || s>=pt2s(PB0, t)) return S_OUTRANGE;
		if (s<=t4ss1(t)) return AREA1;
		if (s>=t4ss2(t)) return AREA2;
		return AREA4;
	}

	if (t<TC)
	{
		if (s<pt3s(PB100, t) || s>=pt2s(PB0, t)) return S_OUTRANGE;
		if (s>=pt2s(p_B23_t(t),t)) return AREA2;
		if (s<=t4ss1(t)||s>=t4ss2(t)) return AREA3;
		return AREA4;
	}

	if (t<TB23)
	{
		if (s<pt3s(PB100, t) || s>=pt2s(PB0, t)) return S_OUTRANGE;
		if (s>=pt2s(p_B23_t(t),t)) return AREA2;
		return AREA3;
	}
	
	if (t<=TB3)
	{
		if (s<pt2s(PB100, t) || s>=pt2s(PB0, t)) return S_OUTRANGE;
		return AREA2;
	}

	if (t<=TB4)
	{
		if (s<pt5s(PB50, t)) return AREAX;
		if (s>=pt5s(PB0, t)) return S_OUTRANGE;
		return AREA5;
	}

	return T_OUTRANGE;
}

/*
DllExport int WINAPI  region_hs(double h, double s)
{
	double p100[5]={97.65457,274.55289,33.34537,3.46198,0.02065};
	double p0[5]={-76359.13632, 27976.05061, -3677.87679,
		208.88107, -4.16865};
	double t800u[5]={83970.73764, -51203.2078, 11916.52765,
		-1205.68817, 45.04378};
	double p10[5]={3596.95945, -844.36753, 185.65815, -18.41495, 1.35144};
	double t2000[2][3]={
		{9.42372, 9.66083, 13.9046},
		{7374.31, 7375.63, 7376.98}};

	double b13[4]={25328.61982, -18178.1734, 4511.41216, -359.27007};
/*	double b23[2][6]={
		{2563.59, 2609.59, 2612.18, 2655.54, 2725.59, 2812.94},	
		{5.21089, 5.1344, 5.06319, 5.04847, 5.06512, 5.09797}};*
	
	
	if (s<SMIN || s>13.9406) return S_OUTRANGE;   /* out S: <TMIN | >TB4 *
	if (h<0 || h>7374.31) return H_OUTRANGE; /* out H: <TMIN & >TB4 *

	if (s<3.3978) /* TB2, PB100 *
	{
		if (h>polynom (s, 4, p100)) return H_OUTRANGE;  /* >PB100 *
		if (h<h_B4_s(s)) return AREA4;
		if (h<(s*2500.9/9.1555)) return H_OUTRANGE; /* area4, <PB0 *
		return AREA1;
	}

	if (s<3.7783) /* t4ss1(TB2) *
	{
		if (h>polynom (s, 4, p100)) return H_OUTRANGE;  /* >PB100 *
		if (h>polynom(s, 3, b13)) return AREA3;
		if (h<h_B4_s(s)) return AREA4;
		if (h<(s*2500.9/9.1555)) return H_OUTRANGE; /* area4, <PB0 *
		return AREA1;
	}

	if (s<6.04) /* PB100, TB3 *
	{
		if (h>polynom (s, 4, p100)) return H_OUTRANGE;  /* >PB100 *
		if (h<h_B4_s(s)) return AREA4;
		if (h>=2812.94) return AREA2;
		else
		{
//			if (s>interpolate2(h, 6, &b23[0][0])&&hs3t(h,s)<0)
			if (hs3t(h,s)<0)
				return AREA2;
		}
		if (h<(s*2500.9/9.1555)) return H_OUTRANGE; /* area4, <PB0 *
		return AREA3;
	}

	if (s<7.40863) /* PB10, TB3 *
	{
		if (h>polynom (s, 4, t800u)) return AREAX;  /* TB3, >PB10 *
		if (h<h_B4_s(s)) return AREA4;
		if (h<(s*2500.9/9.1555)) return H_OUTRANGE;  /*area4, <PB0 *
		return AREA2;
	}

	if (s<9.1555) /* t4ss2(PB0) *
	{
		if (h<=(s*2500.9/9.1555)) return H_OUTRANGE; /* PB10up, >TB3 *
		if (h>polynom (s, 4, p10)) return AREAX;
		if (h>(4160.99186-46.37851*exp(-(pow(s,0.1)-1.22173)/0.00766)))
			return AREA5;
		if (h<h_B4_s(s)) return AREA4;
		return AREA2;
	}

	if (s<9.42372)  /* PB10, TB4 *
	{
		if (h>polynom (s, 4, p10)) return AREAX;  /* >PB10, TB3*
		if (h<=polynom (s, 4, p0)) return H_OUTRANGE; /* <PB0 *
		if (h>(4160.99186-46.37851*exp(-(pow(s,0.1)-1.22173)/0.00766)))
			return AREA5;
		return AREA2;
	}

	if (s<=13.9406)
	{
		if (h<=polynom (s, 4, p0)) return H_OUTRANGE; /* <PB0 *
		if (h>polynom (s, 4, p10)) return AREAX;  /* >PB10, TB3*
		if (h>(4160.99186-46.37851*exp(-(pow(s,0.1)-1.22173)/0.00766)))
			return AREA5;
		return AREA2;
	}

	return AREAX;
} */

/* Implementation with backward formulas */
DllExport int WINAPI  region_hs(double h, double s) /* implemented Aug29,2005 */
{
	double h_TB3_s1[5]={127793.2914, -78557.42435, 18316.36317, -1870.837868, 70.95556447};
	/*double p10[5]={3596.95945, -844.36753, 185.65815, -18.41495, 1.35144}; v.97*/
	double h_PB50_s[5]={-2670.42611976705, 1998.57331301341, -298.395156751131, 20.4450539722304, 0.321543274514163};
	double h_TB3_s2[5]={-41936.89848, 20928.10394, -3568.633582, 270.7620581, -7.709489713};

	double p0[5]={-1526.6994, 913.38289, -30.92447,
		-8.21447, 0.6493};
	double h_TB4_s[5] = {4046.003752, 1148.813257, -147.9685734, 8.43417948, -0.179487245};

	if (s<SMIN || s>13.9406) return S_OUTRANGE;   /* out S: <TMIN | >TB4 */

	if (s<=S1_B14)
	{
		if (h>ps1h(PB100,s) || h<ps1h(PMIN,s)) return H_OUTRANGE;
		else return AREA1;
	}

	if (s<3.397782995) /* p=100, t=TB2=623.15k */
	{
		if (h>ps1h(PB100,s) || h< (273.15*s)) return H_OUTRANGE;
		if (h<h_B4_s(s)) return AREA4;
		else return AREA1;
	}
	
	if (s<3.778281340) /* s'(t=623.15) */
	{
		if (h>ps3h(PB100, s) || h< (273.15*s)) return H_OUTRANGE;
		if (h>h_B13_s(s)) return AREA3;
		if (h<h_B4_s(s)) return AREA4;
		else return AREA1;
	}

	if (s<5.048096828)   /* Smin_B23 */
	{
		if (h>ps3h(PB100, s) || h< (273.15*s)) return H_OUTRANGE;
		if (h<h_B4_s(s)) return AREA4;
		else return AREA3;
	}

	if (s<5.260578707)  /* Smax_B23 */
	{
		if (h>ps3h(PB100, s) || h< (273.15*s)) return H_OUTRANGE;
		if (h<h_B4_s(s)) return AREA4;
		else
		{
			if (h<2.563592004e3) return AREA3;   /* Hmin_B23 */
			if (h>2.812942061e3) return AREA2;
			else
			{
				if (p_2c_hs(h,s) > p_B23_t(t_B23_hs(h,s))) return AREA3;
				else return AREA2;
			}
		}
	}

	if (s<6.0404836717)  /* p=100, t=TB3=800C */
	{
		if (h>ps2h(PB100, s) || h< (273.15*s)) return H_OUTRANGE;
		if (h<h_B4_s(s)) return AREA4;
		else return AREA2;
	}

	/*if (s<7.4086748935)  p =10, t=TB3=800C v.97*/
	if (s<6.52264231228405) /* p =50, t=TB3=800C */
	{
		if (h>polynom (s, 4, h_TB3_s1) || h< (273.15*s)) return AREAX;  /* TB3, >PB50 */
		if (h<h_B4_s(s)) return AREA4;
		else return AREA2;
	}

	if (s<8.66886626515753) /* PB50, TB4 */
	{
		if (h>polynom (s, 4, h_PB50_s) || h< (273.15*s)) return AREAX;
		if (h>polynom(s,4,h_TB3_s2)) return AREA5;
		if (h<h_B4_s(s)) return AREA4;
		return AREA2;
	}

	
	if (s<9.1557593956) /* t"(PB0) */
	{
		if (h>polynom (s, 4, h_TB4_s) || h< (273.15*s)) return AREAX;
		/* TB3 <PB10 s= (7.4086748935, 11.921055069) -(PB10, PMIN) */
/*		if (h>(4160.99186-46.37851*exp(-(pow(s,0.1)-1.22173)/0.00766)))*/ 
/*		if (h>polynom(h, 4, n_h_TB3_s)) return AREA5;*/
		if (h>4160.73787-240597239.7226*exp(-2.08773*s)) return AREA5;
		if (h<h_B4_s(s)) return AREA4;
		return AREA2;
	}

	 /*if (s<9.42372)  PB10, TB4 obsolete v.97
	{
		if (h>polynom (s, 4, p10)) return AREAX;  /* >(PB50, TB3)*
		if (h<=polynom (s, 4, p0)) return H_OUTRANGE; /* <PB0 *
		if (h>4160.73787-240597239.7226*exp(-2.08773*s)) return AREA5;
		return AREA2;
	}*/

	if (s<11.921055069) /* 11.921055069 (PMIN, TB3) area5 has problem when s>11.5 */
	{
		if (h>polynom (s, 4, h_TB4_s)) return AREAX;  /* >(PB50, TB3)*/
		if (h<=polynom (s, 4, p0)) return H_OUTRANGE; /* <PB0 */
		if (h>4160.73787-240597239.7226*exp(-2.08773*s)) return AREA5;
		return AREA2;
	}

	if (s<SMAX) /* 13.904956091 (PMIN, TB4) area5 has problem when s>11.5 */
	{
		if (h>polynom (s, 4, h_TB4_s)) return AREAX;  /* >(PB50, TB3)*/
		if (h<=polynom (s, 4, p0)) return H_OUTRANGE; /* <PB0 */
		return AREA5;
	}


	else return S_OUTRANGE;


	return 1;
}


DllExport int WINAPI  region_tv(double t, double v)
{
	if (t<TMIN) return T_OUTRANGE;

	if (t<=TB2)
	{
		if (v<pt1v(PB100, t)) return V_OUTRANGE;
		if (v<=t4vs1(t)) return AREA1;
		if (v<t4vs2(t)) return AREA4;
		if (v<pt2v(PB0,t)) return AREA2;
		return V_OUTRANGE;
	}

	if (t<TC)
	{

		if (v<1/pt3r(PB100, t)) return V_OUTRANGE;
		if (v<=t4vs1(t)) return AREA3;
		if (v<t4vs2(t)) return AREA4;
		if (v<pt2v(p_B23_t(t),t)) return AREA3;
		if (v<pt2v(PB0,t)) return AREA2;
		return V_OUTRANGE;
	}

	if (t<TB23)
	{
		if (v<1/pt3r(PB100, t)) return V_OUTRANGE;
		if (v<pt2v(p_B23_t(t),t)) return AREA3;
		if (v<pt2v(PB0,t)) return AREA2;
		return V_OUTRANGE;
	}

	if (t<=TB3)
	{
		if (v<pt2v(PB100, t)) return V_OUTRANGE;
		if (v<pt2v(PB0,t)) return AREA2;
		return V_OUTRANGE;
	}

	if (t<=TB4)
	{
		if (v<pt5v(PB50, t)) return AREAX;
		if (v<pt5v(PB0, t)) return AREA5;
		return V_OUTRANGE;
	}

	return T_OUTRANGE;
}

DllExport int WINAPI  region_pv(double p, double v)
{
	if (p<=PB0) return P_OUTRANGE; 

	if (p<=t4p(TB2))
	{
		if (v<pt1v(p,TMIN)) return V_OUTRANGE;
		if (v<=p4vs1(p)) return AREA1;
		if (v<p4vs2(p)) return AREA4;
		if (v<=pt2v(p, TB3)) return AREA2;
		if (v<=pt5v(p, TB4)) return AREA5;
		return AREAX;
	}

	if (p<PC)
	{
		if (v<pt1v(p,TMIN)) return V_OUTRANGE;
		if (v<=pt1v(p,TB2)) return AREA1;
		if (v<=p4vs1(p)) return AREA3;
		if (v<p4vs2(p)) return AREA4;
		if (v<pt2v(p,t_B23_p(p))) return AREA3;
		if (v<=pt2v(p, TB3)) return AREA2;
		if (v<=pt5v(p, TB4)) return AREA5;
		return AREAX;
	}

	if (p<=PB50)
	{
		if (v<pt1v(p,TMIN)) return V_OUTRANGE;
		if (v<=pt1v(p,TB2)) return AREA1;
		if (v<pt2v(p,t_B23_p(p))) return AREA3;
		if (v<=pt2v(p, TB3)) return AREA2;
		if (v<=pt5v(p, TB4)) return AREA5;
		return V_OUTRANGE;
	}

	if (p<=PB100)
	{
		if (v<pt1v(p,TMIN)) return V_OUTRANGE;
		if (v<=pt1v(p,TB2)) return AREA1;
		if (v<pt2v(p,t_B23_p(p))) return AREA3;
		if (v<=pt2v(p, TB3)) return AREA2;
		return AREAX;
	}

	return P_OUTRANGE;
}

DllExport int WINAPI region_vh(double v, double h)
{
	double a1[3]={2013.70605,2.2321,5.26456E-4};
	double a2a[4]={-3.2989E7,1.93324E8,-3.7774E8,2.46098E8};
	double a2b[4]={-767689.44123,4.12249E6,-7.40051E6,4.45837E6};
	/*double a4[4]={2281.28947,24411.71647,314719.29864,-833519.21069}; */
	double h_PB50_v[5] = {2151.803527, 1612.913137, 329.0537991, 74.93169686, -23.78819055};
	double a5[4]={-26566.48525,123952.90936,-168025.99726,76421.89719};
	double a6[3]={1865.72712,-634764.16058,3.01359E8};
/*	double a7[2][11]={{7.13531E-03,6.02744E-03,5.94560E-03,4.38581E-03,
		3.87250E-03,3.50072E-03,3.22366E-03,3.01006E-03,2.84071E-03,
		2.70224E-03,2.58552E-03},
		{2.62637E+03,2.61980E+03,2.60965E+03,2.60553E+03,2.61209E+03,
		2.62937E+03,2.65532E+03,2.68773E+03,2.72529E+03,2.76728E+03,
		2.81259E+03}
		}; */
	double a8[4]={-48247.60047,170290.95347,-184850.7348,67019.35309};

	if (v<=0.000957) return V_OUTRANGE; /* t=0 C */
	if (h<HMIN) return H_OUTRANGE;

	if (v<0.001) /* p0, area1 */
	{
		if (h<-0.04) return H_OUTRANGE;  /* area1, <TMIN */
		return AREA1;
	}

	if (v<=1.312E-3)   /* area1, t<TB2 */
	{
		if (h>polynom(pow(v,0.1),3,a2a)) return H_OUTRANGE; /* area1, >PB100 */
		if (h<=(-0.01197+v*11.96652)) return H_OUTRANGE; /* <PB0 */
		if (h<t4hs1(v4ts(v))) return AREA4;
		return AREA1;
	}

	if (v<=0.00316) /* <VC */
	{
		if (h>polynom(pow(v,0.1),3,a2b)) return H_OUTRANGE; /* >PB100 */
		if (h<=(-0.01197+v*11.96652)) return H_OUTRANGE; /* <PB0 */
		if (h<t4hs1(v4ts(v))) return AREA4;
		if (h<=polynom(v,2,a6)) return AREA1;
		return AREA3;
	}

	if (v<4.336E-3) /* v(p=PB100, t=TB3) */
	{
		if (h>polynom(pow(v,0.1),3,a2b)) return H_OUTRANGE; /* >PB100 */
		if (h<=(-0.01197+v*11.96652)) return H_OUTRANGE; /* <PB0 */
		if (h<t4hs2(v4ts(v))) return AREA4;
		if (h>2812.59) return AREA2;
		else
		{
//			if (interpolate2(v, 11, &a7[0][0]) && rh3t(1/v,h)<0)
			if (rh3t(1/v,h)<0)
				return AREA2;
		}
		return AREA3;
	}
		
	if (v<0.00907413302906931) /* v(PB50, TB3) */
	{
		if (h>polynom(pow(v,0.1),3,a5)) return AREAX; /* >TB3,area2 */
		if (h<=(-0.01197+v*11.96652)) return H_OUTRANGE; /* <PB0 */
		if (h<t4hs2(v4ts(v))) return AREA4;
		if (h>2812.59) return AREA2;
		else
		{
//			if (interpolate2(v, 11, &a7[0][0]) && rh3t(1/v,h)<0)
			if (rh3t(1/v,h)<0)
				return AREA2;
		}
		return AREA3;
	}
		
	if (v<0.0214633996483061) /* v(PB50, TB4) */
	{
		if (h>polynom(v*100,4,h_PB50_v)) return AREAX;  /* >PB50,area5 */
		if (h<=(-0.01197+v*11.96652)) return H_OUTRANGE; /* <PB0 */
		if (h>polynom(pow(v,0.05),3,a8)) return AREA5;
		if (h>=t4hs2(v4ts(v))) return AREA2;
		return AREA4;
	}

	if (v<209.7779) /* v4s2, PB0 */
	{
		if (h<=(-0.01197+v*11.96652)) return H_OUTRANGE; /* <PB0 area4*/
		if (h>7374.75) return H_OUTRANGE;   /* >TB4 */
		if (h<t4hs2(v4ts(v))) return AREA4;
		if (v>=9.9 && h>4160.5) return AREA5;
		if (v<=9.9 &&  
			h>(4160.40518-45.60305*exp(-(pow(v,0.05)-0.85967)/0.04398))) return AREA5;
		return AREA2;
	}

	if (v<809.745) /* PB0, area2,5 crosspoint */
	{
		if (h>7376.75) return H_OUTRANGE;   /* >TB4 */
		if (h<=polynom(v,2,a1)) return H_OUTRANGE; /* <PB0 area2,5*/
		if (h>(4159.87887-45.30513*exp(-(pow(v,0.05)-0.85967)/0.04265))) return AREA5;
		return AREA2;
	}

	if (v<1715.2066)
	{
		if (h>7376.75) return H_OUTRANGE;   /* >TB4 */
		if (h<=polynom(v,2,a1)) return H_OUTRANGE; /* <PB0 area2,5*/
		if (h>(4159.87887-45.30513*exp(-(pow(v,0.05)-0.85967)/0.04265))) return AREA5;
		return AREA5;
	}

	return AREAX;
}

DllExport int WINAPI  region_vs(double v, double s)
{
	double a1[5]={7.29842,0.01115,-1.06067E-5,5.79591E-9,-1.22187E-12};
	double a2[5]={-455.38366,483391.16241,-2.05412E7,8.09063E8};
	double s_TB4_v[5]={10.46775638, 0.461326522, -0.000583286, 0.000224336, -1.93992E-05};
	double s_PB50_v[5]={4.120963826, 3.213218726, -0.733980191, 0.133208984, -0.013414858};

	double a5[4]={-27.844,130.36461,-168.53277,76.69427};
	double a6[3]={0.56541,3123.29036,-734367.09901};
	double a8[3]={-4.57423,17.44734,-4.04929};			
		
	if (v<=0.000957) return V_OUTRANGE; /* t=0 C */
	if (s<SMIN) return S_OUTRANGE;

	if (v<0.001) /* p0, area1 */
	{
		if (s<-0.0002) return S_OUTRANGE;  /* area1, <TMIN */
		else return AREA1;
	}

	if (v<=1.312E-3)   /* area1, t<=TB2 */
	{
//		if (s>(3.87651-2.61707*exp(-(v-1E-3)/1.87591E-4))) return S_OUTRANGE; /* area1, >PB100 */
		if (s>polynom(pow(v,0.1),4,a2)) return S_OUTRANGE; /* >PB100 */
		if (s<=(-4.4445E-5+v*0.04445)) return S_OUTRANGE; /* <PB0 */
//		if (s<=0.0005) return S_OUTRANGE; /* <PB0 */
		if (s<t4ss1(v4ts(v))) return AREA4;
		return AREA1;
	}

	if (v<=0.00316) /* <VC */
	{
		if (s>polynom(pow(v,0.1),4,a2)) return S_OUTRANGE-10; /* >PB100 */
		if (s<=(-4.4445E-5+v*0.04445)) return S_OUTRANGE-20; /* <PB0 */
		if (s<t4ss1(v4ts(v))) return AREA4;
		if (s<=polynom(v,2,a6)) return AREA1;
		return AREA3;
	}

	if (v<4.336E-3) /* v(p=PB100, t=TB3) */
	{
		if (s>polynom(pow(v,0.1),4,a2)) return S_OUTRANGE-10; /* >PB100 */
		if (s<=(-4.4445E-5+v*0.04445)) return S_OUTRANGE-20; /* <PB0 */
		if (s<t4ss2(v4ts(v))) return AREA4;
//		if (s> vbd8s(v) && rs3p(1/v,s)<0) return AREA2;
		if (rs3p(1/v,s)<0) return AREA2;
		return AREA3;
	}
	
	if (v<0.00907413302906931) /* v(PB50, TB3) */
	{
		if (s>polynom(pow(v,0.1),3,a5)) return AREAX; /* >TB3,area2 */
		if (s<=(-4.4445E-5+v*0.04445)) return S_OUTRANGE; /* <PB0 */
		if (s<t4ss2(v4ts(v))) return AREA4;
		if (v>0.0088066) return AREA2;
		else 
		{
//			if (s> vbd8s(v) && rs3p(1/v,s)<0) return AREA2;
			if (rs3p(1/v,s)<0) return AREA2;
		}
		return AREA3;
	}
		
	if (v<0.0214633996483061) /* v(PB50, TB4) */
	{
		if (s>polynom(v*100,4,s_PB50_v)) return AREAX;  /* >PB50,area5 */
		if (s<=(-4.4445E-5+v*0.04445)) return S_OUTRANGE; /* <PB0 */
		if (s<t4ss2(v4ts(v))) return AREA2;
		if (s>polynom(pow(v,0.05),2,a8)) return AREA5;
		return AREA4+430;
	}

	if (v<209.7779) /* point at area4-s2*PB0 */
	{
		if (s<polynom(v,4,a1)) return S_OUTRANGE; /* <PB0 area4*/
		if (s>polynom(log(v),4,s_TB4_v)) return S_OUTRANGE;   /* >TB4 */
		if (s<t4ss2(v4ts(v))) return AREA4+410;
		if (s>polynom(pow(v,0.05),2,a8)) return AREA5;
		return AREA2;
	}

	if (v<809.745) /* PB0, area2,5 crosspoint */
	{
		if (s<polynom(v,4,a1)) return S_OUTRANGE; /* <PB0 area4*/
		if (s>polynom(log(v),4,s_TB4_v)) return S_OUTRANGE;   /* >TB4 */
		if (s>polynom(pow(v,0.05),2,a8)) return AREA5;
		return AREA2;
	}

	if (v<1715.2066)
	{
		if (s<polynom(v,4,a1)) return S_OUTRANGE; /* <PB0 area4*/
		if (s>polynom(log(v),4,s_TB4_v)) return S_OUTRANGE;   /* >TB4 */
		return AREA5;
	}

	return AREAX;
}

DllExport int WINAPI  region_px(double p, double x)
{
	if (x<0 || x>1) return X_OUTRANGE;
	if (p<Pt || p>PC) return P_OUTRANGE;
	return AREA4;
}

DllExport int WINAPI  region_tx(double t, double x)
{
	if (x<0 || x>1) return X_OUTRANGE;
	if (t<Tt || t>TC) return T_OUTRANGE;
	return AREA4;
}

DllExport int WINAPI  region_vx(double v, double x)
{
	if (x<0 || x>1) return X_OUTRANGE;
	if (v<0.001 || v>206.14) return V_OUTRANGE;
	return AREA4;
}

DllExport int WINAPI  region_hx(double h, double x)
{
	if (x<0 || x>1) return X_OUTRANGE;
	if (h<-0.042 || h>2803.28) return H_OUTRANGE;
	return AREA4;
}

DllExport int WINAPI  region_sx(double s, double x)
{
	if (x<0 || x>1) return X_OUTRANGE;
	if (s<-0.0002 || s>9.1558) return S_OUTRANGE;
	return AREA4;
}


