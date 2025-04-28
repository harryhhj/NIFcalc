/*
Filename: generic.cpp
generic function regarding to areas.
*/
#include "stdafx.h"

#include "generic.h"
#include "areas.h"
#include "area1.h"
#include "area2.h"
#include "area3.h"
#include "area4.h"
#include "area5.h"
#include "region.h"

/* Given p, t */
DllExport double WINAPI pt_v(double p, double t, int area)
{	
	double ret;

	if (area==0) area = region_pt(p, t);

	switch (area)
	{
		case AREA1: ret = pt1v(p,t);	break;
		case AREA2: ret = pt2v(p,t);	break;
		case AREA3: ret = 1/pt3r(p,t);	break;
		case AREA5: ret = pt5v(p,t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1)
		return (ret+area); 
	return ret;
}

DllExport double WINAPI pt_u(double p, double t, int area)
{	
	double ret;

	if (area==0) area = region_pt(p, t);

	switch (area)
	{
		case AREA1: ret = pt1u(p,t);	break;
		case AREA2: ret = pt2u(p,t);	break;
		case AREA3: ret = rt3u(pt3r(p,t),t);	break;
		case AREA5: ret = pt5u(p,t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pt_h(double p, double t, int area)
{	
	double ret;

	if (area==0) area = region_pt(p, t);

	switch (area)
	{
		case AREA1: ret = pt1h(p,t);	break;
		case AREA2: ret = pt2h(p,t);	break;
		case AREA3: ret = rt3h(pt3r(p,t),t);	break;
		case AREA5: ret = pt5h(p,t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pt_s(double p, double t, int area)
{	
	double ret;

	if (area==0) area = region_pt(p, t);

	switch (area)
	{
		case AREA1: ret = pt1s(p,t);	break;
		case AREA2: ret = pt2s(p,t);	break;
		case AREA3: ret = rt3s(pt3r(p,t),t);	break;
		case AREA5: ret = pt5s(p,t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pt_cp(double p, double t, int area)
{	
	double ret;

	if (area==0) area = region_pt(p, t);

	switch (area)
	{
		case AREA1: ret = pt1cp(p,t);	break;
		case AREA2: ret = pt2cp(p,t);	break;
		case AREA3: ret = rt3cp(pt3r(p,t),t);	break;
		case AREA5: ret = pt5cp(p,t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pt_cv(double p, double t, int area)
{	
	double ret;

	if (area==0) area = region_pt(p, t);

	switch (area)
	{
		case AREA1: ret = pt1cv(p,t);	break;
		case AREA2: ret = pt2cv(p,t);	break;
		case AREA3: ret = rt3cv(pt3r(p,t),t);	break;
		case AREA5: ret = pt5cv(p,t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pt_w(double p, double t, int area)
{	
	double ret;

	if (area==0) area = region_pt(p, t);

	switch (area)
	{
		case AREA1: ret = pt1w(p,t);	break;
		case AREA2: ret = pt2w(p,t);	break;
		case AREA3: ret = rt3w(pt3r(p,t),t);	break;
		case AREA5: ret = pt5w(p,t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}


/* Given p, h */
DllExport double WINAPI ph_t(double p, double h,  int area)
{	
	double ret;

	if (area==0) area = region_ph(p, h);

	switch (area)
	{
		case AREA1: ret = ph1t(p,h);	break;
		case AREA2: ret = ph2t(p,h);	break;
		case AREA3: ret = ph3t(p,h);	break;
		case AREA4: ret = ts(p);		break;
		case AREA5: ret = ph5t(p,h);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ph_v(double p, double h, int area)
{	
	double ret, t;

	if (area==0) area = region_ph(p, h);
	
	switch (area)
	{
		case AREA4:
			t = ts(p);
			ret = tx4v(t, th4x(t,h));			break;
		case AREA1: ret = pt1v(p,ph1t(p,h));	break;
		case AREA2: ret = pt2v(p,ph2t(p,h));	break;
		case AREA3: ret = ph3v(p,h);			break;
		case AREA5: ret = pt5v(p,ph5t(p,h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ph_u(double p, double h,  int area)
{	
	double ret, t;

	if (area==0) area = region_ph(p, h);
	
	switch (area)
	{
		case AREA4:
			t = ts(p);
			ret = tx4u(t, th4x(t,h));			break;
		case AREA1: ret = pt1u(p,ph1t(p,h));	break;
		case AREA2: ret = pt2u(p,ph2t(p,h));	break;
		case AREA3: ret = rt3u(1/ph3v(p,h), ph3t(p,h));	break;
		case AREA5: ret = pt5u(p,ph5t(p,h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ph_s(double p, double h,  int area)
{	
	double ret, t;

	if (area==0) area = region_ph(p, h);
	
	switch (area)
	{
		case AREA4:
			t = ts(p);
			ret = tx4s(t, th4x(t,h));			break;
		case AREA1: ret = pt1s(p,ph1t(p,h));	break;
		case AREA2: ret = pt2s(p,ph2t(p,h));	break;
		case AREA3: ret = rt3s(1/ph3v(p,h), ph3t(p,h));	break;
		case AREA5: ret = pt5s(p,ph5t(p,h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ph_cp(double p, double h,  int area)
{	
	double ret, t;

	if (area==0) area = region_ph(p, h);
	
	switch (area)
	{
		case AREA4:
			t = ts(p);
			ret = tx4cp(t, th4x(t,h));			break;
		case AREA1: ret = pt1cp(p,ph1t(p,h));	break;
		case AREA2: ret = pt2cp(p,ph2t(p,h));	break;
		case AREA3: ret = rt3cp(1/ph3v(p,h), ph3t(p,h));	break;
		case AREA5: ret = pt5cp(p,ph5t(p,h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ph_cv(double p, double h,  int area)
{	
	double ret, t;

	if (area==0) area = region_ph(p, h);
	
	switch (area)
	{
		case AREA4:
			t = ts(p);
			ret = tx4cv(t, th4x(t,h));			break;
		case AREA1: ret = pt1cv(p,ph1t(p,h));	break;
		case AREA2: ret = pt2cv(p,ph2t(p,h));	break;
		case AREA3: ret = rt3cv(1/ph3v(p,h), ph3t(p,h));	break;
		case AREA5: ret = pt5cv(p,ph5t(p,h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ph_w(double p, double h,  int area)
{	
	double ret, t;

	if (area==0) area = region_ph(p, h);
	
	switch (area)
	{
		case AREA4:
			t = ts(p);
			ret = tx4w(t, th4x(t,h));			break;
		case AREA1: ret = pt1w(p,ph1t(p,h));	break;
		case AREA2: ret = pt2w(p,ph2t(p,h));	break;
		case AREA3: ret = rt3w(1/ph3v(p,h), ph3t(p,h));	break;
		case AREA5: ret = pt5w(p,ph5t(p,h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

/* Given p, s */
DllExport double WINAPI ps_t(double p, double s, int area)
{	
	double ret;

	if (area==0) area = region_ps(p, s);

	switch (area)
	{
		case AREA1: ret = ps1t(p,s);	break;
		case AREA2: ret = ps2t(p,s);	break;
		case AREA3: ret = ps3t(p,s);	break;
		case AREA4: ret = ts(p);		break;
		case AREA5: ret = ps5t(p,s);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ps_v(double p, double s, int area)
{	
	double ret, t;

	if (area==0) area = region_ps(p, s);
	
	switch (area)
	{
		case AREA4: 
			t=ts(p);
			ret = tx4v(t, ts4x(t,s));			break;
		case AREA1: ret = pt1v(p,ps1t(p,s));	break;
		case AREA2: ret = pt2v(p,ps2t(p,s));	break;
		case AREA3: ret = ps3v(p,s);			break;
		case AREA5: ret = pt5v(p,ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ps_u(double p, double s, int area)
{	
	double ret, t;

	if (area==0) area = region_ps(p, s);
	
	switch (area)
	{
		case AREA4:
			t=ts(p);
			ret = tx4u(t, ts4x(t,s));			break;
		case AREA1: ret = pt1u(p,ps1t(p,s));	break;
		case AREA2: ret = pt2u(p,ps2t(p,s));	break;
		case AREA3: ret = rt3u(1/ps3v(p,s), ps3t(p,s));	break;
		case AREA5: ret = pt5u(p,ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ps_h(double p, double s, int area)
{	
	double ret, t;

	if (area==0) area = region_ps(p, s);
	
	switch (area)
	{
		case AREA4:
			t=ts(p);
			ret = tx4h(t, ts4x(t,s));			break;
		case AREA1: ret = pt1h(p,ps1t(p,s));	break;
		case AREA2: ret = pt2h(p,ps2t(p,s));	break;
		case AREA3: ret = rt3h(1/ps3v(p,s), ps3t(p,s));	break;
		case AREA5: ret = pt5h(p,ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ps_cp(double p, double s, int area)
{	
	double ret, t;

	if (area==0) area = region_ps(p, s);
	
	switch (area)
	{
		case AREA4:
			t=ts(p);
			ret = tx4cp(t, ts4x(t,s));			break;
		case AREA1: ret = pt1cp(p,ps1t(p,s));	break;
		case AREA2: ret = pt2cp(p,ps2t(p,s));	break;
		case AREA3: ret = rt3cp(1/ps3v(p,s), ps3t(p,s));	break;
		case AREA5: ret = pt5cp(p,ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ps_cv(double p, double s, int area)
{	
	double ret, t;

	if (area==0) area = region_ps(p, s);
	
	switch (area)
	{
		case AREA4:
			t=ts(p);
			ret = tx4cv(t, ts4x(t,s));			break;
		case AREA1: ret = pt1cv(p,ps1t(p,s));	break;
		case AREA2: ret = pt2cv(p,ps2t(p,s));	break;
		case AREA3: ret = rt3cv(1/ps3v(p,s), ps3t(p,s));	break;
		case AREA5: ret = pt5cv(p,ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ps_w(double p, double s, int area)
{	
	double ret, t;

	if (area==0) area = region_ps(p, s);
	
	switch (area)
	{
		case AREA4:
			t=ts(p);
			ret = tx4w(t, ts4x(t,s));			break;
		case AREA1: ret = pt1w(p,ps1t(p,s));	break;
		case AREA2: ret = pt2w(p,ps2t(p,s));	break;
		case AREA3: ret = rt3w(1/ps3v(p,s), ps3t(p,s));	break;
		case AREA5: ret = pt5w(p,ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

/* Given t, h */
DllExport double WINAPI th_p(double t, double h, int area)
{	
	double ret;

	if (area==0) area = region_th(t, h);
	
	switch (area)
	{
		case AREA4: ret = ps(t);		break;
		case AREA1: ret = th1p(t,h);	break;
		case AREA2: ret = th2p(t,h);	break;
		case AREA3: ret = rt3p(th3r(t,h), t);	break;
		case AREA5: ret = th5p(t,h);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI th_v(double t, double h, int area)
{	
	double ret;

	if (area==0) area = region_th(t, h);
	
	switch (area)
	{
		case AREA4: ret = tx4v(t,th4x(t,h));	break;
		case AREA1: ret = pt1v(th1p(t,h),t);	break;
		case AREA2: ret = pt2v(th2p(t,h),t);	break;
		case AREA3: ret = 1/th3r(t,h);			break;
		case AREA5: ret = pt5v(th5p(t,h),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI th_u(double t, double h, int area)
{	
	double ret;

	if (area==0) area = region_th(t, h);
	
	switch (area)
	{
		case AREA4: ret = tx4u(t,th4x(t,h));	break;
		case AREA1: ret = pt1u(th1p(t,h),t);	break;
		case AREA2: ret = pt2u(th2p(t,h),t);	break;
		case AREA3: ret = rt3u(th3r(t,h),t);	break;
		case AREA5: ret = pt5u(th5p(t,h),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI th_s(double t, double h, int area)
{	
	double ret;

	if (area==0) area = region_th(t, h);
	
	switch (area)
	{
		case AREA4: ret = tx4s(t,th4x(t,h));	break;
		case AREA1: ret = pt1s(th1p(t,h),t);	break;
		case AREA2: ret = pt2s(th2p(t,h),t);	break;
		case AREA3: ret = rt3s(th3r(t,h),t);	break;
		case AREA5: ret = pt5s(th5p(t,h),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI th_cp(double t, double h, int area)
{	
	double ret;

	if (area==0) area = region_th(t, h);
	
	switch (area)
	{
		case AREA4: ret = tx4cp(t,th4x(t,h));	break;
		case AREA1: ret = pt1cp(th1p(t,h),t);	break;
		case AREA2: ret = pt2cp(th2p(t,h),t);	break;
		case AREA3: ret = rt3cp(th3r(t,h),t);	break;
		case AREA5: ret = pt5cp(th5p(t,h),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI th_cv(double t, double h, int area)
{	
	double ret;

	if (area==0) area = region_th(t, h);
	
	switch (area)
	{
		case AREA4: ret = tx4cv(t,th4x(t,h));	break;
		case AREA1: ret = pt1cv(th1p(t,h),t);	break;
		case AREA2: ret = pt2cv(th2p(t,h),t);	break;
		case AREA3: ret = rt3cv(th3r(t,h),t);	break;
		case AREA5: ret = pt5cv(th5p(t,h),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI th_w(double t, double h, int area)
{	
	double ret;

	if (area==0) area = region_th(t, h);
	
	switch (area)
	{
		case AREA4: ret = tx4w(t,th4x(t,h));	break;
		case AREA1: ret = pt1w(th1p(t,h),t);	break;
		case AREA2: ret = pt2w(th2p(t,h),t);	break;
		case AREA3: ret = rt3w(th3r(t,h),t);	break;
		case AREA5: ret = pt5w(th5p(t,h),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

/* Given t, s */
DllExport double WINAPI ts_p(double t, double s, int area)
{	
	double ret;

	if (area==0) area = region_ts(t, s);
	
	switch (area)
	{
		case AREA4: ret = ps(t);		break;
		case AREA1: ret = ts1p(t,s);	break;
		case AREA2: ret = ts2p(t,s);	break;
		case AREA3: ret = rt3p(ts3r(t,s), t);	break;
		case AREA5: ret = ts5p(t,s);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ts_v(double t, double s, int area)
{
	double ret;

	if (area==0) area = region_ts(t, s);
	
	switch (area)
	{
		case AREA4: ret = tx4v(t,ts4x(t,s));	break;
		case AREA1: ret = pt1v(ts1p(t,s),t);	break;
		case AREA2: ret = pt2v(ts2p(t,s),t);	break;
		case AREA3: ret = 1/ts3r(t,s);			break;
		case AREA5: ret = pt5v(ts5p(t,s),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ts_u(double t, double s, int area)
{	
	double ret;

	if (area==0) area = region_ts(t, s);
	
	switch (area)
	{
		case AREA4: ret = tx4u(t,ts4x(t,s));	break;
		case AREA1: ret = pt1u(ts1p(t,s),t);	break;
		case AREA2: ret = pt2u(ts2p(t,s),t);	break;
		case AREA3: ret = rt3u(ts3r(t,s),t);	break;
		case AREA5: ret = pt5u(ts5p(t,s),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ts_h(double t, double s, int area)
{	
	double ret;

	if (area==0) area = region_ts(t, s);
	
	switch (area)
	{
		case AREA4: ret = tx4h(t,ts4x(t,s));	break;
		case AREA1: ret = pt1h(ts1p(t,s),t);	break;
		case AREA2: ret = pt2h(ts2p(t,s),t);	break;
		case AREA3: ret = rt3h(ts3r(t,s),t);	break;
		case AREA5: ret = pt5h(ts5p(t,s),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ts_cp(double t, double s, int area)
{	
	double ret;

	if (area==0) area = region_ts(t, s);
	
	switch (area)
	{
		case AREA4: ret = tx4cp(t,ts4x(t,s));	break;
		case AREA1: ret = pt1cp(ts1p(t,s),t);	break;
		case AREA2: ret = pt2cp(ts2p(t,s),t);	break;
		case AREA3: ret = rt3cp(ts3r(t,s),t);	break;
		case AREA5: ret = pt5cp(ts5p(t,s),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ts_cv(double t, double s, int area)
{	
	double ret;

	if (area==0) area = region_ts(t, s);
	
	switch (area)
	{
		case AREA4: ret = tx4cv(t,ts4x(t,s));	break;
		case AREA1: ret = pt1cv(ts1p(t,s),t);	break;
		case AREA2: ret = pt2cv(ts2p(t,s),t);	break;
		case AREA3: ret = rt3cv(ts3r(t,s),t);	break;
		case AREA5: ret = pt5cv(ts5p(t,s),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ts_w(double t, double s, int area)
{	
	double ret;

	if (area==0) area = region_ts(t, s);
	
	switch (area)
	{
		case AREA4: ret = tx4w(t,ts4x(t,s));	break;
		case AREA1: ret = pt1w(ts1p(t,s),t);	break;
		case AREA2: ret = pt2w(ts2p(t,s),t);	break;
		case AREA3: ret = rt3w(ts3r(t,s),t);	break;
		case AREA5: ret = pt5w(ts5p(t,s),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

/* Given h, s */
DllExport double WINAPI hs_p(double h, double s, int area)
{	
	double ret;

	if (area==0) area = region_hs(h, s);
	
	switch (area)
	{
		case AREA4: ret = ps(hs4t(h,s));	break;
		case AREA1: ret = hs1p(h,s);	break;
		case AREA2: ret = hs2p(h,s);	break;
		case AREA3: ret = hs3p(h,s);	break;
		case AREA5: ret = hs5p(h,s);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI hs_t(double h, double s, int area)
{	
	double ret;

	if (area==0) area = region_hs(h, s);
	
	switch (area)
	{
		case AREA4: ret = hs4t(h,s);	break;
		case AREA1: ret = hs1t(h,s);	break;
		case AREA2: ret = hs2t(h,s);	break;
		case AREA3: ret = hs3t(h,s);	break;
		case AREA5: ret = hs5t(h,s);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI hs_v(double h, double s, int area)
{	
	double ret;

	if (area==0) area = region_hs(h, s);
	
	switch (area)
	{
		case AREA4: ret = tx4v(hs4t(h,s), hs4x(h,s));	break;
		case AREA1: ret = pt1v(hs1p(h,s), hs1t(h,s));	break;
		case AREA2: ret = pt2v(hs2p(h,s), hs2t(h,s));	break;
		case AREA3: ret = hs3v(h,s);					break;
		case AREA5: ret = pt5v(hs5p(h,s), hs5t(h,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI hs_u(double h, double s, int area)
{	
	double ret;

	if (area==0) area = region_hs(h, s);
	
	switch (area)
	{
		case AREA4: ret = tx4u(hs4t(h,s), hs4x(h,s));	break;
		case AREA1: ret = pt1u(hs1p(h,s), hs1t(h,s));	break;
		case AREA2: ret = pt2u(hs2p(h,s), hs2t(h,s));	break;
		case AREA3: ret = rt3u(1/hs3v(h,s), hs3t(h,s));	break;
		case AREA5: ret = pt5u(hs5p(h,s), hs5t(h,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI hs_cp(double h, double s, int area)
{	
	double ret;

	if (area==0) area = region_hs(h, s);
	
	switch (area)
	{
		case AREA4: ret = tx4cp(hs4t(h,s), hs4x(h,s));	break;
		case AREA1: ret = pt1cp(hs1p(h,s), hs1t(h,s));	break;
		case AREA2: ret = pt2cp(hs2p(h,s), hs2t(h,s));	break;
		case AREA3: ret = rt3cp(1/hs3v(h,s), hs3t(h,s));	break;
		case AREA5: ret = pt5cp(hs5p(h,s), hs5t(h,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI hs_cv(double h, double s, int area)
{	
	double ret;

	if (area==0) area = region_hs(h, s);
	
	switch (area)
	{
		case AREA4: ret = tx4cv(hs4t(h,s), hs4x(h,s));	break;
		case AREA1: ret = pt1cv(hs1p(h,s), hs1t(h,s));	break;
		case AREA2: ret = pt2cv(hs2p(h,s), hs2t(h,s));	break;
		case AREA3: ret = rt3cv(1/hs3v(h,s), hs3t(h,s));	break;
		case AREA5: ret = pt5cv(hs5p(h,s), hs5t(h,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI hs_w(double h, double s, int area)
{	
	double ret;

	if (area==0) area = region_hs(h, s);
	
	switch (area)
	{
		case AREA4: ret = tx4w(hs4t(h,s), hs4x(h,s));	break;
		case AREA1: ret = pt1w(hs1p(h,s), hs1t(h,s));	break;
		case AREA2: ret = pt2w(hs2p(h,s), hs2t(h,s));	break;
		case AREA3: ret = rt3w(1/hs3v(h,s), hs3t(h,s));	break;
		case AREA5: ret = pt5w(hs5p(h,s), hs5t(h,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

/* Given t, v */
DllExport double WINAPI tv_p(double t, double v, int area)
{	
	double ret;

	if (area==0) area = region_tv(t, v);
	
	switch (area)
	{
		case AREA4: ret = ps(t);		break;
		case AREA1: ret = tv1p(t,v);	break;
		case AREA2: ret = tv2p(t,v);	break;
		case AREA3: ret = rt3p(1/v, t);	break;
		case AREA5: ret = tv5p(t,v);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tv_u(double t, double v, int area)
{	
	double ret;

	if (area==0) area = region_tv(t, v);
	
	switch (area)
	{
		case AREA4: ret = tx4u(t,tv4x(t,v));	break;
		case AREA1: ret = pt1u(tv1p(t,v),t);	break;
		case AREA2: ret = pt2u(tv2p(t,v),t);	break;
		case AREA3: ret = rt3u(1/v, t);			break;
		case AREA5: ret = pt5u(tv5p(t,v),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tv_h(double t, double v, int area)
{	
	double ret;

	if (area==0) area = region_tv(t, v);
	
	switch (area)
	{
		case AREA4: ret = tx4h(t,tv4x(t,v));	break;
		case AREA1: ret = pt1h(tv1p(t,v),t);	break;
		case AREA2: ret = pt2h(tv2p(t,v),t);	break;
		case AREA3: ret = rt3h(1/v, t);			break;
		case AREA5: ret = pt5h(tv5p(t,v),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tv_s(double t, double v, int area)
{	
	double ret;

	if (area==0) area = region_tv(t, v);
	
	switch (area)
	{
		case AREA4: ret = tx4s(t,tv4x(t,v));	break;
		case AREA1: ret = pt1s(tv1p(t,v),t);	break;
		case AREA2: ret = pt2s(tv2p(t,v),t);	break;
		case AREA3: ret = rt3s(1/v, t);			break;
		case AREA5: ret = pt5s(tv5p(t,v),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tv_cp(double t, double v, int area)
{	
	double ret;

	if (area==0) area = region_tv(t, v);
	
	switch (area)
	{
		case AREA4: ret = tx4cp(t,tv4x(t,v));	break;
		case AREA1: ret = pt1cp(tv1p(t,v),t);	break;
		case AREA2: ret = pt2cp(tv2p(t,v),t);	break;
		case AREA3: ret = rt3cp(1/v, t);		break;
		case AREA5: ret = pt5cp(tv5p(t,v),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tv_cv(double t, double v, int area)
{	
	double ret;

	if (area==0) area = region_tv(t, v);
	
	switch (area)
	{
		case AREA4: ret = tx4cv(t,tv4x(t,v));	break;
		case AREA1: ret = pt1cv(tv1p(t,v),t);	break;
		case AREA2: ret = pt2cv(tv2p(t,v),t);	break;
		case AREA3: ret = rt3cv(1/v, t);		break;
		case AREA5: ret = pt5cv(tv5p(t,v),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tv_w(double t, double v, int area)
{	
	double ret;

	if (area==0) area = region_tv(t, v);
	
	switch (area)
	{
		case AREA4: ret = tx4w(t,tv4x(t,v));	break;
		case AREA1: ret = pt1w(tv1p(t,v),t);	break;
		case AREA2: ret = pt2w(tv2p(t,v),t);	break;
		case AREA3: ret = rt3w(1/v, t);			break;
		case AREA5: ret = pt5w(tv5p(t,v),t);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}


/* Given p, v */
DllExport double WINAPI pv_t(double p, double v, int area)
{	
	double ret;

	if (area==0) area = region_pv(p, v);

	switch (area)
	{
		case AREA1: ret = pv1t(p,v);	break;
		case AREA2: ret = pv2t(p,v);	break;
		case AREA3: ret = rp3t(1/v,p);	break;
		case AREA4: ret = ts(p);		break;
		case AREA5: ret = pv5t(p,v);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pv_u(double p, double v, int area)
{	
	double ret, t;

	if (area==0) area = region_pv(p, v);

	switch (area)
	{
		case AREA1: ret = pt1u(p, pv1t(p,v));	break;
		case AREA2: ret = pt2u(p, pv2t(p,v));	break;
		case AREA3: ret = rt3u(1/v, rp3t(1/v,p));	break;
		case AREA4:
			t = ts(p);
			ret = tx4u(t, tv4x(t,v));			break;
		case AREA5: ret = pt5u(p, pv5t(p,v));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pv_h(double p, double v, int area)
{	
	double ret, t;

	if (area==0) area = region_pv(p, v);

	switch (area)
	{
		case AREA1: ret = pt1h(p, pv1t(p,v));	break;
		case AREA2: ret = pt2h(p, pv2t(p,v));	break;
		case AREA3: ret = rt3h(1/v, rp3t(1/v,p));	break;
		case AREA4:
			t = ts(p);
			ret = tx4h(t, tv4x(t,v));			break;
		case AREA5: ret = pt5h(p, pv5t(p,v));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pv_s(double p, double v, int area)
{	
	double ret, t;

	if (area==0) area = region_pv(p, v);

	switch (area)
	{
		case AREA1: ret = pt1s(p, pv1t(p,v));	break;
		case AREA2: ret = pt2s(p, pv2t(p,v));	break;
		case AREA3: ret = rt3s(1/v, rp3t(1/v,p));	break;
		case AREA4:
			t = ts(p);
			ret = tx4s(t, tv4x(t,v));			break;
		case AREA5: ret = pt5s(p, pv5t(p,v));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pv_cp(double p, double v, int area)
{	
	double ret, t;

	if (area==0) area = region_pv(p, v);

	switch (area)
	{
		case AREA1: ret = pt1cp(p, pv1t(p,v));	break;
		case AREA2: ret = pt2cp(p, pv2t(p,v));	break;
		case AREA3: ret = rt3cp(1/v, rp3t(1/v,p));	break;
		case AREA4:
			t = ts(p);
			ret = tx4cp(t, tv4x(t,v));			break;
		case AREA5: ret = pt5cp(p, pv5t(p,v));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pv_cv(double p, double v, int area)
{	
	double ret, t;

	if (area==0) area = region_pv(p, v);

	switch (area)
	{
		case AREA1: ret = pt1cv(p, pv1t(p,v));	break;
		case AREA2: ret = pt2cv(p, pv2t(p,v));	break;
		case AREA3: ret = rt3cv(1/v, rp3t(1/v,p));	break;
		case AREA4:
			t = ts(p);
			ret = tx4cv(t, tv4x(t,v));			break;
		case AREA5: ret = pt5cv(p, pv5t(p,v));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI pv_w(double p, double v, int area)
{	
	double ret, t;

	if (area==0) area = region_pv(p, v);

	switch (area)
	{
		case AREA1: ret = pt1w(p, pv1t(p,v));	break;
		case AREA2: ret = pt2w(p, pv2t(p,v));	break;
		case AREA3: ret = rt3w(1/v, rp3t(1/v,p));	break;
		case AREA4:
			t = ts(p);
			ret = tx4w(t, tv4x(t,v));			break;
		case AREA5: ret = pt5w(p, pv5t(p,v));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}


/* Given v, h */
DllExport double WINAPI vh_p(double v, double h, int area)
{	
	double ret;

	if (area==0) area = region_vh(v, h);
	
	switch (area)
	{
		case AREA4: ret = ps(vh4t(v,h));	break;
		case AREA1: ret = vh1p(v,h);	break;
		case AREA2: ret = vh2p(v,h);	break;
		case AREA3: ret = rt3p(1/v,rh3t(1/v, h));	break;
		case AREA5: ret = vh5p(v,h);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vh_t(double v, double h, int area)
{	
	double ret;

	if (area==0) area = region_vh(v, h);
	
	switch (area)
	{
		case AREA4: ret = vh4t(v,h);	break;
		case AREA1: ret = ph1t(vh1p(v,h),h);	break;
		case AREA2: ret = ph2t(vh2p(v,h),h);	break;
		case AREA3: ret = rh3t(1/v, h);	break;
		case AREA5: ret = ph5t(vh5p(v,h),h);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vh_u(double v, double h, int area)
{	
	double ret, t, p;

	if (area==0) area = region_vh(v, h);
	
	switch (area)
	{
		case AREA4:
			t = vh4t(v,h);
			ret = tx4u(t, vh4x(v,h));	break;
		case AREA1:
			p = vh1p(v, h);
			ret = pt1u(p, ph1t(p, h));	break;
		case AREA2:
			p = vh2p(v, h);
			ret = pt2u(p, ph2t(p, h));	break;
		case AREA3: ret = rt3u(1/v,rh3t(1/v, h));	break;
		case AREA5:
			p = vh5p(v, h);
			ret = pt5u(p, ph5t(p, h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vh_s(double v, double h, int area)
{	
	double ret, t, p;

	if (area==0) area = region_vh(v, h);
	
	switch (area)
	{
		case AREA4:
			t = vh4t(v,h);
			ret = tx4s(t, vh4x(v,h));	break;
		case AREA1:
			p = vh1p(v, h);
			ret = pt1s(p, ph1t(p, h));	break;
		case AREA2:
			p = vh2p(v, h);
			ret = pt2s(p, ph2t(p, h));	break;
		case AREA3: ret = rt3s(1/v,rh3t(1/v, h));	break;
		case AREA5:
			p = vh5p(v, h);
			ret = pt5s(p, ph5t(p, h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vh_cp(double v, double h, int area)
{	
	double ret, t, p;

	if (area==0) area = region_vh(v, h);
	
	switch (area)
	{
		case AREA4:
			t = vh4t(v,h);
			ret = tx4cp(t, vh4x(v,h));	break;
		case AREA1:
			p = vh1p(v, h);
			ret = pt1cp(p, ph1t(p, h));	break;
		case AREA2:
			p = vh2p(v, h);
			ret = pt2cp(p, ph2t(p, h));	break;
		case AREA3: ret = rt3cp(1/v,rh3t(1/v, h));	break;
		case AREA5:
			p = vh5p(v, h);
			ret = pt5cp(p, ph5t(p, h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vh_cv(double v, double h, int area)
{	
	double ret, t, p;

	if (area==0) area = region_vh(v, h);
	
	switch (area)
	{
		case AREA4:
			t = vh4t(v,h);
			ret = tx4cv(t, vh4x(v,h));	break;
		case AREA1:
			p = vh1p(v, h);
			ret = pt1cv(p, ph1t(p, h));	break;
		case AREA2:
			p = vh2p(v, h);
			ret = pt2cv(p, ph2t(p, h));	break;
		case AREA3: ret = rt3cv(1/v,rh3t(1/v, h));	break;
		case AREA5:
			p = vh5p(v, h);
			ret = pt5cv(p, ph5t(p, h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vh_w(double v, double h, int area)
{	
	double ret, t, p;

	if (area==0) area = region_vh(v, h);
	
	switch (area)
	{
		case AREA4:
			t = vh4t(v,h);
			ret = tx4w(t, vh4x(v,h));	break;
		case AREA1:
			p = vh1p(v, h);
			ret = pt1w(p, ph1t(p, h));	break;
		case AREA2:
			p = vh2p(v, h);
			ret = pt2w(p, ph2t(p, h));	break;
		case AREA3: ret = rt3w(1/v,rh3t(1/v, h));	break;
		case AREA5:
			p = vh5p(v, h);
			ret = pt5w(p, ph5t(p, h));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}


/* Given v, s */
DllExport double WINAPI vs_p(double v, double s, int area)
{	
	double ret;

	if (area==0) if (area==0) area = region_vs(v, s);
	
	switch (area)
	{
		case AREA4: ret = ps(vs4t(v,s));	break;
		case AREA1: ret = vs1p(v,s);	break;
		case AREA2: ret = vs2p(v,s);	break;
		case AREA3: ret = rt3p(1/v,rs3t(1/v, s));	break;
		case AREA5: ret = vs5p(v,s);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vs_t(double v, double s, int area)
{	
	double ret;

	if (area==0) area = region_vs(v, s);
	
	switch (area)
	{
		case AREA4: ret = vs4t(v,s);	break;
		case AREA1: ret = ps1t(vs1p(v,s),s);	break;
		case AREA2: ret = ps2t(vs2p(v,s),s);	break;
		case AREA3: ret = rs3t(1/v, s);			break;
		case AREA5: ret = ps5t(vs5p(v,s),s);	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vs_u(double v, double s, int area)
{	
	double ret, p, t;

	if (area==0) area = region_vs(v, s);
	
	switch (area)
	{
		case AREA4:
			t = vs4t(v,s);
			ret =tx4u(t, vs4x(v,s));			break;
		case AREA1:
			p = vs1p(v,s);
			ret = pt1u(p, ps1t(p,s));	break;
		case AREA2:
			p = vs2p(v,s);
			ret = pt2u(p, ps2t(p,s));	break;
		case AREA3: ret = rt3u(1/v,rs3t(1/v, s));	break;
		case AREA5:
			p = vs5p(v,s);
			ret = pt5u(p, ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vs_h(double v, double s, int area)
{	
	double ret, p, t;

	if (area==0) area = region_vs(v, s);
	
	switch (area)
	{
		case AREA4:
			t = vs4t(v,s);
			ret =tx4h(t, vs4x(v,s));			break;
		case AREA1:
			p = vs1p(v,s);
			ret = pt1h(p, ps1t(p,s));	break;
		case AREA2:
			p = vs2p(v,s);
			ret = pt2h(p, ps2t(p,s));	break;
		case AREA3: ret = rt3h(1/v,rs3t(1/v, s));	break;
		case AREA5:
			p = vs5p(v,s);
			ret = pt5h(p, ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vs_cp(double v, double s, int area)
{
	double ret, p, t;

	if (area==0) area = region_vs(v, s);
	
	switch (area)
	{
		case AREA4:
			t = vs4t(v,s);
			ret =tx4cp(t, vs4x(v,s));			break;
		case AREA1:
			p = vs1p(v,s);
			ret = pt1cp(p, ps1t(p,s));	break;
		case AREA2:
			p = vs2p(v,s);
			ret = pt2cp(p, ps2t(p,s));	break;
		case AREA3: ret = rt3cp(1/v,rs3t(1/v, s));	break;
		case AREA5:
			p = vs5p(v,s);
			ret = pt5cp(p, ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vs_cv(double v, double s, int area)
{
	
	double ret, p, t;

	if (area==0) area = region_vs(v, s);
	
	switch (area)
	{
		case AREA4:
			t = vs4t(v,s);
			ret =tx4cv(t, vs4x(v,s));			break;
		case AREA1:
			p = vs1p(v,s);
			ret = pt1cv(p, ps1t(p,s));	break;
		case AREA2:
			p = vs2p(v,s);
			ret = pt2cv(p, ps2t(p,s));	break;
		case AREA3: ret = rt3cv(1/v,rs3t(1/v, s));	break;
		case AREA5:
			p = vs5p(v,s);
			ret = pt5cv(p, ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vs_w(double v, double s, int area)
{
	
	double ret, p, t;

	if (area==0) area = region_vs(v, s);
	
	switch (area)
	{
		case AREA4:
			t = vs4t(v,s);
			ret =tx4w(t, vs4x(v,s));			break;
		case AREA1:
			p = vs1p(v,s);
			ret = pt1w(p, ps1t(p,s));	break;
		case AREA2:
			p = vs2p(v,s);
			ret = pt2w(p, ps2t(p,s));	break;
		case AREA3: ret = rt3w(1/v,rs3t(1/v, s));	break;
		case AREA5:
			p = vs5p(v,s);
			ret = pt5w(p, ps5t(p,s));	break;
		default: ret = (double)area;
	};
	
	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}


// Saturation line
DllExport double WINAPI ps(double t)
{
	double ps;

	if (t<Tt || t>TC) return T_OUTRANGE;
	ps=t4p(t);

	return ps;
}

DllExport double WINAPI ts(double p)
{
	double ts;

	if (p<Pt || p>PC) return P_OUTRANGE;
	ts=p4t(p);

	return ts;
}

DllExport double WINAPI ph_x(double p, double h, int area)
{
	double ret;

	if (area==0) area = region_ph(p, h);
	
	if (area != AREA4) return AREAX;
	ret = ph4x(p,h);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ps_x(double p, double s, int area)
{
	double ret;

	if (area==0) area = region_ps(p, s);
	
	if (area != AREA4) return AREAX;
	ret = ts4x(p4t(p),s);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI th_x(double t, double h, int area)
{
	double ret;

	if (area==0) area = region_th(t, h);
	
	if (area != AREA4) return AREAX;
	ret = th4x(t,h);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI ts_x(double t, double s, int area)
{
	double ret;

	if (area==0) area = region_ts(t, s);
	
	if (area != AREA4) return AREAX;
	ret = ts4x(t,s);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI hs_x(double h, double s, int area)
{
	double ret;

	if (area==0) area = region_hs(h, s);
	
	if (area != AREA4) return AREAX;
	ret = hs4x(h,s);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tv_x(double t, double v, int area)
{
	double ret;

	if (area==0) area = region_tv(t, v);
	
	if (area != AREA4) return AREAX;
	ret = tv4x(t,v);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vh_x(double v, double h, int area)
{
	double ret;

	if (area==0) area = region_vh(v, h);
	
	if (area != AREA4) return AREAX;
	ret = vh4x(v,h);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vs_x(double v, double s, int area)
{
	double ret;

	if (area==0) area = region_vs(v, s);
	
	if (area != AREA4) return AREAX;
	ret = vs4x(v,s);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}


DllExport double WINAPI tx_v(double t, double x, int area)
{
	double ret;

	if (area==0) area = region_tx(t, x);
	
	if (area != AREA4) return AREAX;
	ret = tx4v(t,x);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tx_h(double t, double x, int area)
{
	double ret;

	if (area==0) area = region_tx(t, x);
	
	if (area != AREA4) return AREAX;
	ret = tx4h(t,x);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tx_s(double t, double x, int area)
{
	double ret;

	if (area==0) area = region_tx(t, x);
	
	if (area != AREA4) return AREAX;
	ret = tx4s(t,x);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI vx_t(double v, double x, int area)
{
	double ret;

	if (area==0) area = region_vx(v, x);
	
	if (area != AREA4) return AREAX;
	ret = vx4t(v,x);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI hx_t(double h, double x, int area)
{
	double ret;

	if (area==0) area = region_hx(h, x);
	
	if (area != AREA4) return AREAX;
	ret = hx4t(h,x);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI sx_t(double s, double x, int area)
{
	double ret;

	if (area==0) area = region_sx(s, x);
	
	if (area != AREA4) return AREAX;
	ret = sx4t(s,x);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tx_u(double t, double x, int area)
{
	double ret;

	if (area==0) area = region_tx(t, x);
	
	if (area != AREA4) return AREAX;
	ret = tx4u(t,x);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tx_cp(double t, double x, int area)
{
	double ret;

	if (area==0) area = region_tx(t, x);
	
	if (area != AREA4) return AREAX;
	ret = tx4cp(t,x);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tx_cv(double t, double x, int area)
{
	double ret;

	if (area==0) area = region_tx(t, x);
	
	if (area != AREA4) return AREAX;
	ret = tx4cv(t,x);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}

DllExport double WINAPI tx_w(double t, double x, int area)
{
	double ret;

	if (area==0) area = region_tx(t, x);
	
	if (area != AREA4) return AREAX;
	ret = tx4w(t,x);

	if (ret>-10 && ret<-1) return (ret+area); 
	return ret;
}


