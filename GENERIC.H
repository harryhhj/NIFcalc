/*
Filename: generic.cpp
generic function regarding to areas.
*/
#include <windows.h>

#define DllImport	__declspec( dllimport )
#define DllExport	__declspec( dllexport )

/* Given p, t */
DllExport double WINAPI pt_v(double p, double t, int area);
DllExport double WINAPI pt_u(double p, double t, int area);
DllExport double WINAPI pt_h(double p, double t, int area);
DllExport double WINAPI pt_s(double p, double t, int area);
DllExport double WINAPI pt_cp(double p, double t, int area);
DllExport double WINAPI pt_cv(double p, double t, int area);
DllExport double WINAPI pt_w(double p, double t, int area);

/* Given p, h */
DllExport double WINAPI ph_t(double p, double h, int area);
DllExport double WINAPI ph_v(double p, double h, int area);
DllExport double WINAPI ph_u(double p, double h, int area);
DllExport double WINAPI ph_s(double p, double h, int area);
DllExport double WINAPI ph_cp(double p, double h, int area);
DllExport double WINAPI ph_cv(double p, double h, int area);
DllExport double WINAPI ph_w(double p, double h, int area);

/* Given p, s */
DllExport double WINAPI ps_t(double p, double s, int area);
DllExport double WINAPI ps_v(double p, double s, int area);
DllExport double WINAPI ps_u(double p, double s, int area);
DllExport double WINAPI ps_h(double p, double s, int area);
DllExport double WINAPI ps_cp(double p, double s, int area);
DllExport double WINAPI ps_cv(double p, double s, int area);
DllExport double WINAPI ps_w(double p, double s, int area);

/* Given t, h */
DllExport double WINAPI th_p(double t, double h, int area);
DllExport double WINAPI th_v(double t, double h, int area);
DllExport double WINAPI th_u(double t, double h, int area);
DllExport double WINAPI th_s(double t, double h, int area);
DllExport double WINAPI th_cp(double t, double h, int area);
DllExport double WINAPI th_cv(double t, double h, int area);
DllExport double WINAPI th_w(double t, double h, int area);

/* Given t, s */
DllExport double WINAPI ts_p(double t, double s, int area);
DllExport double WINAPI ts_v(double t, double s, int area);
DllExport double WINAPI ts_u(double t, double s, int area);
DllExport double WINAPI ts_h(double t, double s, int area);
DllExport double WINAPI ts_cp(double t, double s, int area);
DllExport double WINAPI ts_cv(double t, double s, int area);
DllExport double WINAPI ts_w(double t, double s, int area);

/* Given h, s */
DllExport double WINAPI hs_p(double h, double s, int area);
DllExport double WINAPI hs_t(double h, double s, int area);
DllExport double WINAPI hs_v(double h, double s, int area);
DllExport double WINAPI hs_u(double h, double s, int area);
DllExport double WINAPI hs_cp(double h, double s, int area);
DllExport double WINAPI hs_cv(double h, double s, int area);
DllExport double WINAPI hs_w(double h, double s, int area);

/* Given t, v */
DllExport double WINAPI tv_p(double t, double v, int area);
DllExport double WINAPI tv_u(double t, double v, int area);
DllExport double WINAPI tv_h(double t, double v, int area);
DllExport double WINAPI tv_s(double t, double v, int area);
DllExport double WINAPI tv_cp(double t, double v, int area);
DllExport double WINAPI tv_cv(double t, double v, int area);
DllExport double WINAPI tv_w(double t, double v, int area);

/* Given p, v */
DllExport double WINAPI pv_t(double p, double v, int area);
DllExport double WINAPI pv_u(double p, double v, int area);
DllExport double WINAPI pv_h(double p, double v, int area);
DllExport double WINAPI pv_s(double p, double v, int area);
DllExport double WINAPI pv_cp(double p, double v, int area);
DllExport double WINAPI pv_cv(double p, double v, int area);
DllExport double WINAPI pv_w(double p, double v, int area);

/* Given v, h */
DllExport double WINAPI vh_t(double v, double h, int area);
DllExport double WINAPI vh_p(double v, double h, int area);
DllExport double WINAPI vh_u(double v, double h, int area);
DllExport double WINAPI vh_s(double v, double h, int area);
DllExport double WINAPI vh_cp(double v, double h, int area);
DllExport double WINAPI vh_cv(double v, double h, int area);
DllExport double WINAPI vh_w(double v, double h, int area);

/* Given v, s */
DllExport double WINAPI vs_t(double v, double s, int area);
DllExport double WINAPI vs_p(double v, double s, int area);
DllExport double WINAPI vs_u(double v, double s, int area);
DllExport double WINAPI vs_h(double v, double s, int area);
DllExport double WINAPI vs_cp(double v, double s, int area);
DllExport double WINAPI vs_cv(double v, double s, int area);
DllExport double WINAPI vs_w(double v, double s, int area);

/* In saturation area */
DllExport double WINAPI ps(double t);
DllExport double WINAPI ts(double p);

DllExport double WINAPI ph_x(double p, double h, int area);
DllExport double WINAPI ps_x(double p, double s, int area);
DllExport double WINAPI th_x(double t, double h, int area);
DllExport double WINAPI ts_x(double t, double s, int area);
DllExport double WINAPI hs_x(double h, double s, int area);
DllExport double WINAPI tv_x(double t, double v, int area);
DllExport double WINAPI vh_x(double v, double h, int area);
DllExport double WINAPI vs_x(double v, double s, int area);

DllExport double WINAPI tx_v(double t, double x, int area);
DllExport double WINAPI tx_h(double t, double x, int area);
DllExport double WINAPI tx_s(double t, double x, int area);
DllExport double WINAPI vx_t(double v, double x, int area);
DllExport double WINAPI hx_t(double h, double x, int area);
DllExport double WINAPI sx_t(double s, double x, int area);
DllExport double WINAPI tx_u(double t, double x, int area);
DllExport double WINAPI tx_cp(double t, double x, int area);
DllExport double WINAPI tx_cv(double t, double x, int area);
DllExport double WINAPI tx_w(double t, double x, int area);

