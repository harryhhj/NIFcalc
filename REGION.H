/* Filename: region */
#include <windows.h>

#define DllImport	__declspec( dllimport )
#define DllExport	__declspec( dllexport )

#define P_OUTRANGE -10
#define T_OUTRANGE -20
#define V_OUTRANGE -30
#define H_OUTRANGE -40
#define S_OUTRANGE -50
#define X_OUTRANGE -60
#define AREAX -90 // invalid range
#define AREA0  0    // unidentified range
#define UNKNOWN 0
#define AREA1  1    // Compressed Water
#define WATER  1
#define AREA2  2    // Superheated Steam
#define STEAM  2
#define AREA3  3    // Supercritical area
#define CRITICAL 3
#define AREA4  4    // Saturation area (include 8,9)
#define SATURATE 4
#define AREA5  5    // Vast Superheated Steam
#define ULTRA  5
#define SEMI_SATURATE  6   // less saturation water along saturaed water
#define SEMI_SUPERHEAT 7   // less superheated steam along saturated steam
#define SATUR_WATER  8   // less saturation water along saturaed water
#define SATUR_VAPOR 9   // less superheated steam along saturated steam

DllExport int WINAPI region_pt(double p, double t);
DllExport int WINAPI region_ph(double p, double h);
DllExport int WINAPI region_ps(double p, double s);
DllExport int WINAPI region_th(double t, double h);
DllExport int WINAPI region_ts(double t, double s);
DllExport int WINAPI region_hs(double h, double s);
DllExport int WINAPI region_tv(double t, double v);
DllExport int WINAPI region_pv(double p, double v);
DllExport int WINAPI region_vh(double v, double h);
DllExport int WINAPI region_vs(double v, double s);
DllExport int WINAPI region_px(double p, double x);
DllExport int WINAPI region_tx(double t, double x);
DllExport int WINAPI region_vx(double v, double x);
DllExport int WINAPI region_hx(double h, double x);
DllExport int WINAPI region_sx(double s, double x);

