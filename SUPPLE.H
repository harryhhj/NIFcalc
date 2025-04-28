/* 
   Filename: supple.h
   Supplementary Equations for other Porperties
*/
#include <windows.h>

#define DllImport	__declspec( dllimport )
#define DllExport	__declspec( dllexport )

/* Dynamic Viscosity */
DllExport double WINAPI dyn_visc(double p, double t, double r, int area);
/* Pressure along the Melting and Sublimation Curves of Ordinary Water Substance (September Czech Republic 2011) */
/* required by Thermal Conductivity */
double pm_Ih_t(double t);
double pm_III_t(double t);
double pm_V_t(double t);
double pm_VI_t(double t);
double pm_VII_t(double t);
double pm_t(double t);
double tm_p(double p);
/* Thermal Conductivity */
DllExport double WINAPI thm_cond(double p, double t, double r, int area);
/* Surface Tension */
DllExport double WINAPI surf_tens(double t, int area);
/* Static Dielectric Constant */
DllExport double WINAPI sta_dielec(double p, double t, double r, int area);
/* Refractive Index */
DllExport double WINAPI Ref_Index(double p, double t, double r, double rmd, int area);