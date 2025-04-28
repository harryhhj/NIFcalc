/**********************************************************************/
/*	      Filename: area2.h                                           */
/*		  Function: Define Functions for Region 2 in NIF              */
/**********************************************************************/



// region 2  forward equations
double Gmo2(double Pi, double To);
double Gmo2_p(double Pi, double To);
double Gmo2_pp(double Pi, double To);
double Gmo2_t(double Pi, double To);
double Gmo2_tt(double Pi, double To);
double Gmo2_pt(double Pi, double To);
double Gmt2(double Pi, double To);
double Gmt2_p(double Pi, double To);
double Gmt2_pp(double Pi, double To);
double Gmt2_t(double Pi, double To);
double Gmt2_tt(double Pi, double To);
double Gmt2_pt(double Pi, double To);
/* fundamental equation */
double pt2fv(double p, double t);
double pt2fu(double p, double t);
double pt2fs(double p, double t);
double pt2fh(double p, double t);
double pt2fcp(double p, double t);
double pt2fcv(double p, double t);
double pt2fw(double p, double t);
// Supplementary Equation for Metastable-Vapor Region
double Gmt2s(double Pi, double To);
double Gmt2s_p(double Pi, double To);
double Gmt2s_pp(double Pi, double To);
double Gmt2s_t(double Pi, double To);
double Gmt2s_tt(double Pi, double To);
double Gmt2s_pt(double Pi, double To);
/* specific volume */
double pt2sv(double p, double t);
double pt2su(double p, double t);
double pt2ss(double p, double t);
double pt2sh(double p, double t);
double pt2scp(double p, double t);
double pt2scv(double p, double t);
double pt2sw(double p, double t);
// generic area 2 equation
int metastable(double p, double t);
double pt2v(double p, double t);
double pt2u(double p, double t);
double pt2s(double p, double t);
double pt2h(double p, double t);
double pt2cp(double p, double t);
double pt2cv(double p, double t);
double pt2w(double p, double t);

  // 2a  p,h -> t
double t_2a_ph(double p,double h);
  // 2b  p,h -> t 
double t_2b_ph(double p,double h);
  // 2c  p,h -> t 
double t_2c_ph(double p,double h);

  // 2a  p,s -> t 
double t_2a_ps(double p,double s);
  // 2b  p,s -> t 
double t_2b_ps(double p,double s);
  // 2c  p,s -> t 
double t_2c_ps(double p,double s);

// generic area 2 backward equation
double p_B2bc_h(double h);
double h_B2bc_p(double p);
double ph2t(double p,double h);
double ps2t(double p,double s);
double ps2h(double p, double s);

double ts2p(double t, double s);
double ts2h(double t, double s);

double th2p(double t, double h);
double th2s(double t, double h);

double hs2t(double h, double s);
double h_B2ab_s(double s); /* boundary formula of 2ab */
double p_2a_hs(double h, double s);
double p_2b_hs(double h, double s);
double p_2c_hs(double h, double s);
double hs2p(double h, double s);

double pv2t(double p, double v);
double tv2p(double t, double v);
double ts2v(double t, double s);
double ph2v(double p, double h);
double vh2p(double v, double h);
double ps2v(double p, double s);
double vs2p(double v, double s);

