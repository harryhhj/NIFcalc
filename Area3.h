/**********************************************************************/
/*	     Filename: area3.h                                            */
/*		 Function: NIF area 3 calculations                            */
/**********************************************************************/

double Fi(double rdr,double rdt);
double Fi_d(double rdr,double rdt);
double Fi_t(double rdr,double rdt);
double Fi_dd(double rdr,double rdt);
double Fi_tt(double rdr,double rdt);
double Fi_dt(double rdr,double rdt);
double rt3p(double ro,double  t);
double rt3u(double ro,double  t);
double rt3s(double ro,double  t);
double rt3h(double ro,double  t);
double rt3cv(double ro,double  t);
double rt3cp(double ro,double  t);
double rt3w(double ro,double  t);
// backward function
double pt3v(double p, double t); /* Jul.2005, Santorni, Greece v3(p,t) */

/* pt3r replaced by pt3v(p,t) 2005 release but reserved as accuracy is required */
double pt3r(double p,double  t);
double ts3r(double t, double s);
double th3r(double t, double h);
/*double pt3h(double p, double t);
double pt3s(double p, double t);*/
double ts3h(double t, double s);
double hs3t(double h, double s);
/* double hs3r(double h, double s); obsolete. -> hs3v(h,s) */
double hs3v(double h, double s);
double p_3a_hs(double h, double s); /* Kyoto, Sep. 2004 */
double p_3b_hs(double h, double s); /* Kyoto, Sep. 2004 */
double hs3p(double h, double s);
double rs3t(double r, double s);
double rh3t(double r, double h);
double rh3p(double r, double h);
/*double ph3r(double p, double h); obsolete implementation */
double h_B3ab_p(double p); /*implemented Aug21-2005 */
double v_3a_ph(double p, double h); /*implemented Aug21-2005*/
double v_3b_ph(double p, double h); /*implemented Aug21-2005*/
double ph3v(double p, double h); /*implemented Aug21-2005*/
double t_3a_ph(double p, double h); /*implemented Aug21-2005*/
double t_3b_ph(double p, double h); /*implemented Aug21-2005*/
double ph3t(double p, double h); /*implemented Aug21-2005*/
double rs3p(double r, double s);
/*double ps3r(double p, double s); obsolete */
double v_3a_ps(double p, double s);/* implemented Aug18-2005 */
double v_3b_ps(double p, double s);/* implemented Aug18-2005 */
double ps3v(double p, double s); /* implemented Aug18-2005 */
double t_3a_ps(double p, double s);/* implemented Aug18-2005 */
double t_3b_ps(double p, double s);/* implemented Aug18-2005 */
double ps3t(double p, double s);/* implemented Aug18-2005 */
double rp3t(double r, double p);
double pt3h(double p, double t);
double pt3s(double p, double t);
/* Jul.2005, Santorni, Greece v3(p,t) */
//boundaries
double t_B3ab_p(double p); 
double t_B3cd_p(double p);
double t_B3ef_p(double p); 
double t_B3gh_p(double p); 
double t_B3ij_p(double p);
double t_B3jk_p(double p);
double t_B3mn_p(double p);
double t_B3op_p(double p);
double t_B3qu_p(double p);
double t_B3rx_p(double p);
/* T_B3uv(p) & T_B2wx(P) are very close to critical */
double t_B3uv_p(double p);
double t_B3wx_p(double p);
//formula
double v_3a_pt(double p, double t);
double v_3b_pt(double p, double t);
double v_3c_pt(double p, double t);
double v_3d_pt(double p, double t);
double v_3e_pt(double p, double t);
double v_3f_pt(double p, double t);
double v_3g_pt(double p, double t);
double v_3h_pt(double p, double t);
double v_3i_pt(double p, double t);
double v_3j_pt(double p, double t);
double v_3k_pt(double p, double t);
double v_3l_pt(double p, double t);
double v_3m_pt(double p, double t);
double v_3n_pt(double p, double t);
double v_3o_pt(double p, double t);
double v_3p_pt(double p, double t);
double v_3q_pt(double p, double t);
double v_3r_pt(double p, double t);
double v_3s_pt(double p, double t);
double v_3t_pt(double p, double t);
/* 3u - 3z very close to critical */
double v_3aux_pt(double p, double t); /* very close to critical  */
double v_3u_pt(double p, double t); 
double v_3v_pt(double p, double t);
double v_3w_pt(double p, double t); 
double v_3x_pt(double p, double t);
double v_3y_pt(double p, double t);
double v_3z_pt(double p, double t); 



/* derived ps3h() */
double ps3h(double p, double s);
