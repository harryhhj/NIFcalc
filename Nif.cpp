/*  filename: nif.c                           */
/*  function: new formula for calculating W&S */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <conio.h>
#include "nif.h"
#include "tools.h"
#include "area5.h"
#include "area4.h"
#include "area3.h"
#include "area2.h"
#include "area1.h"
#include "areas.h"
//#include "message.h"
#include "supple.h"
#include "generic.h"
#include "region.h"




void main()
{
	int area;
	int i,j;
	double p,t,r,h,s, v,x,other, out;
	double a[3][3], b[3]={0,1,2};
	//ofstream fp("e:\hhj\nif\data\demo.dat");
		

/* area 1 forward
printf("area 1.\n");
printf("v %16.12lf, %16.12lf, %16.12lf\n", pt1v(3, 300.0), pt1v(80, 300), pt1v(3, 500.0));
printf("u %16.12lf, %16.12lf, %16.12lf\n", pt1u(3, 300.0), pt1u(80, 300), pt1u(3, 500.0));
printf("s %16.12lf, %16.12lf, %16.12lf\n", pt1s(3, 300.0), pt1s(80, 300), pt1s(3, 500.0));
printf("h %16.12lf, %16.12lf, %16.12lf\n", pt1h(3, 300.0), pt1h(80, 300), pt1h(3, 500.0));
printf("cp %16.12lf, %16.12lf, %16.12lf\n", pt1cp(3, 300.0), pt1cp(80, 300), pt1cp(3, 500.0));
printf("cv %16.12lf, %16.12lf, %16.12lf\n", pt1cv(3, 300.0), pt1cv(80, 300), pt1cv(3, 500.0));
printf("w %16.12lf, %16.12lf, %16.12lf\n", pt1w(3, 300.0), pt1w(80, 300), pt1w(3, 500.0));
*/
/* area 1 backward 
printf("t %16.12lf, %16.12lf, %16.12lf\n", ph1t(3, 500.0), ph1t(80, 500), ph1t(80, 1500.0));
printf("t %16.12lf, %16.12lf, %16.12lf\n", ps1t(3, 0.50), ps1t(80, 0.500), ps1t(80, 3.0));
*/
/* area 1 others backward 
t = 172+T0;
s = 2.0;	
p = ts1p(t,s);
printf("t=%16.12lf, s=%16.12lf, p=%16.12lf\n", t,s,p);
printf("p= %16.12lf, t=%16.12lf, s=%16.12lf\n", p,t,pt1s(p,t));

t = 285+T0;
h = 1300.41;	
p = th1p(t,h);
printf("t=%16.12lf, h=%16.12lf, p=%16.12lf\n", t,h,p);
printf("p= %16.12lf, t=%16.12lf, h=%16.12lf\n", p,t,pt1h(p,t));
*//*  
t = 200+T0;
s = 2.3254;	
h = ts1h(t,s);
printf("t=%16.12lf, s=%16.12lf, h=%16.12lf\n", t,s,h);
printf("t=%16.12lf, h=%16.12lf, s=%16.12lf\n", t,h,th1s(t,h));
*//*
h = 1842.65;
s = 4.0333;	
t = hs1t(h,s);
printf("h=%16.12lf, s=%16.12lf, t=%16.12lf\n", h,s,t);
printf("t=%16.12lf, s=%16.12lf, h=%16.12lf\n", t,s,ts1h(t,s));
p = hs1p(h,s);
printf("h=%16.12lf, s=%16.12lf, p=%16.12lf\n", h,s,p);
*//*
p = 50;
v = 0.001115;	
t = pv1t(p,v);
printf("p=%16.12lf, v=%16.12lf, t=%16.12lf\n", p,v,t);
printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", p,t,pt1v(p,t));
*/
/* 
t=25+T0;
v=0.001003;
p=tv1p(t,v,0);
printf("t=%16.12lf, v=%16.12lf, p=%16.12lf, v=%16.12lf\n", t,v,p,pt1v(p,t));
*/
/* 
p=3;
h=115.331;
v=ph1v(p,h);
printf("p=%16.12lf, h=%16.12lf, v=%16.12lf, p=%16.12lf\n", p,h,v,vh1p(v,h,0));
*//* 
p=80;
s=0.36856385;
v=ps1v(p,s);
printf("p=%16.12lf, s=%16.12lf, v=%16.12lf, p=%16.12lf\n", p,s,v,vs1p(v,s));

	t=Tt;
printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", p,t,v,pt1v(p,ps(t)));*/



/* area 2 generic for f
printf("area 2.\n");
printf("v %16.12lf, %16.12lf, %16.12lf\n", pt2v(0.0035, 300.0), pt2v(0.0035, 700.0), pt2v(30, 700.0));
printf("u %16.12lf, %16.12lf, %16.12lf\n", pt2u(0.0035, 300.0), pt2u(0.0035, 700.0), pt2u(30, 700.0));
printf("s %16.12lf, %16.12lf, %16.12lf\n", pt2s(0.0035, 300.0), pt2s(0.0035, 700.0), pt2s(30, 700.0));
printf("h %16.12lf, %16.12lf, %16.12lf\n", pt2h(1, 450.0), pt2h(0.0035, 700.0), pt2h(30, 700.0));
printf("cp %16.12lf, %16.12lf, %16.12lf\n", pt2cp(0.0035, 300.0), pt2cp(0.0035, 700.0), pt2cp(30, 700.0));
printf("cv %16.12lf, %16.12lf, %16.12lf\n", pt2cv(0.0035, 300.0), pt2cv(0.0035, 700.0), pt2cv(30, 700.0));
printf("w %16.12lf, %16.12lf, %16.12lf\n", pt2w(0.0035, 300.0), pt2w(0.0035, 700.0), pt2w(30, 700.0));
*/
/* area 2  generic for s 
printf("area 2s.\n");
printf("v %16.12lf, %16.12lf, %16.12lf\n", pt2v(1, 450.0), pt2v(1, 440.0), pt2v(1.5, 450.0));
printf("u %16.12lf, %16.12lf, %16.12lf\n", pt2u(1, 450.0), pt2u(1, 440.0), pt2u(1.5, 450.0));
printf("s %16.12lf, %16.12lf, %16.12lf\n", pt2s(1, 450.0), pt2s(1, 440.0), pt2s(1.5, 450.0));
printf("h %16.12lf, %16.12lf, %16.12lf\n", pt2h(1.0, 450.0), pt2h(1, 440.0), pt2h(1.5, 450.0));
printf("cp %16.12lf, %16.12lf, %16.12lf\n", pt2cp(1, 450.0), pt2cp(1, 440.0), pt2cp(1.5, 450.0));
printf("cv %16.12lf, %16.12lf, %16.12lf\n", pt2cv(1, 450.0), pt2cv(1, 440.0), pt2cv(1.5, 450.0));
printf("w %16.12lf, %16.12lf, %16.12lf\n", pt2w(1, 450.0), pt2w(1, 440.0), pt2w(1.5, 450.0));
*/
/* area 2s 
printf("area 2s.\n");
printf("v %16.12lf, %16.12lf, %16.12lf\n", pt2sv(1, 450.0), pt2sv(1, 440.0), pt2sv(1.5, 450.0));
printf("u %16.12lf, %16.12lf, %16.12lf\n", pt2su(1, 450.0), pt2su(1, 440.0), pt2su(1.5, 450.0));
printf("s %16.12lf, %16.12lf, %16.12lf\n", pt2ss(1, 450.0), pt2ss(1, 440.0), pt2ss(1.5, 450.0));
printf("h %16.12lf, %16.12lf, %16.12lf\n", pt2sh(1.0, 450.0), pt2sh(1, 440.0), pt2sh(1.5, 450.0));
printf("cp %16.12lf, %16.12lf, %16.12lf\n", pt2scp(1, 450.0), pt2scp(1, 440.0), pt2scp(1.5, 450.0));
printf("cv %16.12lf, %16.12lf, %16.12lf\n", pt2scv(1, 450.0), pt2scv(1, 440.0), pt2scv(1.5, 450.0));
printf("w %16.12lf, %16.12lf, %16.12lf\n", pt2sw(1, 450.0), pt2sw(1, 440.0), pt2sw(1.5, 450.0));
*/
/* area 2 backward equation p,h -> t 
printf("area 2s backward equation p,h -> t.\n");
printf("h2a t %16.12lf, %16.12lf, %16.12lf\n", ph2at(0.001, 3000.0), ph2at(3, 3000.0), ph2at(3, 4000.0));
printf("h2b t %16.12lf, %16.12lf, %16.12lf\n", ph2bt(5, 3500.0), ph2bt(5, 4000.0), ph2bt(25, 3500.0));
printf("h2c t %16.12lf, %16.12lf, %16.12lf\n", ph2ct(40, 2700.0), ph2ct(60, 2700.0), ph2ct(60, 3200.0));
*/
/* area 2 backward equation p,s -> t 
printf("area 2s backward equation p,s -> t.\n");
printf("s2a t %16.12lf, %16.12lf, %16.12lf\n", ps2at(0.1, 7.5), ps2at(0.1, 8), ps2at(2.5, 8));
printf("s2b t %16.12lf, %16.12lf, %16.12lf\n", ps2bt(8, 6), ps2bt(8, 7.5), ps2bt(90, 6));
printf("s2c t %16.12lf, %16.12lf, %16.12lf\n", ps2ct(20, 5.75), ps2ct(80, 5.25), ps2ct(80, 5.75));
*/
/* area 2 backward equation p,h -> t 
printf("area 2s backward equation p,h -> t.\n");
printf("h2a t %16.12lf, %16.12lf, %16.12lf\n", ph2t(0.001, 3000.0), ph2t(3, 3000.0), ph2t(3, 4000.0));
printf("h2b t %16.12lf, %16.12lf, %16.12lf\n", ph2t(5, 3500.0), ph2t(5, 4000.0), ph2t(25, 3500.0));
printf("h2c t %16.12lf, %16.12lf, %16.12lf\n", ph2t(40, 2700.0), ph2t(60, 2700.0), ph2t(60, 3200.0));
*/
/*
printf("area 2 general backward equation p,s -> t.\n");
printf("s2a t %16.12lf, %16.12lf, %16.12lf\n", ps2t(0.1, 7.5), ps2t(0.1, 8), ps2t(2.5, 8));
printf("s2b t %16.12lf, %16.12lf, %16.12lf\n", ps2t(8, 6), ps2t(8, 7.5), ps2t(90, 6));
printf("s2c t %16.12lf, %16.12lf, %16.12lf\n", ps2t(20, 5.75), ps2t(80, 5.25), ps2t(80, 5.75));
*/
/* area 2 	B2bc-equation. 
printf("s2a t %16.12lf, %16.12lf, %16.12lf\n", pB2bch(100.0), hB2bcp(pB2bch(100.0)),pB2bch(100.0));
*/
/* t s -> p ,h
t = 600+T0;
s = 5.2806;	
p = ts2p(t,s);
h = ts2h(t,s);
printf("t=%16.12lf, s=%16.12lf, p=%16.12lf\n", t,s,p);
printf("t=%16.12lf, p=%16.12lf, s=%16.12lf\n", t, p, pt2s(p,t));
printf("t=%16.12lf, s=%16.12lf, h=%16.12lf\n", t,s,h);
*/
/* t h -> p,s
t = 200+T0;
h = 2800;	
p = th2p(t,h);
printf("t=%16.12lf, s=%16.12lf, h=%16.12lf\n", t,h,p);
printf("t=%16.12lf, p=%16.12lf, s=%16.12lf\n", t, p, pt2h(p,t));
*/
/*h = 3462.59;
s = 7.3251;	
t = hs2t(h,s);
printf("h=%16.12lf, s=%16.12lf, t=%16.12lf\n", h,s,t);
printf("t=%16.12lf, s=%16.12lf, h=%16.12lf\n", t, s, ts2h(t,s));
*/
/*h = 3705.96;
s = 9.42;	
t = hs2t(h,s);
p = hs2p(h,s);
printf("h=%16.12lf, s=%16.12lf, p=%16.12lf\n", h,s,p);
printf("p=%16.12lf, s=%16.12lf, h=%16.12lf\n", p,s,ps2h(p,s));
printf("t=%16.12lf, p=%16.12lf, h=%16.12lf\n", t, p, pt2h(p,t));
*/
/* 
p =15;
v=0.001143;
t=pv2t(p,v);
printf("p=%16.12lf, v=%16.12lf, t=%16.12lf, v=%16.12lf\n", p,v,t,pt2v(p,t));
*/
/* 
t=600+T0;
v=201.488;
p=tv2p(t,v,0);
printf("t=%16.12lf, v=%16.12lf, p=%16.12lf, v=%16.12lf\n", t,v,p,pt2v(p,t));
*/
/* 
p=0.0035;
h=2549.91145;
v=ph2v(p,h);
printf("p=%16.12lf, h=%16.12lf, v=%16.12lf, p=%16.12lf\n", p,h,v,vh2p(v,h,0));
*//* 
p=30;
s=5.1754;
v=ps2v(p,s);
printf("p=%16.12lf, s=%16.12lf, v=%16.12lf, p=%16.12lf\n", p,s,v,vs2p(v,s));
*/


/* areas.cpp | area 2,3 boundary 
printf("area 2,3 boundary t <--> p.\n");
printf("boundary 2,3 %16.12lf, %16.12lf, %16.12lf\n", tB23p(0.62315e3), pB23t(tB23p(0.62315e3)), 0.62315e3);
printf("boundary 2,3 %16.12lf, %16.12lf, %16.12lf\n", pB23t(60), pB23t(100), pB23t(80));

j=0;
  for (i=TB2;i<=TB23;i+=5)
  {
t=(double)i;
p=tB23p(t);
printf("%16.12lf, %16.12lf, %16.12lf\n", pt2v(p,t), pt2h(p,t),pt2s(p,t));
if (j++>20){j=0;getch();}
}*/





/* area 5 
printf("v %16.12lf, %16.12lf, %16.12lf\n", pt5v(0.5, 1500.0), pt5v(8.0,1500.0), pt5v(8, 2000.0));
printf("u %16.12lf, %16.12lf, %16.12lf\n", pt5u(0.5, 1500.0), pt5u(8.0,1500.0), pt5u(8, 2000.0));
printf("s %16.12lf, %16.12lf, %16.12lf\n", pt5s(0.5, 1500.0), pt5s(8.0,1500.0), pt5s(8, 2000.0));
printf("h %16.12lf, %16.12lf, %16.12lf\n", pt5h(0.5, 1500.0), pt5h(8.0,1500.0), pt5h(8, 2000.0));
printf("cp %16.12lf, %16.12lf, %16.12lf\n", pt5cp(0.5, 1500.0), pt5cp(8.0,1500.0), pt5cp(8, 2000.0));
printf("cv %16.12lf, %16.12lf, %16.12lf\n", pt5cv(0.5, 1500.0), pt5cv(8.0,1500.0), pt5cv(8, 2000.0));
printf("w %16.12lf, %16.12lf, %16.12lf\n", pt5w(0.5, 1500.0), pt5w(8.0,1500.0), pt5w(8, 2000.0));
*/
/*p=0.001;
h=4398.42;
t=ph5t(p, h);
printf("p=%16.12lf, h=%16.12lf, t=%16.12lf\n", p,h,t);
printf("p=%16.12lf, t=%16.12lf, h=%16.12lf\n", p,t,pt5h(p,t));
*/
/*p=10;
s=8.2123;
t=ps5t(p, s);
printf("p=%16.12lf, s=%16.12lf, t=%16.12lf\n", p,s,t);
printf("p=%16.12lf, t=%16.12lf, s=%16.12lf\n", p,t,pt5s(p,t));
*/
/*t=1000+T0;
s=8.5935;
p=ts5p(t, s);
printf("t=%16.12lf, s=%16.12lf, p=%16.12lf\n", t,s,p);
printf("p=%16.12lf, t=%16.12lf, s=%16.12lf\n", p,t,pt5s(p,t));
*/
/*t=900+T0;
h=4361.92;
p=th5p(t, h);
printf("t=%16.12lf, h=%16.12lf, p=%16.12lf\n", t,h,p);
printf("p=%16.12lf, t=%16.12lf, h=%16.12lf\n", p,t,pt5h(p,t));
*/
/*
t=800+T0;
s=9.2479;
h=ts5h(t, s);
p=ts5p(t,s);
printf("t=%16.12lf, s=%16.12lf, h=%16.12lf\n", t,s,h);
printf("p=%16.12lf, t=%16.12lf, s=%16.12lf\n", p,t,pt5s(p,t));
*/
/*h=5600.52;
s=8.5131;
p=hs5p(h,s);
t=hs5t(h,s);
printf("h=%16.12lf, s=%16.12lf, p=%16.12lf, t=%16.12lf\n", h,s,p,t);
printf("h=%16.12lf, s=%16.12lf\n", pt5h(p,t), pt5s(p,t));
*/
/*
p=10;
v=0.0984;
t=pv5t(p,v);
printf("p=%16.12lf, v=%16.12lf, t=%16.12lf\n", p,v,t);
printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", p,t,pt5v(p,t));
*/
/*
t=1000+T0;
v=0.1171;
p=tv5p(t,v,0);
printf("t=%16.12lf, v=%16.12lf, p=%16.12lf\n", t,v,p);
printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", p,t,pt5v(p,t));
*/
/* 
p=8;
h=6583.8029;
v=ph5v(p,h);
printf("p=%16.12lf, h=%16.12lf, v=%16.12lf, p=%16.12lf\n", p,h,v,vh5p(v,h,0));
*//* 
p=8;
s=9.15671;
v=ps5v(p,s);
printf("p=%16.12lf, s=%16.12lf, v=%16.12lf, p=%16.12lf\n", p,s,v,vs5p(v,s));
*/



/* area 3 
printf("p %16.12lf, %16.12lf, %16.12lf\n", rt3p(500, 650), rt3p(200,650), rt3p(500,750));
printf("u %16.12lf, %16.12lf, %16.12lf\n", rt3u(500, 650), rt3u(200,650), rt3u(500,750));
printf("s %16.12lf, %16.12lf, %16.12lf\n", rt3s(500, 650), rt3s(200,650), rt3s(500,750));
printf("h %16.12lf, %16.12lf, %16.12lf\n", rt3h(500, 650), rt3h(200,650), rt3h(500,750));
printf("cp %16.12lf, %16.12lf, %16.12lf\n", rt3cp(500, 650), rt3cp(200,650), rt3cp(500,750));
printf("cv %16.12lf, %16.12lf, %16.12lf\n", rt3cv(500, 650), rt3cv(200,650), rt3cv(500,750));
printf("w %16.12lf, %16.12lf, %16.12lf\n", rt3w(500, 650), rt3w(200,650), rt3w(500,750));
*/
/*p=78.3096;
t=750;
v = 1/pt3r(p,t);
printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", p,t,v);
printf("r=%16.12lf, t=%16.12lf, p=%16.12lf\n", 1/v,t,rt3p(1/v,t));
*/
/*t=750;
s=4.4697;
v = 1/ts3r(t,s);
printf("t=%16.12lf, s=%16.12lf, v=%16.12lf\n", t,s,v);
printf("r=%16.12lf, t=%16.12lf, s=%16.12lf\n", 1/v,t,rt3s(1/v,t));
*/
/*t=650;
h=2375.124;
v = 1/th3r(t,h);
printf("t=%16.12lf, h=%16.12lf, v=%16.12lf\n", t,h,v);
printf("r=%16.12lf, t=%16.12lf, h=%16.12lf\n", 1/v,t,rt3h(1/v,t));
*/
/*p=25.5837;
t=650;
h = pt3h(p,t);
printf("p=%16.12lf, t=%16.12lf, h=%16.12lf\n", p,t,h);
*/
/*p=25.5837;
t=650;
s = pt3s(p,t);
printf("p=%16.12lf, t=%16.12lf, s=%16.12lf\n", p,t,s);
*/
/*h=2258.69;
s=4.4697;
t = hs3t(h, s);
printf("h=%16.12lf, s=%16.12lf, t=%16.12lf\n", h,s,t);
*/
/*h=1863.43;
s=4.05427;
r = hs3r(h, s);
printf("h=%16.12lf, s=%16.12lf, r=%16.12lf\n", h,s,r);
*/
/*r=500;
s=4.46971906;
t=rs3t(r,s);
printf("r=%16.12lf, s=%16.12lf, t=%16.12lf\n", r,s,t);
*/
/*r=250;
h=2240.0;
t=rh3t(r,h);
printf("r=%16.12lf, h=%16.12lf, t=%16.12lf\n", r,h,t);
*/





/* area 4
printf("area 4.\n");
printf("t %16.12lf, %16.12lf, %16.12lf\n", p4t(0.1), p4t(1), p4t(10));
printf("P %16.12lf, %16.12lf, %16.12lf\n", t4p(p4t(0.1)), t4p(p4t(1)), t4p(p4t(10)));
*/
/*
t=357.0+T0;
printf("t=%16.12lf, vs1=%16.12lf, vs2=%16.12lf\n", t,t4vs1(t),t4vs2(t));
*/
/*t=360+T0;
printf("t=%16.12lf, hs1=%16.12lf, hs2=%16.12lf\n", t,t4hs1(t),t4hs2(t));
*/
/*t=320+T0;
printf("t=%16.12lf, ss1=%16.12lf, ss2=%16.12lf\n", t,t4ss1(t),t4ss2(t));
*/
/*
t=360+T0;
printf("t=%16.12lf, cps1=%16.12lf, cps2=%16.12lf\n", t,t4cps1(t),t4cps2(t));
*//*
t=360+T0;
printf("t=%16.12lf, ws1=%16.12lf, ws2=%16.12lf\n", t,t4ws1(t),t4ws2(t));
*//*
t=360+T0;
printf("t=%16.12lf, vs1=%16.12lf, vs2=%16.12lf\n", i,tx4v(i+T0,0.9),tx4v(i+T0,0.8));
for (i=350;i<360;i++)
printf("t=%16.12lf, vs1=%16.12lf, vs2=%16.12lf\n", (double)i,t4vs1((double)i+T0),t4vs2((double)i+T0));
*//*v=0.003;
x=0.99;
t=vx4t(v,x);
printf("v=%16.12lf, x=%16.12lf, t=%16.12lf\n", v,x,t);
*/
/*h=2700;
x=0.95;
t=hx4t(h,x);
printf("h=%16.12lf, x=%16.12lf, t=%16.12lf\n", h,x,t);
*/
/*s=7;
x=0.7;
t=sx4t(s,x);
printf("s=%16.12lf, x=%16.12lf, t=%16.12lf\n", s,x,t);
*/
/*v=0.005;
h=2000;
t=vh4t(v,h,0);
printf("v=%16.12lf, h=%16.12lf, t=%16.12lf, h=%16.12lf\n", v,h,t,tv4h(t,v,0));
printf("v=%16.12lf, h=%16.12lf, x=%16.12lf\n", v,s,vh4x(v,h,0));
*/
/*v=0.002;
s=4.0;
t=vs4t(v,s);
printf("v=%16.12lf, s=%16.12lf, t=%16.12lf, s=%16.12lf\n", v,s,t,tv4s(t,v,0));
printf("v=%16.12lf, s=%16.12lf, x=%16.12lf\n", v,s,vs4x(v,s));
*//*h=2500;
s=7;
v=hs4v(h,s);
printf("h=%16.12lf, s=%16.12lf, v=%16.12lf, h=%16.12lf\n", h,s,v,vs4h(v,s));
*/





/* CHECK region_AREA*/
/* area 1 */  /* in area1 t,h may be confused by area4 */
/*p=3;
t=300;
v=0.100215168e-2;
h=0.115331273e3;
s=0.392294792;*/
/*p=80;
t=300;
v=0.971180894e-3;
h=0.184142828e3;
s=0.368563852;*/
/*p=3;
t=500;
v=0.120241800e-2;
h=0.975542239e3;
s=0.258041912e1;*/

/* area 2 */
/*p=30;
t=700;
v=0.542946619e-2;
h=0.263149474e4;
s=0.517540298e1;*/
/*p=0.0035;
t=700;
v=0.923015898e2;
h=0.333568375e4;
s=0.101749996e2;*/
/*p=0.0035;
t=300;
v=0.394913866e2;
h=0.254991145e4;
s=0.852238967e1;*/

/* area 3 */
/*p=0.255837018e2;
t=650;
v=0.002;
h=0.186343019e4;
s=0.405427273e1;*/
/*p=0.222930643e2;
t=650;
v=0.005;
h=0.237512401e4;
s=0.485438792e1;*/
/*p=0.783095639e2;
t=750;
v=0.002;
h=0.225868845e4;
s=0.446971906e1;*/

/* area 5 */
/*p=0.5;
t=1500;
v=0.138455354e1;
h=0.52197632e4;
s=0.965408431e1;*/
/*p=8;
t=1500;
v=0.865156616e-1;
h=0.520609634e4;
s=0.836546724e1;
p=34.4738 ;
t=2000;
v=0.115743146;
h=0.658380291e4;
s=5.08698 ;
//i=region_pt(p, t);
//printf("p, t, area=%2d\n", i);
//i= region_ph( p,  h);
//printf("p, h, area=%2d\n", i);
i= region_ps( p,  s);
printf("p, s, area=%2d\n", i);*/
/*i= region_th( t,  h);
printf("t, h, area=%2d\n", i);
i= region_ts( t,  s);
printf("t, s, area=%2d\n", i);
i= region_hs( h,  s);
printf("h, s, area=%2d\n", i);
i= region_tv( t,  v);
printf("t, v, area=%2d\n", i);
i= region_pv( p,  v);
printf("p, v, area=%2d\n", i);
i= region_vh( v,  h);
printf("v, h, area=%2d\n", i);
i= region_vs( v,  s);
printf("v, s, area=%2d\n", i);*/

/*----------------------------------------------------------*/
/* INTEGRAL TEST */
/* area 1 */
/*p=3;
t=300;
v=0.100215168e-2;
h=0.115331273e3;
s=0.392294792;*/
/*p=80;
t=300;
v=0.971180894e-3;
h=0.184142828e3;
s=0.368563852;*/
/*p=3;
t=500;
v=0.120241800e-2;
h=0.975542239e3;
s=0.258041912e1;*/

/* area 2 */
/*p=30;
t=700;
v=0.542946619e-2;
h=0.263149474e4;
s=0.517540298e1;*/
/*p=0.0035;
t=700;
v=0.923015898e2;
h=0.333568375e4;
s=0.101749996e2;*/
/*p=0.0035;
t=300;
v=0.394913866e2;
h=0.254991145e4;
s=0.852238967e1;*/

/* area 3 */
/*p=0.255837018e2;
t=650;
v=0.002;
h=0.186343019e4;
s=0.405427273e1;*/
/*p=0.222930643e2;
t=650;
v=0.005;
h=0.237512401e4;
s=0.485438792e1;*/
/*p=0.783095639e2;
t=750;
v=0.002;
h=0.225868845e4;
s=0.446971906e1;*/

/* area 5 */
/*p=0.5;
t=1500;
v=0.138455354e1;
h=0.52197632e4;
s=0.965408431e1;*/
/*p=8;
t=1500;
v=0.865156616e-1;
h=0.520609634e4;
s=0.836546724e1;*/
/*p=8;
t=2000;
v=0.115743146;
h=0.658380291e4;
s=0.915671044e1;
i=region_pt(p, t);*/

/*


printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", p,t,pt_v(p,t,0));
printf("h=%16.12lf, u=%16.12lf, s=%16.12lf\n", pt_h(p,t,0),pt_u(p,t,0),pt_s(p,t,0));
printf("cp=%16.12lf, cv=%16.12lf, w=%16.12lf\n", pt_cp(p,t,0),pt_cv(p,t,0),pt_w(p,t,0));
printf(" Given p, t; area=%2d \n", area);
 
printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", p,ph_t(p,h,0),ph_v(p,h,0));
printf("h=%16.12lf, u=%16.12lf, s=%16.12lf\n", h,ph_u(p,h,0),ph_s(p,h,0));
printf("cp=%16.12lf, cv=%16.12lf, w=%16.12lf\n", ph_cp(p,h,0),ph_cv(p,h,0),ph_w(p,h,0));
printf(" Given p, h; area=%2d \n", area);

printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", p,ps_t(p,s,0),ps_v(p,s,0));
printf("h=%16.12lf, u=%16.12lf, s=%16.12lf\n", ps_h(p,s,0),ps_u(p,s,0),s);
printf("cp=%16.12lf, cv=%16.12lf, w=%16.12lf\n", ps_cp(p,s,0),ps_cv(p,s,0),ps_w(p,s,0));
printf(" Given p, s; area=%2d \n", area);

printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", th_p(t,h,0),t,th_v(t,h,0));
printf("h=%16.12lf, u=%16.12lf, s=%16.12lf\n", h,th_u(t,h,0),th_s(t,h,0));
printf("cp=%16.12lf, cv=%16.12lf, w=%16.12lf\n", th_cp(t,h,0),th_cv(t,h,0),th_w(t,h,0));
printf(" Given t, h; area=%2d \n", area);

printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", ts_p(t,s,0),t,ts_v(t,s,0));
printf("h=%16.12lf, u=%16.12lf, s=%16.12lf\n", ts_h(t,s,0),ts_u(t,s,0),s);
printf("cp=%16.12lf, cv=%16.12lf, w=%16.12lf\n", ts_cp(t,s,0),ts_cv(t,s,0),ts_w(t,s,0));
printf(" Given t, s; area=%2d \n", area);

getch();
printf(" \n");

printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", hs_p(h,s,0),hs_t(h,s,0),hs_v(h,s,0));
printf("h=%16.12lf, u=%16.12lf, s=%16.12lf\n", h,hs_u(h,s,0),s);
printf("cp=%16.12lf, cv=%16.12lf, w=%16.12lf\n", hs_cp(h,s,0),hs_cv(h,s,0),hs_w(h,s,0));
printf(" Given h, s; area=%2d \n", area);

printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", tv_p(t,v,0),t,v);
printf("h=%16.12lf, u=%16.12lf, s=%16.12lf\n", tv_h(t,v,0),tv_u(t,v,0),tv_s(t,v,0));
printf("cp=%16.12lf, cv=%16.12lf, w=%16.12lf\n", tv_cp(t,v,0),tv_cv(t,v,0),tv_w(t,v,0));
printf(" Given t, v; area=%2d \n", area);

printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", p,pv_t(p,v,0),v);
printf("h=%16.12lf, u=%16.12lf, s=%16.12lf\n", pv_h(p,v,0),pv_u(p,v,0),pv_s(p,v,0));
printf("cp=%16.12lf, cv=%16.12lf, w=%16.12lf\n", pv_cp(p,v,0),pv_cv(p,v,0),pv_w(p,v,0));
printf(" Given p, v; area=%2d \n", area);

printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", vh_p(v,h,0),vh_t(v,h,0),v);
printf("h=%16.12lf, u=%16.12lf, s=%16.12lf\n", h,vh_u(v,h,0),vh_s(v,h,0));
printf("cp=%16.12lf, cv=%16.12lf, w=%16.12lf\n", vh_cp(v,h,0),vh_cv(v,h,0),vh_w(v,h,0));
printf(" Given v, h; area=%2d \n", area);

printf("p=%16.12lf, t=%16.12lf, v=%16.12lf\n", vs_p(v,s,0),vs_t(v,s,0),v);
printf("h=%16.12lf, u=%16.12lf, s=%16.12lf\n", vs_h(v,s,0),vs_u(v,s,0),s);
printf("cp=%16.12lf, cv=%16.12lf, w=%16.12lf\n", vs_cp(v,s,0),vs_cv(v,s,0),vs_w(v,s,0));
printf(" Given v, s; area=%2d \n", area);
*/

/*for (i=0;i<TB2;i+=10)
{
t=(double)i+T0;
p=ps(t);
printf("v=%16.12lf, h=%16.12lf, s=%16.12lf\n", pt1v(p,t),pt1h(p,t),pt1s(p,t));
if (j++>20){j=0;getch();}
}*/
/*s=6.0405;
for (i=300;i<=800;i+=50)
{
t=double(i)+T0;
p=ts2p(t,s);
printf("t=%16.12lf, p=%16.12lf, h=%16.12lf, v=%16.12lf\n", t,p,pt2h(p,t),pt2v(p,t));
}*/

/*for (i=5;i<=9;i++)
{
s=double(i);
t=sx4t(s,1);
printf("s=%16.12lf, t=%16.12lf\n", s,t);
}
j=0;
for (i=350;i<=374;i+=1)
{
t=double(i)+T0;
s=tx4v(t,0);
printf("t=%16.12lf, s=%16.12lf\n", t,s);
if (j++>10) {getch();j=0;}
}

t=TC;
s=tx4s(t,1);
printf("t=%16.12lf, s=%16.12lf\n", t,s);
*/

//p=22.15;
//t=647.5;
//h=2800;
//s= 6.5;
//area= region_ps( p,  s);
//printf("%d\n",area);
//t=ps3v(p,s);
//t = t_3a_ps(p,s);
//p= sB4h(s);
//t=hs4t(h, s);
//printf("%16.12lf\n", t);
//v=ps3v(p,s);
//p = sB4h(s);
//p= p_2c_hs(h,s);
//h=h_B3ab_p(p);
//out=ph3v(p,h);
//out=ph3t(p,h);
//out=p_B34_s( s);
//out=p_B34_h(h);
//out=hs4x(h,s);
//out = rt3h(1/hs3v(h,s), hs3t(h,s));
//out = pt1h(hs1p(h,s),hs1t(h,s));
//out = pt2s(hs2p(h,s),hs2t(h,s));
//out = t_B3wx_p(p);
//out = v_3y_pt(p, t);
//out=pt3v(p,t);
//out = 1/pt3r(p,t);
//printf("%16.12lf\n", out);

t = 773.15;
p = 100;
r = 1.0 / pt_v(p, t,0);
other = thm_cond(p, t,r, UNKNOWN);
//printf("area= %3d\n", area);
printf("r= %16.12lf, %16.12lf, %16.12lf, %16.12lf\n", p,t,r,other);

/* other= surf_tens(t); 
other = surf_tens(t);*/

/* Static Dielectric Constant 
t = 300;
p=10;
other=sta_dielec(p, t);
printf("r= %16.12lf, %16.12lf, %16.12lf, %16.12lf\n", p,t,r,other);
*/

/* Refractive Index 
t = 273.15;
p=0.1;
other = Ref_Index(p, t, 0.2265);
printf("r= %16.12lf, %16.12lf, %16.12lf, %16.12lf\n", p,t,r,other);
*/
// fp << p << "  " << t << endl;

//fp.close;
//getch();
//return 1;
}