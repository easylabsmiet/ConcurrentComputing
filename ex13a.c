//  Linear parabolic equation:
//
//  du/dt = (d/dx)(k(x)du/dx) + f(x,t), xa < x < xb, t>0
//
//  u(x,0) = g0(x), u(xa,t) = g1(t), u(xb,t) = g2(t) 
//
//  Explicit scheme
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mycom.h"
#include "mynet.h"
#include "myio.h"

 int np, mp, nl, ier, lp;
 int mp_l, mp_r;
 char pname[MPI_MAX_PROCESSOR_NAME];
 char vname[10] = "ex13a";
 char sname[20];
 MPI_Status status;
 union_t buf;
 double tick, t1, t2, t3;

 FILE *Fi = NULL;
 FILE *Fo = NULL;

 int nx, ntp, ntm, ntv;
 double xa, xb, xk, x0, r0, q0, u0, u1;
 double k1, k2, tau0, tau1, tmax, epst;
 double tv, u10, omg0, omg1, gt;

double k(double x);
double k(double x) {
  if (x<xk) return k1;
  else      return k2;
}

double f(double x, double t);
double f(double x, double t) {
  double s1 = (x-x0) / r0;
  double s2 = omg0 * t;
  return q0*exp(-s1*s1)*(1.0-exp(-s2));
}

double g0(double x);
double g0(double x) {
  return u0;
}

double g1(double t);
double g1(double t) {
  double s1 = omg1 * t;
  return u0 + u10 * (1.0-exp(-s1));
}

double g2(double t);
double g2(double t) {
  return u0;
}

int main(int argc, char *argv[])
{
  int i, j, ii, i1, i2, nc, ncm;
  double hx, hx2, tau, gam, s0, s1, s2, s3, y0m, y0p;
  double *xx, *aa, *bb, *y0, *y1;

  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  sprintf(sname,"%s.p%02d",vname,mp);
  ier = fopen_m(&Fo,sname,"wt");
  if (ier!=0) mpierr("Protocol file not opened",1);

  if (mp==0) {
    sprintf(sname,"%s.d",vname);
    ier = fopen_m(&Fi,sname,"rt");
    if (ier!=0) mpierr("Data file not opened",2);
    fscanf(Fi,"xa=%le\n",&xa);
    fscanf(Fi,"xb=%le\n",&xb);
    fscanf(Fi,"xk=%le\n",&xk);
    fscanf(Fi,"x0=%le\n",&x0);
    fscanf(Fi,"r0=%le\n",&r0);
    fscanf(Fi,"q0=%le\n",&q0);
    fscanf(Fi,"u0=%le\n",&u0);
    fscanf(Fi,"u1=%le\n",&u1);
    fscanf(Fi,"k1=%le\n",&k1);
    fscanf(Fi,"k2=%le\n",&k2);
    fscanf(Fi,"tau0=%le\n",&tau0);
    fscanf(Fi,"tau1=%le\n",&tau1);
    fscanf(Fi,"tmax=%le\n",&tmax);
    fscanf(Fi,"epst=%le\n",&epst);
    fscanf(Fi,"nx=%d\n",&nx);
    fscanf(Fi,"ntp=%d\n",&ntp);
    fscanf(Fi,"ntm=%d\n",&ntm);
    fscanf(Fi,"lp=%d\n",&lp);
    fclose_m(&Fi);
    if (argc>1) sscanf(argv[1],"%d",&nx);
    if (argc>2) sscanf(argv[2],"%d",&ntp);
    if (argc>3) sscanf(argv[3],"%d",&ntm);
  }

  if (np>1) {
    if (mp==0) {
      buf.ddata[0]  = xa;
      buf.ddata[1]  = xb;
      buf.ddata[2]  = xk;
      buf.ddata[3]  = x0;
      buf.ddata[4]  = r0;
      buf.ddata[5]  = q0;
      buf.ddata[6]  = u0;
      buf.ddata[7]  = u1;
      buf.ddata[8]  = k1;
      buf.ddata[9]  = k2;
      buf.ddata[10] = tau0;
      buf.ddata[11] = tau1;
      buf.ddata[12] = tmax;
      buf.ddata[13] = epst;
      buf.idata[28] = nx;
      buf.idata[29] = ntp;
      buf.idata[30] = ntm;
      buf.idata[31] = lp;
    }
    MPI_Bcast(buf.ddata,16,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (mp>0) {
      xa   = buf.ddata[0];
      xb   = buf.ddata[1];
      xk   = buf.ddata[2];
      x0   = buf.ddata[3];
      r0   = buf.ddata[4];
      q0   = buf.ddata[5];
      u0   = buf.ddata[6];
      u1   = buf.ddata[7];
      k1   = buf.ddata[8];
      k2   = buf.ddata[9];
      tau0 = buf.ddata[10];
      tau1 = buf.ddata[11];
      tmax = buf.ddata[12];
      epst = buf.ddata[13];
      nx   = buf.idata[28];
      ntp  = buf.idata[29];
      ntm  = buf.idata[30];
      lp   = buf.idata[31];
    }
  }

  fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",
    np,mp,pname,tick);

  fprintf(Fo,"xa=%le xb=%le xk=%le x0=%le r0=%le\n",xa,xb,xk,x0,r0);
  fprintf(Fo,"q0=%le u0=%le u1=%le k1=%le k2=%le\n",q0,u0,u1,k1,k2);
  fprintf(Fo,"tau0=%le tau1=%le tmax=%le epst=%le\n",tau0,tau1,tmax,epst);
  fprintf(Fo,"nx=%d ntp=%d ntm=%d lp=%d\n",nx,ntp,ntm,lp);

  t1 = MPI_Wtime();

  u10 = u1 - u0; omg0 = 1.0 / tau0; omg1 = 1.0 / tau1;
  hx = (xb-xa)/nx; hx2 = hx * hx;
  tau = 0.5 * hx2 / dmax(k1,k2); tau = dmin(tau,1.0/q0);
  gam = tau / hx2;
  s0 = dmin(tmax/tau,1000000000.0); ntm = imin(ntm,(int)s0);

  fprintf(Fo,"u10=%le omg0=%le omg1=%le\n",u10,omg0,omg1);
  fprintf(Fo,"hx=%le tau=%le ntm=%d\n",hx,tau,ntm);
  if (mp == 0) fprintf(stderr,"nx=%d hx=%le tau=%le ntm=%d\n",nx,hx,tau,ntm);

  if (mp ==    0) mp_l = -1; else mp_l = mp - 1;
  if (mp == np-1) mp_r = -1; else mp_r = mp + 1;

  MyRange(np,mp,0,nx,&i1,&i2,&nc);
  ncm = nc-1;

  fprintf(Fo,"i1=%d i2=%d nc=%d\n",i1,i2,nc);

  xx = (double*)(malloc(sizeof(double)*nc));
  y0 = (double*)(malloc(sizeof(double)*nc));
  y1 = (double*)(malloc(sizeof(double)*nc));

  aa = (double*)(malloc(sizeof(double)*nc));
  bb = (double*)(malloc(sizeof(double)*nc));

  for (i=0; i<nc; i++) xx[i] = xa + hx * (i1 + i); // grid

  for (i=0; i<nc; i++) {
    ii = i1 + i;

    if ((ii==0) || (ii==nx)) {
      aa[i] = 0.0; bb[i] = 0.0;
    }
    else {
      s0 = k(xx[i]); s1 = k(xx[i]-hx); s2 = k(xx[i]+hx);
      aa[i] = gam * 2.0 * s0 * s1 / (s0 + s1);
      bb[i] = gam * 2.0 * s0 * s2 / (s0 + s2);
    }
  }

  ntv = 0; tv = 0.0; gt = 1.0;

  for (i=0; i<nc; i++) y1[i] = g0(xx[i]); // initial profile

// Time loop:

  do {
    ntv++; tv += tau;

    for (i=0; i<nc; i++) y0[i] = y1[i];

    if (np>1)
      BndAExch1D(mp_l,1,y0+  0,&y0m,
                 mp_r,1,y0+ncm,&y0p);
    else
      { y0m = 0.0; y0p = 0.0; }

    for (i=0; i<nc; i++) {
      ii = i1 + i;

      if (ii==0)
        y1[i] = g1(tv);
      else if (ii==nx)
        y1[i] = g2(tv);
      else {
        if (i==0) s1 = aa[i] * (y0[i] - y0m);
        else      s1 = aa[i] * (y0[i] - y0[i-1]);

        if (i==ncm) s2 = bb[i] * (y0p - y0[i]);
        else        s2 = bb[i] * (y0[i+1] - y0[i]);

        s3 = tau * f(xx[i],tv-tau);

        y1[i] += (s2 - s1 + s3);
      }
    }

    if (ntv % ntp == 0) {
      gt = 0.0;
      for (i=0; i<nc; i++) {
        s0 = (y1[i]/y0[i]-1.0);
        gt = dmax(gt,dabs(s0));
      }
      gt = gt / tau;

      if (np>1) {
        s0 = gt; MPI_Allreduce(&s0,&gt,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      }

      if (mp == 0) {
        t2 = MPI_Wtime() - t1;
        fprintf(stderr,"ntv=%d tv=%le gt=%le tcpu=%le\n",ntv,tv,gt,t2);
      }
    }

    if (lp>0) {
      fprintf(Fo,"ntv=%d tv=%le gt=%le\n",ntv,tv,gt);
      for (i=0; i<nc; i++)
        fprintf(Fo,"i=%8d x=%12le y1=%12le\n",(i1+i),xx[i],y1[i]);
    }

  } while ((ntv<ntm) && (gt>epst));

  t1 = MPI_Wtime() - t1;

  sprintf(sname,"%s_%02d.dat",vname,np);
  OutFun1DP(sname,np,mp,nc,xx,y1);

  fprintf(Fo,"ntv=%d tv=%le gt=%le time=%le\n",ntv,tv,gt,t1);
  if (mp == 0) fprintf(stderr,"ntv=%d tv=%le gt=%le tcpu=%le\n",ntv,tv,gt,t1);

  ier = fclose_m(&Fo);

  MPI_Finalize();
  return 0;
}
