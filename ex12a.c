//  Nonlinear boundary problem:
//
//  (k(u)u')' - q(u)u = - f(u), xa < x < xb
//
//  u(xa) = ua, u(xb) = ub
//
//  Simple iterations
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mycom.h"
#include "mynet.h"
#include "myprog.h"

int np, mp, nl, ier, lp;
char pname[MPI_MAX_PROCESSOR_NAME];
char sname[10] = "ex12a.p00";
MPI_Status status;
union_t buf;
double tick, t1, t2, t3;

FILE *Fi = NULL;
FILE *Fo = NULL;

int nx, it, itm;
double xa, xb, ua, ub, alf, eps, pa, pa2, uc, rka;

double k(double u);
double k(double u) {
  return 1.0 + u*u;
}

double k1(double u);
double k1(double u) {
  return 2.0*u;
}

double q(double u);
double q(double u) {
  return 1.0 / (1.0 + u*u);
}

double u(double x);
double u(double x) {
  return ua*exp(pa*(x-xa));
}

double u1(double u);
double u1(double u) {
  return pa*u;
}

double u2(double u);
double u2(double u) {
  return pa2*u;
}

double f(double u);
double f(double u) {
  return q(u)*u - k1(u)*u1(u)*u1(u) - k(u)*u2(u);
}

int main(int argc, char *argv[])
{
  int i, j, i1, i2, nc, ncm, ncp, ncx;
  double hx, hx2, s0, s1, s2, sm, sp, y0m, y0p;
  double *xx, *aa, *bb, *cc, *dd, *ee, *ff, *al, *y0, *y1, *y2, *y3, *y4;

  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  sprintf(sname+7,"%02d",mp);
  ier = fopen_m(&Fo,sname,"wt");
  if (ier!=0) mpierr("Protocol file not opened",1);

  if (mp==0) {
    ier = fopen_m(&Fi,"ex12a.d","rt");
    if (ier!=0) mpierr("Data file not opened",2);
    fscanf(Fi,"xa=%le\n",&xa);
    fscanf(Fi,"xb=%le\n",&xb);
    fscanf(Fi,"ua=%le\n",&ua);
    fscanf(Fi,"ub=%le\n",&ub);
    fscanf(Fi,"alf=%le\n",&alf);
    fscanf(Fi,"eps=%le\n",&eps);
    fscanf(Fi,"itm=%d\n",&itm);
    fscanf(Fi,"nx=%d\n",&nx);
    fscanf(Fi,"lp=%d\n",&lp);
    fclose_m(&Fi);
    if (argc>1) sscanf(argv[1],"%d",&nx);
  }

  if (np>1) {
    if (mp==0) {
      buf.ddata[0] = xa; buf.ddata[1] = xb;
      buf.ddata[2] = ua; buf.ddata[3] = ub;
      buf.ddata[4] = alf; buf.ddata[5] = eps;
      buf.idata[12] = itm; buf.idata[13] = nx;
      buf.idata[14] = lp;
    }
    MPI_Bcast(buf.ddata,8,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (mp>0) {
      xa = buf.ddata[0]; xb = buf.ddata[1];
      ua = buf.ddata[2]; ub = buf.ddata[3];
      alf = buf.ddata[4]; eps = buf.ddata[5];
      itm = buf.idata[12]; nx = buf.idata[13];
      lp = buf.idata[14];
    }
  }

  ub = ua * exp(alf);

  fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",
    np,mp,pname,tick);

  fprintf(Fo,"xa=%le xb=%le ua=%le ub=%le alf=%le eps=%le itm=%d nx=%d lp=%d\n",
    xa,xb,ua,ub,alf,eps,itm,nx,lp);

  t1 = MPI_Wtime();

  pa = alf / (xb-xa); pa2 = pa*pa; uc = (ub-ua) / (xb-xa);

  hx = (xb-xa)/nx; hx2 = hx * hx;

  MyRange(np,mp,0,nx,&i1,&i2,&nc);
  ncm = nc-1; ncp = 2*(np-1); ncx = imax(nc,ncp);

  fprintf(Fo,"i1=%d i2=%d nc=%d\n",i1,i2,nc);

  xx = (double*)(malloc(sizeof(double)*nc));
  y0 = (double*)(malloc(sizeof(double)*nc));
  y1 = (double*)(malloc(sizeof(double)*nc));

  aa = (double*)(malloc(sizeof(double)*ncx));
  bb = (double*)(malloc(sizeof(double)*ncx));
  cc = (double*)(malloc(sizeof(double)*ncx));
  ff = (double*)(malloc(sizeof(double)*ncx));
  al = (double*)(malloc(sizeof(double)*ncx));

  if (np>1) {
    y2 = (double*)(malloc(sizeof(double)*nc));
    y3 = (double*)(malloc(sizeof(double)*nc));
    y4 = (double*)(malloc(sizeof(double)*ncp));
    dd = (double*)(malloc(sizeof(double)*4*ncp));
    ee = (double*)(malloc(sizeof(double)*4*ncp));
  }

  for (i=0; i<nc; i++) xx[i] = xa + hx * (i1 + i); // grid

  it = 0;

  for (i=0; i<nc; i++) y1[i] = ua + uc * (xx[i]-xa); // start solution

// Iterations:

  do {
    for (i=0; i<nc; i++) y0[i] = y1[i];

    if (np>1)
      BndExch1D(np,mp,1,1,1,1,&(y0[0]),&y0m,&(y0[ncm]),&y0p);
    else {
      y0m = 0.0; y0p = 0.0;
    }

    if (mp==0) {
      aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = ua;
    }
    else {
      s0 = k(y0[0]); s1 = k(y0m); s2 = k(y0[1]);
      aa[0] = 0.5 * (s0 + s1);
      bb[0] = 0.5 * (s0 + s2);
      cc[0] = hx2 * q(y0[0]) + aa[0] + bb[0];
      ff[0] = hx2 * f(y0[0]);
    }

    for (i=1; i<ncm; i++) {
      s0 = k(y0[i]); s1 = k(y0[i-1]); s2 = k(y0[i+1]);
      aa[i] = 0.5 * (s0 + s1);
      bb[i] = 0.5 * (s0 + s2);
      cc[i] = hx2 * q(y0[i]) + aa[i] + bb[i];
      ff[i] = hx2 * f(y0[i]);
    }

    if (mp==np-1) {
      aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = ub;
    }
    else {
      s0 = k(y0[ncm]); s1 = k(y0[ncm-1]); s2 = k(y0p);
      aa[ncm] = 0.5 * (s0 + s1);
      bb[ncm] = 0.5 * (s0 + s2);
      cc[ncm] = hx2 * q(y0[ncm]) + aa[ncm] + bb[ncm];
      ff[ncm] = hx2 * f(y0[ncm]);
    }

    rka = 0.0;

    for (i=0; i<nc; i++) {
      s0 = y0[i];
      if (i==  0) { if (mp>   0) sm = y0m; else sm = 0.0; } else sm = y0[i-1];
      if (i==ncm) { if (mp<np-1) sp = y0p; else sp = 0.0; } else sp = y0[i+1];
      s1 = ff[i] + aa[i] * sm + bb[i] * sp - cc[i] * s0;
      rka = dmax(rka,dabs(s1));
    }

    if (np>1) {
      s0 = rka; MPI_Allreduce(&s0,&rka,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    }

    if (lp>0) {
      if (mp==0) fprintf(stderr,"it=%d rka=%le\n",it,rka);
      fprintf(Fo,"it=%d rka=%le\n",it,rka);
    }

    if (rka<=eps) break;

    it++;

    if (lp>0)
      for (i=0; i<nc; i++)
        fprintf(Fo,"i=%8d a=%12le b=%12le c=%12le f=%12le\n",
          (i1+i),aa[i],bb[i],cc[i],ff[i]);

    if (np<2) ier = prog_right(nc,aa,bb,cc,ff,al,y1);
    else      ier = prog_rightp(np,mp,nc,aa,bb,cc,ff,al,y1,y2,y3,y4,dd,ee);

    if (ier!=0) mpierr("Bad solution",1);

    if (lp>0)
      for (i=0; i<nc; i++)
        fprintf(Fo,"i=%8d x=%12le y0=%12le y1=%12le\n",
          (i1+i),xx[i],y0[i],y1[i]);

  } while (it<=itm);

  t1 = MPI_Wtime() - t1;

  s0 = 0.0;
  for (i=0; i<nc; i++) {
    s1 = u(xx[i]); s2 = dabs(s1-y1[i]); s0 = dmax(s0,s2);
    if (lp>0)
      fprintf(Fo,"i=%8d x=%12le y=%12le u=%12le d=%12le\n",
        (i1+i),xx[i],y1[i],s1,s2);
  }

  if (np>1) {
    s1 = s0; MPI_Allreduce(&s1,&s0,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  }

  if (mp==0) fprintf(stderr,"nx=%d it=%d time=%le dmax=%le\n",nx,it,t1,s0);
  fprintf(Fo,"nx=%d it=%d time=%le dmax=%le\n",nx,it,t1,s0);

  ier = fclose_m(&Fo);

  MPI_Finalize();
  return 0;
}
