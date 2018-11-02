#include "stdlib.h"
#include "mpi.h"
#include "mynet.h"
#include "myprog.h"

//
//  n=1 => c[0]*y[0]=f[0]
//  n=2 => c[0]*y[0]-b[0]*y[1]=f[0], -a[1]*y[0]+c[1]*y[1]=f[1]
//  n>2 => classic system:
//               c[i]*y[i]-b[i]*y[i+1]=f[i], i=0
//  -a[i]*y[i-1]+c[i]*y[i]-b[i]*y[i+1]=f[i], 0<i<n-1
//  -a[i]*y[i-1]+c[i]*y[i]            =f[i], i=n-1
//

int prog_right(int n, double* a, double* b, double* c,
                      double* f, double* al, double* y)
{
  int i; double s;

  if (n<1) return -1;

  if (n==1) {
    if (c[0]==0.0) return -2;
    y[0] = f[0]/c[0];
  }
  else if (n==2) {
    s = c[0]*c[1]-b[0]*a[1];
    if (s==0.0) return -2;
    y[0]=(c[1]*f[0]+b[0]*f[1])/s;
    y[1]=(a[1]*f[0]+c[0]*f[1])/s;
  }
  else {
    if (c[0]==0.0) return -2;
    al[0] = b[0]/c[0];
    y[0] = f[0]/c[0];
    for (i=1; i<n; i++) {
      s = c[i]-a[i]*al[i-1];
      if (s==0.0) return -2;
      al[i] = b[i]/s;
      y[i] = (f[i]+a[i]*y[i-1])/s;
    }
    for (i=n-2; i>=0; i--)
      y[i] = al[i]*y[i+1]+y[i];
  }
  return 0;
}

//
// alternative memory model:
// a[4*i],a[4*i+1],a[4*i+2],a[4*i+3] <=> a[i],b[i],c[i],f[i]
//

int prog_rightm(int n, double* a, double* al, double* y)
{
  int i,m; double s;

  if (n<1) return -1;

  if (n==1) {
    if (a[2]==0.0) return -2;
    y[0] = a[3]/a[2];
  }
  else if (n==2) {
    s = a[2]*a[6]-a[1]*a[4];
    if (s==0.0) return -2;
    y[0]=(a[6]*a[3]+a[1]*a[7])/s;
    y[1]=(a[4]*a[3]+a[2]*a[7])/s;
  }
  else {
    if (a[2]==0.0) return -2;
    al[0] = a[1]/a[2];
    y[0] = a[3]/a[2];
    for (i=1; i<n; i++) {
      m = 4*i;
      s = a[m+2]-a[m]*al[i-1];
      if (s==0.0) return -2;
      al[i] = a[m+1]/s;
      y[i] = (a[m+3]+a[m]*y[i-1])/s;
    }
    for (i=n-2; i>=0; i--)
      y[i] = al[i]*y[i+1]+y[i];
  }
  return 0;
}

// Parallel variant:

int prog_rightp(int np, int mp, int nc,
                double *aa, double *bb, double *cc, double *ff,
                double *al, double *y1, double *y2, double *y3,
                double *y4, double *dd, double *ee)
{
  int i, j, ncm=nc-1, ncp=2*np-2;
  double a0, b0, c0, f0, a1, b1, c1, f1;

  a0 = aa[0]; b0 = bb[0]; c0 = cc[0]; f0 = ff[0];
  a1 = aa[ncm]; b1 = bb[ncm]; c1 = cc[ncm]; f1 = ff[ncm];

  if (mp==0) {
    aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y1); if (i!=0) return i;

    for (i=0; i<ncm; i++) ff[i] = 0.0; ff[ncm] = 1.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y2); if (i!=0) return i;
  }
  else if (mp<np-1) {
    aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = 0.0;
    aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y1); if (i!=0) return i;

    for (i=0; i<ncm; i++) ff[i] = 0.0; ff[ncm] = 1.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y2); if (i!=0) return i;

    ff[0] = 1.0; for (i=1; i<=ncm; i++) ff[i] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y3); if (i!=0) return i;
  }
  else {
    aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y1); if (i!=0) return i;

    ff[0] = 1.0; for (i=1; i<=ncm; i++) ff[i] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y3); if (i!=0) return i;
  }

  aa[0] = a0; bb[0] = b0; cc[0] = c0; ff[0] = f0;
  aa[ncm] = a1; bb[ncm] = b1; cc[ncm] = c1; ff[ncm] = f1;

  for (i=0; i<4*ncp; i++) dd[i] = 0;
  for (i=0; i<4*ncp; i++) ee[i] = 0;

  if (mp==0) {
    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = 0.0;
    ee[0] = a1;
    ee[1] = b1;
    ee[2] = c1;
    ee[3] = f1;
  }
  else if (mp<np-1) {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = b0 * y2[1];
    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = a1 * y3[ncm-1];
    i = mp * 8 - 4;
    ee[i]   = a0;
    ee[i+1] = b0;
    ee[i+2] = c0;
    ee[i+3] = f0;
    ee[i+4] = a1;
    ee[i+5] = b1;
    ee[i+6] = c1;
    ee[i+7] = f1;
  }
  else {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = 0.0;
    i = mp * 8 - 4;
    ee[i]   = a0;
    ee[i+1] = b0;
    ee[i+2] = c0;
    ee[i+3] = f0;
  }

  MPI_Allreduce(ee,dd,4*ncp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

//  i = prog_rightm(ncp,dd,al,y4); if (i!=0) return i;

  for (i=0; i<ncp; i++) {
    j = 4*i;
    aa[i] = dd[j];   bb[i] = dd[j+1];
    cc[i] = dd[j+2]; ff[i] = dd[j+3];
  }

  i = prog_right(ncp,aa,bb,cc,ff,al,y4); if (i!=0) return i;

  if (mp==0){
    b1 = y4[0];
    for (i=0; i<nc; i++) y1[i] = y1[i] + b1 * y2[i];
  }
  else if (mp<np-1) {
    a1 = y4[2*mp-1]; b1 = y4[2*mp];
    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i] + b1 * y2[i];
  }
  else {
    a1 = y4[2*mp-1];
    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i];
  }

  return 0;
}

int prog_rightpm(int np, int mp, int nc, int ip,
                 double *aa, double *bb, double *cc, double *ff,
                 double *al, double *y1, double *y2, double *y3, double *y4)
{
  int i, j, ncm=nc-1, ncp=2*np-2;
  double a0, b0, c0, f0, a1, b1, c1, f1;
  double *dd=y4+ncp, *ee=dd+4*ncp;

  a0 = aa[0]; b0 = bb[0]; c0 = cc[0]; f0 = ff[0];
  a1 = aa[ncm]; b1 = bb[ncm]; c1 = cc[ncm]; f1 = ff[ncm];

  if (mp==0) {
    aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y1); if (i!=0) return i;

    if (ip<2) {
      for (i=0; i<ncm; i++) ff[i] = 0.0; ff[ncm] = 1.0;
      i = prog_right(nc,aa,bb,cc,ff,al,y2); if (i!=0) return i;
    }
  }
  else if (mp<np-1) {
    aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = 0.0;
    aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y1); if (i!=0) return i;

    if (ip<2) {
      for (i=0; i<ncm; i++) ff[i] = 0.0; ff[ncm] = 1.0;
      i = prog_right(nc,aa,bb,cc,ff,al,y2); if (i!=0) return i;

      ff[0] = 1.0; for (i=1; i<=ncm; i++) ff[i] = 0.0;
      i = prog_right(nc,aa,bb,cc,ff,al,y3); if (i!=0) return i;
    }
  }
  else {
    aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y1); if (i!=0) return i;

    if (ip<2) {
      ff[0] = 1.0; for (i=1; i<=ncm; i++) ff[i] = 0.0;
      i = prog_right(nc,aa,bb,cc,ff,al,y3); if (i!=0) return i;
    }
  }

  aa[0] = a0; bb[0] = b0; cc[0] = c0; ff[0] = f0;
  aa[ncm] = a1; bb[ncm] = b1; cc[ncm] = c1; ff[ncm] = f1;

  for (i=0; i<4*ncp; i++) dd[i] = 0;
  for (i=0; i<4*ncp; i++) ee[i] = 0;

  if (mp==0) {
    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = 0.0;
    ee[0] = a1;
    ee[1] = b1;
    ee[2] = c1;
    ee[3] = f1;
  }
  else if (mp<np-1) {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = b0 * y2[1];
    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = a1 * y3[ncm-1];
    i = mp * 8 - 4;
    ee[i]   = a0;
    ee[i+1] = b0;
    ee[i+2] = c0;
    ee[i+3] = f0;
    ee[i+4] = a1;
    ee[i+5] = b1;
    ee[i+6] = c1;
    ee[i+7] = f1;
  }
  else {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = 0.0;
    i = mp * 8 - 4;
    ee[i]   = a0;
    ee[i+1] = b0;
    ee[i+2] = c0;
    ee[i+3] = f0;
  }

  MPI_Allreduce(ee,dd,4*ncp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  i = prog_rightm(ncp,dd,al,y4); if (i!=0) return i;

  if (mp==0){
    b1 = y4[0];
    for (i=0; i<nc; i++) y1[i] = y1[i] + b1 * y2[i];
  }
  else if (mp<np-1) {
    a1 = y4[2*mp-1]; b1 = y4[2*mp];
    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i] + b1 * y2[i];
  }
  else {
    a1 = y4[2*mp-1];
    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i];
  }

  return 0;
}

int prog_rightpn(int np, int mp, MPI_Comm cm, int nc, int ip,
                 double *aa, double *bb, double *cc, double *ff,
                 double *al, double *y1, double *y2, double *y3, double *y4)
{
  int i, j, ncm=nc-1, ncp=2*np-2;
  double a0, b0, c0, f0, a1, b1, c1, f1;
  double *dd=y4+ncp, *ee=dd+4*ncp;

  if (np<2) {
    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    return i;
  }

  a0 = aa[0]; b0 = bb[0]; c0 = cc[0]; f0 = ff[0];
  a1 = aa[ncm]; b1 = bb[ncm]; c1 = cc[ncm]; f1 = ff[ncm];

  if (mp==0) {
    aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y1); if (i!=0) return i;

    if (ip<2) {
      for (i=0; i<ncm; i++) ff[i] = 0.0; ff[ncm] = 1.0;
      i = prog_right(nc,aa,bb,cc,ff,al,y2); if (i!=0) return i;
    }
  }
  else if (mp<np-1) {
    aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = 0.0;
    aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y1); if (i!=0) return i;

    if (ip<2) {
      for (i=0; i<ncm; i++) ff[i] = 0.0; ff[ncm] = 1.0;
      i = prog_right(nc,aa,bb,cc,ff,al,y2); if (i!=0) return i;

      ff[0] = 1.0; for (i=1; i<=ncm; i++) ff[i] = 0.0;
      i = prog_right(nc,aa,bb,cc,ff,al,y3); if (i!=0) return i;
    }
  }
  else {
    aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = 0.0;
    i = prog_right(nc,aa,bb,cc,ff,al,y1); if (i!=0) return i;

    if (ip<2) {
      ff[0] = 1.0; for (i=1; i<=ncm; i++) ff[i] = 0.0;
      i = prog_right(nc,aa,bb,cc,ff,al,y3); if (i!=0) return i;
    }
  }

  aa[0] = a0; bb[0] = b0; cc[0] = c0; ff[0] = f0;
  aa[ncm] = a1; bb[ncm] = b1; cc[ncm] = c1; ff[ncm] = f1;

  for (i=0; i<4*ncp; i++) dd[i] = 0;
  for (i=0; i<4*ncp; i++) ee[i] = 0;

  if (mp==0) {
    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = 0.0;
    ee[0] = a1;
    ee[1] = b1;
    ee[2] = c1;
    ee[3] = f1;
  }
  else if (mp<np-1) {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = b0 * y2[1];
    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = a1 * y3[ncm-1];
    i = mp * 8 - 4;
    ee[i]   = a0;
    ee[i+1] = b0;
    ee[i+2] = c0;
    ee[i+3] = f0;
    ee[i+4] = a1;
    ee[i+5] = b1;
    ee[i+6] = c1;
    ee[i+7] = f1;
  }
  else {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = 0.0;
    i = mp * 8 - 4;
    ee[i]   = a0;
    ee[i+1] = b0;
    ee[i+2] = c0;
    ee[i+3] = f0;
  }

  MPI_Allreduce(ee,dd,4*ncp,MPI_DOUBLE,MPI_SUM,cm);

  i = prog_rightm(ncp,dd,al,y4); if (i!=0) return i;

  if (mp==0){
    b1 = y4[0];
    for (i=0; i<nc; i++) y1[i] = y1[i] + b1 * y2[i];
  }
  else if (mp<np-1) {
    a1 = y4[2*mp-1]; b1 = y4[2*mp];
    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i] + b1 * y2[i];
  }
  else {
    a1 = y4[2*mp-1];
    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i];
  }

  return 0;
}
