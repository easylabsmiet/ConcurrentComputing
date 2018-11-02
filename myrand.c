#include <stdlib.h>
#include <math.h>
#include "myrand.h"

double myrand()
{
  int i; double s;
  i = rand(); s = (1.0 * i) / 65536.0;
  return s;
}

int mycomp(const void* a, const void* b)
{
  if (*((double*)a)>*((double*)b))
    return 1;
  else {
    if (*((double*)a)<*((double*)b))
      return -1;
    else
      return 0;
  }
}

void myindexes(const int np, const int ns, int* ind)
{
  if (np==1) { ind[0] = 0; ind[1] = ns-1; ind[2] = ns; }
  else {
    int ii, i1, i2, ni, mi;
    ni = ns / np; mi = ns - ni * np;
    for(ii=0; ii<np; ii++) {
      if (ii+1<=mi) { i1 = ii * (ni+1); i2 = i1 + ni; }
      else { i1 = mi * (ni+1) + (ii-mi) * ni; i2 = i1 + ni - 1; }
      ind[3*ii  ] = i1;
      ind[3*ii+1] = i2;
      ind[3*ii+2] = i2-i1+1;
    }
  }
  return;
}

void mychebset(int mc, int nc, double *rr)
{
  int i, j, m, *ii;
  double h, *dd;

  if (mc<1) { rr[0] = 1.0; return; }

  ii = (int*)(malloc(sizeof(int)*nc));
  dd = (double*)(malloc(sizeof(double)*2*nc));

  h = M_PI / (2.0 * nc);

  for (i=0; i<nc; i++) dd[i] = 1.0 + cos(h*(2*i+1));

  qsort(dd,nc,sizeof(double),mycomp);

  for (i=0; i<nc; i++) dd[nc+i] = 1.0 / dd[i];

  ii[0] = 1; ii[1] = 0;

  if (nc>2) {
    m = 1;
    for (i=1; i<mc; i++) {
      m = m * 2;
      for (j=m; j>0; j--) ii[2*j-1] = ii[j-1];
      for (j=0; j<m; j++) ii[2*j] = 2*m - 1 - ii[2*j+1];
    }
  }

  for (i=0; i<nc; i++) { j = nc + ii[i]; dd[i] = dd[j]; }

  for (i=0; i<nc; i++) rr[i] = dd[i];

  free(dd); free(ii);

  return;
}
