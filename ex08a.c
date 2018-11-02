#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "mynet.h"

double f1(double x);
double f1(double x) {
  double e=exp(x);
  return 0.5*(e+1.0/e);
}

double f2(double x);
double f2(double x) {
  double e=exp(x);
  return 0.5*(e-1.0/e);
}

double f(double x, double y);
double f(double x, double y) {
  return f1(x)*f2(y);
}

int np, mp, nl, ier, lp;
int np1, np2, mp1, mp2, iflag;
char pname[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;
double tick, t1, t2, t3;

double a=0.0, b=1.0, c=0.0, d=1.0;
int nx=10000, ny=10000;

int main(int argc, char *argv[])
{
  int i, i1, i2, j, j1, j2;
  double s, p, hx, hy, hxy, x, y;

  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  if (argc<3) mpierr("Too small arguments",1);

  sscanf(argv[1],"%d",&np1);
  sscanf(argv[2],"%d",&np2);
  sscanf(argv[3],"%d",&iflag);

  if (np1*np2!=np)
    mpierr("Bad processor grid",2);

  mp2 = mp / np1;
  mp1 = mp % np1;

  fprintf(stderr,"mp=%d grid=%dx%d coord=(%d,%d)\n",mp,np1,np2,mp1,mp2);

  t1 = MPI_Wtime();

  hx=(b-a)/nx; hy=(d-c)/ny; hxy=hx*hy;

  if (np1==1) {
    i1 = 0; i2 = nx-1;
  }
  else {
    i1 = mp1 * (nx/np1); if (mp1<np1-1) i2 = i1+(nx/np1)-1; else i2 = nx-1;
  }

  if (np2==1) {
    j1 = 0; j2 = ny-1;
  }
  else {
    j1 = mp2 * (ny/np2); if (mp2<np2-1) j2 = j1+(ny/np2)-1; else j2 = ny-1;
  }

  fprintf(stderr,"mp=%d i1=%d i2=%d j1=%d j2=%d\n",mp,i1,i2,j1,j2);

  s = 0;

  for (j=j1; j<=j2; j++) {
    y = c + (j*1.0+.5)*hy;
    for (i=i1; i<=i2; i++) {
      x = a + (i*1.0+.5)*hx;
      s = s + hxy * f(x,y);
    }
  }

  t2 = MPI_Wtime();

  t1 = t2 - t1;

  if (np==1)
    t2 = 0;
  else {
    p = s;
    if (iflag==0)
      MPI_Reduce(&p, &s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    else	
      MPI_Allreduce(&p, &s, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    t2 = MPI_Wtime()-t2;
  }

  t3 = t1 + t2;

  fprintf(stderr,"mp=%d t1=%le t2=%le t3=%le int=%le\n",mp,t1,t2,t3,s);

  MPI_Finalize();
  return 0;
}
