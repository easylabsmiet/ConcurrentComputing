#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "mycom.h"
#include "mynet.h"

int np, mp, nl, ier, lp;
char pname[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;
double tick, t1, t2, t3;

double a = 0;
double b = 1;
int ni = 1000000000;
double sum = 0;

double f1(double x);
double f1(double x) { return 4.0/(1.0+x*x); }

double myjob(int mp);
double myjob(int mp)
{
  int n1; double a1, b1, h1, s1;

  n1 = ni / np;
  h1 = (b - a) / np;
  a1 = a + h1 * mp;
  if (mp<np-1) b1 = a1 + h1; else b1 = b;

  s1 = integrate(f1,a1,b1,n1);

  return s1;
}

int main(int argc, char* argv[])
{
  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  if (np<2) {
    t1 = MPI_Wtime();
    sum = integrate(f1,a,b,ni);
    t2 = MPI_Wtime();
    t3 = t2;
  }
  else {
    int i; double p;

    t1 = MPI_Wtime();
    sum = myjob(mp);
    t2 = MPI_Wtime();

    if (mp==0)
      for (i=1; i<np; i++) {
        MPI_Recv(&p, 1, MPI_DOUBLE, i, MY_TAG, MPI_COMM_WORLD, &status);
        sum = sum + p;
      }
    else
      MPI_Send(&sum, 1, MPI_DOUBLE, 0, MY_TAG, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    t3 = MPI_Wtime();
  }

  t1 = t2 - t1;
  t2 = t3 - t2;
  t3 = t1 + t2;

  fprintf(stderr,"mp=%d t1=%lf t2=%lf t3=%lf int=%22.15le\n",mp,t1,t2,t3,sum);

  MPI_Finalize();
  return 0;
}
