#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mycom.h"
#include <omp.h>

int nt, mt;

double a = 0;
double b = 1;
int ni = 1000000000;
int num = 0;
double sum = 0;

double f1(double x);
double f1(double x) { return 4.0/(1.0+x*x); }

void myjobt(int mt);
void myjobt(int mt)
{
  int n1; double a1, b1, h1, s;

  n1 = ni / nt;

  h1 = (b - a) / nt;

  a1 = a + h1 * mt;

  if (mt<nt-1) b1 = a1 + h1; else b1 = b;

  s = integrate(f1,a1,b1,n1);

  printf("mt=%d a1=%le b1=%le n1=%d s1=%le\n",mt,a1,b1,n1,s);

  #pragma omp critical
  {
    sum += s;
    num++;

  } // end critical

  return;
}

int main(int argc, char *argv[])
{
  double t;

  if (argc<2){
    printf("Usage: %s <thread number>\n",argv[0]);
    return 1;
  }

  sscanf(argv[1],"%d",&nt);
  if (nt<1) nt = 1; mt = 0;

  t = mytime(0);

  if (nt<2)
    sum = integrate(f1,a,b,ni);
  else {

    omp_set_num_threads(nt);

    #pragma omp parallel
    {
      int mt = omp_get_thread_num();

      myjobt(mt);

    } // end parallel

    while (num<nt);
  }

  t = mytime(1);

  printf("time=%lf sum=%22.15le\n",t,sum);

  return 0;
}
