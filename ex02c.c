#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include "mycom.h"

double f1(double x);
double f1(double x)
{
  return (4.0 / (1.0 + x*x));
}

double f2(double x);
double f2(double x)
{
  return (4.0 / (1.0 + x*x/3.0));
}

double f3(double x);
double f3(double x)
{
  return (4.0 / (1.1 + exp(x) + sin(x)));
}

double f4(double x);
double f4(double x)
{
  return (4.0 / (1.1 + 2.2*exp(x) + 3.3*sin(x)));
}

int main(int argc, char *argv[])
{	
  int nc=1000000000;
  double t1, t2, dt, sum;

  t1 = mytime(0); sum = integrate(f1,0.0,1.0,nc); t1 = mytime(1);
  t2 = mytime(0); sum = integrate(f2,0.0,1.0,nc); t2 = mytime(1);
  dt = 1.0/dabs(t2-t1);
  printf("Time: %lf %lf sec Div. perf.: %le GFlops\n",t1,t2,dt);

  t1 = mytime(0); sum = integrate(f3,0.0,1.0,nc); t1 = mytime(1);
  t2 = mytime(0); sum = integrate(f4,0.0,1.0,nc); t2 = mytime(1);
  dt = 2.0/dabs(t2-t1);
  printf("Time: %lf %lf sec Mult. perf.: %le GFlops\n",t1,t2,dt);

  return 0;
}
