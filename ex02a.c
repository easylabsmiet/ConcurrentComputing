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
  return (4.0 / (1.0 + x + x*x));
}

double f3(double x);
double f3(double x)
{
  return (4.0 / (1.0 + 3.0*x*x));
}

double f4(double x);
double f4(double x)
{
  return (4.0 / (1.0 + x*x/3.0));
}

int main(int argc, char *argv[])
{	
  int nc=1000000000;
  double t1, t2, sum;

  t1 = mytime(0); sum = integrate(f1,0.0,1.0,nc); t1 = mytime(1);
  printf("Time: %lf sec sum = %14.12lf\n",t1,sum);

  t2 = mytime(0); sum = integrate(f2,0.0,1.0,nc); t2 = mytime(1);
  printf("Time: %lf sec sum = %14.12lf diff: %lf\n",t2,sum,t2-t1);

  t2 = mytime(0); sum = integrate(f3,0.0,1.0,nc); t2 = mytime(1);
  printf("Time: %lf sec sum = %14.12lf diff: %lf\n",t2,sum,t2-t1);

  t2 = mytime(0); sum = integrate(f4,0.0,1.0,nc); t2 = mytime(1);
  printf("Time: %lf sec sum = %14.12lf diff: %lf\n",t2,sum,t2-t1);

  return 0;
}
