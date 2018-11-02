#include <stdio.h>
#include <time.h>
#include <unistd.h>

double f(double x);
double f(double x)
{
  return (4.0 / (1.0 + x*x));
}

double pi_calculate(const int n);
double pi_calculate(const int n)
{
  int i;
  double x,h,sum;
  h = 1.0 / (double) n;
  sum = 0.0;
  for (i=1; i<=n; i++){
    x = h * ((double)i - 0.5);
    sum += f(x);
  }
  return h*sum;
}

int main(int argc, char *argv[])
{	
  long tick;
  double t, mypi;

  tick = sysconf(_SC_CLK_TCK);

  t = (double)(clock()*0.0001/tick);
  mypi = pi_calculate(1000000000);
  t = (double)(clock()*0.0001/tick)-t;
  printf("Time: %lf sec Pi = %14.12lf\n",t,mypi);

  return 0;
}
