#include <stdio.h>

#define read_cpu_tick(high,low)   \
        asm("rdtsc                            \n\
             mov %%eax,%0                     \n\
             mov %%edx,%1": "=m" (low),"=m" (high))

static double freq;

double my_cpu_time();
double my_cpu_time()
{
  unsigned long h=0, l=0;
  double t;
  read_cpu_tick(h,l);
  t = 1e-6*l;
  t += 4294.967296*h;
  return (t/freq);
}

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
  double t, t1, t2, mypi;
  freq = 2000.0;
  sleep(1);
  t1 = my_cpu_time();
  printf("%lf\n",t1);
  mypi = pi_calculate(1000000000);
  t2 = my_cpu_time();
  printf("%lf\n",t2);
  t = t2 - t1;
  printf("Time: %lf sec Pi = %14.12lf\n",t,mypi);
  return 0;
}
