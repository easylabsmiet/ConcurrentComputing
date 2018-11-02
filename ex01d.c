#include <stdio.h>
#include <sys/times.h>
#include <unistd.h>

struct tms tmsBegin, tmsEnd;
long tick;

void time_tick();
void time_tick()
{
  tick = sysconf(_SC_CLK_TCK);
  return;
}

void time_start();
void time_start()
{
  times(&tmsBegin);
  return;
}

double time_stop();
double time_stop()
{
  times(&tmsEnd);
  return ((tmsEnd.tms_utime-tmsBegin.tms_utime)+
          (tmsEnd.tms_stime-tmsBegin.tms_stime))*1.0/tick;
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
  double t, mypi;

  time_tick();
  time_start();
  mypi = pi_calculate(1000000000);
  t = time_stop();
  printf("Time: %lf sec Pi = %14.12lf\n",t,mypi);

  return 0;
}
