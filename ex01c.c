#include <stdio.h>
#include <sys/time.h>

struct timeval tv1,tv2,dtv;
struct timezone tz;

void time_start();
void time_start()
{
  gettimeofday(&tv1,&tz);
}

double time_stop();
double time_stop()
{
  gettimeofday(&tv2, &tz);
  dtv.tv_sec  = tv2.tv_sec  - tv1.tv_sec;
  dtv.tv_usec = tv2.tv_usec - tv1.tv_usec;
  if(dtv.tv_usec<0){
    dtv.tv_sec--;
    dtv.tv_usec+=1000000;
  }
  return (double)(dtv.tv_sec) + (double)(dtv.tv_usec)/1000000.0;
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

  time_start();
  mypi = pi_calculate(1000000000);
  t = time_stop();
  printf("Time: %lf sec Pi = %14.12lf\n",t,mypi);

  return 0;
}
