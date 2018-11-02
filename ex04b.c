#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mycom.h"

#include <pthread.h>
#define LBUF     10

int nt, mt;
pthread_t* threads;
double sums[LBUF];

double a = 0;
double b = 1;
int ni = 1000000000;
double sum = 0;

double f1(double x);
double f1(double x) { return 4.0/(1.0+x*x); }

void* myjobt(void* m);
void* myjobt(void* m)
{
  int n1, mt=(int)(((long long int)m)%1024); double a1, b1, h1;

  n1 = ni / nt;
  h1 = (b - a) / nt;
  a1 = a + h1 * mt;
  if (mt<nt-1) b1 = a1 + h1; else b1 = b;

  sums[mt] = integrate(f1,a1,b1,n1);

  printf("mt=%d a1=%le b1=%le n1=%d s1=%le\n",mt,a1,b1,n1,sums[mt]);

  return 0;
}

int main(int argc, char *argv[])
{
  int i; double t;

  if (argc<2){
    printf("Usage: %s <thread number>\n",argv[0]);
    return 1;
  }

  sscanf(argv[1],"%d",&nt); mt = 0;

  t = mytime(0);

  if (nt<2)
    sum = integrate(f1,a,b,ni);
  else {
    if (!(threads = (pthread_t*) malloc(nt*sizeof(pthread_t))))
      myerr("server: not enough memory",1);

    for (i=0; i<nt; i++)
      if (pthread_create(threads+i,0,myjobt,(void*)((long long int)i)))
        myerr("server: cannot create thread",2);
    
    for (i=0; i<nt; i++)
      if (pthread_join(threads[i],0))
        myerr("server: cannot wait thread",3);
      else
        sum = sum + sums[i];

    free(threads);
  }

  t = mytime(1);

  printf("time=%lf sum=%22.15le\n",t,sum);

  return 0;
}
