#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mycom.h"
#include <pthread.h>

typedef struct tag_data_t {
  int n, nt, mt;
  double a, b, s, *sum;
} data_t;

int nt, mt;
pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
pthread_t* threads;
data_t* data;

double a = 0;
double b = 1;
int ni = 1000000000;
double sum = 0;

double f1(double x);
double f1(double x) { return 4.0/(1.0+x*x); }

void* myjobt(void* d);
void* myjobt(void* d)
{
  int n1; double a1, b1, h1;
  data_t* dd = (data_t *)d; 

  n1 = dd->n / dd->nt;

  h1 = (dd->b - dd->a) / dd->nt;

  a1 = dd->a + h1 * dd->mt;

  if (dd->mt<dd->nt-1) b1 = a1 + h1; else b1 = dd->b;

  dd->s = integrate(f1,a1,b1,n1);

  printf("mt=%d a1=%le b1=%le n1=%d s1=%le\n",dd->mt,a1,b1,n1,dd->s);

  pthread_mutex_lock(&mut); // lock

  *dd->sum += dd->s;

  pthread_mutex_unlock(&mut); // unlock

  return 0;
}

void ThreadInit();
void ThreadInit()
{
  int i;

  if (!(threads = (pthread_t*) malloc(nt*sizeof(pthread_t))))
    myerr("Not enough memory",1);

  if (!(data = (data_t*) malloc(nt*sizeof(data_t))))
    myerr("Not enough memory",1);

  for (i=0; i<nt; i++){
    (data+i)->a=a;
    (data+i)->b=b;
    (data+i)->n=ni;
    (data+i)->nt=nt;
    (data+i)->mt=i;
    (data+i)->sum = &sum;
    if (pthread_create(threads+i,0,myjobt,(void*)(data+i)))
      myerr("Can not create thread",2);
  }

  return;
}

void ThreadDone();
void ThreadDone()
{
  int i;

  for (i=0; i<nt; i++)
    if (pthread_join(threads[i],0))
      myerr("Can not wait thread",3);

  free(data);

  free(threads);

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
    ThreadInit();

    ThreadDone();
  }

  t = mytime(1);

  printf("time=%lf sum=%22.15le\n",t,sum);

  return 0;
}
