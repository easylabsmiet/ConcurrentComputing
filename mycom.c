#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/time.h>
#include <math.h>
#include "mycom.h"

int iabs(int a)
{
  if (a>=0)
    return a;
  else
    return (-a);
}

int imax(int a, int b)
{
  if (a>=b)
    return a;
  else
    return b;
}

int imin(int a, int b)
{
  if (a<=b)
    return a;
  else
    return b;
}

double dabs(double a)
{
  if (a>=0)
    return a;
  else
    return (-a);
}

double dmax(double a, double b)
{
  if (a>=b)
    return a;
  else
    return b;
}

double dmin(double a, double b)
{
  if (a<=b)
    return a;
  else
    return b;
}

double dsin(double x) {
  double s = sin(x);
  if (dabs(s)>1e-15)
    return s;
  else
    return 0.0;
}

double dcos(double x) {
  double s = cos(x);
  if (dabs(s)>1e-15)
    return s;
  else
    return 0.0;
}

double dexp(double x) {
  double s = exp(x);
  if (dabs(s)>1e-15)
    return s;
  else
    return 0.0;
}

double mytime(const int n)
{
  static struct timeval tv1,tv2,dtv;
  static struct timezone tz;
  if (n==0){
    gettimeofday(&tv1,&tz);
    return 0;
  }
  else{
    gettimeofday(&tv2, &tz);
    dtv.tv_sec  = tv2.tv_sec  - tv1.tv_sec;
    dtv.tv_usec = tv2.tv_usec - tv1.tv_usec;
    if(dtv.tv_usec<0){
      dtv.tv_sec--;
      dtv.tv_usec+=1000000;
    }
    return (double)(dtv.tv_sec) + (double)(dtv.tv_usec)*1e-6;
  }
}

void myerr(char *msg, const int n)
{
  puts(msg); exit(n);
}

double integrate(double f(double x), double a, double b, int n)
{
  int i; double h, s;

  h = (b - a) / n;
  s = 0.5 * f(a);
  for (i=1; i<n-1; i++)
    s = s + f(a+h*i);
  s = s + 0.5 * f(b);

  return h*s;
}

int fopen_m(FILE** F, const char* name, const char* mode)
{
  int i;

  *F = fopen(name, mode);

  if (*F == NULL) i = ENOENT; else i = ferror(*F);

  if (i != 0)
   switch (i)
    {
     case ENOENT:
      {
       fprintf(stderr,"Input file not found !\n");
       return EC_Bad_file_name;
      }
     case EMFILE:
      {
       fprintf(stderr,"Too many opened files !\n");
       return EC_Many_open_files;
      }
     case EACCES:
      {
       fprintf(stderr,"Access to input file is not correct !\n");
       return EC_Bad_file_access;
      }
     default:
      {
       fprintf(stderr,"Fatal i/o error !\n");
       return EC_Fatal_error;
      }
    }

  return 0;
}

int fclose_m(FILE **F)
{
  if (fclose(*F) != 0)
  {
   fprintf(stderr,"File close error !\n");
   return EC_Close_file_error;
  }
  return 0;
}
