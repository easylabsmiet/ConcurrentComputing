#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mynet.h"
#include "myrand.h"

int np, mp, nl, ier, lp;
char pname[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;
double tick, t1;

int main(int argc, char *argv[])
{
  int ns, ncm, nc, nn, i1, i2, i, ii;
  int* ind = NULL;
  double* a1 = NULL;
  double* a2 = NULL;
  double amin, amax;

  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  if (argc<2) ns = 10000;
  else        sscanf(argv[1],"%d",&ns);

  if (ns<np) ns = 10000;

  ind = (int*)(malloc(sizeof(int)*3*np));

  myindexes(np,ns,ind);

  i1 = ind[3*mp];
  i2 = ind[3*mp+1];
  nc = ind[3*mp+2];

  fprintf(stderr,"mp=%d ns=%d i1=%d i2=%d nc=%d\n",mp,ns,i1,i2,nc);

  t1 = MPI_Wtime();

  if (np==1) {
    a1 = (double*)(malloc(sizeof(double)*nc));

    for (i=0; i<nc; i++)
      a1[i] = myrand();

    qsort(a1, nc, sizeof(double), mycomp);
  }
  else {
    ncm = ns / np + 1;

    a1 = (double*)(malloc(sizeof(double)*ncm));
    a2 = (double*)(malloc(sizeof(double)*ncm));

    if (mp==0) {
      for (i=0; i<nc; i++)
        a1[i] = myrand();

      for (ii=1; ii<np; ii++) {
        nn = ind[3*ii+2];
        for (i=0; i<nn; i++)
          a2[i] = myrand();
        MPI_Send(a2,nn,MPI_DOUBLE,ii,MY_TAG,MPI_COMM_WORLD);
      }
    }
    else
      MPI_Recv(a1,nc,MPI_DOUBLE,0,MY_TAG,MPI_COMM_WORLD,&status);

    qsort(a1,nc,sizeof(double),mycomp);

    for (ii=0; ii<1; ii++){
      if (mp%2==0){
        if (mp+1<np){
	  MPI_Sendrecv(a1,ncm,MPI_DOUBLE,mp+1,0,
	               a2,ncm,MPI_DOUBLE,mp+1,0,
                       MPI_COMM_WORLD,&status);
          //..
        }
        if (mp>0){
          MPI_Sendrecv(a1,ncm,MPI_DOUBLE,mp-1,0,
                       a2,ncm,MPI_DOUBLE,mp-1,0,
                       MPI_COMM_WORLD,&status);
          //..
        }
      }
      else{
        if (mp>0){
          MPI_Sendrecv(a1,ncm,MPI_DOUBLE,mp-1,0,
		       a2,ncm,MPI_DOUBLE,mp-1,0,
		       MPI_COMM_WORLD,&status);
	  //..
        }
        if (mp+1<np){
          MPI_Sendrecv(a1,ncm,MPI_DOUBLE,mp+1,0,
                       a2,ncm,MPI_DOUBLE,mp+1,0,
                       MPI_COMM_WORLD,&status);
          //..
        }
      }
    }
  }

  amin = a1[0];
  amax = a1[nc-1];

  t1 = MPI_Wtime() - t1;

  fprintf(stderr,"mp=%d amin=%le amax=%le t1=%le\n",mp,amin,amax,t1);

  MPI_Finalize();
  return 0;
}
