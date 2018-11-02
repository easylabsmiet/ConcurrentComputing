#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mycom.h"
#include "mynet.h"
#include "myio.h"

void OutFun1D(char *name, int mode, int nx, double *x, double *f)
{
  int i; FILE *F;

  if (mode<1) i = fopen_m(&F,name,"wt"); else i = fopen_m(&F,name,"at");

  if (i!=0) return;

  if (mode<1) {
    fprintf(F,"TITLE = \"FUNC\"\n");
    fprintf(F,"VARIABLES = \"X\", \"F\"\n");
  }

  fprintf(F,"ZONE I=%d, F=POINT\n",nx);
  for (i=0; i<nx; i++) fprintf(F,"%20.13le %20.13le\n",x[i],f[i]);

  i = fclose_m(&F);

  return;
}

void OutFun1DP(char *name, int np, int mp, int nx, double *x, double *f)
{
  int n, m;
  MPI_Status status;

  if (mp==0) {
    OutFun1D(name,0,nx,x,f);
    if (np>1) {
      for (m=1; m<np; m++) {
        n = 0;
        MPI_Send(&n,1,MPI_INT,m,MY_TAG,MPI_COMM_WORLD);
        MPI_Recv(&n,1,MPI_INT,m,MY_TAG,MPI_COMM_WORLD,&status);
      }
    }
  }
  else {
    MPI_Recv(&n,1,MPI_INT,0,MY_TAG,MPI_COMM_WORLD,&status);
    OutFun1D(name,1,nx,x,f);
    MPI_Send(&n,1,MPI_INT,0,MY_TAG,MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  return;
}

void OutFun2D(char *name, int mode, int nx, int ny,
              double *x, double *y, double *f)
{
  int i,j,m; FILE *F;

  if (mode<1) i = fopen_m(&F,name,"wt"); else i = fopen_m(&F,name,"at");

  if (i!=0) return;

  if (mode<1) {
    fprintf(F,"TITLE = \"FUNC\"\n");
    fprintf(F,"VARIABLES = \"X\", \"Y\", \"F\"\n");
  }

  fprintf(F,"ZONE I=%d, J=%d, F=POINT\n",nx,ny);

  m = 0;
  for (j=0; j<ny; j++)
    for (i=0; i<nx; i++) {
      fprintf(F,"%20.13le %20.13le %20.13le\n",x[i],y[j],f[m]); m++;
    }

  i = fclose_m(&F);

  return;
}

void OutFun2DP(char *name, int np, int mp, int nx, int ny,
               double *x, double *y, double *f)
{
  int n, m;
  MPI_Status status;

  if (mp==0) {
    OutFun2D(name,0,nx,ny,x,y,f);
    if (np>1) {
      for (m=1; m<np; m++) {
        n = 0;
        MPI_Send(&n,1,MPI_INT,m,MY_TAG,MPI_COMM_WORLD);
        MPI_Recv(&n,1,MPI_INT,m,MY_TAG,MPI_COMM_WORLD,&status);
      }
    }
  }
  else {
    MPI_Recv(&n,1,MPI_INT,0,MY_TAG,MPI_COMM_WORLD,&status);
    OutFun2D(name,1,nx,ny,x,y,f);
    MPI_Send(&n,1,MPI_INT,0,MY_TAG,MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  return;
}
