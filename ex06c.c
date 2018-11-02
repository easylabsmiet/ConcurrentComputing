#include <stdio.h>
#include <unistd.h>
#include "mpi.h"
#define MY_TAG 333

int main(int argc, char *argv[])
{
  int np, mp, nl;
  char pname[MPI_MAX_PROCESSOR_NAME];
  MPI_Status status;
  double t0, t1, a;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&mp);
  MPI_Get_processor_name(pname,&nl);
  t0 = MPI_Wtick();

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,t0);
  sleep(1);

  if (np<2) {
    fprintf(stderr,"Too small network\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  t1 = MPI_Wtime();

  a = 0;

  if (mp==0) {
    a = a + 1.0;
    MPI_Send(&a, 1, MPI_DOUBLE, mp+1, MY_TAG, MPI_COMM_WORLD);
    MPI_Recv(&a, 1, MPI_DOUBLE, np-1, MY_TAG, MPI_COMM_WORLD, &status);
  }
  else {
    MPI_Recv(&a, 1, MPI_DOUBLE, mp-1, MY_TAG, MPI_COMM_WORLD, &status);
    a = a + 1.0;
    MPI_Send(&a, 1, MPI_DOUBLE, (mp+1) % np, MY_TAG, MPI_COMM_WORLD);
  }

  t1 = MPI_Wtime()-t1;

  sleep(1);

  fprintf(stderr,"mp=%d time=%12le res=%12le\n",mp,t1,a);

  MPI_Finalize();
  return 0;
}
