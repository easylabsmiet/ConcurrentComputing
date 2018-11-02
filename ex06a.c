#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
  int i, ii, np, mp, nl, nc=10000;
  char pname[MPI_MAX_PROCESSOR_NAME];
  double t0, t1, a;

  MPI_Initialized(&ii);
  fprintf(stderr,"Before MPI_Init ii=%d\n",ii);

  MPI_Init(&argc,&argv);

  MPI_Initialized(&ii);
  fprintf(stderr,"Between MPI_Init & MPI_Finalize ii=%d\n",ii);

  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&mp);
  MPI_Get_processor_name(pname,&nl);

  t0 = MPI_Wtick();

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",
          np,mp,pname,t0);

  t1 = MPI_Wtime();

  a = 0;
  for (i=0; i<nc; i++)
   {a = a + 1.23*i; if (a>1000) a = a / 1000;}

  t1 = MPI_Wtime()-t1;

  fprintf(stderr,"mp=%d, time=%12le res=%12le\n",mp,t1,a);

  MPI_Finalize();

  MPI_Initialized(&ii);
  fprintf(stderr,"After MPI_Finalize ii=%d\n",ii);

  return 0;
}
