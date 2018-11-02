#include <stdio.h>
#include <unistd.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
  int i, np, mp, nl;
  char pname[MPI_MAX_PROCESSOR_NAME];
  double t0, a1, a2, a3;
  MPI_Status status;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&mp);
  MPI_Get_processor_name(pname,&nl);
  t0 = MPI_Wtick();

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,t0);
  sleep(1);

  a1 = -1;
  a2 = mp * 3.14;
  a3 = -1;

  if (np>1){

    sleep(1);

    if (mp%2==0){

      if (mp+1<np){
        MPI_Sendrecv(&a2,1,MPI_DOUBLE,mp+1,0,
	             &a3,1,MPI_DOUBLE,mp+1,0,
                     MPI_COMM_WORLD,&status);
        fprintf(stderr,"%02d <-> %02d\n",mp,mp+1);
      }

      if (mp>0){
        MPI_Sendrecv(&a2,1,MPI_DOUBLE,mp-1,0,
                     &a1,1,MPI_DOUBLE,mp-1,0,
                     MPI_COMM_WORLD,&status);
        fprintf(stderr,"%02d <-> %02d\n",mp,mp-1);
      }

    }
	    
    if (mp%2==1){

      if (mp>0){
	MPI_Sendrecv(&a2,1,MPI_DOUBLE,mp-1,0,
	             &a1,1,MPI_DOUBLE,mp-1,0,
                     MPI_COMM_WORLD,&status);
	fprintf(stderr,"%02d <-> %02d\n",mp,mp-1);
      }

      if (mp+1<np){
        MPI_Sendrecv(&a2,1,MPI_DOUBLE,mp+1,0,
                     &a3,1,MPI_DOUBLE,mp+1,0,
                     MPI_COMM_WORLD,&status);
        fprintf(stderr,"%02d <-> %02d\n",mp,mp+1);
      }

    }

    sleep(1);
  }

  fprintf(stderr,"mp=%d a1=%le a2=%le a3=%le\n",mp,a1,a2,a3);

  MPI_Finalize();
  return 0;
}
