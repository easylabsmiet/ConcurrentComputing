#include "stdio.h"
#include <unistd.h>
#include <math.h>
#include "mynet.h"

int np, mp, nl, ier, lp;
char pname[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;
double tick, t1, t2, t3;

int main(int argc, char *argv[])
{
  MPI_Group gr, gr1, gr2;
  MPI_Comm cm1, cm2, cm3;
  int i, np1, np2, mp1, mp2, ranks[128];
  double s, p;

  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  if (np<2) mpierr("Too small processes for communications",1);

  MPI_Comm_group(MPI_COMM_WORLD, &gr);

  np1 = np/2; np2 = np - np1;

  for (i=0; i<np1; i++)
    ranks[i] = i;

  MPI_Group_incl(gr, np1, ranks, &gr1);
  MPI_Group_excl(gr, np1, ranks, &gr2);

  MPI_Group_rank(gr1, &mp1);
  MPI_Group_rank(gr2, &mp2);

  MPI_Comm_create(MPI_COMM_WORLD, gr1, &cm1);
  MPI_Comm_create(MPI_COMM_WORLD, gr2, &cm2); 

  if (cm1>0) cm3 = cm1; else if (cm2>0) cm3 = cm2; else cm3 = 0;

  s = 5.0 * (mp+1); p = 0.0;

  MPI_Reduce(&s, &p, 1, MPI_DOUBLE, MPI_SUM, 0, cm3);

  fprintf(stderr,"mp=%d np1=%d np2=%d mp1=%d mp2=%d cm1=%d cm2=%d cm3=%d s=%le p=%le\n",
    mp,np1,np2,mp1,mp2,cm1,cm2,cm3,s,p);

  MPI_Group_free(&gr);
  MPI_Group_free(&gr1);
  MPI_Group_free(&gr2);
  MPI_Finalize();

  return 0;
}
