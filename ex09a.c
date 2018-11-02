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
  MPI_Group group, group1, group2;
  int np1, np2, mp1, mp2, mp3, ranks[128];
  int i, a[4], b[4];

  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  if (np<2) mpierr("Too small processes for communications",1);

  MPI_Comm_group(MPI_COMM_WORLD, &group);

  np1 = np/2; np2 = np - np1;

  for (i=0; i<np1; i++)
    ranks[i] = i;

  MPI_Group_incl(group, np1, ranks, &group1);
  MPI_Group_excl(group, np1, ranks, &group2);

  MPI_Group_rank(group1, &mp1);
  MPI_Group_rank(group2, &mp2);

  if (mp1 == MPI_UNDEFINED) {
    if (mp2<np1)
      MPI_Group_translate_ranks(group1, 1, &mp2, group, &mp3);
    else
      mp3 = MPI_UNDEFINED;
  }
  else
    MPI_Group_translate_ranks(group2, 1, &mp1, group, &mp3);

  a[0] = mp;
  a[1] = mp1;
  a[2] = mp2;
  a[3] = mp3;

  if (mp3 != MPI_UNDEFINED) {
    MPI_Sendrecv(a, 4, MPI_INT, mp3, MY_TAG,
                 b, 4, MPI_INT, mp3, MY_TAG,
                 MPI_COMM_WORLD, &status);
    fprintf(stderr,"%02d <-> %02d\n",mp,mp3);
  }

  for (i=0; i<4; i++)
    fprintf(stderr,"mp=%d i=%d a=%d b=%d\n",mp,i,a[i],b[i]);

  MPI_Group_free(&group);
  MPI_Group_free(&group1);
  MPI_Group_free(&group2);
  MPI_Finalize();

  return 0;
}
