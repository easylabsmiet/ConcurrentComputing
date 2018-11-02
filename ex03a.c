#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

void prt_proc(int np, int mp, pid_t pid, pid_t* pids);
void prt_proc(int np, int mp, pid_t pid, pid_t* pids)
{
  int i;

  if (mp==0)
    printf("pid=%d -> I am master, mp=%d\n",pid,mp);
  else
    printf("pid=%d -> I am slave, mp=%d\n",pid,mp);

  for (i=0; i<np; i++)
    printf("pid=%d -> i=%d pids=%d\n",pid,i,pids[i]);

  return;
}

int main(int argc, char *argv[])
{
  pid_t pid, spid=0, pids[10];
  int i, np=5, mp=0;

  pid = getpid();
  printf("pid=%d -> Hello!\n",pid);

  for (i=0; i<np; i++)
    pids[i] = 0;

  i = 0; pids[i] = pid; i++;

  while (i<np){
    if (spid>0 || i==1) spid = fork();
    if (spid>0) pids[i] = spid;
    i++;
    if (spid>0) sleep(1);
  }

  if (spid == 0){
    pid = getpid();
    mp = np;
    while (pids[mp-1] == 0)
      mp--;
  }

  prt_proc(np,mp,pid,pids);

  return 0;
}
