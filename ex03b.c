#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/msg.h>

#define MSG_ID   7777
#define MSG_PERM 00666
#define LBUF     100

typedef struct tag_msg_t
{
 int type; pid_t buf[LBUF];
}
msg_t;

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
  pid_t pid, spid=0, pids[LBUF];
  int i, np=5, mp=0;
  int msgid; msg_t msg;

  pid = getpid();
  printf("pid=%d -> Hello!\n",pid);

  for (i=0; i<np; i++)
    pids[i] = 0;

  i = 0; pids[i] = pid; i++;

  while (i<np){
    if (spid>0 || i==1) spid = fork();
    if (spid>0) pids[i] = spid;
    i++;
  }

  if (spid == 0){
    pid = getpid();
    mp = np;
    while (pids[mp-1] == 0)
      mp--;
  }

  if (mp == 0){
    msgid = msgget(MSG_ID, MSG_PERM | IPC_CREAT);
    msg.type = 888;
    for(i=0;i<np;i++)
      msg.buf[i] = pids[i];
    for(i=1;i<np;i++)
      msgsnd(msgid, &msg, sizeof(msg_t), 0);
    sleep(np+1);
    msgctl(msgid, IPC_RMID, (struct msqid_ds *) 0); 
    prt_proc(np,mp,pid,pids);
  }
  else{
    while((msgid = msgget(MSG_ID, 0)) < 0);
    msgrcv(msgid, &msg, sizeof(msg_t), 0, 0);
    for(i=0;i<np;i++)
      pids[i] = msg.buf[i];
    sleep(mp);
    prt_proc(np,mp,pid,pids);
  }

  return 0;
}
