#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mycom.h"

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>

#define MSG_ID   7777
#define MSG_PERM 0600
#define LBUF     10

typedef struct tag_msg_t
{ int n; double buf[LBUF]; } msg_t;

int np, mp;
int msgid;
msg_t msg; 

double a = 0;
double b = 1;
int ni = 1000000000;
double sum = 0;

double f1(double x);
double f1(double x) { return 4.0/(1.0+x*x); }

double myjobp(int mp);
double myjobp(int mp)
{
  int n1; double a1, b1, h1, s1;

  n1 = ni / np;
  h1 = (b - a) / np;
  a1 = a + h1 * mp;
  if (mp<np-1) b1 = a1 + h1; else b1 = b;

  s1 = integrate(f1,a1,b1,n1);

  printf("mp=%d a1=%le b1=%le n1=%d s1=%le\n",mp,a1,b1,n1,s1);

  return s1;
}

void NetInit(int np, int* mp);
void NetInit(int np, int* mp)
{
  int i; pid_t spid = 0;

  if (np>1){
    i = 1;
    while (i<np){
      if (spid>0 || i==1){ *mp=i; spid = fork();}
      if (spid==0) return;
      i++;
    }
  }

  *mp = 0;
  return;
}

int main(int argc, char *argv[])
{
  int i; double t;

  if (argc<2){
    printf("Usage: %s <process number>\n",argv[0]);
    return 1;
  }

  sscanf(argv[1],"%d",&np); mp = 0;

  t = mytime(0);

  if (np<2)
    sum = integrate(f1,a,b,ni);
  else {
    NetInit(np,&mp);

    if (mp == 0)
      msgid = msgget(MSG_ID, MSG_PERM | IPC_CREAT);
    else
      while((msgid = msgget(MSG_ID, MSG_PERM))<0);

    sum = myjobp(mp);

    if (mp>0){
      msg.n = 1; msg.buf[0] = sum;
      msgsnd(msgid, &msg, sizeof(msg_t), 0);
    }
    else{
      for(i=1;i<np;i++){
        msgrcv(msgid,&msg,sizeof(msg_t),0,0);
        sum = sum + msg.buf[0];
      }
      msgctl(msgid,IPC_RMID,(struct msqid_ds *)0);
    }
  }

  t = mytime(1);

  if (mp == 0)
    printf("time=%lf sum=%22.15le\n",t,sum);

  return 0;
}
