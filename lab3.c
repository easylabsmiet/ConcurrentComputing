#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h> 
#include "mycom.h" 
#include <sys/types.h> 
#include <sys/ipc.h> 
#include <sys/sem.h> 
#include <sys/shm.h> 
#include <pthread.h> 
#define SEM_ID 2007 
#define SHM_ID 2008 
#define PERMS 0600

typedef struct tag_msg_t { int n; double s; } msg_t; typedef struct tag_data_t {
  int n, nt, mt;
  double a, b, s, *sum;
} data_t;

#ifdef _SEM_SEMUN_UNDEFINED
union semun {
  int val;
  struct semid_ds *buf;
  unsigned short int *array;
  struct seminfo *__buf;
};
#endif

int np, mp, nt, mt; 
int semid; 
int shmid; 
msg_t* msg;
struct sembuf sem_loc = {0,-1,0};
struct sembuf sem_unloc = {0, 1,0}; 
union semun s_un; pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER; 
pthread_t* threads;
data_t* data; 
double a = 0; 
double b = 1; 
int ni = 1000000000; 
double sum = 0;

double f1(double x);
double f1(double x) { return 4.0/(1.0+x*x); } 

void* myjobt(void* d); 
void* myjobt(void* d) {
  int n1; 
  double a1, b1, h1;
  data_t* dd = (data_t *)d;
  n1 = dd->n / dd->nt;
  h1 = (dd->b - dd->a) / dd->nt;
  a1 = dd->a + h1 * dd->mt;
  if (dd->mt<dd->nt-1) b1 = a1 + h1; else b1 = dd->b;
  dd->s = integrate(f1,a1,b1,n1);
  pthread_mutex_lock(&mut); // lock
  printf("mt_1=%d a1=%le b1=%le n1=%d s1=%le\n",dd->mt,a1,b1,n1,dd->s);
  *dd->sum += dd->s;
  pthread_mutex_unlock(&mut); // unlock
  return 0;
}

void ThreadInit(); void ThreadInit() {
  if (!(threads = (pthread_t*) malloc(nt*sizeof(pthread_t))))
    myerr("Not enough memory",1);
  if (!(data = (data_t*) malloc(nt*sizeof(data_t))))
    myerr("Not enough memory",1);
  int i;
  for (i=0; i<nt; i++){
    (data+i)->a=a;
    (data+i)->b=b/(nt+np)*nt;
    (data+i)->n=ni/(nt+np)*nt;
    (data+i)->nt=nt;
    (data+i)->mt=i;
    (data+i)->sum = &sum;
    if (pthread_create(threads+i,0,myjobt,(void*)(data+i)))
      myerr("Can not create thread",2);
  }
  return;
}

void ThreadDone(); void ThreadDone() {
  int i;
  for (i=0; i<nt; i++)
    if (pthread_join(threads[i],0))
      myerr("Can not wait thread",3);
  free(data);
  free(threads);
  return;
}

void myjobp(); void myjobp() {
  int n1; 
  double a1, b1, h1, s1;
  if (mp == 0)
     ThreadInit();
  n1 = ni / (np+nt);
  h1 = (b - a) / (np + nt);
  a1 = a + h1 * (nt + mp);
  if (mp<np-1)
	  b1 = a1 + h1;
  else b1 = b;
  s1 = integrate(f1,a1,b1,n1);
  printf("mp=%d a1=%le b1=%le n1=%d s1=%le\n",mp,a1,b1,n1,s1);
  while(semop(semid,&sem_loc,1)<0); // wait + lock
 	 msg->n++; msg->s += s1;
  while(semop(semid,&sem_unloc,1)<0); // wait + unlock

  if (mp == 0){
    ThreadDone();
    double res = sum;
    n1 = 0;
    while(n1<np){
      while(semop(semid,&sem_loc,1)<0); // wait + lock
      n1 = msg->n;
      sum = msg->s;
      while(semop(semid,&sem_unloc,1)<0); // wait + unlock
    }
    sum+=res;
  }
  return;
}

void NetInit(int np, int* mp); void NetInit(int np, int* mp) {
  int i; 
  pid_t spid = 0;
  if (np>1){
    *mp=1; spid = fork();
    if (spid > 0 && np > 2)
      for (i=2;i<np;i++)
        if (spid > 0){ *mp=i; spid = fork();}
    if (spid > 0) *mp = 0;
  }
  else *mp = 0;
  return;
}

void ShrMemInit(); void ShrMemInit() {
  if (mp == 0) {
    if ((shmid = shmget(SHM_ID,sizeof(msg_t),PERMS | IPC_CREAT)) < 0)
      if ((shmid = shmget(SHM_ID,sizeof(msg_t),PERMS)) < 0)
        myerr("Can not find shared memory segment",1);
    if ((msg = (msg_t*) shmat(shmid, 0, 0)) == NULL)
      myerr("Can not attach shared memory segment",2);
    msg->n = 0; msg->s = 0; // initialization of shared memory
    if ((semid = semget(SEM_ID, 1, PERMS | IPC_CREAT )) < 0)
      if ((semid = semget(SEM_ID, 1, PERMS)) < 0)
        myerr("Can not find semaphore",3);
    s_un.val = 1; semctl(semid,0,SETVAL,s_un); // unlock
  }
  else {
    while((shmid = shmget(SHM_ID,sizeof(msg_t),PERMS)) < 0);
    if ((msg = (msg_t*) shmat(shmid, 0, 0)) == NULL)
      myerr("Can not attach shared memory segment",4);
    while((semid = semget(SEM_ID, 1, PERMS)) < 0);
  }
  return;
}

void ShrMemDone(); void ShrMemDone() {
  if (mp == 0) {
    if (semctl(semid, 0, IPC_RMID, (struct semid_ds *) 0) < 0)
      myerr("Can not remove semaphore",5);
    if (shmdt(msg)<0)
      myerr("Can not dettach shared memory segment",6);
    if (shmctl(shmid, IPC_RMID, (struct shmid_ds *) 0) < 0)
      myerr("Can not remove shared memory segment",7);
  }
  else {
    if (shmdt(msg)<0)
      myerr("Can not dettach shared memory segment",8);
    exit(0);
  }
  return;
}

int main(int argc, char *argv[]) {
  double t;
  if (argc<3){
    printf("Usage: %s <process number> <thread number>\n",argv[0]);
    return 1;
  }
  sscanf(argv[1],"%d",&np);
  sscanf(argv[2],"%d",&nt);
  
  if (np<1) np = 1; mp = 0;
  if (nt<1) nt = 1; mt = 0;
  
  t = mytime(0);
  if (np == 1 && nt == 0)
    sum = integrate(f1,a,b,ni);
  else
  {
    NetInit(np,&mp);
    if (np>1) ShrMemInit();
    myjobp();
    if (np>1) ShrMemDone();
  }
  t = mytime(1);
  printf("time=%lf sum=%22.15le\n",t,sum);
  return 0;
}
