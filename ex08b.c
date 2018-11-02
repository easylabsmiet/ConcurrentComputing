#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mycom.h"
#include "mynet.h"

double f1(double x);
double f1(double x) {
  double e=exp(x);
  return 0.5*(e+1.0/e);
}

double f2(double x);
double f2(double x) {
  double e=exp(x);
  return 0.5*(e-1.0/e);
}

double f(double x, double y);
double f(double x, double y) {
  return f1(x)*f2(y);
}

int np, mp, nl, ier, lp;
int np1, np2, mp1, mp2;
char pname[MPI_MAX_PROCESSOR_NAME];
char sname[10] = "ex08b.p00";
MPI_Status Sta[100];
MPI_Request Req[100];
union_t buf;
double tick, t1, t2, t3;

FILE *Fi = NULL;
FILE *Fo = NULL;
int nx, ny, fl=1;
double xa, xb, ya, yb;

int main(int argc, char *argv[])
{
  int i, i1, i2, j, j1, j2;
  double s, p, x, y, hx, hy, hxy;

  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  sprintf(sname+7,"%02d",mp);
  ier = fopen_m(&Fo,sname,"wt");
  if (ier!=0) mpierr("Protocol file not opened",1);

  if (mp==0) {
    ier = fopen_m(&Fi,"ex08b.d","rt");
    if (ier!=0) mpierr("Data file not opened",2);
    i = fscanf(Fi,"xa=%le\n",&xa);
    i = fscanf(Fi,"xb=%le\n",&xb);
    i = fscanf(Fi,"ya=%le\n",&ya);
    i = fscanf(Fi,"yb=%le\n",&yb);
    i = fscanf(Fi,"nx=%d\n",&nx);
    i = fscanf(Fi,"ny=%d\n",&ny);
    i = fscanf(Fi,"fl=%d\n",&fl);
    fclose_m(&Fi);
    fprintf(stderr,"read is OK\n");
  }

  MPI_Bcast(&fl,1,MPI_INT,0,MPI_COMM_WORLD);

  if (np>1) {
    if (fl<1){ // variant 1:
      if (mp==0) {
        buf.ddata[0] = xa; buf.ddata[1] = xb;
        buf.ddata[2] = ya; buf.ddata[3] = yb;
        buf.idata[8] = nx; buf.idata[9] = ny;
        for (i=1; i<np; i++) {
          MPI_Isend(buf.ddata,5,MPI_DOUBLE,i,MY_TAG,MPI_COMM_WORLD,Req+i);
        }
        MPI_Waitall(np-1,Req+1,Sta+1);
      }
      else {
        MPI_Recv(buf.ddata,5,MPI_DOUBLE,0,MY_TAG,MPI_COMM_WORLD,Sta+0);
        xa = buf.ddata[0]; xb = buf.ddata[1];
        ya = buf.ddata[2]; yb = buf.ddata[3];
        nx = buf.idata[8]; ny = buf.idata[9];
      }
    }
    else { // variant 2:
      if (mp==0) {
        buf.ddata[0] = xa; buf.ddata[1] = xb;
        buf.ddata[2] = ya; buf.ddata[3] = yb;
        buf.idata[8] = nx; buf.idata[9] = ny;
      }
      MPI_Bcast(buf.ddata,5,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if (mp>0) {
        xa = buf.ddata[0]; xb = buf.ddata[1];
        ya = buf.ddata[2]; yb = buf.ddata[3];
        nx = buf.idata[8]; ny = buf.idata[9];
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

// Next code:

  if (np==1) { np1 = 1; np2 = 1; }
  else {
    s = sqrt(((double)np)) * ((double)nx) / ((double)ny);
    np1 = floor(s); if (s>0.5+((double)np1)) np1++;
    np2 = np / np1;
    if (np1*np2!=np) {
      if (nx>ny) {np1 = np; np2 = 1;} else {np1 = 1; np2 = np;}
    }
  }

  mp2 = mp / np1;
  mp1 = mp % np1;

  if (mp==0) fprintf(stderr,"Grid=%dx%d\n",np1,np2);

  fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  fprintf(Fo,"Grid=%dx%d coord=(%d,%d)\n",np1,np2,mp1,mp2);
  fprintf(Fo,"xa=%le xb=%le ya=%le yb=%le nx=%d ny=%d fl=%d\n",xa,xb,ya,yb,nx,ny,fl);

  t1 = MPI_Wtime();

  hx = (xb-xa)/nx; hy = (yb-ya)/ny; hxy = hx*hy;

  if (np1==1) { i1 = 0; i2 = nx-1; }
  else {
    i1 = mp1 * (nx/np1);
    if (mp1<np1-1) i2 = i1+(nx/np1)-1; else i2 = nx-1;
  }

  if (np2==1) { j1 = 0; j2 = ny-1; }
  else {
    j1 = mp2 * (ny/np2); 
    if (mp2<np2-1) j2 = j1+(ny/np2)-1; else j2 = ny-1;
  }

  s = 0;

  for (j=j1; j<=j2; j++) {
    y = ya + (j*1.0+.5)*hy;
    for (i=i1; i<=i2; i++) {
      x = xa + (i*1.0+.5)*hx;
      s = s + hxy * f(x,y);
    }
  }

  t2 = MPI_Wtime();

  t1 = t2 - t1;

  if (np==1)
    t2 = 0;
  else {
    p = s;
    MPI_Reduce(&p, &s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    t2 = MPI_Wtime()-t2;
  }

  t3 = t1 + t2;

  if (mp==0) fprintf(stderr,"t1=%le t2=%le t3=%le int=%le\n",t1,t2,t3,s);

  fprintf(Fo,"i1=%d i2=%d j1=%d j2=%d\n",i1,i2,j1,j2);
  fprintf(Fo,"t1=%le t2=%le t3=%le int=%le\n",t1,t2,t3,s);

  ier = fclose_m(&Fo);

  MPI_Finalize();
  return 0;
}
