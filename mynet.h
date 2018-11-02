#include "mpi.h"

#define MY_TAG 777

#define MSG_MAX_CHR      65536
#define MSG_MAX_INT      16384
#define MSG_MAX_DBL       8192

typedef union tag_union_t {
  char   cdata[MSG_MAX_CHR];
  int    idata[MSG_MAX_INT];
  double ddata[MSG_MAX_DBL];
} union_t;

void mpierr(char *msg, const int n);

int MyNetInit(int* argc, char*** argv, int* np, int* mp,
              int* nl, char* pname, double* tick);

void MyRange(int np, int mp, int ia, int ib, int *i1, int *i2, int *nc);

void My2DGrid(int np, int mp, int n1, int n2,
              int *np1, int *np2, int *mp1, int *mp2);

void BndExch1D(int np, int mp, int nsl, int nrl, int nsr, int nrr,
               double *bsl, double *brl, double *bsr, double *brr);

void BndAExch1D(int mp_l, int nn_l, double* ss_l, double* rr_l,
                int mp_r, int nn_r, double* ss_r, double* rr_r);

void BndAExch2D(int mp_l, int nn_l, double* ss_l, double* rr_l,
                int mp_r, int nn_r, double* ss_r, double* rr_r,
                int mp_b, int nn_b, double* ss_b, double* rr_b,
                int mp_t, int nn_t, double* ss_t, double* rr_t);
