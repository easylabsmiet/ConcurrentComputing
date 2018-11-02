int prog_right(int n, double* a, double* b, double* c,
                      double* f, double* al, double* y);

int prog_rightm(int n, double* a, double* al, double* y);

int prog_rightp(int np, int mp, int nc,
                double *aa, double *bb, double *cc, double *ff,
                double *al, double *y1, double *y2, double *y3,
                double *y4, double *dd, double *ee);

int prog_rightpm(int np, int mp, int nc, int ip,
                 double *aa, double *bb, double *cc, double *ff,
                 double *al, double *y1, double *y2, double *y3, double *y4);

int prog_rightpn(int np, int mp, MPI_Comm cm, int nc, int ip,
                 double *aa, double *bb, double *cc, double *ff,
                 double *al, double *y1, double *y2, double *y3, double *y4);
