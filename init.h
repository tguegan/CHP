#ifndef __INIT__
#define __INIT__

#include "type.h"

void init_mat(struct matrice_diag A ,int Nx ,int Ny ,double D ,double *Omg, double Beta, double Alpha, int ProcId, int ProcNo, double dx, double dy, double dt, int Num_prob);
double f(double x, double y, double t, double *Omg, int Num_prob);
double g(double x, double y, int Num_prob);
double h(double x, double y, int Num_prob);
void vectb(double *b, int Nx, int Ny, double *Omg, int Num_prob, double T, int i1, int iN, int ProcId, int ProcNo, double dx, double dy, double D);
void charge(int *i1, int *iN, int N, int ProcId, int ProcNo, int ri, int rs);
void comptage(int x, int Nx, int *i, int *j);
void update_b(double *b, int Nx, int Ny, int N, int ProcId, int ProcNo);
  
#endif
