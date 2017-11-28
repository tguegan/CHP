#ifndef __INIT__
#define __INIT__

#include "type.h"

void init_mat(struct matrice_diag A ,int Nx ,int Ny ,double D ,double *Omg, double dx, double dy);
double f(double x, double y, double t, double *Omg, int Num_prob);
double g(double x, double y, int Num_prob);
double h(double x, double y, int Num_prob);
void vectb(double *b, int Nx, int Ny, double *Omg, int Num_prob, double T, int j1, int i1, int ProcId, int ProcNo, double dx, double dy, int Nbdomv, int Nbdomh);
void charge(int *i1, int *iN, int N, int ProcId, int ProcNo, int ri, int rs);
void comptage(int x, int Nx, int *i, int *j);
  
#endif
