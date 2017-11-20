#ifndef __SOLVER__
#define __SOLVER__

#include "type.h"
#include "init.h"

void prod_mat_vect(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag,int ProcId, int i1, int iN, int ProcNo, double *x_etendu, double *x_p, double *x_s);
void grad_conju(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag, int i1, int iN, int ProcId, int ProcNo, double *x_etendu, double *x_p, double *x_s, double *r, double *p, double *Ap);
double erreur_max(double *x, int Nx, int Ny, double *Omg, int Num_prob, int i1, int iN);
void ecriture_visit(double *x, int Nx, int Ny, double *Omg, char *Name, int ProcId, int ProcNo, int i1, int iN);

#endif

