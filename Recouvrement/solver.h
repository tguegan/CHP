#ifndef __SOLVER__
#define __SOLVER__

#include "type.h"

void prod_mat_vect(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag);
double prod_scal(double *x, double *y, int Nx, int Ny);
double norme_vect(double *x, int Nx, int Ny);
void grad_conju(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag);
double erreur_max(double *x, int Nx, int Ny, double *Omg, int Num_prob, int i1, int iN);
void ecriture_visit(double *x, int Nx, int Ny, double *Omg, char *Name, int ProcId, int ProcNo, int i1, int iN, int rb);

#endif
