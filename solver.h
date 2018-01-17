#ifndef __SOLVER__
#define __SOLVER__

#include "type.h"

void precondi(struct matrice_diag A, double *b, int Nx, int Ny, int Nb_diag, int ProcId, int ProcNo);
double erreur_max(double *x, int Nx, int Ny, double *Omg, int Num_prob, int i1, int iN, double D);
void ecriture_visit(double *x, int Nx, int Ny, double *Omg, char *Name, int ProcId, int ProcNo, int i1, int iN, int rb);
void bigradstab(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag, double *r, double *rtilde, double *p, double *v, double *s, double *t);

#endif
