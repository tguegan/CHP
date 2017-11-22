#ifndef __SCHWARZ__
#define __SCHWARZ__

#include "type.h"

void schwarz(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag, int i1, int iN, int ProcId, int ProcNo, double *Omg);

#endif
