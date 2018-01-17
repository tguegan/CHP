#ifndef __SCHWARZ__
#define __SCHWARZ__

#include "type.h"

void schwarz(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag, int ProcId, int ProcNo, double *Omg,int ri, int rs, int Init, double Beta, double Alpha, double dy);
void critere_schwarz(double *x, double *x_np, double *x_ns, double *xhaut, double *xbas, int Nx, int Ny, int ri, int rs, int ProcId, int ProcNo, double *normexbord);

#endif
