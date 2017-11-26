#include "schwarz.h"
#include "solver.h"

void schwarz(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag, int i1, int iN, int ProcId, int ProcNo, double *Omg, int ri, int rs, double dt)
{
  int i, N, cpt = 0;
  double eps, normex = 0, normexn = 0, *x_s, *x_p, *xcom, *btmp, dy;
  MPI_Status Status;

  eps = 1e-6;
  N = (iN - i1 + 1) * Nx;
  dy = (Omg[3] - Omg[2]) / (Ny + 1.0);
  x_s = (double*) calloc(Nx, sizeof(double));
  x_p = (double*) calloc(Nx, sizeof(double));
  xcom = (double*) calloc(N, sizeof(double));
  btmp = (double*) calloc(N, sizeof(double));
  
  grad_conju(A, x, b, Nx, iN - i1 + 1, Nb_diag);
  
  normex = norme_vect(x, Nx, iN - i1 + 1);

  while (fabs(normex - normexn) > eps)
    {
      cpt++;
      normex = normexn;
      for (i = 0; i < Nx; i++)
	{
	  x_s[i] = x[i + (iN - i1 - rs - ri) * Nx];
	  x_p[i] = x[i + (rs + ri) * Nx];
	}
      
      if (ProcId == 0)
	MPI_Sendrecv(x_s, Nx, MPI_DOUBLE, ProcId + 1, 0, x_p, Nx, MPI_DOUBLE, ProcId + 1, 0, MPI_COMM_WORLD, &Status);
      else if (ProcId < ProcNo - 1)
	{
	  MPI_Sendrecv(x_p, Nx, MPI_DOUBLE, ProcId - 1, 0, xcom, Nx, MPI_DOUBLE, ProcId - 1, 0, MPI_COMM_WORLD, &Status);
	  MPI_Sendrecv(x_s, Nx, MPI_DOUBLE, ProcId + 1, 0, x_p, Nx, MPI_DOUBLE, ProcId + 1, 0, MPI_COMM_WORLD, &Status);
	}
      else 
	MPI_Sendrecv(x_p, Nx, MPI_DOUBLE, ProcId - 1, 0, xcom, Nx, MPI_DOUBLE, ProcId - 1, 0, MPI_COMM_WORLD, &Status);
      
      if (ProcId < ProcNo - 1)
	for (i = 0; i < Nx; i++)
	  xcom[N - Nx + i] = x_p[i];

      for (i = 0; i < N; i++)
	btmp[i] = b[i] + dt * xcom[i]/(dy * dy);

      grad_conju(A, x, btmp, Nx, iN - i1 + 1, Nb_diag);
      normexn = prod_scal(x, x, Nx, iN - i1 + 1);
      MPI_Allreduce(&normexn, &normexn, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      normexn = sqrt(normexn);

      if (ProcId == 0)
	printf("\nCompteur de Schwarz : %d,  Résidu  : %.10lf\n", cpt, fabs(normex - normexn));
    }

  free(btmp);
  free(x_s);
  free(x_p);
  free(xcom);
}
