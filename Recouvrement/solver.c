#include "solver.h"
#include "init.h"

void prod_mat_vect(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag)
{
  int i, k, l, cpt = 0;

  for (i = 0; i < Nx * Ny; i++)
    {
      b[i] = 0;
      for (k = 0; k < Nb_diag; k++)
	{
	  l = i + A.distance[k];
	  if ((l < Nx * Ny) && (l >= 0))
	    b[i] = b[i] + A.valeur[i][k] * x[l];
	  cpt++;
	}
    }
}

double prod_scal(double *x, double *y, int Nx, int Ny)
{
  double somme = 0;
  int i;

  for (i = 0; i < Nx * Ny; i++)
    somme += x[i] * y[i];

  return somme;
}

double norme_vect(double *x, int Nx, int Ny)
{
  return sqrt(prod_scal(x, x, Nx, Ny));
}

void grad_conju(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag)
{
  double *r, *Ap, *p, alp, beta, tmpa, tmpb, epsi, normer = 1, normeb;
  int i, N, cpt = 0;

  N = Nx * Ny;
  epsi = 1e-8;

  r = (double*) calloc(N, sizeof(double));
  p = (double*) calloc(N, sizeof(double));
  Ap = (double*) calloc(N, sizeof(double));

  prod_mat_vect(A, x, r, Nx, Ny, Nb_diag);
  for (i = 0; i < N; i++)
    {
      r[i] = b[i] - r[i];
      p[i] = r[i];
    }

  normeb = norme_vect(b, Nx, Ny);
  
  while ((normer > epsi) && (cpt < 20000))
    {
      prod_mat_vect(A, p, Ap, Nx, Ny, Nb_diag);

      tmpa = prod_scal(r, p, Nx, Ny);
      tmpb = prod_scal(Ap, p, Nx, Ny);

      alp = tmpa / tmpb;

      for (i = 0; i < N; i++)
	{
	  x[i] += alp * p[i];
	  Ap[i] = r[i] - alp * Ap[i];
	}
      tmpb = prod_scal(Ap, Ap, Nx, Ny);

      beta = tmpb / tmpa;

      for (i = 0; i < N; i++)
	{
	p[i] = Ap[i] + beta * p[i];
	r[i] = Ap[i];
	}
      
      normer = sqrt(tmpb) / normeb;
      cpt ++;
    }
  
  free(r);
  free(p);
  free(Ap);
}

void ecriture_visit(double *x, int Nx, int Ny, double *Omg, char *Name, int ProcId, int ProcNo, int i1, int iN, int rb)
{
  FILE *OutFile;
  double dx, dy;
  int i, k, m, n;

  dx = (Omg[1] - Omg[0]) / (Nx + 1.0);
  dy = (Omg[3] - Omg[2]) / (Ny + 1.0);

  for (k = 0; k < ProcNo; k++)
    {
      if (ProcId == k && k == 0)
	{
	  OutFile = fopen(Name, "w");
	  fprintf(OutFile,"TITLE = \"PROBLEME SOLUTION\" \n");
	  fprintf(OutFile,"VARIABLES = \"X\", \"Y\", \"U\", \"Dom\"\n");
	  fprintf(OutFile,"ZONE T=\"SQUARE\", I=%d, J=%d, F=POINT \n", Nx, Ny);
	  for (i = i1; i < iN + 1; i++)
	    {
	      comptage(i, Nx, &m, &n);
	      fprintf(OutFile,"%lf %lf %lf %d\n", (m + 1) * dx, (n + 1) * dy, x[i - i1], k);
	    }
	  fclose(OutFile);
	}
      else if (ProcId == k)
	{
	  OutFile = fopen(Name, "a");
	  for (i = i1; i < iN + 1; i++)
	    {
	      comptage(i, Nx, &m, &n);
	      fprintf(OutFile,"%lf %lf %lf %d\n", (m + 1) * dx, (n + 1) * dy, x[i + rb*Nx - i1], k);
	    }
	  fclose(OutFile);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
}

double erreur_max(double *x, int Nx, int Ny, double *Omg, int Num_prob , int i1 , int iN) 
{
  double tmp, max = 0, dx, dy;
  int i, j;

  dx = (Omg[1] - Omg[0]) / (Nx + 1.0);
  dy = (Omg[3] - Omg[2]) / (Ny + 1.0);
  
  for (i = 0; i < Nx; i++)
    {
      for (j = 0; j < iN - i1 + 1; j++)
	{
	  switch (Num_prob)
	    {
	    case 1:
	      tmp = fabs(x[i + j * Nx] - (- (i +1) * dx * (j + i1 + 1) * dy * (Omg[1] - (i + 1) * dx) * (Omg[3] - (j + i1 + 1) * dy)));
	      if (tmp > max)
		max = tmp;
	      break;
	    case 2:
	      tmp = fabs(x[i + j * Nx] - (sin((i + 1) * dx) + cos((j + i1 + 1) * dy)));
	      if (tmp > max)
		max = tmp;
	      break;
	    }
	}
    }

  MPI_Allreduce(&max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return max;
}
