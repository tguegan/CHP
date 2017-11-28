#include "init.h"

void init_mat(struct matrice_diag A ,int Nx ,int Ny ,double D ,double *Omg, double dx, double dy)
{
  double val[5];
  int Nb_diag = 5, i, j, N, k;

  N = Nx * Ny; /*Nombre total d'inconnues*/

  /*Distance des sous-diag  // diag principale */
  A.distance[0] = - Nx;
  A.distance[1] = - 1;
  A.distance[2] = 0;
  A.distance[3] = 1;
  A.distance[4] = Nx;

  /*Valeur que prend la matrice a chaque ligne*/
  val[0] = - D / (dy * dy);
  val[1] = - D / (dx * dx);
  val[2] = (2 * D) / (dx * dx) + (2 * D) / (dy * dy);
  val[3] = - D / (dx * dx);
  val[4] = - D / (dy * dy);

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < Nb_diag; j++)
  	{
  	  k = i + A.distance[j];
  	  if((k > -1) && (k < N))
  	    A.valeur[i][j] = val[j];/*Affectation des valeurs a la matrice*/
  	}
    }

  for (i = N - Nx - 1; i > Nx - 2; i -= Nx)
    A.valeur[i][3] = 0; /* 0 correspondant aux intersections des blocs de matrices*/

  for (i = Nx; i < N - Nx + 1; i+=Nx)
    A.valeur[i][1] = 0; /*Correspondant aux intersections des blocs de matrices*/
}

double f(double x, double y, double t, double *Omg, int Num_prob)
{
  switch (Num_prob)
    {
    case 1:
      return 2 * ((x * x) - x * Omg[1] + (y * y) - y * Omg[3]);
    case 2:
      return sin(x) + cos(y);
    case 3:
      return exp(- pow(x - Omg[1] / 2, 2.0)) * exp(- pow(y - Omg[3] / 2, 2.0)) * cos(M_PI * t);
    default:
      exit(0);
    }
}

double g(double x, double y, int Num_prob)
{
  switch (Num_prob)
    {
    case 1:
    case 3:
      return 0.0;
    case 2:
      return sin(x) + cos(y);
    default:
      exit(0);
    }
}

double h(double x, double y, int Num_prob)
{
  switch (Num_prob)
    {
    case 1:
      return 0.0;
    case 2:
      return sin(x) + cos(y);
    case 3:
      return 1.0;
    default:
      exit(0);
    }
}

void vectb(double *b, int Nx, int Ny, double *Omg, int Num_prob, double T, int j1, int i1, int ProcId, int ProcNo, double dx, double dy, int Nbdomv, int Nbdomh)
{
  int i, j;

  for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++)
      b[i + j * Nx] = f((i + j1 + 1) * dx, (j + i1 + 1) * dy , T , Omg, Num_prob);

  if (Nbdomv == 1)
    {
      if (ProcId == 0)
	for (i = 0; i < Nx; i++)
	  b[i] += g((i + j1 + 1) * dx, Omg[2], Num_prob) / (dy * dy);
      else if (ProcId == ProcNo - 1)
	for (i = 0; i < Nx; i++)
	  b[i + (Ny - 1) * Nx] += g((i + j1 + 1) * dx, Omg[3], Num_prob) / (dy * dy); 
      for(i = 0; i < Ny; i++)
	{
	  b[i * Nx] += h(Omg[0], (i + i1 + 1) * dy, Num_prob) / (dx * dx);
	  b[Nx - 1 + i * Nx] += h(Omg[1], (i + i1 + 1) * dy, Num_prob) / (dx * dx);
	}
    }

  if (Nbdomh == 1)
    {
      for (i = 0; i < Nx; i++)
	{
	  b[i] += g((i + j1 + 1) * dx, Omg[2], Num_prob) / (dy * dy);
	  b[i + (Ny - 1) * Nx] += g((i + j1 + 1) * dx, Omg[3], Num_prob) / (dy * dy);
	}
      if (ProcId == 0)
	for (i = 0; i < Ny; i++)
	  b[i * Nx] += h(Omg[0], (i + i1 + 1) * dy, Num_prob) / (dx * dx);
      else if (ProcId == ProcNo - 1)
	for (i = 0; i < Ny; i++)
	  b[i + (Ny - 1) * Nx] += g((i + j1 + 1) * dx, Omg[3], Num_prob) / (dy * dy);
    }
}

void charge(int *i1, int *iN, int N, int ProcId, int ProcNo,int ri,int rs)
{
  if (ProcId > 0 && ProcId < ProcNo - 1)
    {
      *i1 = round(ProcId * ((double) N) / ProcNo ) - ri;
      *iN = round((ProcId + 1) * (((double) N) / ProcNo)) - 1 + rs;
    }
  else if (ProcId == 0 && ProcNo != 1)
    {
      *i1 = round(ProcId * ((double) N) / ProcNo );
      *iN = round((ProcId + 1) * (((double) N) / ProcNo)) - 1 + rs;
    }
  else if (ProcId == 0)
    {
      *i1 = round(ProcId * ((double) N) / ProcNo );
      *iN = round((ProcId + 1) * (((double) N) / ProcNo)) - 1;
    }
  else
    {
      *i1 = round(ProcId * ((double) N) / ProcNo) - ri ;
      *iN = N - 1;
    }
}
  
void comptage(int x, int Nx, int *i, int *j)
{
  *j = x / Nx;
  *i = x % Nx;
}
