#define PI 4 * atan(1.0)
#include "init.h"

void init_mat(struct matrice_diag A ,int Nx ,int Ny ,double D ,double *Omg, int i1, int iN)

{

  double dx, dy, val[5];
  int Nb_diag = 5, i, j, N, k;

  N = Nx * (iN - i1 + 1); /*Nombre total d'inconnues*/
  dx = (Omg[1] - Omg[0]) / (Nx + 1.0); /*Le pas suivant l'axe x*/
  dy = (Omg[3] - Omg[2]) / (Ny + 1.0); /* Le pas suivant l'axe y*/

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

void vectb(double *b, int Nx, int Ny, double *Omg, int Num_prob, double T, int i1, int iN)
{
  double dx, dy;
  int i, j;
  
  dx = (Omg[1] - Omg[0]) / (Nx + 1.0); /*Le pas suivant l'axe x*/
  dy = (Omg[3] - Omg[2]) / (Ny + 1.0); /*Le pas suivant l'axe y*/

  for (i = 0; i < Nx; i++)
    {
      for (j = 0; j < iN - i1 + 1; j++)
	b[i + j * Nx] = f((i + 1) * dx, (j + 1) * dy , T , Omg, Num_prob); /*Valeur intérieure du domaine */
    }

  for (i = 0; i < Nx; i++)
    {
      b[i + (Ny - 1) * Nx] += g((i + 1) * dx, Omg[3], Num_prob) / (dy * dy); /* Valeur sur le bord Gamma 0 */
      b[i] += g((i + 1) * dx, Omg[2], Num_prob) / (dy * dy);
    }
  
  for(i = 0; i < iN - i1 + 1; i++)
    {
      b[i * Nx] += h(Omg[0], (i + 1) * dy, Num_prob) / (dx * dx); /*Valeur sur le bord Gamma 1 */
      b[Nx - 1 + i * Nx] += h(Omg[1], (i + 1) * dy, Num_prob) / (dx * dx);
    } 
}

void charge(int *i1, int *iN, int N, int ProcId, int ProcNo)
{
  if (ProcId < ProcNo - 1)
    {
      *i1 = round(ProcId * ((double) N) / ProcNo);
      *iN = round((ProcId + 1) * ((double) N) / ProcNo) - 1;
    }
  else
    {
      *i1 = round(ProcId * ((double) N) / ProcNo);
      *iN = N - 1;
    }
}
  
void comptage(int x, int Nx, int *i, int *j)
{
  *j = x / Nx;
  *i = x % Nx;
}
