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
  


