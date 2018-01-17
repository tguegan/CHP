#include "init.h"

void init_mat(struct matrice_diag A ,int Nx ,int Ny ,double D ,double *Omg, double Beta, double Alpha, int ProcId, int ProcNo, double dx, double dy, double dt, int Num_prob)

{
  double val[5], Id = 0;
  int Nb_diag = 5, i, j, N, k;

  if (Num_prob > 2)
    Id = 1.0;

  N = Nx * Ny; /* Nombre total d'inconnues du processeur courant */

  /*Distance des sous-diag par rapport à la diag principale */
  A.distance[0] = - Nx;
  A.distance[1] = - 1;
  A.distance[2] = 0;
  A.distance[3] = 1;
  A.distance[4] = Nx;

  /* Valeur que prend la matrice a chaque ligne */
  val[0] = - dt / (dy * dy);
  val[1] = - dt / (dx * dx);
  val[2] = (2 / (dx * dx) + 2 / (dy * dy)) * dt + Id / D ;
  val[3] = - dt / (dx * dx);
  val[4] = - dt / (dy * dy); 

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < Nb_diag; j++)
  	{
  	  k = i + A.distance[j];
  	  if((k > -1) && (k < N))
  	    A.valeur[i][j] = val[j]; /*Affectation des valeurs a la matrice*/
  	}
    }

  for (i = N - Nx - 1; i > Nx - 2; i -= Nx)
    A.valeur[i][3] = 0; /*0 correspondant aux intersections des blocs de matrices */

  for (i = Nx; i < N - Nx + 1; i+=Nx)
    A.valeur[i][1] = 0; /*0 correspondant aux intersections des blocs de matrices */
  
  if (ProcNo > 1)
    {
      if (ProcId == ProcNo - 1) // Nx premières lignes pour transmissions Robin pour l'avant dernier proc
	{
	  for (i = 0; i < Nx; i++)
	    {
	      for (j = 0; j < Nb_diag; j++)
		A.valeur[i][j] = 0;
	      A.valeur[i][2] = Beta + Alpha / dy;
	      A.valeur[i][4] = - Alpha / dy;
	    }
	}

      else if (ProcId == 0)  // Nx dernières lignes pour transmissions Robin pour le premier proc
	{
	  for (i = N - Nx; i < N; i++)
	    {
	      for (j = 0; j < Nb_diag; j++)
		A.valeur[i][j] = 0;
	      A.valeur[i][2] = Beta + Alpha / dy;
	      A.valeur[i][0] =  - Alpha / dy;
	    }
	}

      else // Nx dernières et premières lignes pour transmissions Robin pour les autres processeurs
	{
	  for (i = 0; i < Nx; i++)
	    {
	      for (j = 0; j < Nb_diag; j++)
		{
		  A.valeur[i][j] = 0;
		  A.valeur[i + N - Nx][j] = 0;
		}
	      A.valeur[i][2] = Beta + Alpha / dy;
	      A.valeur[i + N - Nx][2] = Beta + Alpha / dy;
	      A.valeur[i][4] = - Alpha / dy;
	      A.valeur[i + N - Nx][0] = - Alpha / dy;
	    }
	}
    }
}

double f(double x, double y, double t, double *Omg, int Num_prob) // Valeur intérieure du domaine
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

double g(double x, double y, int Num_prob) // Conditions de bords horizontales
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

double h(double x, double y, int Num_prob)  // Conditions de bords verticales
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

//**** Second membre adapté pour transmissions Robin ****//
void vectb(double *b, int Nx, int Ny, double *Omg, int Num_prob, double T, int i1, int iN, int ProcId, int ProcNo, double dx, double dy, double D)
{
  int i, j, ni=1, ns=1;
  
  if (ProcId == 0)
    ni = 0;
  
  if ((ProcId == ProcNo - 1) || (ProcNo == 1))
    ns = 0;
  
  for (i = 0; i < Nx; i++)
    {
      for (j = ni; j < iN - i1 + 1 - ns; j++)
	b[i + j * Nx] = f((i + 1) * dx, (j + i1 + 1) * dy , T , Omg, Num_prob) / D; /* Valeur intérieure du domaine */
    }

  for (i = ni; i < iN - i1 + 1 - ns; i++)
    {
      b[i * Nx] +=  h(Omg[0], (i + i1 + 1) * dy, Num_prob) / (D*(dx * dx)); /* Valeur sur les bords verticaux Gamma 1 */
      b[Nx - 1 + i * Nx] +=  h(Omg[1], (i + i1 + 1) * dy, Num_prob) / (D*(dx * dx));
    } 

  if (ProcId == 0)
    for (i = 0; i < Nx; i++)
      b[i] +=  g((i + 1) * dx, Omg[2], Num_prob) / (D*(dy * dy)); /* Valeur sur le bord horizontal bas Gamma 0 */
  
  if ((ProcId == ProcNo - 1) || (ProcNo == 1))
    for (i = 0; i < Nx; i++)
      b[i + (iN - i1) * Nx] +=  g((i + 1) * dx, Omg[3], Num_prob) / (D*(dy * dy)); /* Valeur sur le bord horizontal haut Gamma 0 */

  
}

void update_b(double *b, int Nx, int Ny, int N, int ProcId, int ProcNo)
{
  int i;
  /* Remet à 0 les Nx premiers/derniers éléments de b si besoin */
  for (i = 0; i < Nx; i++)
    {
      if (ProcId == 0)
	b[i + N - Nx] = 0;
      else if (ProcId == ProcNo - 1)
	b[i] = 0;
      else
	{
	  b[i + N - Nx] = 0;
	  b[i] = 0;
	}
    }
}

void charge(int *i1, int *iN, int N, int ProcId, int ProcNo,int ri,int rs)
{
  /* Distribution de la charge avec la méthodre round incluant le recouvrement */
  if ((ProcId > 0) && (ProcId < ProcNo - 1))
    {
      *i1 = round(ProcId * ((double) N) / ProcNo ) - ri;
      *iN = round((ProcId + 1) * (((double) N) / ProcNo)) - 1 + rs;
    }
  else if (ProcId == 0)
    {
      *i1 = round(ProcId * ((double) N) / ProcNo );
      *iN = round((ProcId + 1) * (((double) N) / ProcNo)) - 1 + rs;
    }
  else
    {
      *i1 = round(ProcId * ((double) N) / ProcNo) - ri;
      *iN = N - 1;
    }
}
  
void comptage(int x, int Nx, int *i, int *j)
{
  *j = x / Nx;
  *i = x % Nx;
}
