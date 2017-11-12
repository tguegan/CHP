#include "solver.h"
#include "init.h"

int main(int argc, char *argv[])
{
  struct matrice_diag A;
  double *x, *b, Lx, Ly, D, Omg[4], Tfinal, dt;
  int Nx, Ny, Nb_diag = 5, i, j, Num_prob,  Nt;
  char *Name;

  printf("Entrez Nx :\n");
  scanf("%d",&Nx);

  printf("Entrez Ny :\n");
  scanf("%d",&Ny);

  printf("Entrez Lx :\n");
  scanf("%lf",&Lx);

  printf("Entrez Ly :\n");
  scanf("%lf",&Ly);

  printf("Entrez D :\n");
  scanf("%lf",&D);

  printf("\nQuel probleme voulez-vous resoudre ? :\n");
  printf("1) f = 2*(x^2 - x + y^2 - y) // g = h = 0\n");
  printf("2) f = g = h = sin(x) + cos(y)\n");
  printf("3) Probleme instastionnaire\n");
  scanf("%d", &Num_prob);

  if (Num_prob == 3)
    {
      printf("\nEntrez le temps final :\n");
      scanf("%lf", &Tfinal);
      printf("Entrez Nt\n");
      scanf("%d", &Nt);
    }

  Omg[0] = 0.0;
  Omg[1] = Lx;
  Omg[2] = 0.0;
  Omg[3] = Ly;

  printf("\nInitialisation et calcul de la solution en cours...\n");
  
  b = (double*) calloc(Nx * Ny, sizeof(double));
  x = (double*) calloc(Nx * Ny, sizeof(double));
  A.distance = (int*) calloc(Nb_diag, sizeof(int));
  A.valeur = (double**) calloc(Nx * Ny, sizeof(double*));
  for (i = 0; i < Nx * Ny; i++)
    A.valeur[i] = (double*) calloc(Nx * Ny, sizeof(double));

  switch (Num_prob)
    {
    case 1:
    case 2:
      Name = "solution.plt";
      init_mat(A, Nx, Ny, D, Omg);
      vectb(b, Nx, Ny, Omg, Num_prob, 0);
      grad_conju(A, x, b, Nx, Ny, Nb_diag);
      printf("Erreur max : %.10lf\n", erreur_max(x, Nx, Ny, Omg, Num_prob));
      ecriture_visit(x, Nx, Ny, Omg, Name);
      break;
    case 3:
      dt = Tfinal / Nt;
      init_mat(A, Nx, Ny, D, Omg);
      for (i = 0; i < Nx * Ny; i++)
	{
	  for (j = 0; j < Nb_diag; j++)
	    A.valeur[i][j] *= dt;
	  A.valeur[i][2] += 1;
	}
      for (i = 0; i < Nt; i++)
	{
	  sprintf(Name, "solution_%d.plt", i);
	  vectb(b, Nx, Ny, Omg, Num_prob, i * dt);
	  for (j = 0; j < Nx * Ny; j++)
	    b[j] = dt * b[j] + x[j];
	  grad_conju(A, x, b, Nx, Ny, Nb_diag);
	  ecriture_visit(x, Nx, Ny, Omg, Name);
	}
      break;
    default:
      printf("Erreur lors du choix du probleme.\n");
    }
  
  free(b);
  free(x);
  free(A.distance);
  for (i = 0; i < Nx * Ny; i++)
    free(A.valeur[i]);
  free(A.valeur);
  
  return 0;
}

