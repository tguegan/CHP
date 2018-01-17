#include "init.h"
#include "solver.h"
#include "schwarz.h"
#include <time.h>
#include "string.h"

// Resolution parallèle de l'équation de la chaleur 2D
//  d(u)/dt -D*Laplacien(u) = f avec u =| g au bord gauche/droit
//                                      | h au bord haut/bas
// sur le rectangle [0,Lx]x[0,Ly] en différence finie
// par une méthode de décomposition de domaine avec recouvrement
// et conditions de transimissions mixtes


int main(int argc, char *argv[])
{
  struct matrice_diag A;  // Stockage diagonal de la matrice de discrétisation
  double *Omg, Lx, Ly, D, *b, *x,  Tfinal, max, dt, dx, dy, Alpha, Beta, t1, t2;
  int Init = 1, Export, ProcId, ProcNo, i1, iN, j1, jN, N, Nx, Ny, Nb_diag = 5, i, j, Num_prob, Nt, ri, rs, Precond;
  char Name[50];
  FILE *file;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ProcId);
  MPI_Comm_size(MPI_COMM_WORLD, &ProcNo);


  //*****Lecture et partage des parametres*****//
  ///////////////////////////////////////////////
  ///////////////////////////////////////////////
  if (ProcId == 0)
    {
      file = fopen("input.txt", "r");

      if (file != NULL)
	{
	  fgets(Name, 50, file);
	  fscanf(file, "%d", &Nx);
	  fscanf(file, "%d", &Ny);
	  fscanf(file, "%lf", &Lx);
	  fscanf(file, "%lf", &Ly);
	  fscanf(file, "%lf", &D);
	  fscanf(file, "%lf", &Tfinal);
	  fscanf(file, "%d", &Nt);
	  printf("\nQuel probleme voulez-vous resoudre ? :\n");
	  printf("1) f = 2*(x^2 - x + y^2 - y) // g = h = 0\n");
	  printf("2) f = g = h = sin(x) + cos(y)\n");
	  printf("3) Probleme instastionnaire\n");
	  scanf("%d", &Num_prob);
	  printf("\nPreconditionnement (oui = 1, non = 0)\n");
	  scanf("%d", &Precond);
	  printf("\nTaille recouvrement inferieur (minimal = 1)\n");
	  scanf("%d", &ri);
	  printf("\nTaille recouvrement superieur (minimal = 1)\n");
	  scanf("%d", &rs);
	  printf("\nValeur d'Alpha (Alpha = 1 recommande)\n");
	  scanf("%s", Name);
	  if (strstr(Name, "pi") != NULL)
	    Alpha = M_PI;
	  else
	    Alpha = atof(Name);
	  printf("\nValeur de Beta (Beta = pi recommande)\n");
	  scanf("%s", Name);
	  if (strstr(Name, "pi") != NULL)
	    Beta = M_PI;
	  else
	    Beta = atof(Name);
	  printf("\nExport de la solution pour Visit (oui = 1, non = 0)\n");
	  scanf("%d", &Export);

	  if (ProcNo == 1)
	    ri = rs = 0;
	  fclose(file);
	}

      else
	{
	  printf("\nFichier input.txt manquant.\n");
	  printf("\nGeneration d'un fichier par defaut.\n");
	  file = fopen("input.txt", "w");
	  fprintf(file, "Nx Ny Lx Ly D Tfinal Nt\n");
	  fprintf(file, "500 500 1 1 1 5 50");
	  fclose(file);
	  MPI_Abort(MPI_COMM_WORLD, 0);
	}
    }

  MPI_Bcast(&Nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Lx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Ly, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&D, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Num_prob, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ri, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rs, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Tfinal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Nt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Precond, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Export, 1, MPI_INT, 0, MPI_COMM_WORLD);

  ///////////////////////////////////////////////
  ///////////////////////////////////////////////

  Omg = (double*) calloc(4, sizeof(double)); // Sommets du domaine [0, Lx, 0, Ly]
  Omg[1] = Lx;
  Omg[3] = Ly;

  dx = (Omg[1] - Omg[0]) / (Nx + 1.0); /*Le pas suivant l'axe x*/
  dy = (Omg[3] - Omg[2]) / (Ny + 1.0); /*Le pas suivant l'axe y*/

  if (ProcId == 0)
    printf("\nInitialisation et calcul de la solution en cours...\n");

  charge( &i1, &iN, Ny, ProcId, ProcNo, ri, rs);// Charge en fonction du recouvrement inferieur et superieur
  N = Nx * (iN - i1 + 1);// Nombre total d'inconnues par processeur avec recouvrement

  A.valeur = (double**) calloc(N, sizeof(double*));
  for ( i = 0 ; i < N ; i++)
    A.valeur[i] = (double*) calloc(Nb_diag, sizeof(double));
  A.distance = (int*) calloc(Nb_diag, sizeof(int));
  b = (double*) calloc(N, sizeof(double));
  x = (double*) calloc(N, sizeof(double));

  t1 = MPI_Wtime();
  
  switch (Num_prob)
    {
    case 1:
    case 2:
      sprintf(Name, "solution.plt");
      init_mat(A, Nx, iN - i1 + 1, D, Omg, Beta, Alpha, ProcId, ProcNo, dx, dy, 1, Num_prob); // Initialisation de la matrice de discrétisation
      vectb(b, Nx, Ny, Omg, Num_prob, 0, i1, iN, ProcId, ProcNo, dx, dy, D); // Initialisation du seconb membre du système
      if (Precond == 1)
	precondi(A, b, Nx, iN - i1 + 1, Nb_diag, ProcId, ProcNo); // Preconditionnement de Jacobi - diagonal
      schwarz(A, x, b, Nx, iN - i1 + 1, Nb_diag, ProcId, ProcNo, Omg, ri, rs, Init, Beta, Alpha, dy); // Algorithme de schwarz
      
      max = erreur_max(x, Nx, Ny, Omg, Num_prob, i1, iN, D); // Calcul de l'erreur en norme inf
      MPI_Allreduce(&max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); // Max des max de l'erreur
      if (ProcId == 0)
	printf("Erreur max : %.10lf\n", max);

      if (Export == 1) // Export de la solution pour VIsit
	{
	  charge(&i1, &iN, Ny, ProcId, ProcNo, 0, 0);
	  ecriture_visit(x, Nx, Ny, Omg, Name, ProcId, ProcNo, i1*Nx, iN*Nx + Nx - 1, ri);
	}
      break;
    case 3:
      dt = Tfinal / Nt;
      init_mat(A, Nx, iN - i1 + 1, D, Omg, Beta, Alpha, ProcId, ProcNo, dx, dy, dt, Num_prob);
      charge(&j1, &jN, Ny, ProcId, ProcNo, 0, 0);
      for (i = 0; i < Nt; i++) // Boucle en temps pour un schema d'euler implicite
	{
	  vectb(b, Nx, Ny, Omg, Num_prob, i * dt, i1, iN, ProcId, ProcNo, dx, dy, D); 
	  for (j = 0; j < N; j++)
	    b[j] = x[j] / D + dt * b[j];
	  update_b(b, Nx, Ny, N, ProcId, ProcNo);// Second membre au temps n + 1
	  schwarz(A, x, b, Nx, iN - i1 + 1, Nb_diag, ProcId, ProcNo, Omg, ri, rs, Init, Beta, Alpha, dy);
	  Init = 0;
	  if (Export == 1)
	    {
	      sprintf(Name, "solution_%d.plt", i);
	      ecriture_visit(x, Nx, Ny, Omg, Name, ProcId, ProcNo, j1*Nx, jN*Nx + Nx - 1, ri);
	    }
	}
      break;
    default:
      if (ProcId == 0)
	printf("Erreur lors du choix du probleme.\n");
    }

  t2 = MPI_Wtime();
  t1 = t2 - t1;
  MPI_Allreduce(&t1, &t1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (ProcId == 0)
    printf("\nTemps de resolution : %.10lf s\n", t1);

  for ( i = 0 ; i < N ; i++)
    free(A.valeur[i]);
  free(A.valeur);
  free(A.distance);
  free(Omg);
  free(b);
  free(x);

  MPI_Finalize();
  return 0;
}
