#include "init.h"
#include "solver.h"
#include "schwarz.h"
#include <time.h>

int main(int argc, char *argv[])
{
  struct matrice_diag A;
  double *Omg, Lx, Ly, D, *b, *x,  Tfinal, max, dx, dy, dt;
  int ProcId, ProcNo, j1, jN, i1, iN, N, Nx, Ny, Nb_diag = 5, i, j, Num_prob, Nt, rb, rh, rg, rd, Nbdomv, Nbdomh;
  unsigned long duration;
  char Name[50];
  FILE *file;
  MPI_Status Status;
 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ProcId);
  MPI_Comm_size(MPI_COMM_WORLD, &ProcNo);

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
	  printf("\nNombre de sous-domaine vertical ?\n");
	  scanf("%d", &Nbdomv);
	  printf("\nNombre de sous-domaine horizontal ?\n");
	  scanf("%d", &Nbdomh);
	  if (ProcNo != Nbdomv * Nbdomh)
	    {
	      printf("\nNombre de sous-domaine incompatible avec nombre de processeur.\n");
	      MPI_Abort(MPI_COMM_WORLD, 2);
	    }
	  printf("\nQuel probleme voulez-vous resoudre ? :\n");
	  printf("1) f = 2*(x^2 - x + y^2 - y) // g = h = 0\n");
	    printf("2) f = g = h = sin(x) + cos(y)\n");
	    printf("3) Probleme instastionnaire\n");
	    scanf("%d", &Num_prob);

	  if (Num_prob == 3)
	    {
	      fscanf(file, "%lf", &Tfinal);
	      fscanf(file, "%d", &Nt);
	    }

	  printf("\nTaille recouvrement inferieur (0 = minimal)\n");
	  scanf("%d", &rb);
	  printf("\nTaille recouvrement superieur (0 = minimal)\n");
	  scanf("%d", &rh);
	  printf("\nTaille recouvrement lateral gauche (0 = minimal)\n");
	  scanf("%d", &rg);
	  printf("\nTaille recouvrement lateral droit (0 = minimal)\n");
	  scanf("%d", &rd);
	  
	  fclose(file);
	}
      
      else
	{
	  printf("\nFichier input.txt manquant.\n");
	  MPI_Abort(MPI_COMM_WORLD, 0);
	}
    }
    
  MPI_Bcast(&Nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Lx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Ly, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&D, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Num_prob, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rb, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rh, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rd, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Nbdomv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Nbdomh, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(Num_prob == 3)
    {
      MPI_Bcast(&Tfinal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&Nt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

  for (i = 0; i < Nbdomv; i++)
    {
      if (ProcId == i)
	{
	  charge(&j1, &jN, Nx, ProcId, Nbdomv, rg, rd);
	  for (j = 1; j < Nbdomh; j++)
	    {
	      MPI_Send(&j1, 1, MPI_INT, ProcId + j * Nbdomv, 100, MPI_COMM_WORLD);
	      MPI_Send(&jN, 1, MPI_INT, ProcId + j * Nbdomv, 100, MPI_COMM_WORLD);
	    }
	}
      else if ((ProcId > Nbdomv - 1) && (ProcId % Nbdomv == i))
	{
	  MPI_Recv(&j1, 1, MPI_INT, i, 100, MPI_COMM_WORLD, &Status);
	  MPI_Recv(&jN, 1, MPI_INT, i, 100, MPI_COMM_WORLD, &Status);
	}
    }

  for (i = 0; i < Nbdomh; i++)
    {
      if ((ProcId % Nbdomv == 0) && (ProcId / Nbdomv == i))
	{
	  charge(&i1, &iN, Ny, i, Nbdomh, rb, rh);
	  for (j = 1; j < Nbdomv; j++)
	    {
	      MPI_Send(&i1, 1, MPI_INT, ProcId + j, 100, MPI_COMM_WORLD);
	      MPI_Send(&iN, 1, MPI_INT, ProcId + j, 100, MPI_COMM_WORLD);
	    }
	}
      else if (ProcId / Nbdomv == i)
      	{
      	  MPI_Recv(&i1, 1, MPI_INT, i * Nbdomv, 100, MPI_COMM_WORLD, &Status);
      	  MPI_Recv(&iN, 1, MPI_INT, i * Nbdomv, 100, MPI_COMM_WORLD, &Status);
      	}
    }
  
  clock_t start = clock();

  Omg = (double*) calloc(4 , sizeof(double));
  Omg[1] = Lx;
  Omg[3] = Ly;

  dx = (Omg[1] - Omg[0]) / (Nx + 1.0);
  dy = (Omg[3] - Omg[2]) / (Ny + 1.0);

  if (ProcId == 0)
    printf("\nInitialisation et calcul de la solution en cours...\n");

  N = (iN - i1 + 1) * (jN - j1 + 1);
  
  A.valeur = (double**) malloc(N*sizeof(double*));
  for ( i = 0 ; i < N ; i++)
    A.valeur[i] = (double*)calloc(Nb_diag,sizeof(double));
  A.distance = (int*)calloc(Nb_diag,sizeof(int));
  b = (double*) calloc(N, sizeof(double));
  x = (double*) calloc(N, sizeof(double));


  init_mat(A, jN - j1 + 1, iN - i1 + 1 , D, Omg, dx, dy);

  /* for (rs = 0; rs < ProcNo; rs++) */
  /*   { */
  /*     if (ProcId == rs) */
  /* 	{ */
  /* 	  printf("\nProcId : %d  j1 = %d  jN = %d  i1 = %d  iN = %d\n", ProcId, j1, jN, i1, iN); */
  /* 	  for (i = 0; i < N; i++) */
  /* 	    { */
  /* 	      for (j = 0; j < Nb_diag; j++) */
  /* 		printf("%lf   ", A.valeur[i][j]); */
  /* 	      printf("\n"); */
  /* 	    } */
  /* 	} */
  /*     MPI_Barrier(MPI_COMM_WORLD); */
  /*   } */
  
  vectb(b, jN - j1 + 1, iN - i1 + 1, Omg, Num_prob, 0, j1, i1, ProcId, ProcNo, dx, dy, Nbdomv, Nbdomh);
  schwarz(A, x, b, Nx, Ny, Nb_diag, i1, iN, ProcId, ProcNo, Omg, rb, rh, 1);
  max = erreur_max(x, Nx, Ny, Omg, Num_prob, i1, iN);
  if (ProcId == 0)
    printf("Erreur max : %.10lf\n", max);
  
  charge(&i1, &iN, Ny, ProcId, ProcNo, 0, 0);
  ecriture_visit(x, Nx, Ny, Omg, "sol.plt", ProcId, ProcNo, i1*(jN-j1+1), iN*(jN-j1+1) + jN-j1+1 - 1, rb);
  
  clock_t end = clock();
  duration = (end - start) * 1000 / CLOCKS_PER_SEC;
  MPI_Allreduce(&duration, &duration, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  if (ProcId == 0)
    printf("\nTemps de resolution : %ld ms\n", duration);
  
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
