#include "init.h"
#include "solver.h"
#include <time.h>

int main(int argc, char *argv[])
{
    struct matrice_diag A;
    double Omg[4], Lx, Ly, D, *b, *x, *x_etendu, *x_p, *x_s, *r, *p, *Ap, Tfinal;
    int ProcId, ProcNo, i1, iN, Nx, Ny, Nb_diag = 5, i, Num_prob, Nt;
    unsigned long duration;
    char Name[50];
    FILE *file;

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

    if(Num_prob == 3)
      {
	MPI_Bcast(&Tfinal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Nt, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
    
    clock_t start = clock();
    
    Omg[0] = 0.0;
    Omg[1] = Lx;
    Omg[2] = 0.0;
    Omg[3] = Ly;

    if (ProcId == 0)
      printf("\nInitialisation et calcul de la solution en cours...\n");

    charge(&i1, &iN, Nx * Ny, ProcId, ProcNo);

    x = (double*) calloc(iN - i1 + 1, sizeof(double));
    b = (double*) calloc(iN - i1 + 1, sizeof(double));
    x_etendu = (double*) calloc(iN - i1 + 1 + 2 * Nx, sizeof(double));
    x_p = (double*) calloc(Nx, sizeof(double));
    x_s = (double*) calloc(Nx, sizeof(double));
    r = (double*) calloc(iN - i1 + 1, sizeof(double));
    p = (double*) calloc(iN - i1 + 1, sizeof(double));
    Ap = (double*) calloc(iN - i1 + 1, sizeof(double));
    A.distance = (int*) calloc(Nb_diag, sizeof(int));
    A.valeur = (double**) calloc(iN - i1 + 1, sizeof(double*));
    for (i = 0; i < iN - i1 + 1; i++)
      A.valeur[i] = (double*) calloc(Nb_diag, sizeof(double));

    switch (Num_prob)
      {
      case 1:
      case 2:
	init_mat(A, Nx, Ny, D, Omg, i1, iN);
	vectb(b, Nx, Ny, Omg, Num_prob, 0, i1, iN);
	grad_conju(A, x, b, Nx, Ny, Nb_diag, i1, iN, ProcId, ProcNo, x_etendu, x_p, x_s, r, p, Ap);
	if (ProcId == 0)
	  printf("Erreur max : %.10lf\n", erreur_max(x, Nx, Ny, Omg, Num_prob, i1, iN));
	break;
      default:
	if (ProcId == 0)
	  printf("Erreur lors du choix du probleme.\n");
      }
    
    clock_t end = clock();
    duration = (end - start ) * 1000 / CLOCKS_PER_SEC;
    MPI_Allreduce(&duration, &duration, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    if (ProcId == 0)
      printf("\nTemps de resolution : %ld ms\n", duration);

    free(x);
    free(x_etendu);
    free(x_p);
    free(x_s);
    free(r);
    free(p);
    free(Ap);
    free(b);
    free(A.distance);
    for (i = 0; i < iN - i1 + 1; i++)
        free(A.valeur[i]);
    free(A.valeur);

    MPI_Finalize();
    return 0;
}

