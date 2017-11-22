#include "init.h"
#include "solver.h"
#include "schwarz.h"

int main(int argc , char **argv) {

  int Nx  ,  Ny ,  i1 , iN, ProcId , ProcNo,i , j , k;
  int N, Nb_diag;
  double D = 1.0 , Lx = 1.0 , Ly = 1.0;
  double *Omg, *b, *x;
  FILE* fichier_entree;
  char tmp[5];
  struct matrice_diag A;
 
  fichier_entree = fopen("parameter.txt","r");

  fscanf(fichier_entree,"%s %s %d",tmp,tmp,&Nx);

  fscanf(fichier_entree,"%s %s %d",tmp,tmp,&Ny);

  fscanf(fichier_entree,"%s %s %lf",tmp,tmp,&Lx);

  fscanf(fichier_entree,"%s %s %lf",tmp,tmp,&Ly);

  fscanf(fichier_entree,"%s %s %lf",tmp,tmp,&D);

  fscanf(fichier_entree,"%s %s %d",tmp,tmp,&Nb_diag);

  fclose(fichier_entree);

  Omg = (double*)calloc(4 , sizeof(double));
  Omg[1] = Lx;
  Omg[3] = Ly;
 
  
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&ProcId);
  MPI_Comm_size(MPI_COMM_WORLD,&ProcNo);

  charge( &i1, &iN, Ny, ProcId, ProcNo);
  
  N = Nx * (iN - i1 + 1);
  
  A.valeur = (double**) malloc(N*sizeof(double*));
  for ( i = 0 ; i < N ; i++)
    A.valeur[i] = (double*)calloc(Nb_diag,sizeof(double));
  A.distance = (int*)calloc(Nb_diag,sizeof(int));
  b = (double*) calloc(N, sizeof(double));
  x = (double*) calloc(N, sizeof(double));

  A.distance[0] = -Nx;
  A.distance[1] = -1;
  A.distance[2] = 0;
  A.distance[3] = 1;
  A.distance[4] = Nx;


  init_mat( A ,Nx ,Ny , D ,Omg,  i1, iN);
  vectb(b, Nx, Ny, Omg, 1, 1, i1, iN);
  schwarz(A, x, b, Nx, Ny, Nb_diag, i1, iN, ProcId, ProcNo, Omg);

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
