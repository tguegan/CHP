#include "init.h"

int main(int argc , char **argv) {

  int Nx  ,  Ny ,  i1 , iN, ProcID , ProcNo,i , j , k;
  int N;
  int Nb_diag = 5;
  double D = 1.0 , Lx = 1.0 , Ly = 1.0;
  double* Omg;
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
  MPI_Comm_rank(MPI_COMM_WORLD,&ProcID);
  MPI_Comm_size(MPI_COMM_WORLD,&ProcNo);

  charge( &i1, &iN, Ny, ProcID, ProcNo);
  
  N = Nx * (iN - i1 + 1);
  
  A.valeur = (double**)malloc(N*sizeof(double));
  for ( i = 0 ; i < N ; i++)
    A.valeur[i] = (double*)calloc(Nb_diag,sizeof(double));
  
  A.distance = (int*)calloc(Nb_diag,sizeof(int));
  A.distance[0] = -Nx;
  A.distance[1] = -1;
  A.distance[2] = 0;
  A.distance[3] = 1;
  A.distance[4] = Nx;


  init_mat( A ,Nx ,Ny , D ,Omg,  i1, iN);

    for ( k = 0 ; k < ProcNo ; k++) {

      if ( k == ProcID ) {
	printf(" \n Voici la matrice pour le processeur numero %i i1 = %i  iN = %i \n",ProcID,i1,iN);

	for ( i = 0 ; i < N ; i++) {

	  for( j = 0 ; j < Nb_diag ; j++) {

	    printf("%lf  " , A.valeur[i][j]);
	  }
	  printf("\n");
	}

      }
      MPI_Barrier(MPI_COMM_WORLD);
    
    }

    
  for ( i = 0 ; i < N ; i++)
    free(A.valeur[i]);
  free(A.valeur);
  free(A.distance);
  free(Omg);


  MPI_Finalize();

}

  
