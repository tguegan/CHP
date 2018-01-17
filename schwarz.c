#include "schwarz.h"
#include "solver.h"

/////////////////////////////////////////////////
//******* Méthode de Schwarz additive ********//
///////////////////////////////////////////////
 
void schwarz(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag, int ProcId, int ProcNo, double *Omg, int ri, int rs, int Init, double Beta, double Alpha, double dy)
{
  int i, N, cpt = 0;
  double epsi, normebord = 1, normebordpre = 0, *xcom, *btmp, *x_ns, *x_np, *xhaut, *xbas, *r, *rtilde, *p, *ptilde, *Ap, *Aptilde;
  MPI_Status Status;

  epsi = 1e-6;
  N = Nx * Ny;

  x_np = (double*) calloc(Nx, sizeof(double));  //
  x_ns = (double*) calloc(Nx, sizeof(double));  //
  xhaut = (double*) calloc(Nx, sizeof(double)); // Allocations pour les commununications 
  xbas = (double*) calloc(Nx, sizeof(double));  //
  btmp = (double*) calloc(N, sizeof(double));   //
  xcom = (double*) calloc(N, sizeof(double));   //
  r = (double*) calloc(N, sizeof(double));        //
  rtilde = (double*) calloc(N, sizeof(double));   //
  p = (double*) calloc(N, sizeof(double));        // Allocations pour le BiCGstab
  ptilde = (double*) calloc(N, sizeof(double));   //
  Ap = (double*) calloc(N, sizeof(double));       //
  Aptilde = (double*) calloc(N, sizeof(double));  //

  if (Init == 1)
    bigradstab(A, x, b, Nx, Ny, Nb_diag, r, rtilde, p, ptilde, Ap, Aptilde); // Initialisation d'une solution dans chaque sous domaine

  while ((normebord > epsi) && (fabs(normebord - normebordpre) > epsi) && ProcNo > 1)  // Critére d'arrêt sur la continuité des interfaces
    {
      cpt++;
      normebordpre = normebord;
      for (i = 0; i < Nx; i++)       // Stockage des 2 lignes à envoyer (au lieu de 4) en ayant construit la dérivée avant l'envoi
	{
	  x_ns[i] = A.valeur[i + N - Nx][0] * x[i + (Ny - rs - ri + 1) * Nx] + A.valeur[i + N - Nx][2] * x[i + (Ny - rs - ri) * Nx] ; // pour proc suivant
	  x_np[i] = A.valeur[i][4] * x[i + (rs + ri - 2) * Nx] + A.valeur[i][2] * x[i + (rs + ri - 1) * Nx];  // pour proc precédant
	}

      if (ProcId == 0)                                                                                                    // Si ProcId = 0
	MPI_Sendrecv(x_ns, Nx, MPI_DOUBLE, ProcId + 1, 0, x_np, Nx, MPI_DOUBLE, ProcId + 1, 0, MPI_COMM_WORLD, &Status);  // Communications avec le proc suivant uniquement

      else if (ProcId < ProcNo - 1)
	{
	  MPI_Sendrecv(x_np, Nx, MPI_DOUBLE, ProcId - 1, 0, xcom, Nx, MPI_DOUBLE, ProcId - 1, 0, MPI_COMM_WORLD, &Status);// Communications avec le proc précédant
	  MPI_Sendrecv(x_ns, Nx, MPI_DOUBLE, ProcId + 1, 0, x_np, Nx, MPI_DOUBLE, ProcId + 1, 0, MPI_COMM_WORLD, &Status);// "           " avec le proc suivant
	}

      else                                                                                                                // Si ProcId = ProcNo -1
	MPI_Sendrecv(x_np, Nx, MPI_DOUBLE, ProcId - 1, 0, xcom, Nx, MPI_DOUBLE, ProcId - 1, 0, MPI_COMM_WORLD, &Status);  // Communications avec le proc précédant uniquement
	    
      if (ProcId < ProcNo - 1)
	for (i = 0; i < Nx; i++)
	  xcom[N - Nx + i] = x_np[i];    

      for (i = 0; i < N; i++)
	btmp[i] = b[i] + xcom[i];           // second membre temporaire sans modifier le second membre initial
      
      bigradstab(A, x, btmp, Nx, Ny, Nb_diag, r, rtilde, p, ptilde, Ap, Aptilde);

      critere_schwarz(x, x_np, x_ns, xhaut, xbas, Nx, Ny, ri, rs, ProcId, ProcNo, &normebord);   // Critére d'arrêt

      if (ProcId == 0)
	printf("\nCompteur de Schwarz : %d, norme max des interfaces : %.15lf\n", cpt, normebord);
    }

  free(btmp);
  free(xcom);
  free(x_ns);
  free(x_np);
  free(xhaut);
  free(xbas);
  free(r);
  free(rtilde);
  free(p);
  free(ptilde);
  free(Ap);
  free(Aptilde);
}

void critere_schwarz(double *x, double *x_np, double *x_ns, double *xhaut, double *xbas, int Nx, int Ny, int ri, int rs, int ProcId, int ProcNo, double *normexbord)
{

  // Critére d'arrêt sur la continuité des interfaces
  // Nécessite donc le même type de communication pour échanger les valeurs aux interfaces
  
  int i;
  double normexb, normexh, tmpa, tmpb;
  MPI_Status Status;

  for (i = 0; i < Nx; i++)
    {
      x_ns[i] = x[i + (Ny - rs - ri) * Nx]; // ligne correspondant au bord du domaine suivant
      x_np[i] = x[i + (rs + ri - 1) * Nx];  // ligne correspondant au bord du domaine précédant
    }
  
  if (ProcId == 0)
    MPI_Sendrecv(x_ns, Nx, MPI_DOUBLE, ProcId + 1, 0, xhaut, Nx, MPI_DOUBLE, ProcId + 1, 0, MPI_COMM_WORLD, &Status);
  else if (ProcId < ProcNo - 1)
    {
      MPI_Sendrecv(x_np, Nx, MPI_DOUBLE, ProcId - 1, 0, xbas, Nx, MPI_DOUBLE, ProcId - 1, 0, MPI_COMM_WORLD, &Status);
      MPI_Sendrecv(x_ns, Nx, MPI_DOUBLE, ProcId + 1, 0, xhaut, Nx, MPI_DOUBLE, ProcId + 1, 0, MPI_COMM_WORLD, &Status);
    }
  else
    MPI_Sendrecv(x_np, Nx, MPI_DOUBLE, ProcId - 1, 0, xbas, Nx, MPI_DOUBLE, ProcId - 1, 0, MPI_COMM_WORLD, &Status);

  normexh = 0;
  normexb = 0;
  for (i = 0; i < Nx; i++)
    {
      tmpa = xhaut[i] - x[Nx * (Ny - 1) + i];
      tmpb = xbas[i] - x[i];
      normexh += tmpa * tmpa;
      normexb += tmpb * tmpb;
    }

  normexh = sqrt(normexh);
  normexb = sqrt(normexb);

  if (ProcId == 0)  // Pas d'interface inférieure à vérifier
    normexb = 0;
  if (ProcId == ProcNo - 1) // Pas d'interface supérieure à vérifier
    normexh = 0;

  if (normexh < normexb)
    normexh = normexb;

  MPI_Allreduce(&normexh, &normexh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); // critére global

  *normexbord = normexh;
}
