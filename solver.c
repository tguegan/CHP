#include "solver.h"
#include "init.h"

void precondi(struct matrice_diag A, double *b, int Nx, int Ny, int Nb_diag, int ProcId, int ProcNo)
{
  double P;
  int i, k;
  /* Preconditionnement diagonal par la gauche aussi dit de Jacobi */
  /* Cela revient a diviser chaque ligne de A et b par le terme diagonal de A */
  for (i = 0; i < Nx * Ny; i++)
    {
      P = 1 / A.valeur[i][2];
      b[i] *= P;
      for (k = 0; k < Nb_diag; k++)
	A.valeur[i][k] *= P ;
    }
}

void ecriture_visit(double *x, int Nx, int Ny, double *Omg, char *Name, int ProcId, int ProcNo, int i1, int iN, int rb)
{
  FILE *OutFile;
  double dx, dy;
  int i, k, m, n;
  /* Ecriture d'un fichier *.plt permettant la visualisation de la solution sous visit */
  /* Chaque processeur ecrit l'un apres l'autre dans le fichier */
  /* MPI_Barrier permet de bloquer les processeurs et force une ecriture l'un apres l'autre */
  dx = (Omg[1] - Omg[0]) / (Nx + 1.0);
  dy = (Omg[3] - Omg[2]) / (Ny + 1.0);

  for (k = 0; k < ProcNo; k++)
    {
      if (ProcId == k && k == 0)
	{
	  OutFile = fopen(Name, "w");
	  fprintf(OutFile,"TITLE = \"PROBLEME SOLUTION\" \n");
	  fprintf(OutFile,"VARIABLES = \"X\", \"Y\", \"U\", \"Dom\"\n");
	  fprintf(OutFile,"ZONE T=\"SQUARE\", I=%d, J=%d, F=POINT \n", Nx, Ny);
	  for (i = i1; i < iN + 1; i++)
	    {
	      comptage(i, Nx, &m, &n);
	      fprintf(OutFile,"%lf %lf %lf %d\n", (m + 1) * dx, (n + 1) * dy, x[i - i1], k);
	    }
	  fclose(OutFile);
	}
      else if (ProcId == k)
	{
	  OutFile = fopen(Name, "a");
	  for (i = i1; i < iN + 1; i++)
	    {
	      comptage(i, Nx, &m, &n);
	      fprintf(OutFile,"%lf %lf %lf %d\n", (m + 1) * dx, (n + 1) * dy, x[i + rb*Nx - i1], k);
	    }
	  fclose(OutFile);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
}

double erreur_max(double *x, int Nx, int Ny, double *Omg, int Num_prob , int i1 , int iN, double D)
{
  /* Calcul de l'erreur max entre la solution exacte et approchee sur chaque processeur */
  /* Puis Allreduce du max des max */
  double tmp, max = 0, dx, dy;
  int i, j;

  dx = (Omg[1] - Omg[0]) / (Nx + 1.0);
  dy = (Omg[3] - Omg[2]) / (Ny + 1.0);
  
  for (i = 0; i < Nx; i++)
    {
      for (j = 0; j < iN - i1 + 1; j++)
	{
	  switch (Num_prob)
	    {
	    case 1:
	      tmp = fabs(x[i + j * Nx] - (- (i +1) * dx * (j + i1 + 1) * dy * (Omg[1] - (i + 1) * dx) * (Omg[3] - (j + i1 + 1) * dy)) / D);
	      if (tmp > max)
		max = tmp;
	      break;
	    case 2:
	      tmp = fabs(x[i + j * Nx] - (sin((i + 1) * dx) + cos((j + i1 + 1) * dy)) / D);
	      if (tmp > max)
		max = tmp;
	      break;
	    }
	}
    }

  MPI_Allreduce(&max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return max;
}

void bigradstab(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag, double *r, double *rtilde, double *p, double *v, double *s, double *t)
{
  double normeb = 0, rho1, rho2, alpha, beta, omega, epsi, tmpa, tmpb, tmpc;
  int i, k, l, N, cpt = 0;
  
  /* Methode stabilisee du gradient biconjugue */
  /* x -- solution approchee de Ax = b */
  
  N = Nx * Ny; /*Nombre total d'inconnues*/
  epsi = 1e-8; /* critere d'arret : (|r| / |b| < epsi) */

  for (i = 0; i < N; i++)
    {
      tmpa = 0;
      for (k = 0; k < Nb_diag; k++)
	{
	  l = i + A.distance[k];
	  if ((l < N) && (l >= 0))
	    tmpa += A.valeur[i][k] * x[l] ;
	}
      tmpb = b[i];
      r[i] = tmpb - tmpa; /* r = b - Ax */
      rtilde[i] = tmpb - tmpa; /* rtilde = r */
      normeb += tmpb * tmpb; /* (b, b) */
    }

  normeb = sqrt(normeb); /* norme(b) */

  while (cpt < 20000)
    {
      rho1 = 0;
      for (i = 0; i < N; i++)
	rho1 += rtilde[i] * r[i]; /* rho1 = (rtilde, r) */
      
      if (cpt == 0) /* Si premiere iteration */
	for (i = 0; i < N; i++)
	  p[i] = r[i]; /* p = r */
      else /* Sinon */
	{
	  beta = (rho1 / rho2) * (alpha / omega); /* beta = (rho1 / rho2) * (alpha / omega) */
	  for (i = 0; i < N; i++)
	    p[i] = r[i] + beta * (p[i] - omega * v[i]); /* p = r + beta * (p - omega * v) */
	}

      tmpb = 0;
      for (i = 0; i < N; i++)
	{
	  tmpa = 0;
	  for (k = 0; k < Nb_diag; k++)
	    {
	      l = i + A.distance[k];
	      if ((l < Nx * Ny) && (l >= 0))
		tmpa += A.valeur[i][k] * p[l] ;
	    }
	  v[i] = tmpa; /* v = A * p */
	  tmpb += rtilde[i] * tmpa; /* (rtilde, v) */
	}
      
      alpha = rho1 / tmpb; /* alpha = rho1 / (rtilde, v) */

      tmpb = 0;
      for (i = 0; i < N; i++)
	{
	  tmpa = r[i] - alpha * v[i];
	  s[i] = tmpa; /* s = r - alpha * v */
	  tmpb += tmpa * tmpa; /* (s, s) */
	}

      if (sqrt(tmpb) / normeb < epsi) /* Si (norme(s) / norme(b) < epsi) on a convergence */
	{
	  for (i = 0; i < N; i++)
	    x[i] += alpha * p[i]; /* x = x + alpha * p */
	  break;
	}

      tmpb = 0;
      tmpc = 0;
      for (i = 0; i < Nx * Ny; i++)
	{
	  tmpa = 0;
	  for (k = 0; k < Nb_diag; k++)
	    {
	      l = i + A.distance[k];
	      if ((l < Nx * Ny) && (l >= 0))
		tmpa += A.valeur[i][k] * s[l] ;
	    }
	  t[i] = tmpa; /* t = A * s */
	  tmpb += tmpa * s[i]; /* (t, s) */
	  tmpc += tmpa * tmpa; /* (t, t) */
	}

      omega = tmpb / tmpc; /* omega = (t, s) / (t, t) */

      tmpb = 0;
      for (i = 0; i < N; i++)
	{
	  x[i] += alpha * p[i] + omega * s[i]; /* x = x + alpha * p + omega * s */
	  tmpa = s[i] - omega * t[i];
	  r[i] = tmpa; /* r = s - omega * t */
	  tmpb += tmpa * tmpa; /* (r, r) */
	}

      rho2 = rho1;
      if (sqrt(tmpb) / normeb < epsi) /* Si (norme(r) / norme(b) < epsi) on a convergence */
	break;
      cpt++;
    }
}
