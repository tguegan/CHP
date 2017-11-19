#include "solver.h"

void prod_mat_vect(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag,int ProcId, int i1, int iN, int ProcNo, double *x_etendu, double *x_p, double *x_s)
{
    int i, k;
    MPI_Request Request;

    for (i = 0 ; i< iN - i1 + 1; i++)
        x_etendu[i + Nx]=x[i];

    for (i = 0; i < Nx; i++)
    {
        x_p[i]=x[i];
        x_s[i]=x[i + iN - i1 - Nx + 1];
    }

    if (ProcNo > 1)
    {
        if (ProcId == 0)
        {
            MPI_Isend(x_s, Nx, MPI_DOUBLE, ProcId + 1, 102, MPI_COMM_WORLD, &Request);
            MPI_Irecv(x_p, Nx, MPI_DOUBLE, ProcId + 1, 101, MPI_COMM_WORLD, &Request);
        }
        else if (ProcId == ProcNo - 1)
        {
            MPI_Isend(x_p, Nx, MPI_DOUBLE, ProcId - 1, 101, MPI_COMM_WORLD, &Request);
            MPI_Irecv(x_etendu, Nx, MPI_DOUBLE, ProcId - 1, 102, MPI_COMM_WORLD, &Request);
        }
        else
        {
            MPI_Isend(x_p, Nx, MPI_DOUBLE, ProcId - 1, 101, MPI_COMM_WORLD, &Request);
            MPI_Isend(x_s, Nx, MPI_DOUBLE, ProcId + 1, 102, MPI_COMM_WORLD, &Request);
            MPI_Irecv(x_etendu, Nx, MPI_DOUBLE, ProcId - 1, 102, MPI_COMM_WORLD, &Request);
            MPI_Irecv(x_p, Nx, MPI_DOUBLE, ProcId + 1, 101, MPI_COMM_WORLD, &Request);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (ProcId == 0)
            for (i = iN - i1 + 1; i < iN - i1 + 1 + Nx; i++)
                x_etendu[i + Nx] = x_p[i - iN + i1 - 1];
        else if (ProcId < ProcNo - 1)
            for (i = iN - i1 + 1; i < iN - i1 + 1 + Nx; i++)
                x_etendu[i + Nx]=x_p[i - iN + i1 - 1];
    }

    for (i = 0; i < iN - i1 + 1; i++)
    {
        b[i] = 0;
        for (k = 0; k < Nb_diag; k++)
        {
            b[i] += A.valeur[i][k] * x_etendu[i + A.distance[k] + Nx];
        }
    }
}

double prod_scal(double *x, double *y, int Nx, int Ny, int i1, int iN)
{
    double somme = 0;
    int i;

    for (i = 0; i < iN - i1 + 1; i++)
        somme += x[i] * y[i];

    return somme;
}

double norme_vect(double *x, int Nx, int Ny, int i1, int iN)
{
    return sqrt(prod_scal(x, x, Nx, Ny,i1, iN));
}

void grad_conju(struct matrice_diag A, double *x, double *b, int Nx, int Ny, int Nb_diag, int i1, int iN, int ProcId, int ProcNo, double *x_etendu, double *x_p, double *x_s, double *r, double *p, double *Ap)
{
    double alp, beta, tmpa, tmpb, epsi, normer = 1, normeb;
    int i, N, cpt = 0;

    N = iN - i1 + 1;
    epsi = 1e-8;

    prod_mat_vect(A, x, r, Nx, Ny, Nb_diag, ProcId, i1, iN, ProcNo, x_etendu, x_p, x_s);

    for (i = 0; i < N; i++)
    {
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }

    normeb = prod_scal(b, b, Nx, Ny, i1, iN);
    MPI_Allreduce(&normeb, &normeb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    normeb = sqrt(normeb);

    while ((normer > epsi) && (cpt < 20000))
    {
        prod_mat_vect(A, p, Ap, Nx, Ny, Nb_diag, ProcId, i1, iN, ProcNo, x_etendu, x_p, x_s);

        tmpa = prod_scal(r, p, Nx, Ny, i1, iN);
        tmpb = prod_scal(Ap, p, Nx, Ny, i1, iN);
        MPI_Allreduce(&tmpa, &tmpa, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&tmpb, &tmpb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alp = tmpa / tmpb;

        for (i = 0; i < N; i++)
        {
            x[i] += alp * p[i];
            Ap[i] = r[i] - alp * Ap[i];
        }

        tmpb = prod_scal(Ap, Ap, Nx, Ny, i1, iN);
        MPI_Allreduce(&tmpb, &tmpb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta = tmpb / tmpa;

        for (i = 0; i < N; i++)
        {
            p[i] = Ap[i] + beta * p[i];
            r[i] = Ap[i];
        }
        normer = sqrt(tmpb) / normeb;
        cpt ++;
    }

    if (ProcId == 0)
        printf("Nombre d'iteration : %d\n", cpt);
}

double erreur_max(double *x, int Nx, int Ny, double *Omg, int Num_prob, int i1, int iN)
{
    double tmp, max = 0, dx, dy;
    int i, j, k;

    dx = (Omg[1] - Omg[0]) / (Nx + 1.0);
    dy = (Omg[3] - Omg[2]) / (Ny + 1.0);

    for (i = i1; i < iN + 1; i++)
    {
        comptage(i, Nx, &j, &k);
        switch (Num_prob)
        {
        case 1:
            tmp = fabs(x[i - i1] - (- (j + 1) * dx * (k + 1) * dy * (Omg[1] - (j + 1) * dx) * (Omg[3] - (k + 1) * dy)));
            if (tmp > max)
                max = tmp;
            break;
        case 2:
            tmp = fabs(x[i - i1] - (sin((j + 1) * dx) + cos((k + 1) * dy)));
            if ((tmp > max))
                max = tmp;
            break;
        }
    }

    MPI_Allreduce(&max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return max;
}

void ecriture_visit(double *x, int Nx, int Ny, double *Omg, char *Name)
{
    FILE *OutFile;
    double dx, dy;
    int i, j;

    dx = (Omg[1] - Omg[0]) / (Nx + 1.0);
    dy = (Omg[3] - Omg[2]) / (Ny + 1.0);

    OutFile = fopen(Name, "w");

    fprintf(OutFile,"TITLE = \"PROBLEME SOLUTION\" \n");
    fprintf(OutFile,"VARIABLES = \"X\", \"Y\", \"U\" \n");
    fprintf(OutFile,"ZONE T=\"SQUARE\", I=%d, J=%d, F=POINT \n", Nx, Ny);

    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
            fprintf(OutFile,"%lf %lf %lf\n", (i + 1) * dx, (j + 1) * dy, x[i + j * Nx]);
    }
    fclose(OutFile);
}
