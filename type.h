#ifndef __TYPE__
#define __TYPE__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

struct matrice_diag {
  double **valeur;
  int *distance;
};

#endif
