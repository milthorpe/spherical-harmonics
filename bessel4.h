#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>

double ftrunc(double x);
double fmax2(double x, double y);

int imin2(int x, int y);

#define nsig_BESS	16
#define ensig_BESS	1e16
#define rtnsig_BESS	1e-4
#define enmten_BESS	8.9e-308
#define enten_BESS	1e308

#define PI 3.141592653589793

void J_bessel_HalfInt(double x, int nb, double *b, int *ncalc);

