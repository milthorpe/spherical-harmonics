#include <cmath>
#include <climits>
#include <cfloat>
#include <cstdio>

#define nsig_BESS	16
#define ensig_BESS	1e16
#define rtnsig_BESS	1e-4
#define enmten_BESS	8.9e-308
#define enten_BESS	1e308

void J_bessel_HalfInt(double x, int nb, double *b, int *ncalc);

