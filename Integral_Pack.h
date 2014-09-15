#ifndef __INTEGRAL_PACK_H_
#define __INTEGRAL_PACK_H_

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "SphericalHarmonics.h"

#ifdef __X10_HEADERS
#include <x10aux/RTT.h>
#endif

typedef struct {int x,y,z;} Point;

namespace au {
namespace edu {
namespace anu {
namespace qm {
namespace ro {
class Integral_Pack {
 public:
#ifdef __X10_HEADERS
RTT_H_DECLS_CLASS;
#endif
    static Integral_Pack* _make(int N, int L, double Type, double roThresh, double rad, double roZ);
    Integral_Pack(int N, int L,double Type,double roThresh, double rad, double roZ);
    ~Integral_Pack();
    void GenclassY(const double *A, const double *B, const double *zetaA, const double *zetaB, int dconA, int dconB, int Ln, double *Ylm);
    void Genclass(int angA, int angB, const double *A, const double *B, const double *zetaA, const double *zetaB, const double *conA, const double *conB, int dconA, int dconB, int n, int Ln, double *Ylm, int maxL, double* aux);
    void Genclass(int angA, int angB, const double *A, const double *B, const double *zetaA, const double *zetaB, const double *conA, const double *conB, int dconA, int dconB, int n, int Ln, double *Ylm, int maxL, int off, double* aux);
    void getNL(int *n_l);

 private:
    int N, L;
    double Type, roZ;
    int Ncal, Nprime;
    double thresh, rad,omega;
    double *arrV;

    SphericalHarmonics* sh;

    // BRA
    #define MAX_BRA_L 10 //for hh
    #define MAX_TOTAL_BRA_L (MAX_BRA_L+1)*(MAX_BRA_L+2)*(MAX_BRA_L+3)/6
    // #define MAX_C 10 // Contraction - for 5Z
    int map3[MAX_BRA_L+1][MAX_BRA_L+1][MAX_BRA_L+1];
    Point inverseMap3[MAX_TOTAL_BRA_L];
    int buildMap[MAX_TOTAL_BRA_L];
    int totalBraL[MAX_BRA_L+2],noOfBra[MAX_BRA_L+1];
    Point *HRRMAP[MAX_BRA_L+1][MAX_BRA_L+1],*HRRMAP2[MAX_BRA_L+1][MAX_BRA_L+1];

    // KET
    #define MAX_KET_L 200
    #define MAX_KET_LM (MAX_KET_L+1)*(MAX_KET_L+1)
    #define Koffset(k,a) ((k)*(K)+(a))
    double cxminus[MAX_KET_LM],cxplus[MAX_KET_LM],cyminus[MAX_KET_LM],cyplus[MAX_KET_LM],cz[MAX_KET_LM];

    double *lambda, *q;

    void initialize();
    void initializeCoulomb(int N);
    void initializeEwald(int N, int L, double Omega, double thresh, double rad);
    void GenJ(double *B, double x, int L);

    void Genclass2(int angA, int angB, const double *A, const double *B, const double *zetaA, const double *zetaB, const double *conA, const double *conB, int dconA, int dconB, int n, int Ln, double *Ylm, int maxL, double* aux);
    void Genclass3(int angA, int angB, const double *A, const double *B, const double *zetaA, const double *zetaB, const double *conA, const double *conB, int dconA, int dconB, int n, int Ln, double *Ylm, int maxL, double* aux);
};
}  // namespace ro
}  // namespace qm
}  // namespace anu
}  // namespace edu
}  // namespace au

#endif  // __INTEGRAL_PACK_H_
