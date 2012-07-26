#ifndef __INTEGRAL_PACK
#define __INTEGRAL_PACK

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define sqr(x) ((x)*(x))
#define SQRT2 1.4142135623730951

typedef struct {int x,y,z;} Point;

namespace au {
    namespace edu {
        namespace anu {
            namespace qm {
                namespace ro {
                    class Integral_Pack {
                        public:
                            static Integral_Pack* _make(int N,int L);
                            Integral_Pack(int N, int L);
                            ~Integral_Pack();
                            void GenclassY(double *A, double *B, double *zetaA, double *zetaB, int dconA, int dconB, int Ln, double *Ylm);
                            void Genclass(int a, int b, double *A, double *B, double *zetaA, double *zetaB, double *conA, double *conB, int dconA, int dconB, double* temp, int N, int Ln, double *Ylm, int maxL);

                        private:
                            int N,L;
                            // BRA
                            #define MAX_BRA_L 10 //for hh
                            #define MAX_TOTAL_BRA_L (MAX_BRA_L+1)*(MAX_BRA_L+2)*(MAX_BRA_L+3)/6
                            // #define MAX_C 10 // Contraction - for 5Z
                            int map3[MAX_BRA_L+1][MAX_BRA_L+1][MAX_BRA_L+1];
                            Point inverseMap3[MAX_TOTAL_BRA_L];
                            int buildMap[MAX_TOTAL_BRA_L];
                            int totalBraL[MAX_BRA_L+2],noOfBra[MAX_BRA_L+1];
                            Point *HRRMAP[MAX_BRA_L+1][MAX_BRA_L+1];

                            // KET
                            #define MAX_KET_L 200
                            #define MAX_KET_LM (MAX_KET_L+1)*(MAX_KET_L+1)
                            #define lm2k(l,m) ((l)*(l)+(l)+(m))
                            #define Koffset(k,a) ((k)*(K)+(a))
                            double cxminus[MAX_KET_LM],cxplus[MAX_KET_LM],cyminus[MAX_KET_LM],cyplus[MAX_KET_LM],cz[MAX_KET_LM];

                            double *lambda, *q;

                            void initialize();
                            void initializeCoulomb(int N);
                            void GenJ(double *B, double x, int L);
                            void GenY(double *Y, double X, double phi, int L);

                    };
                }
            }
        }
    }
}

#endif // __INTEGRAL_PACK
