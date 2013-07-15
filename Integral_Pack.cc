#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "bessel4.h"
#include "Integral_Pack.h"

#ifdef __X10_HEADERS
#include <x10aux/config.h>
#include <x10aux/RTT.h>
using namespace x10aux;
#endif

#define sqr(x) ((x)*(x))
#define SQRT2 1.4142135623730951

// Papers
// LHG2012 (RO#6): 10.1063/1.3691829
// LMRG2013 (RO#7): 10.1021/ct301110y

int delta[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double JpY00[11]={0.28209479177387814, 0.09403159725795937, 0.018806319451591877,
       		0.0026866170645131254, 0.000298513007168125,0.000027137546106193183, 
                2.087503546630245e-6,1.39166903108683e-7, 8.186288418157823e-9,
                4.308572851662012e-10, 2.0517013579342914e-11}; // LMRG2013 eq 19b
// 0-10 accomodate (hh|phi) - need more values in the array for higer angular momentum 

namespace au {
    namespace edu {
        namespace anu {
            namespace qm {
                namespace ro {

#ifdef __X10_HEADERS
RTT_CC_DECLS0(Integral_Pack, "Integral_Pack", RuntimeType::class_kind)
#endif

    Integral_Pack* Integral_Pack::_make(int N, int L,double Type, double roThresh, double mrad, double roZ) {
      return new Integral_Pack(N, L, Type, roThresh, mrad, roZ);
    }

    Integral_Pack::Integral_Pack(int N, int L, double Type, double roThresh, double mrad, double roZ) {
        this->N = N; this->L=L; this->Type=Type;
        initialize();
        thresh=roThresh; rad=mrad; omega=Type; this->roZ=roZ;
        if (Type<=0.) initializeCoulomb(N);
        else initializeEwald(N, L, Type, roThresh, mrad);         
        arrV=(double *)malloc(totalBraL[MAX_BRA_L+1]*(L+1)*(L+1)*sizeof(double)*2);
        // printf("Integral_Pack.cc N=%d L=%d type=%f\n",N,L,Type);
    }

    Integral_Pack::~Integral_Pack() {
        free(lambda); free(q); free(arrV); 
        for (int i=1; i<=MAX_BRA_L; i++) for (int j=1; j<=i; j++) {
            free(HRRMAP[i][j]); free(HRRMAP2[j][i]);
        }
    }

    void Integral_Pack::initialize() {
        int x,y,z,l,m=0;
        totalBraL[0]=0;
        totalBraL[1]=1;
        for (l=0; l<=MAX_BRA_L; l++) {
            noOfBra[l]=((l+1)*(l+2))/2;
            if (l>0) totalBraL[l+1]=totalBraL[l]+noOfBra[l];    
            for (x=0; x<=l; x++) for (y=0; y<=l-x; y++) { // must be consistent with PowerList.x10
                map3[x][y][z=l-x-y] = m;
                inverseMap3[m].x = x;
                inverseMap3[m].y = y;
                inverseMap3[m].z = z;
                // LHG2012 Table II & boolean before (36)
                if (z == 1) buildMap[m]=2; //z
                else if (y == 1) buildMap[m]=1; //y
                else if (x == 1) buildMap[m]=0; //x
                else if (z >= 2) buildMap[m]=2; //z
                else if (y >= 2) buildMap[m]=1; //y
                else buildMap[m]=0; //x
                m++;
            }
        }
        // HRR
        for (int i=1; i<=MAX_BRA_L; i++) for (int j=1; j<=i; j++) {
            HRRMAP[i][j] = (Point *) malloc(sizeof(Point)*noOfBra[i]*noOfBra[j]);
            if (HRRMAP[i][j]==NULL) {printf("Integral_Pack.cc malloc failed at ln75\n"); exit(1);}
            for (int a=0; a<noOfBra[i]; a++) for (int b=0; b<noOfBra[j]; b++) {
             	int aInt=totalBraL[i]+a,bInt=totalBraL[j]+b;
             	int ax=inverseMap3[aInt].x, ay=inverseMap3[aInt].y, az=inverseMap3[aInt].z, bx=inverseMap3[bInt].x, by=inverseMap3[bInt].y, bz=inverseMap3[bInt].z;
             	int increment;
             	if (bx) increment=0; else if (by) increment=1; else /*if (bz)*/ increment=2;
             	int bm1 = map3[bx-delta[0][increment]][by-delta[1][increment]][bz-delta[2][increment]]-totalBraL[j-1];
             	int ap1 = map3[ax+delta[0][increment]][ay+delta[1][increment]][az+delta[2][increment]]-totalBraL[i+1];
             	int leftindex=noOfBra[j]*a+b;
             	int rightindexA=noOfBra[j-1]*ap1+bm1;
             	int rightindexB=noOfBra[j-1]*a+bm1;
                HRRMAP[i][j][leftindex].x=rightindexA;
                HRRMAP[i][j][leftindex].y=rightindexB;
                HRRMAP[i][j][leftindex].z=increment;
            }
        }
        // HRR2 
        
        for (int i=1; i<=MAX_BRA_L; i++) for (int j=1; j<=i; j++) { // ok
            HRRMAP2[j][i] = (Point *) malloc(sizeof(Point)*noOfBra[i]*noOfBra[j]); // a<b
            if (HRRMAP2[j][i]==NULL) {printf("Integral_Pack.cc malloc failed at ln104\n"); exit(1);} // a<b
            for (int b=0; b<noOfBra[i]; b++) for (int a=0; a<noOfBra[j]; a++) { // a<b
             	int aInt=totalBraL[j]+a, bInt=totalBraL[i]+b, increment; // a<b
             	int ax=inverseMap3[aInt].x, ay=inverseMap3[aInt].y, az=inverseMap3[aInt].z, bx=inverseMap3[bInt].x, by=inverseMap3[bInt].y, bz=inverseMap3[bInt].z;
             	if (ax) increment=0; else if (ay) increment=1; else increment=2; // a<b
             	int am1 = map3[ax-delta[0][increment]][ay-delta[1][increment]][az-delta[2][increment]]-totalBraL[j-1]; // ok
             	int bp1 = map3[bx+delta[0][increment]][by+delta[1][increment]][bz+delta[2][increment]]-totalBraL[i+1]; // ok
             	int leftindex=noOfBra[i]*a+b; // a<b
             	int rightindexA=noOfBra[i-1]*am1+bp1; // a<b
             	int rightindexB=noOfBra[i-1]*am1+b; // a<b
                HRRMAP2[j][i][leftindex].x=rightindexA;
                HRRMAP2[j][i][leftindex].y=rightindexB;
                HRRMAP2[j][i][leftindex].z=increment;
            }
        } 

        // LMRG2013 eq 22a-22f
        for (l=0; l<=MAX_KET_L; l++) for (m=-l; m<=l; m++) {
            double cminus=.5*sqrt((l-m-1.0)*(l-m)*(2.0*l+1.0)/(2.0*l-1.0)), cplus=.5*sqrt((l+m-1.0)*(l+m)*(2.0*l+1.0)/(2.0*l-1.0));
            int index=lm2k(l,m);
            if (m<=-2) cxplus[index]=-cminus;
                else if (m>=1) cxplus[index]=cminus;
                else if (m==-1) cxplus[index]=0;
                else /*if (m==0)*/ cxplus[index]=SQRT2*cminus;
            if (m<=-1) cxminus[index]=cplus;
                else if (m>=2) cxminus[index]=-cplus;
                else if (m==0) cxminus[index]=0;
                else /*if (m==1)*/ cxminus[index]=-SQRT2*cplus;
            if (m<=-2) cyplus[index]=-cminus;
                else if (m>=1) cyplus[index]=cminus;
                else if (m==-1) cyplus[index]=-SQRT2*cminus;
                else /*if (m==0)*/ cyplus[index]=SQRT2*cminus;
            if (m<=-1) cyminus[index]=-cplus;
                else if (m>=2) cyminus[index]=cplus;
                else /*if (m==0 || m==1)*/ cyminus[index]=0;
            cz[index] = sqrt((l*l-m*m)*(2.0*l+1.0)/(2.0*l-1.0));
        }
    }

    void Integral_Pack::GenJ(double *B, double x, int L) {
        int l,nb=L+1,ncalc;
        if (x<1e-20) { // To guard division by zero later in "fac" - 1e-20 should be justified by checking with Mathematica.
            for (l=0; l<=L; l++) B[l]=0.;
            B[0]=1.;
            return;
        }
        J_bessel_HalfInt(x, nb, B, &ncalc);
        //if (ncalc!=L+1) printf("bessel %d %f\n",ncalc,x); // Check what happens? Zero out the rest?
        double fac = sqrt(PI*.5/x); 
        for (l=0; l<ncalc; l++) B[l]*=fac;
        for (; l<=L; l++) B[l]=0.;
    }

    void Integral_Pack::GenY(double *Y, double X, double phi, int L) {
        int l,m;
        // Plm calculation according to (6.7.9)-(6.7.10) NR 3rd ed
        double Plm[L+1][L+1],sintheta = sqrt(1-X*X);
        Plm[0][0] = 0.5/sqrt(PI);
        for (l=1; l<=L; l++) Plm[l][l]=-sqrt(1.0+0.5/l)*sintheta*Plm[l-1][l-1];
        for (m=0; m<L; m++) Plm[m+1][m]=X*sqrt(2.0*m+3.0)*Plm[m][m];
        for (l=2; l<=L; l++) {
            double ls=l*l, lm1s = (l-1)*(l-1);
            for (m=0; m <= l-2; m++) {
                double ms=m*m;
                Plm[l][m] = sqrt((4.0*ls-1.0)/(ls-ms))*(X*Plm[l-1][m]-sqrt((lm1s-ms)/(4.0*lm1s-1.0))*Plm[l-2][m]);
            }
        }
        // Real Ylm
        for (l=0; l<=L; l++) {
            Y[lm2k(l,0)]  = Plm[l][0];
            for (m=1; m<=l; m++) {
                Y[lm2k(l,m)]  = Plm[l][m]*SQRT2*cos(m*phi);
                Y[lm2k(l,-m)] = Plm[l][m]*SQRT2*sin(m*phi);
            }
        }
    }

    void Integral_Pack::GenclassY(double *A, double *B, double *zetaA, double *zetaB, int dconA, int dconB, int Ln, double *Ylm){
        for (int ii=0; ii<dconA; ii++) for (int jj=0; jj<dconB; jj++) {
            double zeta=zetaA[ii]+zetaB[jj];
            double P[3]={(zetaA[ii]*A[0]+zetaB[jj]*B[0])/zeta,(zetaA[ii]*A[1]+zetaB[jj]*B[1])/zeta,(zetaA[ii]*A[2]+zetaB[jj]*B[2])/zeta};
            double r=sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
            if (r<1e-15/roZ) continue; // CHECK 
            double *Y=&Ylm[(ii*dconB+jj)*(Ln+1)*(Ln+1)],phi=atan2(P[1],P[0]),X=P[2]/r;
            GenY(Y,X,phi,Ln); 
        }
    } 

    void Integral_Pack::Genclass(int a, int b, double *A, double *B, double *zetaA, double *zetaB, double *conA, double *conB, int dconA, int dconB, double* temp, int n, int Ln, double *Ylm, int maxL){ // This function is for a>=b
        if (a<b) { Genclass2(a, b, A, B, zetaA, zetaB, conA, conB, dconA, dconB, temp, n, Ln, Ylm, maxL); return;}
        int K=(Ln+1)*(Ln+1); 
        double ldn=lambda[n];
        bool lzero=(ldn<1e-15/roZ); //  
        double onelambda = lzero? 0.0 : -1.0/ldn;
        double* V1 = arrV; //malloc(totalBraL[a+b+1]*K*sizeof(double));
        double* V2 = arrV+totalBraL[a+b+1]*K;
        double (*HRR[a+b+1][b+1])[K];
        for (int i=a; i<=a+b; i++) {
            if (i==a && b==0) HRR[a][0]=(double (*)[K])temp; 
            else HRR[i][0] = (double (*)[K])(malloc(sizeof(double)*K*noOfBra[i]));
            if (HRR[i][0]==NULL) {printf("Integral_Pack.cc HRR[%d][0] allocation failed sized=%d*sizeof(double)\n",i,K*noOfBra[i]); exit(1);}
            memset(HRR[i][0],0,sizeof(double)*K*noOfBra[i]);
        }
        double J[Ln+a+b+1], *Y=NULL, rAB2=sqr(A[0]-B[0])+sqr(A[1]-B[1])+sqr(A[2]-B[2]);
        if (Ylm==NULL) Y=(double *)malloc(sizeof(double)*K); 
        // printf("A %e %e %e B %e %e %e\n",A[0],A[1],A[2],B[0],B[1],B[2]); 
        for (int ii=0; ii<dconA; ii++) for (int jj=0; jj<dconB; jj++) {
            double zeta=zetaA[ii]+zetaB[jj];
            double P[3]={(zetaA[ii]*A[0]+zetaB[jj]*B[0])/zeta,(zetaA[ii]*A[1]+zetaB[jj]*B[1])/zeta,(zetaA[ii]*A[2]+zetaB[jj]*B[2])/zeta};
            double gAB=exp(-zetaA[ii]*zetaB[jj]/zeta*rAB2)*pow(PI/zeta,1.5)*conA[ii]*conB[jj];
            double one2zeta=.5/zeta;
            double r=sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
            bool rzero=(r<1e-15/roZ); // scaled so that this codition is insensitive to roZ
            //printf("ii=%d jj=%d zetaA=%e zetaB=%e zeta=%e conA=%e conB=%e rAB2=%e gAB=%e\n",ii,jj,zetaA[ii],zetaB[jj],zeta,conA[ii],conB[jj],rAB2,gAB);
            if (!lzero && !rzero) { 
                // if (r<1e-14) printf("small r\n");
                if (Ylm!=NULL)
                    Y=&Ylm[(ii*dconB+jj)*(maxL+1)*(maxL+1)];
                else {                
                    double phi=atan2(P[1],P[0]),X=P[2]/r;
                    GenY(Y,X,phi,Ln);                 
                } // if (Ln<10) printf("IntegralPack ln204 J(%f,%d)\n",r*ldn,Ln+a+b); 
                GenJ(J,r*ldn,Ln+a+b); // if (Ln<10) printf("OK\n");
            } 
            bool swap = false;
            for (int p=a+b; p>=0; p--) {
                swap = !swap;
                double *Va = swap?V2:V1, *Vb = swap?V1:V2;
                memset(Va,0,sizeof(double)*K*totalBraL[a+b+1]);
                // ldn=0.0 is taken care of by separately. (3-term RR & trivial initial conditions) only p=0 contributes! LMRG2013 eqs 19c and 24
                if (lzero && p==0) {
                    Va[0]=JpY00[0]*q[0]*gAB;
                    for (int e=1; e<a+b+1; e++) for (int i=0; i<noOfBra[e]; i++) {
                        int aplusIndex = totalBraL[e]+i;
                        int j=buildMap[aplusIndex]; // printf("j=%d\n",j);
                        int x=inverseMap3[aplusIndex].x,y=inverseMap3[aplusIndex].y,z=inverseMap3[aplusIndex].z;
                        int aIndex = map3[x-delta[0][j]][y-delta[1][j]][z-delta[2][j]];
                        int aminusIndex = map3[abs(x-2*delta[0][j])][abs(y-2*delta[1][j])][abs(z-2*delta[2][j])]; // Be careful
                        int aj = delta[0][j]*(x-1) + delta[1][j]*(y-1) + delta[2][j]*(z-1);
                        Va[Koffset(aplusIndex,0)]=(P[j]-A[j])*Va[Koffset(aIndex,0)];
                        if (aj>0) Va[Koffset(aplusIndex,0)] += aj*one2zeta*Va[Koffset(aminusIndex,0)];
                    }
                }
                // Fill e=0
                if (!rzero && !lzero) {
                    double nfactor=q[n]*gAB*exp(-.25*sqr(ldn)/zeta)*pow(-.5*ldn/zeta/r,p);
                    for (int l=0; l<=Ln; l++) {
                        int ll=l*l+l; 
                        for (int m=-l; m<=l; m++) Va[Koffset(0,ll+m)] = nfactor*J[l+p]*Y[ll+m]; // LHG2012 eq 23 LMRG2013 eq 19a
                    }
                }
                else if (rzero && !lzero) 
                    Va[0] = q[n]*gAB*exp(-.25*sqr(ldn)/zeta)*pow(-.5*sqr(ldn)/zeta,p)*JpY00[p]; // LMRG2013 l=m=0 only - eq 19b  
                // Fill higher e
                if (!lzero) {
                    for (int e=1; e<a+b+1-p; e++) for (int i=0; i<noOfBra[e]; i++)  {
                        int aplusIndex = totalBraL[e]+i;
                        int j=buildMap[aplusIndex]; // printf("j=%d\n",j);
                        int x=inverseMap3[aplusIndex].x,y=inverseMap3[aplusIndex].y,z=inverseMap3[aplusIndex].z;
                        int aIndex = map3[x-delta[0][j]][y-delta[1][j]][z-delta[2][j]];
                        int aminusIndex = map3[abs(x-2*delta[0][j])][abs(y-2*delta[1][j])][abs(z-2*delta[2][j])];
                        int aj = delta[0][j]*(x-1) + delta[1][j]*(y-1) + delta[2][j]*(z-1);
                        double paj = P[j]-A[j];
                        /*cblas_dcopy(K, Va[aIndex], 1, Va[aplusIndex], 1); // Va[aplusIndex][...]=Va[aIndex][...]
                        cblas_dscal(K, paj, Va[aplusIndex], 1); // Va[aplusIndex][...]=paj*Va[aplusIndex][...]
                        cblas_daxpy(K, P[j], Vb[aIndex], 1, Va[aplusIndex], 1); // Va[aplusIndex][...] = P[j]*Vb[aIndex][...]+Va[aplusIndex][...]
                        if (aj>0) {
                            cblas_daxpy(K, aj*one2zeta, Va[aminusIndex], 1, Va[aplusIndex], 1);
                            cblas_daxpy(K, aj*one2zeta, Vb[aminusIndex], 1, Va[aplusIndex], 1);
                        }*/
                        if (aj) for (int k=0; k<K; k++)
                            Va[Koffset(aplusIndex,k)]=P[j]*Vb[Koffset(aIndex,k)]+paj*Va[Koffset(aIndex,k)]+aj*one2zeta*(Va[Koffset(aminusIndex,k)]+Vb[Koffset(aminusIndex,k)]);
                        else for (int k=0; k<K; k++)
                            Va[Koffset(aplusIndex,k)] = P[j]*Vb[Koffset(aIndex,k)]+paj*Va[Koffset(aIndex,k)];

                        switch (j) {
                        case 2: //z
                            for (int l=1; l<=Ln; l++) {
                                int ll=l*l+l; int ll1=l*l-l;
                                for (int m=-l+1; m<l; m++) Va[Koffset(aplusIndex,ll+m)] += onelambda*cz[ll+m]*Vb[Koffset(aIndex,ll1+m)];
                            }
                        break;
                        case 1: //y
                            if (Ln>=1) {
                                Va[Koffset(aplusIndex,1)] += onelambda*cyplus[1]*Vb[Koffset(aIndex,0)];
                                Va[Koffset(aplusIndex,3)] += onelambda*cyminus[3]*Vb[Koffset(aIndex,0)];
                            }
                            for (int l=2; l<=Ln; l++) {
                                int ll=l*l+l; int ll1=l*l-l;
                                Va[Koffset(aplusIndex,ll-l)] += onelambda*cyplus[ll-l]*Vb[Koffset(aIndex,ll1+l-1)];
                                Va[Koffset(aplusIndex,ll+l)] += onelambda*cyminus[ll+l]*Vb[Koffset(aIndex,ll1-l+1)];
                                Va[Koffset(aplusIndex,ll-l+1)] += onelambda*cyplus[ll-l+1]*Vb[Koffset(aIndex,ll1+l-2)];
                                Va[Koffset(aplusIndex,ll+l-1)] += onelambda*cyminus[ll+l-1]*Vb[Koffset(aIndex,ll1-l+2)];
                                for (int m=-l+2; m<l-1; m++) Va[Koffset(aplusIndex,ll+m)] += onelambda*(cyplus[ll+m]*Vb[Koffset(aIndex,ll1-m-1)]+cyminus[ll+m]*Vb[Koffset(aIndex,ll1-m+1)]);
                            }
                        break;
                        //case 0:
                        default: //x
                            if (Ln>=1) {
                                Va[Koffset(aplusIndex,1)] += onelambda*cxplus[1]*Vb[Koffset(aIndex,0)];
                                Va[Koffset(aplusIndex,3)] += onelambda*cxminus[3]*Vb[Koffset(aIndex,0)];
                            }
                            for (int l=2; l<=Ln; l++) {
                                int ll=l*l+l; int ll1=l*l-l;
                                Va[Koffset(aplusIndex,ll-l)] += onelambda*cxplus[ll-l]*Vb[Koffset(aIndex,ll1-l+1)];
                                Va[Koffset(aplusIndex,ll+l)] += onelambda*cxminus[ll+l]*Vb[Koffset(aIndex,ll1+l-1)];
                                Va[Koffset(aplusIndex,ll-l+1)] += onelambda*cxplus[ll-l+1]*Vb[Koffset(aIndex,ll1-l+2)];
                                Va[Koffset(aplusIndex,ll+l-1)] += onelambda*cxminus[ll+l-1]*Vb[Koffset(aIndex,ll1+l-2)];                                
                                for (int m=-l+2; m<l-1; m++) Va[Koffset(aplusIndex,ll+m)] += onelambda*(cxplus[ll+m]*Vb[Koffset(aIndex,ll1+m+1)]+cxminus[ll+m]*Vb[Koffset(aIndex,ll1+m-1)]);

                            }
                        }
                    }
                }

            }
            double* Vx = swap ? V2 : V1;
            for (int i=a; i<=a+b; i++) for (int bra=0; bra<noOfBra[i]; bra++) for (int k=0; k<K; k++)
                HRR[i][0][bra][k] += Vx[Koffset(bra+totalBraL[i],k)]; // cblas_daxpy(K, 1.0, Va[bra+totalBraL[i]], 1, HRR[i][0][bra], 1); 
        }
        if (Ylm==NULL) free(Y);
        double dd[3]={A[0]-B[0],A[1]-B[1],A[2]-B[2]};   
        for (int j=1; j<=b; j++) for (int i=a; i<=a+b-j; i++)  {
            if (i==a && j==b) HRR[a][b]=(double (*)[K])temp;
            else HRR[i][j]=(double (*)[K])(malloc(sizeof(double)*K*noOfBra[i]*noOfBra[j]));
            if (HRR[i][j]==NULL) {printf("Integral_Pack.cc HRR[%d][%d] size=%d*sizeof(double)\n",i,j,K*noOfBra[i]*noOfBra[j]); exit(1);}
            for (int ii=0; ii<noOfBra[i]; ii++) for (int jj=0; jj<noOfBra[j]; jj++) {
            	int lindex = ii*noOfBra[j] + jj;
            	int rindex1=HRRMAP[i][j][lindex].x;
            	int rindex2=HRRMAP[i][j][lindex].y;
            	double factor = dd[HRRMAP[i][j][lindex].z];
                for (int k=0; k<K; k++) HRR[i][j][lindex][k]=factor*HRR[i][j-1][rindex2][k]+HRR[i+1][j-1][rindex1][k];
                /*double* lhs = HRR[i][j][lindex], rhs1 = HRR[i+1][j-1][rindex1], rhs2 = HRR[i][j-1][rindex2];
                // lhs[...]=factor*rhs2[...]+rhs1[...]
                cblas_dcopy(K, rhs1, 1, lhs, 1); // lhs[...] = rhs1[...]
                cblas_daxpy(K, factor, rhs2, 1, lhs, 1); // lhs[...] = factor*rhs2[...] + lhs[...]*/
            }
            free(HRR[i][j-1]);
            if (i==a+b-j) free(HRR[i+1][j-1]);
        }
        //memcpy(temp, HRR[a][b], noOfBra[a]*noOfBra[b]*K*sizeof(double)); free(HRR[a][b]);
    }

    void Integral_Pack::Genclass2(int a, int b, double *A, double *B, double *zetaA, double *zetaB, double *conA, double *conB, int dconA, int dconB, double* temp, int n, int Ln, double *Ylm, int maxL){ //a<b
        int K=(Ln+1)*(Ln+1); 
        double ldn=lambda[n];
        bool lzero=(ldn<1e-15/roZ); //  
        double onelambda = lzero? 0.0 : -1.0/ldn;
        double* V1 = arrV; 
        double* V2 = arrV+totalBraL[a+b+1]*K;
        double (*HRR[a+1][a+b+1])[K]; /*a<b*/
        for (int i=b; i<=a+b; i++) { /*a<b*/
            if (i==b && a==0) HRR[0][b]=(double (*)[K])temp; /*a<b*/
            else HRR[0][i] = (double (*)[K])(malloc(sizeof(double)*K*noOfBra[i])); /*a<b*/
            if (HRR[0][i]==NULL) {printf("Integral_Pack.cc HRR[0][%d] allocation failed sized=%d*sizeof(double)\n",i,K*noOfBra[i]); exit(1);} /*a<b*/
            memset(HRR[0][i],0,sizeof(double)*K*noOfBra[i]);
        }
        double J[Ln+a+b+1], *Y=NULL, rAB2=sqr(A[0]-B[0])+sqr(A[1]-B[1])+sqr(A[2]-B[2]);
        if (Ylm==NULL) Y=(double *)malloc(sizeof(double)*K); 
        for (int ii=0; ii<dconA; ii++) for (int jj=0; jj<dconB; jj++) {
            double zeta=zetaA[ii]+zetaB[jj];
            double P[3]={(zetaA[ii]*A[0]+zetaB[jj]*B[0])/zeta,(zetaA[ii]*A[1]+zetaB[jj]*B[1])/zeta,(zetaA[ii]*A[2]+zetaB[jj]*B[2])/zeta};
            double gAB=exp(-zetaA[ii]*zetaB[jj]/zeta*rAB2)*pow(PI/zeta,1.5)*conA[ii]*conB[jj];
            double one2zeta=.5/zeta;
            double r=sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
            bool rzero=(r<1e-15/roZ); 
            if (!lzero && !rzero) { 
                if (Ylm!=NULL)
                    Y=&Ylm[(ii*dconB+jj)*(maxL+1)*(maxL+1)];
                else {                
                    double phi=atan2(P[1],P[0]),X=P[2]/r;
                    GenY(Y,X,phi,Ln);                 
                } 
                GenJ(J,r*ldn,Ln+a+b); 
            } 
            bool swap = false;
            for (int p=a+b; p>=0; p--) {
                swap = !swap;
                double *Va = swap?V2:V1, *Vb = swap?V1:V2;
                memset(Va,0,sizeof(double)*K*totalBraL[a+b+1]);
                if (lzero && p==0) {
                    Va[0]=JpY00[0]*q[0]*gAB;
                    for (int e=1; e<a+b+1; e++) for (int i=0; i<noOfBra[e]; i++) {
                        int aplusIndex = totalBraL[e]+i;
                        int j=buildMap[aplusIndex]; 
                        int x=inverseMap3[aplusIndex].x,y=inverseMap3[aplusIndex].y,z=inverseMap3[aplusIndex].z;
                        int aIndex = map3[x-delta[0][j]][y-delta[1][j]][z-delta[2][j]];
                        int aminusIndex = map3[abs(x-2*delta[0][j])][abs(y-2*delta[1][j])][abs(z-2*delta[2][j])];
                        int aj = delta[0][j]*(x-1) + delta[1][j]*(y-1) + delta[2][j]*(z-1);
                        Va[Koffset(aplusIndex,0)]=(P[j]-B[j])*Va[Koffset(aIndex,0)]; /*a<b*/
                        if (aj>0) Va[Koffset(aplusIndex,0)] += aj*one2zeta*Va[Koffset(aminusIndex,0)];
                    }
                }
                if (!rzero && !lzero) {
                    double nfactor=q[n]*gAB*exp(-.25*sqr(ldn)/zeta)*pow(-.5*ldn/zeta/r,p);
                    for (int l=0; l<=Ln; l++) {
                        int ll=l*l+l; 
                        for (int m=-l; m<=l; m++) Va[Koffset(0,ll+m)] = nfactor*J[l+p]*Y[ll+m];
                    }
                }
                else if (rzero && !lzero) 
                    Va[0] = q[n]*gAB*exp(-.25*sqr(ldn)/zeta)*pow(-.5*sqr(ldn)/zeta,p)*JpY00[p]; 
                if (!lzero) {
                    for (int e=1; e<a+b+1-p; e++) for (int i=0; i<noOfBra[e]; i++)  {
                        int aplusIndex = totalBraL[e]+i;
                        int j=buildMap[aplusIndex]; 
                        int x=inverseMap3[aplusIndex].x,y=inverseMap3[aplusIndex].y,z=inverseMap3[aplusIndex].z;
                        int aIndex = map3[x-delta[0][j]][y-delta[1][j]][z-delta[2][j]];
                        int aminusIndex = map3[abs(x-2*delta[0][j])][abs(y-2*delta[1][j])][abs(z-2*delta[2][j])];
                        int aj = delta[0][j]*(x-1) + delta[1][j]*(y-1) + delta[2][j]*(z-1);
                        double pbj = P[j]-B[j]; /*a<b*/
                        if (aj) for (int k=0; k<K; k++)
                            Va[Koffset(aplusIndex,k)]=P[j]*Vb[Koffset(aIndex,k)]+pbj*Va[Koffset(aIndex,k)]+aj*one2zeta*(Va[Koffset(aminusIndex,k)]+Vb[Koffset(aminusIndex,k)]); /*a<b*/
                        else for (int k=0; k<K; k++)
                            Va[Koffset(aplusIndex,k)] = P[j]*Vb[Koffset(aIndex,k)]+pbj*Va[Koffset(aIndex,k)]; /*a<b*/
                        switch (j) {
                        case 2: 
                            for (int l=1; l<=Ln; l++) {
                                int ll=l*l+l; int ll1=l*l-l;
                                for (int m=-l+1; m<l; m++) Va[Koffset(aplusIndex,ll+m)] += onelambda*cz[ll+m]*Vb[Koffset(aIndex,ll1+m)];
                            }
                        break;
                        case 1: 
                            if (Ln>=1) {
                                Va[Koffset(aplusIndex,1)] += onelambda*cyplus[1]*Vb[Koffset(aIndex,0)];
                                Va[Koffset(aplusIndex,3)] += onelambda*cyminus[3]*Vb[Koffset(aIndex,0)];
                            }
                            for (int l=2; l<=Ln; l++) {
                                int ll=l*l+l; int ll1=l*l-l;
                                Va[Koffset(aplusIndex,ll-l)] += onelambda*cyplus[ll-l]*Vb[Koffset(aIndex,ll1+l-1)];
                                Va[Koffset(aplusIndex,ll+l)] += onelambda*cyminus[ll+l]*Vb[Koffset(aIndex,ll1-l+1)];
                                Va[Koffset(aplusIndex,ll-l+1)] += onelambda*cyplus[ll-l+1]*Vb[Koffset(aIndex,ll1+l-2)];
                                Va[Koffset(aplusIndex,ll+l-1)] += onelambda*cyminus[ll+l-1]*Vb[Koffset(aIndex,ll1-l+2)];
                                for (int m=-l+2; m<l-1; m++) Va[Koffset(aplusIndex,ll+m)] += onelambda*(cyplus[ll+m]*Vb[Koffset(aIndex,ll1-m-1)]+cyminus[ll+m]*Vb[Koffset(aIndex,ll1-m+1)]);
                            }
                        break;
                        default: 
                            if (Ln>=1) {
                                Va[Koffset(aplusIndex,1)] += onelambda*cxplus[1]*Vb[Koffset(aIndex,0)];
                                Va[Koffset(aplusIndex,3)] += onelambda*cxminus[3]*Vb[Koffset(aIndex,0)];
                            }
                            for (int l=2; l<=Ln; l++) {
                                int ll=l*l+l; int ll1=l*l-l;
                                Va[Koffset(aplusIndex,ll-l)] += onelambda*cxplus[ll-l]*Vb[Koffset(aIndex,ll1-l+1)];
                                Va[Koffset(aplusIndex,ll+l)] += onelambda*cxminus[ll+l]*Vb[Koffset(aIndex,ll1+l-1)];
                                Va[Koffset(aplusIndex,ll-l+1)] += onelambda*cxplus[ll-l+1]*Vb[Koffset(aIndex,ll1-l+2)];
                                Va[Koffset(aplusIndex,ll+l-1)] += onelambda*cxminus[ll+l-1]*Vb[Koffset(aIndex,ll1+l-2)];                                
                                for (int m=-l+2; m<l-1; m++) Va[Koffset(aplusIndex,ll+m)] += onelambda*(cxplus[ll+m]*Vb[Koffset(aIndex,ll1+m+1)]+cxminus[ll+m]*Vb[Koffset(aIndex,ll1+m-1)]);

                            }
                        }
                    }
                }

            }
            double* Vx = swap ? V2 : V1;
            for (int i=b; i<=a+b; i++) for (int bra=0; bra<noOfBra[i]; bra++) for (int k=0; k<K; k++) /*a<b*/
                HRR[0][i][bra][k] += Vx[Koffset(bra+totalBraL[i],k)]; /*a<b*/
        }
        if (Ylm==NULL) free(Y);
        double dd[3]={B[0]-A[0], B[1]-A[1], B[2]-A[2]}; /*a<b*/
        for (int i=1; i<=a; i++) for (int j=b; j<=a+b-i; j++) { /*a<b*/
            if (i==a && j==b) HRR[a][b]=(double (*)[K])temp;
            else HRR[i][j]=(double (*)[K])(malloc(sizeof(double)*K*noOfBra[i]*noOfBra[j]));
            if (HRR[i][j]==NULL) {printf("Integral_Pack.cc HRR[%d][%d] size=%d*sizeof(double)\n",i,j,K*noOfBra[i]*noOfBra[j]); exit(1);}
            for (int ii=0; ii<noOfBra[i]; ii++) for (int jj=0; jj<noOfBra[j]; jj++) {
            	int lindex = ii*noOfBra[j] + jj; /*OK*/
            	int rindex1=HRRMAP2[i][j][lindex].x; /*OK*/
            	int rindex2=HRRMAP2[i][j][lindex].y; /*OK*/
            	double factor = dd[HRRMAP2[i][j][lindex].z]; /*a<b*/
                for (int k=0; k<K; k++) HRR[i][j][lindex][k]=factor*HRR[i-1][j][rindex2][k]+HRR[i-1][j+1][rindex1][k]; /*a<b*/
            }
            free(HRR[i-1][j]); /*a<b*/
            if (i==a+b-j) free(HRR[i-1][j+1]); /*a<b*/
        }
    }

    void Integral_Pack::initializeCoulomb(int N){
        lambda = (double *) malloc(sizeof(double)*(N+1));
        q = (double *) malloc(sizeof(double)*(N+1));
        int i;
        for (i=0; i<=N; i++) {
            lambda[i]=i;
            q[i]=2.*sqrt(2.);
        }
        q[0]=2.;
    }

    void Integral_Pack::initializeEwald(int N, int L, double Omega, double thresh, double rad){
        // TODO: Take into account eqn11 and eqn12 in RO#5
        //int Ncal,Nprime;
        double R=rad*Omega*2.,th=thresh*PI/4/Omega;
        Ncal=(int) ceil(R*R/4.+(sqrt(-log10(th))-1)*R+2.); //RO Thesis Eq (5.11) and RO#5 Eq (11)
        if (Ncal>N) Ncal=N; // Limit by user input (if any)
        Nprime = (int) ceil(2./PI*sqrt(-(Ncal+1)*log(th))-1.); // RO Thesis Eq (5.12) and RO#5 Eq (12)
        if (Nprime>Ncal) Nprime=Ncal;        

        lambda = (double *) malloc(sizeof(double)*(Nprime+1));
        q = (double *) malloc(sizeof(double)*(Nprime+1));
        if (!lambda || !q) {printf("Allocation failed ln393 IntPack\n"); exit(1);}
        int n;

        // TODO: Replaced external file read-in by on-the-fly generation of roots and weights        
        FILE *fptr1,*fptr2;
        char fname1[255],fname2[255];
        sprintf(fname1,"Ewald/roots%d.txt",Ncal);
        sprintf(fname2,"Ewald/weights%d.txt",Ncal);
        fptr1=(FILE *)fopen(fname1,"r");
        fptr2=(FILE *)fopen(fname2,"r");
        if (!fptr1 || !fptr2) {printf("Integral_Pack.cc can't find Hermite root/weight for Ewald calculation.\n"); exit(1);}
        for (n=0; n<=Nprime; n++) {
            if (fscanf(fptr1,"%lf",&lambda[n]) != 1) {
                fprintf(stderr, "error: invalid lambda[%d]\n", n);
                exit(1);
            }
            lambda[n]*=2.*Omega;
            if (fscanf(fptr2,"%lf",&q[n]) != 1) {
                fprintf(stderr, "error: invalid q[%d]\n", n);
                exit(1);
            }
            q[n]=4.*sqrt(q[n]*Omega);
        }
        fclose(fptr1); fclose(fptr2);

    }

    void Integral_Pack::getNL(int *n_l) {
        printf("*** Integral_Pack::getNL ****\n");
        int n,l,maxl=-1;
        double th=thresh/Nprime; // this thresh has been scaled
        //printf("th=%e thresh=%e roZ=%e\n",th,thresh,roZ);
        for (n=0; n<=Nprime; n++) {
            double J[L+1]; //printf("lambda[%d]*rad=%e lambda=%e\n",n,lambda[n]*rad,lambda[n]);
            GenJ(J,lambda[n]*rad,L); // 2 omega beta[n] r1           
            for (l=L; l>=0 && fabs(sqr(q[n]*J[l])/4./PI)<th;) {/*printf("%e\n",fabs(sqr(q[n]*J[l])/4./PI));*/ l--;} //q[n] and th has roZ factor     
            
            printf("n=%d Ltest=%d\n",n,l);
            n_l[n+1]=l;
            if (l>maxl) maxl=l;
        }
        n_l[0]=Nprime;
        n_l[Nprime+2]=maxl;
        printf("*** Override N and L ***\n");
        printf("Ncal=%d Nprime=%d L=%d\n", Ncal,Nprime,maxl);
    }

// end Integral_Pack

                }
            }
        }
    }
}

/*
int main() {

    // Test Genclass
    int a=2,b=0,N=10,L=10,dconA=6,dconB=3;
    double A[3],B[3];
    double conA[dconA],conB[dconB],zetaA[dconA],zetaB[dconB];

    au::edu::anu::qm::ro::Integral_Pack* ip = new au::edu::anu::qm::ro::Integral_Pack(N, L, 0.3, 1e-6, 2.0, 1.0);

    A[0]=.1; B[0]=.4;
    A[1]=.2; B[1]=.3;
    A[2]=.3; B[2]=.2;

    zetaA[0]=1.;  conA[0]=6.;
    zetaA[1]=2.;  conA[1]=5.;
    zetaA[2]=3.;  conA[2]=4.;
    zetaA[3]=4.;  conA[3]=3.;
    zetaA[4]=5.;  conA[4]=2.;
    zetaA[5]=6.;  conA[5]=1.;

    zetaB[0]=1.;  conB[0]=1.5;
    zetaB[1]=2.;  conB[1]=2.5;
    zetaB[2]=4.5; conB[2]=3.5;

    double temp[(L+1)*(L+1)*(a+1)*(a+2)/2*(b+1)*(b+2)/2],Y[(L+1)*(L+1)*dconA*dconB];
    ip->GenclassY(A,B,zetaA,zetaB,dconA,dconB,L,Y);
    for (int i=0; i<100000; i++)
        ip->Genclass(a,b,A,B,zetaA,zetaB,conA,conB,dconA,dconB,temp,N,L,Y,L);

    delete ip;
}
*/
