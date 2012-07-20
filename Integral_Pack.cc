#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "bessel4.h"
#include "Integral_Pack.h"
#include "cblas.h"

#define sqr(x) ((x)*(x))
#define SQRT2 1.4142135623730951

int delta[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double JpY00[11]={0.28209479177387814, 0.09403159725795937, 0.018806319451591877,
       		0.0026866170645131254, 0.000298513007168125,0.000027137546106193183, 
                2.087503546630245e-6,1.39166903108683e-7, 8.186288418157823e-9,
                4.308572851662012e-10, 2.0517013579342914e-11}; 
double *arrV;

// 0-10 accomodate (hh|phi) - need more for higer angular momentum - see RO#7 eqn 6b 
// double *YCon[MAX_C*MAX_C];

namespace au {
    namespace edu {
        namespace anu {
            namespace qm {
                namespace ro {

    Integral_Pack* Integral_Pack::_make(int N, int L) {
      return new Integral_Pack(N,L);
    }

    Integral_Pack::Integral_Pack(int N,int L) {
        this->N = N; this->L=L;
        initialize();
        initializeCoulomb(N);

        arrV=(double *)malloc(totalBraL[MAX_BRA_L+1]*(L+1)*(L+1)*sizeof(double)*2);
        printf("Integral_Pack.cc N=%d L=%d\n",N,L);
    }

    Integral_Pack::~Integral_Pack() {
        free(lambda); free(q); free(arrV); 
        for (int i=1; i<=MAX_BRA_L; i++) for (int j=1; j<=i; j++) {
            free(HRRMAP[i][j]);
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

                // LHG 2012 Table II & boolean before (36)
                if (z == 1) buildMap[m]=2; //z
                else if (y == 1) buildMap[m]=1; //y
                else if (x == 1) buildMap[m]=0; //x
                else if (z >= 2) buildMap[m]=2; //z
                else if (y >= 2) buildMap[m]=1; //y
                else buildMap[m]=0; //x

                //printf("%d : %d %d %d j=%d\n", m, x,y,z, buildMap[m]);
                m++;
            }
        }

        // HRR
        int i,j,a,b;
        for (i=1; i<=MAX_BRA_L; i++) for (j=1; j<=i; j++) {
            HRRMAP[i][j] = (Point *) malloc(sizeof(Point)*noOfBra[i]*noOfBra[j]);
            if (HRRMAP[i][j]==NULL) {printf("Integral_Pack.cc malloc failed at ln75\n"); exit(1);}
            for (a=0; a<noOfBra[i]; a++) for (b=0; b<noOfBra[j]; b++) {
             	int aInt=totalBraL[i]+a,bInt=totalBraL[j]+b;
             	int ax=inverseMap3[aInt].x, ay=inverseMap3[aInt].y, az=inverseMap3[aInt].z,
             		bx=inverseMap3[bInt].x,	by=inverseMap3[bInt].y,	bz=inverseMap3[bInt].z;

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

        //printf("%d\n",MAX_TOTAL_BRA_L);
        // LMR 2012 manuscript in preparation - similar to LHG 2012 eqn (27)
        for (l=0; l<=MAX_KET_L; l++) for (m=-l; m<=l; m++) {
            double cminus=.5*sqrt((l-m-1.0)*(l-m)*(2.0*l+1.0)/(2.0*l-1.0)),
                cplus=.5*sqrt((l+m-1.0)*(l+m)*(2.0*l+1.0)/(2.0*l-1.0));
            int index=lm2k(l,m);

            if (m<=-2) cxplus[index]=-cminus;
                else if (m>=1) cxplus[index]=cminus;
                else if (m==-1) cxplus[index]=0;
                else /*if (m==0)*/ cxplus[index]=SQRT2*cminus;

            if (m<=-1) cxminus[index]=cplus;
                else if (m>=2) cxminus[index]=-cplus;
                else if (m==0) cxminus[index]=0;
                else /*if (m==1)*/ cxminus[index]=-SQRT2*cplus;

            if (m<=-2) cyplus[index] = -cminus;
                else if (m>=1) cyplus[index] = cminus;
                else if (m==-1) cyplus[index] = -SQRT2*cminus;
                else /*if (m==0)*/ cyplus[index] = SQRT2*cminus;

            if (m<=-1) cyminus[index]=-cplus;
                else if (m>=2) cyminus[index]=cplus;
                else /*if (m==0)*/ cyminus[index]=0;
                //else if (m==1) cyminus[index]=0;

            cz[index] = sqrt((l*l-m*m)*(2.0*l+1.0)/(2.0*l-1.0));
            //printf("%2d %2d : %15.6e %15.6e %15.6e %15.6e %15.6e : %15.6e %15.6e %d\n",l,m,cxplus[index],cxminus[index],cyplus[index],cyminus[index],cz[index],cplus,cminus,index);
        }

/*        for (i=0; i<MAX_C*MAX_C; i++) {
            YCon[i]=(double *)malloc((L+1)*(L+1)*sizeof(double));
            if (YCon[i]==NULL) {printf("Integral_Pack.cc allocation failed at ln128\n"); exit(1);}
        }*/
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
        // Plm calculation according to (6.7.9)-(6.7.10) NR 3rd ed
        int l,m;
        double Plm[L+1][L+1],sintheta = sqrt(1-X*X);

        Plm[0][0] = 0.5/sqrt(PI);
        for (l=1;l<=L; l++)
            Plm[l][l]=-sqrt(1.0+0.5/l)*sintheta*Plm[l-1][l-1];

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
        int ii,jj;
        for (ii=0; ii<dconA; ii++) for (jj=0; jj<dconB; jj++) {
            double zeta=zetaA[ii]+zetaB[jj];
            double P[3]={(zetaA[ii]*A[0]+zetaB[jj]*B[0])/zeta,(zetaA[ii]*A[1]+zetaB[jj]*B[1])/zeta,(zetaA[ii]*A[2]+zetaB[jj]*B[2])/zeta};
            double r=sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
            if (r<=1e-14) continue; // CHECK 
            double *Y=&Ylm[(ii*dconB+jj)*(Ln+1)*(Ln+1)],phi=atan2(P[1],P[0]),X=P[2]/r;
            GenY(Y,X,phi,Ln); 
        }
    } 

    void Integral_Pack::Genclass(int a, int b, double *A, double *B, double *zetaA, double *zetaB, double *conA, double *conB, int dconA, int dconB, double* temp, int n, int Ln, double *Ylm, int maxL){
        int bra,K=(Ln+1)*(Ln+1),p,e,i,ii,j,jj,k,l,m,ll,ll1; 
        double ldn=lambda[n],onelambda;
        bool swap,lzero=(ldn==0.); // should be replace by ldn<1e-y 
        onelambda = lzero? 0.0 : -1.0/ldn;
        double (*V1)[K]=(double (*)[K])arrV;///malloc(totalBraL[a+b+1]*K*sizeof(double));
        double (*V2)[K]=(double (*)[K])(arrV+totalBraL[a+b+1]*K);//malloc(totalBraL[a+b+1]*K*sizeof(double));
        //if (V1==NULL || V2==NULL) {printf("Integral_Pack.cc V1/V2 allocation failed size=%d*sizeof(double)\n",totalBraL[a+b+1]*K); exit(1);}
        double (*HRR[a+b+1][b+1])[K];
        for (i=a; i<=a+b; i++) {
            if (i==a && b==0) HRR[a][0]=(double (*)[K])temp; 
            else HRR[i][0] = (double (*)[K])(malloc(sizeof(double)*K*noOfBra[i]));
            if (HRR[i][0]==NULL) {printf("Integral_Pack.cc HRR[%d][0] allocation failed sized=%d*sizeof(double)\n",i,K*noOfBra[i]); exit(1);}
            memset(HRR[i][0],0,sizeof(double)*K*noOfBra[i]);
        }
        double J[Ln+a+b+1], *Y=NULL, rAB2=sqr(A[0]-B[0])+sqr(A[1]-B[1])+sqr(A[2]-B[2]);
        if (Ylm==NULL) Y=(double *)malloc(sizeof(double)*K); 
        // printf("A %e %e %e B %e %e %e\n",A[0],A[1],A[2],B[0],B[1],B[2]); 
        for (ii=0; ii<dconA; ii++) for (jj=0; jj<dconB; jj++) {
            double zeta=zetaA[ii]+zetaB[jj];
            double P[3]={(zetaA[ii]*A[0]+zetaB[jj]*B[0])/zeta,(zetaA[ii]*A[1]+zetaB[jj]*B[1])/zeta,(zetaA[ii]*A[2]+zetaB[jj]*B[2])/zeta};
            double gAB=exp(-zetaA[ii]*zetaB[jj]/zeta*rAB2)*pow(PI/zeta,1.5)*conA[ii]*conB[jj];
            double one2zeta=.5/zeta;

            double r=sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
            bool rzero=(r<=1e-14); // CHECK

            //printf("ii=%d jj=%d zetaA=%e zetaB=%e zeta=%e conA=%e conB=%e rAB2=%e gAB=%e\n",ii,jj,zetaA[ii],zetaB[jj],zeta,conA[ii],conB[jj],rAB2,gAB);

            if (!lzero && !rzero) { 
                // if (r<1e-14) printf("small r\n");
                if (Ylm!=NULL)
                    Y=&Ylm[(ii*dconB+jj)*(maxL+1)*(maxL+1)];
                else {                
                    double phi=atan2(P[1],P[0]),X=P[2]/r;
                    GenY(Y,X,phi,Ln);                 
                }
                // if (Ln<10) printf("IntegralPack ln204 J(%f,%d)\n",r*ldn,Ln+a+b); 
                GenJ(J,r*ldn,Ln+a+b); 
                //if (Ln<10) printf("OK\n");
            } 

            swap = false;
            for (p=a+b; p>=0; p--) {
                swap = !swap;
                double (* Va)[K] = swap ? V2 : V1;
                double (* Vb)[K] = swap ? V1 : V2;
                memset(Va,0,sizeof(double)*K*totalBraL[a+b+1]);

                // ldn=0.0 is taken care of by separately. (3-term RR & trivial initial conditions) only p=0 contributes! RO#7 eqn11
                if (lzero && p==0) {
                    Va[0][0]=JpY00[0]*q[0]*gAB;
                    for (e=1; e<a+b+1; e++) for (i=0; i<noOfBra[e]; i++) {
                        int aplusIndex = totalBraL[e]+i;
                        int j=buildMap[aplusIndex]; // printf("j=%d\n",j);
                        int x=inverseMap3[aplusIndex].x,y=inverseMap3[aplusIndex].y,z=inverseMap3[aplusIndex].z;
                        int aIndex = map3[x-delta[0][j]][y-delta[1][j]][z-delta[2][j]];
                        int aminusIndex = map3[abs(x-2*delta[0][j])][abs(y-2*delta[1][j])][abs(z-2*delta[2][j])]; // Be careful
                        int aj = delta[0][j]*(x-1) + delta[1][j]*(y-1) + delta[2][j]*(z-1);
                        Va[aplusIndex][0]=(P[j]-A[j])*Va[aIndex][0];
                        if (aj>0) Va[aplusIndex][0] += aj*one2zeta*Va[aminusIndex][0];
                    }
                }

                // Fill e=0
                if (!rzero && !lzero) {
                    double nfactor=q[n]*gAB*exp(-.25*sqr(ldn)/zeta)*pow(-.5*ldn/zeta/r,p);
                    for (l=0; l<=Ln; l++) {
                        ll=l*l+l; for (m=-l; m<=l; m++) Va[0][ll+m] = nfactor*J[l+p]*Y[ll+m]; // eqn (23) //6a
                    }
                }
                else if (rzero && !lzero) 
                    Va[0][0] = q[n]*gAB*exp(-.25*sqr(ldn)/zeta)*pow(-.5*sqr(ldn)/zeta,p)*JpY00[p]; // l=m=0 only //6b   

                // Fill higher e
                if (!lzero) {
                    for (e=1; e<a+b+1-p; e++) for (i=0; i<noOfBra[e]; i++)  {
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

                        if (aj) for (k=0; k<K; k++)
                            Va[aplusIndex][k]=P[j]*Vb[aIndex][k]+paj*Va[aIndex][k]+aj*one2zeta*(Va[aminusIndex][k]+Vb[aminusIndex][k]);
                        else for (k=0; k<K; k++)
                            Va[aplusIndex][k] = P[j]*Vb[aIndex][k]+paj*Va[aIndex][k];

                        switch (j) {
                        case 2: //z
                            for (l=1; l<=Ln; l++) {
                                ll=l*l+l; ll1=l*l-l;
                                for (m=-l+1; m<l; m++) Va[aplusIndex][ll+m] += onelambda*cz[ll+m]*Vb[aIndex][ll1+m];
                            }
                        break;
                        case 1: //y
                            if (Ln>=1) {
                                Va[aplusIndex][1] += onelambda*cyplus[1]*Vb[aIndex][0];
                                Va[aplusIndex][3] += onelambda*cyminus[3]*Vb[aIndex][0];
                            }
                            for (l=2; l<=Ln; l++) {
                                ll=l*l+l; ll1=l*l-l;
                                Va[aplusIndex][ll-l] += onelambda*cyplus[ll-l]*Vb[aIndex][ll1+l-1];
                                Va[aplusIndex][ll+l] += onelambda*cyminus[ll+l]*Vb[aIndex][ll1-l+1];
                                Va[aplusIndex][ll-l+1] += onelambda*cyplus[ll-l+1]*Vb[aIndex][ll1+l-2];
                                Va[aplusIndex][ll+l-1] += onelambda*cyminus[ll+l-1]*Vb[aIndex][ll1-l+2];
                                for (m=-l+2; m<l-1; m++) 
                                    Va[aplusIndex][ll+m] += onelambda*(cyplus[ll+m]*Vb[aIndex][ll1-m-1]+cyminus[ll+m]*Vb[aIndex][ll1-m+1]);
                            }
                        break;
                        //case 0:
                        default: //x
                            if (Ln>=1) {
                                Va[aplusIndex][1] += onelambda*cxplus[1]*Vb[aIndex][0];
                                Va[aplusIndex][3] += onelambda*cxminus[3]*Vb[aIndex][0];
                            }
                            for (l=2; l<=Ln; l++) {
                                ll=l*l+l; ll1=l*l-l;
                                Va[aplusIndex][ll-l] += onelambda*cxplus[ll-l]*Vb[aIndex][ll1-l+1];
                                Va[aplusIndex][ll+l] += onelambda*cxminus[ll+l]*Vb[aIndex][ll1+l-1];
                                Va[aplusIndex][ll-l+1] += onelambda*cxplus[ll-l+1]*Vb[aIndex][ll1-l+2];
                                Va[aplusIndex][ll+l-1] += onelambda*cxminus[ll+l-1]*Vb[aIndex][ll1+l-2];                                
                                for (m=-l+2; m<l-1; m++) 
                                    Va[aplusIndex][ll+m] += onelambda*(cxplus[ll+m]*Vb[aIndex][ll1+m+1]+cxminus[ll+m]*Vb[aIndex][ll1+m-1]);

                            }
                        }
                    }
                }

            }
            double (* Va)[K] = swap ? V2 : V1;
            for (i=a; i<=a+b; i++) for (bra=0; bra<noOfBra[i]; bra++) for (k=0; k<K; k++)
                HRR[i][0][bra][k] = Va[bra+totalBraL[i]][k]+HRR[i][0][bra][k]; 
                // cblas_daxpy(K, 1.0, Va[bra+totalBraL[i]], 1, HRR[i][0][bra], 1); 
                // HRR[i][0][bra][...] = Va[bra+totalBraL[i]][...]+HRR[i][0][bra][...] 
        }

        //free(V1);free(V2);
        if (Ylm==NULL) free(Y);

        double dd[3]={A[0]-B[0],A[1]-B[1],A[2]-B[2]};
   
        for (j=1; j<=b; j++) for (i=a; i<=a+b-j; i++)  {
            if (i==a && j==b) HRR[a][b]=(double (*)[K])temp;
            else HRR[i][j]=(double (*)[K])(malloc(sizeof(double)*K*noOfBra[i]*noOfBra[j]));
            if (HRR[i][j]==NULL) {printf("Integral_Pack.cc HRR[%d][%d] size=%d*sizeof(double)\n",i,j,K*noOfBra[i]*noOfBra[j]); exit(1);}
            for (ii=0; ii<noOfBra[i]; ii++) for (jj=0; jj<noOfBra[j]; jj++) {
            	int lindex = ii*noOfBra[j] + jj;
            	int rindex1=HRRMAP[i][j][lindex].x;
            	int rindex2=HRRMAP[i][j][lindex].y;
            	double factor = dd[HRRMAP[i][j][lindex].z];
                for (k=0; k<K; k++) HRR[i][j][lindex][k]=factor*HRR[i][j-1][rindex2][k]+HRR[i+1][j-1][rindex1][k];
                /*double* lhs = HRR[i][j][lindex];
                double* rhs1 = HRR[i+1][j-1][rindex1];
                double* rhs2 = HRR[i][j-1][rindex2];
                // lhs[...]=factor*rhs2[...]+rhs1[...]
                cblas_dcopy(K, rhs1, 1, lhs, 1); // lhs[...] = rhs1[...]
                cblas_daxpy(K, factor, rhs2, 1, lhs, 1); // lhs[...] = factor*rhs2[...] + lhs[...]*/
            }
            free(HRR[i][j-1]);
            if (i==a+b-j) free(HRR[i+1][j-1]);
        }
        //memcpy(temp, HRR[a][b], noOfBra[a]*noOfBra[b]*K*sizeof(double));
        //free(HRR[a][b]);
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

// end Integral_Pack

                }
            }
        }
    }
}

