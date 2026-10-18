#include "mex.h"
#include <math.h>
#include <float.h>
/* Optional implementation of the identical real Biot-Savart kernels.
 * Mode 0: independent finite segments. Mode 1: derivative of the shared
 * vertex between two collinear segments. MATLAB implementations remain the
 * reference and handle complex-step evaluations. No time history is cut. */
static double dot(const double *a,const double *b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static void cross(const double *a,const double *b,double *c) {c[0]=a[1]*b[2]-a[2]*b[1];c[1]=a[2]*b[0]-a[0]*b[2];c[2]=a[0]*b[1]-a[1]*b[0];}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    if(nrhs!=5) mexErrMsgIdAndTxt("uvlm:Native","Five arguments required.");
    int mode=(int)mxGetScalar(prhs[0]);
    mwSize np=mxGetM(prhs[1]),ns=mxGetM(prhs[2]);
    if((mode!=0 && mode!=1)||nlhs!=(mode==0?1:3)||mxGetN(prhs[1])!=3||mxGetN(prhs[2])!=3||mxGetN(prhs[3])!=3||mxGetM(prhs[3])!=(mode==0?ns:1)) mexErrMsgIdAndTxt("uvlm:Native","Invalid kernel dimensions.");
    for(int k=1;k<5;k++) if(!mxIsDouble(prhs[k])||mxIsComplex(prhs[k])) mexErrMsgIdAndTxt("uvlm:Native","Real double inputs required.");
    const double *P=mxGetPr(prhs[1]),*A=mxGetPr(prhs[2]),*B=mxGetPr(prhs[3]);
    double core=mxGetScalar(prhs[4]),e2=core*core,cut=core*1e-3,inv4pi=1/(4*acos(-1.0));
    mwSize dims[3]={np,3,ns};double *out[3];
    for(int k=0;k<nlhs;k++){plhs[k]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);out[k]=mxGetPr(plhs[k]);}
    for(mwSize j=0;j<ns;j++){
        double a[3],d[3];for(int k=0;k<3;k++){a[k]=A[j+ns*k];d[k]=mode==0?B[j+ns*k]-a[k]:B[k];}
        double dd=dot(d,d);
        for(mwSize i=0;i<np;i++){
            double r[3],rp[3],rm[3],c[3],n0,nplus,nminus,D;
            for(int k=0;k<3;k++){r[k]=P[i+np*k]-a[k];rp[k]=r[k]+d[k];rm[k]=r[k]-d[k];}
            n0=sqrt(dot(r,r));nminus=sqrt(dot(rm,rm));
            if(mode==0){
                cross(r,rm,c);D=dot(c,c)+e2*dd;
                if(n0<=cut||nminus<=cut||D<=DBL_MIN)continue;
                double factor=0;for(int k=0;k<3;k++)factor+=(r[k]/n0-rm[k]/nminus)*d[k];factor*=inv4pi/D;
                for(int k=0;k<3;k++)out[0][i+np*(k+3*j)]=c[k]*factor;
            } else {
                nplus=sqrt(dot(rp,rp));cross(d,r,c);D=dot(c,c)+e2*dd;
                if(n0<=cut||nminus<=cut||nplus<=cut||D<=DBL_MIN)continue;
                double sp=0,sm=0,R[3],cxR[3],G[3];
                for(int k=0;k<3;k++){sp+=d[k]*(rp[k]/nplus-r[k]/n0);sm+=d[k]*(r[k]/n0-rm[k]/nminus);}
                for(int k=0;k<3;k++)R[k]=-sp*rp[k]+sm*rm[k];cross(c,R,cxR);
                for(int k=0;k<3;k++)G[k]=(rp[k]/nplus-2*r[k]/n0+rm[k]/nminus)/D-(2*cxR[k]+2*e2*d[k]*(sp-sm))/(D*D);
                double rotation[3][3]={{0,R[2],-R[1]},{-R[2],0,R[0]},{R[1],-R[0],0}};
                for(int axis=0;axis<3;axis++)for(int k=0;k<3;k++)out[axis][i+np*(k+3*j)]=(rotation[axis][k]/D+c[k]*G[axis])*inv4pi;
            }
        }
    }
}
