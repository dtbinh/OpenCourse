/* l1tf.c
 *
 * Copyright (C) 2007 Kwangmoo Koh, Seung-Jean Kim and Stephen Boyd.
 *
 * See the file "COPYING.TXT" for copyright and warranty information.
 *
 * Author: Kwangmoo Koh (deneb1@stanford.edu)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "blas.h"
#include "lapack.h"

#ifndef F77_CALL
#define F77_CALL(x) x ## _
#endif

/* macro functions */
#define  max(x,y)       ((x)>(y)?(x):(y))
#define  min(x,y)       ((x)<(y)?(x):(y))

#define  TRUE           (1)
#define  FALSE          (0)

/* constant variables for fortran call */
const double done       = 1.0;
const double dtwo       = 2.0;
const double dminusone  =-1.0;
const int    ione       = 1;
const int    itwo       = 2;
const int    ithree     = 3;
const int    iseven     = 7;

/* function declaration */
void   Dx(const int n, const double *x, double *y); /* y = D*x */
void   DTx(const int n, const double *x, double *y); /* y = D'*x */
void   yainvx(int n, const double a, const double *x, double *y); /* y = a./x */
void   vecset(int n, const double val, double *dst);
double vecsum(int n, const double *src);

void   print_dvec(int n, const double *x); /* for debug */
void   print_ivec(int n, const int *x); /* for debug */

/****************************************************************************
 *                                                                          *
 *              l1tf : main routine for l1 trend filtering                  *
 *                                                                          *
 ****************************************************************************/
int l1tf(const int n, const double *y, const double lambda, double *x)
{
    /* parameters */
    const double ALPHA      = 0.01; /* linesearch parameter (0,0.5] */
    const double BETA       = 0.5;  /* linesearch parameter (0,1) */
    const double MU         = 2;    /* IPM parameter: t update */
    const double MAXITER    = 40;   /* IPM parameter: max iter. of IPM */
    const double MAXLSITER  = 20;   /* IPM parameter: max iter. of linesearch */
    const double TOL        = 1e-4; /* IPM parameter: tolerance */

    /* dimension */
    const int    m          = n-2;  /* length of Dx */

    /* matrix of size 3xm */
    double *S, *DDTF;

    /* vector of size n */
    double *DTz;

    /* vector of size m */
    double *z;                      /* dual variable */
    double *mu1, *mu2;              /* dual of dual variables */
    double *f1, *f2;                /* constraints */
    double *dz, *dmu1, *dmu2;       /* search directions */
    double *w, *rz, *tmp_m2, *tmp_m1, *Dy, *DDTz;
    double norm2_res, norm2_newres;

    double pobj1, pobj2, pobj, dobj, gap;
    double t;                       /* barrier update parameter */
    double step, ratio;
    double *dptr;

    int iters, lsiters;             /* IPM and linesearch iterations */
    int i, info;
    int *IPIV;
    int ddtf_chol;

    /* memory allocation */
    S       = malloc(sizeof(double)*m*3);
    DDTF    = malloc(sizeof(double)*m*7);

    DTz     = malloc(sizeof(double)*n);
    Dy      = malloc(sizeof(double)*m);
    DDTz    = malloc(sizeof(double)*m);

    z       = malloc(sizeof(double)*m);
    mu1     = malloc(sizeof(double)*m);
    mu2     = malloc(sizeof(double)*m);
    f1      = malloc(sizeof(double)*m);
    f2      = malloc(sizeof(double)*m);
    dz      = malloc(sizeof(double)*m);
    dmu1    = malloc(sizeof(double)*m);
    dmu2    = malloc(sizeof(double)*m);
    w       = malloc(sizeof(double)*m);
    rz      = malloc(sizeof(double)*m);
    tmp_m1  = malloc(sizeof(double)*m);
    tmp_m2  = malloc(sizeof(double)*m);

    IPIV    = malloc(sizeof(int)*m);

    /* INITIALIZATION */

    /* DDT : packed representation of symmetric banded matrix
     *
     * fortran style storage (should be transposed in C)
     *  6  6  6  ...  6  6  6
     * -4 -4 -4  ... -4 -4  *
     *  1  1  1  ...  1  *  *
     */

    ddtf_chol = TRUE;
    /* DDTF stores Cholesky factor of DDT in packed symm.-band representation */
    dptr = DDTF;
    for (i = 0; i < m; i++)
    {
        *dptr++ = 6.0;
        *dptr++ =-4.0;
        *dptr++ = 1.0;
    }
    F77_CALL(dpbtrf)("L",&m,&itwo,DDTF,&ithree,&info);
    if (info > 0) /* if Cholesky factorization fails, try LU factorization */
    {
        fprintf(stderr,"Changing to LU factorization\n");
        ddtf_chol = FALSE;
        dptr = DDTF;
        for (i = 0; i < m; i++)
        {
            dptr++;
            dptr++;
            *dptr++ = 1.0;
            *dptr++ =-4.0;
            *dptr++ = 6.0;
            *dptr++ =-4.0;
            *dptr++ = 1.0;
        }
        F77_CALL(dgbtrf)(&m,&m,&itwo,&itwo,DDTF,&iseven,IPIV,&info);
    }

    Dx(n,y,Dy);

    /* variable initialization */
    t       = -1;
    step    =  1;

    vecset(m,0.0,z);
    vecset(m,1.0,mu1);
    vecset(m,1.0,mu2);
    vecset(m,-lambda,f1);
    vecset(m,-lambda,f2);

    DTx(m,z,DTz); /* DTz = D'*z */
    Dx(n,DTz,DDTz); /* DDTz = D*D'*z */

    fprintf(stderr,"%s %13s %12s %8s\n","Iteration","Primal obj.", \
            "Dual obj.","Gap");

    /*---------------------------------------------------------------------*
     *                          MAIN LOOP                                  *
     *---------------------------------------------------------------------*/
    for (iters = 0; iters <= MAXITER; iters++)
    {
        double zTDDTz;

        /* COMPUTE DUALITY GAP */
        zTDDTz = F77_CALL(ddot)(&n,DTz,&ione,DTz,&ione);

        F77_CALL(dcopy)(&m,Dy,&ione,w,&ione); /* w = D*y-(mu1-mu2)  */
        F77_CALL(daxpy)(&m,&dminusone,mu1,&ione,w,&ione);
        F77_CALL(daxpy)(&m,&done,mu2,&ione,w,&ione);

        F77_CALL(dcopy)(&m,Dy,&ione,tmp_m2,&ione); /* tmp_m2 = D*y-D*D'*z */
        F77_CALL(daxpy)(&m,&dminusone,DDTz,&ione,tmp_m2,&ione);

        F77_CALL(dcopy)(&m,w,&ione,tmp_m1,&ione); /* tmp_m1 = w*(D*D')^-1*w */
        if (ddtf_chol == TRUE) {
            F77_CALL(dpbtrs)("L",&m,&itwo,&ione,DDTF,&ithree,tmp_m1,&m,&info);
        } else {
            F77_CALL(dgbtrs)("N",&m,&itwo,&itwo,&ione,DDTF,&iseven,IPIV,tmp_m1,&m,&info);
        }

        pobj1 = 0.5*F77_CALL(ddot)(&m,w,&ione,tmp_m1,&ione)
                +lambda*(vecsum(m,mu1)+vecsum(m,mu2));
        pobj2 = 0.5*zTDDTz + lambda*F77_CALL(dasum)(&m,tmp_m2,&ione);
        pobj  = min(pobj1,pobj2);
        dobj  =-0.5*zTDDTz+F77_CALL(ddot)(&m,Dy,&ione,z,&ione);

        gap   = pobj - dobj;

        fprintf(stderr,"%6d %15.4e %13.5e %10.2e\n",iters,pobj,dobj,gap);

        /* STOPPING CRITERION */

        if (gap <= TOL)
        {
            fprintf(stderr,"Solved\n");
            F77_CALL(dcopy)(&n,y,&ione,x,&ione);
            F77_CALL(daxpy)(&n,&dminusone,DTz,&ione,x,&ione);
            return(0);
        }

        if (step >= 0.2)
        {
            t = max(2*m*MU/gap, 1.2*t);
        }

        /* CALCULATE NEWTON STEP */

        F77_CALL(dcopy)(&m,DDTz,&ione,rz,&ione); /* rz = D*D'*z-w */
        F77_CALL(daxpy)(&m,&dminusone,w,&ione,rz,&ione);

        yainvx(m,+1.0/t,f1,dz); /* dz = r = D*y-D*D'*z+(1/t)./f1-(1/t)./f2 */
        yainvx(m,-1.0/t,f2,tmp_m1);
        F77_CALL(daxpy)(&m,&done,tmp_m1,&ione,dz,&ione);
        F77_CALL(daxpy)(&m,&done,tmp_m2,&ione,dz,&ione);

        dptr = S; /* S = D*D'-diag(mu1./f1-mu2./f2) */
        for (i = 0; i < m; i++)
        {
            *dptr++ = 6-mu1[i]/f1[i]-mu2[i]/f2[i];
            *dptr++ =-4.0;
            *dptr++ = 1.0;
        }

        F77_CALL(dpbsv)("L",&m,&itwo,&ione,S,&ithree,dz,&m,&info); /* dz=S\r */

        norm2_res = F77_CALL(ddot)(&m,rz,&ione,rz,&ione);
        for (i = 0; i < m; i++)
        {
            double tmp1, tmp2;
            tmp1 = -mu1[i]*f1[i]-(1/t);
            tmp2 = -mu2[i]*f2[i]-(1/t);
            norm2_res += tmp1*tmp1+tmp2*tmp2;

            dmu1[i] = -(mu1[i]+((1/t)+dz[i]*mu1[i])/f1[i]);
            dmu2[i] = -(mu2[i]+((1/t)-dz[i]*mu2[i])/f2[i]);
        }
        norm2_res = sqrt(norm2_res);
        
        /* BACKTRACKING LINESEARCH */

        ratio = 2;   /* any number larger than 1/0.99 */
        for (i = 0; i < m; i++)
        {
            if (dmu1[i]<0 && -mu1[i]/dmu1[i]<ratio) ratio = -mu1[i]/dmu1[i];
            if (dmu2[i]<0 && -mu2[i]/dmu2[i]<ratio) ratio = -mu2[i]/dmu2[i];
        }
        step = min(1,0.99*ratio);

        /* compute new values of z, dmu1, dmu2, f1, f2 */
        F77_CALL(daxpy)(&m,&step,dz  ,&ione,z  ,&ione);
        F77_CALL(daxpy)(&m,&step,dmu1,&ione,mu1,&ione);
        F77_CALL(daxpy)(&m,&step,dmu2,&ione,mu2,&ione);

        for (lsiters = 0; lsiters < MAXLSITER; lsiters++)
        {
            int linesearch_skip;
            double diff_step;

            linesearch_skip = 0;
            for (i = 0; i < m; i++)
            {
                f1[i] =  z[i]-lambda;
                f2[i] = -z[i]-lambda;
                if (f1[i] > 0 || f2[i] > 0) linesearch_skip = 1;
            }
            if (linesearch_skip != 1)
            {
                DTx(m,z,DTz); /* rz = D*D'*z-D*y-(mu1-mu2) */
                Dx(n,DTz,DDTz);
                F77_CALL(dcopy)(&m,DDTz,&ione,rz,&ione);
                F77_CALL(daxpy)(&m,&dminusone,Dy,&ione,rz,&ione);
                F77_CALL(daxpy)(&m,&done,mu1,&ione,rz,&ione);
                F77_CALL(daxpy)(&m,&dminusone,mu2,&ione,rz,&ione);

                /* UPDATE RESIDUAL */

                /* compute  norm([rz; -mu1.*f1-1/t; -mu2.*f2-1/t]) */
                norm2_newres = F77_CALL(ddot)(&m,rz,&ione,rz,&ione);
                for (i = 0; i < m; i++)
                {
                    double tmp1, tmp2;
                    tmp1 = -mu1[i]*f1[i]-(1/t);
                    tmp2 = -mu2[i]*f2[i]-(1/t);
                    norm2_newres += tmp1*tmp1+tmp2*tmp2;
                }
                norm2_newres = sqrt(norm2_newres);

                if (norm2_newres <= (1-ALPHA*step)*norm2_res)
                    break;
            }
            diff_step = -step*(1.0-BETA);
            F77_CALL(daxpy)(&m,&diff_step,dz  ,&ione,z  ,&ione);
            F77_CALL(daxpy)(&m,&diff_step,dmu1,&ione,mu1,&ione);
            F77_CALL(daxpy)(&m,&diff_step,dmu2,&ione,mu2,&ione);
            step *= BETA;
        }
    }
    fprintf(stderr,"Maxiter exceeded\n");
    F77_CALL(dcopy)(&n,y,&ione,x,&ione);
    F77_CALL(daxpy)(&n,&dminusone,DTz,&ione,x,&ione);
    return(0);
}

double l1tf_lambdamax(const int n, double *y)
{
    int i, m, info;
    double maxval;
    double *vec, *mat, *dptr;
    int *piv;

    m = n-2;
    vec = malloc(sizeof(double)*m);
    mat = malloc(sizeof(double)*7*m);
    piv = malloc(sizeof(int)*m);

    Dx(n,y,vec);
    dptr = mat;
    for (i = 0; i < m; i++)
    {
        *dptr++ = 6;
        *dptr++ =-4.0;
        *dptr++ = 1.0;
    }

    F77_CALL(dpbsv)("L",&m,&itwo,&ione,mat,&ithree,vec,&m,&info);
    if (info > 0) /* if Cholesky factorization fails, try LU factorization */
    {
        fprintf(stderr,"Changing to LU factorization\n");
        dptr = mat;
        for (i = 0; i < m; i++)
        {
            dptr++;
            dptr++;
            *dptr++ = 1.0;
            *dptr++ =-4.0;
            *dptr++ = 6.0;
            *dptr++ =-4.0;
            *dptr++ = 1.0;
        }

        F77_CALL(dgbsv)(&m,&itwo,&itwo,&ione,mat,&iseven,piv,vec,&m,&info);
        if (info > 0) return -1.0;  /* if LU fails, return -1 */
    }
    maxval = 0;
    for (i = 0; i < m; i++)
    {
        if (fabs(vec[i]) > maxval) maxval = fabs(vec[i]);
    }

    return maxval;
}

/* Computes y = D*x, where x has length n
 *
 *     | 1 -2  1  0  0 |
 * y = | 0  1 -2  1  0 |*x
 *     | 0  0  1 -2  1 |
 */
void Dx(const int n, const double *x, double *y)
{
    int i;
    for (i = 0; i < n-2; i++,x++)
        *y++ = *x-*(x+1)-*(x+1)+*(x+2); /* y[0..n-3]*/
}

/* Computes y = D^T*x, where x has length n
 *
 *     | 1  0  0 |
 *     |-2  1  0 |
 * y = | 1 -2  1 |*x
 *     | 0  1 -2 |
 *     | 0  0  1 |
 */
void DTx(const int n, const double *x, double *y)
{
    int i;
    *y++ = *x;                          /* y[0]     */
    *y++ = -*x-*x+*(x+1);               /* y[1]     */
    for (i = 2; i < n; i++,x++)
        *y++ = *x-*(x+1)-*(x+1)+*(x+2); /* y[2..n-1]*/
    *y++ = *x-*(x+1)-*(x+1); x++;       /* y[n]     */
    *y = *x;                            /* y[n+1]   */
} 

/* Computes y = a./x, where x has length n */
void yainvx(int n, const double a, const double *x, double *y)
{
    while (n-- != 0)
        *y++ = a/ *x++;
}

/* Set dst = val, where dst has length n */
void vecset(int n, const double val, double *dst)
{
    while (n-- != 0)
        *dst++ = val;
}

/* Computes sum(x) */
double vecsum(int n, const double *x)
{
    double ret = 0.0;
    while (n-- != 0)
        ret += *x++;
    return ret;
}


/* Print vector of double type */
void print_dvec(int n, const double *x)
{
    int i;
    fprintf(stdout,"\n");
    for (i = 0; i < n; i++)
        fprintf(stdout,"%e\n",x[i]);
    fprintf(stdout,"\n");
}

/* Print vector of int type */
void print_ivec(int n, const int *x)
{
    int i;
    fprintf(stdout,"\n");
    for (i = 0; i < n; i++)
        fprintf(stdout,"%d  ",x[i]);
    fprintf(stdout,"\n");
}
