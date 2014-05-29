#include <iostream>
#include <stdio.h>
#include <stack>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include "lapacke.h"
#include "cblas.h"
#include "matrix.h"
#include <assert.h>

template <class T, class Q>
class BB {
    public:
        Q C;
        T container;
        int branch(T ps){
            unsigned int pop;
            this->C.push(ps);
            pop = this->C.top();
            this->C.pop();
            printf("Got: %d\r\n",pop);
            return(0);
        }
};


class PartialSolution {
    public:
        virtual int eval() = 0;
        virtual int branch() = 0;
};

class ILP : public PartialSolution{
    public:
        int eval(){ return 0;}
        int branch(){ return 0;}

        /** \brief Minimize c'x st. Ax=b. Using Simplex Algorithm.
         * \fn Vector& solve_simplex(Matrix &A, Vector &b, Vector &c);
         * \param[in] A constraint matrix.
         * \param[in] b constraint vector.
         * \param[in] c cost vector.
         * \return Minimizing vector.
         *
         */

        Vector<double>& solve_simplex(Matrix<double> &A, Vector<double> &b, Vector<double> &c);
        int solve_simplex(const double *A, double *b, double *c, unsigned int M, unsigned int N);
    private:
        int simplex_core(const double *A, const double *b, double *c, unsigned int *bv, int M, int N);
        int price(double *r, unsigned int *dv, const double *A, const double *c, const unsigned int *bv, unsigned int M,unsigned int N);
        void print_matrix(const double *A,unsigned int M, unsigned int N);
        int findBasisToEnter(const double *r,unsigned int N);
        int findBasisToLeave(const double *A,const double *b, const unsigned int j, unsigned int M);
        int getBasicSolution(const double *A,const unsigned int *bv, double *b,unsigned int M);
        int copybycols(double *dest, const double *src,const unsigned int *col, unsigned int M, unsigned int numcols);

};


int ILP::getBasicSolution(const double *A,const unsigned int *bv,double *b,unsigned int M){
    double *B,*btemp;
    int i,*pivot;

    B = (double *)malloc(M*M*sizeof(double)); // basis
    btemp = (double *)malloc(M*sizeof(double)); // basis
    pivot = (int *)malloc(M*sizeof(int));

    memcpy(btemp,b,M*sizeof(double));

    copybycols(B,A,bv,M,M);
    // compute LU factorization of B.
    LAPACKE_dgetrf(LAPACK_COL_MAJOR,M,M,B,M,pivot);

    // solve this system. cbv is now lambda.
    LAPACKE_dgetrs(LAPACK_COL_MAJOR,'N',M,1,B,M,pivot,btemp,M);

    for(i=0;i<M;i++){
        b[i] = btemp[i];
    }

    

    free(B);
    free(pivot);
}


int ILP::findBasisToEnter(const double *r,unsigned int N){
    int i;
    for(i=0;i<N;i++){
        if(r[i]<0)
            return(i);
    }
    return(-1);
}

int ILP::findBasisToLeave(const double *A,const double *b, const unsigned int j,unsigned int M){
    double ratio;
    double min_ratio;
    int i,i_out = -1;
    min_ratio = INFINITY;
    for(i=0;i<M;i++){
        if((ratio=b[i]/A[(j*M)+i])>0){
            if(ratio<min_ratio){
                min_ratio=ratio;
                i_out = i;
            }
        }
    }
    return(i_out);
}

int ILP::price(double *r, unsigned int *dv, const double *A, const double *c, const unsigned int *bv, unsigned int M,unsigned int N){
    double *lambda, *B,*cbv,*r_temp;
    int *pivot;
    bool flag;
    unsigned int i,j,k,*dv_temp;
    lambda = (double *)malloc(M*sizeof(double));
    pivot = (int *)malloc(M*sizeof(int));
    B = (double *)malloc(M*M*sizeof(double)); // basis
    cbv = (double *)malloc(M*sizeof(double));
    r_temp = (double *)malloc((N-M)*sizeof(double));
    dv_temp = (unsigned int *)malloc((N-M)*sizeof(unsigned int));

    copybycols(B,A,bv,M,M);
    copybycols(cbv,c,bv,1,M);

    // compute LU factorization of B.
    LAPACKE_dgetrf(LAPACK_COL_MAJOR,M,M,B,M,pivot);

    // solve this system. cbv is now lambda.
    LAPACKE_dgetrs(LAPACK_COL_MAJOR,'T',M,1,B,M,pivot,cbv,M);
    for(i=0;i<M;i++){
        lambda[i] = cbv[i];
    }


    k=0;
    for(i=0;i<N;i++){
        flag = 0;

        // is i in the basis?
        for(j=0;j<M;j++){
            if(i==bv[j]){
                flag = true;
                break;
            }
        }

        // if not compute price. 
        if(!flag){
          r[k]=c[i]-cblas_ddot(M,lambda,1,&A[i*M],1);
          dv[k] = i;
          k++;
        }

        
        
    }


    //loop through non-basis vectors and compute price. 

   free(pivot);
   free(lambda);
   free(cbv);
   free(B);

}


void ILP::print_matrix(const double *A,unsigned int M, unsigned int N){
    unsigned int i,j;
    for(i=0;i<M;i++){
        for(j=0;j<N;j++){
            printf("%f\t",A[(j*M)+i]);
        }
        printf("\r\n");
    }
}

int ILP::copybycols(double *dest, const double *src,const unsigned int *col, unsigned int M, unsigned int numcols){
    unsigned int i;
    for(i=0;i<numcols;i++){
        memcpy(&dest[i*M],&src[M*col[i]],M*sizeof(double));
    }
}

int ILP::simplex_core(const double *A, const double *b, double *c, unsigned int *bv, int M, int N) {
    int i,j,k;
    double *r,*bcpy;
    unsigned int *dv;
    int dvidx;

    if(M>N){
        assert(M>N);
        printf("Number of rows must be larger than number of columns\r\n");
        return(1);
    }
    

    r = (double *)malloc((N-M)*sizeof(double));
    dv = (unsigned int *)malloc((N-M)*sizeof(unsigned int));
    bcpy = (double *)malloc(M*sizeof(double));

    memcpy(bcpy,b,M*sizeof(double));
    
    price(r, dv, A, c, bv, M, N);

    printf("r: \r\n");
    print_matrix(r,N-M,1);

    while(0<=(dvidx=findBasisToEnter(r,N-M))){
        j=dv[dvidx];
        if((k=findBasisToLeave(A,bcpy,j,M))<0){
            printf("Problem is unbounded from below.");
            exit(1);
        }

        printf("Basis to enter: %d Basis to leave: %d \r\n",j,bv[k]);

        dv[dvidx]=bv[k];
        bv[k]=j;

        // b = linsolve(A(:,bv),b)

        getBasicSolution(A,bv,bcpy,M);

        price(r, dv, A, c, bv, M, N);
        
        printf("r: \r\n");
        print_matrix(r,N-M,1);

    }

    
    memset(c,0,N*sizeof(double));
    for(i=0;i<M;i++){
        c[bv[i]]=bcpy[i];
    }

    free(r);
    free(dv);

    return(0);
}

int ILP::solve_simplex(const double *A, double *b, double *c, unsigned int M, unsigned int N){
    unsigned int i,j,k,*bv;
    double *aa,*ca;
    aa = (double *) calloc(M*(M+N),sizeof(double));
    ca = (double *) calloc((M+N),sizeof(double));
    bv = (unsigned int *) calloc(M,sizeof(unsigned int));
    memcpy(aa,A,M*N*sizeof(double));


    /* Lets now augment aa with the identity matrix, though these should be chased out of the basis eventually */
    for(i=0;i<M;i++){
        aa[(M*N)+(i*M)+i]=1;
    }

    k=0;
    for(i=N;i<N+M;i++){
        ca[i]=1;
        bv[k++]=i;
    }

    simplex_core(aa,b,ca,bv,M,M+N);

    printf("Solution to part I: c\r\n");

    print_matrix(ca,N+M,1);

    printf("Solution to part I: bv\r\n");

    for(i=0;i<M;i++){
        printf("%d ",bv[i]);
    }


    simplex_core(A,b,c,bv,M,N);

    printf("\r\nSolution to part II: \r\n");

    print_matrix(c,N,1);

    printf("Solution to part II: bv\r\n");

    for(i=0;i<M;i++){
        printf("%d ",bv[i]);
    }



    free(aa);
    free(ca);
    free(bv);
    return 0;
}


Vector<double>& ILP::solve_simplex(Matrix<double> &A, Vector<double> &b, Vector<double> &c){
    solve_simplex(A.getData(),b.getData(),c.getData(),A.getNumRows(),A.getNumCols());
    return b;
}

using namespace std;

int main()
{
    ILP *p = new ILP;
    double a[9] = {2,1,2,3,3,1};
    double b[4] = {4,3};
    double c[3] = {4,1,1};


    Matrix<double> A(a,2,3);
    Vector<double> B(b,2);
    Vector<double> C(c,3);

    p->solve_simplex(A,B,C);


    delete p;
}

