#include <iostream>
#include <stdio.h>
#include <stack>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
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
        int solve_simplex(const double *A, const double *b, const double *c, unsigned int M, unsigned int N);
    private:
        int simplex_core(const double *A, const double *b, const double *c, unsigned int *bv, int M, int N);
        int price(double *r, const double *A, const double *c, unsigned int *bv, unsigned int M,unsigned int N);
        void print_matrix(const double *A,unsigned int M, unsigned int N);
        int in(double val,const double *A, unsigned int M, unsigned int N);
        int copybycols(double *dest, const double *src,unsigned int *col, unsigned int M, unsigned int numcols);

};


int ILP::price(double *r, const double *A, const double *c, unsigned int *bv, unsigned int M,unsigned int N){
    double *lambda, *B,*cbv;
    int *pivot;
    unsigned int i,k;
    lambda = (double *)malloc(M*sizeof(double));
    pivot = (int *)malloc(M*sizeof(int));
    B = (double *)malloc(M*M*sizeof(double)); // basis
    cbv = (double *)malloc(M*sizeof(double));

    copybycols(B,A,bv,M,M);
    copybycols(cbv,c,bv,1,M);

    // compute LU factorization of B.
    LAPACKE_dgetrf(LAPACK_COL_MAJOR,M,M,B,M,pivot);

    // solve this system. cbv is now lambda.
    LAPACKE_dgetrs(LAPACK_COL_MAJOR,'T',M,1,B,M,pivot,cbv,M);

    printf("cbv: \r\n");
    print_matrix(cbv,M,1);

    printf("A: \r\n");
    print_matrix(A,M,N);

    k=0;
    for(i=0;i<N;i++){
        if(!in((double)i,(double *)bv,M,1)){
            printf("not in basis: %d \r\n",i);
            r[k++]=c[i]-cblas_ddot(M,cbv,1,&A[i*M],1);
        }
    }

    //loop through non-basis vectors and compute price. 



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

int ILP::copybycols(double *dest, const double *src,unsigned int *col, unsigned int M, unsigned int numcols){
    unsigned int i;
    for(i=0;i<numcols;i++){
        memcpy(&dest[i*M],&src[M*col[i]],M*sizeof(double));
    }
}



int ILP::in(double val,const double *A, unsigned int M, unsigned int N){
    unsigned int i;
    for(i=0;i<M*N;i++){
        if(A[i]==val)
            return(1);
    }
    return(0);
}


int ILP::simplex_core(const double *A, const double *b, const double *c, unsigned int *bv, int M, int N) {
    int i,j,numRHS=1,*pivot;
    unsigned int *dv;
    double *Bbv,*D,*cbv,*cdv,*lambda;

    if(M>N){
        assert(M>N);
        printf("Number of rows must be larger than number of columns\r\n");
        return(1);
    }
    

    Bbv = (double *)malloc(M*M*sizeof(double)); // basis
    D = (double *)malloc(M*(N-M)*sizeof(double)); // non-basis
    cbv = (double *)malloc(M*sizeof(double));
    cdv = (double *)malloc((N-M)*sizeof(double));
    dv = (unsigned int *)malloc((N-M)*sizeof(unsigned int));
    
    lambda = (double *)malloc(M*sizeof(double));
    pivot = (int *)malloc(M*sizeof(int));

    price(cdv, A, c, bv, M, N);
    // copy vector in basis to its own matrix.

    // copybycols(Bbv,A,bv,M,M);
    // copybycols(cbv,c,bv,1,M);
    // copybycols(D,A,dv,M,N-M);
    // copybycols(cdv,A,dv,1,N-M);

    //LAPACKE_dgesv(LAPACK_COL_MAJOR,M,numRHS,Bbv,M,pivot,cbv,M);
    // cbv now has the value lambda.

    printf("Printing matrix\r\n\r\n");

    //print_matrix(Bbv,M,M);
    print_matrix(cdv,M,1);


    free(Bbv);
    free(cbv);
    return(0);
}

int ILP::solve_simplex(const double *A, const double *b, const double *c, unsigned int M, unsigned int N){
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

    free(aa);
    free(ca);
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

    // A.print();
    // B.print();
    // C.print();





    p->solve_simplex(A,B,C);


    //Matrix A(a,3,3);
    //A.print();
    //A.pivot(0,0).pivot(1,1).appendRow(c).appendColumn(b).print();
    //A.print();
    /* double b[3] = {9,1,35};
    ILP PS;
    BB< ILP,std::stack<ILP> > ILPSolver;
    Matrix M;
    M.printm(A,3,3);
    M.pivot(A,b,3,3,0,0);

    M.printm(A,3,3);
    M.printm(b,3,1);
    return 0;*/
}

