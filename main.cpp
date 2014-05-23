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
        int solve_simplex(const double *A, double *b, const double *c, unsigned int M, unsigned int N);
    private:
        int simplex_core(const double *A, double *b, const double *c, unsigned int *bv, int M, int N);
        int price(double *r, unsigned int *dv, const double *A, const double *c, const unsigned int *bv, unsigned int M,unsigned int N);
        void print_matrix(const double *A,unsigned int M, unsigned int N);
        int findBasisToEnter(const double *r,unsigned int N);
        int findBasisToLeave(const double *A,const double *b, const unsigned int j, unsigned int M);
        int getBasicSolution(const double *A,const unsigned int *bv, double *b,unsigned int M);
        int copybycols(double *dest, const double *src,const unsigned int *col, unsigned int M, unsigned int numcols);

};


int ILP::getBasicSolution(const double *A,const unsigned int *bv,double *b,unsigned int M){
    double *B;
    int *pivot;

    B = (double *)malloc(M*M*sizeof(double)); // basis
    pivot = (int *)malloc(M*sizeof(int));

    copybycols(B,A,bv,M,M);
    // compute LU factorization of B.
    LAPACKE_dgetrf(LAPACK_COL_MAJOR,M,M,B,M,pivot);

    // solve this system. cbv is now lambda.
    LAPACKE_dgetrs(LAPACK_COL_MAJOR,'N',M,1,B,M,pivot,b,M);
    

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
    double *lambda, *B,*cbv;
    int *pivot;
    bool flag;
    unsigned int i,j,k;
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

    printf("bv: \r\n");
    for(i=0;i<M;i++)
        printf("%d ",bv[i]);

    printf("\r\n");

    printf("A: \r\n");
    print_matrix(A,M,N);

    printf("N: %d \r\n",N);

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
          printf("not in basis: %d \r\n",i);
          r[k]=c[i]-cblas_ddot(M,cbv,1,&A[i*M],1);
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

int ILP::simplex_core(const double *A, double *b, const double *c, unsigned int *bv, int M, int N) {
    int i,j,k;
    double *r;
    unsigned int *dv;
    int dvidx;

    if(M>N){
        assert(M>N);
        printf("Number of rows must be larger than number of columns\r\n");
        return(1);
    }
    

    r = (double *)malloc((N-M)*sizeof(double));
    dv = (unsigned int *)malloc((N-M)*sizeof(unsigned int));
    

    price(r, dv, A, c, bv, M, N);
    while(0<=(dvidx=findBasisToEnter(r,N-M))){
        printf("dvidx: %d\r\n",dvidx);
        j=dv[dvidx];
        printf("Basis to leave: %d \r\n",j);
        if((k=findBasisToLeave(A,b,j,M))<0){
            printf("Problem is unbounded from below.");
            exit(1);
        }

        dv[dvidx]=bv[k];
        bv[k]=j;

        // b = linsolve(A(:,bv),b)

        getBasicSolution(A,bv,b,M);

        price(r, dv, A, c, bv, M, N);

    }

    printf("Printing matrix\r\n\r\n");

    print_matrix(b,M,1);

    printf("Basis vectors: \r\n\r\n");

    for(i=0;i<M;i++){
        printf("%d ",bv[i]);
    }

    free(r);
    free(dv);

    return(0);
}

int ILP::solve_simplex(const double *A, double *b, const double *c, unsigned int M, unsigned int N){
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

    delete p;
}

