#include <iostream>
#include <stdio.h>
#include <stack>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "matrix.h"

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
        int simplex_core(Matrix<double> &A, Vector<double> &b, Vector<double> &c, Vector<unsigned int> &bv);
        Matrix<double>& submatrix(const Matrix<double> &A, const Vector<unsigned int> &bv) const;

};

Matrix<double>& ILP::submatrix(const Matrix<double> &A, const Vector<unsigned int> &bv) const{
    unsigned int i,len;


    len = bv.length();

    printf("inside submatrix\r\n");
    for(i=0;i<len;i++){
        printf("getting cols: %d\r\n",bv[i]);
    }
    return *new Matrix<double>(A);
}


int ILP::simplex_core(Matrix<double> &A, Vector<double> &b, Vector<double> &c, Vector<unsigned int> &bv) {
    submatrix(A,bv);
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

    Matrix<double> *Am = new Matrix<double>(aa,M,M+N);
    Vector<double> *bm = new Vector<double>(ca,M);
    Vector<double> *cm = new Vector<double>(b,M+N);
    Vector<unsigned int> *bvm = new Vector<unsigned int>(bv,M);

    simplex_core(*Am,*bm,*cm,*bvm);



    free(aa);
    free(ca);
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

    A.print();
    B.print();
    C.print();



    (p->solve_simplex(A,B,C)).print();


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

