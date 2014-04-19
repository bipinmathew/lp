#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stack>

class Matrix {
    public:

        /** \brief Constructor to create Matrix object from C-style double array
         *
         * \param[in] in pointer to C-style matrix.
         * \param[in] n Number of rows.
         * \param[in] m Number of columns.
         * \return none
         *
         */
        Matrix(const double *in, unsigned int n, unsigned int m);

        /** \fn Matrix& pivot(unsigned int pi,unsigned int pj);
        * \brief Pivot matrix on a given element.
        * \param[in] pi Row index of element to pivot on.
        * \param[in] pj Column index of element to pivot on.
        * \return 0 on success, non-zero on failure.
        *
        */

        /** \brief Destructor frees allocated memory and cleans up.
         * \return none
         *
         */

        ~Matrix();
        Matrix& pivot(unsigned int pi,unsigned int pj);

        /** \fn print()
        *  \brief Output matrix to standard out.
        *  \return void
        */
        void print();


    private:
        double *A;
        unsigned int N,  M;
};

Matrix::Matrix(const double *in, unsigned int n, unsigned int m){
    N = n;
    M = m;
    A = (double *)calloc(n*m,sizeof(double));
    memcpy(A,in,M*N*sizeof(double));
}

Matrix::~Matrix(){
    free(A);
}

Matrix& Matrix::pivot(unsigned int pi,unsigned int pj){
    unsigned int i,j;
    double p,multiplier;
    p = 1/A[(N*pi)+pj];
    for(j=0;j<N;j++){
        A[(N*pi)+j] *= p;
    }

    for(i=0;i<M;i++){
        if(i==pi)
            continue;
        multiplier = (A[(N*i)+pj]);
        for(j=0;j<N;j++){
            //printf("(%d,%d)\r\n",i,j);
            A[(N*i)+j] = A[(N*i)+j]-(multiplier*A[(N*pi)+j]);
        }
    }
    return *this;
}



void Matrix::print(){
    unsigned int i, j;
    for(j=0;j<M;j++){
        for(i=0;i<N;i++){
            printf("%f\t",A[(N*j)+i]);
        }
        printf("\r\n");
    }
    printf("\r\n");
}


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
};

using namespace std;

int main()
{
    double a[9] = {1,3,1,1,1,-1,3,11,5};
    Matrix A(a,3,3);
    A.print();
    A.pivot(0,0).print();
    A.print();
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
