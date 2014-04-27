#include <iostream>
#include <stdio.h>
#include <stack>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
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
        Matrix& solve_simplex(Matrix &A, Vector &b, Vector &c);
    private:
        int find_bfs(Matrix &A, Matrix &b, Matrix &c);
        Matrix& build_tableau(const Matrix &A, const Vector &b, const Vector &c);
};

Matrix& ILP::build_tableau(const Matrix &A, const Vector &b, const Vector &c){
    Matrix *Tbl = new Matrix(A);
    return Tbl->appendColumn(b);
}

Matrix& ILP::solve_simplex(Matrix &A, Vector &b, Vector &c){
    return build_tableau(A,b,c);
}

using namespace std;

int main()
{
    ILP *p = new ILP;
    double a[9] = {2,1,2,3,3,1};
    double b[4] = {4,3};
    double c[3] = {4,1,1};

    Matrix A(a,2,3);
    Vector B(b,2);
    Vector C(c,3);

    A.appendRow(C).print();
    B.print();
    C.print();



    // (p->solve_simplex(A,B,C)).print();


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
