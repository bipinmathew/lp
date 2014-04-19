#include <iostream>
#include <stdio.h>
#include <stack>


int pivot(double *A,double *b,unsigned int N, unsigned int M,unsigned int pi,unsigned int pj){
    unsigned int i,j;
    double p,multiplier;
    p = 1/A[(N*pi)+pj];
    for(j=0;j<N;j++){
        A[(N*pi)+j] *= p;
    }
    b[pi] *= p;


    for(i=0;i<M;i++){
        if(i==pi)
            continue;
        multiplier = (A[(N*i)+pj]);
        for(j=0;j<N;j++){
            //printf("(%d,%d)\r\n",i,j);
            A[(N*i)+j] = A[(N*i)+j]-(multiplier*A[(N*pi)+j]);
        }
        b[i]=b[i]-(multiplier*b[pi]);
    }
    return 0;
}

void printm(double *A,unsigned int N, unsigned int M){
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
    double A[9] = {1,3,1,1,1,-1,3,11,5};
    double b[3] = {9,1,35};
    ILP PS;
    BB< ILP,std::stack<ILP> > ILPSolver;
    printm(A,3,3);
    pivot(A,b,3,3,0,0);

    printm(A,3,3);
    printm(b,3,1);
    return 0;
}
