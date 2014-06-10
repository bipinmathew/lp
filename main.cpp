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
#include <unordered_map>
#include <assert.h>

using namespace std;

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
        void print_matrix(const double *A,unsigned int M, unsigned int N);

        /** \brief Minimize c'x st. Ax=b. Using Simplex Algorithm.
         * \fn Vector& solve_simplex(Matrix &A, Vector &b, Vector &c);
         * \param[in] A constraint matrix.
         * \param[in] b constraint vector.
         * \param[in] c cost vector.
         * \return Minimizing vector.
         *
         */

        Vector<double>& solve_simplex(Matrix<double> &A, Vector<double> &b, Vector<double> &c);
        int solve_simplex(const double *A, const double *b, double *c, unsigned int M, unsigned int N);
    private:
        int simplex_core(const double *A, const double *b, double *c, unsigned int *bv, int M, int N);
        int price(double *r, unsigned int *dv, const double *A, const double *c, const unsigned int *bv, unsigned int M,unsigned int N);
        int findBasisToEnter(const double *r,unsigned int N);
        int findBasisToLeave(const double *A,const double *b, const unsigned int j, unsigned int M);
        int getBasicSolution(const double *A,const unsigned int *bv, double *b,unsigned int M);
        int copybycols(double *dest, const double *src,const unsigned int *col, unsigned int M, unsigned int numcols);

};


int ILP::getBasicSolution(const double *A,const unsigned int *bv,double *b,unsigned int M){
    double *B;
    int i,*pivot;

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

    while(0<=(dvidx=findBasisToEnter(r,N-M))){
        j=dv[dvidx];
        if((k=findBasisToLeave(A,bcpy,j,M))<0){
            printf("Problem is unbounded from below.");
            exit(1);
        }

        dv[dvidx]=bv[k];
        bv[k]=j;

        // b = linsolve(A(:,bv),b)

        getBasicSolution(A,bv,bcpy,M);

        price(r, dv, A, c, bv, M, N);
        

    }

    
    memset(c,0,N*sizeof(double));
    for(i=0;i<M;i++){
        c[bv[i]]=bcpy[i];
    }

    free(r);
    free(dv);
    free(bcpy);

    return(0);
}

int ILP::solve_simplex(const double *A, const double *b, double *c, unsigned int M, unsigned int N){
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

    // Solve the phase I problem. to obtail and BFS
    simplex_core(aa,b,ca,bv,M,M+N);

    // Take that BFS and solve the original problem.
    simplex_core(A,b,c,bv,M,N);

    print_matrix(c,N,1);





    free(aa);
    free(ca);
    free(bv);
    return 0;
}


Vector<double>& ILP::solve_simplex(Matrix<double> &A, Vector<double> &b, Vector<double> &c){
    solve_simplex(A.getData(),b.getData(),c.getData(),A.getNumRows(),A.getNumCols());
    return b;
}


#define NAMECARD 0
#define ROWCARD 1
#define COLCARD 2
#define RHSCARD 3
#define BOUNDCARD 4


class Constraint {
    public:
        int setName(char *n);
        int setType(const char t);
        int setNum(const unsigned int num);
        string& getName();
        char getType();
        unsigned int getNum();
        Constraint(char *n, const char t,const unsigned int num);
    private:
        string name;
        char type;
        unsigned int num;
};

Constraint::Constraint(char *n, const char t,const unsigned int num){
    setName(n);
    setType(t);
    setNum(num);
}

int Constraint::setName(char *n){
    name=n;
}

int Constraint::setType(const char t){
    type = t;
}

int Constraint::setNum(const unsigned int n){
    num=n;
}

string& Constraint::getName(){
    return name;
}

char Constraint::getType(){
    return type;
}

unsigned int Constraint::getNum(){
    return num;
}

class Variable {
};


unsigned int stris(const char *s1,const char *s2){
    return(strncmp(s1,s2,strlen(s2))==0);
}

int main()
{
    FILE *fp;
    ILP *p = new ILP;
    //double a[9] = {2,1,2,3,3,1};
    //double b[4] = {4,3};
    //double c[3] = {4,1,1};
    char str[1024];
    char chardata[1024],name[1024],type;
    char varname[1024], constraint1[1024], constraint2[1024];
    double f1,f2;
    unsigned int card,numVars=0,numConstraints=0,numvars,k;
    double *A,*b,*c;
    A = NULL;
    unordered_map<string,Constraint * > ConstraintSet;
    unordered_map<string,Variable * > VariableSet;

    /* Matrix<double> A(a,2,3);
    Vector<double> B(b,2);
    Vector<double> C(c,3);

    p->solve_simplex(A,B,C); */

    fp = fopen("./data/test.mps","r");

    while(!feof(fp)){
        fgets(str,1024,fp);
        if(stris(str,"ENDATA")){
            break;
        }
        else if(str[0]=='*'){
            continue;
        }
        else if(str[0]!=' '){
            if(stris(str,"NAME")){
                card=NAMECARD;
                sscanf(str,"%s %s",chardata,name);
                printf("GOT: %s , %s\r\n",chardata,name);
            }
            else if(stris(str,"ROWS")){
                card=ROWCARD;
            }
            else if(stris(str,"COLUMNS")){
                card=COLCARD;
            }
            else if(stris(str,"RHS")){
                card=RHSCARD;
            }
            else if(stris(str,"BOUNDS")){
                card=BOUNDCARD;
            }
        }
        else if(str[0]==' '){
            switch(card){
                case NAMECARD:
                    printf("  manage name card\r\n");
                break;
                case ROWCARD:
                    if(sscanf(str," %c %s",&type,chardata) != 2){
                        printf("Parse eror in ROWCARD");
                        exit(1);
                    }
                    if(type!='N'){
                        ConstraintSet[chardata] = new Constraint(chardata,type,numConstraints++);
                    }
                    if(type=='N'){
                        printf("Got N Type!\r\n");
                    }
                break;
                case COLCARD:
                    numvars=sscanf(str," %s %s %lf %s %lf",varname, constraint1, &f1, constraint2, &f2);

                    printf("Got numvars: %d\r\n",numvars);


                    /* Need to add logic to dynamically enlarge the A matrix as necessary. */
                    if(3<=numvars){
                        if(ConstraintSet.end() != ConstraintSet.find(constraint1)){
                            if(VariableSet.end() == VariableSet.find(varname)){
                                numVars++;
                                VariableSet[varname] = new Variable();
                                A = (double *)realloc(A,numVars*numConstraints*sizeof(double));
                            }
                            A[((numVars-1)*numConstraints)+ConstraintSet[constraint1]->getNum()] = f1;
                        }
                    }
                    if(5==numvars){
                        if(ConstraintSet.end() != ConstraintSet.find(constraint2)){
                            if(VariableSet.end() == VariableSet.find(varname)){
                                numVars++;
                                VariableSet[varname] = new Variable();
                                A = (double *)realloc(A,numVars*numConstraints*sizeof(double));
                            }
                            A[((numVars-1)*numConstraints)+ConstraintSet[constraint2]->getNum()] = f2;
                        }
                    }
                    if((numvars!=3) && (numvars!=5)){
                        printf("Parse error in Column Card.");
                        exit(1);
                    }

                break;
                case RHSCARD:
                    printf("  manage rhs card\r\n");
                break;
                case BOUNDCARD:
                    printf("  manage bound card\r\n");
                break;
                default:
                    printf("  unknown card\r\n");
                break;
            }
        }
        else{
            printf("Parse error\r\n");
            return(0);
        }
        // printf("%s",str);
    }

   printf("Printing matrix of %d rows and %d cols...\r\n",numConstraints,numVars); 
    p->print_matrix(A,numConstraints,numVars);
    // for(auto& x: ConstraintSet){
    //     cout << "key: " << x.first << " name: " << x.second->getName() << " type: " << x.second->getType() << " num: " << x.second->getNum() << endl;

    // }

    fclose(fp);
    delete p;
}

