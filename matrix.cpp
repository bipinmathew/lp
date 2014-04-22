#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

Matrix::Matrix(const double *in, unsigned int m,unsigned int n){
    M = m;
    N = n;
    A = (double *)calloc(n*m,sizeof(double));
    memcpy(A,in,M*N*sizeof(double));
}


Matrix::Matrix( const Matrix& other){
    N = other.N;
    M = other.M;
    A = (double *)malloc(N*M*sizeof(double));
    memcpy(A,other.A,N*M*sizeof(double));
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
            A[(N*i)+j] = A[(N*i)+j]-(multiplier*A[(N*pi)+j]);
        }
    }
    return *this;
}

Matrix& Matrix::appendRow(const Matrix &in){
    A=(double *)realloc(A,N*(M+1)*sizeof(double));
    memcpy(&A[N*M],in.A,N*sizeof(double));
    M+=1;
    return *this;
}

Matrix& Matrix::appendColumn(const Matrix &in){
    unsigned int i,numRows,numCols;
    unsigned int startFrom,newSize;


    numRows = getNumRows();
    numCols = getNumCols();
    newSize= numRows*(numCols+1);

    A=(double *)realloc(A,newSize*sizeof(double));

    startFrom = numCols;


    for(i=0;i<numRows;i++){
        memmove(&A[startFrom+1],&A[startFrom],(newSize-startFrom)*sizeof(double));
        A[startFrom++] = in.A[i];
        startFrom+=numCols;
    }
    setNumCols(numCols+1);
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

