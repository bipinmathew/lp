#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

template <class T>
Matrix<T>::Matrix(const T *in, unsigned int m,unsigned int n){
    unsigned int i,j,k;
    M = m;
    N = n;
    A = (T *)calloc(n*m,sizeof(T));
    //memcpy(A,in,M*N*sizeof(double)); // This is for row major order, we are looking for column major order.
    k=0;
    for(j=0;j<N;j++){
        for(i = 0;i<M;i++){
            A[k++]=in[(N*i)+j];
        }
    }
}

template <class T>
Matrix<T>::Matrix( const Matrix& other){

    N = other.N;
    M = other.M;
    A = (T *)malloc(N*M*sizeof(T));
    memcpy(A,other.A,N*M*sizeof(T));

}

template <class T>
Matrix<T>::~Matrix(){
    free(A);
}

template <class T>
Matrix<T>& Matrix<T>::pivot(unsigned int pi,unsigned int pj){
    unsigned int i,j;
    double p,multiplier;
    p = 1/A[(M*pj)+pi];
    for(j=0;j<N;j++){
        A[(M*j)+pi] *= p;
    }

    for(i=0;i<M;i++){
        if(i==pi)
            continue;
        multiplier = (A[(M*pj)+i]);
        for(j=0;j<N;j++){
            A[(M*j)+i] = A[(M*j)+i]-(multiplier*A[(M*j)+pi]);
        }
    }
    return *this;
}


template <class T>
Matrix<T>& Matrix<T>::appendColumn(const Matrix &in){
    A=(T *)realloc(A,M*(N+1)*sizeof(T));
    memcpy(&A[N*M],in.A,M*sizeof(T));
    N+=1;
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::appendRow(const Matrix &in){
    unsigned int i,numRows,numCols;
    unsigned int startFrom,newSize;


    numRows = getNumRows();
    numCols = getNumCols();
    newSize= numRows*(numCols+1);

    A=(T *)realloc(A,newSize*sizeof(T));

    startFrom = numRows;


    for(i=0;i<numCols;i++){
        memmove(&A[startFrom+1],&A[startFrom],(newSize-startFrom)*sizeof(T));
        A[startFrom++] = in.A[i];
        startFrom+=numRows;
    }
    setNumRows(numRows+1);
    return *this;
}


template < >
void Matrix<double>::print(){
    unsigned int i, j;
    for(j=0;j<M;j++){
        for(i=0;i<N;i++){
            printf("%f\t",A[(M*i)+j]);
        }
        printf("\r\n");
    }
    printf("\r\n");
}


template <>
void Matrix< unsigned int >::print(){
    unsigned int i, j;
    for(j=0;j<M;j++){
        for(i=0;i<N;i++){
            printf("%d\t",A[(M*i)+j]);
        }
        printf("\r\n");
    }
    printf("\r\n");
}

template class Matrix<double>;
template class Matrix<unsigned int>;
