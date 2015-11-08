BLAS_LIBDIR = /usr/local/OpenBLAS/lib
LAPACK_LIBDIR = /usr/local/lapack/lapack-3.5.0
LAPACKE_INCLUDEDIR = /usr/local/lapack/lapack-3.5.0/lapacke/include

Debug: 
	g++ -pedantic -Wall -O0 -g -o simplex -std=c++0x -I$(LAPACKE_INCLUDEDIR) -L$(LAPACK_LIBDIR) -L$(BLAS_LIBDIR) main.cpp matrix.cpp -llapacke -llapack -lopenblas  
