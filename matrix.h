#ifndef MATRIX_H
#define MATRIX_H
class Matrix {
    friend class Vector;
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

        /** \brief Copy constructor to do a deep copy of Matrix object.
         * \fn Matrix( const Matrix& other)
         * \param[in] other Matrix to copy.
         * \return none
         *
         */
        Matrix( const Matrix& other);


        /** \brief Destructor frees allocated memory and cleans up.
         * \return none
         *
         */
        ~Matrix();



        /** \fn Matrix& pivot(unsigned int pi,unsigned int pj);
        * \brief Pivot matrix on a given element.
        * \param[in] pi Row index of element to pivot on.
        * \param[in] pj Column index of element to pivot on.
        * \return 0 on success, non-zero on failure.
        *
        */
        Matrix& pivot(unsigned int pi,unsigned int pj);


        /** \brief Append row to matrix.
         * \fn Matrix& appendRow(const Matrix& in)
         * \param[in] Matrix to append.
         * \return reference to updated Matrix object
         *
         */
        Matrix& appendRow(const Matrix& in);


        /** \brief Append column to matrix.
         * \fn Matrix& appendColumn(const Matrix& in)
         * \param[in] Matrix to append.
         * \return reference to updated Matrix object
         *
         */
        Matrix& appendColumn(const Matrix& in);

        unsigned int getNumRows() const {
            return M;
        }

        unsigned int getNumCols() const {
            return N;
        }


        /** \fn print()
        *  \brief Output matrix to standard out.
        *  \return void
        */
        void print();


    private:
        double *A;
        unsigned int N,  M;/**< M rows by N columns */
        void setNumRows (unsigned int numrows){
            M = numrows;
        }

        void setNumCols(unsigned int numcols){
            N = numcols;
        }
};

class Vector : public Matrix{
    public:
        /** \brief Constructor to create Vector from C-style double array
         * \fn Vector(const double *in, unsigned int n)
         * \param[in] in Pointer to array.
         * \param[in] n Number of elements in array.
         * \return none
         *
         */

        Vector(const double *in, unsigned int n) : Matrix(in,n,1){};
        Vector& transpose(){
            unsigned int temp;
            temp = getNumRows();
            setNumRows(getNumCols());
            setNumCols(temp);
            return *this;
        }
};

#endif // MATRIX_H
