//try to implement Matrix class, based on 
// http://www.parashift.com/c++-faq-lite/freestore-mgmt.html#faq-16.16

#ifndef MATRIX_H
#define MATRIX_H
  template<typename T>  // See section on templates for more
 class Matrix {
 public:
   Matrix(unsigned nrows, unsigned ncols);
   // Throws a BadSize object if either size is zero
   class BadSize { };
 
   // Based on the Law Of The Big Three:
  ~Matrix();
   Matrix(const Matrix<T>& m);
   Matrix<T>& operator= (const Matrix<T>& m);
 
   // Access methods to get the (i,j) element:
   T&       operator() (unsigned i, unsigned j);        // subscript operators often come in pairs
   const T& operator() (unsigned i, unsigned j) const;  // subscript operators often come in pairs
   // These throw a BoundsViolation object if i or j is too big
   class BoundsViolation { };
 
 private:
   unsigned nrows_, ncols_;
   T* data_;
 };
 
 template<typename T>
 inline T& Matrix<T>::operator() (unsigned row, unsigned col)
 {
   if (row >= nrows_ || col >= ncols_) throw BoundsViolation();
   return data_[row*ncols_ + col];
 }
 
 template<typename T>
 inline const T& Matrix<T>::operator() (unsigned row, unsigned col) const
 {
   if (row >= nrows_ || col >= ncols_)
     throw BoundsViolation();
   return data_[row*ncols_ + col];
 }
 
 template<typename T>
 inline Matrix<T>::Matrix(unsigned nrows, unsigned ncols)
   : nrows_ (nrows)
   , ncols_ (ncols)
 //, data_ <--initialized below (after the 'if/throw' statement)
 {
   if (nrows == 0 || ncols == 0)
     throw BadSize();
   data_ = new T[nrows * ncols];
 }
 
 template<typename T>
 inline Matrix<T>::~Matrix()
 {
   delete[] data_;
 } 


#endif


