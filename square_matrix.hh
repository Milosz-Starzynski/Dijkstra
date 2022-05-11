#pragma once
#include "dynamic_array.hh"

template<typename T>
class Matrix{
  DynamicArray<DynamicArray<T>*> matrix;
  unsigned _dimension = 0;
  public:
  void resize(unsigned n){
    matrix.resize(n);
    
    for(unsigned i=0; i < n-1; ++i)
      matrix[i]->resize(n);
    matrix[n-1] = new DynamicArray<T>;
    matrix[n-1]->resize(n);
    
    _dimension = n;
  }

  DynamicArray<T>* operator[] (unsigned index){
    return matrix[index];
  }
  
  T& operator() (unsigned i, unsigned j){
    return (*matrix[i])[j];
  }

  unsigned dimension() {return _dimension;}
};
