#pragma once
#include <stdexcept>

template<typename T, unsigned int initial_capacity=1>
class DynamicArray{
  T *_array;
  int _size, _capacity;

  void realloc(unsigned new_capacity= 0){
    if(new_capacity < 2*_capacity)
      new_capacity = 2*_capacity;
    
    T *new_array = new T[new_capacity];
      
      for(unsigned i=0; i<_size; ++i)
	new_array[i] = _array[i];

      delete[] _array;
      _array = new_array;
      _capacity = new_capacity;
  }

  public:
  DynamicArray(){
    _array = new T[initial_capacity];
    _capacity = initial_capacity;
    _size = 0;
  }

  ~DynamicArray() {delete[] _array;}

  /*const T& operator[] (unsigned int index) const{
    if(index >= _size)
      throw std::out_of_range("DynamicArray: index out of range");

    return _array[index];
    }*/

  T& operator[] (unsigned int index){
    if(index >= _capacity)
      throw std::out_of_range("DynamicArray: index out of range");
    if(index >= _size)
      _size = index+1;
      
    return _array[index];
  }

  void push_back(const T& new_element) {
    if(_size == _capacity) realloc();

    ++_size;
    _array[_size-1] = new_element;
  }

  void pop_back() {
    --_size;
  }

  int size() const {return _size;}
  int idx(const T *ptr) const {return ptr - _array;}
  int capacity() const {return _capacity;}
  void resize(unsigned n){
    if(n > _capacity) realloc(n);
    _size = n;
  }
};
