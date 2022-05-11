#pragma once

template<typename T>
class List{
  struct Node;
  
  Node *header, *tailer;
  int _size;

  Node* createNode(const T& element) const;
  
  public:
  class Iterator;

  List(){
    header = new Node; tailer = new Node;
    header->previous = nullptr; header->next = tailer;
    tailer->previous = header;  tailer->next = nullptr;
    _size = 0;
  }
  /*~List() {
    for(Iterator iter = begin(); iter != end(); ++iter)
      delete iter.node->previous;
    delete tailer;
    }*/

  int size() const {return _size;}
  bool empty() const {return _size == 0;}
  Iterator begin() {return Iterator(header->next);}
  Iterator end() {return Iterator(tailer);}
  void insertFront(const T& new_element);
  void insertBack(const T& new_element);
  void insertBefore(const Iterator& next, const T& new_element);
  void eraseFront();
  void eraseBack();
  void erase(Iterator& to_erase);
};

template<typename T>
struct List<T>::Node{
  T element;
  Node *next, *previous;
  bool operator==(const Node& node) {
    return next == node.next && previous == node.previous;
  }
};

template<typename T>
class List<T>::Iterator{
  Node *node;

  Iterator(Node* new_node) : node(new_node) {}

  public:
  Iterator(const Iterator& iter) : node(iter.node) {}
  Iterator() : node(nullptr) {}
  T& operator*() {return node->element;}
  bool operator==(const Iterator& iterator) {return node == iterator.node;}
  bool operator!=(const Iterator& iterator) {return !(*this == iterator);}
  
  Iterator& operator++(){
    node = node->next;
    return *this;
  }
  Iterator& operator--(){
    node = node->previous;
    return *this;
  }
  Iterator operator+(unsigned i){
    if(!i) return *this;
    
    Iterator iter = *this;
    do ++iter; while(--i);
    return iter;
  }
  friend List;
};

template<typename T>
typename List<T>::Node* List<T>::createNode(const T &element) const{
  Node *node = new Node;
  node->element= element;

  return node;
}

template<typename T>
void  List<T>::insertFront(const T& element){
  insertBefore(begin(), element);
}

template<typename T>
void  List<T>::insertBack(const T& element){
  insertBefore(end(), element);
}

template<typename T>
void  List<T>::insertBefore(const Iterator &next, const T& element){
  Node *new_node = createNode(element);
  
  new_node->previous = next.node->previous; new_node->next = next.node;
  next.node->previous->next = new_node;
  next.node->previous = new_node;
  ++_size;
}

template<typename T>
void  List<T>::eraseFront(){
  Iterator iter(begin());
  erase(iter);
}

template<typename T>
void  List<T>::eraseBack(){
  Iterator iter(--end());
  erase(iter);
}

template<typename T>
void  List<T>::erase(Iterator &to_erase){
  Node *node_to_erase = to_erase.node;
  node_to_erase->previous->next = node_to_erase->next;
  node_to_erase->next->previous = node_to_erase->previous;
  delete node_to_erase;
  --_size;
}
