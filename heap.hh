#pragma once
#include "dynamic_array.hh"

template<typename T>
class Heap{
  public:
  struct Node;

  Heap() {V.push_back(*(new Node));} // V[1] is "empty" (placeholder)
  int size() const {return V.size()-1;}
  Node* left(const Node& node) {return &V[2*V.idx(&node)];}
  Node* right(const Node& node) {return &V[2*V.idx(&node)+1];}
  Node* parent(const Node& node) {return &V[V.idx(&node)/2];}
  bool hasLeft(const Node& node) const {return 2*V.idx(&node) <= size();}
  bool hasRight(const Node& node) const {return 2*V.idx(&node)+1 <= size();}
  bool isRoot(const Node& node) const {return V.idx(&node) == 1;}
  Node* root() {return &V[1];}
  Node* last() {return &V[size()];}

  void upHeapBubling(){
    Node *node = last();
    while(!isRoot(*node)){
      if(node->key < parent(*node)->key)
	swap(*node, *parent(*node));
      else break;
      node = parent(*node);
    }
  }

  void downHeapBubling(){
    Node *parent, *actual= root();
    while(true){
      parent = actual;
      
      if(hasLeft(*actual))
	if(left(*actual)->key < actual->key)
	actual = left(*actual);
      if(hasRight(*parent))
	if(right(*parent)->key < actual->key)
	actual = right(*parent);
    
      if(actual == parent)
	break;
    
      swap(*actual, *parent);
    }
  }
  
  void addLast(const T& new_element, int new_key) {
    V.push_back(Node({new_element, new_key}));
  }
  void removeLast() {V.pop_back();}
  void swap(Node& node1, Node& node2){
    Node tmp = node1; node1 = node2; node2 = tmp;
  }

  private:
  DynamicArray<Node, 10> V;
};

template<typename T>
struct Heap<T>::Node{
  T element;
  int key;
};
