#pragma once
#include "list.hh"
#include "square_matrix.hh"
#include "priority_queue.hh"
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <limits>
#include <iostream>

template<typename T> // forward declarations
class Edge;
template<typename T>
class Graph;

template<typename T>
class Vertex{
  protected:
  typename List<Vertex*>::Iterator position_in_list;
  T element;
  
  public:
  unsigned label;
  bool visited;
  
  Vertex(const T& new_element) : element(new_element), label(0), visited(false) {}
  typename List<Vertex*>::Iterator& get_position_in_list()
  {return position_in_list;}
  T& get_element() {return element;}
  T operator*() {return element;}
  virtual void fun() = 0;
  friend class Edge<T>;
  friend class Graph<T>;
};

template<typename T>
class Edge{
  protected:
  T weight;
  typename List<Edge*>::Iterator position_in_list;
  Vertex<T> *_endVertices[2];
  public:
  Edge(const T& new_weight, Vertex<T> **endVertices) : weight(new_weight){
    _endVertices[0] = endVertices[0]; _endVertices[1] = endVertices[1];}
  typename List<Edge*>::Iterator& get_position_in_list()
  {return position_in_list;}
  T& get_element() {return weight;}
  T operator*() {return weight;}
  Vertex<T>** endVertices() {return _endVertices;}
  Vertex<T>* opposite(Vertex<T> *v){
    if(!v) return v;
    if(*(v->position_in_list) == _endVertices[0]) return _endVertices[1];
    else if(*(v->position_in_list) == _endVertices[1]) return _endVertices[0];
    else throw std::invalid_argument("Edge::opposite : given vertex is not incident on the edge");
  }
  bool isAdjacent(const Edge *e){
    return e->_endVertices[0] == _endVertices[0]
      || e->_endVertices[0] == _endVertices[1]
      || e->_endVertices[1] == _endVertices[0]
      || e->_endVertices[1] == _endVertices[1];
  }
  bool isIncidentOn(Vertex<T> *v){
    return *(v->get_position_in_list()) == _endVertices[0] ||
      *(v->get_position_in_list()) == _endVertices[1];
  }
  virtual void fun() {};
  friend class Vertex<T>;
  friend class Graph<T>;
};

template<typename T>
class Graph{
  protected:
  List<Vertex<T>*> list_of_vertices;
  List<Edge<T>*>   list_of_edges;
  public:
  List<Vertex<T>*>& vertices() {return list_of_vertices;}
  List<Edge<T>*>& edges() {return list_of_edges;}
  virtual Vertex<T>* insertVertex(T new_element) = 0;
  virtual Edge<T>* insertEdge(Vertex<T> *v, Vertex<T> *u, const T& new_element) = 0;
  virtual void eraseVertex(Vertex<T> *v) = 0;
  virtual void eraseEdge(Edge<T> *e) = 0;

  virtual List<Edge<T>*> incidentEdges(Vertex<T> *v) = 0;
  virtual bool isAdjacent(Vertex<T> *v, Vertex<T> *u) = 0;
};

namespace AdjacencyList{
  template<typename T> // forward declaration
  class Graph;
  
  template<typename T>
  class Vertex : public ::Vertex<T>{
    List<::Edge<T>*> incidence_list;
  public:
    Vertex(T& new_element) : ::Vertex<T>(new_element) { }
    List<::Edge<T>*>& incidentEdges() {return incidence_list;}
    bool isAdjacent(::Vertex<T> *v){
      for(typename List<::Edge<T>*>::Iterator iter = incidence_list.begin(); iter != incidence_list.end(); ++iter)
	if(((*iter)->endVertices())[0] == v || ((*iter)->endVertices())[1] == v) return true;
      return false;
    }
    void fun() override {}
    friend Graph<T>;
  };
  
  template<typename T>
  class Edge : public ::Edge<T>{
    typename List<::Edge<T>*>::Iterator position_in_incidence_lists[2];
    public:
    Edge(const T& new_weight, ::Vertex<T> **endVertices) : ::Edge<T>(new_weight, endVertices){ 
      for(unsigned j = 0; j < 2; ++j){
	Vertex<T> *v = dynamic_cast<Vertex<T>*> (endVertices[j]);
	for(typename List<::Edge<T>*>::Iterator i = v->incidentEdges().begin(); i != v->incidentEdges().end(); ++i)
	  if(*i == this)
	    position_in_incidence_lists[j] = i;
      }
    }
    typename List<::Edge<T>*>::Iterator* get_position_in_incidence_lists() {return position_in_incidence_lists;}
    public:
    void fun() override {}
  };
  
  template<typename T>
  class Graph : public ::Graph<T>{
  public:
    ::Vertex<T>* insertVertex(T new_element) override{
      Vertex<T> *new_vertex = new Vertex<T>(new_element);
      Graph<T>::list_of_vertices.insertBack(new_vertex);
      new_vertex->position_in_list = --Graph<T>::list_of_vertices.end();
      return dynamic_cast<::Vertex<T>*> (*new_vertex->position_in_list);
    }
    /*Vertex<T>* insertVertex(T&& new_element){
      T tmp = new_element;
      Vertex<T> *new_vertex = new Vertex<T>(tmp);
      Graph<T>::list_of_vertices.insertBack(new_vertex);
      new_vertex->position_in_list = --Graph<T>::list_of_vertices.end();
      return dynamic_cast<Vertex<T>*> (*new_vertex->position_in_list);
      }*/
    ::Edge<T>* insertEdge(::Vertex<T> *v, ::Vertex<T> *u, const T& new_weight) override{
      ::Vertex<T> *endv[2]={*v->get_position_in_list(), *u->get_position_in_list()};
      Edge<T> *new_edge = new Edge<T>(new_weight, endv);
      Graph<T>::list_of_edges.insertBack(new_edge);
      new_edge->get_position_in_list() = --Graph<T>::list_of_edges.end();
      dynamic_cast<Vertex<T>*>(v)->incidence_list.insertBack(*new_edge->get_position_in_list());
      dynamic_cast<Vertex<T>*>(u)->incidence_list.insertBack(*new_edge->get_position_in_list());
      return dynamic_cast<::Edge<T>*> (*(new_edge->get_position_in_list()));
    }
    void eraseVertex(::Vertex<T> *vertex) override{
      Vertex<T> *v = dynamic_cast<Vertex<T>*> (vertex);
      Graph<T>::list_of_vertices.erase(v->get_position_in_list());
      delete v;
      for(typename List<::Edge<T>*>::Iterator i = v->incidence_list.begin(); i != v->incidence_list.end(); ++i)
	eraseEdge(*i);
    }
    void eraseEdge(::Edge<T> *edge) override{
      Edge<T> *e = dynamic_cast<Edge<T>*> (edge);
      ::Vertex<T> **endVertices = e->endVertices();
      dynamic_cast<Vertex<T>*>(endVertices[0])->incidence_list.erase(e->get_position_in_incidence_lists()[0]);
      dynamic_cast<Vertex<T>*>(endVertices[1])->incidence_list.erase(e->get_position_in_incidence_lists()[1]);
      Graph<T>::list_of_edges.erase(e->get_position_in_list());
      delete e;
    }

    virtual List<::Edge<T>*> incidentEdges(::Vertex<T> *v) override{
      return dynamic_cast<Vertex<T>*> (v)->incidentEdges();
    }
    virtual bool isAdjacent(::Vertex<T> *v, ::Vertex<T> *u) override{
      return dynamic_cast<Vertex<T>*> (v)->isAdjacent(u);
    }
  };
}

namespace AdjacencyMatrix{
  template<typename T> // forward declaration
  class Graph;
  
  template<typename T>
  class Vertex : public ::Vertex<T>{
    int key;
  public:
    Vertex(T& new_element, int new_key) : ::Vertex<T>(new_element),  key(new_key)
    { }
    void fun() override {}
    friend Graph<T>;
  };
  
  template<typename T>
  class Graph : public ::Graph<T>{
    ::Matrix<::Edge<T>*> adjacency_matrix;
    unsigned actual_key = 0;
  public:
    Vertex<T>* insertVertex(T new_element) override{
      Vertex<T> *new_vertex = new Vertex<T>(new_element, actual_key);
      
      Graph<T>::list_of_vertices.insertBack(new_vertex);
      new_vertex->position_in_list = --Graph<T>::list_of_vertices.end();
      
      adjacency_matrix.resize(actual_key+1);
      for(int i = 0; i < adjacency_matrix.dimension(); ++i)
	adjacency_matrix(i, actual_key) = adjacency_matrix(actual_key, i) = nullptr;

      ++actual_key;
      
      return dynamic_cast<Vertex<T>*> (*new_vertex->position_in_list);
    }
    /*Vertex<T>* insertVertex(T&& new_element){
      T tmp = new_element;
      insertVertex(tmp);
      }*/
    
    ::Edge<T>* insertEdge(::Vertex<T> *v, ::Vertex<T> *u, const T& new_weight) override{
      ::Vertex<T> *endv[2]={*v->get_position_in_list(), *u->get_position_in_list()};
      ::Edge<T> *new_edge = new ::Edge<T>(new_weight, endv);
      Graph<T>::list_of_edges.insertBack(new_edge);
      new_edge->get_position_in_list() = --Graph<T>::list_of_edges.end();
      Vertex<T> *v2 = dynamic_cast<Vertex<T>*> (v),  *u2 = dynamic_cast<Vertex<T>*> (u);
      adjacency_matrix(v2->key, u2->key) = adjacency_matrix(u2->key, v2->key) = *(new_edge->get_position_in_list());
      return *(new_edge->get_position_in_list());
    }
    
    void eraseVertex(::Vertex<T> *v) override{
      Graph<T>::list_of_vertices.erase(v->get_position_in_list());
      Vertex<T> *v2 = dynamic_cast<Vertex<T>*> (v);
      
      for(int i = 0; i < adjacency_matrix.dimension(); ++i){
	if(adjacency_matrix(v2->key, i))
	  eraseEdge(adjacency_matrix(v2->key, i));
      }
      
      delete v;
    }
    
    void eraseEdge(::Edge<T> *e) override{
      ::Vertex<T> **endVertices = e->endVertices();
      Vertex<T> *v = dynamic_cast<Vertex<T>*> (endVertices[0]),  *u = dynamic_cast<Vertex<T>*> (endVertices[1]);
      adjacency_matrix(v->key, u->key) = adjacency_matrix(u->key, v->key) = nullptr;
      Graph<T>::list_of_edges.erase(e->get_position_in_list());
      delete e;
    }

    List<::Edge<T>*> incidentEdges(::Vertex<T> *v){
      List<::Edge<T>*> incident_edges;
      Vertex<T> *v2 = dynamic_cast<Vertex<T>*> (v);
      for(int i = 0; i < adjacency_matrix.dimension(); ++i){
	if(adjacency_matrix(v2->key, i))
	  incident_edges.insertBack(adjacency_matrix(v2->key, i));
      }
      
      return incident_edges;
    }
    
    bool isAdjacent(::Vertex<T> *v, ::Vertex<T> *u){
      Vertex<T> *v2 = dynamic_cast<Vertex<T>*> (v),  *u2 = dynamic_cast<Vertex<T>*> (u);
      return adjacency_matrix(v2->key, u2->key) != nullptr;
    }
    void printMatrix(){
      for(int i = 0; i<adjacency_matrix.dimension(); ++i){
	for(int j = 0; j<adjacency_matrix.dimension(); ++j){
	  if(adjacency_matrix(i,j))
	    std::cout << adjacency_matrix(i,j)->get_element() << ' ';
	  else std::cout << 0  << ' ';
	}
	 std::cout << '\n';
      }
    }
  };
}

template<typename T> // unsigned int have to be convertable to T
void generate_random_minimal_connected_graph(Graph<T> *G, unsigned size){
  Vertex<T> *vertex1 = G->insertVertex(0), *vertex2;
  
  srand(time(NULL));
  
  for(unsigned i = 1; i < size; ++i){
    vertex2 = G->insertVertex(i);
    G->insertEdge(vertex1, vertex2, rand() % 300 + 1);
    vertex1 = vertex2;
  }
}

template<typename T> // unsigned int have to be convertable to T
void generate_random_complete_graph(Graph<T> *G, unsigned size){  
  srand(time(NULL));
  
  for(unsigned i = 0; i < size; ++i){
    G->insertVertex(i);
  }

  for(typename List<Vertex<T>*>::Iterator i = G->vertices().begin(); i != G->vertices().end(); ++i){
    typename List<Vertex<T>*>::Iterator j = ++i; --i;
    
    for(; j != G->vertices().end(); ++j)
      G->insertEdge(*i, *j, rand() % 300 + 1);
  }
}

template<typename T>
void generate_random_connected_graph(Graph<T> *G, unsigned size, double density){
  if(density > 1 || density < 0)
    throw std::invalid_argument("generate_random_connected_graph: density must be in [0,1]");
  
  if(density == 1){
    generate_random_complete_graph(G, size);
    return;
  }
  
  generate_random_minimal_connected_graph(G, size);

  srand(time(NULL));

  typename List<Vertex<T>*>::Iterator v1 = G->vertices().begin();
  typename List<Vertex<T>*>::Iterator v2 = G->vertices().begin();
  unsigned l1, l2;
  
  for(unsigned i = size-1; i < size*(size-1)*density/2; ++i){
    do{
      l1 = rand() % (size-1);
      l2 = rand() % (size-1);
      v1 = G->vertices().begin() + l1;
      v2 = G->vertices().begin() + l2;
    }while(G->isAdjacent(*v1, *v2) || (*v1)->get_element() == (*v2)->get_element());
    G->insertEdge(*v1, *v2, rand() % 300 + 1);
  }

}

template<typename T>
void export_to_file(Graph<T> *G, const char *file_name){
  std::ofstream file(file_name, std::ofstream::out);

  file << G->edges().size() << ' ' << G->vertices().size() << '\n';

  for(typename List<Edge<T>*>::Iterator i = G->edges().begin(); i != G->edges().end(); ++i)
    file << ((*i)->endVertices()[0])->get_element() << ' ' << ((*i)->endVertices()[1])->get_element() << ' ' << (*i)->get_element() << '\n';

  file.close();
}

template<typename T>
void import_from_file(Graph<T> *G, const char *file_name){
  std::ifstream file(file_name, std::ifstream::in);
  unsigned number_of_edges, number_of_vertices;
  
  file >> number_of_edges >> number_of_vertices;

  if(file.bad()) throw std::runtime_error("import_to_file: file reading error");

  for(unsigned i = 0; i < number_of_vertices; ++i)
    G->insertVertex(i);

  typename List<Vertex<T>*>::Iterator v1;
  typename List<Vertex<T>*>::Iterator v2;
  unsigned l1, l2, weight;
  for(unsigned i =0; i < number_of_edges; ++i){
    file >> l1 >> l2 >> weight;
    
    if(file.bad()) throw std::runtime_error("import_to_file: file reading error");

    v1 = G->vertices().begin() + l1;
    v2 = G->vertices().begin() + l2;

    G->insertEdge(*v1, *v2, weight);
  }
  
  file.close();
}

template<typename T>
void Dijkstra(Graph<T> *G, Vertex<T> *v, /*Graph<T> *result*/ const char *file_name=nullptr){
  
  for(typename List<Vertex<T>*>::Iterator i = G->vertices().begin(); i != G->vertices().end(); ++i)
    (*i)->label = std::numeric_limits<int>::max();
  
  v->label = 0;

  PriorityQueue<Vertex<T>*> Q;
  
  unsigned *previous_of = new unsigned[G->vertices().size()];
  for(unsigned i=0; i < G->vertices().size(); ++i)
    previous_of[i] = v->get_element();

  for(typename List<Vertex<T>*>::Iterator i = G->vertices().begin(); i != G->vertices().end(); ++i)
    Q.insert(*i, (*i)->label);

  Vertex<T> *u = nullptr;
  while(!Q.empty()){
    Vertex<T> *u = Q.min(); Q.removeMin();
    /* results->insertVertex(u->get_element()); */

    if(!u->visited){
      auto incidentEdges = G->incidentEdges(u);
    
      for(typename List<Edge<T>*>::Iterator e = incidentEdges.begin(); e != incidentEdges.end(); ++e){      
	Vertex<T> *z = (*e)->opposite(u);

	if(z->label < u->label) continue; // out of queue
	if(z->visited) continue; // already done
      
      
	if(u->label + (*e)->get_element() < z->label){ // relaxation
	  z->label = u->label + (*e)->get_element();
	  Q.insert(z, z->label);
	  previous_of[z->get_element()] = u->get_element();
	}
      }
      u->visited = true;
    }
    
  }

  /* auto v0 = G->vertices().begin();
  for(unsigned i = 0; i < G->vertices().size(); ++i){ // saving result to graph
    auto v1 = *(v0+i), v2 = *(v0+previous[i]);
    result->insertEdge(v1, v2, v1->label-v2->label);}*/

  if(file_name == nullptr) return; // below saving to file

  std::ofstream file(file_name, std::ofstream::out);

  file << "Droga (lista kolejnych wierzcho³ków) " << "Koszt\n";

  for(unsigned i=0; i < G->vertices().size(); ++i){
    file << i << ' ';
    unsigned j = i;
    do{
      j = previous_of[j];
      file << j << ' ';
    }while(j != v->get_element());
   
    file << (*(G->vertices().begin()+i))->label << '\n';   
  }
  
}

template<typename T>
void Dijkstra(Graph<T> *G, unsigned n = 0, /*Graph<T> *result*/ const char *file_name=nullptr){
  Vertex<T> *v = *(G->vertices().begin() + n);
  Dijkstra(G, v, file_name);
}

