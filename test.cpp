#include "graph.hh"
#include <iostream>

int main(){
  const char *file_name = "graph", *f1 = "matrix", *f2 = "list";

  
  Graph<int> *G_matrix = new AdjacencyMatrix::Graph<int>;
  Graph<int> *G_list = new AdjacencyList::Graph<int>;

	
  generate_random_connected_graph(G_list, 7, 0.5); // genereating graph
  export_to_file(G_list, file_name);                      // copying graph
  import_from_file(G_matrix, file_name);

  Dijkstra(G_matrix, 0, f1);

  Dijkstra(G_list, 0, f2);

  
  return 0;
}
