#include "graph.hh"
#include <iostream>
#include <iomanip>
#include <chrono>

int main(){
  srand(time(nullptr));

  int number_of_graphs = 100;
  int sizes_of_graphs[] = {10, 50, 100, 500, 1000};
  double densities_of_graphs[] = {0.25, 0.5, 0.75, 1};
  
  auto begin = std::chrono::high_resolution_clock::now(), end  = std::chrono::high_resolution_clock::now();
  
  std::ofstream matrix_f("dat/Adjacency_matrix.csv"), list_f("dat/Adjacency_list.csv");
  matrix_f << "size density[%] time[ns]\n";
  list_f << "size density[%] time[ns]\n";  

  const char *file_name = "graph.dat";
  
  for(int s = 0; s < sizeof(sizes_of_graphs)/sizeof(*sizes_of_graphs); ++s)
    for(int d = 0; d < sizeof(densities_of_graphs)/sizeof(*densities_of_graphs); ++d)
      for(int n = 0; n < number_of_graphs; ++n){
	Graph<int> *G_matrix = new AdjacencyMatrix::Graph<int>;
	Graph<int> *G_list = new AdjacencyList::Graph<int>;
	int size = sizes_of_graphs[s];
	double density = densities_of_graphs[d];
	
	generate_random_connected_graph(G_list, size, density); // genereating graph
	export_to_file(G_list, file_name);                      // copying graph
	import_from_file(G_matrix, file_name);

	begin = std::chrono::high_resolution_clock::now();
	Dijkstra(G_matrix);
	end = std::chrono::high_resolution_clock::now();
	auto time = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
	
	matrix_f << size << ' ' << density*100 << ' ' << time.count() << '\n';

	
	begin = std::chrono::high_resolution_clock::now();
	Dijkstra(G_list);
	end = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
	
	list_f << size << ' ' << density*100 << ' ' << time.count() << '\n';

	
	delete G_matrix; delete G_list;
      }
  
  return 0;
}
