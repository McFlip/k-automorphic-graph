#include <iostream>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/graphviz.hpp>
using namespace std;
using namespace boost;


int main(int argc, char* argv[]) {

 typedef adjacency_list<setS, vecS, undirectedS> graph_type;
 int fromNode;
 int toNode;
 ifstream inputFile;
 inputFile.open(argv[1]);
 char trashBuffer[256];
 const int MAX_NODES = 325729;


 // Build graph1
 int num_vertices1 = 8; graph_type graph1(num_vertices1);
 add_edge(0, 6, graph1); add_edge(0, 7, graph1);
 add_edge(1, 5, graph1); add_edge(1, 7, graph1);
 add_edge(2, 4, graph1); add_edge(2, 5, graph1); add_edge(2, 6, graph1);
 add_edge(3, 4, graph1);

 // Build graph2
 graph_type graph2(MAX_NODES);

 inputFile.getline(trashBuffer, 256);
 inputFile.getline(trashBuffer, 256);
 inputFile.getline(trashBuffer, 256);
 inputFile.getline(trashBuffer, 256);

 while(inputFile >> fromNode >> toNode)
 {
  add_edge(fromNode,  toNode, graph2);
 }

 // Create callback to print mappings
 vf2_print_callback<graph_type, graph_type> callback(graph1, graph2);

 // Create ordering of small graph
 vf2_order_by_mult<graph_type> vos(graph1);

 // Print out all subgraph isomorphism mappings between graph1 and graph2.
 // Vertices and edges are assumed to be always equivalent.
  vf2_subgraph_iso(graph1, graph2, callback, vos, vertices_equivalent([&graph_small, &graph_large](Graph::vertex_descriptor small_vd, Graph::vertex_descriptor large_vd) {
    return out_degree(small_vd, graph_small) == out_degree(large_vd, graph_large);})
   );

  // print out the graph
//   write_graphviz(cout, graph2);

 return 0;
}