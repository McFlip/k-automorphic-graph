#include <iostream>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;
using namespace boost;

int main(int argc, char* argv[]) {

  typedef adjacency_list<vecS, vecS, undirectedS> graph_type;
  graph_type graph1(0);
  graph_type graph2(0);
  ifstream inputFile;
  inputFile.open(argv[1]);
  dynamic_properties dp(ignore_other_properties);
  read_graphviz(inputFile, graph1, dp);
//   read_graphviz(inputFile,graph2, dp);
  graph2 = graph1;
  inputFile.close();

  // Create callback to print mappings
  vf2_print_callback<graph_type, graph_type> callback(graph1, graph2);
  vf2_subgraph_iso(graph1, graph2, callback);

}
