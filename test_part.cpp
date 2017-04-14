#include <iostream>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/graphviz.hpp>
#include "metis.h"

using namespace std;
using namespace boost;

int main(int argc, char* argv[]) {

  typedef adjacency_list<vecS, vecS, undirectedS> graph_type;
  graph_type graph1(0);
  //  graph_type graph2(0);
  ifstream inputFile;
  inputFile.open(argv[1]);
  dynamic_properties dp(ignore_other_properties);
  read_graphviz(inputFile, graph1, dp);

  //  graph2 = graph1;
  inputFile.close();

  // Create callback to print mappings
  //  vf2_print_callback<graph_type, graph_type> callback(graph1, graph2);
  //  vf2_subgraph_iso(graph1, graph2, callback);

  int nvert = num_vertices(graph1);
  int nedge = num_edges(graph1);
  idx_t *xadj = new idx_t[nvert];
  idx_t *adjncy = new idx_t[nedge];

  xadj[0] = 0;
  typedef graph_traits<graph_type>::vertex_iterator vertex_iter;
  typedef std::pair<vertex_iter, vertex_iter> vrange_t;
  typedef graph_traits<graph_type>::adjacency_iterator adj_iter;
  typedef std::pair<adj_iter, adj_iter> adjrange_t;

  vrange_t vpair = vertices(graph1);
  int i=0, j=0;
  for(vpair.first; vpair.first != vpair.second; ++vpair.first){
    xadj[i+1] = xadj[i] + out_degree(*vpair.first, graph1);
    adjrange_t adjpair = adjacent_vertices(*vpair.first, graph1);
    j = xadj[i];
    for(adjpair.first; adjpair.first != adjpair.second; ++adjpair.first){
      adjncy[j] = *adjpair.first;
      ++j;
    }
    ++i;
  }
  /*  for(i = 0; i < nvert + 1; ++i){
    cout << xadj[i] << ' ';
  }
  cout << endl;
  for(i = 0; i < nedge * 2; ++i){
    cout << adjncy[i] << ' ';
  }
  cout << endl;
  */
  delete[] xadj;
  delete[] adjncy;
  return 0;
}
