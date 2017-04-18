#include <iostream>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_utility.hpp>
#include "metis.h"

using namespace std;
using namespace boost;

const int K = 2;
idx_t ncon = 1;
idx_t nparts = K;
int main(int argc, char* argv[]) {
  typedef subgraph< adjacency_list<vecS, vecS, undirectedS, no_property,
    property< edge_index_t, int > > > graph_type;
  typedef graph_traits<graph_type>::vertex_iterator vertex_iter;
  typedef std::pair<vertex_iter, vertex_iter> vrange_t;
  typedef graph_traits<graph_type>::adjacency_iterator adj_iter;
  typedef std::pair<adj_iter, adj_iter> adjrange_t;
  typedef property_map<graph_type, vertex_index_t>::type IndexMap;

  graph_type graph1(0);
  ifstream inputFile;
  inputFile.open(argv[1]);
  dynamic_properties dp(ignore_other_properties);
  read_graphviz(inputFile, graph1, dp);

  inputFile.close();

  // Create callback to print mappings
  //  vf2_print_callback<graph_type, graph_type> callback(graph1, graph2);
  //  vf2_subgraph_iso(graph1, graph2, callback);

  // Partition the graph
  idx_t nvert = num_vertices(graph1);
  int nedge = num_edges(graph1);
  idx_t *xadj = new idx_t[nvert + 1];
  idx_t *adjncy = new idx_t[nedge * 2];
  idx_t objval = 0;
  idx_t *part = new idx_t[nvert];
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);


  vrange_t vpair = vertices(graph1);
  int i=0, j=0;
  xadj[0] = 0;

  // Transform the adjacency list into the Compressed Storage Format
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

  // Print out to verify
  cout << "METIS partitioning: " << endl;
  for(i = 0; i < nvert + 1; ++i){
    cout << xadj[i] << ' ';
  }
  cout << endl;
  for(i = 0; i < nedge * 2; ++i){
    cout << adjncy[i] << ' ';
  }
  cout << endl;


  // Use k-way graph partition algorithm
  METIS_PartGraphKway(&nvert, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts, NULL, NULL, NULL, &objval, part);

  // Print out the partition vector
  for(i=0; i < nvert; ++i){
    cout << part[i] << ' ';
  }
  cout << endl;

  // Create induced subgraphs based on partition vector
  graph_type *subgraph_vect = new graph_type[K];
  IndexMap index = get(vertex_index, graph1);
  for (i=0; i<K; ++i) {
    subgraph_vect[i] = graph1.create_subgraph();
  }
  graph_type::children_iterator ci, ci_end;
  boost::tie(ci, ci_end) = graph1.children();
  j = 0;
  for (i=0; i < nvert; ++i){
    cout << "Mapping " << index[i] << " to subgraph " << part[i] << endl;
    add_vertex(index[i], subgraph_vect[part[i]]);
    if (part[i] > j){
      ++j;
      ++ci;
    }
    add_vertex(index[i], *ci);
  }

  // print for testing
  cout << "root:" << endl;
  print_graph(graph1, get(vertex_index, graph1));
  cout << endl;
  for (i=0; i < K; ++i){
    cout << "subgraph " << i << ":" << endl;
    print_graph(subgraph_vect[i], get(vertex_index, subgraph_vect[i]));
  }
  cout << endl;
  int num = 1;
  for (boost::tie(ci, ci_end) = graph1.children(); ci != ci_end; ++ci){
    cout << "G" << num++ << ":" << endl;
    print_graph(*ci, get(vertex_index, *ci));
    cout << endl;
  }

  //cleanup
  delete[] xadj;
  delete[] adjncy;
  delete[] part;
  delete[] subgraph_vect;
  xadj = NULL;
  adjncy = NULL;
  part = NULL;
  subgraph_vect = NULL;


  return 0;
}
