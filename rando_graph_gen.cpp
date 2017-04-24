#include <iostream>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/random.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
using namespace std;
using namespace boost;


int main(int argc, char* argv[]) {

 typedef adjacency_list<setS, vecS, undirectedS> graph_type;
 size_t numNodes = atoi(argv[1]);
 size_t numEdges = atoi(argv[2]);
 ofstream outputFile;
 outputFile.open(argv[3]);
 boost::mt19937 rng;
 rng.seed(uint32_t(time(0)));

// Build graph
 graph_type graph(0);

 generate_random_graph(graph, numNodes, numEdges, rng, false, false);

 write_graphviz(outputFile, graph);
 outputFile.close();

 return 0;
}