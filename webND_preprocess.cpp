#include <iostream>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
using namespace std;
using namespace boost;


int main(int argc, char* argv[]) {

 typedef adjacency_list<setS, vecS, undirectedS> graph_type;
 int fromNode;
 int toNode;
 ifstream inputFile;
 inputFile.open(argv[1]);
 ofstream outputFile;
 outputFile.open(argv[2]);
 char trashBuffer[256];
 const int MAX_NODES = 325729;

// Build graph
 graph_type graph(MAX_NODES);

 inputFile.getline(trashBuffer, 256);
 inputFile.getline(trashBuffer, 256);
 inputFile.getline(trashBuffer, 256);
 inputFile.getline(trashBuffer, 256);

 while(inputFile >> fromNode >> toNode)
 {
  add_edge(fromNode,  toNode, graph);
 }

 write_graphviz(outputFile, graph);
 inputFile.close();
 outputFile.close();

 return 0;
}