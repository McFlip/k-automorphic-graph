#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_utility.hpp>
#include "metis.h"

using namespace std;
using namespace boost;

const int K = 2;
idx_t ncon = 1;
idx_t nparts = K;
int main(int argc, char* argv[])
{	
  typedef subgraph< adjacency_list<vecS, vecS, undirectedS, no_property,
    property< edge_index_t, int > > > graph_type;
  typedef graph_traits<graph_type>::vertex_iterator vertex_iter;
  typedef graph_traits<graph_type>::edge_iterator edge_iter;
  typedef std::pair<vertex_iter, vertex_iter> vrange_t;
  typedef graph_traits<graph_type>::adjacency_iterator adj_iter;
  typedef std::pair<adj_iter, adj_iter> adjrange_t;
  typedef property_map<graph_type, vertex_index_t>::type IndexMap;
  typedef graph_traits<graph_type>::vertex_descriptor v_descriptor;
  typedef std::vector<v_descriptor> vert_vec;
  //  typedef std::vector<vert_vec> avt_type;
  typedef std::queue<v_descriptor> vert_que;
  typedef std::vector<bool> colormap;


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
  for(vpair.first; vpair.first != vpair.second; ++vpair.first)
  {
    xadj[i+1] = xadj[i] + out_degree(*vpair.first, graph1);
    adjrange_t adjpair = adjacent_vertices(*vpair.first, graph1);
    j = xadj[i];
    for(adjpair.first; adjpair.first != adjpair.second; ++adjpair.first)
    {
      adjncy[j] = *adjpair.first;
      ++j;
    }
    ++i;
  }

  // Print out to verify
  cout << "METIS partitioning: " << endl;
  for(i = 0; i < nvert + 1; ++i)
  {
    cout << xadj[i] << ' ';
  }
  cout << endl;
  for(i = 0; i < nedge * 2; ++i)
  {
    cout << adjncy[i] << ' ';
  }
  cout << endl;


  // Use k-way graph partition algorithm
  METIS_PartGraphKway(&nvert, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts, NULL, NULL, NULL, &objval, part);

  // Print out the partition vector
  for(i=0; i < nvert; ++i)
  {
    cout << part[i] << ' ';
  }
  cout << endl;

  // Create induced subgraphs based on partition vector
  graph_type *subgraph_vect = new graph_type[K];
  IndexMap index = get(vertex_index, graph1);
  for (i=0; i<K; ++i) 
  {
    subgraph_vect[i] = graph1.create_subgraph();
  }
  graph_type::children_iterator ci, ci_end;
  for (i=0; i < nvert; ++i)
  {
    cout << "Mapping " << index[i] << " to subgraph " << part[i] << endl;
    add_vertex(index[i], subgraph_vect[part[i]]);
  }
  /*
  j = 0;
  for (  boost::tie(ci, ci_end) = graph1.children(); ci != ci_end; ++ci)
  {
    for (i=0; i < nvert; ++i)
    {
      if(part[i] == j)
      {
        add_vertex(index[i], *ci);
      }
    }
    ++j;
  }
  */

  // print for testing

  //This will print out the global vertex IDs from the root graph
  vertex_iter v, v_end;
  edge_iter e, e_end;

  cout << "root:" << endl;
  print_graph(graph1, get(vertex_index, graph1));
  cout << endl;
  for (i=0; i < K; ++i)
  {
    cout << "subgraph " << i << ":" << endl;
    cout << "vertices = ";
    for (boost::tie(v, v_end) = vertices(subgraph_vect[i]); v != v_end; ++v)
    {
      cout << subgraph_vect[i].local_to_global(*v) << ", ";
    }
    cout << endl;
    cout << "edges = ";
    for (boost::tie(e, e_end) = edges(subgraph_vect[i]); e != e_end; ++e)
    {
      cout << subgraph_vect[i].local_to_global(*e) << ", ";
    }
    cout << endl;
  }

  // This does the same thing but will print out the local vertex
  // IDs for each subgraph
/*
  for (i=0; i < K; ++i)
  {
    cout << "subgraph " << i << ":" << endl;
    print_graph(subgraph_vect[i], get(vertex_index, subgraph_vect[i]));
  }
  cout << endl;
  int num = 1;
  vertex_iter v, v_end;
  for (boost::tie(ci, ci_end) = graph1.children(); ci != ci_end; ++ci)
  {
    cout << "G" << num++ << ": " << endl;
    cout << "vertices = ";
    for (boost::tie(v, v_end) = vertices(*ci); v != v_end; ++v)
    {
      cout << ci->local_to_global(*v) << ", ";
    }
    cout << endl;
    print_graph(*ci, get(vertex_index, *ci));
    cout << endl;
  }
*/

  // Start to build the AVT
  // Each vector will represent a column in the table
  // Each column in the table represents a subgraph
  cout << "cp1" << endl;
  vert_vec *avt = new vert_vec[K];
  vert_vec *avt_unmatched = new vert_vec[K];  //to store intermediate results
  
  // Build colormap
  cout << "cp2" << endl;

  colormap *clr_arr = new colormap[K];

  // First step is to create the initial row
  cout << "cp3" << endl;
  v_descriptor *avtrow = new v_descriptor[K];
  for (i=0; i < K; ++i)
  {
    int maxdegree = 0;
    v_descriptor vertID = 0;
    cout << "cp4" << endl;

    for (boost::tie(v, v_end) = vertices(subgraph_vect[i]); v != v_end; ++v)
    {
      int temp = out_degree(*v, subgraph_vect[i]);
      if (temp > maxdegree)
      {
        maxdegree = temp;
	    vertID = subgraph_vect[i].local_to_global(*v);
	    avtrow[i] = vertID;
      }
      cout << "maxdegree = " << maxdegree << endl;
      // while we are processing every vertex in each subgraph
      // initialise the colormap
      clr_arr[i].push_back(false) ;
    }
  }
    
  //Print out the color array
  for(i=0; i<K; ++i)
  {
       for(j=0; j < clr_arr[i].size(); ++j)
       {
            cout << clr_arr[i][j] << ' ';
       }
  }
  cout << "cp5" << endl;

  // Print out the starting point of AVT
  cout << endl << "AVT:" <<endl;
  for (i=0; i<K; ++i)
  {
    cout << avtrow[i] << ' ';
  }
  
  // Load in the first row
  // TODO: do this lol
  for(i=0; i<K; ++i)
    {
      avt[i].push_back(avtrow[i]);
    }

  // build ques for BFS
  vert_que *vque_arr = new vert_que[K];
  for(i=0; i<K; ++i)
  {
    vque_arr[i].push(avtrow[i]);
  }

  // Process using BFS
  cout << "cp6" << endl;
  // Mark all starting nodes as visited
  for (i=0; i < K; ++i)
    {
      clr_arr[i][0] = true;
    }
  // Process the que
  for(i=0; i < K; ++i)
    {
      while (vque_arr[i].empty()==false)
	{
	  v_descriptor vertID = vque_arr[i].front();
	  vque_arr[i].pop();
	  adj_iter vi, vi_end;
	  for(boost::tie(vi, vi_end) = adjacent_vertices(vertID, graph1); vi != vi_end; ++v)
	    {
	      if(clr_arr[i][*vi] == false)
		{
		  vque_arr[i].push(*vi);
		  clr_arr[i][*vi] = true;
		}
	    }
	}
    }

  // Print out the AVT
  int avtRows = avt[0].size();
  
  cout << endl << "AVT:" << endl;
  for(i=0; i < avtRows; ++i)
  {
    for(j=0; j<K; ++j)
    {
      cout << avt[j][i] << ' '; 
    }
    cout << endl;
  }

  cout << endl;


  //cleanup
  delete[] xadj;
  delete[] adjncy;
  delete[] part;
  delete[] subgraph_vect;
  delete[] avt;
  delete[] avtrow;
  delete[] clr_arr;
  delete[] vque_arr;
  delete[] avt_unmatched;
  xadj = NULL;
  adjncy = NULL;
  part = NULL;
  subgraph_vect = NULL;
  avt = NULL;
  avtrow = NULL;
  clr_arr = NULL;
  vque_arr = NULL;
  avt_unmatched = NULL;


  return 0;
}
