#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <utility>
#include <stdlib.h>
#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_utility.hpp>
#include "metis.h"

using namespace std;
using namespace boost;

// global declarations
const int K = 2;                                            // TODO: Change this to main param
idx_t ncon = 1;                                             // default METIS tuning param
idx_t nparts = K;                                           // set num partition to K


int main(int argc, char* argv[])
{
    //********************************** Declarations **********************************
	
    // typedef declarations for templates
    typedef subgraph< adjacency_list<vecS, vecS, undirectedS, uint32_t,
        property< edge_index_t, int > > > graph_type;
    typedef graph_traits<graph_type>::vertex_iterator vertex_iter;
    typedef graph_traits<graph_type>::edge_iterator edge_iter;
    typedef std::pair<vertex_iter, vertex_iter> vrange_t;
    typedef graph_traits<graph_type>::adjacency_iterator adj_iter;
    typedef std::pair<adj_iter, adj_iter> adjrange_t;
    typedef property_map<graph_type, vertex_index_t>::type IndexMap;
    typedef graph_traits<graph_type>::vertex_descriptor v_descriptor;
    typedef std::vector<v_descriptor> vert_vec;
    typedef std::pair< v_descriptor, int > hop_pair_t;
    typedef std::vector<hop_pair_t> avt_vector_t;
    typedef std::queue<hop_pair_t> vert_que;
    typedef std::vector<bool> colormap;
    typedef std::vector<int> score_vec;

    //*** variable declarations ***

    // The root graph
    graph_type graph1(0);

    // subgraphs from the root
    graph_type *subgraph_vect = new graph_type[K];

    // METIS library args
    idx_t nvert;
    int nedge;
    idx_t *xadj;
    idx_t *adjncy;
    idx_t objval;
    idx_t *part;
    idx_t options[METIS_NOPTIONS];

    // The graphviz dot file
    ifstream inputFile;

    // iterators
    vrange_t vpair;
    int i = 0;
    int j = 0;
    graph_type::children_iterator ci;
    graph_type::children_iterator ci_end;
    vertex_iter v;
    vertex_iter v_end;
    //    vertex_iter v_next_col;
    edge_iter e;
    edge_iter e_end;
    adj_iter vi;
    adj_iter vi_end;

    // index for random access
    IndexMap index;

    // Alignment Vertex table
    vert_vec *avt = new vert_vec[K];
    avt_vector_t *avt_unmatched = new avt_vector_t[K];      // to store intermediate results
    v_descriptor *avtrow = new v_descriptor[K];             // the initial row
    int avtRows;                                            // number of rows in table
    int hopcount;
    int degrees;
    int score;
    int bestscore;
    int pos_right;
    int best_position;
    int maxlength;
    score_vec *score_table = new score_vec[K];
    v_descriptor bestV;
    hop_pair_t hopPair;
    
    // colormap & que used in BFS
    colormap *clr_arr = new colormap[K];
    vert_que *vque_arr = new vert_que[K];

    // descriptors
    v_descriptor vLocalID;
    v_descriptor vGlobalID;

	
    //********************************** Reading and Building Graph **********************************
    // Read in the graph from disk
    inputFile.open(argv[1]);
    dynamic_properties dp(ignore_other_properties);
    read_graphviz(inputFile, graph1, dp);
    inputFile.close();

    //Initialize Graph-Dependent Variables
    nvert = num_vertices(graph1);
    nedge = num_edges(graph1);
    xadj = new idx_t[nvert + 1];
    adjncy = new idx_t[nedge * 2];
    part = new idx_t[nvert];
    index = get(vertex_index, graph1);
    
	
    //********************************** Partition the graph with METIS algorithm **********************************
    // initialise
    METIS_SetDefaultOptions(options);
    xadj[0] = 0;

    // Transform the adjacency list into the Compressed Storage Format
    // TODO: consider using tie function - syntactic sugar
    for(vpair = vertices(graph1); vpair.first != vpair.second; ++vpair.first)
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
    objval = 0;
    METIS_PartGraphKway(&nvert, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts, NULL, NULL, NULL, &objval, part);

    // Print out the partition vector
    for(i=0; i < nvert; ++i)
    {
        cout << part[i] << ' ';
    }
    cout << endl;

	//********************************** Create Subgraphs based on METIS result **********************************
    // Create induced subgraphs based on partition vector
    for (i=0; i<K; ++i)
    {
        subgraph_vect[i] = graph1.create_subgraph();
    }
    for (i=0; i < nvert; ++i)
    {
        cout << "Mapping " << index[i] << " to subgraph " << part[i] << endl;
        add_vertex(index[i], subgraph_vect[part[i]]);
    }
	
	// print for testing

    //This will print out the global vertex IDs from the root graph
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

    
    //********************************** Building the AVT **********************************
    // Start to build the AVT
    // Each vector will represent a column in the table
    // Each column in the table represents a subgraph

    // Build colormap
    // First step is to create the initial row
    for (i=0; i < K; ++i)
    {
        int maxdegree = 0;
        v_descriptor vertID = 0;
        // cout << "cp4" << endl;

        for (boost::tie(v, v_end) = vertices(subgraph_vect[i]); v != v_end; ++v)
        {
            int temp = out_degree(*v, subgraph_vect[i]);
            if (temp > maxdegree)
            {
                maxdegree = temp;
                vertID = *v;
                avtrow[i] = vertID;
            }
            // cout << "maxdegree = " << maxdegree << endl;
            // while we are processing every vertex in each subgraph
            // initialise the colormap
            clr_arr[i].push_back(false) ;
        }
    }

    //Print out the color array
    cout << std::boolalpha;
    for(i=0; i<K; ++i)
    {
        for(j=0; j < clr_arr[i].size(); ++j)
        {
            cout << clr_arr[i][j] << ' ';
        }
        cout << endl;
    }

    // Print out the starting point of AVT
    cout << endl << "AVT:" <<endl;
    for (i=0; i<K; ++i)
    {
        cout << subgraph_vect[i].local_to_global(avtrow[i]) << ' ';
    }
    cout << endl;
    
    // Load in the first row
    for(i=0; i<K; ++i)
    {
        avt_unmatched[i].push_back(std::make_pair(avtrow[i], 0));
    }

    // build ques for BFS
    for(i=0; i<K; ++i)
    {
        vque_arr[i].push(std::make_pair(avtrow[i], 0));
    }

    // Process using BFS
    // Mark all starting nodes as visited
    for (i=0; i < K; ++i)
    {
        index = get(vertex_index, subgraph_vect[i]);

        clr_arr[i][index[avtrow[i]]] = true;
    }

    // Print  out the color array
    cout << std::boolalpha;
    for(i=0; i<K; ++i)
    {
        for(j=0; j < clr_arr[i].size(); ++j)
        {
            cout << clr_arr[i][j] << ' ';
        }
        cout << endl;
    }

    // Process the que
	//*** Breadth First Search ***
    for(i=0; i < K; ++i)
    {
        while (vque_arr[i].empty()==false)
        {
            //  v_descriptor vGlobalID = vque_arr[i].front();
            //  v_descriptor vLocalID = subgraph_vect[i].global_to_local(vGlobalID);
            hopPair = vque_arr[i].front();

            vLocalID = hopPair.first;
            hopcount = hopPair.second;
            //  cout << "Processing BFS node: " << vLocalID << endl;
            cout << "Processing BFS node: " << index[vLocalID] << endl;
            index = get(vertex_index, subgraph_vect[i]);
            vque_arr[i].pop();
            for(boost::tie(vi, vi_end) = adjacent_vertices(vLocalID, subgraph_vect[i]); vi != vi_end; ++vi)
            {
                //  cout << subgraph_vect[i][*vi] << endl;
                //  cout << graph1[*vi] << endl;
                if(clr_arr[i][index[*vi]] == false)
                {
                    cout << "Processing adjacent node: " << index[*vi] << endl;
                    vque_arr[i].push(std::make_pair(*vi, hopcount+1));
                    avt_unmatched[i].push_back(std::make_pair(*vi, hopcount+1));
                    clr_arr[i][index[*vi]] = true;
                    // Print out the color map as updated
                    for(int z=0; z<K; ++z)
                    {
                        for(j=0; j < clr_arr[z].size(); ++j)
                        {
                            cout << clr_arr[z][j] << ' ';
                        }
                        cout << endl;
                    }
                }
            }
        }
    }
    
    // Print out the AVT
    cout << endl << "AVT (unmatched - UNbalanced) <read this sideways>:" << endl;
    for(i=0; i < K; ++i)
    {
        avtRows = avt_unmatched[i].size();
        for(j=0; j < avtRows; ++j)
        {
            cout << subgraph_vect[i].local_to_global(avt_unmatched[i][j].first) << '(' << avt_unmatched[i][j].second << ") ";
        }
        cout << endl;
    }
    cout << endl;

    // balance out the table if odd number of vertices
    // find the max length column
    maxlength = 0;
    for (i=0; i<K; ++i)
    {
        if (avt_unmatched[i].size() > maxlength)
        {
            maxlength = avt_unmatched[i].size();
        }
    }
	
    //Add vertices until all columns are maxlength
    for (i=0; i<K; ++i)
    {
		while (avt_unmatched[i].size() < maxlength)
    	{
	    	vLocalID = add_vertex(subgraph_vect[i]);
			avt_unmatched[i].push_back(std::make_pair(vLocalID, 10000));
			clr_arr[i].push_back(false);
    	}
    }

	
    // Print out the AVT
    cout << endl << "AVT (unmatched - balanced) <read this sideways>:" << endl;
    for(i=0; i < K; ++i)
    {
        avtRows = avt_unmatched[i].size();
        for(j=0; j < avtRows; ++j)
        {
            cout << subgraph_vect[i].local_to_global(avt_unmatched[i][j].first) << '(' << avt_unmatched[i][j].second << ") ";
        }
        cout << endl;
    }
    cout << endl;
    
    
    //********************************** using local scores **********************************
    
    // use colormap for tracking matches
    // initialise
    for(i=0; i<K; ++i)
    {
        for(j=0; j < clr_arr[i].size(); ++j)
        {
            if(j==0)
            {
                clr_arr[i][j] = true;
            }
            else
            {
                clr_arr[i][j] = false;
            }
        }
        cout << endl;
    }
    
    // Print  out the color array
    for(i=0; i<K; ++i)
    {
        for(j=0; j < clr_arr[i].size(); ++j)
        {
            cout << clr_arr[i][j] << ' ';
        }
        cout << endl;
    }
    
    // copy the first column as is
    for (i=0; i < avt_unmatched[0].size(); ++i)
    {
        avt[0].push_back(avt_unmatched[0][i].first);
    }
    // match from left to right by score
    // score of 0 is considered a perfect match
    // as in same degree and same distance from hub
    for (i=0; i < K-1; ++i)
    {
        // copy the first row as is
        avt[i+1].push_back(avt_unmatched[i+1][0].first);
        
        // skip over first iteration so j=1 and begin()+1
        for(j=1; j < avt[i].size(); ++j)
        {
            cout << "avt size=" << avt[i].size() << endl;
            bestscore = 10000;
            pos_right = 0;
            best_position = 0;
            for(auto v_next_col = avt_unmatched[i+1].begin(); v_next_col != avt_unmatched[i+1].end(); ++v_next_col)
            {
                if(clr_arr[i+1][pos_right] == false)
                {
                    int dleft = out_degree(avt_unmatched[i][j].first, subgraph_vect[i]);
                    int dright = out_degree(v_next_col->first, subgraph_vect[i+1]);
                    cout << "comparing " << subgraph_vect[i].local_to_global(avt_unmatched[i][j].first) << " to " << subgraph_vect[i+1].local_to_global(v_next_col->first) << endl;
                    cout << "degree left: " << dleft
                        << " degree right: " << dright
                        << endl;
                    degrees = std::abs(dleft - dright);
                    hopcount = std::abs(avt_unmatched[i][j].second - v_next_col->second);
                    score = degrees + hopcount;
                    cout << "degree diff: " << degrees << " hopcount diff: " << hopcount << " score: " << score << endl;
                    if(score < bestscore)
                    {
                        bestscore = score;
                        vGlobalID = v_next_col->first;
                        best_position = pos_right;
                    }
                }
                ++pos_right;
            }
            clr_arr[i+1][best_position] = true;
            avt[i+1].push_back(vGlobalID);
        }
        
    }

    // Print  out the color array
    for(i=0; i<K; ++i)
    {
        for(j=0; j < clr_arr[i].size(); ++j)
        {
            cout << clr_arr[i][j] << ' ';
        }
        cout << endl;
    }

    
    // Print out the AVT

    cout << endl << "AVT <read this sideways>:" << endl;
    for(i=0; i < K; ++i)
    {
        avtRows = avt[i].size();
        for(j=0; j < avtRows; ++j)
        {
            cout << subgraph_vect[i].local_to_global(avt[i][j]) << ' ';
        }
        cout << endl;
    }
    cout << endl;
	
	//********************************** Perform Block Alignment **********************************
	
	
	
	
	//********************************** Perform Edge Copy **********************************
    
    
    // *******End***********
    //cleanup
    // TODO: move these higher
    delete[] xadj;
    delete[] adjncy;
    delete[] part;
    delete[] subgraph_vect;
    delete[] avt;
    delete[] avtrow;
    delete[] clr_arr;
    delete[] vque_arr;
    delete[] avt_unmatched;
    delete[] score_table;
    xadj = NULL;
    adjncy = NULL;
    part = NULL;
    subgraph_vect = NULL;
    avt = NULL;
    avtrow = NULL;
    clr_arr = NULL;
    vque_arr = NULL;
    avt_unmatched = NULL;
    score_table = NULL;


    return 0;
}
