#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <utility>
#include <stdlib.h>
#include <unordered_set>
#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_utility.hpp>
#include "metis.h"

// Custom hash function for hash table
template <class Edge>
struct myhash
{
    std::size_t operator()(Edge const& e) const
    {
        return _h(e.idx);
    }
};


using namespace std;
using namespace boost;


int main(int argc, char* argv[])
{
    // error check
    if (argc != 3)
    {
        cout << "usage: main.exe inputFile outputFile" << endl;
        cout << "inputFile must be in graphviz dot file format." << endl;
        return 1;
    }
    
    //********************************** Declarations **********************************
	
	// global declarations
    const int K = 2;
    idx_t ncon = 1;                                             // default METIS tuning param
    idx_t nparts = K;                                           // set num partition to K
    
    
    // typedef declarations for templates
    typedef subgraph< adjacency_list<vecS, vecS, undirectedS, uint32_t,
        property< edge_index_t, int > > > graph_type;
    typedef graph_traits<graph_type>::vertex_iterator vertex_iter;
    typedef graph_traits<graph_type>::edge_iterator edge_iter;
    typedef graph_traits<graph_type>::out_edge_iterator out_edge_iter;
    typedef std::pair<vertex_iter, vertex_iter> vrange_t;
    typedef graph_traits<graph_type>::adjacency_iterator adj_iter;
    typedef std::pair<adj_iter, adj_iter> adjrange_t;
    typedef property_map<graph_type, vertex_index_t>::type IndexMap;
    typedef graph_traits<graph_type>::vertex_descriptor v_descriptor;
    typedef graph_traits<graph_type>::edge_descriptor e_descriptor;
    typedef std::vector<v_descriptor> vert_vec;
    typedef std::pair< v_descriptor, int > hop_pair_t;
    typedef std::vector<hop_pair_t> avt_vector_t;
    typedef std::queue<hop_pair_t> vert_que;
    typedef std::vector<bool> colormap;
    typedef std::vector<int> score_vec;
    typedef std::pair<v_descriptor, bool> found_t;
    //typedef std::unordered_set<e_descriptor, typename myhash> u_set_t;
    typedef std::set<e_descriptor> u_set_t;

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
    
    // Output graphviz dot file
    ofstream outputFile;
    // ostream& os = cout;

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
    out_edge_iter child_e;
    out_edge_iter child_e_end;
    out_edge_iter parent_e;
    out_edge_iter parent_e_end;
    adj_iter vi;
    adj_iter vi_end;

    // index for random access
    IndexMap index;

    // Alignment Vertex table
    vert_vec *avt = new vert_vec[K];
    v_descriptor *avt_lookup;
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
    
    // edge copy
    u_set_t edge_set;
    u_set_t added_edge_set;
    found_t found;
    int from;
    int to;

    // descriptors
    v_descriptor vLocalID;
    v_descriptor vGlobalID;
    v_descriptor orig_source;
    v_descriptor orig_target;
    v_descriptor copy_source;
    v_descriptor copy_target;
    e_descriptor eLocalID;
    e_descriptor eGlobalID;
    
    
    
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



    // Use k-way graph partition algorithm
    objval = 0;
    METIS_PartGraphKway(&nvert, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts, NULL, NULL, NULL, &objval, part);

    //********************************** Create Subgraphs based on METIS result **********************************
    // Create induced subgraphs based on partition vector
    for (i=0; i<K; ++i)
    {
        subgraph_vect[i] = graph1.create_subgraph();
    }
    for (i=0; i < nvert; ++i)
    {
        add_vertex(index[i], subgraph_vect[part[i]]);
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

        for (boost::tie(v, v_end) = vertices(subgraph_vect[i]); v != v_end; ++v)
        {
            int temp = out_degree(*v, subgraph_vect[i]);
            if (temp > maxdegree)
            {
                maxdegree = temp;
                vertID = *v;
                avtrow[i] = vertID;
            }
            // while we are processing every vertex in each subgraph
            // initialise the colormap
            clr_arr[i].push_back(false) ;
        }
    }


    
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
            index = get(vertex_index, subgraph_vect[i]);
            vque_arr[i].pop();
            for(boost::tie(vi, vi_end) = adjacent_vertices(vLocalID, subgraph_vect[i]); vi != vi_end; ++vi)
            {
                if(clr_arr[i][index[*vi]] == false)
                {
                    vque_arr[i].push(std::make_pair(*vi, hopcount+1));
                    avt_unmatched[i].push_back(std::make_pair(*vi, hopcount+1));
                    clr_arr[i][index[*vi]] = true;
                    
                }
            }
        }
    }
    
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
			avt_unmatched[i].push_back(std::make_pair(vLocalID, 99));
			clr_arr[i].push_back(false);
    	}
    }

	
 
    
    
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
    }
    

    
    // create the avt_lookup table
    // We used a 2d array flattened to 1d
    // The values in this table are the index positions in the avt
    // Every time we push into avt we update avt_lookup
    avt_lookup = new v_descriptor[K * avt_unmatched[0].size() * K];
    
    // copy the first column as is
    
    //use index to index into avt_lookup
    index = get(vertex_index, subgraph_vect[0]);
    for (i=0; i < avt_unmatched[0].size(); ++i)
    {
        avt[0].push_back(avt_unmatched[0][i].first);
        avt_lookup[index[avt_unmatched[0][i].first]] = i;
    }
    // match from left to right by score
    // score of 0 is considered a perfect match
    // as in same degree and same distance from hub
    
    // for each column, fill in next column row by row
    for (i=0; i < K-1; ++i)
    {
        // copy the first row as is
        avt[i+1].push_back(avt_unmatched[i+1][0].first);
        // get new index for each subgraph
        index = get(vertex_index, subgraph_vect[i+1]);
        avt_lookup[(i+1) * avt_unmatched[0].size() + (index[avt_unmatched[i+1][0].first])] = i;
        
        // for each item in left column, match against items in right column
        // skip over first iteration so j=1
        for(j=1; j < avt[i].size(); ++j)
        {
            bestscore = 10000;  // initialise to silly high number; track best score so far
            pos_right = 0;      // current position of item on right
            best_position = 0;  // save position of best match
            // for each item in right column see if that item is a good match
            for(auto v_next_col = avt_unmatched[i+1].begin(); v_next_col != avt_unmatched[i+1].end(); ++v_next_col)
            {
                if(clr_arr[i+1][pos_right] == false)
                {
                    int dleft = out_degree(avt_unmatched[i][j].first, subgraph_vect[i]);
                    int dright = out_degree(v_next_col->first, subgraph_vect[i+1]);
                    degrees = std::abs(dleft - dright);
                    hopcount = std::abs(avt_unmatched[i][j].second - v_next_col->second);
                    score = degrees + hopcount;  // can change to weighted score
                    if(score < bestscore)
                    {
                        bestscore = score;
                        vLocalID = v_next_col->first;
                        best_position = pos_right;
                    }
                }
                ++pos_right;
            }
            
            // update the colormap, push in to avt, update lookup table
            clr_arr[i+1][best_position] = true;
            avt[i+1].push_back(vLocalID);
            avt_lookup[(i + 1) * avt_unmatched[i].size() + (index[vLocalID])] = j;
        }
        
    }


	
	//********************************** Perform Block Alignment **********************************

	//*** First add all edges to the first column ***
	//*** first column will become the rolemodel ***
	adj_iter one;
	adj_iter one_end;
	adj_iter two;
	adj_iter two_end;
	bool matchFound;
	//*** Make the first column the rolemodel with all edges copied ***
	for(i = 0; i < avt[0].size(); ++i)
	{
		for(j = 1; j < K; ++j)
		{
			index = get(vertex_index ,subgraph_vect[j]);
			for(boost::tie(two, two_end) = adjacent_vertices(avt[j][i], subgraph_vect[j]); two != two_end; ++two)
			{
				//find pair to vertex in column one
				v_descriptor pair_vertex = index[avt[0][avt_lookup[j * avt[0].size() + (*two)]]];
				matchFound = false;
				for(boost::tie(one, one_end) = adjacent_vertices(avt[0][i], subgraph_vect[0]); one != one_end; ++one)
				{
					if(pair_vertex == *one)
					{
						matchFound = true;
					}
				}
				if(matchFound == false)
					add_edge(avt[0][i], pair_vertex, subgraph_vect[0]);
			}
		}
	}
	//*** compare all other columns to the rolemodel column and add edges ***
	for(i = 0; i < avt[0].size(); ++i)
	{
		for(j = 1; j < K; ++j)
		{
			index = get(vertex_index ,subgraph_vect[0]);
			for(boost::tie(one, one_end) = adjacent_vertices(avt[0][i], subgraph_vect[0]); one != one_end; ++one)
			{
				v_descriptor pair_vertex = index[avt[j][avt_lookup[*one]]];
				matchFound = false;
				for(boost::tie(two, two_end) = adjacent_vertices(avt[j][i], subgraph_vect[j]); two != two_end; ++two)
				{
					if(pair_vertex == *two)
					{
						matchFound = true;
					}
				}
				if(matchFound == false)
					add_edge(avt[j][i], pair_vertex, subgraph_vect[j]);
			}
		}
	}
	

	
	//********************************** Perform Edge Copy **********************************
    // for each child subgraph
    for (i=0; i < K; ++i)
    {
        // for each vertex in subgraph, compare degree to parent vertex
        for(boost::tie(v, v_end) = vertices(subgraph_vect[i]); v != v_end; ++v)
        {
            // if degrees don't match then there is at least one crossing edge
            if (out_degree(*v, subgraph_vect[i]) != out_degree(subgraph_vect[i].local_to_global(*v), graph1))
            {
                // for each child edge, load into hash table
                edge_set.clear();
                for(boost::tie(child_e, child_e_end) = out_edges(*v, subgraph_vect[i]); child_e != child_e_end; ++child_e)
                {
                    eGlobalID = subgraph_vect[i].local_to_global(*child_e);
                    edge_set.insert(eGlobalID);
                }
                // for each edge of parent vertex, check if it's in hash table
                for(boost::tie(parent_e, parent_e_end) = out_edges(subgraph_vect[i].local_to_global(*v), graph1); parent_e != parent_e_end; ++ parent_e)
                {
                    if(edge_set.find(*parent_e) == edge_set.end()
                        && added_edge_set.find(*parent_e) == added_edge_set.end())
                    {
                        // do edge copy
                        // get the source & target of the edge
                        orig_source = source(*parent_e, graph1);
                        orig_target = target(*parent_e, graph1);
                        // find them in the avt
                        // use avt_lookup to find the index of the source
                        index = get(vertex_index, subgraph_vect[i]);
                        vLocalID = subgraph_vect[i].global_to_local(orig_source);
                        from = avt_lookup[i * avt[0].size() + index[vLocalID]];
                        // find the subgraph of the target & get it's local ID
                        for (j=0; j < K; ++j)
                        {
                            found = subgraph_vect[j].find_vertex(orig_target);
                            if(found.second)
                            {
                                // do lookup of the target
                                index = get(vertex_index, subgraph_vect[j]);
                                vLocalID = found.first;
                                to = avt_lookup[j * avt[0].size() + index[vLocalID]];
                                copy_source = subgraph_vect[j].local_to_global(avt[j][from]);
                                copy_target = subgraph_vect[i].local_to_global(avt[i][to]);
                                if(edge(copy_source, copy_target, graph1).second == false)
                                {
                                    add_edge(copy_source, copy_target, graph1);
                                    added_edge_set.insert(eGlobalID);
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // *******End***********
    
    
    // Dump the output to file
    outputFile.open(argv[2]);
    write_graphviz_dp(outputFile, graph1, dp.property("node_id", get(boost::vertex_index, graph1)));
    // uncomment this to print dot file to std out
    //write_graphviz_dp(cout, graph1, dp.property("node_id", get(boost::vertex_index, graph1)));
    outputFile.close();
    
    
    //cleanup
    // TODO: move these higher
    if(xadj != nullptr)
    {
        delete[] xadj;
    }
    if(adjncy != nullptr)
    {
        delete[] adjncy;
    }
    if(part != nullptr)
    {
        delete[] part;
    }
    if(subgraph_vect != nullptr)
    {
        delete[] subgraph_vect;
    }
    if(avt != nullptr)
    {
        delete[] avt;
    }
    if(avt_lookup != nullptr)
    {
        delete[] avt_lookup;
    }
    if(avtrow != nullptr)
    {
        delete[] avtrow;
    }
    if(clr_arr != nullptr)
    {
        delete[] clr_arr;
    }
    if(vque_arr != nullptr)
    {
        delete[] vque_arr;
    }
    if (avt_unmatched != nullptr)
    {
        delete[] avt_unmatched;
    }
    if(score_table != nullptr)
    {
        delete[] score_table;
    }
    
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
    avt_lookup = NULL;

    return 0;
}
