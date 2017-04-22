// Create callback to print mappings
//  vf2_print_callback<graph_type, graph_type> callback(graph1, graph2);
//  vf2_subgraph_iso(graph1, graph2, callback);

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
 
 //Old if else for adding vertices
 /*
 if (avt_unmatched[i].size() < avt_unmatched[i+1].size())
 {
     vLocalID = add_vertex(subgraph_vect[i]);
     avt_unmatched[i].push_back(subgraph_vect[i].local_to_global(vLocalID));
 }
 else if (avt_unmatched[i].size() > avt_unmatched[i+1].size())
 {
     vLocalID = add_vertex(subgraph_vect[i+1]);
     avt_unmatched[i+1].push_back(subgraph_vect[i+1].local_to_global(vLocalID));
 }*/
 
 //*** Add Vertices at each "wave" from starting point ***
 /*typedef std::vector< int > score_vec;

 //Insert necessary noise vertexs into the AVT table
 int current_hopcount;
 std::vector<hop_pair_t>::iterator it;
 std::vector<hop_pair_t>::iterator it2;
 it = avt_unmatched[0].begin();
 
 //for each item in the first AVT column
 for(i = 0; it+i != avt_unmatched[0].end(); ++i)
 {
     //set the hopcount equal to the items hopcount
     current_hopcount = (it+i)->second;
     //for each entry in that items row
     for(int j = 1; j < K; ++j)
     {
         it2 = avt_unmatched[j].begin()+i;
         if( it2->second > current_hopcount)
         {
             avt_unmatched[j].insert(it2, std::make_pair(boost::add_vertex(subgraph_vect[i]), hopcount));
         }
         else( it2->second < current_hopcount)
         {
             current_hopcount = it2->second;
             j = 0;
         }
     }
 }*/
 
 //*** Assign a global score to each vertex ***
 //global score optimization
 /*score_vec *score_table = new score_vec[K];
 int global_score;
 for(i=0; i < K; ++i)
 {
     avtRows = avt_unmatched[i].size();
     for(j = 1; j < avtRows; ++j)
     {
         global_score = out_degree(subgraph_vect[i].global_to_local(avt_unmatched[i][j].first), subgraph_vect[i]);
         global_score += (1000 * avt_unmatched[i][j].second);
         score_table[i].push_back(global_score);
     }
 }
 
 //*** Pair vertices based on global score
 std::vector< int >::iterator it;
 for(i = 0; i < K; ++i)
 {
     avt[i].push_back(avt_unmatched[i].first);
     avt_unmatched[i].erase(avt_unmatched[i].begin());
 }
 int bestMatch;
 int bestDiff;
 for(i = 0; i < K; ++i)
 {
     while(avt_unmatched[i].empty() == false)
     {
         for(j = 0; j < K; ++j)
         {
             if(j != i)
             {
                 bestDiff = 10000;
                 for(int zy = 0; zy < score_table[j].size(); ++zy)
                 {
                     if(bestDiff > std::abs(score_table[j][zy] - score_table[i][0]))
                     {
                         bestDiff = std::abs(score_table[j][zy] - score_table[i][0]);
                         bestMatch = zy;
                     }
                 }
                 if(bestDiff == 10000)
                 {
                     avt[j].push_back(boost::add_vertex(subgraph_vect[j]));
                 }
                 else
                 {
                     avt[j].push_back((avt_unmatched[j].begin())+bestMatch);
                     avt_unmatched[j].erase((avt_unmatched[j].begin())+bestMatch);
                     score_table[j].erase((score_table[j].begin())+bestMatch);
                 }
             }
         }
         avt[i].push_back(avt_unmatched[i].first);
         avt_unmatched[i].erase(avt_unmatched[i].begin());
         score_table[i].erase(score_table[i].begin());
     }
 }*/
