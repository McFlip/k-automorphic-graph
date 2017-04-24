# k-automorphic-graph
Advanced DB project consisting of privacy preserving graphs and queries.

1.
WebND_preprocess.exe path-to-web-NotreDame.txt
-preprocesses the text file format that the DB is written to
-outputs a dot file that is used by the other programs

2.
automorph.exe path-to-dot-file
-prints out the automorphic functions of given graph

3.
main.exe path-to-input-dot-file path-to-output-dot-file
-the main program
-produces the k-anonymized version of the graph
-saves the new graph to a dot file

4.
run automorph again and compare results to previous run


Compiling flags required:

anything with boost -I path/to/boost i.e. /usr/local/boost_1_63_0/boost/

The Boost C++ Libraries were successfully built!

The following directory should be added to compiler include paths:

/usr/local/boost_1_63_0

The following directory should be added to linker library paths:

/usr/local/boost_1_63_0/stage/lib

g++ -I /usr/local/boost_1_63_0/ -std=c++11 source.cpp -L/usr/local/boost_1_63_0/stage/lib/ -lboost_graph

METIS library:

g++ -I /usr/local/boost_1_63_0/ /usr/local/metis-5.1.0/include/ -o test_part.o -c test_part.cpp
g++ -o part_test.exe test_part.o -L /usr/local/lib/ -lboost_grap
g++ -o part_test.exe test_part.o -L /usr/local/lib/ -lboost_graph \/usr/local/lib/libmetis.so


notes:
web-NotreDame has Nodes: 325,729 Edges: 1,497,134

to test run time and space use the time command with full path
$ /usr/bin/time -v ./test.exe web-NotreDame.txt
