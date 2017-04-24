all: automorph.exe main.exe webND_preprocess.exe rando_graph_gen.exe
automorph.exe: automorph.o
	g++ -o automorph.exe automorph.o -L /usr/local/lib/ -lboost_graph
automorph.o: test_k.cpp
	g++ -std=c++11 -ggdb -o automorph.o -c test_k.cpp
main.exe: main.o
	g++ -o main.exe main.o -L /usr/local/lib/ -lboost_graph \/usr/local/lib/libmetis.so
main.o: main.cpp
	g++ -std=c++14 -ggdb -o main.o -c main.cpp
webND_preprocess.exe: webND_preprocess.o
	g++ -o webND_preprocess.exe webND_preprocess.o -L /usr/local/lib/ -lboost_graph
webND_preprocess.o: webND_preprocess.cpp
	g++ -std=c++11 -ggdb -o webND_preprocess.o -c webND_preprocess.cpp
rando_graph_gen.exe: rando_graph_gen.o
	g++ -o rando_graph_gen.exe rando_graph_gen.o
rando_graph_gen.o: rando_graph_gen.cpp
	g++ -std=c++11 -ggdb -o rando_graph_gen.o -c rando_graph_gen.cpp
clean:
	rm -f *.o main.exe automorph.exe webND_preprocess.exe
