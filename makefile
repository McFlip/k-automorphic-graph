all: automorph.exe main.exe
automorph.exe: automorph.o
	g++ -o automorph.exe automorph.o -L /usr/local/lib/ -lboost_graph
automorph.o: test_k.cpp
	g++ -std=c++11 -ggdb -o automorph.o -c test_k.cpp
main.exe: main.o
	g++ -o main.exe main.o -L /usr/local/lib/ -lboost_graph \/usr/local/lib/libmetis.so
main.o: test_part.cpp
	g++ -std=c++14 -ggdb -o main.o -c test_part.cpp
clean:
	rm -f *.o main.exe automorph.exe
