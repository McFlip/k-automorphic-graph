part_test.exe: test_part.o
	g++ -o part_test.exe test_part.o -L /usr/local/lib/ -lboost_graph \/usr/local/lib/libmetis.so
test_part.o: test_part.cpp
	g++ -std=c++11 -ggdb -o test_part.o -c test_part.cpp
clean:
	rm -f *.o part_test.exe

