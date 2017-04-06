# k-automorphic-graph
Advanced DB project consisting of privacy preserving graphs and queries.
Compiling flags required:
main.cpp -mcmodel=medium
anything with boost -I path/to/boost i.e. /usr/local/boost_1_63_0/boost/

notes:
web-NotreDame has Nodes: 325,729 Edges: 1,497,134
Relying on compiler to do initialization of the matrix dropped execution time by 30 seconds.
Going to test boost adjacency list to measure memmory difference.
$ /usr/bin/time -v ./test.exe web-NotreDame.txt
        Command being timed: "./test.exe web-NotreDame.txt"
        User time (seconds): 0.19
        System time (seconds): 1.36
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.56
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 12953704
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 13289
        Voluntary context switches: 11
        Involuntary context switches: 3
        Swaps: 0
        File system inputs: 0
        File system outputs: 0
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0

Time results from boost adjacency list:
$ /usr/bin/time -v ./adjlist.exe web-NotreDame.txt
        Command being timed: "./adjlist.exe web-NotreDame.txt"
        User time (seconds): 2.36
        System time (seconds): 0.06
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.43
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 293312
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 65059
        Voluntary context switches: 11
        Involuntary context switches: 4
        Swaps: 0
        File system inputs: 0
        File system outputs: 0
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0