#include <iostream>
#include <fstream>
using namespace std;

char graphMatrix[325729][40717];

void parse(char* argv[]);
void printMatrix();

int main(int argc, char* argv[])
{
	//Zero out the graph matrix
	for(int i = 0; i<325729; ++i)
	{
		for(int j = 0; j<40717; ++j)
		{
			graphMatrix[i][j] = '\0';
		}
	}
	
	//build graph matrix from input file
	parseWebNotre(argv);
	
	printMatrix();
	
	return 0;
}

void parseWebNotre(char* argv[])
{
	char trashBuffer[256];
	int fromNode;
	int toNode;
	char maskArray[8];
	
	//initialize mask array
	maskArray[0] = 0x80;
	maskArray[1] = 0x40;
	maskArray[2] = 0x20;
	maskArray[3] = 0x10;
	maskArray[4] = 0x08;
	maskArray[5] = 0x04;
	maskArray[6] = 0x02;
	maskArray[7] = 0x01;
	
	//open input file
	ifstream inputFile;
	inputFile.open(argv[1]);
	
	inputFile.getline(trashBuffer, 256);
	inputFile.getline(trashBuffer, 256);
	inputFile.getline(trashBuffer, 256);
	inputFile.getline(trashBuffer, 256);
	
	while(inputFile >> fromNode >> toNode)
	{
		char mask = maskArray[toNode % 8];
		graphMatrix[fromNode][toNode/8] |= mask;
		mask = maskArray[fromNode % 8];
		graphMatrix[toNode][fromNode/8] |= mask;
	}
	
}

void printMatrix()
{
	char maskArray[8];
	maskArray[0] = 0x80;
	maskArray[1] = 0x40;
	maskArray[2] = 0x20;
	maskArray[3] = 0x10;
	maskArray[4] = 0x08;
	maskArray[5] = 0x04;
	maskArray[6] = 0x02;
	maskArray[7] = 0x01;
	
	for(int i = 0; i < 16; ++i)
	{
		for(int j = 0; j < 16; ++j)
		{
			if((graphMatrix[j][i/8] & maskArray[i%8]) == 0)
				cout << "0";
			else
				cout << "1";
			cout << ' ';
		}
		cout << '\n';
	}
}

