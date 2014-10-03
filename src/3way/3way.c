#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sonLib.h"
#include "3_Absorb3edge2x.h"

static int LINE_MAX = 5000;

// Parsing
static void parseLine(char * line, stList * vertices) {
	char * token;
	int64_t val;
	stList * vertex = stList_construct();

	token = strtok(line, "\t");
	while (token != NULL) {
		sscanf(token, "%i", &val);
		stList_append(vertex, stIntTuple_constructN(1, &val));
		token = strtok(NULL, "\t");
	}

	stList_append(vertices, vertex);
}

static stList * parseFile(char* filename) {
	FILE * file = fopen(filename, "r");
	if (file == NULL)
		abort();
	
	stList * vertices = stList_construct();
	char line[LINE_MAX];
	while(fgets(line, LINE_MAX, file)) 
		parseLine(line, vertices);	

	close(file);
	return vertices;
}

// Writing
static void exportNet(FILE * file, stList * net) {
	int val;
	if (stList_length(net) == 0)
		return;

	fprintf(file, "%i", stIntTuple_get(stList_pop(net),0));
	while(stList_length(net) > 0)
		fprintf(file, "\t%i", stIntTuple_get(stList_pop(net),0));
	
	fprintf(file, "\n");
}

static void exportNets(char* filename, stList* nets) {
	FILE * file = fopen(filename, "w");
	if (file == NULL)
		abort();
	
	while(stList_length(nets) > 0)
		exportNet(file, stList_pop(nets));
	
	close(file);
}

// Main function
int main(int argc, char** argv) {
	if (argc!=3)
		abort();
	stList* vertices = parseFile(argv[1]);
	stList* nets = computeThreeEdgeConnectedComponents(vertices);
	exportNets(argv[2], nets);
	return 0;
}
