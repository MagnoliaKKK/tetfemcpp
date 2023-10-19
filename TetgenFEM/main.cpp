#include <iostream>
#include <vector>
#include <cstring>  // for std::strcpy
#include "tetgen.h"  // Include the TetGen header file
#include <fstream>
// Function to read STL file
void readSTL(const std::string& filename, tetgenio& in) {
	std::ifstream file(filename, std::ios::binary);
	if (!file) {
		std::cerr << "Error: Could not open file: " << filename << std::endl;
		return;
	}

	char header[80];
	file.read(header, 80);  // Read header

	uint32_t numTriangles;
	file.read(reinterpret_cast<char*>(&numTriangles), 4);  // Read number of triangles

	in.numberofpoints = numTriangles * 3;  // Each triangle has 3 vertices
	in.pointlist = new REAL[in.numberofpoints * 3];
	in.numberoffacets = numTriangles;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for (uint32_t i = 0; i < numTriangles; ++i) {
		float normal[3];
		file.read(reinterpret_cast<char*>(normal), 12);  // Read normal vector

		tetgenio::facet& f = in.facetlist[i];
		f.numberofpolygons = 1;
		f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
		f.numberofholes = 0;
		f.holelist = NULL;
		tetgenio::polygon& p = f.polygonlist[0];
		p.numberofvertices = 3;
		p.vertexlist = new int[p.numberofvertices];

		for (int j = 0; j < 3; ++j) {  // Read vertices
			float vertex[3];
			file.read(reinterpret_cast<char*>(vertex), 12);
			int index = i * 3 + j;
			in.pointlist[index * 3] = vertex[0];
			in.pointlist[index * 3 + 1] = vertex[1];
			in.pointlist[index * 3 + 2] = vertex[2];
			p.vertexlist[j] = index + 1;  // Indices in TetGen start from 1
		}

		uint16_t attributeByteCount;
		file.read(reinterpret_cast<char*>(&attributeByteCount), 2);  // Read attribute byte count

		in.facetmarkerlist[i] = 0;  // No facet markers
	}

	file.close();
}


int main() {

	tetgenio in, out;
	in.firstnumber = 1;  // All indices start from 1
	readSTL("C:/Users/XYX/Documents/VSProgramming/2.stl", in);

	// Configure TetGen behavior
	tetgenbehavior behavior;
	char args[] = "pq1.414a0.1";
	behavior.parse_commandline(args);

	// Call TetGen to tetrahedralize the geometry
	tetrahedralize(&behavior, &in, &out);

	// ... rest of your code ...

	return 0;
}
