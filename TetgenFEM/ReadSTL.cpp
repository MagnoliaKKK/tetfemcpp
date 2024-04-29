#include "ReadSTL.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


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



// Assuming REAL is defined as float or double somewhere, and tetgenio structures are defined appropriately
void readOBJ(const std::string& filename, tetgenio& in) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::vector<std::vector<float>> vertices;
    std::vector<std::vector<int>> faces;

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::string type;
        ss >> type;

        if (type == "v") {
            std::vector<float> vertex(3);
            ss >> vertex[0] >> vertex[1] >> vertex[2];
            vertices.push_back(vertex);
        }
        else if (type == "f") {
            std::vector<int> face;
            std::string vertex_index;
            while (ss >> vertex_index) {
                size_t slash_pos = vertex_index.find('/');
                int index = std::stoi(vertex_index.substr(0, slash_pos)) - 1;  // OBJ indices are 1-based
                face.push_back(index);
            }
            faces.push_back(face);
        }
    }

    file.close();

    // Allocate tetgenio data
    in.numberofpoints = vertices.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    for (size_t i = 0; i < vertices.size(); ++i) {
        in.pointlist[i * 3] = vertices[i][0];
        in.pointlist[i * 3 + 1] = vertices[i][1];
        in.pointlist[i * 3 + 2] = vertices[i][2];
    }

    in.numberoffacets = faces.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    for (size_t i = 0; i < faces.size(); ++i) {
        tetgenio::facet& f = in.facetlist[i];
        f.numberofpolygons = 1;
        f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
        f.numberofholes = 0;
        f.holelist = NULL;
        tetgenio::polygon& p = f.polygonlist[0];
        p.numberofvertices = faces[i].size();
        p.vertexlist = new int[p.numberofvertices];

        for (size_t j = 0; j < faces[i].size(); ++j) {
            p.vertexlist[j] = faces[i][j] + 1;  // Adjust for 1-based indexing used in TetGen
        }

        in.facetmarkerlist[i] = 0;  // No facet markers
    }
}