#pragma once
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "tetgen.h"
#include <cstring> 
#include <fstream>
#include <iostream>
#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry>
class Vertex {
public:
	double x, y, z;
	int index;  // Add an index field to help identify vertices

	Vertex(double x, double y, double z, int index) : x(x), y(y), z(z), index(index) {}
};

class Edge {
public:
	Vertex* vertices[2];
	bool isBoundary;

	Edge(Vertex* v1, Vertex* v2) : isBoundary(false) {
		vertices[0] = v1;
		vertices[1] = v2;
	}
};

class Tetrahedron {
public:
	Vertex* vertices[4];
	Edge* edges[6];  // Each tetrahedron has six edges

	Tetrahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4) {
		vertices[0] = v1;
		vertices[1] = v2;
		vertices[2] = v3;
		vertices[3] = v4;
	}
	Eigen::MatrixXd createElementK(double E, double nu);
};


class Group {
public:
	std::vector<Tetrahedron*> tetrahedra;
	std::unordered_map<int, Vertex*> verticesMap;

	void addTetrahedron(Tetrahedron* tet);
	std::vector<Vertex*> getUniqueVertices();
};

class Object {
public:
	Group groups[3];

	Group& getGroup(int index);
};

void findBoundaryEdges(tetgenio& out);
void divideIntoGroups(tetgenio& out, Object& object, int numGroups);