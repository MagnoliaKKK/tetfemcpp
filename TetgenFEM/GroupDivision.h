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
	const double initx, inity, initz; //初始化后不可更改
	int index;  // Add an index field to help identify vertices
	double vertexMass; // mass of vertices

	//Vertex(double x, double y, double z, int index) : x(x), y(y), z(z), index(index) {}
	Vertex(double x, double y, double z, int index)
		: initx(x), inity(y), initz(z), // 首先初始化const成员
		x(x), y(y), z(z), // 然后是可变成员
		index(index), vertexMass(1) // 其他成员可以直接赋值
	{}
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
	double massTetra;
	double volumeTetra;

	Tetrahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4) {
		vertices[0] = v1;
		vertices[1] = v2;
		vertices[2] = v3;
		vertices[3] = v4;
	}
	Eigen::MatrixXd createElementK(double E, double nu, const Eigen::Vector3d& groupCenterOfMass);
	double calMassTetra(double den);

};


class Group {
public:
	std::vector<Tetrahedron*> tetrahedra;
	std::unordered_map<int, Vertex*> verticesMap;
	Eigen::Vector3d centerofMass;
	double groupMass;//每组的质量
	Eigen::MatrixXd massMatrix;//group mass matrix
	Eigen::MatrixXd massDistribution;
	Eigen::MatrixXd groupK;//group stiffness matrix
	Eigen::VectorXd primeVec;
	Eigen::VectorXd groupVelocity;
	Eigen::VectorXd groupExf;
	Eigen::MatrixXd rotationMatrix;
	Eigen::VectorXd gravity;
	Eigen::MatrixXd dampingMatrix;

	void addTetrahedron(Tetrahedron* tet);
	std::vector<Vertex*> getUniqueVertices();
	void calCenterofMass();
	void calMassGroup();
	Eigen::MatrixXd calMassMatrix(double den);
	void setVertexMassesFromMassMatrix();
	void calMassDistributionMatrix();
	void calGroupK(double E, double nu);
	void calPrimeVec();
	void calDampingMatrix();
};

class Object {
public:
	Group groups[3];

	Group& getGroup(int index);
};

void findBoundaryEdges(tetgenio& out);
void divideIntoGroups(tetgenio& out, Object& object, int numGroups);