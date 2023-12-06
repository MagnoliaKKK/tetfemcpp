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
#include "Eigen/Sparse"
#include "GMRES.h"
#include "Eigen/IterativeLinearSolvers"



class Vertex {
public:
	double x, y, z;
	const double initx, inity, initz; //初始化后不可更改
	int index;  // global index
	int localIndex; // 组内的本地索引
	double vertexMass; // mass of vertices
	double velx, vely, velz;//速度的三个分量
	bool isFixed;

	//Vertex(double x, double y, double z, int index) : x(x), y(y), z(z), index(index) {}
	Vertex(double x, double y, double z, int index)
		: initx(x), inity(y), initz(z), // 首先初始化const成员
		x(x), y(y), z(z), // 然后是可变成员
		index(index), vertexMass(1), // 其他成员可以直接赋值
		velx(0), vely(0), velz(0), // 初始化速度分量为0
		isFixed(false) // 默认不是固定点
	{}
	void setFixedIfBelowThreshold() {
		if (initx < -0.619) {
			isFixed = true;
		}
	}
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
	Eigen::MatrixXd elementK;

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
	Eigen::Vector3d initCOM;//initial center of mass
	Eigen::VectorXd initLocalPos;//initial position - center of mass
	Eigen::MatrixXd FEMLHS;
	Eigen::VectorXd FEMRHS;
	Eigen::VectorXd Fbind;
	Eigen::VectorXd deltaX;
	Eigen::SparseMatrix<double> rotationSparse;
	Eigen::SparseMatrix<double> rotationTransSparse;//R^t sparse
	Eigen::SparseMatrix<double> kSparse;
	Eigen::SparseMatrix<double> massSparse;
	Eigen::SparseMatrix<double> dampingSparse;
	Eigen::SparseMatrix<double> massDistributionSparse;
	Eigen::SparseMatrix<double> massDampingSparseInv; //(M+C').inv sparse
	Eigen::VectorXd currentPosition;//计算bindf用的位置信息，不用做位置更新


	
	void addTetrahedron(Tetrahedron* tet);
	std::vector<Vertex*> getUniqueVertices();
	void calCenterofMass();
	void calMassGroup();
	Eigen::MatrixXd calMassMatrix(double den);
	void setVertexMassesFromMassMatrix();
	void calMassDistributionMatrix();
	void calGroupK(double E, double nu);
	void calPrimeVec(int w);
	void calDampingMatrix();
	void calInitCOM();
	void calRotationMatrix();
	void calLocalPos();//calculate initial local position
	Eigen::Vector3d axlAPD(Eigen::Matrix3d a);
	Eigen::Vector3d clamp2(Eigen::Vector3d x, double y, double z);
	Eigen::Quaterniond Exp2(Eigen::Vector3d a);
	void calLHS();
	void calRHS();
	void calDeltaX();
	void calculateCurrentPositions();
	void calFbind(const std::vector<Vertex*>& commonVerticesGroup1,
		const std::vector<Vertex*>& commonVerticesGroup2,
		const Eigen::VectorXd& currentPositionGroup1,
		const Eigen::VectorXd& currentPositionGroup2,
		double k);
	void updatePosition();
	void updateVelocity();
	void initialize();
	//void updateVertexPositions();


	Group()
		: centerofMass(Eigen::Vector3d::Zero()),  // Initialize Eigen vector
		groupMass(0.0),                         // Initialize double
		verticesMap()
		
	{
		// Additional initialization logic, if needed
	}
};


class Object {
public:
	Group groups[1]; // change this
	std::pair<std::vector<Vertex*>, std::vector<Vertex*>> commonPoints;
	std::pair<std::vector<Vertex*>, std::vector<Vertex*>> commonPoints1;

	Group& getGroup(int index);
	std::pair<std::vector<Vertex*>, std::vector<Vertex*>> findCommonVertices(const Group& group1, const Group& group2);// find common vertex
	void assignLocalIndicesToAllGroups(); // local Index
	void updateIndices();
	void generateUniqueVertices();//generate unique vertices
	void PBDLOOP(int looptime);
};

void findBoundaryEdges(tetgenio& out);
void divideIntoGroups(tetgenio& out, Object& object, int numGroups);
