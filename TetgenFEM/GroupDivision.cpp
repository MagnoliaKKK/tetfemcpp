#include "GroupDivision.h"
double timeStep = 0.01;
double alpha = 0.1;

void Group::addTetrahedron(Tetrahedron* tet) {
	tetrahedra.push_back(tet);
	for (int i = 0; i < 4; ++i) {
		verticesMap[tet->vertices[i]->index] = tet->vertices[i];
	}
}
std::vector<Vertex*> Group::getUniqueVertices() {
	std::vector<Vertex*> uniqueVertices;
	for (auto& pair : verticesMap) {
		uniqueVertices.push_back(pair.second);
	}
	return uniqueVertices;
}

void Group::calCenterofMass() {
	double totalMass = 0.0;
	Eigen::Vector3d weightedSum(0.0, 0.0, 0.0);

	for (const auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		totalMass += vertex->vertexMass;
		weightedSum += vertex->vertexMass * Eigen::Vector3d(vertex->x, vertex->y, vertex->z);
	}

	if (totalMass > 0.0) {
		centerofMass = weightedSum / totalMass;
	}
	else {
		// Handle the case where totalMass is zero to avoid division by zero.
		// You can set centerofMass to a default value or handle it according to your requirements.
		centerofMass = Eigen::Vector3d(0.0, 0.0, 0.0);
	}
}


//
Eigen::MatrixXd Tetrahedron::createElementK(double E, double nu, const Eigen::Vector3d& groupCenterOfMass) {

	// Initialize the B matrix
	Eigen::MatrixXd B(6, 12);
	B.setZero();

	double N_x[4], N_y[4], N_z[4];

	double p1x = vertices[0]->x;
	double p1y = vertices[0]->y;
	double p1z = vertices[0]->z;
	double p2x = vertices[1]->x;
	double p2y = vertices[1]->y;
	double p2z = vertices[1]->z;
	double p3x = vertices[2]->x;
	double p3y = vertices[2]->y;
	double p3z = vertices[2]->z;
	double p4x = vertices[3]->x;
	double p4y = vertices[3]->y;
	double p4z = vertices[3]->z;

	N_x[0] = (-p3y * p4z + p3z * p4y + p2y * p4z - p2y * p3z - p2z * p4y + p2z * p3y);
	N_y[0] = (p3x * p4z - p3z * p4x - p2x * p4z + p2x * p3z + p2z * p4x - p2z * p3x);
	N_z[0] = (-p3x * p4y + p3y * p4x + p2x * p4y - p2x * p3y - p2y * p4x + p2y * p3x);

	N_x[1] = (p3y * p4z - p3z * p4y - p1y * p4z + p1y * p3z + p1z * p4y - p1z * p3y);
	N_y[1] = (-p3x * p4z + p3z * p4x + p1x * p4z - p1x * p3z - p1z * p4x + p1z * p3x);
	N_z[1] = (p3x * p4y - p3y * p4x - p1x * p4y + p1x * p3y + p1y * p4x - p1y * p3x);

	N_x[2] = (-p2y * p4z + p2z * p4y + p1y * p4z - p1y * p2z - p1z * p4y + p1z * p2y);
	N_y[2] = (p2x * p4z - p2z * p4x - p1x * p4z + p1x * p2z + p1z * p4x - p1z * p2x);
	N_z[2] = (-p2x * p4y + p2y * p4x + p1x * p4y - p1x * p2y - p1y * p4x + p1y * p2x);

	N_x[3] = (p2y * p3z - p2z * p3y - p1y * p3z + p1y * p2z + p1z * p3y - p1z * p2y);
	N_y[3] = (-p2x * p3z + p2z * p3x + p1x * p3z - p1x * p2z - p1z * p3x + p1z * p2x);
	N_z[3] = (p2x * p3y - p2y * p3x - p1x * p3y + p1x * p2y + p1y * p3x - p1y * p2x);

	for (unsigned int i = 0; i < 4; i++) {
		B(0, 3 * i) = N_x[i];
		B(1, 3 * i + 1) = N_y[i];
		B(2, 3 * i + 2) = N_z[i];
		B(3, 3 * i) = N_y[i]; 
		B(3, 3 * i + 1) = N_x[i];
		B(4, 3 * i) = N_z[i]; 
		B(4, 3 * i + 2) = N_x[i];
		B(5, 3 * i + 1) = N_z[i]; 
		B(5, 3 * i + 2) = N_y[i];
	}
	// Calculate the material's elasticity matrix D
	Eigen::MatrixXd D(6, 6);
	double factor = E / ((1 + nu) * (1 - 2 * nu));
	D <<
		(1 - nu) * factor, nu* factor, nu* factor, 0, 0, 0,
		nu* factor, (1 - nu)* factor, nu* factor, 0, 0, 0,
		nu* factor, nu* factor, (1 - nu)* factor, 0, 0, 0,
		0, 0, 0, ((1 - 2 * nu) / 2)* factor, 0, 0,
		0, 0, 0, 0, ((1 - 2 * nu) / 2)* factor, 0,
		0, 0, 0, 0, 0, ((1 - 2 * nu) / 2)* factor;

	// Calculate the element stiffness matrix K
	Eigen::MatrixXd K = (B / (6 * volumeTetra)).transpose() * D * (B/(6 * volumeTetra)) * volumeTetra;
	return K;
}
void Group::calGroupK(double E, double nu) {
	// Initialize groupK to the right size. Assuming 3 degrees of freedom per vertex
	int dof = verticesMap.size() * 3;
	groupK.resize(dof, dof);
	groupK.setZero();

	// Iterate over each tetrahedron to assemble the global stiffness matrix
	for (auto& tetra : tetrahedra) {
		// Get the local stiffness matrix for the current tetrahedron
		Eigen::MatrixXd localK = tetra->createElementK(E, nu, centerofMass);

		// Determine where to add the local stiffness matrix in the global stiffness matrix
		for (int i = 0; i < 4; ++i) { // Each tetrahedron has 4 vertices
			Vertex* vertex = tetra->vertices[i];
			int globalIndex = vertex->index * 3; // Assuming the index is set up correctly in verticesMap

			for (int j = 0; j < 4; ++j) {
				Vertex* otherVertex = tetra->vertices[j];
				int otherGlobalIndex = otherVertex->index * 3;

				// Add the 3x3 submatrix of localK to the correct place in groupK
				groupK.block<3, 3>(globalIndex, otherGlobalIndex) += localK.block<3, 3>(i * 3, j * 3);
			}
		}
	}
}


void Group::calMassGroup() {
	groupMass = 0.0; // 初始化组的质量为0
	for (auto& tet : tetrahedra) { // 遍历每一个四面体
		groupMass += tet->massTetra; // 累加每一个四面体的质量到组的质量
	}
}

Eigen::MatrixXd Group::calMassMatrix(double den) {
	int N = verticesMap.size();  // Number of unique vertices
	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(3 * N, 3 * N);  // Initialize mass matrix

	for (auto& tet : tetrahedra) {
		double tetMass = tet->calMassTetra(den);
		double vertexMass = tetMass / 4.0;  // Assume uniform distribution of mass among vertices

		for (int i = 0; i < 4; i++) {
			int idx = tet->vertices[i]->index;
			M(3 * idx, 3 * idx) += vertexMass;     // x-direction
			M(3 * idx + 1, 3 * idx + 1) += vertexMass; // y-direction
			M(3 * idx + 2, 3 * idx + 2) += vertexMass; // z-direction
		}
	}
	massMatrix = M;
	return M;
}

void Group::calDampingMatrix() {
	dampingMatrix.setZero();
	dampingMatrix = alpha * massMatrix;
}

void Group::calMassDistributionMatrix() {

	massDistribution = Eigen::MatrixXd::Zero(3 * verticesMap.size(), 3 * verticesMap.size());
	Eigen::MatrixXd SUMsub_M_Matrix = Eigen::MatrixXd::Zero(3, 3 * verticesMap.size());
	for (const auto& vertexEntry : verticesMap) {
		Vertex* vertex = vertexEntry.second;
		int vertexIndex = vertex->index;

		// Create a 3x3 identity matrix scaled by the vertex's mass
		SUMsub_M_Matrix.block(0, 3 * vertexIndex, 3, 3) = (vertex->vertexMass / groupMass) * Eigen::Matrix3d::Identity();
	}
	
	for (const auto& vertexEntry : verticesMap) {
		Vertex* vertex = vertexEntry.second;
		int vertexIndex = vertex->index;

		massDistribution.block(3 * vertexIndex, 0, 3, 3 * vertexIndex) = SUMsub_M_Matrix;
	}
	// Reset SUM_M to a sparse view
	//massDistributionSparse.setZero();
	//massDistributionSparse = massDistribution.sparseView();

	
}
void Group::setVertexMassesFromMassMatrix() {
	int N = verticesMap.size();  // Number of unique vertices

	for (int i = 0; i < N; i++) {
		double mass = massMatrix(3 * i, 3 * i);
		// Assuming that the vertices are stored sequentially in the verticesMap
		// And the index in the massMatrix corresponds to the index in the Vertex object
		verticesMap[i]->vertexMass = mass;
	}
}

void Group::calPrimeVec() {
	Eigen::VectorXd gravity(3 * verticesMap.size());
	Eigen::MatrixXd inverseTerm = (massMatrix + dampingMatrix * timeStep).inverse();
	gravity.setZero();
	for (int i = 1; i < 3 * verticesMap.size(); i += 3) {
		gravity(i) = 9.8;
	}
	primeVec = Eigen::VectorXd::Zero(3 * verticesMap.size());
	for (int pi = 0; pi < verticesMap.size(); pi++)
	{
		Eigen::Vector3d currentPosition = Eigen::Vector3d(verticesMap[pi]->x, verticesMap[pi]->y, verticesMap[pi]->z);
		groupVelocity.block(3 * pi, 0, 3, 1) = groupVelocity.block(3 * pi, 0, 3, 1)  + gravity.block(3 * pi, 0, 3, 1) * timeStep;
		Eigen::Vector3d velocityUpdate = inverseTerm * massMatrix * groupVelocity.block(3 * pi, 0, 3, 1) * timeStep;
		Eigen::Vector3d newPosition = currentPosition + velocityUpdate;
		// 更新primeVec为当前位置加上速度更新
		
		verticesMap[pi]->x = newPosition.x();
		verticesMap[pi]->y = newPosition.y();
		verticesMap[pi]->z = newPosition.z();
		primeVec.block(3 * pi, 0, 3, 1) = newPosition;
	}
	
	
}

Group& Object::getGroup(int index) {
	return groups[index];
}

std::unordered_set<std::string> boundaryEdgesSet;  // Set to store boundary edges

void findBoundaryEdges(tetgenio& out) {
	int indexOffset = out.firstnumber;  // Get the index offset (0 or 1)
	for (int i = 0; i < out.numberoftrifaces; ++i) {
		for (int j = 0; j < 3; ++j) {
			int vertexIndex1 = out.trifacelist[i * 3 + j] - indexOffset;
			int vertexIndex2 = out.trifacelist[i * 3 + ((j + 1) % 3)] - indexOffset;
			std::string edgeKey = vertexIndex1 < vertexIndex2 ?
				std::to_string(vertexIndex1) + "-" + std::to_string(vertexIndex2) :
				std::to_string(vertexIndex2) + "-" + std::to_string(vertexIndex1);
			boundaryEdgesSet.insert(edgeKey);
		}
	}
}

void divideIntoGroups(tetgenio& out, Object& object, int numGroups) {

	findBoundaryEdges(out);  // Populate the boundaryEdgesSet

	double minX = out.pointlist[0];
	double maxX = out.pointlist[0];

	// Find min and max X coordinates
	for (int i = 0; i < out.numberofpoints; ++i) {
		double x = out.pointlist[i * 3];
		if (x < minX) minX = x;
		if (x > maxX) maxX = x;
	}

	double range = maxX - minX;
	double groupRange = range / numGroups;

	// Create vertices
	std::vector<Vertex*> vertices;
	for (int i = 0; i < out.numberofpoints; ++i) {
		double x = out.pointlist[i * 3];
		double y = out.pointlist[i * 3 + 1];
		double z = out.pointlist[i * 3 + 2];
		vertices.push_back(new Vertex(x, y, z, i));
	}

	// Create tetrahedra and assign to groups
	for (int i = 0; i < out.numberoftetrahedra; ++i) {
		Vertex* v1 = vertices[out.tetrahedronlist[i * 4] - 1];
		Vertex* v2 = vertices[out.tetrahedronlist[i * 4 + 1] - 1];
		Vertex* v3 = vertices[out.tetrahedronlist[i * 4 + 2] - 1];
		Vertex* v4 = vertices[out.tetrahedronlist[i * 4 + 3] - 1];

		Tetrahedron* tet = new Tetrahedron(v1, v2, v3, v4);


		// Determine group based on average X coordinate
		double avgX = (v1->x + v2->x + v3->x + v4->x) / 4;
		int groupIndex = (avgX - minX) / groupRange;
		if (groupIndex >= numGroups) groupIndex = numGroups - 1;  // Update here

		object.getGroup(groupIndex).addTetrahedron(tet);

		// Set up edges for each tetrahedron
		static int edgeIndices[6][2] = { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3} };
		for (int j = 0; j < 6; ++j) {
			Vertex* vertex1 = tet->vertices[edgeIndices[j][0]];
			Vertex* vertex2 = tet->vertices[edgeIndices[j][1]];
			Edge* edge = new Edge(vertex1, vertex2);

			std::string edgeKey = vertex1->index < vertex2->index ?
				std::to_string(vertex1->index) + "-" + std::to_string(vertex2->index) :
				std::to_string(vertex2->index) + "-" + std::to_string(vertex1->index);

			edge->isBoundary = boundaryEdgesSet.count(edgeKey) > 0;
			tet->edges[j] = edge;
		}
	}
}

double Tetrahedron::calMassTetra(double den) {
	
	//double volume;
	Eigen::Vector3d AB(vertices[1]->x - vertices[0]->x, vertices[1]->y - vertices[0]->y, vertices[1]->z - vertices[0]->z);
	Eigen::Vector3d AC(vertices[2]->x - vertices[0]->x, vertices[2]->y - vertices[0]->y, vertices[2]->z - vertices[0]->z);
	Eigen::Vector3d AD(vertices[3]->x - vertices[0]->x, vertices[3]->y - vertices[0]->y, vertices[3]->z - vertices[0]->z);

	// Calculate volume using the formula
	volumeTetra = (AB.cross(AC)).dot(AD) / 6.0;
	volumeTetra = std::abs(volumeTetra);
	massTetra = volumeTetra * den;
	return massTetra;

	
}