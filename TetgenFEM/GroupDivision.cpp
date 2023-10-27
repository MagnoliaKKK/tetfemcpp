#include "GroupDivision.h"


void Group::addTetrahedron(Tetrahedron* tet) {
	tetrahedra.push_back(tet);
	for (int i = 0; i < 4; ++i) {
		verticesMap[tet->vertices[i]->index] = tet->vertices[i];
	}
}

Eigen::MatrixXd Tetrahedron::createElementK(double E, double nu) {
	Eigen::MatrixXd J(3, 3);
	Eigen::MatrixXd K(12, 12);
	J(0, 0) = vertices[0]->x - vertices[3]->x;
	J(0, 1) = vertices[1]->x - vertices[3]->x;
	J(0, 2) = vertices[2]->x - vertices[3]->x;

	J(1, 0) = vertices[0]->y - vertices[3]->y;
	J(1, 1) = vertices[1]->y - vertices[3]->y;
	J(1, 2) = vertices[2]->y - vertices[3]->y;

	J(2, 0) = vertices[0]->z - vertices[3]->z;
	J(2, 1) = vertices[1]->z - vertices[3]->z;
	J(2, 2) = vertices[2]->z - vertices[3]->z;
	Eigen::MatrixXd J_inv = J.inverse();

	Eigen::MatrixXd B(6, 12);
	B.setZero();
	// Node 0
	Eigen::Vector3d dN0_dxieta(1, 0, 0);
	Eigen::Vector3d dN0_dx = J_inv * dN0_dxieta;
	B(0, 0) = dN0_dx(0);
	B(1, 1) = dN0_dx(1);
	B(2, 2) = dN0_dx(2);
	B(3, 0) = dN0_dx(1);
	B(3, 1) = dN0_dx(0);
	B(4, 1) = dN0_dx(2);
	B(4, 2) = dN0_dx(1);
	B(5, 0) = dN0_dx(2);
	B(5, 2) = dN0_dx(0);

	// Node 1
	Eigen::Vector3d dN1_dxieta(0, 1, 0);
	Eigen::Vector3d dN1_dx = J_inv * dN1_dxieta;
	B(0, 3) = dN1_dx(0);
	B(1, 4) = dN1_dx(1);
	B(2, 5) = dN1_dx(2);
	B(3, 3) = dN1_dx(1);
	B(3, 4) = dN1_dx(0);
	B(4, 4) = dN1_dx(2);
	B(4, 5) = dN1_dx(1);
	B(5, 3) = dN1_dx(2);
	B(5, 5) = dN1_dx(0);

	// Node 2
	Eigen::Vector3d dN2_dxieta(0, 0, 1);
	Eigen::Vector3d dN2_dx = J_inv * dN2_dxieta;
	B(0, 6) = dN2_dx(0);
	B(1, 7) = dN2_dx(1);
	B(2, 8) = dN2_dx(2);
	B(3, 6) = dN2_dx(1);
	B(3, 7) = dN2_dx(0);
	B(4, 7) = dN2_dx(2);
	B(4, 8) = dN2_dx(1);
	B(5, 6) = dN2_dx(2);
	B(5, 8) = dN2_dx(0);

	// Node 3
	Eigen::Vector3d dN3_dxieta(-1, -1, -1);
	Eigen::Vector3d dN3_dx = J_inv * dN3_dxieta;
	B(0, 9) = dN3_dx(0);
	B(1, 10) = dN3_dx(1);
	B(2, 11) = dN3_dx(2);
	B(3, 9) = dN3_dx(1);
	B(3, 10) = dN3_dx(0);
	B(4, 10) = dN3_dx(2);
	B(4, 11) = dN3_dx(1);
	B(5, 9) = dN3_dx(2);
	B(5, 11) = dN3_dx(0);
	Eigen::MatrixXd D(6, 6);
	double factor = E / ((1 + nu) * (1 - 2 * nu));
	D <<
		(1 - nu) * factor, nu* factor, nu* factor, 0, 0, 0,
		nu* factor, (1 - nu)* factor, nu* factor, 0, 0, 0,
		nu* factor, nu* factor, (1 - nu)* factor, 0, 0, 0,
		0, 0, 0, ((1 - 2 * nu) / 2)* factor, 0, 0,
		0, 0, 0, 0, ((1 - 2 * nu) / 2)* factor, 0,
		0, 0, 0, 0, 0, ((1 - 2 * nu) / 2)* factor;
	K = B.transpose() * D * B;
	return K;
}

std::vector<Vertex*> Group::getUniqueVertices() {
	std::vector<Vertex*> uniqueVertices;
	for (auto& pair : verticesMap) {
		uniqueVertices.push_back(pair.second);
	}
	return uniqueVertices;
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
