#include "Group.h"

void Group::initialize() {
	groupVelocity = Eigen::VectorXf::Zero(3 * verticesMap.size());
	groupVelocityFEM = Eigen::VectorXf::Zero(3 * verticesMap.size());
	Fbind = Eigen::VectorXf::Zero(3 * verticesMap.size());
	prevFbind = Eigen::VectorXf::Zero(3 * verticesMap.size());
	currentPosition = Eigen::VectorXf::Zero(3 * verticesMap.size());
	currentPositionFEM = Eigen::VectorXf::Zero(3 * verticesMap.size());
	gravity = Eigen::VectorXf::Zero(3 * verticesMap.size());

	verticesVector.reserve(verticesMap.size());
	for (const auto& pair : verticesMap) {
		verticesVector.push_back(pair.second);
	}
}
void Group::addTetrahedron(Tetrahedron* tet) {
	tetrahedra.push_back(tet);
	//for (int i = 0; i < 4; ++i) {
	//	verticesMap[tet->vertices[i]->index] = tet->vertices[i]; //�Y���l�ʑ̓I��?�C�c�l�ʑ̓I?�_����verticesMap
	//}
}
std::vector<Vertex*> Group::getUniqueVertices() { //?��?�����v�I�C�������chashmap??��vertexGroup
	std::vector<Vertex*> uniqueVertices;
	for (auto& pair : verticesMap) {
		uniqueVertices.push_back(pair.second);
	}
	return uniqueVertices;
}

void Group::calCenterofMass() {
	float totalMass = 0.0;
	Eigen::Vector3f weightedSum(0.0, 0.0, 0.0);

	for (const auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		totalMass += vertex->vertexMass;
		weightedSum += vertex->vertexMass * Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
	}

	if (totalMass > 0.0) {
		centerofMass = weightedSum / totalMass;
	}
	else {
		// Handle the case where totalMass is zero to avoid division by zero.
		// You can set centerofMass to a default value or handle it according to your requirements.
		centerofMass = Eigen::Vector3f(0.0, 0.0, 0.0);
	}
}
void Group::calLocalPos() {
	// ?��initLocalPos�L��?�I��?����?���L�I�Ǖ��ʒu
	initLocalPos.resize(3 * verticesMap.size());

	for (const auto& vertexPair : verticesMap) {
		const Vertex* vertex = vertexPair.second;
		Eigen::Vector3f initial_position(vertex->initx, vertex->inity, vertex->initz);
		// ?�Z���n�ʒu�^���n�d�S�I��?
		Eigen::Vector3f local_position = initial_position - initCOM;

		// ���Ǖ��ʒu��?��initLocalPos���C����index���v����3��??��?�_�L3����??
		initLocalPos.segment<3>(vertex->localIndex * 3) = local_position;
	}
}
Eigen::Vector3f Group::axlAPD(Eigen::Matrix3f a) {
	Eigen::Vector3f g = Eigen::Vector3f::Zero();
	g[0] = a(1, 2) - a(2, 1);
	g[1] = a(2, 0) - a(0, 2);
	g[2] = a(0, 1) - a(1, 0);
	return g;
}
Eigen::Vector3f Group::clamp2(Eigen::Vector3f x, float y, float z) {
	if (x.norm() < y) {
		return (y / x.norm()) * x;
	}
	else if (x.norm() > z) {
		return (z / x.norm()) * x;
	}
	else {
		return x;
	}
}
Eigen::Quaternionf Group::Exp2(Eigen::Vector3f a) {
	float s = sin((a * 0.5).norm());
	float x = s * a.x() / a.norm();
	float y = s * a.y() / a.norm();
	float z = s * a.z() / a.norm();
	Eigen::Quaternionf qq = Eigen::Quaternionf(cos((a * 0.5).norm()), x, y, z);
	return  qq;
}
void Group::calRotationMatrix() {
	Eigen::MatrixXf Apq = Eigen::MatrixXf::Zero(3, 3);
	Eigen::Matrix3f tempA = Eigen::Matrix3f::Zero();
	Eigen::Vector3f center_grid = Eigen::Vector3f::Zero();

	Eigen::Vector3f tempCenterGrid = Eigen::Vector3f::Zero();
#pragma omp parallel num_threads(4)
	{
		Eigen::Vector3f localCenterGrid = Eigen::Vector3f::Zero();
#pragma omp for
		for (int i = 0; i < static_cast<int>(verticesMap.size()); ++i) {
			auto it = verticesMap.begin();
			std::advance(it, i);
			const std::pair<const int, Vertex*>& vertexEntry = *it;

			localCenterGrid.noalias() += massDistribution(0, 3 * vertexEntry.second->localIndex) * primeVec.segment<3>(3 * vertexEntry.second->localIndex);
		}
#pragma omp critical
		{
			tempCenterGrid += localCenterGrid;
		}
	}
	center_grid = tempCenterGrid;


	// ?�ZApq��?
	Eigen::MatrixXf tempApq = Eigen::MatrixXf::Zero(3, 3);
#pragma omp parallel num_threads(4)
	{
		Eigen::MatrixXf localApq = Eigen::MatrixXf::Zero(3, 3);
#pragma omp for
		for (int i = 0; i < static_cast<int>(verticesMap.size()); ++i) {
			auto it = verticesMap.begin();
			std::advance(it, i);
			const std::pair<const int, Vertex*>& vertexEntry = *it;

			tempA.noalias() = (primeVec.segment<3>(3 * vertexEntry.second->localIndex) - center_grid) * initLocalPos.segment<3>(3 * vertexEntry.second->localIndex).transpose();
			localApq.noalias() += massMatrix(3 * vertexEntry.second->localIndex, 3 * vertexEntry.second->localIndex) * tempA;
		}
#pragma omp critical
		{
			tempApq += localApq;
		}
	}
	Apq = tempApq;

	// ���n���l�����a��?��?
	Eigen::Vector3f omega = Eigen::Vector3f::Identity();
	Eigen::Quaternionf quaternion(Eigen::Quaternionf::Identity());
	rotate_matrix = Eigen::Matrix3f::Identity();
	Eigen::Vector3f gradR = Eigen::Vector3f::Zero();
	Eigen::Matrix3f HesseR = Eigen::Matrix3f::Zero();
	Eigen::Matrix3f S = Eigen::Matrix3f::Zero();

	// �R��?�Q�ŉ���?
	for (unsigned int ci = 0; ci < 20; ci++) {
		Eigen::Matrix3f R = quaternion.matrix();
		Eigen::Matrix3f S = R.transpose() * Apq;
		Eigen::Vector3f gradR = axlAPD(S);
		Eigen::Matrix3f HesseR = S.trace() * Eigen::Matrix3f::Identity() - (S + S.transpose()) * 0.5;
		Eigen::Vector3f omega = -HesseR.inverse() * gradR;
		float w = omega.norm();
		if (w < 1.0e-9) {
			break;
		}
		omega = clamp2(omega, -PI, PI);
		Eigen::Quaternionf temp2 = Exp2(omega);
		quaternion = quaternion * temp2;
	}

	rotate_matrix = quaternion.matrix();

	// ?����?��?�I3N x 3N�Ŗ{
	rotationMatrix = Eigen::MatrixXf::Zero(3 * verticesMap.size(), 3 * verticesMap.size());

#pragma omp parallel for
	for (int pi = 0; pi < static_cast<int>(verticesMap.size()); pi++) {
		rotationMatrix.block<3, 3>(3 * pi, 3 * pi) = rotate_matrix;
	}
}
//

void Group::calGroupK(float E, float nu) {
	// Initialize groupK to the right size. Assuming 3 degrees of freedom per vertex
	int dof = verticesMap.size() * 3;
	groupK.resize(dof, dof);
	groupK.setZero();

	// Iterate over each tetrahedron to assemble the global stiffness matrix
	for (auto& tetra : tetrahedra) {
		// Get the local stiffness matrix for the current tetrahedron
		Eigen::MatrixXf localK = tetra->createElementK(E, nu, initCOM);

		// Determine where to add the local stiffness matrix in the global stiffness matrix
		for (int i = 0; i < 4; ++i) { // Each tetrahedron has 4 vertices
			Vertex* vertex = tetra->vertices[i];
			int localIndex = vertex->localIndex * 3; // Use local index

			for (int j = 0; j < 4; ++j) {
				Vertex* otherVertex = tetra->vertices[j];
				int otherLocalIndex = otherVertex->localIndex * 3;

				// Add the 3x3 submatrix of localK to the correct place in groupK
				groupK.block<3, 3>(localIndex, otherLocalIndex) += localK.block<3, 3>(i * 3, j * 3);
			}
		}
	}

	kSparse = groupK.sparseView();
}
void Group::calGroupKAni(float E1, float E2, float E3, float nu) {
	// Initialize groupK to the right size. Assuming 3 degrees of freedom per vertex
	int dof = verticesMap.size() * 3;
	groupK.resize(dof, dof);
	groupK.setZero();

	// Iterate over each tetrahedron to assemble the global stiffness matrix
	for (auto& tetra : tetrahedra) {
		// Get the local stiffness matrix for the current tetrahedron
		Eigen::MatrixXf localK = tetra->createElementKAni(E1, E2, E3, nu, initCOM);

		// Determine where to add the local stiffness matrix in the global stiffness matrix
		for (int i = 0; i < 4; ++i) { // Each tetrahedron has 4 vertices
			Vertex* vertex = tetra->vertices[i];
			int localIndex = vertex->localIndex * 3; // Use local index

			for (int j = 0; j < 4; ++j) {
				Vertex* otherVertex = tetra->vertices[j];
				int otherLocalIndex = otherVertex->localIndex * 3;

				// Add the 3x3 submatrix of localK to the correct place in groupK
				groupK.block<3, 3>(localIndex, otherLocalIndex) += localK.block<3, 3>(i * 3, j * 3);
			}
		}
	}

	kSparse = groupK.sparseView();
}
void Group::calGroupKFEM(float E, float nu) {
	// Initialize groupK to the right size. Assuming 3 degrees of freedom per vertex
	int dof = verticesMap.size() * 3;
	groupKFEM.resize(dof, dof);
	groupKFEM.setZero();

	// Iterate over each tetrahedron to assemble the global stiffness matrix
	for (auto& tetra : tetrahedra) {
		// Get the local stiffness matrix for the current tetrahedron
		Eigen::MatrixXf localK = tetra->createElementKFEM(E, nu);

		// Determine where to add the local stiffness matrix in the global stiffness matrix
		for (int i = 0; i < 4; ++i) { // Each tetrahedron has 4 vertices
			Vertex* vertex = tetra->vertices[i];
			int localIndex = vertex->localIndex * 3; // Use local index

			for (int j = 0; j < 4; ++j) {
				Vertex* otherVertex = tetra->vertices[j];
				int otherLocalIndex = otherVertex->localIndex * 3;

				// Add the 3x3 submatrix of localK to the correct place in groupK
				groupKFEM.block<3, 3>(localIndex, otherLocalIndex) += localK.block<3, 3>(i * 3, j * 3);
			}
		}
	}

	kSparseFEM = groupKFEM.sparseView();
}


void Group::calMassGroup() {
	groupMass = 0.0; // ���n��?�I?��?0
	for (auto& tet : tetrahedra) { // ��??�꘢�l�ʑ�
		groupMass += tet->massTetra; // �݉�?�꘢�l�ʑ̓I?�ʓ�?�I?��
	}
}

Eigen::MatrixXf Group::calMassMatrix(float den) {
	int N = verticesMap.size();  // Number of unique vertices
	Eigen::MatrixXf M = Eigen::MatrixXf::Zero(3 * N, 3 * N);  // Initialize mass matrix

	for (auto& tet : tetrahedra) {
		float tetMass = tet->calMassTetra(den);
		float vertexMass = tetMass / 4.0;  // Assume uniform distribution of mass among vertices

		for (int i = 0; i < 4; i++) {
			int localIdx = tet->vertices[i]->localIndex;  // Use local index
			M(3 * localIdx, 3 * localIdx) += vertexMass;     // x-direction
			M(3 * localIdx + 1, 3 * localIdx + 1) += vertexMass; // y-direction
			M(3 * localIdx + 2, 3 * localIdx + 2) += vertexMass; // z-direction
		}
	}
	massMatrix = M;
	return M;
	massSparse = massMatrix.sparseView();  // Commented out as it's not part of the provided code snippet
}


void Group::calDampingMatrix() {
	dampingMatrix.setZero();
	dampingMatrix = dampingConst * massMatrix;
	dampingSparse = dampingMatrix.sparseView();
}

void Group::calMassDistributionMatrix() {
	int N = verticesMap.size();  // Number of unique vertices
	massDistribution = Eigen::MatrixXf::Zero(3 * N, 3 * N);
	Eigen::MatrixXf SUMsub_M_Matrix = Eigen::MatrixXf::Zero(3, 3 * N);

	for (const auto& vertexEntry : verticesMap) {
		Vertex* vertex = vertexEntry.second;
		int localIdx = vertex->localIndex;  // Use local index

		// Create a 3x3 identity matrix scaled by the vertex's mass
		SUMsub_M_Matrix.block(0, 3 * localIdx, 3, 3) = (vertex->vertexMass / groupMass) * Eigen::Matrix3f::Identity();
	}

	for (const auto& vertexEntry : verticesMap) {
		Vertex* vertex = vertexEntry.second;
		int localIdx = vertex->localIndex;  // Use local index

		massDistribution.block(3 * localIdx, 0, 3, 3 * N) = SUMsub_M_Matrix;
	}

	// Reset SUM_M to a sparse view
	massDistributionSparse = massDistribution.sparseView();
}

void Group::setVertexMassesFromMassMatrix() {
	int N = verticesMap.size();  // Number of unique vertices in the hash map

	for (auto& vertexPair : verticesMap) {
		int localIdx = vertexPair.second->localIndex;  // Access local index from the second value of the pair
		float mass = massMatrix(3 * localIdx, 3 * localIdx);
		vertexPair.second->vertexMass = mass;
	}
}



void Group::calInitCOM() {
	float totalMass = 0.0;
	Eigen::Vector3f weightedSum(0.0, 0.0, 0.0);

	for (const auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		totalMass += vertex->vertexMass;
		weightedSum += vertex->vertexMass * Eigen::Vector3f(vertex->initx, vertex->inity, vertex->initz);
	}

	if (totalMass > 0.0) {
		initCOM = weightedSum / totalMass;
	}
	else {
		// Handle the case where totalMass is zero to avoid division by zero.
		// You can set centerofMass to a default value or handle it according to your requirements.
		initCOM = Eigen::Vector3f(0.0, 0.0, 0.0);
	}
}
void Group::calPrimeVec1(int w) {
	// Assuming Gravity and timeStep are defined elsewhere
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());
	gravity = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// Find min and max x-values
	float minX = std::numeric_limits<float>::max();
	float maxX = std::numeric_limits<float>::lowest();
	for (const auto& vertex : verticesVector) {
		if (vertex->x < minX) minX = vertex->x;
		if (vertex->x > maxX) maxX = vertex->x;
	}

	// Calculate the midpoint of the x-axis
	float midX = (minX + maxX) / 2;

	// Apply forces based on position relative to the midpoint
	for (size_t i = 0; i < verticesVector.size(); ++i) {
		float x = verticesVector[i]->x;
		int yIndex = 3 * i + 1; // Index for the y-component in the gravity vector

		if (x <= midX) {
			// Front half, apply +y force
			gravity(yIndex) = Gravity;
		}
		else {
			// Back half, apply -y force
			gravity(yIndex) = -Gravity;
		}
	}

	// Update groupVelocity
	groupVelocity += gravity * timeStep;

	// Compute velocityUpdate
	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;

	// Update primeVec and vertex positions
	for (auto& vertex : verticesVector) {
		int localPi = vertex->localIndex; // Use local index

		// Get the current vertex's velocity update part
		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);

		// Calculate new position
		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;

		// Update primeVec
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
	}
}
//void Group::calPrimeVecT(int w) {
//	// ?��primeVec��?���n����?�u?��?�I�ڐ�
//	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());
//
//	// ���n����?�͌���
//	Eigen::VectorXf twistForce = Eigen::VectorXf::Zero(3 * verticesVector.size());
//
//	// ??���ې���O?�C�@�ʐ��C??����?�_�{����?��
//	if (this->groupIndex == 2) {
//		// ?�u?��?�_��?�͓I�召�C?��?��??��召�C���?���v����??��v?��
//		float forceMagnitude = 10.0f; // ����͓I�召
//
//		// ?�_�Ǖ�������?
//		std::vector<int> indices = { 116, 91, 24, 22 };
//
//		// ??��?�_�{����?��
//		for (int index : indices) {
//			// ?��͓I�����F����?�_��yz�ʏ�I???��?�����{��
//			if (index == 116) {
//				// ?��116�C��?�͌����iy??�����j
//				twistForce(3 * index + 1) = -forceMagnitude;
//			}
//			else if (index == 91) {
//				// ?��91�C��?�͌��E�iz?�������j
//				twistForce(3 * index + 2) = forceMagnitude;
//			}
//			else if (index == 24) {
//				// ?��24�C��?�͌���iy?�������j
//				twistForce(3 * index + 1) = forceMagnitude;
//			}
//			else if (index == 22) {
//				// ?��22�C��?�͌����iz??�����j
//				twistForce(3 * index + 2) = -forceMagnitude;
//			}
//		}
//	}
//
//	// �X�VgroupVelocity
//	groupVelocity += twistForce * timeStep;
//
//	// �g�p������??�ZvelocityUpdate
//	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;
//
//	// �X�VprimeVec�a?�_�ʒu
//	for (auto& vertexPair : verticesVector) {
//		Vertex* vertex = vertexPair;
//		int localPi = vertex->localIndex; // �g�p�Ǖ�����
//
//		// ?�擖�O?�_�I���x�X�V����
//		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);
//
//		// ?�Z�V�I�ʒu
//		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;
//
//		// �X�VprimeVec
//		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
//	}
//}
//void Group::calPrimeVecS(int wKey) {
//	// ?��primeVec��?���n����?�u?��?�I�ڐ�
//	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());
//
//	// ���n���꘢�͌���
//	Eigen::VectorXf appliedForce = Eigen::VectorXf::Zero(3 * verticesVector.size());
//
//	// ����?�_�Ǖ�������?
//	std::vector<int> indices = { 55, 35, 57, 70 };
//
//	// �{���͓I�召�C�ȍ���??��v?��
//	float forceMagnitude = 1000.0f;
//
//	// �@�ʐ���O?�C����wKey�I??����?�_�{����
//	if (this->groupIndex == 2) {
//		for (int localPi : indices) {
//			int indexX = 3 * localPi; // x�����I����
//
//			// ����???��wKey�I?��x������{����
//			if (wKey == 3) {
//				appliedForce(indexX) -= forceMagnitude; // �����{����
//			}
//			else if (wKey == 4) {
//				appliedForce(indexX) += forceMagnitude; // ���E�{����
//			}
//			// ?��wKey == 1��2�C��?�s��x������{����
//		}
//	}
//
//	// �X�VgroupVelocity
//	groupVelocity += appliedForce * timeStep;
//
//	// �g�p������??�ZvelocityUpdate
//	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;
//
//	// �X�VprimeVec�a?�_�ʒu
//	for (auto& vertexPair : verticesVector) {
//		Vertex* vertex = vertexPair;
//		int localPi = vertex->localIndex; // �g�p�Ǖ�����
//
//		// ?�擖�O?�_�I���x�X�V����
//		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);
//
//		// ?�Z�V�I�ʒu
//		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;
//
//		// �X�VprimeVec
//		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
//	}
//}
void Group::calPrimeVecS(int wKey) {
	// ?��primeVec��?���n����?�u?��?�I�ڐ�
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// ���n���꘢�͌���
	Eigen::VectorXf appliedForce = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// �{���͓I�召�C�ȍ���??��v?��
	float forceMagnitude = 10.0f;

	// �@�ʐ���O?�C����wKey�I???�����L?�_�{����
	if (this->groupIndex == 2) {
		for (size_t i = 0; i < verticesVector.size(); ++i) {
			int indexX = 3 * i; // x�����I����

			// ����???��wKey�I?��x������{����
			if (wKey == 3) {
				appliedForce(indexX) -= forceMagnitude; // �����{����
			}
			else if (wKey == 4) {
				appliedForce(indexX) += forceMagnitude; // ���E�{����
			}
			// ?��wKey == 1��2�C��?�s��x������{����
		}
	}

	// �X�VgroupVelocity
	groupVelocity += appliedForce * timeStep;

	// �g�p������??�ZvelocityUpdate
	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;

	// �X�VprimeVec�a?�_�ʒu
	for (auto& vertexPair : verticesVector) {
		Vertex* vertex = vertexPair;
		int localPi = vertex->localIndex; // �g�p�Ǖ�����

		// ?�擖�O?�_�I���x�X�V����
		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);

		// ?�Z�V�I�ʒu
		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;

		// �X�VprimeVec
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
	}
}


void Group::calPrimeVecT(int w) {
	// ?��primeVec��?���n����?�u?��?�I�ڐ�
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// ���n����?�͌���
	Eigen::VectorXf twistForce = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// ����?�_�Ǖ�������?
	std::vector<int> indices = { 22, 17, 9, 15 };
	std::vector<int> indices1 = { 24, 18, 5, 17 };
	if (this->groupIndex == 0) {
		// ��?�Z���ψʒu
		Eigen::Vector3f avgPosition = Eigen::Vector3f::Zero();
		for (int index : indices1) {
			Vertex* vertex = verticesVector[index];
			avgPosition += Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
		}
		avgPosition /= indices1.size();

		// ??��?�_�{����?��
		for (int index : indices1) {
			Vertex* vertex = verticesVector[index];
			Eigen::Vector3f pos = Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
			Eigen::Vector3f toCenter = pos - avgPosition;

			// ��yz���ʏ�?�Z���ؓI�͓I����
			Eigen::Vector3f tangentForce = Eigen::Vector3f(0, -toCenter.z(), toCenter.y()).normalized();

			// �X�VtwistForce�C?���g�ptangentForce���Ȉ꘢?�ʗ��\���͓I�召
			float forceMagnitude = 150.0f; // ����͓I�召
			twistForce(3 * index) += 0; // X������s�{����
			twistForce(3 * index + 1) += forceMagnitude * tangentForce.y(); // Y����
			twistForce(3 * index + 2) += forceMagnitude * tangentForce.z(); // Z����
		}
	}

	// ??���ې���O?�C�@�ʐ��C??����?�_�{����?��
	if (this->groupIndex == 3) {
		// ��?�Z���ψʒu
		Eigen::Vector3f avgPosition = Eigen::Vector3f::Zero();
		for (int index : indices) {
			Vertex* vertex = verticesVector[index];
			avgPosition += Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
		}
		avgPosition /= indices.size();

		// ??��?�_�{����?��
		for (int index : indices) {
			Vertex* vertex = verticesVector[index];
			Eigen::Vector3f pos = Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
			Eigen::Vector3f toCenter = pos - avgPosition;

			// ��yz���ʏ�?�Z���ؓI�͓I����
			Eigen::Vector3f tangentForce = Eigen::Vector3f(0, toCenter.z(), -toCenter.y()).normalized();

			// �X�VtwistForce�C?���g�ptangentForce���Ȉ꘢?�ʗ��\���͓I�召
			float forceMagnitude = 150.0f; // ����͓I�召
			twistForce(3 * index) += 0; // X������s�{����
			twistForce(3 * index + 1) += forceMagnitude * tangentForce.y(); // Y����
			twistForce(3 * index + 2) += forceMagnitude * tangentForce.z(); // Z����
		}
	}

	// �X�VgroupVelocity
	groupVelocity += twistForce * timeStep;

	// �g�p������??�ZvelocityUpdate
	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;

	// �X�VprimeVec�a?�_�ʒu
	for (auto& vertexPair : verticesVector) {
		Vertex* vertex = vertexPair;
		int localPi = vertex->localIndex; // �g�p�Ǖ�����

		// ?�擖�O?�_�I���x�X�V����
		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);

		// ?�Z�V�I�ʒu
		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;

		// �X�VprimeVec
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
	}
}


void Group::calPrimeVec2(int w) {
	// ... [?�L��?] ...
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());


	gravity = Eigen::VectorXf::Zero(3 * verticesVector.size());
	// ?�u�w��_�I gravity
	int localPi = 2; // �w��I localIndex
	int globalPi = 69;
	if (globalPi < verticesVector.size()) {
		Vertex* v = verticesVector[globalPi];
		int gravityIndex = 3 * (v->localIndex) + 1; // ��?�� y ������{����
		if (w == 1 || w == 2) {
			gravity(gravityIndex) = (w == 1) ? Gravity : -Gravity; // ���� w �I?�r��͓I����
		}
	}
	groupVelocity += gravity * timeStep;

	// �g�p������??�ZvelocityUpdate
	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;

	// �X�VprimeVec�a?�_�ʒu
	for (auto& vertexPair : verticesVector) {
		Vertex* vertex = vertexPair;
		int localPi = vertex->localIndex; // �g�p�Ǖ�����

		// ?�擖�O?�_�I���x�X�V����
		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);

		// ?�Z�V�I�ʒu
		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;

		// �X�VprimeVec
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;

	}
	// ... [�X�V groupVelocity �a primeVec �I��?] ...
}
void Group::calPrimeVec(int w) {
	// ?��groupVelocity��?���n����?�u?��?�I�ڐ�
	//groupVelocity = Eigen::VectorXf::Zero(3 * verticesMap.size());
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());


	gravity = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// ���n��gravity����
	if (w == 4) {
		for (int i = 0; i < 3 * verticesVector.size(); i += 3) {
			gravity(i) = -Gravity; // y������?�u�d�� �E
		}
	}
	else if (w == 2) {
		for (int i = 1; i < 3 * verticesVector.size(); i += 3) {
			gravity(i) = Gravity; // y������?�u�d�� ��
		}
	}
	else if (w == 1) {
		for (int i = 1; i < 3 * verticesVector.size(); i += 3) {
			gravity(i) = -Gravity; // y������?�u�d�� ��
		}
	}
	else if (w == 3) {
		for (int i = 0; i < 3 * verticesVector.size(); i += 3) {
			gravity(i) = Gravity; // y������?�u�d�� ��
		}
	}


	// �X�VgroupVelocity
	//groupVelocityFEM += gravity * timeStep;
	Eigen::VectorXf exfUpdate = timeStep * timeStep * inverseTerm * massMatrix * gravity;
	Eigen::VectorXf velocityUpdate = inverseTerm * massMatrix * groupVelocity * timeStep;

	// �X�VprimeVec�a?�_�ʒu
	for (auto& vertexPair : verticesVector) {
		Vertex* vertex = vertexPair;
		int localPi = vertex->localIndex; // �g�p�Ǖ�����

		// ?�擖�O?�_�I���x�X�V����
		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);
		Eigen::Vector3f currentExfUpdate = exfUpdate.segment<3>(3 * localPi);
		// ?�Z�V�I�ʒu
		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate + currentExfUpdate;

		// �X�VprimeVec
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
	}

}
void Group::calPrimeVec() {
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());

	if (!gravityApplied) {


		// ���n��gravity���ʁC����y�����{���d��
		for (int i = 1; i < 3 * verticesVector.size(); i += 3) {
			gravity(i) = Gravity; // y������?�u�d��
		}


		gravityApplied = true; // 
	}
	//groupVelocity += gravity * timeStep;
	
	//Eigen::VectorXf exfUpdate = timeStep * timeStep * massMatrix * gravity;
	//Eigen::VectorXf exfUpdate = timeStep * timeStep *inverseTerm * massMatrix * gravity;
	Eigen::VectorXf exfUpdate = timeStep * timeStep * inverseTerm * massMatrix * gravity;
	Eigen::VectorXf velocityUpdate = inverseTerm * massMatrix * groupVelocity * timeStep;

	// �X�VprimeVec�a?�_�ʒu
	for (auto& vertexPair : verticesVector) {
		Vertex* vertex = vertexPair;
		int localPi = vertex->localIndex; // �g�p�Ǖ�����
		if (vertexPair->isFixed == true)
		{
			primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = Eigen::Vector3f(vertex->initx, vertex->inity, vertex->initz);
		}
		else
		{
			Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);
			Eigen::Vector3f currentExfUpdate = exfUpdate.segment<3>(3 * localPi);
			// ?�Z�V�I�ʒu
			Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate + currentExfUpdate;

			// �X�VprimeVec
			primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
		}
			// ?�擖�O?�_�I���x�X�V����
			
	}
		
}

void Group::calLHS() {
	//A = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK;
	//B = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * massDistribution;

	float reference = 0.0f; // float?�^�I�Q�l?
	float epsilon = std::numeric_limits<float>::epsilon(); // float?�^�Iepsilon
	massDampingSparseInv = (massMatrix + timeStep * dampingMatrix).inverse().sparseView(reference, epsilon);
	LHS_A = timeStep * timeStep * massDampingSparseInv * kSparse;
	LHS_B = LHS_A * massDistributionSparse;

	// ?�Z�t��? ���ʕs����LHS�C?�֎Z
	inverseTerm = (massMatrix + dampingMatrix * timeStep).inverse(); //??�֔c?���Z��
	inverseTermSparse = inverseTerm.sparseView();
	RHS_E = timeStep * timeStep * massDampingSparseInv * kSparse;
	RHS_A = RHS_E * initLocalPos;

	FEMLHS = LHS_I + LHS_A - LHS_B;
	FEMLHS_Inv = FEMLHS.inverse();
}
void Group::calLHSFEM() {
	//A = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK;
	//B = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * massDistribution;

	float reference = 0.0f; // float?�^�I�Q�l?
	float epsilon = std::numeric_limits<float>::epsilon(); // float?�^�Iepsilon
	LHSFEM = (massMatrix.sparseView() + timeStep * dampingMatrix.sparseView() + timeStep * timeStep * kSparseFEM);
}
void Group::calRHS() {

	//Fbind = Eigen::VectorXf::Zero(3 * verticesMap.size());
	//A = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * initLocalPos;
	//B = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * rotationMatrix.transpose() * primeVec;
	//C = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * rotationMatrix.transpose() * massDistribution * primeVec;
	//D = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * rotationMatrix.inverse() * Fbind;
	//rotationTransSparse = rotationMatrix.transpose().sparseView();




	RHS_D = RHS_G * Fbind;
	FEMRHS = RHS_AsubBplusC + RHS_D;

}
void Group::calRHSFEM()
{
	RHSFEM = (massMatrix.sparseView() + timeStep * dampingMatrix.sparseView()) * primeVec + timeStep * timeStep * kSparseFEM * initLocalPos;
}
void Group::calRInvLocalPos() {
	curLocalPos.resize(3 * verticesMap.size());
	RInvPos.resize(3 * verticesMap.size());

	for (const auto& vertexPair : verticesMap) {
		const Vertex* vertex = vertexPair.second;
		Eigen::Vector3f local_position(vertex->x, vertex->y, vertex->z);
		// ?�Z���n�ʒu�^���n�d�S�I��?
		Eigen::Vector3f positiondifference = local_position - centerofMass;

		// ���Ǖ��ʒu��?��initLocalPos���C����index���v����3��??��?�_�L3����??
		curLocalPos.segment<3>(vertex->localIndex * 3) = positiondifference;
	}
	RInvPos = rotationMatrix.inverse() * curLocalPos;
	std::cout << RInvPos << std::endl;
}




void Group::calDeltaX() {

	// ��?������Ax = b
	deltaX = FEMLHS_Inv * FEMRHS;
	//deltaX = FEMLHS.colPivHouseholderQr().solve(FEMRHS);

	// �� FEMLHS ???�H�`��?
	//float threshold = 1e-18;
	//Eigen::SparseMatrix<float> sparseFEMLHS = FEMLHS.sparseView(threshold);

	//Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> solver;
	//solver.compute(sparseFEMLHS);
	//deltaX = solver.solve(FEMRHS);

	deltaX = rotationMatrix * deltaX;
}
void Group::calDeltaXFEM() {

	// ��?������Ax = b
	deltaXFEM = LHSFEM.inverse() * RHSFEM;
}
void Group::calculateCurrentPositions() {
	// ��?���L?�_
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localidx;
		localidx = vertex->localIndex;

		// ?��primeVec��???�_�I�ʒu
		Eigen::Vector3f primePosition = primeVec.segment<3>(3 * localidx);

		// ?��deltaX��???�_�I�ʈ�
		Eigen::Vector3f displacement = deltaX.segment<3>(3 * localidx);

		// ?�Z���O�ʒu
		Eigen::Vector3f currentPos = primePosition + displacement;
		currentPosition.segment<3>(3 * localidx) = currentPos;

	}
}
void Group::calculateCurrentPositionsFEM() {
	// ��?���L?�_
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localidx;
		localidx = vertex->localIndex;


		// ?�Z���O�ʒu
		Eigen::Vector3f currentPos = deltaXFEM.segment<3>(3 * localidx);
		currentPositionFEM.segment<3>(3 * localidx) = currentPos;

	}
}
//void Group::updateVertexPositions() {
//	for (auto& vertexPair : verticesMap) {
//		Vertex* vertex = vertexPair.second;
//
//		// �g�p�Ǖ�������?�搳?�I��??�aprimeVec����
//		int localIndex = vertex->localIndex;
//		Eigen::Matrix3f rotationBlock = rotationMatrix.block<3, 3>(3 * localIndex, 3 * localIndex);
//		Eigen::Vector3f positionInPrimeVec = primeVec.segment<3>(3 * localIndex);
//
//		if (vertex->isFixed) {
//			// ?���Œ�_�C���ʒu?�u?���n�ʒu
//			vertex->x = vertex->initx;
//			vertex->y = vertex->inity;
//			vertex->z = vertex->initz;
//		}
//		else {
//			// �g�p��?��??����primeVec���I�ʒu
//			Eigen::Vector3f newPosition = rotationBlock * positionInPrimeVec;
//
//			// �X�V?�_�ʒu
//			vertex->x = newPosition.x();
//			vertex->y = newPosition.y();
//			vertex->z = newPosition.z();
//		}
//	}
//}

void Group::calFbind1(const std::vector<Vertex*>& commonVerticesGroup1,
	const std::vector<Vertex*>& commonVerticesGroup2,
	const Eigen::VectorXf& currentPositionGroup1,
	const Eigen::VectorXf& currentPositionGroup2,
	float k) {
	// Initialize Fbind, with a length three times the number of vertices in the group
	//Fbind = Eigen::VectorXf::Zero(verticesMap.size() * 3);
	Eigen::Vector3f posThisGroup;
	Eigen::Vector3f posOtherGroup;
	Eigen::Vector3f avgPosition;
	Eigen::Vector3f posDifference;
	Eigen::Vector3f force;
	Eigen::VectorXf distances = Eigen::VectorXf::Zero(commonVerticesGroup1.size() * 3);
	posThisGroup = Eigen::Vector3f::Zero();
	posOtherGroup = Eigen::Vector3f::Zero();
	avgPosition = Eigen::Vector3f::Zero();
	posDifference = Eigen::Vector3f::Zero();
	force = Eigen::Vector3f::Zero();
	// Iterate through all common vertices
	for (size_t i = 0; i < commonVerticesGroup1.size(); ++i) {
		// Get the common vertex in this group and the other group
		Vertex* vertexThisGroup = commonVerticesGroup1[i];
		Vertex* vertexOtherGroup = commonVerticesGroup2[i];

		// Directly use x, y, z from the Vertex objects for position
		/*Eigen::Vector3f posThisGroup(vertexThisGroup->x, vertexThisGroup->y, vertexThisGroup->z);
		Eigen::Vector3f posOtherGroup(vertexOtherGroup->x, vertexOtherGroup->y, vertexOtherGroup->z);*/
		posThisGroup = currentPositionGroup1.segment<3>(3 * vertexThisGroup->localIndex);
		posOtherGroup = currentPositionGroup2.segment<3>(3 * vertexOtherGroup->localIndex);
		avgPosition = (posThisGroup + posOtherGroup) / 2;
		// Compute the position difference between the current vertex and the other group's vertex
		posDifference = posThisGroup - avgPosition;
		// Compute the constraint force
		force = k * posDifference;
		float maxForce = 100000;
		if (abs(force.x()) > maxForce)
		{
			force.x() = force.x() / abs(force.x()) * maxForce;
		}
		if (abs(force.y()) > maxForce)
		{
			force.y() = force.y() / abs(force.y()) * maxForce;
		}
		if (abs(force.z()) > maxForce)
		{
			force.z() = force.z() / abs(force.z()) * maxForce;
		}
		// Place the constraint force in Fbind at the appropriate position using the local index
		Fbind.segment<3>(3 * vertexThisGroup->localIndex) += force;
		distances.segment<3>(3 * i) = (posThisGroup - posOtherGroup).cwiseAbs();
	}

}

void Group::calFbind(const Eigen::VectorXf& currentPositionThisGroup, const std::vector<Eigen::VectorXf>& allCurrentPositionsOtherGroups, float k) {
	Eigen::Vector3f posThisGroup;
	Eigen::Vector3f posOtherGroup;
	Eigen::Vector3f avgPosition;
	Eigen::Vector3f posDifference;
	Eigen::Vector3f force;
	// ��?verticesMap����?��?��?�_�����f�˓I?��
	Fbind = Eigen::VectorXf::Zero(currentPositionThisGroup.size()); // ��??��?�_��p3���ʒu

	for (int direction = 0; direction < 6; ++direction) {
		int adjacentGroupIdx = adjacentGroupIDs[direction];

		if (adjacentGroupIdx != -1) {
			const auto& currentPositionOtherGroup = allCurrentPositionsOtherGroups[adjacentGroupIdx];
			const auto& commonVerticesPair = commonVerticesInDirections[direction];

			for (size_t i = 0; i < commonVerticesPair.first.size(); ++i) {
				Vertex* vertexThisGroup = commonVerticesPair.first[i];
				Vertex* vertexOtherGroup = commonVerticesPair.second[i];

				// ���ڎg�p???�_�I���O�ʒu
				posThisGroup = currentPositionThisGroup.segment<3>(3 * vertexThisGroup->localIndex);
				posOtherGroup = currentPositionOtherGroup.segment<3>(3 * vertexOtherGroup->localIndex);
				avgPosition = (posThisGroup + posOtherGroup) / 2;
				posDifference = posThisGroup - avgPosition;
				force = k * posDifference;

				Fbind.segment<3>(3 * vertexThisGroup->localIndex) = force;
			}
		}
	}
}


void Group::updatePosition() {
	static float frameTime = 0;
	frameTime += timeStep;
	Eigen::Vector3f pos = Eigen::Vector3f::Zero();
	// ��?���L?�_
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localIndex = vertex->localIndex;

		// ��currentPosition��?��???�_�I�ʒu
		pos = currentPosition.segment<3>(3 * localIndex);
		/*vertex->x = pos.x();
		vertex->y = pos.y();
		vertex->z = pos.z();*/
		if (vertex->isFixed == true) {
			// ?���Œ�_�C���ʒu?�u?���n�ʒu
			/*vertex->x = vertex->x = vertex->initx - 0.25 * sin(0.4 * frameTime);
			vertex->y = vertex->y = vertex->inity + 0.1 * sin(0.7 * frameTime);*/
			vertex->x = vertex->initx;
			vertex->y = vertex->inity;
			vertex->z = vertex->initz;
		}
		else {
			// �g�p��?��??����primeVec���I�ʒu


			vertex->x = pos.x();
			vertex->y = pos.y();
			vertex->z = pos.z();
		}

		/* �X�V?�_�I�ʒu
		if (vertex->isFixed == true)
		{
			vertex->x += 0.01;
		}*/
	}
}
void Group::updatePositionFEM() {
	Eigen::Vector3f pos = Eigen::Vector3f::Zero();
	// ��?���L?�_
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localIndex = vertex->localIndex;

		// ��currentPosition��?��???�_�I�ʒu
		pos = deltaXFEM.segment<3>(3 * localIndex);
		/*vertex->x = pos.x();
		vertex->y = pos.y();
		vertex->z = pos.z();*/
		if (vertex->isFixed) {
			// ?���Œ�_�C���ʒu?�u?���n�ʒu
			vertex->x = vertex->initx;
			vertex->y = vertex->inity;
			vertex->z = vertex->initz;
		}
		else {
			// �g�p��?��??����primeVec���I�ʒu


			vertex->x = pos.x();
			vertex->y = pos.y();
			vertex->z = pos.z();
		}

		// �X�V?�_�I�ʒu

	}
}
void Group::updateVelocity() {
	Eigen::Vector3f previousPos = Eigen::Vector3f::Zero();
	Eigen::Vector3f currentPos = Eigen::Vector3f::Zero();
	Eigen::Vector3f velocity = Eigen::Vector3f::Zero();
	Kinematics = 0.0;

	
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localIndex = vertex->localIndex;
		if (vertex->isFixed == false)
		{
			previousPos.x() = vertex->x;
			previousPos.y() = vertex->y;
			previousPos.z() = vertex->z;

			
			currentPos = currentPosition.segment<3>(3 * localIndex);

		
			velocity = (currentPos - previousPos) / timeStep;
			groupVelocity.segment<3>(3 * localIndex) = velocity;

			//Kinematics += 0.5 * vertex->vertexMass * velocity.norm();
		}
		/*std::cout << "\r";
		std::cout << groupVelocity << std::endl;*/
		

		// �X�V vertex �I���x
		// ��@�Fvertex->velocity = velocity;

		// �ۑ����O�ʒu��?����?�I�g���?�ʒu�h
		//previousPosition.segment<3>(3 * localIndex) = currentPos;
	}

	//// �X�VcurrentPosition?�{?�ō@�I�ʒu
	//for (auto& vertexPair : verticesMap) {
	//	Vertex* vertex = vertexPair.second;
	//	int localIndex = vertex->localIndex;

	//	currentPosition.segment<3>(3 * localIndex) = Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
	//}
}
void Group::updateVelocityFEM() {
	Eigen::Vector3f previousPos = Eigen::Vector3f::Zero();
	Eigen::Vector3f currentPos = Eigen::Vector3f::Zero();
	Eigen::Vector3f velocity = Eigen::Vector3f::Zero();

	// ��?���L?�_�C�X�V���x��ۑ����O�ʒu
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localIndex = vertex->localIndex;

		// ?�擖�O�ʒu
		previousPos.x() = vertex->x;
		previousPos.y() = vertex->y;
		previousPos.z() = vertex->z;

		// �� previousPosition ?����?�I�ʒu
		currentPos = currentPositionFEM.segment<3>(3 * localIndex);

		// ?�Z���x
		velocity = (currentPos - previousPos) / timeStep;
		groupVelocityFEM.segment<3>(3 * localIndex) = velocity;
		// �X�V vertex �I���x
		// ��@�Fvertex->velocity = velocity;

		// �ۑ����O�ʒu��?����?�I�g���?�ʒu�h
		//previousPosition.segment<3>(3 * localIndex) = currentPos;
	}

	//// �X�VcurrentPosition?�{?�ō@�I�ʒu
	//for (auto& vertexPair : verticesMap) {
	//	Vertex* vertex = vertexPair.second;
	//	int localIndex = vertex->localIndex;

	//	currentPosition.segment<3>(3 * localIndex) = Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
	//}
}