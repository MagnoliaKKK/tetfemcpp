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
	//	verticesMap[tet->vertices[i]->index] = tet->vertices[i]; //揧壛巐柺懱揑摨?丆攃巐柺懱揑?揰壛擖verticesMap
	//}
}
std::vector<Vertex*> Group::getUniqueVertices() { //?槩?惀廀梫揑丆憡摉槹攃hashmap??惉vertexGroup
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
	// ?曐initLocalPos桳懌?揑嬻?棃懚?強桳揑嬊晹埵抲
	initLocalPos.resize(3 * verticesMap.size());

	for (const auto& vertexPair : verticesMap) {
		const Vertex* vertex = vertexPair.second;
		Eigen::Vector3f initial_position(vertex->initx, vertex->inity, vertex->initz);
		// ?嶼弶巒埵抲梌弶巒廳怱揑嵎?
		Eigen::Vector3f local_position = initial_position - initCOM;

		// 彨嬊晹埵抲懚?嵼initLocalPos拞丆拲堄index廀梫槱埲3場??槩?揰桳3槩嵖??
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


	// ?嶼Apq嬮?
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

	// 弶巒壔巐尦悢榓慁?嬮?
	Eigen::Vector3f omega = Eigen::Vector3f::Identity();
	Eigen::Quaternionf quaternion(Eigen::Quaternionf::Identity());
	rotate_matrix = Eigen::Matrix3f::Identity();
	Eigen::Vector3f gradR = Eigen::Vector3f::Zero();
	Eigen::Matrix3f HesseR = Eigen::Matrix3f::Zero();
	Eigen::Matrix3f S = Eigen::Matrix3f::Zero();

	// 揜戙?漄嵟壚慁?
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

	// ?寶慁?嬮?揑3N x 3N斉杮
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
	groupMass = 0.0; // 弶巒壔?揑?検?0
	for (auto& tet : tetrahedra) { // 曊??堦槩巐柺懱
		groupMass += tet->massTetra; // 椵壛?堦槩巐柺懱揑?検摓?揑?検
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
//	// ?曐primeVec涍?弶巒壔妿?抲?惓?揑広悺
//	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());
//
//	// 弶巒壔慁?椡岦検
//	Eigen::VectorXf twistForce = Eigen::VectorXf::Zero(3 * verticesVector.size());
//
//	// ??惀斲惀戞嶰?丆擛壥惀丆??摿掕?揰巤壛慁?椡
//	if (this->groupIndex == 2) {
//		// ?抲?槩?揰慁?椡揑戝彫丆?棦?壔??堦戝彫丆嬶懱?廀梫崻悩??忣檝?惍
//		float forceMagnitude = 10.0f; // 帵椺椡揑戝彫
//
//		// ?揰嬊晹嶕堷悢?
//		std::vector<int> indices = { 116, 91, 24, 22 };
//
//		// ??槩?揰巤壛慁?椡
//		for (int index : indices) {
//			// ?掕椡揑曽岦丗崻悩?揰嵼yz柺忋揑???慁?曽岦巤壛
//			if (index == 116) {
//				// ?槹116丆橈?椡岦壓乮y??曽岦乯
//				twistForce(3 * index + 1) = -forceMagnitude;
//			}
//			else if (index == 91) {
//				// ?槹91丆橈?椡岦塃乮z?惓曽岦乯
//				twistForce(3 * index + 2) = forceMagnitude;
//			}
//			else if (index == 24) {
//				// ?槹24丆橈?椡岦忋乮y?惓曽岦乯
//				twistForce(3 * index + 1) = forceMagnitude;
//			}
//			else if (index == 22) {
//				// ?槹22丆橈?椡岦嵍乮z??曽岦乯
//				twistForce(3 * index + 2) = -forceMagnitude;
//			}
//		}
//	}
//
//	// 峏怴groupVelocity
//	groupVelocity += twistForce * timeStep;
//
//	// 巊梡惍槩嬮??嶼velocityUpdate
//	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;
//
//	// 峏怴primeVec榓?揰埵抲
//	for (auto& vertexPair : verticesVector) {
//		Vertex* vertex = vertexPair;
//		int localPi = vertex->localIndex; // 巊梡嬊晹嶕堷
//
//		// ?庢摉慜?揰揑懍搙峏怴晹暘
//		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);
//
//		// ?嶼怴揑埵抲
//		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;
//
//		// 峏怴primeVec
//		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
//	}
//}
//void Group::calPrimeVecS(int wKey) {
//	// ?曐primeVec涍?弶巒壔妿?抲?惓?揑広悺
//	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());
//
//	// 弶巒壔堦槩椡岦検
//	Eigen::VectorXf appliedForce = Eigen::VectorXf::Zero(3 * verticesVector.size());
//
//	// 摿掕?揰嬊晹嶕堷悢?
//	std::vector<int> indices = { 55, 35, 57, 70 };
//
//	// 巤壛椡揑戝彫丆壜埲崻悩??忣檝?惍
//	float forceMagnitude = 1000.0f;
//
//	// 擛壥惀戞嶰?丆崻悩wKey揑??摿掕?揰巤壛椡
//	if (this->groupIndex == 2) {
//		for (int localPi : indices) {
//			int indexX = 3 * localPi; // x曽岦揑嶕堷
//
//			// 崻悩???擖wKey揑?嵼x曽岦忋巤壛椡
//			if (wKey == 3) {
//				appliedForce(indexX) -= forceMagnitude; // 岦嵍巤壛椡
//			}
//			else if (wKey == 4) {
//				appliedForce(indexX) += forceMagnitude; // 岦塃巤壛椡
//			}
//			// ?槹wKey == 1埥2丆変?晄嵼x曽岦忋巤壛椡
//		}
//	}
//
//	// 峏怴groupVelocity
//	groupVelocity += appliedForce * timeStep;
//
//	// 巊梡惍槩嬮??嶼velocityUpdate
//	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;
//
//	// 峏怴primeVec榓?揰埵抲
//	for (auto& vertexPair : verticesVector) {
//		Vertex* vertex = vertexPair;
//		int localPi = vertex->localIndex; // 巊梡嬊晹嶕堷
//
//		// ?庢摉慜?揰揑懍搙峏怴晹暘
//		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);
//
//		// ?嶼怴揑埵抲
//		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;
//
//		// 峏怴primeVec
//		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
//	}
//}
void Group::calPrimeVecS(int wKey) {
	// ?曐primeVec涍?弶巒壔妿?抲?惓?揑広悺
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// 弶巒壔堦槩椡岦検
	Eigen::VectorXf appliedForce = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// 巤壛椡揑戝彫丆壜埲崻悩??忣檝?惍
	float forceMagnitude = 10.0f;

	// 擛壥惀戞嶰?丆崻悩wKey揑???撪強桳?揰巤壛椡
	if (this->groupIndex == 2) {
		for (size_t i = 0; i < verticesVector.size(); ++i) {
			int indexX = 3 * i; // x曽岦揑嶕堷

			// 崻悩???擖wKey揑?嵼x曽岦忋巤壛椡
			if (wKey == 3) {
				appliedForce(indexX) -= forceMagnitude; // 岦嵍巤壛椡
			}
			else if (wKey == 4) {
				appliedForce(indexX) += forceMagnitude; // 岦塃巤壛椡
			}
			// ?槹wKey == 1埥2丆変?晄嵼x曽岦忋巤壛椡
		}
	}

	// 峏怴groupVelocity
	groupVelocity += appliedForce * timeStep;

	// 巊梡惍槩嬮??嶼velocityUpdate
	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;

	// 峏怴primeVec榓?揰埵抲
	for (auto& vertexPair : verticesVector) {
		Vertex* vertex = vertexPair;
		int localPi = vertex->localIndex; // 巊梡嬊晹嶕堷

		// ?庢摉慜?揰揑懍搙峏怴晹暘
		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);

		// ?嶼怴揑埵抲
		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;

		// 峏怴primeVec
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
	}
}


void Group::calPrimeVecT(int w) {
	// ?曐primeVec涍?弶巒壔妿?抲?惓?揑広悺
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// 弶巒壔慁?椡岦検
	Eigen::VectorXf twistForce = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// 摿掕?揰嬊晹嶕堷悢?
	std::vector<int> indices = { 22, 17, 9, 15 };
	std::vector<int> indices1 = { 24, 18, 5, 17 };
	if (this->groupIndex == 0) {
		// 愭?嶼暯嬒埵抲
		Eigen::Vector3f avgPosition = Eigen::Vector3f::Zero();
		for (int index : indices1) {
			Vertex* vertex = verticesVector[index];
			avgPosition += Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
		}
		avgPosition /= indices1.size();

		// ??槩?揰巤壛慁?椡
		for (int index : indices1) {
			Vertex* vertex = verticesVector[index];
			Eigen::Vector3f pos = Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
			Eigen::Vector3f toCenter = pos - avgPosition;

			// 嵼yz暯柺忋?嶼憡愗揑椡揑曽岦
			Eigen::Vector3f tangentForce = Eigen::Vector3f(0, -toCenter.z(), toCenter.y()).normalized();

			// 峏怴twistForce丆?棦巊梡tangentForce槱埲堦槩?検棃昞帵椡揑戝彫
			float forceMagnitude = 150.0f; // 帵椺椡揑戝彫
			twistForce(3 * index) += 0; // X曽岦忋晄巤壛椡
			twistForce(3 * index + 1) += forceMagnitude * tangentForce.y(); // Y曽岦
			twistForce(3 * index + 2) += forceMagnitude * tangentForce.z(); // Z曽岦
		}
	}

	// ??惀斲惀戞嶰?丆擛壥惀丆??摿掕?揰巤壛慁?椡
	if (this->groupIndex == 3) {
		// 愭?嶼暯嬒埵抲
		Eigen::Vector3f avgPosition = Eigen::Vector3f::Zero();
		for (int index : indices) {
			Vertex* vertex = verticesVector[index];
			avgPosition += Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
		}
		avgPosition /= indices.size();

		// ??槩?揰巤壛慁?椡
		for (int index : indices) {
			Vertex* vertex = verticesVector[index];
			Eigen::Vector3f pos = Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
			Eigen::Vector3f toCenter = pos - avgPosition;

			// 嵼yz暯柺忋?嶼憡愗揑椡揑曽岦
			Eigen::Vector3f tangentForce = Eigen::Vector3f(0, toCenter.z(), -toCenter.y()).normalized();

			// 峏怴twistForce丆?棦巊梡tangentForce槱埲堦槩?検棃昞帵椡揑戝彫
			float forceMagnitude = 150.0f; // 帵椺椡揑戝彫
			twistForce(3 * index) += 0; // X曽岦忋晄巤壛椡
			twistForce(3 * index + 1) += forceMagnitude * tangentForce.y(); // Y曽岦
			twistForce(3 * index + 2) += forceMagnitude * tangentForce.z(); // Z曽岦
		}
	}

	// 峏怴groupVelocity
	groupVelocity += twistForce * timeStep;

	// 巊梡惍槩嬮??嶼velocityUpdate
	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;

	// 峏怴primeVec榓?揰埵抲
	for (auto& vertexPair : verticesVector) {
		Vertex* vertex = vertexPair;
		int localPi = vertex->localIndex; // 巊梡嬊晹嶕堷

		// ?庢摉慜?揰揑懍搙峏怴晹暘
		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);

		// ?嶼怴揑埵抲
		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;

		// 峏怴primeVec
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
	}
}


void Group::calPrimeVec2(int w) {
	// ... [?桳戙?] ...
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());


	gravity = Eigen::VectorXf::Zero(3 * verticesVector.size());
	// ?抲巜掕揰揑 gravity
	int localPi = 2; // 巜掕揑 localIndex
	int globalPi = 69;
	if (globalPi < verticesVector.size()) {
		Vertex* v = verticesVector[globalPi];
		int gravityIndex = 3 * (v->localIndex) + 1; // 橈?嵼 y 曽岦忋巤壛椡
		if (w == 1 || w == 2) {
			gravity(gravityIndex) = (w == 1) ? Gravity : -Gravity; // 崻悩 w 揑?檙掕椡揑曽岦
		}
	}
	groupVelocity += gravity * timeStep;

	// 巊梡惍槩嬮??嶼velocityUpdate
	Eigen::VectorXf velocityUpdate = inverseTermSparse * (massMatrix * groupVelocity) * timeStep;

	// 峏怴primeVec榓?揰埵抲
	for (auto& vertexPair : verticesVector) {
		Vertex* vertex = vertexPair;
		int localPi = vertex->localIndex; // 巊梡嬊晹嶕堷

		// ?庢摉慜?揰揑懍搙峏怴晹暘
		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);

		// ?嶼怴揑埵抲
		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;

		// 峏怴primeVec
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;

	}
	// ... [峏怴 groupVelocity 榓 primeVec 揑戙?] ...
}
void Group::calPrimeVec(int w) {
	// ?曐groupVelocity涍?弶巒壔妿?抲?惓?揑広悺
	//groupVelocity = Eigen::VectorXf::Zero(3 * verticesMap.size());
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());


	gravity = Eigen::VectorXf::Zero(3 * verticesVector.size());

	// 弶巒壔gravity岦検
	if (w == 4) {
		for (int i = 0; i < 3 * verticesVector.size(); i += 3) {
			gravity(i) = -Gravity; // y曽岦忋?抲廳椡 塃
		}
	}
	else if (w == 2) {
		for (int i = 1; i < 3 * verticesVector.size(); i += 3) {
			gravity(i) = Gravity; // y曽岦忋?抲廳椡 壓
		}
	}
	else if (w == 1) {
		for (int i = 1; i < 3 * verticesVector.size(); i += 3) {
			gravity(i) = -Gravity; // y曽岦忋?抲廳椡 忋
		}
	}
	else if (w == 3) {
		for (int i = 0; i < 3 * verticesVector.size(); i += 3) {
			gravity(i) = Gravity; // y曽岦忋?抲廳椡 嵍
		}
	}


	// 峏怴groupVelocity
	//groupVelocityFEM += gravity * timeStep;
	Eigen::VectorXf exfUpdate = timeStep * timeStep * inverseTerm * massMatrix * gravity;
	Eigen::VectorXf velocityUpdate = inverseTerm * massMatrix * groupVelocity * timeStep;

	// 峏怴primeVec榓?揰埵抲
	for (auto& vertexPair : verticesVector) {
		Vertex* vertex = vertexPair;
		int localPi = vertex->localIndex; // 巊梡嬊晹嶕堷

		// ?庢摉慜?揰揑懍搙峏怴晹暘
		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);
		Eigen::Vector3f currentExfUpdate = exfUpdate.segment<3>(3 * localPi);
		// ?嶼怴揑埵抲
		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate + currentExfUpdate;

		// 峏怴primeVec
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
	}

}
void Group::calPrimeVec() {
	primeVec = Eigen::VectorXf::Zero(3 * verticesVector.size());

	if (!gravityApplied) {
		for (int i = 0; i < 3 * verticesVector.size(); i += 3) {
			gravity(i) = -Gravity; 
		}


		gravityApplied = true; // 
	}
	//groupVelocity += gravity * timeStep;
	
	//Eigen::VectorXf exfUpdate = timeStep * timeStep * massMatrix * gravity;
	//Eigen::VectorXf exfUpdate = timeStep * timeStep *inverseTerm * massMatrix * gravity;
	Eigen::VectorXf exfUpdate = timeStep * timeStep * inverseTerm * massMatrix * gravity;
	Eigen::VectorXf velocityUpdate = inverseTerm * massMatrix * groupVelocity * timeStep;

	for (auto& vertexPair : verticesVector) {
		Vertex* vertex = vertexPair;
		int localPi = vertex->localIndex;
	/*	if (vertexPair->isFixed == true)
		{
			primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = Eigen::Vector3f(vertex->initx, vertex->inity, vertex->initz);
		}
		else
		{
			Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);
			Eigen::Vector3f currentExfUpdate = exfUpdate.segment<3>(3 * localPi);
			Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate + currentExfUpdate;
			primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
		}*/
		Eigen::Vector3f currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);
		Eigen::Vector3f currentExfUpdate = exfUpdate.segment<3>(3 * localPi);
		Eigen::Vector3f newPosition = Eigen::Vector3f(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate + currentExfUpdate;
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
			
			
	}
		
}

void Group::calLHS() {
	//A = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK;
	//B = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * massDistribution;

	float reference = 0.0f; // float?宆揑嶲峫?
	float epsilon = std::numeric_limits<float>::epsilon(); // float?宆揑epsilon
	massDampingSparseInv = (massMatrix + timeStep * dampingMatrix).inverse().sparseView(reference, epsilon);
	LHS_A = timeStep * timeStep * massDampingSparseInv * kSparse;
	LHS_B = LHS_A * massDistributionSparse;

	// ?嶼媡嬮? 壓柺晄懏槹LHS丆?曋嶼
	inverseTerm = (massMatrix + dampingMatrix * timeStep).inverse(); //??曋攃?槩嶼椆
	inverseTermSparse = inverseTerm.sparseView();
	RHS_E = timeStep * timeStep * massDampingSparseInv * kSparse;
	RHS_A = RHS_E * initLocalPos;

	FEMLHS = LHS_I + LHS_A - LHS_B;
	FEMLHS_Inv = FEMLHS.inverse();
}
void Group::calLHSFEM() {
	//A = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK;
	//B = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * massDistribution;

	float reference = 0.0f; // float?宆揑嶲峫?
	float epsilon = std::numeric_limits<float>::epsilon(); // float?宆揑epsilon
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
		// ?嶼弶巒埵抲梌弶巒廳怱揑嵎?
		Eigen::Vector3f positiondifference = local_position - centerofMass;

		// 彨嬊晹埵抲懚?嵼initLocalPos拞丆拲堄index廀梫槱埲3場??槩?揰桳3槩嵖??
		curLocalPos.segment<3>(vertex->localIndex * 3) = positiondifference;
	}
	RInvPos = rotationMatrix.inverse() * curLocalPos;
	std::cout << RInvPos << std::endl;
}




void Group::calDeltaX() {

	// 夝?惈曽掱Ax = b
	deltaX = FEMLHS_Inv * FEMRHS;
	//deltaX = FEMLHS.colPivHouseholderQr().solve(FEMRHS);

	// 彨 FEMLHS ???婬慲嬮?
	//float threshold = 1e-18;
	//Eigen::SparseMatrix<float> sparseFEMLHS = FEMLHS.sparseView(threshold);

	//Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> solver;
	//solver.compute(sparseFEMLHS);
	//deltaX = solver.solve(FEMRHS);

	deltaX = rotationMatrix * deltaX;
}
void Group::calDeltaXFEM() {

	// 夝?惈曽掱Ax = b
	deltaXFEM = LHSFEM.inverse() * RHSFEM;
}
void Group::calculateCurrentPositions() {
	// 曊?強桳?揰
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localidx;
		localidx = vertex->localIndex;

		// ?庢primeVec拞???揰揑埵抲
		Eigen::Vector3f primePosition = primeVec.segment<3>(3 * localidx);

		// ?庢deltaX拞???揰揑埵堏
		Eigen::Vector3f displacement = deltaX.segment<3>(3 * localidx);

		// ?嶼摉慜埵抲
		Eigen::Vector3f currentPos = primePosition + displacement;
		currentPosition.segment<3>(3 * localidx) = currentPos;

	}
}
void Group::calculateCurrentPositionsFEM() {
	// 曊?強桳?揰
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localidx;
		localidx = vertex->localIndex;


		// ?嶼摉慜埵抲
		Eigen::Vector3f currentPos = deltaXFEM.segment<3>(3 * localidx);
		currentPositionFEM.segment<3>(3 * localidx) = currentPos;

	}
}
//void Group::updateVertexPositions() {
//	for (auto& vertexPair : verticesMap) {
//		Vertex* vertex = vertexPair.second;
//
//		// 巊梡嬊晹嶕堷棃?庢惓?揑嬮??榓primeVec晹暘
//		int localIndex = vertex->localIndex;
//		Eigen::Matrix3f rotationBlock = rotationMatrix.block<3, 3>(3 * localIndex, 3 * localIndex);
//		Eigen::Vector3f positionInPrimeVec = primeVec.segment<3>(3 * localIndex);
//
//		if (vertex->isFixed) {
//			// ?槹屌掕揰丆彨埵抲?抲?弶巒埵抲
//			vertex->x = vertex->initx;
//			vertex->y = vertex->inity;
//			vertex->z = vertex->initz;
//		}
//		else {
//			// 巊梡慁?嬮??槱埲primeVec拞揑埵抲
//			Eigen::Vector3f newPosition = rotationBlock * positionInPrimeVec;
//
//			// 峏怴?揰埵抲
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

void Group::calBindFixed() {
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		if (vertex->isFixed) {
			Eigen::Vector3f initPos(vertex->initx, vertex->inity, vertex->initz);
			Eigen::Vector3f currentPos = currentPosition.segment<3>(3 * vertex->localIndex);
			Eigen::Vector3f diff = currentPos - initPos;

			// Compute the constraint force for fixed vertices
			Eigen::Vector3f constraintForce = -100 * diff;
			Fbind.segment<3>(3 * vertex->localIndex) += constraintForce;
		}
	}
}

void Group::calFbind(const Eigen::VectorXf& currentPositionThisGroup, const std::vector<Eigen::VectorXf>& allCurrentPositionsOtherGroups, float k) {
	Eigen::Vector3f posThisGroup;
	Eigen::Vector3f posOtherGroup;
	Eigen::Vector3f avgPosition;
	Eigen::Vector3f posDifference;
	Eigen::Vector3f force;
	// 橈?verticesMap惀掕?椆?撪?揰嶕堷塮幩揑?検
	Fbind = Eigen::VectorXf::Zero(currentPositionThisGroup.size()); // 橈??槩?揰愯梡3槩埵抲

	for (int direction = 0; direction < 6; ++direction) {
		int adjacentGroupIdx = adjacentGroupIDs[direction];

		if (adjacentGroupIdx != -1) {
			const auto& currentPositionOtherGroup = allCurrentPositionsOtherGroups[adjacentGroupIdx];
			const auto& commonVerticesPair = commonVerticesInDirections[direction];

			for (size_t i = 0; i < commonVerticesPair.first.size(); ++i) {
				Vertex* vertexThisGroup = commonVerticesPair.first[i];
				Vertex* vertexOtherGroup = commonVerticesPair.second[i];

				// 捈愙巊梡???揰揑摉慜埵抲
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


//

void Group::updatePositionFEM() {
	Eigen::Vector3f pos = Eigen::Vector3f::Zero();
	// 曊?強桳?揰
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localIndex = vertex->localIndex;

		// 樃currentPosition拞?庢???揰揑埵抲
		pos = deltaXFEM.segment<3>(3 * localIndex);
		/*vertex->x = pos.x();
		vertex->y = pos.y();
		vertex->z = pos.z();*/
		if (vertex->isFixed) {
			// ?槹屌掕揰丆彨埵抲?抲?弶巒埵抲
			vertex->x = vertex->initx;
			vertex->y = vertex->inity;
			vertex->z = vertex->initz;
		}
		else {
			// 巊梡慁?嬮??槱埲primeVec拞揑埵抲


			vertex->x = pos.x();
			vertex->y = pos.y();
			vertex->z = pos.z();
		}

		// 峏怴?揰揑埵抲

	}
}
void Group::updatePosition() {
	static float frameTime = 0;
	frameTime += timeStep;
	Eigen::Vector3f pos = Eigen::Vector3f::Zero();
	// 遍?所有?点
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localIndex = vertex->localIndex;

		// 从currentPosition中?取???点的位置
		pos = currentPosition.segment<3>(3 * localIndex);

		vertex->x = pos.x();
		vertex->y = pos.y();
		vertex->z = pos.z();
		/*vertex->x = pos.x();
		vertex->y = pos.y();
		vertex->z = pos.z();*/
		//if (vertex->isFixed == true) {
		//	// ?于固定点，将位置?置?初始位置
		//	/*vertex->x = vertex->x = vertex->initx - 0.25 * sin(0.4 * frameTime);
		//	vertex->y = vertex->y = vertex->inity + 0.1 * sin(0.7 * frameTime);*/
		//	vertex->x = vertex->initx;
		//	vertex->y = vertex->inity;
		//	vertex->z = vertex->initz;
		//}
		//else {
		//	// 使用旋?矩??乘以primeVec中的位置


		//	vertex->x = pos.x();
		//	vertex->y = pos.y();
		//	vertex->z = pos.z();
		//}

		/* 更新?点的位置
		if (vertex->isFixed == true)
		{
			vertex->x += 0.01;
		}*/
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
		

		// 峏怴 vertex 揑懍搙
		// 椺擛丗vertex->velocity = velocity;

		// 曐懚摉慜埵抲嶌?壓堦?揑乬忋堦?埵抲乭
		//previousPosition.segment<3>(3 * localIndex) = currentPos;
	}

	//// 峏怴currentPosition?杮?嵟岪揑埵抲
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

	// 曊?強桳?揰丆峏怴懍搙涹曐懚摉慜埵抲
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localIndex = vertex->localIndex;

		// ?庢摉慜埵抲
		previousPos.x() = vertex->x;
		previousPos.y() = vertex->y;
		previousPos.z() = vertex->z;

		// 樃 previousPosition ?庢忋堦?揑埵抲
		currentPos = currentPositionFEM.segment<3>(3 * localIndex);

		// ?嶼懍搙
		velocity = (currentPos - previousPos) / timeStep;
		groupVelocityFEM.segment<3>(3 * localIndex) = velocity;
		// 峏怴 vertex 揑懍搙
		// 椺擛丗vertex->velocity = velocity;

		// 曐懚摉慜埵抲嶌?壓堦?揑乬忋堦?埵抲乭
		//previousPosition.segment<3>(3 * localIndex) = currentPos;
	}

	//// 峏怴currentPosition?杮?嵟岪揑埵抲
	//for (auto& vertexPair : verticesMap) {
	//	Vertex* vertex = vertexPair.second;
	//	int localIndex = vertex->localIndex;

	//	currentPosition.segment<3>(3 * localIndex) = Eigen::Vector3f(vertex->x, vertex->y, vertex->z);
	//}
}