#include "GroupDivision.h"

double timeStep = 0.01;
double dampingConst = 4;
const  double PI = 3.14159265358979265358979;
const double Gravity = -9.8;

void Object::assignLocalIndicesToAllGroups() { // local index generation
	for (Group& group : groups) {
		int currentLocalIndex = 0;
		std::unordered_set<Vertex*> processedVertices; // 用于跟踪已处理的顶点

		for (Tetrahedron* tetra : group.tetrahedra) {
			for (int i = 0; i < 4; ++i) {
				Vertex* vertex = tetra->vertices[i];

				// 检查顶点是否已经处理过
				if (processedVertices.find(vertex) == processedVertices.end()) {
					vertex->localIndex = currentLocalIndex++; // 分配本地索引
					processedVertices.insert(vertex); // 标记为已处理
				}
			}
		}
	}
}

void Object::updateIndices() {
	std::unordered_set<int> globalIndices;
	std::unordered_map<int, Vertex*> indexToVertexMap; // 旧索引到新顶点的映射
	int maxIndex = 0;

	// 首先遍历所有顶点以找到最大索引值
	for (Group& group : groups) {
		for (Tetrahedron* tetra : group.tetrahedra) {
			for (int i = 0; i < 4; ++i) {
				Vertex* vertex = tetra->vertices[i];
				maxIndex = std::max(maxIndex, vertex->index);
			}
		}
	}

	int nextAvailableIndex = maxIndex + 1;

	// 再次遍历所有顶点以更新索引
	for (Group& group : groups) {
		std::unordered_set<int> localIndices; // 每个组内的本地索引集合

		for (Tetrahedron* tetra : group.tetrahedra) {
			for (int i = 0; i < 4; ++i) {
				Vertex* vertex = tetra->vertices[i];

				if (localIndices.find(vertex->index) == localIndices.end()) {
					localIndices.insert(vertex->index);

					if (globalIndices.find(vertex->index) != globalIndices.end()) {
						// 如果索引已在全局集合中，创建新顶点并更新映射
						Vertex* newVertex = new Vertex(vertex->x, vertex->y, vertex->z, nextAvailableIndex++);
						indexToVertexMap[vertex->index] = newVertex;
						tetra->vertices[i] = newVertex;
						vertex = newVertex;
					}
					globalIndices.insert(vertex->index);
				}
				else if (indexToVertexMap.find(vertex->index) != indexToVertexMap.end()) {
					// 更新为新的顶点引用
					tetra->vertices[i] = indexToVertexMap[vertex->index];
				}
			}
		}
	}
}

void Object::generateUniqueVertices() { //执行这个函数以后，verticesMap就会装满这个组的vertices， 不重复
	//std::vector<Vertex*> uniqueVertices;
	
	for (Group& group : groups) {
		group.verticesMap.clear(); // 清空现有的映射

		for (Tetrahedron* tetra : group.tetrahedra) {
			for (int i = 0; i < 4; ++i) {
				Vertex* vertex = tetra->vertices[i];

				// 如果顶点尚未在verticesMap中，则添加
				if (group.verticesMap.find(vertex->index) == group.verticesMap.end()) {
					group.verticesMap[vertex->index] = vertex;
				}
			}
		}
		group.initialize();
	}

}

std::pair<std::vector<Vertex*>, std::vector<Vertex*>> Object::findCommonVertices(const Group& group1, const Group& group2) { //寻找共同点
	std::vector<Vertex*> commonVerticesGroup1;
	std::vector<Vertex*> commonVerticesGroup2;

	// 遍历group1的verticesMap中的所有顶点
	for (auto& mapEntry1 : group1.verticesMap) {
		Vertex* vertex1 = mapEntry1.second;

		// 遍历group2的verticesMap中的所有顶点
		for (auto& mapEntry2 : group2.verticesMap) {
			Vertex* vertex2 = mapEntry2.second;

			// 检查坐标是否相同
			if (vertex1->x == vertex2->x && vertex1->y == vertex2->y && vertex1->z == vertex2->z) {
				commonVerticesGroup1.push_back(vertex1);
				commonVerticesGroup2.push_back(vertex2);
			}
		}
	}

	return { commonVerticesGroup1, commonVerticesGroup2 };
}

void Group::initialize() {
	groupVelocity = Eigen::VectorXd::Zero(3 * verticesMap.size());
	Fbind = Eigen::VectorXd(3 * verticesMap.size());
	currentPosition = Eigen::VectorXd::Zero(3 * verticesMap.size());;
	gravity = Eigen::VectorXd::Zero(3 * verticesMap.size());
}

void Group::addTetrahedron(Tetrahedron* tet) {
	tetrahedra.push_back(tet);
	//for (int i = 0; i < 4; ++i) {
	//	verticesMap[tet->vertices[i]->index] = tet->vertices[i]; //添加四面体的同时，把四面体的顶点加入verticesMap
	//}
}
std::vector<Vertex*> Group::getUniqueVertices() { //这个还是需要的，相当于把hashmap转换成vertexGroup
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
void Group::calLocalPos() {
	// 确保initLocalPos有足够的空间来存储所有的局部位置
	initLocalPos.resize(3 * verticesMap.size());

	for (const auto& vertexPair : verticesMap) {
		const Vertex* vertex = vertexPair.second;
		Eigen::Vector3d initial_position(vertex->initx, vertex->inity, vertex->initz);
		// 计算初始位置与初始重心的差值
		Eigen::Vector3d local_position = initial_position - initCOM;

		// 将局部位置存储在initLocalPos中，注意index需要乘以3因为每个顶点有3个坐标值
		initLocalPos.segment<3>(vertex->localIndex * 3) = local_position;
	}
}
Eigen::Vector3d Group::axlAPD(Eigen::Matrix3d a) {
	Eigen::Vector3d g = Eigen::Vector3d::Zero();
	g[0] = a(1, 2) - a(2, 1);
	g[1] = a(2, 0) - a(0, 2);
	g[2] = a(0, 1) - a(1, 0);
	return g;
}
Eigen::Vector3d Group::clamp2(Eigen::Vector3d x, double y, double z) {
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
Eigen::Quaterniond Group::Exp2(Eigen::Vector3d a) {
	double s = sin((a * 0.5).norm());
	double x = s * a.x() / a.norm();
	double y = s * a.y() / a.norm();
	double z = s * a.z() / a.norm();
	Eigen::Quaterniond qq = Eigen::Quaterniond(cos((a * 0.5).norm()), x, y, z);
	return  qq;
}
void Group::calRotationMatrix() {
	Eigen::Matrix3d Apq = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d tempA = Eigen::Matrix3d::Zero();
	Eigen::Vector3d center_grid = Eigen::Vector3d::Zero();
	
	for (const auto& vertexEntry : verticesMap) {
		center_grid[0] = massDistribution(0, 3 * vertexEntry.second->localIndex) * primeVec[3 * vertexEntry.second->localIndex];
		center_grid[1] = massDistribution(0, 3 * vertexEntry.second->localIndex) * primeVec[3 * vertexEntry.second->localIndex + 1];
		center_grid[2] = massDistribution(0, 3 * vertexEntry.second->localIndex) * primeVec[3 * vertexEntry.second->localIndex + 2];
	}
	// 计算Apq矩阵
	for (const auto& vertexEntry : verticesMap) {
		
		tempA = (primeVec.block<3, 1>(3 * vertexEntry.second->localIndex, 0) - center_grid) * initLocalPos.block<3, 1>(3 * vertexEntry.second->localIndex, 0).transpose();
		Apq += massMatrix(3 * vertexEntry.second->localIndex, 3 * vertexEntry.second->localIndex) * tempA;
	}

	// 初始化四元数和旋转矩阵
	Eigen::Vector3d omega = Eigen::Vector3d::Identity();
	Eigen::Quaterniond quaternion(Eigen::Quaterniond::Identity());
	Eigen::Matrix3d rotate_matrix = Eigen::Matrix3d::Identity();
	Eigen::Vector3d gradR = Eigen::Vector3d::Zero();
	Eigen::Matrix3d HesseR = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d S = Eigen::Matrix3d::Zero();

	// 迭代寻找最佳旋转
	for (unsigned int ci = 0; ci < 20; ci++) {
		Eigen::Matrix3d R = quaternion.matrix();
		Eigen::Matrix3d S = R.transpose() * Apq;
		Eigen::Vector3d gradR = axlAPD(S);
		Eigen::Matrix3d HesseR = S.trace() * Eigen::Matrix3d::Identity() - (S + S.transpose()) * 0.5;
		Eigen::Vector3d omega = -1 * HesseR.inverse() * gradR;

		double w;
		w = omega.norm();
		if (w < 1.0e-9) {
			break;
		}

		omega = clamp2(omega, -1 * PI, PI);
		Eigen::Quaterniond  temp2;
		temp2 = Exp2(omega);
		quaternion = quaternion * temp2;
	}

	rotate_matrix = quaternion.matrix();

	// 构建旋转矩阵的3N x 3N版本
	rotationMatrix = Eigen::MatrixXd::Zero(3 * verticesMap.size(), 3 * verticesMap.size());
	for (unsigned int pi = 0; pi < verticesMap.size(); pi++) {
		rotationMatrix.block<3, 3>(3 * pi, 3 * pi) = rotate_matrix;
	}

	//rotationSparse = rotationMatrix.sparseView();
	//rotationTransSparse = rotationMatrix.transpose().sparseView();

	// 将旋转矩阵转换为稀疏格式（如果需要）
	//Eigen::SparseMatrix<double> Rn_Matrix_Sparse = rotate_matrix3N.sparseView();
	//Eigen::SparseMatrix<double> Rn_MatrixTR_Sparse = rotate_matrix3N.transpose().sparseView();
}
//
Eigen::MatrixXd Tetrahedron::createElementK(double E, double nu, const Eigen::Vector3d& groupCenterOfMass) {
	// 定义节点坐标
	double x1 = vertices[0]->x - groupCenterOfMass.x();
	double y1 = vertices[0]->y - groupCenterOfMass.y();
	double z1 = vertices[0]->z - groupCenterOfMass.z();
	double x2 = vertices[1]->x - groupCenterOfMass.x();
	double y2 = vertices[1]->y - groupCenterOfMass.y();
	double z2 = vertices[1]->z - groupCenterOfMass.z();
	double x3 = vertices[2]->x - groupCenterOfMass.x();
	double y3 = vertices[2]->y - groupCenterOfMass.y();
	double z3 = vertices[2]->z - groupCenterOfMass.z();
	double x4 = vertices[3]->x - groupCenterOfMass.x();
	double y4 = vertices[3]->y - groupCenterOfMass.y();
	double z4 = vertices[3]->z - groupCenterOfMass.z();

		// 构建体积计算矩阵 A
	Eigen::Matrix4d A;
	A << x1, y1, z1, 1,
		 x2, y2, z2, 1,
		 x3, y3, z3, 1,
		 x4, y4, z4, 1;

		// 计算四面体的体积
	double V = std::abs(A.determinant() / 6);

		// 定义 mbeta, mgamma, mdelta 矩阵
	Eigen::Matrix3d mbeta1, mbeta2, mbeta3, mbeta4, mgamma1, mgamma2, mgamma3, mgamma4, mdelta1, mdelta2, mdelta3, mdelta4;


	mbeta1 << 1, y2, z2, 1, y3, z3, 1, y4, z4;
	mbeta2 << 1, y1, z1, 1, y3, z3, 1, y4, z4;
	mbeta3 << 1, y1, z1, 1, y2, z2, 1, y4, z4;
	mbeta4 << 1, y1, z1, 1, y2, z2, 1, y3, z3;

	mgamma1 << 1, x2, z2, 1, x3, z3, 1, x4, z4;
	mgamma2 << 1, x1, z1, 1, x3, z3, 1, x4, z4;
	mgamma3 << 1, x1, z1, 1, x2, z2, 1, x4, z4;
	mgamma4 << 1, x1, z1, 1, x2, z2, 1, x3, z3;

	mdelta1 << 1, x2, y2, 1, x3, y3, 1, x4, y4;
	mdelta2 << 1, x1, y1, 1, x3, y3, 1, x4, y4;
	mdelta3 << 1, x1, y1, 1, x2, y2, 1, x4, y4;
	mdelta4 << 1, x1, y1, 1, x2, y2, 1, x3, y3;

		// 计算 beta, gamma 和 delta 值
	double beta1 = -mbeta1.determinant();
	double beta2 = mbeta2.determinant();
	double beta3 = -mbeta3.determinant();
	double beta4 = mbeta4.determinant();

	double gamma1 = mgamma1.determinant();
	double gamma2 = -mgamma2.determinant();
	double gamma3 = mgamma3.determinant();
	double gamma4 = -mgamma4.determinant();

	double delta1 = -mdelta1.determinant();
	double delta2 = mdelta2.determinant();
	double delta3 = -mdelta3.determinant();
	double delta4 = mdelta4.determinant();

	// 定义 B 矩阵
	Eigen::MatrixXd B(6, 12);
	
	B << beta1, 0, 0, beta2, 0, 0, beta3, 0, 0, beta4, 0, 0,
		0, gamma1, 0, 0, gamma2, 0, 0, gamma3, 0, 0, gamma4, 0,
		0, 0, delta1, 0, 0, delta2, 0, 0, delta3, 0, 0, delta4,
		gamma1, beta1, 0, gamma2, beta2, 0, gamma3, beta3, 0, gamma4, beta4, 0,
		0, delta1, gamma1, 0, delta2, gamma2, 0, delta3, gamma3, 0, delta4, gamma4,
		delta1, 0, beta1, delta2, 0, beta2, delta3, 0, beta3, delta4, 0, beta4;

	B /= (6 * V);

		// 定义材料属性矩阵 D
	    // 泊松比
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(6, 6);
	
	D << 1 - nu, nu, nu, 0, 0, 0,
		nu, 1 - nu, nu, 0, 0, 0,
		nu, nu, 1 - nu, 0, 0, 0,
		0, 0, 0, (1 - 2 * nu) / 2, 0, 0,
		0, 0, 0, 0, (1 - 2 * nu) / 2, 0,
		0, 0, 0, 0, 0, (1 - 2 * nu) / 2;

	D *= (E / ((1 + nu) * (1 - 2 * nu)));

	// 计算刚度矩阵 k
	Eigen::MatrixXd k=V * (B.transpose() * D * B);
	
	elementK = k;
	return k;


}
void Group::calGroupK(double E, double nu) {
	// Initialize groupK to the right size. Assuming 3 degrees of freedom per vertex
	int dof = verticesMap.size() * 3;
	groupK.resize(dof, dof);
	groupK.setZero();

	// Iterate over each tetrahedron to assemble the global stiffness matrix
	for (auto& tetra : tetrahedra) {
		// Get the local stiffness matrix for the current tetrahedron
		Eigen::MatrixXd localK = tetra->createElementK(E, nu, initCOM);

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
			int localIdx = tet->vertices[i]->localIndex;  // Use local index
			M(3 * localIdx, 3 * localIdx) += vertexMass;     // x-direction
			M(3 * localIdx + 1, 3 * localIdx + 1) += vertexMass; // y-direction
			M(3 * localIdx + 2, 3 * localIdx + 2) += vertexMass; // z-direction
		}
	}
	massMatrix = M;
	return M;
	// massSparse = massMatrix.sparseView();  // Commented out as it's not part of the provided code snippet
}


void Group::calDampingMatrix() {
	dampingMatrix.setZero();
	dampingMatrix = dampingConst * massMatrix;
	dampingSparse = dampingMatrix.sparseView();
}

void Group::calMassDistributionMatrix() {
	int N = verticesMap.size();  // Number of unique vertices
	massDistribution = Eigen::MatrixXd::Zero(3 * N, 3 * N);
	Eigen::MatrixXd SUMsub_M_Matrix = Eigen::MatrixXd::Zero(3, 3 * N);

	for (const auto& vertexEntry : verticesMap) {
		Vertex* vertex = vertexEntry.second;
		int localIdx = vertex->localIndex;  // Use local index

		// Create a 3x3 identity matrix scaled by the vertex's mass
		SUMsub_M_Matrix.block(0, 3 * localIdx, 3, 3) = (vertex->vertexMass / groupMass) * Eigen::Matrix3d::Identity();
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
		double mass = massMatrix(3 * localIdx, 3 * localIdx);
		vertexPair.second->vertexMass = mass;
	}
}



void Group::calInitCOM() {
	double totalMass = 0.0;
	Eigen::Vector3d weightedSum(0.0, 0.0, 0.0);

	for (const auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		totalMass += vertex->vertexMass;
		weightedSum += vertex->vertexMass * Eigen::Vector3d(vertex->initx, vertex->inity, vertex->initz);
	}

	if (totalMass > 0.0) {
		initCOM = weightedSum / totalMass;
	}
	else {
		// Handle the case where totalMass is zero to avoid division by zero.
		// You can set centerofMass to a default value or handle it according to your requirements.
		initCOM = Eigen::Vector3d(0.0, 0.0, 0.0);
	}
}

void Group::calPrimeVec(int w) {
	// 确保groupVelocity已经初始化且设置为正确的尺寸
	//groupVelocity = Eigen::VectorXd::Zero(3 * verticesMap.size());
	primeVec = Eigen::VectorXd::Zero(3 * verticesMap.size());
	// 计算逆矩阵
	Eigen::MatrixXd inverseTerm = (massMatrix + dampingMatrix * timeStep).inverse();

	for (int i = 0; i < 3 * verticesMap.size(); i += 1) {
		gravity(i) = 0; // y方向上设置重力
	}

	// 初始化gravity向量
	if (w == 4) {
		for (int i = 0; i < 3 * verticesMap.size(); i += 3) {
			gravity(i) = -Gravity; // y方向上设置重力 右
		}
	}
	else if (w == 2) {
		for (int i = 1; i < 3 * verticesMap.size(); i += 3) {
			gravity(i) = Gravity; // y方向上设置重力 下
		}
	}
	else if (w == 1) {
		for (int i = 1; i < 3 * verticesMap.size(); i += 3) {
			gravity(i) = -Gravity; // y方向上设置重力 上
		}
	}
	else if (w == 3) {
		for (int i = 0; i < 3 * verticesMap.size(); i += 3) {
			gravity(i) = Gravity; // y方向上设置重力 左
		}
	}
	

	// 更新groupVelocity
	groupVelocity += gravity * timeStep;

	// 使用整个矩阵计算velocityUpdate
	Eigen::VectorXd velocityUpdate = inverseTerm * (massMatrix * groupVelocity) * timeStep;

	// 更新primeVec和顶点位置
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localPi = vertex->localIndex; // 使用局部索引

		// 获取当前顶点的速度更新部分
		Eigen::Vector3d currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);

		// 计算新的位置
		Eigen::Vector3d newPosition = Eigen::Vector3d(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;

		// 更新primeVec
		primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;
		//if (vertex->isFixed) {
		//	// 对于固定点，将位置设置为初始位置
		//	primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = Eigen::Vector3d(vertex->initx, vertex->inity, vertex->initz);

		//	// 保持顶点位置不变
		//	vertex->x = vertex->initx;
		//	vertex->y = vertex->inity;
		//	vertex->z = vertex->initz;
		//}
		//else {
		//	// 获取当前顶点的速度更新部分
		//	Eigen::Vector3d currentVelocityUpdate = velocityUpdate.segment<3>(3 * localPi);

		//	// 计算新的位置
		//	Eigen::Vector3d newPosition = Eigen::Vector3d(vertex->x, vertex->y, vertex->z) + currentVelocityUpdate;

		//	// 更新primeVec
		//	primeVec.segment<3>(3 * static_cast<Eigen::Index>(localPi)) = newPosition;

		//	// 更新顶点位置
		//	vertex->x = newPosition.x();
		//	vertex->y = newPosition.y();
		//	vertex->z = newPosition.z();
		//}
	}

}

void Group::calLHS() {
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3 * verticesMap.size(), 3 * verticesMap.size());
	Eigen::MatrixXd A;
	Eigen::MatrixXd B;
	Eigen::MatrixXd C;
	A = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK;
	B = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * massDistribution;
	/*massDampingSparseInv = (massMatrix + timeStep * dampingMatrix).inverse().sparseView();
	A = timeStep * timeStep * massDampingSparseInv * kSparse;
	B = timeStep * timeStep * massDampingSparseInv * kSparse * massDistributionSparse;*/
	FEMLHS = I + A - B;

}

void Group::calRHS() {
	Eigen::VectorXd A;
	Eigen::VectorXd B;
	Eigen::VectorXd C;
	Eigen::VectorXd D;
	//Fbind = Eigen::VectorXd::Zero(3 * verticesMap.size());
	A = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * initLocalPos;
	B = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * rotationMatrix.transpose() * primeVec;
	C = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * groupK * rotationMatrix.transpose() * massDistribution * primeVec;
	D = timeStep * timeStep * (massMatrix + timeStep * dampingMatrix).inverse() * rotationMatrix.inverse() * Fbind;
	/*A = timeStep * timeStep * massDampingSparseInv * kSparse * initLocalPos;
	B = timeStep * timeStep * massDampingSparseInv * kSparse * rotationTransSparse * primeVec;
	C = timeStep * timeStep * massDampingSparseInv * kSparse * rotationTransSparse * massDistributionSparse * primeVec;
	D = timeStep * timeStep * massDampingSparseInv * rotationTransSparse * Fbind;*/
	//RHS = A - B + C + D;
	FEMRHS = A - B + C +D;

}

void Group::calDeltaX() {
	
	//deltaX = LHS.inverse() * RHS;
	if (FEMLHS.determinant() != 0) {
		// 解线性方程Ax = b
		deltaX = FEMLHS.colPivHouseholderQr().solve(FEMRHS);
	}
	else {
		std::cout << "矩阵A是奇异的，无法解决方程组." << std::endl;
	}
	deltaX = rotationMatrix.transpose() * deltaX;
}

void Group::calculateCurrentPositions() {
	// 遍历所有顶点
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localidx;
		localidx = vertex->localIndex;

		// 获取primeVec中对应顶点的位置
		Eigen::Vector3d primePosition = primeVec.segment<3>(3 * localidx);

		// 获取deltaX中对应顶点的位移
		Eigen::Vector3d displacement = deltaX.segment<3>(3 * localidx);

		// 计算当前位置
		Eigen::Vector3d currentPos = primePosition + displacement;
		currentPosition.segment<3>(3 * localidx) = currentPos;

	}
}
//void Group::updateVertexPositions() {
//	for (auto& vertexPair : verticesMap) {
//		Vertex* vertex = vertexPair.second;
//
//		// 使用局部索引来获取正确的矩阵块和primeVec部分
//		int localIndex = vertex->localIndex;
//		Eigen::Matrix3d rotationBlock = rotationMatrix.block<3, 3>(3 * localIndex, 3 * localIndex);
//		Eigen::Vector3d positionInPrimeVec = primeVec.segment<3>(3 * localIndex);
//
//		if (vertex->isFixed) {
//			// 对于固定点，将位置设置为初始位置
//			vertex->x = vertex->initx;
//			vertex->y = vertex->inity;
//			vertex->z = vertex->initz;
//		}
//		else {
//			// 使用旋转矩阵块乘以primeVec中的位置
//			Eigen::Vector3d newPosition = rotationBlock * positionInPrimeVec;
//
//			// 更新顶点位置
//			vertex->x = newPosition.x();
//			vertex->y = newPosition.y();
//			vertex->z = newPosition.z();
//		}
//	}
//}

void Group::calFbind(const std::vector<Vertex*>& commonVerticesGroup1,
	const std::vector<Vertex*>& commonVerticesGroup2,
	const Eigen::VectorXd& currentPositionGroup1,
	const Eigen::VectorXd& currentPositionGroup2,
	double k) {
	// Initialize Fbind, with a length three times the number of vertices in the group
	Fbind = Eigen::VectorXd::Zero(verticesMap.size() * 3);
	Eigen::Vector3d posThisGroup;
	Eigen::Vector3d posOtherGroup;
	Eigen::Vector3d avgPosition;
	Eigen::Vector3d posDifference;
	Eigen::Vector3d force;
	posThisGroup = Eigen::Vector3d::Zero();
	posOtherGroup = Eigen::Vector3d::Zero();
	avgPosition = Eigen::Vector3d::Zero();
	posDifference = Eigen::Vector3d::Zero();
	force = Eigen::Vector3d::Zero();
	// Iterate through all common vertices
	for (size_t i = 0; i < commonVerticesGroup1.size(); ++i) {
		// Get the common vertex in this group and the other group
		Vertex* vertexThisGroup = commonVerticesGroup1[i];
		Vertex* vertexOtherGroup = commonVerticesGroup2[i];

		// Directly use x, y, z from the Vertex objects for position
		/*Eigen::Vector3d posThisGroup(vertexThisGroup->x, vertexThisGroup->y, vertexThisGroup->z);
		Eigen::Vector3d posOtherGroup(vertexOtherGroup->x, vertexOtherGroup->y, vertexOtherGroup->z);*/
		posThisGroup = currentPositionGroup1.segment<3>(3 * vertexThisGroup->localIndex);
		posThisGroup = currentPositionGroup2.segment<3>(3 * vertexOtherGroup->localIndex);
		avgPosition = (posThisGroup + posOtherGroup) / 2;
		// Compute the position difference between the current vertex and the other group's vertex
		posDifference = posThisGroup - avgPosition;
		// Compute the constraint force
		force = k * posDifference;
		// Place the constraint force in Fbind at the appropriate position using the local index
		Fbind.segment<3>(3 * vertexThisGroup->localIndex) = force;
	}
}



void Group::updatePosition() {
	Eigen::Vector3d pos = Eigen::Vector3d::Zero();
	// 遍历所有顶点
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localIndex = vertex->localIndex;

		// 从currentPosition中获取对应顶点的位置
		pos = currentPosition.segment<3>(3 * localIndex);
		/*vertex->x = pos.x();
		vertex->y = pos.y();
		vertex->z = pos.z();*/
		if (vertex->isFixed) {
			// 对于固定点，将位置设置为初始位置
			vertex->x = vertex->initx;
			vertex->y = vertex->inity;
			vertex->z = vertex->initz;
		}
		else {
			// 使用旋转矩阵块乘以primeVec中的位置
			

			vertex->x = pos.x();
			vertex->y = pos.y();
			vertex->z = pos.z();
		}

		// 更新顶点的位置
		
	}
}

void Group::updateVelocity() {
	Eigen::Vector3d previousPos = Eigen::Vector3d::Zero();
	Eigen::Vector3d currentPos = Eigen::Vector3d::Zero();
	Eigen::Vector3d velocity =  Eigen::Vector3d::Zero();
	
	// 遍历所有顶点，更新速度并保存当前位置
	for (auto& vertexPair : verticesMap) {
		Vertex* vertex = vertexPair.second;
		int localIndex = vertex->localIndex;

		// 获取当前位置
		previousPos.x() = vertex->x;
		previousPos.y() = vertex->y;
		previousPos.z() = vertex->z;

		// 从 previousPosition 获取上一帧的位置
		currentPos = currentPosition.segment<3>(3 * localIndex);

		// 计算速度
		velocity = (currentPos - previousPos) / timeStep;
		groupVelocity.segment<3>(3 * localIndex) = velocity;
		// 更新 vertex 的速度
		// 例如：vertex->velocity = velocity;

		// 保存当前位置作为下一帧的“上一帧位置”
		//previousPosition.segment<3>(3 * localIndex) = currentPos;
	}

	//// 更新currentPosition为本帧最后的位置
	//for (auto& vertexPair : verticesMap) {
	//	Vertex* vertex = vertexPair.second;
	//	int localIndex = vertex->localIndex;

	//	currentPosition.segment<3>(3 * localIndex) = Eigen::Vector3d(vertex->x, vertex->y, vertex->z);
	//}
}

Group& Object::getGroup(int index) {
	return groups[index];
}

void Object::PBDLOOP(int looptime) {


	// 1. 初始化：将每个组的 Fbind 置零
	for (auto& g : groups) {
		g.Fbind = Eigen::VectorXd::Zero(3 * g.verticesMap.size()); // 假设 Group 类有一个方法来清除 Fbind
	}

	// 2. 开始迭代
	for (int iter = 0; iter < looptime; ++iter) {
		// 每组计算 RHS
		for (auto& g : groups) {
			g.calRHS(); 
			g.calDeltaX();
			g.calculateCurrentPositions();
						
		}
		//groups[0].calFbind(commonPoints.first, commonPoints.second, groups[0].currentPosition, groups[1].currentPosition, - 10000);
		//groups[1].calFbind(commonPoints.second, commonPoints.first, groups[1].currentPosition,groups[0].currentPosition, -10000);

		//groups[1].calFbind(commonPoints1.first, commonPoints1.second, 1000);
		//groups[2].calFbind(commonPoints1.second, commonPoints1.first,1000);		
	}
	for (auto& g : groups)
	{
		g.updateVelocity();
		g.updatePosition();

	}
	// 迭代完成后更新位置和速度
	//for (int i = 0; i < 3; ++i) {
	//	// 更新位置，这里可能需要一些逻辑来获取最后一次迭代的结果
	//	groups[i].updateFinalPositions(); // 假设这个方法用最后一次迭代的结果更新顶点位置

	//	// 更新速度
	//	groups[i].updateVelocities(timestep); // 假设这个方法用 (现在位置 - 上一帧位置) / timestep 计算速度
	//}

	// ... 现在，所有的组都应该有了更新后的位置和速度，可以传递给绘图功能
	// drawGroups(); // 假设有一个方法来绘制或输出最新的组状态
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

		Tetrahedron* tet = new Tetrahedron(v1, v2, v3, v4); //把vertex打包进四面体


		// Determine group based on average X coordinate
		double avgX = (v1->x + v2->x + v3->x + v4->x) / 4;
		int groupIndex = (avgX - minX) / groupRange;
		if (groupIndex >= numGroups) groupIndex = numGroups - 1;  // Update here

		object.getGroup(groupIndex).addTetrahedron(tet); //把四面体打包进group

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