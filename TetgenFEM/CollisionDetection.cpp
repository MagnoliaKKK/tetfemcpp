#include "CollisionDetection.h"
#include "Group.h"
#include <iostream>
#include <algorithm>

void CollisionDetection::applyGravity(Object& object, float deltaTime) {
    for (int groupIdx = 0; groupIdx < object.groupNum; ++groupIdx) {
        Group& group = object.groups[groupIdx];
        
        for (const auto& vertexPair : group.verticesMap) {
            Vertex* vertex = vertexPair.second;
            
            // Skip fixed vertices
            if (vertex->isFixed) {
                continue;
            }
            
            // Apply gravity acceleration to velocity
            vertex->vely += gravity * deltaTime;
            
            // Apply damping to all velocity components
            vertex->velx *= damping;
            vertex->vely *= damping;
            vertex->velz *= damping;
        }
    }
}

void CollisionDetection::detectCollisions(Object& object) {
    clearConstraints();
    
    for (int groupIdx = 0; groupIdx < object.groupNum; ++groupIdx) {
        Group& group = object.groups[groupIdx];
        
        for (const auto& vertexPair : group.verticesMap) {
            Vertex* vertex = vertexPair.second;
            
            // Skip fixed vertices
            if (vertex->isFixed) {
                continue;
            }
            
            CollisionConstraint constraint;
            if (checkVertexGroundCollision(vertex, constraint)) {
                constraints.push_back(constraint);
            }
        }
    }
}

bool CollisionDetection::checkVertexGroundCollision(Vertex* vertex, CollisionConstraint& constraint) {
    Eigen::Vector3f vertexPos(vertex->x, vertex->y, vertex->z);
    
    for (const auto& ground : grounds) {
        float distance = ground.distanceToPoint(vertexPos);
        
        // Check if vertex is below or touching the ground
        if (distance <= 0.0f) {
            // Project vertex to ground plane
            Eigen::Vector3f projection = projectToGround(vertexPos, ground);
            float penetration = -distance; // Positive penetration depth
            
            constraint = CollisionConstraint(vertex, projection, ground.normal, penetration);
            return true;
        }
    }
    
    return false;
}

Eigen::Vector3f CollisionDetection::projectToGround(const Eigen::Vector3f& point, const Ground& ground) {
    float distance = ground.distanceToPoint(point);
    return point - distance * ground.normal;
}

void CollisionDetection::solveConstraints(Object& object, int iterations) {
    for (int iter = 0; iter < iterations; ++iter) {
        for (auto& constraint : constraints) {
            if (!constraint.isActive) continue;
            
            Vertex* vertex = constraint.vertex;
            
            // Calculate correction vector
            Eigen::Vector3f currentPos(vertex->x, vertex->y, vertex->z);
            Eigen::Vector3f correction = constraint.projection - currentPos;
            
            // Apply positional correction with stiffness factor
            float stiffness = 1.0f; // Full stiffness for hard constraints
            correction *= stiffness;
            
            // Update vertex position
            vertex->x += correction.x();
            vertex->y += correction.y();
            vertex->z += correction.z();
        }
    }
}

void CollisionDetection::applyCollisionResponse(Object& object, float deltaTime) {
    for (auto& constraint : constraints) {
        if (!constraint.isActive) continue;
        
        Vertex* vertex = constraint.vertex;
        Eigen::Vector3f velocity(vertex->velx, vertex->vely, vertex->velz);
        
        // Find the ground involved in collision
        const Ground* collisionGround = nullptr;
        for (const auto& ground : grounds) {
            Eigen::Vector3f vertexPos(vertex->x, vertex->y, vertex->z);
            if (ground.distanceToPoint(vertexPos) <= 0.01f) { // Small threshold
                collisionGround = &ground;
                break;
            }
        }
        
        if (!collisionGround) continue;
        
        // Calculate velocity components relative to surface normal
        float normalVelocity = velocity.dot(constraint.normal);
        Eigen::Vector3f tangentialVelocity = velocity - normalVelocity * constraint.normal;
        
        // Apply restitution to normal component
        if (normalVelocity < 0) { // Moving towards surface
            normalVelocity = -normalVelocity * collisionGround->restitution;
        }
        
        // Apply friction to tangential component
        float tangentialSpeed = tangentialVelocity.norm();
        if (tangentialSpeed > 0.001f) { // Avoid division by zero
            float frictionForce = collisionGround->friction * std::abs(normalVelocity);
            float speedReduction = std::min(tangentialSpeed, frictionForce * deltaTime);
            tangentialVelocity *= (tangentialSpeed - speedReduction) / tangentialSpeed;
        }
        
        // Reconstruct velocity
        Eigen::Vector3f newVelocity = normalVelocity * constraint.normal + tangentialVelocity;
        
        // Update vertex velocity
        vertex->velx = newVelocity.x();
        vertex->vely = newVelocity.y();
        vertex->velz = newVelocity.z();
    }
}