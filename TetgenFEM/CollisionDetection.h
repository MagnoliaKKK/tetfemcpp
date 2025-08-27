#pragma once
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "Vertex.h"
#include "Object.h"

// Ground plane representation
class Ground {
public:
    Eigen::Vector3f position;   // Position on the plane
    Eigen::Vector3f normal;     // Normal vector of the plane
    float restitution;          // Coefficient of restitution (bounce factor)
    float friction;             // Friction coefficient
    
    Ground(const Eigen::Vector3f& pos = Eigen::Vector3f(0, -2.0f, 0), 
           const Eigen::Vector3f& norm = Eigen::Vector3f(0, 1.0f, 0),
           float rest = 0.3f, float fric = 0.5f)
        : position(pos), normal(norm.normalized()), restitution(rest), friction(fric) {}
    
    // Get distance from point to plane (negative if below plane)
    float distanceToPoint(const Eigen::Vector3f& point) const {
        return normal.dot(point - position);
    }
};

// Collision constraint for PBD system
struct CollisionConstraint {
    Vertex* vertex;             // Vertex involved in collision
    Eigen::Vector3f projection; // Projected position to satisfy constraint
    Eigen::Vector3f normal;     // Collision normal
    float penetration;          // Penetration depth
    bool isActive;              // Whether constraint is active
    
    CollisionConstraint() : vertex(nullptr), penetration(0.0f), isActive(false) {}
    
    CollisionConstraint(Vertex* v, const Eigen::Vector3f& proj, 
                       const Eigen::Vector3f& norm, float pen)
        : vertex(v), projection(proj), normal(norm), penetration(pen), isActive(true) {}
};

// Main collision detection class
class CollisionDetection {
public:
    std::vector<Ground> grounds;
    std::vector<CollisionConstraint> constraints;
    float gravity;              // Gravity acceleration
    float damping;              // Velocity damping factor
    
    CollisionDetection(float g = -9.81f, float d = 0.99f) : gravity(g), damping(d) {
        // Add default ground plane
        grounds.emplace_back();
    }
    
    // Add a ground plane
    void addGround(const Ground& ground) {
        grounds.push_back(ground);
    }
    
    // Apply gravity to all vertices
    void applyGravity(Object& object, float deltaTime);
    
    // Detect collisions and generate constraints
    void detectCollisions(Object& object);
    
    // Solve collision constraints using PBD approach
    void solveConstraints(Object& object, int iterations = 3);
    
    // Apply collision response (velocity updates)
    void applyCollisionResponse(Object& object, float deltaTime);
    
    // Clear all constraints
    void clearConstraints() {
        constraints.clear();
    }
    
    // Check if a vertex collides with any ground
    bool checkVertexGroundCollision(Vertex* vertex, CollisionConstraint& constraint);
    
    // Project vertex to nearest point on ground plane
    Eigen::Vector3f projectToGround(const Eigen::Vector3f& point, const Ground& ground);
    
    // Get ground at index (for rendering)
    const std::vector<Ground>& getGrounds() const { return grounds; }
};