#include <iostream>
#include <vector>
#include <Eigen/Core>

// Simple test for collision detection logic
struct SimpleVertex {
    float x, y, z;
    float velx, vely, velz;
    bool isFixed;
    
    SimpleVertex(float x_, float y_, float z_) 
        : x(x_), y(y_), z(z_), velx(0), vely(0), velz(0), isFixed(false) {}
};

// Simple ground plane
struct SimpleGround {
    Eigen::Vector3f position;
    Eigen::Vector3f normal;
    float restitution;
    
    SimpleGround() : position(0, -2.0f, 0), normal(0, 1.0f, 0), restitution(0.3f) {}
    
    float distanceToPoint(const Eigen::Vector3f& point) const {
        return normal.dot(point - position);
    }
};

// Simple collision test
void testCollisionLogic() {
    std::cout << "=== Collision Detection Test ===" << std::endl;
    
    // Create test vertex above ground
    SimpleVertex vertex(0, 0, 0);
    vertex.vely = -5.0f; // Moving downward
    
    // Create ground
    SimpleGround ground;
    
    std::cout << "Initial vertex position: (" << vertex.x << ", " << vertex.y << ", " << vertex.z << ")" << std::endl;
    std::cout << "Initial vertex velocity: (" << vertex.velx << ", " << vertex.vely << ", " << vertex.velz << ")" << std::endl;
    
    // Simulate falling
    float deltaTime = 0.016f; // ~60 FPS
    float gravity = -9.81f;
    
    for (int frame = 0; frame < 150; ++frame) {
        // Apply gravity
        vertex.vely += gravity * deltaTime;
        
        // Update position
        vertex.x += vertex.velx * deltaTime;
        vertex.y += vertex.vely * deltaTime;
        vertex.z += vertex.velz * deltaTime;
        
        // Check collision with ground
        Eigen::Vector3f pos(vertex.x, vertex.y, vertex.z);
        float distance = ground.distanceToPoint(pos);
        
        if (distance <= 0.0f) {
            std::cout << "COLLISION at frame " << frame << "!" << std::endl;
            std::cout << "Position: (" << vertex.x << ", " << vertex.y << ", " << vertex.z << ")" << std::endl;
            std::cout << "Velocity before: (" << vertex.velx << ", " << vertex.vely << ", " << vertex.velz << ")" << std::endl;
            
            // Project to ground plane
            vertex.y = ground.position.y(); // Simple projection for horizontal ground
            
            // Apply restitution
            vertex.vely = -vertex.vely * ground.restitution;
            
            std::cout << "Velocity after: (" << vertex.velx << ", " << vertex.vely << ", " << vertex.velz << ")" << std::endl;
            std::cout << "---" << std::endl;
            
            // Stop if velocity is too small
            if (std::abs(vertex.vely) < 0.1f) {
                std::cout << "Vertex settled on ground." << std::endl;
                break;
            }
        }
        
        // Print every 30 frames
        if (frame % 30 == 0) {
            std::cout << "Frame " << frame << ": pos(" << vertex.x << ", " << vertex.y << ", " << vertex.z 
                      << "), vel(" << vertex.velx << ", " << vertex.vely << ", " << vertex.velz << ")" << std::endl;
        }
    }
    
    std::cout << "Final position: (" << vertex.x << ", " << vertex.y << ", " << vertex.z << ")" << std::endl;
    std::cout << "Final velocity: (" << vertex.velx << ", " << vertex.vely << ", " << vertex.velz << ")" << std::endl;
    std::cout << "=== Test Complete ===" << std::endl;
}

int main() {
    testCollisionLogic();
    return 0;
}