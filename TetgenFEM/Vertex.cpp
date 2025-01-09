#include "Vertex.h"
#include <iostream>

Vertex::Vertex(float x, float y, float z, int index)
    : initx(x), inity(y), initz(z), // Initialize const members first
    x(x), y(y), z(z), // Then initialize non-const members
    index(index), vertexMass(1), // Initialize other members
    velx(0), vely(0), velz(0), // Initialize velocity components to 0
    isFixed(false) // Default to not fixed
{}

void Vertex::setFixedIfBelowThreshold() {
    if (inity > 0.42) {
        isFixed = true;
    }
}

//bunnyfront ymin=-0.61
//cloth (initx < -0.6 && inity < -0.25) || (initx > 0.6 && inity < -0.25)
//armadillo ��r initx < -1.1 && inity > 0.63) || (initx > 1.1 && inity > 0.63) || inity < -0.42
//bunny x<-0.5