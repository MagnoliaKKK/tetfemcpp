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
    if (initx < -0.6/*initx < -0.64*/ /* ||  || initx > -0.15 || initx > 0.62 */) {
        isFixed = true;
    }
}

