#include "Edge.h"

Edge::Edge(Vertex* v1, Vertex* v2) : isBoundary(false) {
    vertices[0] = v1;
    vertices[1] = v2;
}
