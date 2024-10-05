#pragma once
#ifndef EDGE_H
#define EDGE_H

#include "Vertex.h"

class Edge {
public:
    Vertex* vertices[2];
    bool isBoundary;

    Edge(Vertex* v1, Vertex* v2);
};

#endif // EDGE_H
#pragma once
