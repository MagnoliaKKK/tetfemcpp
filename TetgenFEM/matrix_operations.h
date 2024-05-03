#pragma once
#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

// ÉùÃ÷ CUDA ¾ØÕó³Ë·¨º¯Êý
void cudaMatrixMultiply(double* a, double* b, double* c, int width, int common, int height);

#endif // MATRIX_OPERATIONS_H
