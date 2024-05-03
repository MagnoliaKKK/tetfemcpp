#include <cuda_runtime.h>
#include <iostream>

__global__ void matrixMultiplyKernel(double* a, double* b, double* c, int width, int common, int height) {
    // 计算行和列
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < height && col < width) {
        double sum = 0;
        for (int i = 0; i < common; i++) {
            sum += a[row * common + i] * b[i * width + col];
        }
        c[row * width + col] = sum;
    }
}

void cudaMatrixMultiply(double* a, double* b, double* c, int width, int common, int height) {
    double* dev_a, * dev_b, * dev_c;

    // 分配 GPU 内存
    cudaMalloc((void**)&dev_a, width * common * sizeof(double));
    cudaMalloc((void**)&dev_b, common * height * sizeof(double));
    cudaMalloc((void**)&dev_c, width * height * sizeof(double));

    // 复制数据到 GPU
    cudaMemcpy(dev_a, a, width * common * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, common * height * sizeof(double), cudaMemcpyHostToDevice);

    // 设置 CUDA 核函数的维度
    dim3 dimGrid(ceil(width / 16.0), ceil(height / 16.0), 1);
    dim3 dimBlock(16, 16, 1);

    // 启动 CUDA 核函数
    matrixMultiplyKernel << <dimGrid, dimBlock >> > (dev_a, dev_b, dev_c, width, common, height);

    // 将结果复制回主机
    cudaMemcpy(c, dev_c, width * height * sizeof(double), cudaMemcpyDeviceToHost);

    // 释放 GPU 内存
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);
}

