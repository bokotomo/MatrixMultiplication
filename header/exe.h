void calculationMatrixMultiplication(){
    long i, j, k, ib, jb, kb;
    int block = 16;
    double tA, tA2;

    #pragma omp parallel
    {
      #pragma omp for private(i, j, k, ib, jb, kb)
      for(ib = 0; ib < matrixSize; ib += block){
        for(j = 0; j < matrixSize; j += block){
          for(k = 0; k < matrixSize; k += block){
            for(i = ib; i < ib + block; i += 1){

              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + j] * B[matrixSize * j + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + k];
              C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + k];
              
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + j] * B[matrixSize * j + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+1)];
              C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+1)];
              
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + j] * B[matrixSize * j + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+2)];
              C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+2)];
              
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + j] * B[matrixSize * j + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+3)];
              C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+3)];
              
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + j] * B[matrixSize * j + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+4)];
              C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+4)];
              
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + j] * B[matrixSize * j + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+5)];
              C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+5)];
              
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + j] * B[matrixSize * j + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+6)];
              C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+6)];
              
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + j] * B[matrixSize * j + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+7)];
              C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+7)];
              
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + j] * B[matrixSize * j + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+8)];
              C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+8)];
              
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + j] * B[matrixSize * j + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+9)];
              C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+9)];
              
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + j] * B[matrixSize * j + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+10)];
              C[i * matrixSize + (k+10)] = C[i * matrixSize + (k+10)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+10)];
              
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + j] * B[matrixSize * j + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+11)];
              C[i * matrixSize + (k+11)] = C[i * matrixSize + (k+11)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+11)];
              
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + j] * B[matrixSize * j + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+12)];
              C[i * matrixSize + (k+12)] = C[i * matrixSize + (k+12)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+12)];
              
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + j] * B[matrixSize * j + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+13)];
              C[i * matrixSize + (k+13)] = C[i * matrixSize + (k+13)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+13)];
              
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + j] * B[matrixSize * j + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+14)];
              C[i * matrixSize + (k+14)] = C[i * matrixSize + (k+14)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+14)];
              
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + j] * B[matrixSize * j + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+1)] * B[matrixSize * (j+1) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+2)] * B[matrixSize * (j+2) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+3)] * B[matrixSize * (j+3) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+4)] * B[matrixSize * (j+4) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+5)] * B[matrixSize * (j+5) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+6)] * B[matrixSize * (j+6) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+7)] * B[matrixSize * (j+7) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+8)] * B[matrixSize * (j+8) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+9)] * B[matrixSize * (j+9) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+10)] * B[matrixSize * (j+10) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+11)] * B[matrixSize * (j+11) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+12)] * B[matrixSize * (j+12) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+13)] * B[matrixSize * (j+13) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+14)] * B[matrixSize * (j+14) + (k+15)];
              C[i * matrixSize + (k+15)] = C[i * matrixSize + (k+15)] + A[i * matrixSize + (j+15)] * B[matrixSize * (j+15) + (k+15)];

            }
          }
        }
      }
    }
  }