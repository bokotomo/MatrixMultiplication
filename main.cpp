#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <cmath> 
#include <immintrin.h>
#include <omp.h>
#include "matrix_multiplication.h"

using namespace std;

int main(int argc, char *argv[]){
  int size = 1024 * 1;
  int threadsNum;
  double executionTime, GFlops;

  threadsNum = atoi(argv[1]);

  MatrixMultiplication matMul(size, threadsNum);
  executionTime = matMul.measurement();
  //matMul.showMatrixData();
  matMul.memoryFree();

  cout << "Time   is " << fixed << executionTime << endl;

  GFlops = ( ((double)size/1000 * (double)size/1000 * (double)size/1000) * 2 ) / executionTime;
  cout << "GFlops is " << fixed << GFlops << endl;

  return 0;
}
