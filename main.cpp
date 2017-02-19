#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <immintrin.h>
#include <omp.h>
#include "matrix_multiplication.h"

using namespace std;

int main(int argc, char *argv[]){
  long size = 1024;
  int threadsNum;
  double executionTime, GFlops;
  
  threadsNum = atoi(argv[1]);

  MatrixMultiplication matMul(size, threadsNum);
  executionTime = matMul.measurement();

  cout << "Time   is " << fixed << executionTime << endl;

  GFlops = ( ((double)size/1000 * (double)size/1000 * (double)size/1000) * 2 ) / executionTime;
  cout << "GFlops is " << fixed << GFlops << endl;
  
  return 0;
}
