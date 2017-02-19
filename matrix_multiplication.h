using namespace std;

class MatrixMultiplication
{
private:
  double *A, *B, *C, *tmpC;
  int matrixSize;

  void setMatrixData(){
    int i;

    srand( time(NULL) );

    #pragma omp parallel
    {
      #pragma omp for
      for(i = 0 ; i < matrixSize * matrixSize ; i++){
        A[i] = rand () % 5 + 0.1;
        B[i] = rand () % 5 + 0.2;
      }
    }
  }

  void setMatrixDataC(){
    int i;

    #pragma omp parallel
    {
      #pragma omp for
      for(i = 0 ; i < matrixSize * matrixSize ; i++){
        C[i] = 0;
      }
    }
  }

  bool checkMatrixData(){
    bool checkFlag = false;
    int i;

    for(i = 0 ; i < matrixSize * matrixSize ; i++){
      if(fabs(C[i] - tmpC[i]) > 0.001){
        checkFlag = true;
        break;
      }
    }

    return checkFlag;
  }

  void showMatrixData(){
    int i, j;
    int num = matrixSize;

    for(i = 0 ; i < num  ; i++){
      for(j = 0 ; j < num ; j++){
        cout << A[i * num + j] << " ";
      }
      cout << endl;
    }
    
    cout << endl;
    for(i = 0 ; i < num  ; i++){
      for(j = 0 ; j < num ; j++){
        cout << B[i * num + j] << " ";
      }
      cout << endl;
    }
    
    cout << endl;
    for(i = 0 ; i < num  ; i++){
      for(j = 0 ; j < num ; j++){
        cout << C[i * num + j] << " ";
      }
      cout << endl;
    }
    
    cout << endl;
    for(i = 0 ; i < num  ; i++){
      for(j = 0 ; j < num ; j++){
        cout << tmpC[i * num + j] << " ";
      }
      cout << endl;
    }
  }

  void calculationMatrixMultiplication(){
    int i, j, k;
    double tA;

    #pragma omp parallel
    {
      #pragma omp for private(i, j, k, tA)
      for(i = 0; i < matrixSize; i++){
      for(j = 0; j < matrixSize; j++){
      for(k = 0; k < matrixSize; k++){
        C[i * matrixSize + k] = C[i * matrixSize + k] + A[i * matrixSize + j] * B[matrixSize * j + k];
      }}}
    }
  }

  void calculationMatrixMultiplicationSimple(){
    int i, j, k;

    #pragma omp parallel
    {
      #pragma omp for private(i, j, k)
      for(i = 0; i < matrixSize; i++){
      for(j = 0; j < matrixSize; j++){
      for(k = 0; k < matrixSize; k++){
        tmpC[i * matrixSize + k] = tmpC[i * matrixSize + k] + A[i * matrixSize + j] * B[matrixSize * j + k];
      }}}
    }
  }

  void calculationMatrixMultiplicationSIMD(){
    int i, j, k;
  	__m256d rA1, rB1, rA2, rB2, rA3, rB3, rA4, rB4, rA5, rB5, rA6, rB6, rA7, rB7, rA8, rB8, rA21, rA22, rA23, rA24, rA25, rA26, rA27, rA28, rA31, rA32, rA33, rA34, rA35, rA36, rA37, rA38, rA41, rA42, rA43, rA44, rA45, rA46, rA47, rA48, rSUM1, rSUM2, rSUM3, rSUM4;
  	double d[4] = {0};
     
    #pragma omp parallel
    {
      #pragma omp for private(i, j, k, rA1, rB1, rA2, rB2, rA3, rB3, rA4, rB4, rA5, rB5, rA6, rB6, rA7, rB7, rA8, rB8, rA21, rA22, rA23, rA24, rA25, rA26, rA27, rA28, rA31, rA32, rA33, rA34, rA35, rA36, rA37, rA38, rA41, rA42, rA43, rA44, rA45, rA46, rA47, rA48, rSUM1, rSUM2, rSUM3, rSUM4, d)
      for(i = 0; i < matrixSize; i++){
        for(k = 0; k < matrixSize; k+=4){
          rSUM1 = _mm256_setzero_pd();
          rSUM2 = _mm256_setzero_pd();
          rSUM3 = _mm256_setzero_pd();
          rSUM4 = _mm256_setzero_pd();

          for(j = 0; j < matrixSize; j+=32){
            rA1 = _mm256_load_pd(&A[k * matrixSize + j]);
            rA2 = _mm256_load_pd(&A[k * matrixSize + j + 4]);
            rA3 = _mm256_load_pd(&A[k * matrixSize + j + 8]);
            rA4 = _mm256_load_pd(&A[k * matrixSize + j + 12]);
            rA5 = _mm256_load_pd(&A[k * matrixSize + j + 16]);
            rA6 = _mm256_load_pd(&A[k * matrixSize + j + 20]);
            rA7 = _mm256_load_pd(&A[k * matrixSize + j + 24]);
            rA8 = _mm256_load_pd(&A[k * matrixSize + j + 28]);
            
            rB1 = _mm256_load_pd(&B[i * matrixSize + j]);
            rB2 = _mm256_load_pd(&B[i * matrixSize + j + 4]);
            rB3 = _mm256_load_pd(&B[i * matrixSize + j + 8]);
            rB4 = _mm256_load_pd(&B[i * matrixSize + j + 12]);
            rB5 = _mm256_load_pd(&B[i * matrixSize + j + 16]);
            rB6 = _mm256_load_pd(&B[i * matrixSize + j + 20]);
            rB7 = _mm256_load_pd(&B[i * matrixSize + j + 24]);
            rB8 = _mm256_load_pd(&B[i * matrixSize + j + 28]);
            
            rA21 = _mm256_load_pd(&A[(k+1) * matrixSize + j]);
            rA22 = _mm256_load_pd(&A[(k+1) * matrixSize + j + 4]);
            rA23 = _mm256_load_pd(&A[(k+1) * matrixSize + j + 8]);
            rA24 = _mm256_load_pd(&A[(k+1) * matrixSize + j + 12]);
            rA25 = _mm256_load_pd(&A[(k+1) * matrixSize + j + 16]);
            rA26 = _mm256_load_pd(&A[(k+1) * matrixSize + j + 20]);
            rA27 = _mm256_load_pd(&A[(k+1) * matrixSize + j + 24]);
            rA28 = _mm256_load_pd(&A[(k+1) * matrixSize + j + 28]);
            
            rA31 = _mm256_load_pd(&A[(k+2) * matrixSize + j]);
            rA32 = _mm256_load_pd(&A[(k+2) * matrixSize + j + 4]);
            rA33 = _mm256_load_pd(&A[(k+2) * matrixSize + j + 8]);
            rA34 = _mm256_load_pd(&A[(k+2) * matrixSize + j + 12]);
            rA35 = _mm256_load_pd(&A[(k+2) * matrixSize + j + 16]);
            rA36 = _mm256_load_pd(&A[(k+2) * matrixSize + j + 20]);
            rA37 = _mm256_load_pd(&A[(k+2) * matrixSize + j + 24]);
            rA38 = _mm256_load_pd(&A[(k+2) * matrixSize + j + 28]);
            
            rA41 = _mm256_load_pd(&A[(k+3) * matrixSize + j]);
            rA42 = _mm256_load_pd(&A[(k+3) * matrixSize + j + 4]);
            rA43 = _mm256_load_pd(&A[(k+3) * matrixSize + j + 8]);
            rA44 = _mm256_load_pd(&A[(k+3) * matrixSize + j + 12]);
            rA45 = _mm256_load_pd(&A[(k+3) * matrixSize + j + 16]);
            rA46 = _mm256_load_pd(&A[(k+3) * matrixSize + j + 20]);
            rA47 = _mm256_load_pd(&A[(k+3) * matrixSize + j + 24]);
            rA48 = _mm256_load_pd(&A[(k+3) * matrixSize + j + 28]);
            
            rSUM1 = _mm256_fmadd_pd(rA1, rB1, rSUM1);
            rSUM2 = _mm256_fmadd_pd(rA21, rB1, rSUM2);
            rSUM3 = _mm256_fmadd_pd(rA31, rB1, rSUM3);
            rSUM4 = _mm256_fmadd_pd(rA41, rB1, rSUM4);
            rSUM1 = _mm256_fmadd_pd(rA2, rB2, rSUM1);
            rSUM2 = _mm256_fmadd_pd(rA22, rB2, rSUM2);
            rSUM3 = _mm256_fmadd_pd(rA32, rB2, rSUM3);
            rSUM4 = _mm256_fmadd_pd(rA42, rB2, rSUM4);
            rSUM1 = _mm256_fmadd_pd(rA3, rB3, rSUM1);
            rSUM2 = _mm256_fmadd_pd(rA23, rB3, rSUM2);
            rSUM3 = _mm256_fmadd_pd(rA33, rB3, rSUM3);
            rSUM4 = _mm256_fmadd_pd(rA43, rB3, rSUM4);
            rSUM1 = _mm256_fmadd_pd(rA4, rB4, rSUM1);
            rSUM2 = _mm256_fmadd_pd(rA24, rB4, rSUM2);
            rSUM3 = _mm256_fmadd_pd(rA34, rB4, rSUM3);
            rSUM4 = _mm256_fmadd_pd(rA44, rB4, rSUM4);
            rSUM1 = _mm256_fmadd_pd(rA5, rB5, rSUM1);
            rSUM2 = _mm256_fmadd_pd(rA25, rB5, rSUM2);
            rSUM3 = _mm256_fmadd_pd(rA35, rB5, rSUM3);
            rSUM4 = _mm256_fmadd_pd(rA45, rB5, rSUM4);
            rSUM1 = _mm256_fmadd_pd(rA6, rB6, rSUM1);
            rSUM2 = _mm256_fmadd_pd(rA26, rB6, rSUM2);
            rSUM3 = _mm256_fmadd_pd(rA36, rB6, rSUM3);
            rSUM4 = _mm256_fmadd_pd(rA46, rB6, rSUM4);
            rSUM1 = _mm256_fmadd_pd(rA7, rB7, rSUM1);
            rSUM2 = _mm256_fmadd_pd(rA27, rB7, rSUM2);
            rSUM3 = _mm256_fmadd_pd(rA37, rB7, rSUM3);
            rSUM4 = _mm256_fmadd_pd(rA47, rB7, rSUM4);
            rSUM1 = _mm256_fmadd_pd(rA8, rB8, rSUM1);
            rSUM2 = _mm256_fmadd_pd(rA28, rB8, rSUM2);
            rSUM3 = _mm256_fmadd_pd(rA38, rB8, rSUM3);
            rSUM4 = _mm256_fmadd_pd(rA48, rB8, rSUM4);

          }

          _mm256_store_pd(d, _mm256_hadd_pd(rSUM1, rSUM2));
          C[k * matrixSize + i] = d[0] + d[2];
          C[(k+1) * matrixSize + i] = d[1] + d[3];
          _mm256_store_pd(d, _mm256_hadd_pd(rSUM3, rSUM4));
          C[(k+2) * matrixSize + i] = d[0] + d[2];
          C[(k+3) * matrixSize + i] = d[1] + d[3];
        }
      }
    }
  }

  void transpositionMatrix(){
    int i, j;
    double *tmp;
    
    tmp = (double *)_mm_malloc(sizeof(double) * matrixSize * matrixSize, 32);
    #pragma omp parallel
    {
      #pragma omp for private(i)
      for(i = 0; i < matrixSize * matrixSize; i+=4){
        tmp[i] = B[i];
        tmp[i+1] = B[i+1];
        tmp[i+2] = B[i+2];
        tmp[i+3] = B[i+3];
      }
    }
    
    #pragma omp parallel
    {
      #pragma omp for private(i, j)
      for(i = 0; i < matrixSize; i++){
        for(j = 0; j < matrixSize; j+=4){
          B[i* matrixSize + j] = tmp[j* matrixSize + i];
          B[i* matrixSize + (j+1)] = tmp[(j+1)* matrixSize + i];
          B[i* matrixSize + (j+2)] = tmp[(j+2)* matrixSize + i];
          B[i* matrixSize + (j+3)] = tmp[(j+3)* matrixSize + i];
        }
      }
    }
    _mm_free(tmp);
  }

  void setThreadsNum(int threadsNum){
    omp_set_num_threads(threadsNum);
  
    #pragma omp parallel num_threads(threadsNum)
    {
      #pragma omp single
      {
        cout << "ThreadsNum = " << omp_get_num_threads() << endl;
      }
    }
  }
public:
  MatrixMultiplication(int size, int threadsNum){
    cout << "Matrix Multiplication" << endl;

    setThreadsNum(threadsNum);

    matrixSize = size;
    A = (double *)_mm_malloc(sizeof(double) * size * size, 32);
    B = (double *)_mm_malloc(sizeof(double) * size * size, 32);
    C = (double *)_mm_malloc(sizeof(double) * size * size, 32);
    tmpC = (double *)_mm_malloc(sizeof(double) * size * size, 32);
    
    setMatrixData();

    calculationMatrixMultiplicationSimple();

    cout << "------ MatrixData Setting complete" << endl;
  }

  double measurement(){
    double startTime, endTime, resultTime = 0.0;
    int i;

    transpositionMatrix();
    calculationMatrixMultiplicationSIMD();
    
    for(i = 0; i < 5; i++){
      setMatrixDataC();
      startTime = omp_get_wtime();
      calculationMatrixMultiplicationSIMD();
      endTime = omp_get_wtime();
      resultTime += (double)(endTime - startTime);
    }
    resultTime = resultTime / 5;

    //showMatrixData();

    if( checkMatrixData() == true ){
      cout << "MatrixData Error!!!" << endl;
      exit(1);
    }

    _mm_free(A);
    _mm_free(B);
    _mm_free(C);
    _mm_free(tmpC);

    return resultTime;
  }
};
