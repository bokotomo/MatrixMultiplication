using namespace std;

class MatrixMultiplication
{
private:
  double *A, *B, *C, *tmpC;
  long matrixSize;
  void setMatrixData(){
    long i;

    srand( time(NULL) );

    #pragma omp parallel
    {
      #pragma omp for
      for(i = 0 ; i < matrixSize * matrixSize ; i++){
        A[i] = i+1;//rand () % 100 + 0.1;
        B[i] = i+2;//rand () % 100 + 0.1;
      }
    }
  }
  void setMatrixDataC(){
    long i;

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
    long i;

    for(i = 0 ; i < matrixSize * matrixSize ; i++){
      if(C[i] != tmpC[i]){
        checkFlag = true;
        break;
      }
    }

    return checkFlag;
  }
  void showMatrixData(){
    long i, j;
    long num = matrixSize;

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
    long i, j, k;
    double tA;

    #pragma omp parallel
    {
      #pragma omp for private(i, j, k, tA)
      for(i = 0; i < matrixSize; i++){
        for(j = 0; j < matrixSize; j++){
          tA = A[i * matrixSize + j];
          for(k = 0; k < matrixSize; k+=8){
            C[i * matrixSize + k] = C[i * matrixSize + k] + tA * B[matrixSize * j + k];
            C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + tA * B[matrixSize * j + (k+1)];
            C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + tA * B[matrixSize * j + (k+2)];
            C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + tA * B[matrixSize * j + (k+3)];
            C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + tA * B[matrixSize * j + (k+4)];
            C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + tA * B[matrixSize * j + (k+5)];
            C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + tA * B[matrixSize * j + (k+6)];
            C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + tA * B[matrixSize * j + (k+7)];
          }
        }
      }
    }
  }
  void calculationMatrixMultiplicationSimple(){
    long i, j, k;

    #pragma omp parallel
    {
      #pragma omp for private(i, j, k)
      for(i = 0; i < matrixSize; i++){
        for(j = 0; j < matrixSize; j++){
          for(k = 0; k < matrixSize; k++){
            tmpC[i * matrixSize + k] = tmpC[i * matrixSize + k] + A[i * matrixSize + j] * B[matrixSize * j + k];
          }
        }
      }
    }
  }
  void calculationMatrixMultiplicationSIMD(){
    long i, j, k;
  	__m256d reg1, regA, regB, regC, regSUM;
  	double d[4] = {0};
     
    #pragma omp parallel
    {
      #pragma omp for private(i, j, k, regA, regB, regSUM, d)
      for(i = 0; i < matrixSize; i++){
        for(k = 0; k < matrixSize; k++){
          regSUM = _mm256_setzero_pd();
          for(j = 0; j < matrixSize; j+=4){
            
  					regA = _mm256_load_pd(&A[k * matrixSize + j]);
  					regB = _mm256_load_pd(&B[i * matrixSize + j]);
            regSUM = _mm256_fmadd_pd(regA, regB, regSUM);
            
          }
          _mm256_store_pd(d, regSUM);
          C[k * matrixSize + i] += d[0]+d[1]+d[2]+d[3];
        }
      }
    }
  }
  void transpositionMatrix(){
    long i, j;
    double *tmp;
    
    tmp = (double *)_mm_malloc(sizeof(double) * matrixSize * matrixSize, 32);
    for(i = 0; i < matrixSize * matrixSize; i++){
      tmp[i] = B[i];
    }
    for(i = 0; i < matrixSize; i++){
      for(j = 0; j < matrixSize; j++){
        B[i* matrixSize + j] = tmp[j* matrixSize + i];
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
    tmpC = new double[size * size];

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
    return resultTime;
  }
};