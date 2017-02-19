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
  void showMatrixData(long num){
    long i, j;

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
    double tA, tA2;

    #pragma omp parallel
    {
      #pragma omp for private(i, j, k, tA, tA2)
      for(i = 0; i < matrixSize; i++){
        for(j = 0; j < matrixSize; j+=2){
          tA = A[i * matrixSize + j];
          tA2 = A[i * matrixSize + (j+1)];
          for(k = 0; k < matrixSize; k+=10){

            C[i * matrixSize + k] = C[i * matrixSize + k] + tA * B[matrixSize * j + k];
            C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + tA * B[matrixSize * j + (k+1)];
            C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + tA * B[matrixSize * j + (k+2)];
            C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + tA * B[matrixSize * j + (k+3)];
            C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + tA * B[matrixSize * j + (k+4)];
            C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + tA * B[matrixSize * j + (k+5)];
            C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + tA * B[matrixSize * j + (k+6)];
            C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + tA * B[matrixSize * j + (k+7)];
            C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + tA * B[matrixSize * j + (k+8)];
            C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + tA * B[matrixSize * j + (k+9)];
            
            C[i * matrixSize + k] = C[i * matrixSize + k] + tA2 * B[matrixSize * (j+1) + k];
            C[i * matrixSize + (k+1)] = C[i * matrixSize + (k+1)] + tA2 * B[matrixSize * (j+1) + (k+1)];
            C[i * matrixSize + (k+2)] = C[i * matrixSize + (k+2)] + tA2 * B[matrixSize * (j+1) + (k+2)];
            C[i * matrixSize + (k+3)] = C[i * matrixSize + (k+3)] + tA2 * B[matrixSize * (j+1) + (k+3)];
            C[i * matrixSize + (k+4)] = C[i * matrixSize + (k+4)] + tA2 * B[matrixSize * (j+1) + (k+4)];
            C[i * matrixSize + (k+5)] = C[i * matrixSize + (k+5)] + tA2 * B[matrixSize * (j+1) + (k+5)];
            C[i * matrixSize + (k+6)] = C[i * matrixSize + (k+6)] + tA2 * B[matrixSize * (j+1) + (k+6)];
            C[i * matrixSize + (k+7)] = C[i * matrixSize + (k+7)] + tA2 * B[matrixSize * (j+1) + (k+7)];
            C[i * matrixSize + (k+8)] = C[i * matrixSize + (k+8)] + tA2 * B[matrixSize * (j+1) + (k+8)];
            C[i * matrixSize + (k+9)] = C[i * matrixSize + (k+9)] + tA2 * B[matrixSize * (j+1) + (k+9)];
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
    A = new double[size * size];
    B = new double[size * size];
    C = new double[size * size];
    tmpC = new double[size * size];

    setMatrixData();

    calculationMatrixMultiplicationSimple();

    cout << "------ MatrixData Setting complete" << endl;
  }
  double measurement(){
    double startTime, endTime, resultTime = 0.0;
    int i;

    calculationMatrixMultiplication();
    
    for(i = 0; i < 5; i++){
      setMatrixDataC();
      startTime = omp_get_wtime();
      calculationMatrixMultiplication();
      endTime = omp_get_wtime();
      resultTime += (double)(endTime - startTime);
    }
    resultTime = resultTime / 5;

    //showMatrixData(20);

    if( checkMatrixData() == true ){
      cout << "MatrixData Error!!!" << endl;
      exit(1);
    }

    return resultTime;
  }
};