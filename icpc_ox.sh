NAME=$1
echo $NAME
echo "icpc $NAME -o main.out -fopenmp -O2 -march=core-avx2"
icpc $NAME -o main.out -fopenmp -O2 -march=core-avx2