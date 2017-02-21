# MatrixMultiplication
行列行列積を高速に解くプログラム  

# 実行結果

[環境1]  
CPU : Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz  
周波数 : 2.4GHz  
コア数 : 28  
NUMAノード : 2  
命令セット : AVX 2.0  
理論演算性能 : 2.4GHz * 28core * 2(ALU) * 4(SIMD) * 2(FMA) = 1075.2GFlops  
ループ数 : 1024 * 1024 * 1024ループ  
演算時間 : 0.014秒  
演算性能 : 145.3GFlops  

[環境2]  
CPU : Intel(R) Xeon(R) CPU E5-2695 v4 @ 2.10GHz  
周波数 : 2.1GHz  
コア数 : 36  
NUMAノード : 2  
命令セット : AVX 2.0  
理論演算性能 : 2.1GHz * 36core * 2(ALU) * 4(SIMD) * 2(FMA) = 1209.6GFlops  
ループ数 : 1024 * 1024 * 1024ループ  
演算時間 : 0.005秒  
演算性能 : 368.6GFlops  

# ファイル説明
・asm_g++.sh  
g++で作ったアセンブリデータをasmフォルダに保存するシェル  

・icpc.sh  
iccコンパイラを使った最適化オプションをつけたコンパイル  

・g++.sh  
g++コンパイラを使ったコンパイル  

・main.cpp  
エントリーポイント  

・matrix_multiplication.h  
行列行列積解くためのクラス  

# 実行方法
[1] sh g++.sh  
[2] ./main.out スレッド数  

# 使用した高速化手法
・メモリの連続アクセス  
・ループアンローリング  
・tmp変数によるパイプラインハザードの回避  
・openMPによる並列化  
・SIMD命令  
・ループ交換法  
・キャッシュブロック化  
・numactl --locallocコマンド  

