# MatrixMultiplication
行列行列積を高速に解くプログラム  

# 実行結果

[環境1]  
CPU : Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz  
周波数 : 2.4GHz  
コア数 : 28  
NUMAノード : 2  
理論演算性能 : 2.4 * 28 * 2 = 134.4GFlops  
演算時間 : 0.0秒  
演算性能 : 00.0GFlops  

[環境2]  
CPU : Intel(R) Xeon(R) CPU E5-2695 v4 @ 2.10GHz  
周波数 : 2.4GHz  
コア数 : 36  
NUMAノード : 2  
理論演算性能 : 2.1 * 36 * 2 = 151.2GFlops  
演算時間 : 0.0秒  
演算性能 : 00.0GFlops  

# ファイル説明
・asm_g++.sh  
g++で作ったアセンブリデータをasmフォルダに保存するシェル  

・icpc.sh  
iccコンパイラを使ったコンパイル  

・g++.sh  
g++コンパイラを使ったコンパイル  

・main.cpp  
エントリーポイント  

・matrix_multiplication.h  
行列行列積解くためのクラス  


# 使用した高速化手法
・メモリの連続アクセス  
・ループアンローリング  
・tmp変数によるパイプラインハザードの回避  
・openMPによる並列化  
・SIMD命令  
・ループ交換法  
・キャッシュブロック化  
・numactl --locallocコマンド  

