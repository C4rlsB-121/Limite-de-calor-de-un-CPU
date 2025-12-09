#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Matrices
float** crearMatriz(int filasNum, int colNum);
void liberarMatriz(float** matriz, int filasNum);
void prodMatriz(int n, int m, int l, float** matrizA, float** matrizB, float** matrizR);
void copiarMatriz(int filasNum, int colNum, float** matrizNueva, float** matrizCopiada);
void initZerosM(float** matriz, int filasNum, int colNum);
void initZerosV(float* v, int n);
void flat_Interior(float** matriz, int filasNum, int colNum, float* v_flat);

// Fuentes de calor
float calor_3Flat(float x, float y, float t);
float calor_3GaussUp(float x, float y, float t);
float calor_4GaussSparced(float x, float y, float t);
