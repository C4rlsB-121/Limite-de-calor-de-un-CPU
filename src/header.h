#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Matrices
float** crearMatriz(int filasNum, int colNum);
void liberarMatriz(float** matriz, int filasNum);
void prodMatriz(int n, int m, int l, float** matrizA, float** matrizB, float** matrizR);
void copiarMatriz(int filasNum, int colNum, float** matrizNueva, float** matrizCopiada);

// Fuentes de calor
float calor_1Gauss(float x, float y, float t);
float calor_1Flat(float x, float y, float t);
float calor_3Gauss(float x, float y, float t);
float calor_4GaussSparced(float x, float y, float t);
