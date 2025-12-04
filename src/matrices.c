#include<stdlib.h>
#include"header.h"

float** crearMatriz(int filasNum, int colNum){
        float**  matriz = malloc(filasNum * sizeof(float*));
        for (int i = 0; i < filasNum; i++){
                matriz[i]  = malloc(colNum * sizeof(float));
        }
        return matriz;
}

void liberarMatriz(float** matriz, int filasNum){
        for (int i = 0; i < filasNum; i++){
                free(matriz[i]);
        }
        free(matriz);
}



void prodMatriz(int n, int m, int l, float** matrizA, float** matrizB, float** matrizR){
        for (int i = 0; i < n; i++){
                for (int j = 0; j < l; j++){
                        matrizR[i][j] = 0;
                }
        }

        for (int i = 0; i < n; i++){
                for (int j = 0; j < l; j++){
                        for (int k = 0; k < m; k++){
                                matrizR[i][j] += matrizA[i][k] * matrizB[k][j];
                        }
                }
        }

}


void copiarMatriz(int filasNum, int colNum, float** matrizNueva, float** matrizCopiada){
	for (int i = 0; i < filasNum; i++){
		for (int j = 0; j < colNum; j++){
			matrizNueva[i][j] = matrizCopiada[i][j];
		}
	}
}


void initZerosM(float** matriz, int filNum, int colNum){
        for (int i = 0; i < filNum; i++){
                for (int j = 0; j < colNum; j++){
                        matriz[i][j] = 0.0;
                }
        }

}


void initZerosV(float* v, int n){
	for (int i = 0; i < n; i++) v[i] = 0;
}
