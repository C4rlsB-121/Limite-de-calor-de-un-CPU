#include<stdio.h> // Input Output
#include<stdlib.h> //standard library
#include<mpi.h> // Division del dominio
#include"header.h" // Enlazar funciones y objetos declarados en los demás archivos

#define Nx 11 // Numero de puntos interiores  en cada eje
#define Ny Nx
#define xmin 0  // Determinar espacio del dominio
#define xmax 1
#define ymin 0
#define ymax 1

#define Tf 1800
#define dt 30 // Tamaño de paos temporal

// ------------------------------- Main --------------------------------- //
int main(){

int rank, size;
float ancho = xmax - xmin;
float alto = ymax - ymin;
float dx = ancho / ((float)Nx - 1);
float dy = alto / ((float)Ny - 1);
int n_pts = Nx*Ny;
int my_coords[2];
MPI_Comm comm_estd, comm_cart;

// Declaraciones del modelo
float a = 9e-5;
float k = 152.24;
float eta = 180;
float T_out = 25;
float r = a * (float)dt / (dx * dy);

float nt = (float)Tf / (float)dt;


// Condiciones Iniciales
float** u = crearMatriz(Nx + 2, Ny + 2); // Considerando Nodos fantasma
for (int i = 0; i < Nx; i++){
	for (int j = 0; j < Ny; j++){
		u[i][j] = T_out;
	}
}


MPI_Init(NULL, NULL);
comm_estd = MPI_COMM_WORLD; // Comunicador estandar por defecto
MPI_Comm_rank(comm_estd, &rank);
MPI_Comm_size(comm_estd, &size);

// Creamos la topologia cartesiana
int dims[2] = {0,0};
int periods[2] = {0,0};
MPI_Dims_create(size, 2, dims); // Factoriza el numero de procesos y actualiza dims
int nprocs_x = dims[0], nprocs_y = dims[1];

MPI_Cart_create(comm_estd, 2, dims, periods, 1, &comm_cart);
MPI_Cart_coords(comm_cart, rank, 2, my_coords);
int coord_i = my_coords[0], coord_j = my_coords[1]; // Indexadas en 0

// Delimitamos propiedades de cada subdivision
int Nx_pN = Nx / nprocs_x;
int Ny_pN = Ny / nprocs_y;
int local_Nx = Nx_pN;
int local_Ny = Ny_pN;

int inicio_nx = Nx_pN * coord_i;
int inicio_ny = Ny_pN * coord_j;
float inicio_x = (float)inicio_nx * dx;
float inicio_y = (float)inicio_ny * dy;
//printf("Rank %d tiene inicios: (%f,%f)\n", rank, inicio_x, inicio_y);
//printf("Rank %d tiene coordenadas (%d, %d)\n", rank, local_Nx, local_Ny);

	/*
´	 Puede que no a todos los procesos se les asigne la misma cantidad de puntos en la malla
	 Hacer el ajusta para el último comunicador en cada dimension
	*/
if (coord_i == nprocs_x - 1) local_Nx = Nx - inicio_nx + 1;
if (coord_j == nprocs_y - 1) local_Ny = Ny - inicio_ny + 1;
int local_size = local_Nx * local_Ny;

// Vemos vecinos //
int up, down, right, left;
MPI_Cart_shift(comm_cart, 0, 1, &left, &right);
MPI_Cart_shift(comm_cart, 1, 1, &down, &up);

float** u_prev = crearMatriz(Nx+4, Ny+4);
copiarMatriz(Nx+4, Ny+4, u_prev, u);

/*
for (int n = 0; n < nt + 1 ; n++){ // Ciclo Temporal //
// Nodos fantasma
	if (down < 0){ // Frontera inferior
		for (int i = 0; i < local_Nx; i++){
			u[Ny + 3][i+1] = u[Ny+1][i+1] + (2*dy*eta) / k * (T_out - u[Ny + 2][i+1]);
		}
	}
	if (up < 0){ // Frontera superior
		for (int i = 0; i < local_Nx; i++){
			u[0][i+1] = u[2][i+1] + (2*dx*dy) / k * (T_out - u[1][i+1]);
		}
	}
	if (left < 0){ // Frontera Izquierda
		for (int j = 0; j < Ny; j++){
			u[j+1][0] = u[j+1][2] + (2*dx*dy) / k * (T_out - u[j+1][1]);
		}
	}
	if (right < 0){ // Frontera Derecha
		for (int j = 0; j < Ny; j++){
		u[j+1][Nx +3] = u[j+1][Nx + 1] + (2*dx*dy) / k * (T_out - u[j+1][Nx+2]);
		}
	}
	// Los demas nodos usando método de lineas //
}

*/



// Vemos dominio
for (int k =0; k<size;k++){
if (rank == k){
printf("Rank: %d: Coord (%d,%d)\n\n", rank, coord_i,coord_j);
for (int i = 0; i < local_Nx+2; i++){
	for (int j = 0; j < local_Ny+2; j++) printf("(x,y)=(%f,%f)\n num nodos=%d \n", inicio_x + i*dx, inicio_y + j*dy, local_size);
	}
}
MPI_Barrier(MPI_COMM_WORLD);
}


liberarMatriz(u, Nx);
liberarMatriz(u_prev, Nx);
MPI_Finalize();
}




void saveImg(int rank){
	if (rank != 0) return;
}
