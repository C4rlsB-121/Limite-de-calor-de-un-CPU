#include<stdio.h> // Input Output
#include<stdlib.h> //standard library
#include<mpi.h> // Division del dominio
#include"header.h" // Enlazar funciones y objetos declarados en los demás archivos

#define Nx 22 // Numero de puntos en el cerrado por eje
#define Ny Nx
#define xmin 0  // Determinar espacio del dominio
#define xmax 1
#define ymin 0
#define ymax 1

#define Tf 1800
#define dt 30 // Tamaño de paos temporal

// -------------- Funciones ------------------- //
void preparar_envio(float** matrizOrigen, int local_Ny, int local_Nx, int left, int right, int up, int down,
		 float* send_left, float* send_right, float* send_up, float* send_down);
void halo_Update(float** subMatriz, int local_Ny, int local_Nx, int left, int right, int up, int down,
                 float* recv_left, float* recv_right, float* recv_up, float* recv_down);

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
float r = a / (dx * dy);
int tipo_emision = 4;

float nt = (float)Tf / (float)dt;

float (*emision)(float, float, float);
switch (tipo_emision){
	case 1: emision = calor_1Gauss;
	case 2: emision = calor_1Flat;
	case 3: emision = calor_3Gauss;
	case 4: emision = calor_4GaussSparced;
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
if (coord_i == nprocs_x - 1) local_Nx = Nx - inicio_nx;
if (coord_j == nprocs_y - 1) local_Ny = Ny - inicio_ny;
int local_size = local_Nx * local_Ny;

// Vemos vecinos //
int up, down, right, left;
MPI_Cart_shift(comm_cart, 0, 1, &left, &right);
MPI_Cart_shift(comm_cart, 1, 1, &down, &up);

// Definiciones necesarias para el metodo de lineas //


//// Condiciones Iniciales en subcuadriculas con nodos fantasma y halo
float** u = crearMatriz(local_Nx + 2, local_Ny + 2);
for (int i = 0; i <local_Nx + 2; i++){
        for (int j = 0; j < local_Ny + 2; j++){
                u[0][0] = T_out;
        }
}
float** u_prev = crearMatriz(local_Ny+2, local_Nx+2);
copiarMatriz(local_Ny+2, local_Nx+2, u_prev, u);

//// pasos de RK
float** k1 = crearMatriz(local_Ny + 2, local_Nx + 2);
float** k2 = crearMatriz(local_Ny + 2, local_Nx + 2);
float** k3 = crearMatriz(local_Ny + 2, local_Nx + 2);
float** k4 = crearMatriz(local_Ny + 2, local_Nx + 2);


//// vectores a enviar y recibir
float* send_left = malloc(local_Ny * sizeof(float));
float* send_right = malloc(local_Ny * sizeof(float));
float* send_up = malloc(local_Nx * sizeof(float));
float* send_down = malloc(local_Nx * sizeof(float));

float* recv_left = malloc(local_Ny * sizeof(float));
float* recv_right = malloc(local_Ny * sizeof(float));
float* recv_up = malloc(local_Nx * sizeof(float));
float* recv_down = malloc(local_Nx * sizeof(float));

initZerosV(send_left, local_Ny);
initZerosV(send_right, local_Ny);
initZerosV(send_up, local_Nx);
initZerosV(send_down, local_Nx);
initZerosV(recv_left, local_Ny);
initZerosV(recv_right, local_Ny);
initZerosV(recv_up, local_Nx);
initZerosV(recv_down, local_Nx);

int num_requests;

//	Ciclo Temporal	//
for (int n = 0; n <  nt + 1 ; n++){
// Nodos fantasma
	if (down < 0){ // Frontera inferior
		for (int i = 0; i < local_Nx; i++){
			u[local_Ny + 1][i + 1] = u[local_Ny-1][i + 1] + (2*dy*eta) / k * (T_out - u[local_Ny][i +1]);
		}
	}
	if (up < 0){ // Frontera superior
		for (int i = 0; i < local_Nx; i++){
			u[0][i + 1] = u[2][i + 1] + (2*dy*eta) / k * (T_out - u[1][i + 1]);
		}
	}
	if (left < 0){ // Frontera Izquierda
		for (int j = 0; j < local_Ny; j++){
			u[j + 1][0] = u[j + 1][2] + (2*dx*eta) / k * (T_out - u[j + 1][1]);
		}
	}
	if (right < 0){ // Frontera Derecha
		for (int j = 0; j < local_Ny; j++){
			u[j + 1][local_Nx + 1] = u[j + 1][Nx - 1] + (2*dx*eta) / k * (T_out - u[j + 1][local_Nx]);
		}
	}
// Los demas nodos usando método de lineas particularmente RK 4 //
	//  ------------------------ K1 ---------------------//
	for (int j = 1; j < local_Ny; j++){
		for (int i = 1; i < local_Nx+1; i++){
			k1[local_Ny + 1 - j][i] = r*u_prev[local_Ny + 1 - j][i - 1] - 4*r*u_prev[local_Ny + 1 - j][i] + r*u_prev[local_Ny + 1 - j][i+1]
							+ r*u_prev[local_Ny + 1 - j - 1][i] + r*u_prev[local_Ny + 1 - j + 1][i] + emision( inicio_x + (i-1)*dx, inicio_y + (j-1)*dy, n*dt);
		}
	}


	// Armar array a enviar
	preparar_envio(k1, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);


	// Comunicación
        MPI_Request requests_k1[8];
        num_requests = 0;
	if (left > -1){ // Tiene vecino izquierdo
		MPI_Isend(send_left, local_Ny, MPI_FLOAT, left, 0, MPI_COMM_WORLD, &requests_k1[num_requests++] );
		MPI_Irecv(recv_left, local_Ny, MPI_FLOAT, left, 1, MPI_COMM_WORLD, &requests_k1[num_requests++] );
	}
	if (right > -1){ // Tiene vecino derecho
		MPI_Isend(send_right, local_Ny, MPI_FLOAT, right, 1, MPI_COMM_WORLD, &requests_k1[num_requests++] );
                MPI_Irecv(recv_right, local_Ny, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &requests_k1[num_requests++] );
        }
	if (up > -1){ // Tiene vecino superior
		MPI_Isend(send_up, local_Nx, MPI_FLOAT, up, 2, MPI_COMM_WORLD, &requests_k1[num_requests++] );
                MPI_Irecv(recv_up, local_Nx, MPI_FLOAT, up, 3, MPI_COMM_WORLD, &requests_k1[num_requests++] );
        }
	if (down > -1){ // Tiene vecino inferior
		MPI_Isend(send_down, local_Nx, MPI_FLOAT, down, 3, MPI_COMM_WORLD, &requests_k1[num_requests++] );
                MPI_Irecv(recv_down, local_Nx, MPI_FLOAT, down, 2, MPI_COMM_WORLD, &requests_k1[num_requests++] );
        }

	// Esperamos a los procesos
	MPI_Waitall(num_requests, requests_k1, MPI_STATUSES_IGNORE);

	// Actualizamos el halo
	halo_Update(k1, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);

	// -----------------------  K2 -------------------------------------- //
	for (int j = 1; j < local_Ny; j++){
                for (int i = 1; i < local_Nx+1; i++){
                        k2[local_Ny + 1 - j][i] = r*(u_prev[local_Ny + 1 - j][i - 1] + 0.5*dt*k1[local_Ny + 1 - j][i -1]) - 4*r*(u_prev[local_Ny + 1 - j][i] + 0.5*dt*k1[local_Ny + 1 -j][i]) +
			r*(u_prev[local_Ny + 1 - j][i+1] + 0.5*dt*k1[local_Ny + 1 - j][i + 1]) + r*(u_prev[local_Ny + 1 - j - 1][i] + 0.5*dt*k1[local_Ny + 1 - j - 1][i]) +
			r*(u_prev[local_Ny + 1 - j + 1][i] + 0.5*dt*k1[local_Ny + 1 - j +1][i]) + emision( inicio_x + (i-1)*dx, inicio_y + (j-1)*dy, (n+0.5)*dt);
                }
        }

	// Armar array a enviar
        preparar_envio(k2, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);


        // Comunicación
        MPI_Request requests_k2[8];
        num_requests = 0;
        if (left > -1){ // Tiene vecino izquierdo
                MPI_Isend(send_left, local_Ny, MPI_FLOAT, left, 0, MPI_COMM_WORLD, &requests_k2[num_requests++] );
                MPI_Irecv(recv_left, local_Ny, MPI_FLOAT, left, 1, MPI_COMM_WORLD, &requests_k2[num_requests++] );
        }
        if (right > -1){ // Tiene vecino derecho
                MPI_Isend(send_right, local_Ny, MPI_FLOAT, right, 1, MPI_COMM_WORLD, &requests_k2[num_requests++] );
                MPI_Irecv(recv_right, local_Ny, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &requests_k2[num_requests++] );
        }
        if (up > -1){ // Tiene vecino superior
                MPI_Isend(send_up, local_Nx, MPI_FLOAT, up, 2, MPI_COMM_WORLD, &requests_k2[num_requests++] );
                MPI_Irecv(recv_up, local_Nx, MPI_FLOAT, up, 3, MPI_COMM_WORLD, &requests_k2[num_requests++] );
        }
        if (down > -1){ // Tiene vecino inferior
                MPI_Isend(send_down, local_Nx, MPI_FLOAT, down, 3, MPI_COMM_WORLD, &requests_k2[num_requests++] );
                MPI_Irecv(recv_down, local_Nx, MPI_FLOAT, down, 2, MPI_COMM_WORLD, &requests_k2[num_requests++] );
        }

        // Esperamos a los procesos
        MPI_Waitall(num_requests, requests_k2, MPI_STATUSES_IGNORE);

        // Actualizamos el halo
        halo_Update(k2, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);

	// ------------------------------------- K3 ----------------------------- //
        for (int j = 1; j < local_Ny; j++){
                for (int i = 1; i < local_Nx+1; i++){
                        k3[local_Ny + 1 - j][i] = r*(u_prev[local_Ny + 1 - j][i - 1] + 0.5*dt*k2[local_Ny + 1 - j][i -1]) - 4*r*(u_prev[local_Ny + 1 - j][i] + 0.5*dt*k2[local_Ny + 1 -j][i]) +
                        r*(u_prev[local_Ny + 1 - j][i+1] + 0.5*dt*k2[local_Ny + 1 - j][i + 1]) + r*(u_prev[local_Ny + 1 - j - 1][i] +0.5*dt*k2[local_Ny + 1 - j - 1][i]) +
                        r*(u_prev[local_Ny + 1 - j + 1][i] + 0.5*dt*k2[local_Ny + 1 - j +1][i]) + emision( inicio_x + (i-1)*dx, inicio_y + (j-1)*dy, (n+0.5)*dt);
                }
        }

        // Armar array a enviar
        preparar_envio(k3, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);


        // Comunicación
        MPI_Request requests_k3[8];
        num_requests = 0;
        if (left > -1){ // Tiene vecino izquierdo
                MPI_Isend(send_left, local_Ny, MPI_FLOAT, left, 0, MPI_COMM_WORLD, &requests_k3[num_requests++] );
                MPI_Irecv(recv_left, local_Ny, MPI_FLOAT, left, 1, MPI_COMM_WORLD, &requests_k3[num_requests++] );
        }
        if (right > -1){ // Tiene vecino derecho
                MPI_Isend(send_right, local_Ny, MPI_FLOAT, right, 1, MPI_COMM_WORLD, &requests_k3[num_requests++] );
                MPI_Irecv(recv_right, local_Ny, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &requests_k3[num_requests++] );
        }
        if (up > -1){ // Tiene vecino superior
                MPI_Isend(send_up, local_Nx, MPI_FLOAT, up, 2, MPI_COMM_WORLD, &requests_k3[num_requests++] );
                MPI_Irecv(recv_up, local_Nx, MPI_FLOAT, up, 3, MPI_COMM_WORLD, &requests_k3[num_requests++] );
        }
        if (down > -1){ // Tiene vecino inferior
                MPI_Isend(send_down, local_Nx, MPI_FLOAT, down, 3, MPI_COMM_WORLD, &requests_k3[num_requests++] );
                MPI_Irecv(recv_down, local_Nx, MPI_FLOAT, down, 2, MPI_COMM_WORLD, &requests_k3[num_requests++] );
        }

        // Esperamos a los procesos
        MPI_Waitall(num_requests, requests_k3, MPI_STATUSES_IGNORE);

        // Actualizamos el halo
        halo_Update(k3, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);


	// ------------------------------------- K4 ----------------------------- //
        for (int j = 1; j < local_Ny; j++){
                for (int i = 1; i < local_Nx+1; i++){
                        k4[local_Ny + 1 - j][i] = r*(u_prev[local_Ny + 1 - j][i - 1] + dt*k3[local_Ny + 1 - j][i -1]) - 4*r*(u_prev[local_Ny + 1 - j][i] + dt*k3[local_Ny + 1 -j][i]) +
                        r*(u_prev[local_Ny + 1 - j][i+1] + dt*k1[local_Ny + 1 - j][i + 1]) + r*(u_prev[local_Ny + 1 - j - 1][i] + dt*k3[local_Ny + 1 - j - 1][i]) +
                        r*(u_prev[local_Ny + 1 - j + 1][i] + dt*k3[local_Ny + 1 - j +1][i]) + emision( inicio_x + (i-1)*dx, inicio_y + (j-1)*dy, (n+1)*dt);
                }
        }

        // Armar array a enviar
        preparar_envio(k4, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);


        // Comunicación
        MPI_Request requests_k4[8];
        num_requests = 0;

        if (left > -1){ // Tiene vecino izquierdo
                MPI_Isend(send_left, local_Ny, MPI_FLOAT, left, 0, MPI_COMM_WORLD, &requests_k4[num_requests++] );
                MPI_Irecv(recv_left, local_Ny, MPI_FLOAT, left, 1, MPI_COMM_WORLD, &requests_k4[num_requests++] );
        }
        if (right > -1){ // Tiene vecino derecho
                MPI_Isend(send_right, local_Ny, MPI_FLOAT, right, 1, MPI_COMM_WORLD, &requests_k4[num_requests++] );
                MPI_Irecv(recv_right, local_Ny, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &requests_k4[num_requests++] );
        }
        if (up > -1){ // Tiene vecino superior
                MPI_Isend(send_up, local_Nx, MPI_FLOAT, up, 2, MPI_COMM_WORLD, &requests_k4[num_requests++] );
                MPI_Irecv(recv_up, local_Nx, MPI_FLOAT, up, 3, MPI_COMM_WORLD, &requests_k4[num_requests++] );
        }
        if (down > -1){ // Tiene vecino inferior
                MPI_Isend(send_down, local_Nx, MPI_FLOAT, down, 3, MPI_COMM_WORLD, &requests_k4[num_requests++] );
                MPI_Irecv(recv_down, local_Nx, MPI_FLOAT, down, 2, MPI_COMM_WORLD, &requests_k4[num_requests++] );
        }

        // Esperamos a los procesos
        MPI_Waitall(num_requests, requests_k4, MPI_STATUSES_IGNORE);

        // Actualizamos el halo
        halo_Update(k4, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);





}

//printf("%f\n", u[0][1]);


/*
// Vemos dominio
for (int k =0; k<size;k++){
if (rank == k){
printf("Rank: %d: Coord (%d,%d)\n\n", rank, coord_i,coord_j);
for (int i = 0; i < local_Nx; i++){
	for (int j = 0; j < local_Ny; j++) printf("(x,y)=(%f,%f)\n num nodos=%d \n", inicio_x + i*dx, inicio_y + j*dy, local_size);
	}
}
MPI_Barrier(MPI_COMM_WORLD);
}
*/



liberarMatriz(u, local_Ny + 2);
liberarMatriz(u_prev, local_Ny+2);
liberarMatriz(k1, local_Ny+2);
liberarMatriz(k2, local_Ny+2);
liberarMatriz(k3, local_Ny+2);
liberarMatriz(k4, local_Ny+2);
free(send_left);
free(send_right);
free(send_up);
free(send_down);
free(recv_left);
free(recv_right);
free(recv_up);
free(recv_down);

MPI_Finalize();
} // End main





void preparar_envio(float** matrizOrigen, int local_Ny, int local_Nx, int left, int right, int up, int down,
		    float* send_left, float* send_right, float* send_up, float* send_down){
	if (left > -1) for (int j = 0; j < local_Ny; j++) send_left[j] = matrizOrigen[j+1][0];
	if (right > -1) for (int j = 0; j < local_Ny; j++) send_right[j] = matrizOrigen[j+1][local_Nx + 1];
	if (up > -1) for (int i = 0; i < local_Ny; i++) send_up[i] = matrizOrigen[0][i +1];
	if (down > -1) for (int i = 0; i < local_Ny; i++) send_down[i] = matrizOrigen[local_Ny + 1][i + 1];
}

void halo_Update(float** subMatriz, int local_Ny, int local_Nx, int left, int right, int up, int down,
		 float* recv_left, float* recv_right, float* recv_up, float* recv_down){
		if (left > -1){
			for (int j = 0; j < local_Ny; j++) subMatriz[j + 1][0] = recv_left[j];
		}
		if (right > -1){
			for (int j = 0; j < local_Ny; j++) subMatriz[j + 1][local_Nx + 1] = recv_right[j];
		}
		if (up > -1){
			for (int i = 0; i < local_Nx; i++) subMatriz[0][i + 1] = recv_up[i];
		}
		if (down > -1){
			for (int i = 0; i < local_Nx; i++) subMatriz[local_Ny + 1][i + 1] = recv_down[i];
		}
}



void saveImg(int rank){
	if (rank != 0) return;
}
