#include<stdio.h> // Input Output
#include<stdlib.h> //standard library
#include<mpi.h> // Division del dominio
#include"header.h" // Enlazar funciones y objetos declarados en los demás archivos

#define Nx 13 // Numero de puntos en el cerrado por eje
#define Ny Nx
#define xmin 0  // Determinar espacio del dominio
#define xmax 1
#define ymin 0
#define ymax 1

#define Tf 1800
#define dt 30 // Tamaño de paos temporal
#define tipo_calor 1

#define num_frames 100
#define animate 3

// -------------- Funciones ------------------- //
void show_domain(int size, int rank, int local_Nx, int local_Ny, int inicio_nx, int inicio_ny,
                 float dx, float dy);
void set_emision(float (*emision)(float, float, float), int tipo_emision);
void preparar_envio(float** matrizOrigen, int local_Ny, int local_Nx, int left, int right, int up, int down,
		 float* send_left, float* send_right, float* send_up, float* send_down);
void halo_Update(float** subMatriz, int local_Ny, int local_Nx, int left, int right, int up, int down,
                 float* recv_left, float* recv_right, float* recv_up, float* recv_down);
void save_ppm_frame(const char* filename, float** data, int n_filas, int n_columnas);
// ------------------------------- Main --------------------------------- //
int main(){

int rank, size;
float ancho = xmax - xmin;
float alto = ymax - ymin;
float dx = ancho / ((float)Nx - 1);
float dy = alto / ((float)Ny - 1);
int my_coords[2];
MPI_Comm comm_estd, comm_cart;

// Declaraciones del modelo
float a = 9e-5;
float k = 152.24;
float eta = 180;
float T_out = 25;
float r = a / (dx * dy);
int tipo_emision = 4;
float T_limite = 100;
float (*emision)(float, float, float);
set_emision(emision, tipo_calor);


float nt = (float)Tf / (float)dt;
int frame_time = (nt + 1) / (int)num_frames;

MPI_Init(NULL, NULL);
float t1 = MPI_Wtime();
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

// Vemos vecinos //
int up, down, right, left;
MPI_Cart_shift(comm_cart, 0, 1, &left, &right);
MPI_Cart_shift(comm_cart, 1, 1, &down, &up);

// Definiciones necesarias para el metodo de lineas //


//// Condiciones Iniciales en subcuadriculas con nodos fantasma y halo
float** u = crearMatriz(local_Ny + 2, local_Nx + 2);
/*
for (int i = 0; i <local_Nx + 2; i++){
        for (int j = 0; j < local_Ny + 2; j++){
                u[0][0] = T_out;
        }
}
*/
float** u_prev = crearMatriz(local_Ny+2, local_Nx+2);
//copiarMatriz(local_Ny+2, local_Nx+2, u_prev, u);
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

//initZerosV(send_left, local_Ny);
//initZerosV(send_right, local_Ny);
//initZerosV(send_up, local_Nx);
//initZerosV(send_down, local_Nx);
//initZerosV(recv_left, local_Ny);
//initZerosV(recv_right, local_Ny);
//initZerosV(recv_up, local_Nx);
//initZerosV(recv_down, local_Nx);

int num_requests;

//	Ciclo Temporal	//
for (int n = 0; n <  nt + 1 ; n++){
    /*
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
	*/
	/*
// Los demas nodos usando método de lineas particularmente RK 4 //
	//  ------------------------ K1 ---------------------//
	for (int j = 1; j < local_Ny; j++){
		for (int i = 1; i < local_Nx+1; i++){
			k1[local_Ny + 1 - j][i] = r*u_prev[local_Ny + 1 - j][i - 1] - 4*r*u_prev[local_Ny + 1 - j][i] + r*u_prev[local_Ny + 1 - j][i+1]
							+ r*u_prev[local_Ny + 1 - j - 1][i] + r*u_prev[local_Ny + 1 - j + 1][i] + emision( inicio_x + (i-1)*dx, inicio_y + (j-1)*dy, n*dt);
		}
	}
	*/

	// Armar array a enviar
	//preparar_envio(k1, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);


	// Comunicación
	/*
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
        */
	// Actualizamos el halo
	//halo_Update(k1, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);


	// --------------------------------------- Paso temporal ----------------------------------------------

        // Esperamos a los procesos
        //MPI_Waitall(num_requests, requests_k4, MPI_STATUSES_IGNORE);

        /*
        // Actualizamos el halo
        halo_Update(k4, local_Ny, local_Nx, left, right, up, down, recv_left, recv_right, recv_up, recv_down);

        // Damos el paso
        for (int j = 0; j < local_Ny; j++){
                for (int i = 0; i < local_Ny; i++){
                        u[local_Ny - j][i] = u_prev[local_Ny - j][i] + dt/6 * (k1[local_Ny - j][i] + 2*k2[local_Ny - j][i] + 2*k3[local_Ny - j][i] + k4[local_Ny - j][i]);
                        }
                }

        copiarMatriz(local_Ny + 2, local_Nx + 2, u_prev, u);

        // Vemos que no se haya pasado el limite
        for (int j = 2; j < local_Ny; j++){
                for (int i = 2; i < local_Nx; i++){
                        if (u[j][i] > T_limite) printf("%f\n", u[j][i]);// printf("Peligro: se ha alcanzado la temperatura máxima permitida en %d\n", dt*n);
                }
        }
        */
        /*
        // Animación
        if (size && animate && rank == 0){
                if (n % frame_time == 0){
                        char filename[100];
                        sprintf(filename, "frames/frame_%06d", n/frame_time);
                        save_ppm_frame(filename, u, local_Ny + 2, local_Nx + 2);
                }
        }
        */
}

//show_domain(size, rank, local_Nx, local_Ny, inicio_nx, inicio_ny, dx, dy);

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

MPI_Barrier(MPI_COMM_WORLD);
float t2 = MPI_Wtime();
float time = t2 - t1;
if (rank == 0) printf("----------------------------------------------------\nTiempo de ejecución con %d nucleos: %f\n-------------------------------------------------------\n", size, time);

MPI_Finalize();
} // End main


void show_domain(int size, int rank,  int local_Nx, int local_Ny, int inicio_nx, int inicio_ny,
                 float dx, float dy){
        // Vemos dominio
        for (int k = 0; k<size; k++){
                if (rank == k){
                        for (int i = 0; i < local_Nx; i++){
                                for (int j = 0; j < local_Ny; j++) {
                                        printf("(x,y)=(%f,%f)\n num nodos=%d \n", (inicio_nx+i) * dx, (inicio_ny+j) * dy, local_Nx*local_Ny);
                                }
                        }
                printf("\n");
                }
                MPI_Barrier(MPI_COMM_WORLD);
        }
}

void set_emision(float (*emision)(float, float, float), int tipo_emision){
        switch (tipo_emision){
        case 1: emision = calor_3Flat;
        case 2: emision = calor_3GaussUp;
        case 3: emision = calor_4GaussSparced;
        }
}

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


void save_ppm_frame(const char* filename, float** data, int n_filas, int n_columnas){
        FILE *fp = fopen(filename, "wb");
        fprintf(fp, "P6\n");
        fprintf(fp, "%d %d\n", n_columnas, n_filas);  // Ancho, Alto
        fprintf(fp, "255\n");

        float min_val = data[0][0];
        float max_val = data[0][0];

        for (int j = 0; j < n_filas; j++) {
                for (int i = 0; i < n_columnas; i++) {
                        if (data[j][i] < min_val) min_val = data[j][i];
                        if (data[j][i] > max_val) max_val = data[j][i];
                }
        }

    float range = max_val - min_val;
        if (range == 0) range = 1.0f; // Evitar divisionn por 0

        // Escribir datos de píxeles
        for (int j = 0; j < n_filas; j++) {
                for (int i = 0; i < n_columnas; i++) {
                        float valor = data[j][i];

                        // Normalizar
                        valor = (valor - min_val) / range;
                        valor = fmaxf(0.0f, fminf(1.0f, valor));  // Usar fmaxf/fminf para float

                        unsigned char r, g, b;

                        // Mapa de colores azul -> verde -> rojo
                        if (valor < 0.25f) {
                                r = 0;
                                g = (unsigned char)(valor / 0.25f * 255);
                                b = 255;
                        }
                        else if (valor < 0.5f) {
                                r = 0;
                                g = 255;
                                b = (unsigned char)((0.5f - valor) / 0.25f * 255);
                        }
                        else if (valor < 0.75f) {
                                r = (unsigned char)((valor - 0.5f) / 0.25f * 255);
                                g = 255;
                                b = 0;
                        }
                        else {
                                r = 255;
                                g = (unsigned char)((1.0f - valor) / 0.25f * 255);
                                b = 0;
                        }

                        fwrite(&r, 1, 1, fp);
                        fwrite(&g, 1, 1, fp);
                        fwrite(&b, 1, 1, fp);
                }
        }
        fclose(fp);
}
