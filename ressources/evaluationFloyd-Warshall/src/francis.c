//
// Created by Anas Francis on 04/05/2021.
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <limits.h>

int rank; int numproc;
long ** result_scatter;
int r;
int int_size_of_mat;
int N;


#define INFINI LONG_MAX
#define ALLOC_SIZE 2
struct t_table {
    long **tab;
    int nb_line;
    int nb_colon;
    int nb_elements;
};

typedef struct t_table *Matrix;

 /**
  * prototypes
  * */


Matrix create_matrix(int alloc_size);// O(1) in parallel
void destroy_matrix(Matrix *table);// O(1) in parallel
void print(Matrix mat);// O(N²)
Matrix transpose(Matrix mat);// O(1) in parallel
long* lineariser (Matrix A);// O(1) in parallel
Matrix read_file(char** argv);// O(N) size of the file (one iteration)
long* lineariser_column(long** matrix, int r_x_nb_lines);// O(1) in parallel
void scatter_line(Matrix A,long* save_data,int recv_count);// O(p-1)(l+Bt)
void scatter_colum(Matrix A,long* save_data,int recv_count); // O(p-1)(l+Bt)
void broadcast(int* size_of_mat);// O(p-1)(l+Bt) broadcast the size of matrix to all process
Matrix create_matrix_from_table(long *tab);// O(1) in parallel
void circuler(long** r_line_tmp, int *size); // O(1) each process is going to send/receive
void gather(long* data,int nb_to_send,long* save_into,int nb_to_recv);// O(p-1)(l+Bt)
void calcule_r_line_r_colon(long *save_line,long * save_column,
                            long ** N_r_matrix,int r_x_nb_lines, int r_colone, int start);
long minimum(long a, long b); // O(1)
Matrix compute_and_get_w_from(Matrix mat);// O(1) see specification on https://lms.univ-cotedazur.fr/course/view.php?id=1302
int get_start_line(int r_lines);// O(1) knowing the rank of process and the value of r we can guess from where to start
void initialization(char** argv,Matrix* address_A, Matrix* address_w, int* add_r_x_nb_lines); // read files and
                                                                                            // initilize the matrix w
                                                                                            // O(N)
int mod(int,int);// O(1) gives a positive number when making a%b (even if a < 0)


Matrix create_matrix(int alloc_size) {
    Matrix t         = malloc(sizeof(struct t_table));
    long **tab = (long**)malloc(alloc_size * sizeof(long*));
    #pragma omp parallel for
    for(int i=0; i<alloc_size;i++)
        tab[i] = (long*)malloc(alloc_size * sizeof (long ));
    if (!t || !tab) {
        fprintf(stderr, "cannot allocate memory\n");
        return 0;
    }
    t->tab = tab;
    t->nb_line = t->nb_colon = alloc_size;
    t->nb_elements = 0;
    return t;
}


Matrix transpose(Matrix mat){
    Matrix to_ret = create_matrix(mat->nb_line);
    #pragma omp parallel for
    for(int i = 0; i < mat->nb_line;i++){
        #pragma omp parallel for
        for(int j = 0; j < mat->nb_line;j++){
            to_ret->tab[i][j] = mat->tab[j][i];
        }
    }
    return to_ret;
}


void destroy_matrix(Matrix *table){
    Matrix T = *table;
    #pragma omp parallel for
    for(int i = 0;i<T->nb_line;i++)
        free(T->tab[i]);
    free(T);
    *table = NULL; // good programmer habit
}


long* lineariser (Matrix A){

    long* table = malloc(sizeof(long)*A->nb_line*A->nb_colon);
    #pragma omp parallel for
    for(int i = 0;i<A->nb_line;i++){
        #pragma omp parallel for
        for(int j = 0; j<A->nb_colon; j++){
            table[i*A->nb_colon+j] = A->tab[i][j];
        }
    }
    return table;
}


static int add_word(Matrix *matrix, int index_line,int index_colon, int my_nb) {
    Matrix T           = *matrix;
    long **tab = T->tab;
    int nb_lign        = T->nb_line;
    int nb_colon         = T->nb_colon;
    if (index_colon>=nb_colon) {
        // La table est pleine, la "rallonger" avant d'essayer d'insérer my_nb
        nb_lign += 1;
        nb_colon += 1;
        tab = realloc(tab, nb_lign*sizeof(long*));
        if (!tab) {
            fprintf(stderr, "cannot reallocate memory\n");
            return 0;
        }
        int i=0;
        for(i=0; i<nb_lign-1;i++){
            tab[i] = realloc(tab[i] ,nb_colon * sizeof (long ));
        }
        tab[i] = (long*)malloc(nb_colon * sizeof (long ));
        // conserver les nouvelles valeurs dans la table
        T->tab = tab;
        T->nb_line = nb_lign;
        T->nb_colon = nb_colon;
    }
    // Insérer le nouveau number à  la position index
    tab[index_line][index_colon]= my_nb;
    // On a un number de plus dans la table
    T->nb_elements += 1;
    return 1;                     // car ce number apparaît une fois
}


Matrix read_file(char** argv){
    FILE *fp = fopen(argv[1], "r");
    if (!fp)
    {
        fprintf(stderr, "couldn't open file name %s...\n\nProgramme exit with code error -1.\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    else
    {
        char ch, buff[50];
        int j = 0, index_of_colon = 0,index_lign=0;
        Matrix to_ret = create_matrix(ALLOC_SIZE);
        do
        {
            ch = fgetc(fp);
            buff[j++]=ch;
            if(ch==' ' || ch=='\t' || ch == '\n'){
                buff[j] = '\0';
                int value = atoi(buff);
                add_word(&to_ret,index_lign,index_of_colon++, value);
                if(ch=='\n'){
                    index_lign ++;
                    index_of_colon=0;
                }
                j = 0;
            }
        }while(ch!=EOF);
        fclose(fp);
        return to_ret;
    }
}


long* lineariser_column(long** matrix, int r_x_nb_lines){
    long* to_ret = malloc(sizeof(long)* r_x_nb_lines);
    int r_line = r_x_nb_lines/N;
    #pragma omp parallel for
    for(int i = 0;i<r_line;i++){
        #pragma omp parallel for
        for(int j = 0;j<N;j++){
            to_ret[j+i*N] = matrix[j][i];
        }
    }
    return to_ret;
}


void scatter_line(Matrix A,long* save_data,int recv_count){
    MPI_Status status;
    if(rank == 0){ // if we're the sender

        long* linearize_lines = lineariser(A);
        #pragma omp parallel for
        for(int i_lines = 0 ; i_lines<recv_count;i_lines++){ //save the 1st r lines inside the buffer
            save_data[i_lines] = linearize_lines[i_lines];
        }

        // send the rest of elements to the next process. Starting from the last r x N elements.
        for(int i =1; i<numproc;i++){
            MPI_Send(linearize_lines+(int_size_of_mat-(r*A->nb_colon)*(i)),r*A->nb_colon,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
        }
    }

    else{
        long* recieved;
        for(int i = 0; i<numproc-rank-1;i++){
            recieved = malloc(sizeof(long) * recv_count);
            MPI_Recv(recieved,recv_count,MPI_LONG,((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
            MPI_Send(recieved,recv_count,MPI_LONG,(rank+1)%numproc,99,MPI_COMM_WORLD);
        }
        recieved = malloc(sizeof(long) * recv_count);
        MPI_Recv(recieved,recv_count,MPI_LONG,((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
        #pragma omp parallel for
        for(int i = 0 ; i<recv_count;i++){
            save_data[i] = recieved[i];
        }
    }
}


void scatter_colum(Matrix A,long* save_data,int recv_count){
    MPI_Status status;

    if(rank == 0){
        Matrix t_A = transpose(A);
        long* linearize_colonne = lineariser(t_A);
        #pragma omp parallel for
        for(int i_colone = 0 ; i_colone<recv_count;i_colone++){
            save_data[i_colone] = linearize_colonne[i_colone];
        }

        for(int i =1; i<numproc;i++){
            MPI_Send(linearize_colonne+(int_size_of_mat-(r*A->nb_colon)*(i)),r*A->nb_colon,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
        }
    }

    else{
        long* recieved;
        for(int i = 0; i<numproc-rank-1;i++){
            recieved = malloc(sizeof(long) * recv_count);
            MPI_Recv(recieved,recv_count,MPI_LONG,((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
            MPI_Send(recieved,recv_count,MPI_LONG,(rank+1)%numproc,99,MPI_COMM_WORLD);
        }
        recieved = malloc(sizeof(long) * recv_count);
        MPI_Recv(recieved,recv_count,MPI_LONG,((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
        #pragma omp parallel for
        for(int i = 0 ; i<recv_count;i++){
            save_data[i] = recieved[i];
        }
    }
}


void broadcast(int* size_of_mat){
    MPI_Status status;
    if(rank==0){
        for(int i = 0;i<numproc-1;i++){
            MPI_Send(size_of_mat,1,MPI_INT,rank+1,99,MPI_COMM_WORLD);
        }
    }
    else{
        for(int i = 0; i<numproc-rank-1;i++){
            MPI_Recv(&int_size_of_mat,1,MPI_INT,((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
            MPI_Send(&int_size_of_mat,1,MPI_INT,(rank+1)%numproc,99,MPI_COMM_WORLD);
        }
        MPI_Recv(&int_size_of_mat,1,MPI_INT, ((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
    }
}


Matrix create_matrix_from_table(long *tab){

    Matrix to_ret = create_matrix(N);
    #pragma omp parallel for
    for(int i = 0;i<N;i++){
        #pragma omp parallel for
        for(int j = 0;j<N;j++){
            to_ret->tab[i][j] = tab[i*N+j];
        }
    }

    return transpose(to_ret);
}


void circuler(long** r_line_tmp, int *size){
    MPI_Status status;
    int number;//nb d'entiers que l'on doit faire circuler
    long* r_line = *r_line_tmp;
    long *recv;
    if(rank%2 == 0){
        MPI_Send(r_line,*size,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
        MPI_Probe(((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
        MPI_Get_count(&status, MPI_LONG, &number); // et qu'on découvre dynamiquement
        recv = malloc(sizeof(long) * number);
        MPI_Recv(recv,number,MPI_LONG, ((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
    }
    else{
        MPI_Probe(((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
        MPI_Get_count(&status, MPI_LONG, &number);
        recv = malloc(sizeof(long) * number);
        MPI_Recv(recv,number,MPI_LONG, ((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
        MPI_Send(r_line,*size,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
    }
    *r_line_tmp = recv;
    *size = number;
}

void gather(long* data,int nb_to_send,long* save_into,int nb_to_recv){
    MPI_Status status;

    if(rank == 0){

        int i,k;
        #pragma omp parallel for
        for(k = 0; k<nb_to_send; k++){ // save for process 0
            save_into[k] = data[k];
        }

        // recieve elements starting with the last one, and we put it at the end.
        for(i =0 ; i < numproc-1 ; i++){
            MPI_Recv(save_into+int_size_of_mat-(nb_to_recv*(i+1)),nb_to_recv,MPI_LONG, numproc-1,99,MPI_COMM_WORLD,&status);
        }
    }

    else{
        long* received;
        MPI_Send(data,nb_to_send,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
        for(int i = 0; i<rank-1;i++){
            received = malloc(sizeof(long) * nb_to_recv);
            MPI_Recv(received,nb_to_recv,MPI_LONG,((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
            MPI_Send(received,nb_to_send,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
        }
    }
}


void print(Matrix mat){
    for(int i = 0;i<mat->nb_line;i++){
        for (int j = 0; j < mat->nb_colon; j++) {
            if(mat->tab[i][j]==INFINI)printf("i ");
            else printf("%ld ",mat->tab[i][j]);
        }
        printf("\n");
    }
}


long minimum(long a, long b){
    return a<b?a:b;
}

// compute and put the result, of matrix product of 2 matrix of (r,N) and (N,r), inside the matrix N_r_matrix that has N lines and r column
void calcule_r_line_r_colon(long *save_line,long * save_column,
                            long ** N_r_matrix,int r_x_nb_lines, int r_colone, int start){

    int r_line = r_x_nb_lines/N;

    for(int i = 0;i<r_line;i++){
        for(int j = 0;j<r_colone;j++){
            long cur_min = INFINI;

            for(int k=0; k<N; k++){

                long line_value = save_line[i*N+k];
                long column_value = save_column[j*N+k];

                long total;

                if(line_value == INFINI || column_value == INFINI){
                    total = INFINI;
                }
                else{
                    total = line_value + column_value;
                }
                cur_min = minimum(total, cur_min);
            }

            N_r_matrix[start][j] = cur_min; // we computed start previously so we know where to start from
        }
        start++; 
    }
}


Matrix compute_and_get_w_from(Matrix mat){
    Matrix to_ret = create_matrix(mat->nb_line);

    #pragma omp parallel for
    for(int i = 0;i<mat->nb_line;i++){
        #pragma omp parallel for
        for(int j = 0;j<mat->nb_colon;j++){
            if(i == j){
                to_ret->tab[i][j] = 0;
            }
            else if(mat->tab[i][j] > 0){
                to_ret->tab[i][j] = mat->tab[i][j];
            }
            else{
                to_ret->tab[i][j] = INFINI;
            }
        }
    }

    return to_ret;
}


// this function will return the nb of line at which each process will start knowing it's rank and r, r = nb_lines/numproc
int get_start_line(int r_lines){
    if(rank == 0)return 0;

    else{
        int rest = N%numproc;
        int biggest_r = N/numproc+rest;
        return biggest_r+(rank-1)*r_lines;
    }
}


void initialization(char** argv,Matrix* address_A, Matrix* address_w, int* add_r_x_nb_lines){

    int r_x_nb_lines = *add_r_x_nb_lines;
    Matrix A = *address_A;
    A = read_file(argv);
    int_size_of_mat = A->nb_elements;
    r = A->nb_line/numproc;

    *address_w = compute_and_get_w_from(A);

    if(A->nb_line % numproc == 0){
        r_x_nb_lines = int_size_of_mat/numproc;
    }
    else{
        int rest = A->nb_line%numproc;
        r_x_nb_lines = ( (r+rest)*A->nb_colon );
    }
    *address_A = A;
    *add_r_x_nb_lines = r_x_nb_lines;
}


// found on stackOverFlow
int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

void floyd_workshall(char**argv){

    struct timeval start_ticking, end_ticking;

    Matrix A = NULL;
    Matrix w = NULL;
    int r_x_nb_lines;
    int q,mult_of_lines;
    if(rank == 0){
        initialization(argv,&A,&w,&r_x_nb_lines);
    }

    broadcast(&int_size_of_mat);// all process should know what is the size of the matrix
    N = sqrt(int_size_of_mat);
    q=int_size_of_mat/numproc; // number of elements divided by nb of process so we can know how much elements
                                // each process should have EVENLY except of process 0 whose going to take the rest
                                // in addition.
    mult_of_lines=q/N; // this variable is going to determine how many lines all the other process will have
                        // they should all have the SAME number of lines
                        // only process 0 will have the rest of euclidean division between numproc and number of lines
                        // of the matrix

    if(rank!=0){
        r_x_nb_lines =mult_of_lines * N;
    }

    Matrix w_power_i = NULL;
    int r_colone = r_x_nb_lines/N;

    int save_r_x_nb_lines = r_x_nb_lines; // r x nb_lines is nb of elements in each r line.
                                            // it's different for the process 0 if N is not a multiple of N
                                            // as the process 0 will take the rest of remaining lines
                                            // result of euclidean division between nb of line of the matrix and
                                            // nb of process.
                                            // this variable will be update at each circulation. Therefore the
                                            // other process will have this huge number of lines to compute it with
                                            // their column. So we need to save so we can know at the end of circualtion
                                            // of each iterations what is the exact size of the result table, that
                                            // each process has computed

    long* save_line = malloc(sizeof(long)*r_x_nb_lines);
    long* save_colum = malloc(sizeof(long)*r_x_nb_lines);
    scatter_colum(w, save_colum, r_x_nb_lines); // all process will have the same line as we're going to compute w^i x w

    // loop to compute w^N.
    for(int i = 0;i<N;i++){
        if(rank == 0){
            //we take the current time and store it in start
            gettimeofday(&start_ticking, NULL);
        }

        scatter_line(w, save_line, r_x_nb_lines); // at each iteration we scatter/gather the new matrix w^i
        long** N_r_matrix = malloc(sizeof(long*) * N); // Matrix in which we're going to store the result of
                                                            // computation for each r line and r column

        #pragma omp parallel for
        for(int count_r = 0;count_r<N;count_r++){
            N_r_matrix[count_r] = malloc(sizeof(long) * r_colone);
        }

        int start = get_start_line(save_r_x_nb_lines/N);// get the start line for each process

        // this loop is going to compute for eache process the r lines and r columns and place it
        // inside the N x r matrix (for each process) we will gather those N x r matrix in the gather function
        for(int i = 0;i<numproc;i++){
            calcule_r_line_r_colon(save_line, save_colum, N_r_matrix,r_x_nb_lines, r_colone,start);

            circuler(&save_line,&r_x_nb_lines);// we change the size of matrix
                                                        // therefore we update the variable r_x_nb_lines
                                                        // by giving its address to the function circuler.

            start = mod(start-r_x_nb_lines/N,N); // start at the new line given by the number of r lines we have
                                                    // after each circulation

        }

        int nb_to_send = save_r_x_nb_lines;
        int nb_to_recv;
        if (rank == 0)nb_to_recv = mult_of_lines * N; // for gathering the process 0 must know how many r line
                                                        // each process has. because we're using a ring structure
                                                        // for transferring data between process.
        else nb_to_recv = save_r_x_nb_lines;
        long* recieve_the_mat = malloc(sizeof(long) * int_size_of_mat);

        long* linearisedN_x_rMatrix = lineariser_column(N_r_matrix, save_r_x_nb_lines);
        gather(linearisedN_x_rMatrix,nb_to_send,recieve_the_mat,nb_to_recv); // each process has computed it's N x r column
                                                                                // we're gathering inside p0 as table 1D
        //free the memory at each iteration.
        free(linearisedN_x_rMatrix);
        w_power_i = create_matrix_from_table(recieve_the_mat); // transforming 1D table into a matrix.
        free(recieve_the_mat);
        w = w_power_i;

        if(rank == 0){
            //we  store the current time in end
            gettimeofday(&end_ticking, NULL);

            // Uncomment the bellow printf(...) to see the time of each iteration of the process P0.

            //timeval is a struct with 2 parts for time, one in seconds and the other in
            //microseconds. So we convert everything to microseconds before computing
            //the elapsed time
            //printf("time = %ld\n", ((end_ticking.tv_sec * 1000000 + end_ticking.tv_usec)
            //               - (start_ticking.tv_sec * 1000000 + start_ticking.tv_usec)));

            #pragma omp parallel for
            for(int count_r = 0;count_r<N;count_r++){
                free(N_r_matrix[count_r]);
            }
            free(N_r_matrix);
        }
    }

    if(rank == 0){ // we free the heap

        print(w);
        destroy_matrix(&A);
        destroy_matrix(&w);

        free(save_colum);
        free(save_line);

    }

}


int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    floyd_workshall(argv);

    MPI_Finalize();

    return 0;
}
