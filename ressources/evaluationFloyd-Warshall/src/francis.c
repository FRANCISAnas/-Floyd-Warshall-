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

Matrix create_matrix(int alloc_size);

void destroy_matrix(Matrix *table);

void print(Matrix mat);

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


void scatter_colum(Matrix A,int send_count,long* recv_data_colum,int recv_count){
    MPI_Status status;

    if(rank == 0){
        Matrix t_A = transpose(A);
        long* linearize_colonne = lineariser(t_A);
        #pragma omp parallel for
        for(int i_colone = 0 ; i_colone<recv_count;i_colone++){
            recv_data_colum[i_colone] = linearize_colonne[i_colone];
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
            recv_data_colum[i] = recieved[i];
        }
    }
}

long* lineariser_column(long** table, int recv_count){
    long* to_ret = malloc(sizeof(long)* recv_count);
    int r_local = recv_count/N;
    #pragma omp parallel for
    for(int i = 0;i<r_local;i++){
        #pragma omp parallel for
        for(int j = 0;j<N;j++){
            to_ret[j+i*N] = table[j][i];
        }
    }
    return to_ret;
}

void scatter_line(Matrix A,int send_count,long* recv_data_line,int recv_count){
    MPI_Status status;
    if(rank == 0){ // si on est emetteur

        long* linearize_lines = lineariser(A);
        #pragma omp parallel for
        for(int i_lines = 0 ; i_lines<recv_count;i_lines++){
            recv_data_line[i_lines] = linearize_lines[i_lines];
        }

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
            recv_data_line[i] = recieved[i];
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
    int number;//nb d'entiers que l'on doit faire circuler
    long* r_line = *r_line_tmp;
    MPI_Status status;
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



void gather(long* send_data,int send_count,long* recv_data,int recv_count){
    MPI_Status status;

    if(rank == 0){

        int i,k;
        #pragma omp parallel for
        for(k = 0; k<send_count; k++){ // save for process 0
            recv_data[k] = send_data[k];
        }
        for(i =0 ; i < numproc-1 ; i++){
            MPI_Recv(recv_data+int_size_of_mat-(recv_count*(i+1)),recv_count,MPI_LONG, numproc-1,99,MPI_COMM_WORLD,&status);
        }
    }

    else{

        long* received;
        MPI_Send(send_data,send_count,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
        for(int i = 0; i<rank-1;i++){
            received = malloc(sizeof(long) * recv_count);
            MPI_Recv(received,recv_count,MPI_LONG,((rank-1)+numproc)%numproc,99,MPI_COMM_WORLD,&status);
            MPI_Send(received,send_count,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
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



void calcule_r_line_r_colon(long *save_line,long * save_colum, long ** N_r_matrix,int recv_count, int r_colone, int start){

    int r_local = recv_count/N;

    for(int i = 0;i<r_local;i++){
        for(int j = 0;j<r_colone;j++){
            long cur_min = INFINI;

            for(int k=0; k<N; k++){

                long line_value = save_line[i*N+k];
                long colum_value = save_colum[j*N+k];

                long total;
                if(line_value == INFINI || colum_value == INFINI){
                    total = INFINI;
                }
                else{
                    total = line_value + colum_value;
                }
                cur_min = minimum(total, cur_min);
            }

            N_r_matrix[start][j] = cur_min;
        }
        start++;
    }

}


Matrix compute_from(Matrix mat){
    Matrix to_ret = create_matrix(mat->nb_line);

    #pragma omp parallel for
    for(int i = 0;i<mat->nb_line;i++){
        #pragma omp parallel for
        for(int j = 0;j<mat->nb_colon;j++){
            if(i == j && mat->tab[i][j] == 0){
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


int get_start_line(int r_x_N){
    if(rank == 0)return 0;

    else{
        int rest = N%numproc;
        int biggest_r = N/numproc+rest;
        int sum_to_ret = biggest_r;
        sum_to_ret = sum_to_ret+(rank-1)*r_x_N;
        return sum_to_ret;
    }
}

// found on stackOverFlow
int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

void computation(Matrix* address_w, int* address_reieved_count,int mult_of_lines ){
    struct timeval start_ticking, end_ticking;
    int recieved_count = *address_reieved_count;
    Matrix w = *address_w;

    Matrix w_power_i = NULL;
    int r_colone = recieved_count/N;

    int save_recieved_count = recieved_count;
    long* save_line = malloc(sizeof(long)*recieved_count);
    long* save_colum = malloc(sizeof(long)*recieved_count);
    scatter_colum(w,0,save_colum,recieved_count);
    for(int i = 0;i<N;i++){
        if(rank == 0){
            //we take the current time and store it in start
            gettimeofday(&start_ticking, NULL);
        }
        scatter_line(w,0,save_line,recieved_count);
        long** N_r_matrix = malloc(sizeof(long*) * N);
        #pragma omp parallel for
        for(int count_r = 0;count_r<N;count_r++){
            N_r_matrix[count_r] = malloc(sizeof(long) * r_colone);
        }

        int start = get_start_line(save_recieved_count/N);

        for(int i = 0;i<numproc;i++){
            calcule_r_line_r_colon(save_line, save_colum, N_r_matrix,recieved_count, r_colone,start);

            circuler(&save_line,address_reieved_count);

            recieved_count = *address_reieved_count;

            start = mod(start-recieved_count/N,N);

        }

        recieved_count = save_recieved_count;

        int send_count = recieved_count;
        int recv_count;
        if (rank == 0)recv_count = mult_of_lines * N;
        else recv_count = recieved_count;
        long* recieve_the_mat = malloc(sizeof(long) * int_size_of_mat);

        long* linearisedN_x_rMatrix = lineariser_column(N_r_matrix, recieved_count);
        gather(linearisedN_x_rMatrix,send_count,recieve_the_mat,recv_count);

        w_power_i = create_matrix_from_table(recieve_the_mat);
        //free(recieve_the_mat);
        w = w_power_i;
        if(rank == 0){
            //we  store the current time in end
            gettimeofday(&end_ticking, NULL);


            //timeval is a struct with 2 parts for time, one in seconds and the other in
            //microseconds. So we convert everything to microseconds before computing
            //the elapsed time
            //printf("time = %ld\n", ((end_ticking.tv_sec * 1000000 + end_ticking.tv_usec)
              //               - (start_ticking.tv_sec * 1000000 + start_ticking.tv_usec)));
        }
    }

    *address_w = w;
    if(rank == 0){
        free(save_colum);
        free(save_line);
        /*#pragma omp parallel for
        for(int count_r = 0;count_r<N;count_r++){
             free(N_r_matrix[count_r]);
        }
        free(N_r_matrix);*/
        //destroy_matrix(&w_power_i);
    }

}


void initialization(char** argv,Matrix* address_A, Matrix* address_w, int* add_recieved_count){

    int recieved_count = *add_recieved_count;
    Matrix A = *address_A;
    A = read_file(argv);
    int_size_of_mat = A->nb_elements;
    r = A->nb_line/numproc;

    *address_w = compute_from(A);

    if(A->nb_line % numproc == 0){
        recieved_count = int_size_of_mat/numproc;
    }
    else{
        int rest = A->nb_line%numproc;
        recieved_count = ( (r+rest)*A->nb_colon );
    }
    *address_A = A;
    *add_recieved_count = recieved_count;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Matrix A = NULL;
    Matrix w = NULL;
    int recieved_count;
    int q,mult_of_lines;
    if(rank == 0){
        initialization(argv,&A,&w,&recieved_count);
    }
    broadcast(&int_size_of_mat);
    N = sqrt(int_size_of_mat);
    q=int_size_of_mat/numproc;
    mult_of_lines=q/N;

    if(rank!=0){
        recieved_count =mult_of_lines * N;
    }

    computation(&w,&recieved_count, mult_of_lines);

    if(rank == 0){
        print(w);
        destroy_matrix(&A);
        destroy_matrix(&w);
    }

    MPI_Finalize();
    return 0;
}
