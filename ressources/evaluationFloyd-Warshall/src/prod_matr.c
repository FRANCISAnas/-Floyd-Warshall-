//
// Created by Anas Francis on 04/05/2021.
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include<time.h>
#include <string.h>
#include <limits.h>

int rank; int numproc;
long ** result_scatter;
int r;
int int_size_of_mat;

#define INFINI LONG_MAX
#define ALLOC_SIZE 4
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
    int k = 0;
    long* table = malloc(sizeof(long)*A->nb_line*A->nb_colon);
    for(int i = 0;i<A->nb_line;i++){
        for(int j = 0; j<A->nb_colon; j++){
            table[k++] = A->tab[i][j];
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
        int i = 0;
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
        int k = 0, i_colone = 0;
        for(i_colone = 0 ; i_colone<recv_count;i_colone++){
            recv_data_colum[k++] = linearize_colonne[i_colone];
        }
        if(A->nb_line%numproc != 0){
            for(int i =numproc; i>=2;i--){
                MPI_Send(linearize_colonne+(A->nb_colon*r*(i)),r*A->nb_colon,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
            }
        }
        else{
            for(int i =numproc; i>=2;i--){
                MPI_Send(linearize_colonne+(A->nb_colon*r*(i-1)),r*A->nb_colon,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
            }
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


void scatter_line(Matrix A,int send_count,long* recv_data_line,int recv_count){
    MPI_Status status;
    if(rank == 0){ // si on est emetteur
        /*int rest = A->nb_line%numproc;
        int size_of_recv_data=recv_count + rest*A->nb_colon;
        recv_data = malloc(sizeof(long) * size_of_recv_data);
     */
        long* linearize_lines = lineariser(A);
        int k = 0, i_lines = 0;
        for(i_lines = 0 ; i_lines<recv_count;i_lines++){
            recv_data_line[k++] = linearize_lines[i_lines];
        }
        /*
        printf("linearize_lines = \n");
        for(int i = 0;i<A->nb_line * A->nb_colon; i++){
            printf("%ld ",linearize_lines[i]);
        }
        printf("\n");

        printf("r = %d\n", r);*/
        if(A->nb_line%numproc != 0){
            for(int i =numproc; i>=2;i--){
                MPI_Send(linearize_lines+(A->nb_colon*r*(i)),r*A->nb_colon,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
            }
        }
        else{
            for(int i =numproc; i>=2;i--){
                MPI_Send(linearize_lines+(A->nb_colon*r*(i-1)),r*A->nb_colon,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
            }
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

void reverse_string(char** string, int len){
    char* str = *string;
    char temp;
    for(int i = 0; i<len/2;i++){
        temp = str[i];
        str[i] = str[len-i-1];
        str[len-i-1] = temp;
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

Matrix create_matrix_from_table(long *tab, int size){
    int nb_line = sqrt(size);
    Matrix to_ret = create_matrix(nb_line);
    int k = 0;
    for(int i = 0;i<nb_line;i++){
        for(int j = 0;j<nb_line;j++){
            to_ret->tab[i][j] = tab[k++];
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


void gather(long* send_data,int send_count,long* recv_data,int recv_count,int int_size_of_mat){
    MPI_Status status;
    //printf("recv_count = %d\n",recv_count);
    if(rank == 0){
        int i,k=0;
        for(k = 0; k<send_count; k++){ // save for process 0
            //printf("%ld ",send_data[k]);
            recv_data[k] = send_data[k];
        }
        printf("\n");
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
        printf("\t");
        for (int j = 0; j < mat->nb_colon; j++) {
            printf("%ld ",mat->tab[i][j]);
        }
        printf("\n");
    }
}
long minimum(long a, long b){
    return a<b?a:b;
}

void calcule(long *save_line,long * save_colum, long * save_resulte_for_proc_i, int recv_count, int nb_circulations){

    //printf("nb_circulations = %d\n",nb_circulations);
    int one_line_colum = sqrt(int_size_of_mat);
    int r_local = recv_count/one_line_colum;
    for(int i = 0;i<r_local;i++){
        for(int j = 0;j<r_local;j++){
            long cur_min = INFINI;
            int pos = (((rank-nb_circulations) * r_local)+one_line_colum+j)%one_line_colum;
            //printf("pos = %d, rank = %d\n",pos,rank);
            for(long k=0; k<one_line_colum; k++){
                long line_value = save_line[i*one_line_colum+k];
                long colum_value = save_colum[j*one_line_colum+k];

                long total;
                if(line_value == INFINI || colum_value == INFINI){
                    total = INFINI;
                }
                else{
                    total = line_value + colum_value;
                }
                cur_min = minimum(total, cur_min);

            }
            //printf("%ld ",cur_min);
            save_resulte_for_proc_i[(one_line_colum*i+pos)] = cur_min;
        }
    }
}


Matrix compute_next_rank_matrix(Matrix prev, int rank_of_prev_mat){
    Matrix to_ret = create_matrix(prev->nb_line);
    for(int i = 0;i<prev->nb_line;i++){
        for(int j = 0;j<prev->nb_colon;j++){
            if(i == rank_of_prev_mat + 1 || j == rank_of_prev_mat + 1){
                to_ret->tab[i][j] = prev->tab[i][j];
            }
            else {
                to_ret->tab[i][j] = minimum(prev->tab[i][j], prev->tab[i][rank_of_prev_mat]+prev->tab[rank_of_prev_mat][j]);
            }
        }
    }
    return to_ret;
}


Matrix compute_from(Matrix mat){
    Matrix to_ret = create_matrix(mat->nb_line);

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

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Matrix A = NULL;
    Matrix w = NULL;
    int recieved_count;
    int q,mult_of_lines;
    if(rank == 0){
        A = read_file(argv);
        int_size_of_mat = A->nb_elements;
        r = A->nb_line/numproc;

        w = compute_from(A);
        printf("w = \n");
        print(w);
        if(int_size_of_mat % numproc == 0){
            recieved_count = int_size_of_mat/numproc;
        }
        else{
            int rest = A->nb_line%numproc;
            recieved_count = ( ((A->nb_line/numproc)+rest)*A->nb_colon );
        }
    }
    broadcast(&int_size_of_mat);
    if(rank!=0){
        q=int_size_of_mat/numproc;
        mult_of_lines=q/sqrt(int_size_of_mat);
        recieved_count =mult_of_lines * sqrt(int_size_of_mat);
    }
    //broadcast(0, t_m3x3_B);
    //long* save = malloc(sizeof(long) * r *m4x4_A->nb_colon);
    /* int nb_line_colon = sqrt(int_size_of_mat);
     int rest = (int)nb_line_colon/numproc;
     recieved_count = 2*nb_line_colon * rest;
  }
  printf("recieved_count = %d\n",recieved_count);*/
    long* save_line = malloc(sizeof(long)*recieved_count);
    long* save_colum = malloc(sizeof(long)*recieved_count);
    scatter_line(w,0,save_line,recieved_count);
    scatter_colum(w,0,save_colum,recieved_count);
    long* save_resulte_for_proc_i = malloc(sizeof(long) * recieved_count);
    int save_recieved_count = recieved_count;
    int counter_of_circulation = 0;
    // broadcast line and colums of matrix A and transpos of matrix B
    for(int i = 0;i<numproc;i++){
        calcule(save_line, save_colum, save_resulte_for_proc_i,save_recieved_count,counter_of_circulation);
        int j = 0;
        if(i == 0){
            /*printf("i'm the process %d and i've recieved the following table:\n", rank);
            j = 0;

            while(j<recieved_count){
              printf("%ld ",save_line[j++]);
            }
            printf("\n");
            printf("\n");*/
        }

        circuler(&save_line,&recieved_count);

        counter_of_circulation++;

    }

    printf("i'm the process %d and i've computed the following table:\n", rank);
    for(int k = 0;k<recieved_count;k++){
        printf("%ld ",save_resulte_for_proc_i[k]);
    }
    printf("\n");
    printf("\n");

    int send_count = save_recieved_count;
    int recv_count;
    recv_count = save_recieved_count;
    long* recieve_the_mat = malloc(sizeof(long) * int_size_of_mat);
    gather(save_resulte_for_proc_i,send_count,recieve_the_mat,recv_count,int_size_of_mat);
    Matrix w_i = create_matrix_from_table(recieve_the_mat, int_size_of_mat);
    if(rank == 0){
        printf("w^i = \n");
        print(w_i);
    }
    MPI_Finalize();
    return 0;
}


