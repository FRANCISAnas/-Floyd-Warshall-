//
// Created by Anas Francis on 04/05/2021.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include<time.h>
#include <string.h>

int rank; int numproc;
long ** result_scatter;
int r;
int int_size_of_mat;

/*
int main(int argc, char *argv[]) {

  int rank, numprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int nb_proc = 4;
  return 0;
}*/

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

void scatter(Matrix A,int send_count,long* recv_data,int recv_count){


    MPI_Status status;

    if(rank == 0){ // si on est emetteur
        /*int rest = A->nb_line%numproc;
        int size_of_recv_data=recv_count + rest*A->nb_colon;
        recv_data = malloc(sizeof(long) * size_of_recv_data);
     */
        long* linearize_lines = lineariser(A);
        Matrix t_A = transpose(A);
        long* linearize_colonne = lineariser(t_A);
        int k = 0, i_lines = 0, i_colone = 0;
        for(i_lines = 0 ; i_lines<recv_count/2;i_lines++){
            recv_data[k++] = linearize_lines[i_lines];
        }
        for(i_colone = 0 ; i_colone<recv_count/2;i_colone++){
            recv_data[k++] = linearize_colonne[i_colone];
        }
/*
        printf("linearize_lines = \n");

        for(int i = 0;i<A->nb_line * A->nb_colon; i++){
            printf("%ld ",linearize_lines[i]);
        }

        printf("\n");
        printf("\n");
*/

        int size = (A->nb_line-(recv_count/2/A->nb_line))*A->nb_colon;
        //printf("size = %d\n",size);
        long r_line_r_colonne [size*2];
        int j = 0;
        k = 0;
        /*printf("i_lines = %d\n",i_lines);
        printf("i_colone = %d\n",i_colone);*/
        for(int i = 0;i<numproc-1;i++){
            j = k/2;
            int save_j = j;
            int counter = 0;

            while(counter<r*A->nb_colon){
                r_line_r_colonne[k++] = linearize_lines[i_lines++];
                counter++;
            }
            j = save_j;
            counter = 0;
            while(counter<r*A->nb_colon){
                r_line_r_colonne[k++] = linearize_colonne[i_colone++];
                counter++;
            }
        }

        /* printf("\n");
         printf("r_line_r_colonne = \n");
         for(int i = 0;i<k;i++){
             printf("%ld ",r_line_r_colonne[i]);
         }
         printf("\n");*/
        for(int i =numproc; i>=2;i--){
            MPI_Send(r_line_r_colonne+(2*A->nb_colon*r*(i-2)),r*A->nb_colon*2,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);
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
            recv_data[i] = recieved[i];
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

void broadcast(int size_of_mat){
    MPI_Status status;
    if(rank==0){
        for(int i = 0;i<numproc-1;i++){
            MPI_Send(&size_of_mat,1,MPI_INT,i+1,99,MPI_COMM_WORLD);
        }
    }
    else{

        MPI_Recv(&int_size_of_mat,1,MPI_INT, 0,99,MPI_COMM_WORLD,&status);
    }
}

void circuler(long* r_line_r_colonne, int size){
    MPI_Status status;

    //MPI_Probe(rank, 99, MPI_COMM_WORLD, &status);

    MPI_Send(r_line_r_colonne,size/2,MPI_LONG, (rank+1)%numproc,99,MPI_COMM_WORLD);

    //MPI_Get_count(&status, MPI_INT, &n_items);


    MPI_Recv(r_line_r_colonne,size/2,MPI_LONG, (rank-1)%numproc,99,MPI_COMM_WORLD,&status);
}


void gather(long* send_data,int send_count,long* recv_data,int recv_count,int int_size_of_mat){

    //printf("i'm proecss %d\n", rank);
    MPI_Status status;
    if(rank == 0){

        int i,k=0;
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


int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Matrix A = NULL;

    int recieved_count;

    int q,mult_of_lines;
    if(rank == 0){

        A = read_file(argv);

        int_size_of_mat = A->nb_elements;

        r = A->nb_line/numproc;
        printf("A = \n");
        print(A);
        if(int_size_of_mat % numproc == 0){
            recieved_count = int_size_of_mat/numproc*2;
        }
        else{
            int rest = A->nb_line%numproc;
            recieved_count = 2 *( ((A->nb_line/numproc)+rest)*A->nb_colon );
        }

    }
    broadcast(int_size_of_mat);
    if(rank!=0){
        q=int_size_of_mat/numproc;
        mult_of_lines=q/sqrt(int_size_of_mat);
        recieved_count =mult_of_lines * sqrt(int_size_of_mat) * 2;
    }

    //broadcast(0, t_m3x3_B);*/

    //long* save = malloc(sizeof(long) * r *m4x4_A->nb_colon);


    /* int nb_line_colon = sqrt(int_size_of_mat);
     int rest = (int)nb_line_colon/numproc;
     recieved_count = 2*nb_line_colon * rest;
  }*/

    printf("recieved_count = %d\n",recieved_count);
    long* save_line_colonne = malloc(sizeof(long)*recieved_count);
    scatter(A,0,save_line_colonne,recieved_count);
    // broadcast line and colums of matrix A and transpos of matrix B
    for(int i = 0;i<numproc;i++){
        int j = 0;
        if(i == 2){
            printf("i'm the process %d and i've recieved the following table:\n", rank);
            j = 0;
            while(j<recieved_count){
                printf("%ld ",save_line_colonne[j++]);
            }

            printf("\n");
            printf("\n");
        }
        circuler(save_line_colonne,recieved_count);

    }


    int send_count = recieved_count/2;
    int recv_count;
    if(rank==0)recv_count =  mult_of_lines * sqrt(int_size_of_mat);
    else recv_count = recieved_count/2;
    long* recieve_the_mat = malloc(sizeof(long) * int_size_of_mat);
    //gather(save_line_colonne,send_count,recieve_the_mat,recv_count,int_size_of_mat);

    if(rank == 0){
        for(int i = 0; i < int_size_of_mat; i++){
            printf("%ld ",recieve_the_mat[i]);
        }
        printf("\n");

        Matrix third_rank = compute_next_rank_matrix(A,1);

        print(third_rank);
    }




    MPI_Finalize();



    return 0;

}
