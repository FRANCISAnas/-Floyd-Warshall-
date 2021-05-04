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


void scatter_r_lines_colonne(Matrix A,int send_count,long* recv_data,int recv_count){


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

        printf("linearize_lines = \n");

        for(int i = 0;i<A->nb_line * A->nb_colon; i++){
            printf("%ld ",linearize_lines[i]);
        }

        printf("\n");
        printf("\n");


        int size = (A->nb_line-(recv_count/2/A->nb_line))*A->nb_colon;
        printf("size = %d\n",size);
        long r_line_r_colonne [size*2];
        int j = 0;
        k = 0;
        printf("i_lines = %d\n",i_lines);
        printf("i_colone = %d\n",i_colone);
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


        printf("\n");
        printf("r_line_r_colonne = \n");
        for(int i = 0;i<k;i++){
            printf("%ld ",r_line_r_colonne[i]);
        }
        printf("\n");
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


void print(Matrix mat){
    for(int i = 0;i<mat->nb_line;i++){
        printf("\t");
        for (int j = 0; j < mat->nb_colon; j++) {
            printf("%ld ",mat->tab[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    time_t t;
    srand((unsigned) time(&t));
    char* size_of_matrix = malloc(sizeof(char)*255);
    int len = strlen(argv[5]);
    int i = 0;
    while(argv[5][len-i-1] != '='){
        size_of_matrix[i]=argv[5][len-i-1];
        i++;
    }
    size_of_matrix[i] = '\0';
    reverse_string(&size_of_matrix,i);
    int int_size_of_mat = atoi(size_of_matrix);
    int recieved_count;
    Matrix m4x4_A=NULL,m4x4_B=NULL;

    if(rank == 0){
        m4x4_A = create_matrix(8);
        m4x4_B = create_matrix(8);


        for(int i = 0; i<m4x4_A->nb_line;i++){
            for(int j = 0;j<m4x4_A->nb_colon;j++){
                m4x4_A->tab[i][j]=rand() % 50;
                m4x4_B->tab[i][j]=rand() % 50;
            }
        }
        r  = m4x4_A->nb_line/numproc;
        printf("A = \n");
        print(m4x4_A);
        if(int_size_of_mat % numproc == 0){
            recieved_count = int_size_of_mat/numproc*2;
        }
        else{
            int rest = m4x4_A->nb_line%numproc;
            printf("rest = %d\n", rest);
            recieved_count = 2 *( ((m4x4_A->nb_line/numproc)+rest)*m4x4_A->nb_colon );
        }
    }
    else{
        int q = int_size_of_mat/numproc;
        int mult_of_lines = q/sqrt(int_size_of_mat);
        recieved_count = mult_of_lines * sqrt(int_size_of_mat) * 2;
    }

    //broadcast(0, t_m3x3_B);*/

    //long* save = malloc(sizeof(long) * r *m4x4_A->nb_colon);


    /* int nb_line_colon = sqrt(int_size_of_mat);
     int rest = (int)nb_line_colon/numproc;
     recieved_count = 2*nb_line_colon * rest;
  }*/

    printf("recieved_count = %d\n",recieved_count);
    long* save_line_colonne = malloc(sizeof(long)*recieved_count);
    scatter_r_lines_colonne(m4x4_A,0,save_line_colonne,recieved_count);
    // broadcast line and colums of matrix A and transpos of matrix B
    int j = 0;
    switch(rank){
        case 0:
            printf("i'm the process %d and i've recieved the following table:\n", rank);
            j = 0;
            while(j<recieved_count){
                printf("%ld ",save_line_colonne[j++]);
            }

            printf("\n");
            printf("\n");
            break;
        case 1:
            printf("i'm the process %d and i've recieved the following table:\n", rank);
            j = 0;
            while(j<recieved_count){
                printf("%ld ",save_line_colonne[j++]);
            }

            printf("\n");
            printf("\n");
            break;
        case 2:
            printf("i'm the process %d and i've recieved the following table:\n", rank);
            j = 0;
            while(j<recieved_count){
                printf("%ld ",save_line_colonne[j++]);
            }

            printf("\n");
            printf("\n");
            break;
        default:
            printf("i'm the process %d and i've recieved the following table:\n", rank);
            j = 0;
            while(j<recieved_count){
                printf("%ld ",save_line_colonne[j++]);
            }
            printf("\n");
            printf("\n");
            break;
    }

    /*
      for(int j = 0; j<r*m4x4_A->nb_colon;j++){
        printf("%ld ", result_scatter[0][j]);
      }
  */


    MPI_Finalize();



    return 0;

}
