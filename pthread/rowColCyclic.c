/* 
 * Original author:  UNKNOWN
 *
 * Modified:         NAME  (January 2018)
 */
#include <pthread.h>

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>

/* #define DEBUG */

#define SWAP(a,b)       {double tmp; tmp = a; a = b; b = tmp;}

/* Solve the equation:
 *   matrix * X = R 
 */


int task_num;
pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
int nsize = 0;
int side;


double **matrix, *X, *R;

/* Pre-set solution. */

double *X__;


void printM(){
    for(int i=0;i<nsize;i++){
        for(int j=0;j<nsize;j++){
            printf("%.0f ",matrix[i][j]);
        }
        printf("\n");

    }
}



/* Initialize the matirx. */

int initMatrix(const char *fname)
{
    FILE *file;
    int l1, l2, l3;
    double d;
    int nsize;
    int i, j;
    double *tmp;
    char buffer[1024];

    if ((file = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "The matrix file open error\n");
        exit(-1);
    }
    
    /* Parse the first line to get the matrix size. */
    fgets(buffer, 1024, file);
    sscanf(buffer, "%d %d %d", &l1, &l2, &l3);
    nsize = l1;
#ifdef DEBUG
    fprintf(stdout, "matrix size is %d\n", nsize);
#endif

    /* Initialize the space and set all elements to zero. */
    matrix = (double**)malloc(nsize*sizeof(double*));
    assert(matrix != NULL);
    tmp = (double*)malloc(nsize*nsize*sizeof(double));
    assert(tmp != NULL);    
    for (i = 0; i < nsize; i++) {
        matrix[i] = tmp;
        tmp = tmp + nsize;
    }
    for (i = 0; i < nsize; i++) {
        for (j = 0; j < nsize; j++) {
            matrix[i][j] = 0.0;
        }
    }

    /* Parse the rest of the input file to fill the matrix. */
    for (;;) {
    fgets(buffer, 1024, file);
    sscanf(buffer, "%d %d %lf", &l1, &l2, &d);
    if (l1 == 0) break;

    matrix[l1-1][l2-1] = d;
#ifdef DEBUG
    fprintf(stdout, "row %d column %d of matrix is %e\n", l1-1, l2-1, matrix[l1-1][l2-1]);
#endif
    }

    fclose(file);
    return nsize;
}

/* Initialize the right-hand-side following the pre-set solution. */

void initRHS(int nsize)
{
    int i, j;

    X__ = (double*)malloc(nsize * sizeof(double));
    assert(X__ != NULL);
    for (i = 0; i < nsize; i++) {
    X__[i] = i+1;
    }

    R = (double*)malloc(nsize * sizeof(double));
    assert(R != NULL);
    for (i = 0; i < nsize; i++) {
    R[i] = 0.0;
    for (j = 0; j < nsize; j++) {
        R[i] += matrix[i][j] * X__[j];
    }
    }
}

/* Initialize the results. */

void initResult(int nsize)
{
    int i;

    X = (double*)malloc(nsize * sizeof(double));
    assert(X != NULL);
    for (i = 0; i < nsize; i++) {
    X[i] = 0.0;
    }
}

/* Get the pivot - the element on column with largest absolute value. */

void getPivot(int nsize, int currow)
{
    int i, pivotrow;

    pivotrow = currow;
    for (i = currow+1; i < nsize; i++) {
    if (fabs(matrix[i][currow]) > fabs(matrix[pivotrow][currow])) {
        pivotrow = i;
    }
    }

    if (fabs(matrix[pivotrow][currow]) == 0.0) {
        fprintf(stderr, "The matrix is singular\n");
        exit(-1);
    }
    
    if (pivotrow != currow) {
#ifdef DEBUG
    fprintf(stdout, "pivot row at step %5d is %5d\n", currow, pivotrow);
#endif
        for (i = currow; i < nsize; i++) {
            SWAP(matrix[pivotrow][i],matrix[currow][i]);
        }
        SWAP(R[pivotrow],R[currow]);
    }
}

/* For all the rows, get the pivot and eliminate all rows and columns
 * for that particular pivot row. */

int getBeg(int id, int task_size){
    int beg;
 //   printf("id %d, size %d\n", id, task_size);
    if (task_size >= task_num)
        beg = (id)*task_size/(task_num);
    else {
        if ( id < task_size )
            return id;
        else
            return 0;
    }
    return beg;
}

int getEnd(int id, int task_size){
    int end;
    if (task_size >= task_num)
        end = (id +1)*task_size/(task_num)-1;
    else {
        if ( id < task_size )
            return id;
        else
            return -1;
    }
    return end;
}


int currentRow;
int beg2,end2;
pthread_t *tid;
int *id;
int *endList;



void* task0(void* id){
    int task_id = *(int*)id;
 //   int beg, end;

    double pivotval;
    int j,k, cur;

    cur = currentRow;
    int side1 = side;
    int  sameCol =  task_id%side1 == cur%side1?1:0; 
//    beg = endList[task_id]+1;
 //   end = endList[task_id+1];
    //printM();

 //   if ( end > cur){

//    if (beg < cur+1)
  //      beg = cur+1;
      //  printf("beg: %d  end :  %d   cur:  %d\n",beg,end,cur);
        
          //  printf("   jjj: %d  \n",task_id/side);
    for (j = task_id/side1; j < nsize; j+=side1) {

        if (j > cur){
            pivotval = matrix[j][cur];
           // matrix[j][cur] = 0.0;
            for (k = task_id%side1; k < nsize; k+=side1) {
                if (k >cur){
                    matrix[j][k] -= pivotval * matrix[cur][k];
                }
            }
            if (sameCol){
                R[j] -= pivotval * R[cur];              // printf("%d\n", j );
            }           // if ( task_id == 2 || task_id == 0)
        }
    }

  //  }
 //printf("\n");
}





void computeGauss(int nsize)
{
    int i, j, k;
    double pivotval;
    
    for (i = 0; i < nsize; i++) {
    getPivot(nsize,i);
        
    /* Scale the main row. */
    pivotval = matrix[i][i];
    currentRow = i;
    if (pivotval != 1.0) {
        matrix[i][i] = 1.0;
        for (j = i + 1; j < nsize; j++) {
        matrix[i][j] /= pivotval;
        }
        R[i] /= pivotval;
    }
        



        for (k=0;k<task_num;k++){
            id[k] = k;
            pthread_create(&tid[k], NULL, task0, &id[k]);
          //  beg2 = getBeg(task_id, nsize);
          //  end2 = getEnd(task_id, nsize);    dothis please!~~~
        }

    /* Factorize the rest of the matrix. 
        for (j = i + 1; j < nsize; j++) {
            pivotval = matrix[j][i];
            matrix[j][i] = 0.0;
            for (k = i + 1; k < nsize; k++) {
                matrix[j][k] -= pivotval * matrix[i][k];
            }
            R[j] -= pivotval * R[i];
        }
*/
        for (k=0;k<task_num;k++){
            pthread_join (tid[k], NULL);
        }
    // /printM(); printf("\n");
    } //printf("fail to create thread\n");
}

/* Solve the equation. */

void solveGauss(int nsize)
{
    int i, j;

    X[nsize-1] = R[nsize-1];
    for (i = nsize - 2; i >= 0; i --) {
        X[i] = R[i];
        for (j = nsize - 1; j > i; j--) {
            X[i] -= matrix[i][j] * X[j];
        }
    }

#ifdef DEBUG
    fprintf(stdout, "X = [");
    for (i = 0; i < nsize; i++) {
        fprintf(stdout, "%.6f ", X[i]);
    }
    fprintf(stdout, "];\n");
#endif
}



int main(int argc, char *argv[])
{
    int i;
    struct timeval start, finish;
    double error;
    
    if (argc < 2) {
    fprintf(stderr, "usage: %s <matrixfile>\n", argv[0]);
    exit(-1);
    }

    nsize = initMatrix(argv[1]);
    initRHS(nsize);
    initResult(nsize);



    //pthread num and init

    if (argv[2] != NULL)
        task_num = atoi(argv[2]);
    else
        task_num = 1;

    side = sqrt(task_num);

    id = (int*)malloc(task_num*sizeof(int));
    tid = (pthread_t*)malloc(task_num*sizeof(pthread_t));

    if (!id || !tid)
        printf("fail to create thread\n");

    printf("the num of tasks: %d\n",task_num);

    endList = (int*)malloc(task_num*sizeof(int));
    endList[0] = 0;
    for (i=1;i<=task_num;i++){
        endList[i] = getEnd(i-1, nsize);
        printf("beg: %d", endList[i-1]+1);
        printf("end: %d\n", endList[i]);
        ;
    }


 
  //  printM();
    gettimeofday(&start, 0);
        computeGauss(nsize);
    gettimeofday(&finish, 0);

    solveGauss(nsize);

    //printM();
    fprintf(stdout, "Time:  %f seconds\n", (finish.tv_sec - start.tv_sec) + (finish.tv_usec - start.tv_usec)*0.000001);

    error = 0.0;
    for (i = 0; i < nsize; i++) {
    double error__ = (X__[i]==0.0) ? 1.0 : fabs((X[i]-X__[i])/X__[i]);
    if (error < error__) {
        error = error__;
    }
    }

    printf("li hisshibushi sha\n");
    fprintf(stdout, "Error: %e\n", error);
    free(id);
    free(tid);
    free(endList);

    return 0;
}
