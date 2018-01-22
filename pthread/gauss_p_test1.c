/* 
 * Original author:  UNKNOWN
 *
 * Modified:         XX XX (January 2018)
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


int task_num;
pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
int nsize = 0;


/* Solve the equation:
 *   matrix * X = R 
 */

double **matrix, *X, *R;

/* Pre-set solution. */

double *X__;




void barrier(int expect)
{
    static int arrived = 0;

    pthread_mutex_lock (&mut);  //lock

    arrived++;
    if (arrived < expect)
        pthread_cond_wait (&cond, &mut);
    else {
        arrived = 0;        // reset the barrier before broadcast is important
        pthread_cond_broadcast (&cond);
    }

    pthread_mutex_unlock (&mut);    //unlock
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


/*
void *task1(void *id){
    int task_id = *(int*)id;
    int l,i,j,k;
    double pivotval;

    int beg, end;


    for (l = 0; l < nsize; l++) {
        barrier(task_num); //1

        if (pivotrow != l) {
            beg = getBeg(task_id, nsize-l);
            end = getEnd(task_id, nsize-l);
            
            for (i = beg+l; i <= end+l; i++) {
                SWAP(matrix[pivotrow][i],matrix[l][i]);
            }

            barrier(task_num); //2

        }

        if (matrix[l][l] != 1.0 ) {

            //current row /= pivot
            beg = getBeg(task_id, nsize-l-1);
            end = getEnd(task_id, nsize-l-1);

            for (j = beg+l+1; j <= end+l+1; j++) {  //para 2
                matrix[l][j] /= matrix[l][l];
            }
            barrier(task_num); //3
            barrier(task_num); //4
        }  

        beg = getBeg(task_id, nsize-l-1);
        end = getEnd(task_id, nsize-l-1);
//        printf("111 beg  %d   end  %d\n", beg, end);

        for (j = beg+l+1; j <= end+l+1; j++) {    //para 3
            pivotval = matrix[j][l];
            matrix[j][l] = 0.0;
             //
            for (k = l + 1; k < nsize; k++) {   
                matrix[j][k] -= pivotval * matrix[l][k];
            }
            R[j] -= pivotval * R[l];
        }  


        barrier(task_num);//3
    }

//    for ( i=nsize-1;i>=0;i--){       //solve
        //X[i] = R[i];
        //-= with solved s[i] 
        barrier(task_num); //6

        beg = getBeg(task_id, i);
        end = getEnd(task_id, i);
        for (int j=beg;j<=end;j++){        //para 4
            R[j] -= X[i] * matrix[j][i];  
 //       }
        barrier(task_num); //6

    } 

}

*/

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
int *max;
void printMax(){
    for(int i=0;i<task_num;i++){        
        printf("lai %d ",max[i]);
    }
    printf("\n");
}


void *task0(void *id){
    int task_id = *(int*)id;

    struct timeval start, finish;
    if (task_id == 0) {   
        gettimeofday(&start, 0);
    }
    int maxId,i, j, k, l, beg1, end1,beg2, end2;
    static int pivotrow;

//printf("%d id \n",task_id);
    double pivotval;
    for (l = 0; l < nsize; l++) {


        beg1 = getBeg(task_id, nsize-l);
        end1 = getEnd(task_id, nsize-l);
   //   printf("000 beg  %d   end  %d\n", beg1, end1);

        max[task_id] = beg1+l;
        for (i = beg1+l; i < end1+l; i++) {
            if (fabs(matrix[i][l]) > fabs(matrix[max[task_id]][l])) {
                max[task_id] = i;
            }
        }
        barrier(task_num); //1
        if (task_id == 0){
            maxId = 0;
            for (i=1;i<task_num;i++)
                if (fabs(matrix[max[i]][l]) > fabs(matrix[max[maxId]][l])) {
                    maxId = i;
                }
            pivotrow = max[maxId];
        }
       
        barrier(task_num); //1
        if (fabs(matrix[pivotrow][l]) == 0.0) {
            fprintf(stderr, "The matrix is singular\n");
            exit(-1);
        }
        
      //  barrier(task_num); //1
    
        if (pivotrow != l) {   //para 1

            for (i = beg1+l; i <= end1+l; i++) {
                SWAP(matrix[pivotrow][i],matrix[l][i]);
            }
            if (task_id == 0) {   
                SWAP(R[pivotrow],R[l]);
            }
            barrier(task_num); //2


        }
        beg2 = getBeg(task_id, nsize-l-1);
        end2 = getEnd(task_id, nsize-l-1);
        // Scale the main row. 
        //factor to 1
        //current row /= pivot
        pivotval = matrix[l][l];
        if (matrix[l][l] != 1.0) {
          //  barrier(task_num); //3


            for (j = beg2+l+1; j <= end2+l+1; j++) {  //para 2

                matrix[l][j] /= pivotval;
            }
            if (task_id == 0) {     //seq
                R[l] /= pivotval;
            }
            barrier(task_num); //3
            
            if (task_id == task_num -1) {   //seq
                matrix[l][l] = 1.0;
            }
        } 



        for (j = beg2+l+1; j <= end2+l+1; j++) {    //para 3
            pivotval = matrix[j][l];
            matrix[j][l] = 0.0;
             //
            for (k = l + 1; k < nsize; k++) {   
                matrix[j][k] -= pivotval * matrix[l][k];
            }
            R[j] -= pivotval * R[l];
        } 
      //  printf("\n");
   //   barrier(task_num); //3
    }

 /*   for ( i=nsize-1;i>=0;i--){       //solve
        X[i] = R[i];
        //-= with solved s[i] 
        barrier(task_num); //6
        barrier(task_num); //6

    }  */

    if (task_id == 0) {
        gettimeofday(&finish, 0);
        fprintf(stdout, "Time:  %f seconds\n", (finish.tv_sec - start.tv_sec) + (finish.tv_usec - start.tv_usec)*0.000001);
    }
}
void printM(){
    for(int i=0;i<nsize;i++){
        for(int j=0;j<nsize;j++){
            printf("%.0f ",matrix[i][j]);
        }
        printf("\n");

    }
}
void printX(){
    for(int i=0;i<nsize;i++){        
        printf("%.0f ",X[i]);
    }
    printf("\n");
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

void computeGauss(int nsize)
{
    int i, j, k;
    double pivotval;
    
    for (i = 0; i < nsize; i++) {
    getPivot(nsize,i);
        
    /* Scale the main row. */
        pivotval = matrix[i][i];
    if (pivotval != 1.0) {
        matrix[i][i] = 1.0;
        for (j = i + 1; j < nsize; j++) {
        matrix[i][j] /= pivotval;
        }
        R[i] /= pivotval;
    }
        
    /* Factorize the rest of the matrix. */
        for (j = i + 1; j < nsize; j++) {
            pivotval = matrix[j][i];
            matrix[j][i] = 0.0;
            for (k = i + 1; k < nsize; k++) {
                matrix[j][k] -= pivotval * matrix[i][k];
            }
            R[j] -= pivotval * R[i];
        }
    }
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
    double error;
    
    if (argc < 2) {
    fprintf(stderr, "usage: %s <matrixfile>\n", argv[0]);
    exit(-1);
    }

    //init mX = b
    nsize = initMatrix(argv[1]);
    initRHS(nsize);
    initResult(nsize);
    max = (int*)malloc(nsize*sizeof(int));


//printM();
    //pthread num and init
    if (argv[2] != NULL)
        task_num = atoi(argv[2]);
    else
        task_num = 1;
    pthread_t *tid;
    int *id;

    id = (int*)malloc(task_num*sizeof(int));
    tid = (pthread_t*)malloc(task_num*sizeof(pthread_t));

    //if thread init fail
    if (!id || !tid)
        printf("fail to create thread\n");

    printf("the num of tasks: %d\n",task_num);

    //create threads and assign tasks
//  id[task_num-1] = task_num-1;
   
 //   pthread_create(&tid[task_num-1], NULL, task0, &id[task_num-1]);

    for (i=0;i<=task_num-1;i++){
        id[i] = i;
        pthread_create(&tid[i], NULL, task0, &id[i]);
    }
    for (i=0;i<task_num;i++){
        pthread_join (tid[i], NULL);
    }
   


//    computeGauss(nsize);

    solveGauss(nsize);
 //printM();
 //printX();

   

    error = 0.0;
    for (i = 0; i < nsize; i++) {
    double error__ = (X__[i]==0.0) ? 1.0 : fabs((X[i]-X__[i])/X__[i]);
    if (error < error__) {
        error = error__;
    }
    }
    fprintf(stdout, "Error: %e\n", error);


    free(id);
    free(tid);
    free(max);
    return 0;
}
