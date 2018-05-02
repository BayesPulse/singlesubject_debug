#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

double M;

int main (int argc, char *argv[]) {

    int v, i;
    unsigned long *seed;
    double **S,  **rwish, det, **rwishinv;
    FILE *fwish,  *fseed;

    int rwishart(double **, double **, int, int, unsigned long *, int);
    int cholesky_decomp(double **, int);

    seed  = (unsigned long *)calloc(3, sizeof(unsigned long));
    fseed = fopen("seed.dat", "r");
    fscanf(fseed, "%u %u %u\n", &seed[0], &seed[1], &seed[2]);
    fclose(fseed);

    fwish = fopen("wishart.dat", "w");
    
    S = (double **)calloc(2, sizeof(double *));
    for (i=0;i<2;i++) {
        S[i] = (double *)calloc(2, sizeof(double));
    }

    rwish = (double **)calloc(2, sizeof(double *));
    for (i=0;i<2;i++) {
        rwish[i] = (double *)calloc(2, sizeof(double));
    }
      
    rwishinv = (double **)calloc(2, sizeof(double *));
    
    M = exp(-32*log(2.0));
    for (i=0;i<2;i++) {
      rwishinv[i] = (double *)calloc(2, sizeof(double));
    }
    S[0][0] = S[1][1] = 1.0/5.0;
    S[1][0] = S[0][1] = 0.;
      
      
    if (!cholesky_decomp(S, 2)){
        printf("not PSD matrix\n");
        exit(0);
    }
     
    for (i=0;i<10000;i++) {
        /* printf("first %lf %lf %lf \n", S[0][0], S[1][1], S[0][1]);*/
        if (!rwishart(rwish, S, 2, 4, seed, 1)) printf("1 wishart\n");
        /* printf("after %lf %lf %lf \n", S[0][0], S[1][1], S[0][1]);*/
        det = rwish[0][0]*rwish[1][1] - rwish[0][1]*rwish[1][0];
        rwishinv[0][0] = rwish[1][1]/det;
        rwishinv[1][1] = rwish[0][0]/det;
        rwishinv[0][1] = rwishinv[1][0] = -rwish[0][1]/det;  
        fprintf(fwish, "%lf %lf %lf %lf %lf %lf\n", rwish[0][0], rwish[1][1],
                rwish[0][1], rwishinv[0][0], rwishinv[1][1], rwishinv[0][1]);
     }
     
    /* save the current random number as the seed for the next simulation */
    fseed = fopen("seed.dat", "w");
    fprintf(fseed, "%u %u %u\n", seed[0], seed[1], seed[2]);
    fclose(fseed);
     
    for (i=0;i<2;i++) {
        free(S[i]);
        free(rwish[i]);
        free(rwishinv[i]);
    }
     
    free(S);
    free(rwish);
    free(rwishinv);
    fclose(fwish);

}


    
