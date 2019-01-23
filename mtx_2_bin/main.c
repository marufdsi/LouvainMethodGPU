/* 
 *   Matrix Market I/O example program
 *
 *   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
 *   and copies it to stdout.  This porgram does nothing useful, but
 *   illustrates common usage of the Matrix Matrix I/O routines.
 *   (See http://math.nist.gov/MatrixMarket for details.)
 *
 *   Usage:  a.out [filename] > output
 *
 *       
 *   NOTES:
 *
 *   1) Matrix Market files are always 1-based, i.e. the index of the first
 *      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
 *      OFFSETS ACCORDINGLY offsets accordingly when reading and writing 
 *      to files.
 *
 *   2) ANSI C requires one to use the "l" format modifier when reading
 *      double precision floating point numbers in scanf() and
 *      its variants.  For example, use "%lf", "%lg", or "%le"
 *      when reading doubles, otherwise errors will occur.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include "mmio.h"

int main(int argc, char *argv[]) {
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    int i, *I, *J;
    double *val=NULL;

    printf("\nsizeof (unsigned long) = %lu\n", sizeof (unsigned long));

    if (argc < 2) {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    } else {
        if ((f = fopen(argv[1], "r")) == NULL) {
            printf("\n    Can't open file \n");
            exit(1);
        }
    }

    printf("\nhere\n");
    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode)) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
        exit(1);

    assert(M == N);
    unsigned long *degC = (unsigned long*) calloc(M + 1, sizeof (unsigned long));

    int isPattern = mm_is_pattern(matcode);


    //Prefix Sum
    for (i = 1; i <= M; i++) {
        degC[i] = degC[i] + degC[i - 1];
    }
    printf("\n2*nnz = %d, totDeg= %lu\n", 2 * nz, degC[M]);

    /* reseve memory for matrices */

    I = (int *) malloc(nz * sizeof (int));
    J = (int *) malloc(nz * sizeof (int));
    if (!isPattern)
        val = (double *) malloc(nz * sizeof (double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    size_t nrSelfLoop = 0;
    int firstLine=1;
    int nrValid = 0;
    for (i = 0; i < nz; i++) {
        if (!isPattern) {
            fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        } else {
            if(firstLine){printf("\nReading Pattern\n"); firstLine=0;}
            fscanf(f, "%d %d\n", &I[i], &J[i]);
        }
        if(I[i]>M || J[i]>M){
            printf("\n Invalid edge  %d - %d\n", I[i], J[i]); 
        }
        else if (I[i] != J[i]) { //NOTE: Ignore self loop 
            //printf(" Edge %d %d\n", I[i] , J[i] );
            I[i]--; /* adjust from 1-based to 0-based */
            J[i]--;
            degC[ I[i] + 1 ]++;
            degC[ J[i] + 1 ]++; // degC[a+1] =  deg(v_a); degC[1] =deg(v_0)
            nrValid++;
            //if (I[i] == 1) {
                //printf("\nDegree of %d is %d\n", I[i], degC[I[i]]);
            //}
        }else {
            nrSelfLoop++;

        }
    }
    int nrDegZeroV=0, nrDegOneV=0;
    //Prefix Sum
    for (i = 1; i <= M; i++) {

        int a = degC[i];

        degC[i] = degC[i] + degC[i - 1];

        if (a == 0) nrDegZeroV++;

        if (a == 1) nrDegOneV++;
    }
    //NOTE degC[0] = 0; degC[1] = deg(v_0); degC[2] = deg(v_0) + deg(v_1);


    printf("\n2*nnz = %d, totDeg= %lu, nrValid= %d nrDegZeroV= %d, nrDegOneV=%d \n", 2 * nz, degC[M], nrValid, nrDegZeroV, nrDegOneV);

    // allocate memory for compress neighbor list
    int szLinkArray = degC[M];
    unsigned int *links = (unsigned int*) malloc(szLinkArray * sizeof (unsigned int));

    int nrCopied = 0;

    printf("\n ---  > %d %d\n", szLinkArray, 2*nrCopied);
    for (i = 0; i < nz; i++) {

        if (I[i] >=M ||  J[i] >=M) {
            printf("\n Ignoring  edge  %d--%d\n", I[i], J[i]); 
        }
        else if (I[i] != J[i]) {

            unsigned long *ptrR = &degC[I[i]];
            unsigned long *ptrC = &degC[J[i]];
            if(*ptrR>=szLinkArray||*ptrC>=szLinkArray )
                printf("\ncan't happen\n");
            links[*ptrR ] = J[i];
            links[*ptrC] = I[i];

            *ptrC = *ptrC + 1;
            *ptrR = *ptrR + 1;
            nrCopied++;
        }
    }

    printf("\n ---  > %d %d\n", szLinkArray, 2*nrCopied);
    if (f != stdin) fclose(f);
    assert(szLinkArray == 2 * nrCopied);
    // Write to a binary file

    /*
    for (i = 0; i < szLinkArray; i++) {
        printf("%lu ", links[i]);
    }
    printf("\n Degree counters \n");
    for (i = 0; i <M; i++) {
        printf("%lu ", degC[i]);

    }
     */
    char* rest;
    char* inputFileName = argv[1];
    while (1) {
        rest = strchr(inputFileName, '/');
        if (rest == NULL)break;
        rest++;
        printf("\n%s; len(rest)= %lu\n", rest, strlen(rest));
        inputFileName = rest;
    }


    char* outFileName = (char*) calloc(strlen(inputFileName) + 10, sizeof (char));

    strcat(outFileName, inputFileName);
    strcat(outFileName, ".bin");
    printf("\nWriting to %s, %d vertices and %lu links \n", outFileName, M, szLinkArray);

    FILE* binFile = fopen(outFileName, "wb");

    fwrite(&M, sizeof (int), 1, binFile);
    fwrite(degC, sizeof (unsigned long), M, binFile);
    fwrite(links, sizeof (unsigned int), szLinkArray, binFile);

    fclose(binFile);

    /*
    //now write out matrix 
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    for (i = 0; i < nz; i++)
        fprintf(stdout, "%d %d %20.19g\n", I[i] + 1, J[i] + 1, val[i]);
     */

    //Go Green!

    if (I)free(I);
    if (J)free(J);
    if (val)free(val);
    if (degC)free(degC);
    if (links)free(links);
    if (outFileName) free(outFileName);
    return 0;
}
