#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mesh.h"
#include <unistd.h>
#include "pdb_coordinates.h"

//#include <stddef.h>
//atom_t *readCoordinates(char *fname, int *_nAtom);

atom_t *readCoordinates(char *fname, int *_nAtom) {
    atom_t *head = NULL;
    atom_t *old = NULL;
    head = malloc(sizeof(atom_t));
    atom_t *root = head;
    FILE * fp;
    size_t len = 0;
    ssize_t read;
    char * line = NULL;

    fp = fopen(fname, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    char *p;
    char *end;

    int nAtom = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        /*printf("Retrieved line of length %zu :\n", read);
        printf("%s", line);*/
        p = line;
        double buf[3];
        int i = 0;
        for (double f = strtod(p, &end); p != end; f = strtod(p, &end)) {
            buf[i] = f;
            /*printf("'%.*s' -> ", (int)(end-p), p);*/
            p = end;
            //printf("%f\n", f);
            i++;
        }

        head->x = buf[0];
        head->y = buf[1];
        head->z = buf[2];
        head->nextAtomList = malloc(sizeof(atom_t));
        head = head->nextAtomList;
        head->nextAtomList = NULL;
        nAtom++;
    }
    fclose(fp);
    if (line)
        free(line);
    atom_t *atomArray = malloc (nAtom * sizeof(atom_t));
    head = root;
    int i = 0;
    while (head != NULL) {
        if (i < nAtom) {
            atomArray[i].x = head->x;
            atomArray[i].y = head->y;
            atomArray[i].z = head->z;
        }
        old = head;
        head = head->nextAtomList;
        free(old);
        i++;
    }
    free(head);

    *_nAtom = nAtom;
    return atomArray;
}

void freeBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n) {
    free(x);
    free(y);
    free(z);
    free(chainID);
    for (int i = 0; i < n ; i++) {
        free(resID[i]);
        free(resName[i]);
        free(name[i]);
    }
    free(resID);
    free(resName);
    free(name);
}

// We dont add iCode here, see wheter it is added to resSeq ou resName

int readPdbFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {
    pdbCoordinateContainer_t *pdbCoordinateContainer = pdbFileToContainer(fname);

    int n = pdbCoordinateContainer->atomCount;

    *x = malloc(n * sizeof(double));
    *y = malloc(n * sizeof(double));
    *z = malloc(n * sizeof(double));
    *chainID = malloc(n * sizeof(char));
    *resID = malloc(n * sizeof(char*));
    *resName = malloc(n * sizeof(char*));
    *name = malloc(n * sizeof(char*));
    for (int i = 0 ; i < n ; i++) {

        (*chainID)[i] = pdbCoordinateContainer->atomRecordArray[i].chainID;
        (*x)[i] = pdbCoordinateContainer->atomRecordArray[i].x;
        (*y)[i] = pdbCoordinateContainer->atomRecordArray[i].y;
        (*z)[i] = pdbCoordinateContainer->atomRecordArray[i].z;
        (*resID)[i] = malloc((strlen(pdbCoordinateContainer->atomRecordArray[i].resSeq) + 1) * sizeof(char));
        strcpy((*resID)[i], pdbCoordinateContainer->atomRecordArray[i].resSeq);
        (*resName)[i] = malloc((strlen(pdbCoordinateContainer->atomRecordArray[i].resName) + 1) * sizeof(char));
        strcpy((*resName)[i], pdbCoordinateContainer->atomRecordArray[i].resName);
        (*name)[i] = malloc((strlen(pdbCoordinateContainer->atomRecordArray[i].name) + 1) * sizeof(char));
        strcpy((*name)[i], pdbCoordinateContainer->atomRecordArray[i].name);
    }

    destroyPdbCoordinateContainer(pdbCoordinateContainer);

    return n;
}



int readFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {
    double xBuffer[20000];
    double yBuffer[20000];
    double zBuffer[20000];
    char resIDBuffer[20000][6];
    char chainIDBuffer[20000];
    char resNameBuffer[20000][4];
    char nameBuffer[20000][5];

    FILE *fp;
    fp = fopen(fname, "r");
    char string[100];
    int n = 0;
    char * pch;

    while(fgets(string, 100, fp)) {
        pch = strtok (string,",\n");
        int nF = 0;
        while (pch != NULL) {
            if (nF == 0) chainIDBuffer[n] = pch[0];
            if (nF == 1)
                strcpy(resIDBuffer[n], pch);
            if (nF == 2) xBuffer[n] = atof(pch);
            if (nF == 3) yBuffer[n] = atof(pch);
            if (nF == 4) zBuffer[n] = atof(pch);
            if (nF == 5)
                strcpy(resNameBuffer[n], pch);
            if (nF == 6)
                strcpy(nameBuffer[n], pch);

            pch = strtok (NULL, ",\n");
            nF++;
        }
        n++;
    }

    *x = malloc(n * sizeof(double));
    *y = malloc(n * sizeof(double));
    *z = malloc(n * sizeof(double));
    *chainID = malloc(n * sizeof(char));
    *resID = malloc(n * sizeof(char*));
    *resName = malloc(n * sizeof(char*));
    *name = malloc(n * sizeof(char*));
    for (int i = 0 ; i < n ; i++) {
        (*x)[i] = 2.45;
        (*chainID)[i] = chainIDBuffer[i];
        (*x)[i] = xBuffer[i];
        (*y)[i] = yBuffer[i];
        (*z)[i] = zBuffer[i];
        (*resID)[i] = malloc((strlen(resIDBuffer[i]) + 1) * sizeof(char));
        strcpy((*resID)[i], resIDBuffer[i]);
        (*resName)[i] = malloc((strlen(resNameBuffer[i]) + 1) * sizeof(char));
        strcpy((*resName)[i], resNameBuffer[i]);
        (*name)[i] = malloc((strlen(nameBuffer[i]) + 1) * sizeof(char));
        strcpy((*name)[i], nameBuffer[i]);
    }
    fclose(fp);
    return n;
}


void runSingle( char *fname, float dist, int (*readerFunc)(char*, double**, double**, double**, char**, char***, char***, char***) ) {
    /*  ONE SET OF COORDINATES  */
    double *x;
    double *y;
    double *z;
    char *chainID;
    char **resSeq;
    char **resName;
    char **atomName;
    atom_t *atomList = NULL;
    int nAtom = 0;

    printf("Reading coordinates from %s\n", fname);
    nAtom = (*readerFunc)(fname, &x, &y, &z, &chainID, &resSeq, &resName, &atomName);
    atomList = readFromArrays(nAtom, x, y, z, chainID, resSeq, resName, atomName);
    char *ccmap = residueContactMap(atomList, nAtom, dist);

    atomList = destroyAtomList(atomList, nAtom);

    printf("JSON Single ccmap\n%s\n", ccmap);
    freeBuffers(x, y, z, chainID, resSeq, resName, atomName, nAtom);
    free(ccmap);
}

void runDual( char *iFname, char *jFname, float dist, int (*readerFunc)(char*, double**, double**, double**, char**, char***, char***, char***) ) {
    /*  ONE SET OF COORDINATES  */
    double *x;
    double *y;
    double *z;
    char *chainID;
    char **resSeq;
    char **resName;
    char **atomName;
    atom_t *atomList = NULL;
    int nAtom = 0;
    printf("Reading coordinates from %s\n", iFname);
    nAtom = (*readerFunc)(iFname, &x, &y, &z, &chainID, &resSeq, &resName, &atomName);
    atomList = readFromArrays(nAtom, x, y, z, chainID, resSeq, resName, atomName);

    /*  OPTIONAL SECOND SET OF COORDINATES  */
    double *x_other;
    double *y_other;
    double *z_other;
    char *chainID_other;
    char **resSeq_other;
    char **resName_other;
    char **atomName_other;
    atom_t *atomList_other = NULL;
    int nAtom_other = 0;
    nAtom_other = (*readerFunc)(jFname, &x_other, &y_other, &z_other, &chainID_other, &resSeq_other, &resName_other, &atomName_other);
    atomList_other = readFromArrays(nAtom_other, x_other, y_other, z_other, chainID_other, resSeq_other, resName_other, atomName_other);

    char *ccmap = residueContactMap_DUAL(atomList, nAtom, atomList_other, nAtom_other, dist);

// CLEAR
    atomList = destroyAtomList(atomList, nAtom);
    atomList_other = destroyAtomList(atomList_other, nAtom_other);
    freeBuffers(x, y, z, chainID, resSeq, resName, atomName, nAtom);
    freeBuffers(x_other, y_other, z_other, chainID_other, resSeq_other, resName_other, atomName_other, nAtom_other);

    printf("JSON Dual ccmap\n%s\n", ccmap);
    free(ccmap);
}


int main (int argc, char *argv[]) {

    #ifdef DEBUG
    fprintf(stderr, "*** Debug Mode***\n");
    #endif

    int c;
    char *iFile = NULL;
    char *jFile = NULL;
    char *pdbFile = NULL;

    extern char *optarg;
    extern int optind, optopt, opterr;
    int errflg = 0;
    char *optDist = NULL;
//int readFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {


    int (*readerFunc)(char*, double**, double**, double**, char**, char***, char***, char***) = NULL;
    while ((c = getopt(argc, argv, "a:b:d:p:q:")) != -1) {
        switch(c) {
            case 'p':
                iFile = optarg;
                readerFunc = &readPdbFile;
                break;
            case 'q':
                jFile = optarg;
                break;
            case 'a':
                iFile = optarg;
                readerFunc = &readFile;
                break;
            case 'b':
                jFile = optarg;
                break;
            case 'd':
                optDist = optarg;
                break;
                case ':':/* -f or -o without operand */
                    fprintf(stderr,
                            "Option -%c requires an operand\n", optopt);
                    errflg++;
                    break;
            case '?':
                    fprintf(stderr,
                            "Unrecognized option: -%c\n", optopt);
            errflg++;
        }
    }

    if ( errflg || optDist == NULL || iFile == NULL) {
        fprintf(stderr, "usage: . . . ");
        exit(2);
    }

    double step = atof(optDist);

    if(jFile == NULL)
        runSingle(iFile, step, readerFunc);
    else
        runDual(iFile, jFile, step, readerFunc);

    exit(1);

// Listing neighbouring cells
    //meshDummy(6, 3, 3);

/*
    printf ("--->%s\n", argv[1]);

    atomList = readCoordinates(argv[1], &nAtom);
    printf("Read a %d atoms list\n", nAtom);
    for (int i = 0 ; i < nAtom ; i++) {
        printf("[%d] %g %g %g\n", i, atomList[i].x, atomList[i].y, atomList[i].z);
    }

    mesh(atomList, nAtom, step);
*/
}

