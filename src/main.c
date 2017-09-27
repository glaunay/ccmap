#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mesh.h"
#include <unistd.h>
#include "pdb_coordinates.h"
#include <getopt.h>

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

// PASS Transflag and rotate on the fly
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

void runDual( char *iFname, char *jFname, float dist,
              int (*readerFunc)(char*, double**, double**, double**, char**, char***, char***, char***)
              ) {
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


void stringToThreeFloats(char *input, float (*vector)[3]) {
    char *err, *p = input;
    if (input == NULL) return;
    float val;
    int i = 0;
    while (*p) {
        val = strtod(p, &err);
        if (p == err) p++;
        else if ((err == NULL) || (*err == 0)) {
            (*vector)[i] = val;
            i++;
           // printf("Value: %f\n", val);
            break;
        }
        else {
            //printf("errValue: %f\n", val);
            p = err + 1;
            (*vector)[i] = val;
            i++;
        }
    }
}

void parseTransform(char *eulerString, char *translateString, float (*eulers)[3], float (*translate)[3]){
    if (eulerString != NULL){
        stringToThreeFloats(eulerString, eulers);
        //printf ("Euler's angles values %g %g %g\n", (*eulers)[0], (*eulers)[1], (*eulers)[2]);
    }
    if (translateString != NULL){
        stringToThreeFloats(translateString, translate);
        //printf ("Cartesian translation vector %g %g %g\n", (*translate)[0], (*translate)[1], (*translate)[2]);
    }
}




/*
    Main program to develop and test C library to manipulate PDB structure in Python 2.7
    --rec ReceptorPdbFile
    --lig LigandPdbFile
    --fmt format (default is PDB)
    --eul Euler's angles triplet
    --trs Cartesian translation vector
    --dst Treshold distance to compute contact
*/

int main (int argc, char *argv[]) {

    #ifdef DEBUG
    fprintf(stderr, "*** Debug Mode***\n");
    #endif

    int c;
    char *iFile = NULL;
    char *jFile = NULL;
    char *pdbFile = NULL;
    char *outFile = NULL;
    char *euler = NULL;
    char *translate = NULL;
    extern char *optarg;
    extern int optind, optopt, opterr;
    int errflg = 0;
    int transflg = 0;
    char *optDist = NULL;
    float eulerAngle[3]  = { 0.0, 0.0, 0.0 };
    float translation[3] = { 0.0, 0.0, 0.0 };

    pdbCoordinateContainer_t *pdbCoordinateContainerJ = NULL;
    pdbCoordinateContainer_t *pdbCoordinateContainerI = NULL;
//int readFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {


    int (*readerFunc)(char*, double**, double**, double**, char**, char***, char***, char***) = NULL;

    readerFunc = &readPdbFile;
    const char    *short_opt = "ha:b:e:t:f:d:w:";
    struct option   long_opt[] =
    {
        {"help",               no_argument, NULL, 'h'},
        {"fmt",          required_argument, NULL, 'f'},
        {"rec",          required_argument, NULL, 'a'},
        {"lig",          required_argument, NULL, 'b'},
        {"euler",          required_argument, NULL, 'e'},
        {"trans",          required_argument, NULL, 't'},
        {"dist",          required_argument, NULL, 'd'},
        {"dump",          required_argument, NULL, 'w'},
        {NULL,            0,                NULL, 0  }
    };

    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        switch(c) {
            case -1:       /* no more arguments */
            case 0:        /* long options toggles */
            break;

            case 'a':
                iFile = strdup(optarg);
                printf("you entered REC \"%s\"\n", optarg);
                break;
            case 'b':
                jFile = strdup(optarg);
                printf("you entered LIG \"%s\"\n", jFile);
                //fprintf(stderr, "you entered LIG \"%s\"\n", jFile);
                break;
            case 't':
                translate = strdup(optarg);
                fprintf(stderr, "%s\n", translate);
                break;
            case 'e':
                euler = strdup(optarg);
                fprintf(stderr, "%s\n", euler);
                break;
            case 'd':
                optDist = strdup(optarg);
                break;
            case 'w':
                outFile = strdup(optarg);
                fprintf(stderr, "will dump to %s\n", outFile);
                break;
            case 'h':
                printf("Usage: %s [OPTIONS]\n", argv[0]);
                printf("  -f file                   file\n");
                printf("  -h, --help                print this help and exit\n");
                printf("\n");
                return(0);

            case ':':
            case '?':
                fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
                return(-2);

            default:
                fprintf(stderr, "%s: invalid option -- %c\n", argv[0], c);
                fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
                return(-2);
        }
    }

    if (iFile != NULL)
        pdbCoordinateContainerI = pdbFileToContainer(iFile);

    if (jFile != NULL)
        pdbCoordinateContainerJ = pdbFileToContainer(jFile);


     if (translate != NULL || euler != NULL) {
        transflg++;
        parseTransform(euler, translate, &eulerAngle, &translation);
    }

    if (jFile != NULL && transflg) {
        printf("Transformation to Ligand molecule\n");
        transformPdbCoordinateContainer(pdbCoordinateContainerJ, eulerAngle, translation);
        //pdbCoordinateContainerJ = pdbFileToContainer(jFile);
    }

    if(outFile != NULL) {
        pdbContainerToFile(pdbCoordinateContainerJ, outFile);
    }


// No Distance, no ccmap computations
    if (optDist == NULL) {
    }

    fprintf(stderr,"Exiting\n");
    if (pdbCoordinateContainerJ != NULL)
        destroyPdbCoordinateContainer(pdbCoordinateContainerJ);


    exit(0);
}
/*
    if ( errflg || optDist == NULL || iFile == NULL) {
        fprintf(stderr, "usage: . . . ");
        exit(2);
    }

    double step = atof(optDist);

    if(jFile == NULL)
        runSingle(iFile, step, readerFunc);
    else
        runDual(iFile, jFile, step, readerFunc);
*/

