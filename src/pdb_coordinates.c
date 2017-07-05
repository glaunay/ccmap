#include "pdb_coordinates.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#define MAX_RECORD 200000
/*
COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
*/
pdbCoordinateContainer_t *pdbFileToContainer(char *fileName) {

    pdbCoordinateContainer_t *pdbCoordinateContainer = malloc(sizeof(pdbCoordinateContainer_t));
    fprintf(stderr, "Reading pdb file %s\n", fileName);

    FILE *fp;
    fp = fopen(fileName, "r");



    char lineBuffer[100];

    int atomCount = 0;
    char atomFlag[] = "ATOM  ";

    static char atomRecordBuffer[MAX_RECORD][81];
    char recordType[7];



    while(fgets(lineBuffer, 100, fp)) {
        memcpy( recordType, &lineBuffer[0], 6 );
        recordType[6] = '\0';
        if (strcmp(recordType, atomFlag) != 0) continue;
        strcpy(atomRecordBuffer[atomCount], lineBuffer);
        atomCount++;
        //printf("%s", lineBuffer);
    }
    fclose(fp);

    fprintf(stderr, "%d atoms have been read\n", atomCount);

    pdbCoordinateContainer->atomRecordArray = malloc(atomCount * sizeof(atomRecord_t));
    pdbCoordinateContainer->atomCount = atomCount;
    atomRecord_t *newAtom = NULL;
    for (int i = 0; i < atomCount ; i++) {
        newAtom = &(pdbCoordinateContainer->atomRecordArray[i]);

        createAtomRecord(atomRecordBuffer[i], newAtom);
        stringifyAtomRecord( newAtom, lineBuffer );


     /*   printf("%s", atomRecordBuffer[i]);
        printf("-->\n%s\n\n", lineBuffer);*/
    }


    return pdbCoordinateContainer;
    /*
*/
}

pdbCoordinateContainer_t *destroyPdbCoordinateContainer(pdbCoordinateContainer_t *pdbCoordinateContainer) {
    free(pdbCoordinateContainer->atomRecordArray);
    free(pdbCoordinateContainer);

    return pdbCoordinateContainer;
}

void createAtomRecord(char *recordString, atomRecord_t *newAtom) {
    char buf[10];

    memcpy( newAtom->recordName, &recordString[0], 6 );
    newAtom->recordName[6] = '\0';

    memcpy( buf, &recordString[6], 5 );
    buf[5] = '\0';
    newAtom->serial = atoi(buf);

    memcpy( newAtom->name, &recordString[12], 4 );
    newAtom->name[4] = '\0';

    newAtom->altLoc = recordString[16];

    memcpy( newAtom->resName, &recordString[17], 3 );
    newAtom->resName[3] = '\0';

    newAtom->chainID = recordString[21];

    memcpy( newAtom->resSeq, &recordString[22], 4 );
    newAtom->resSeq[4] = '\0';

    newAtom->iCode = recordString[26];

    memcpy( buf, &recordString[30], 8 );
    buf[8] = '\0';
    newAtom->x = atof(buf);

    memcpy( buf, &recordString[38], 8 );
    buf[8] = '\0';
    newAtom->y = atof(buf);

    memcpy( buf, &recordString[46], 8 );
    buf[8] = '\0';
    newAtom->z = atof(buf);

    memcpy( buf, &recordString[54], 6 );
    buf[6] = '\0';
    newAtom->occupancy = atof(buf);

    memcpy( buf, &recordString[60], 6 );
    buf[6] = '\0';
    newAtom->tempFactor = atof(buf);

    memcpy( newAtom->element, &recordString[76], 2 );
    newAtom->element[3] = '\0';

    memcpy( newAtom->charge, &recordString[78], 2 );
    newAtom->charge[3] = '\0';

    //return newAtom;
}


void stringifyAtomRecord(atomRecord_t *atomRecord, char *atomRecordString) {
    sprintf(atomRecordString, "%6s%5d %4s%C%3s %c%4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %-2s%-2s",\
            atomRecord->recordName, atomRecord->serial, atomRecord->name, atomRecord->altLoc,\
            atomRecord->resName, atomRecord->chainID, atomRecord->resSeq, atomRecord->iCode,\
            atomRecord->x, atomRecord->y, atomRecord->z, atomRecord->occupancy,\
            atomRecord->tempFactor, atomRecord->element, atomRecord->charge);

}
