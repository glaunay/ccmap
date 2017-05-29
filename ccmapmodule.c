#include <Python.h>
#include <stdlib.h>
#include "mesh.h"
#include "mesh_io.h"
//#include "python_logging.h"
/*
    Python C-API, pure-C mesh implementation interface
*/

static int unpackChainID(PyObject *pListChainID, char **buffer) {
#ifdef DEBUG
    PySys_WriteStdout("--->Unpack chainID\n");
#endif
    PyObject *pItem;
    Py_ssize_t n;
    int i;


    n = PyList_Size(pListChainID);
   // PySys_WriteStdout("--->%d\n", n);

    *buffer = PyMem_New(char, n);

    /*PyObject* objectsRepresentation = PyObject_Repr(yourObject);
    const char* s = PyString_AsString(objectsRepresentation);
    */
    PyObject* objectsRepresentation;
    const char* s;
    for (i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListChainID, i);
        objectsRepresentation = PyObject_Repr(pItem);
        s = PyString_AsString(objectsRepresentation);
        (*buffer)[i] = s[1];
    }
    return 1;
}

static int unpackString(PyObject *pListOfStrings, char ***buffer) {
#ifdef DEBUG
    PySys_WriteStdout("--->Unpack string\n");
#endif
    PyObject *pItem;
    Py_ssize_t n;
    int i;

    //char **buffer = *_buffer;
    n = PyList_Size(pListOfStrings);
    *buffer = PyMem_New(char*, n);

    /*PyObject* objectsRepresentation = PyObject_Repr(yourObject);
    const char* s = PyString_AsString(objectsRepresentation);
    */
    PyObject* objectsRepresentation;
    const char* s;
    int sLen;
    for (i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListOfStrings, i);
        objectsRepresentation = PyObject_Repr(pItem);
        s = PyString_AsString(objectsRepresentation);
        sLen =  strlen(s); // This corresponds to the actual string surrounded by \' , ie : 'MYTSRING'
        //PySys_WriteStdout("--->%s[%d]\n", s, strlen(s));
        (*buffer)[i] = PyMem_New(char, sLen - 1);
        for (int j = 1 ; j < sLen - 1 ; j++) {
            (*buffer)[i][j - 1] = s[j];
        }
        (*buffer)[i][sLen - 2] = '\0';
        //PySys_WriteStdout("translated to --->\"%s\"[%d]\n", (*buffer)[i], sLen - 1);
       // PySys_WriteStdout("translated to --->%s[%d]\n", (*buffer)[i]);
    }
    return 1;
}


static int unpackCoordinates(PyObject *pListCoor, double **buffer) {
#ifdef DEBUG
    PySys_WriteStdout("*** Unpack coordinates ***\n");
#endif
    PyObject *pItem;
    Py_ssize_t n;

    int i;
    n = PyList_Size(pListCoor);
#ifdef DEBUG
    PySys_WriteStdout("Allocating for %d coordinates\n", n);
#endif
    *buffer = PyMem_New(double, n);

    double u;
    for (i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListCoor, i);
        if(!PyFloat_Check(pItem)) {
            PyErr_SetString(PyExc_TypeError, "coordinate items must be float.");
            PyMem_Free(*buffer);
            return NULL;
        }

        PyObject_AsDouble(pItem, &u);
        PyObject_AsDouble(pItem, &(*buffer)[i]);
    }
    return 1;
}





static PyObject *
ccmap_compute(PyObject *self, PyObject *args)
{
    char dummy[] = "XY\0";

/*
PyArg_ParseTuple() returns true (nonzero) if all arguments have the right type and its components have been stored in
the variables whose addresses are passed. It returns false (zero) if an invalid argument list was passed.
In the latter case it also raises an appropriate exception so the calling function can return NULL immediately (as we saw in the example).
*/

PyObject *pListX, *pListY, *pListZ, *pListChainID, *pListResName;
float userThreshold;
if (!PyArg_ParseTuple(args, "O!O!O!O!O!f", &PyList_Type, &pListX, &PyList_Type, &pListY, &PyList_Type, &pListZ, &PyList_Type, &pListChainID, &PyList_Type, &pListResSeq, &PyList_Type, &pListResName, &PyList_Type, &pListAtomName,&userThreshold)) {
    PyErr_SetString(PyExc_TypeError, "parameters must be Five and a distance value lists.");
    return NULL;
}
    Py_ssize_t n = PyList_Size(pListX);
#ifdef DEBUG
    PySys_WriteStdout("Unpacking %d atoms [contact distance is %f]\n", n, userThreshold);
#endif
    /* Unpacking atom specification vectors */
    double *coorX, *coorY, *coorZ;
    unpackCoordinates(pListX, &coorX);
    unpackCoordinates(pListY, &coorY);
    unpackCoordinates(pListZ, &coorZ);

    char *chainID;
    unpackChainID(pListChainID, &chainID);

    char **resSeq;
    unpackString(pListResSeq, &resSeq);

    char **resName;
    unpackString(pListResName, &resName);

    char **atomName;
    unpackString(pListAtomName, &atomName);
    /* Create data structures and compute */
    atom_t *atomList = readFromArrays(n, coorX, coorY, coorZ, chainID, resSeq, resName, atomName);

    char *ccmap = residueContactMap(atomList, n, userThreshold);


// We may have to copy string in PythonString Object b4 return ?

    return Py_BuildValue("s", ccmap);
}


int PyObject_AsDouble(PyObject *py_obj, double *x)
{
  PyObject *py_float;

  py_float = PyNumber_Float(py_obj);

  if (py_float == NULL) return -1;

  *x = PyFloat_AsDouble(py_float);

  Py_DECREF(py_float);
  return 0;
}

int PyList_IntoDoubleArray(PyObject *py_list, double *x, int size) {
    int i;

    if (py_list == NULL) return 1;

    if (!PyList_Check(py_list)) return 1;

    if (size != PyList_Size(py_list)) return 1;

    for (i=0; i<size; i++) {
        PyObject *py_float = PyList_GetItem(py_list, i);
        if (py_float == NULL || PyObject_AsDouble(py_float, &(x[i])))
        return 1;
    }

    return 0;
}

static PyMethodDef CcmapMethods[] = {
    //...
    {"compute",  ccmap_compute, METH_VARARGS,
     "Compute a contact map."},
    //...
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initccmap(void)
{
    (void) Py_InitModule("ccmap", CcmapMethods);
}


int
main(int argc, char *argv[])
{
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    initccmap();
}
