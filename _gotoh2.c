#include <Python.h>
#include "numpy/arrayobject.h"
#include <string.h>
#include <stdio.h>

struct align_settings {
    int is_global;
    int gap_open_penalty;
    int gap_extend_penalty;
    int alphabet_length;
    const char * alphabet;
    double * score_matrix;
};

struct align_output {
    const char * aligned_seq1;
    const char * aligned_seq2;
    int alignment_score;
};


struct align_output align(const char * seq1, const char * seq2, struct align_settings my_settings) {
    struct align_output my_output;
    my_output.aligned_seq1 = seq1;
    my_output.aligned_seq2 = seq2;
    my_output.alignment_score = 0;
    return my_output;
}

static PyObject * align_wrapper(PyObject * self, PyObject * args) {
    const char * alphabet;
    const char * seq1;
    const char * seq2;
    int gop, gep, is_global;

    struct align_settings my_settings;
    struct align_output my_output;

    PyObject * obj = NULL;  // variables for parsing array from Python NumPy object
    PyObject * ndarray = NULL;

    if (!PyArg_ParseTuple(args, "ssiiisO", &seq1, &seq2, &gop, &gep, &is_global, &alphabet, &obj)) {
        return NULL;
    }

    // transfer arguments to struct
    my_settings.gap_open_penalty = gop;
    my_settings.gap_extend_penalty = gep;
    my_settings.is_global = (is_global > 0);
    my_settings.alphabet_length = (int) strlen(alphabet);
    my_settings.alphabet = alphabet;

    // parse NumPy array
    ndarray = PyArray_FROM_OTF(obj, NPY_DOUBLE, NPY_IN_ARRAY);
    if (ndarray == NULL) {
        return NULL;
    }
    my_settings.score_matrix = (double *) PyArray_DATA(ndarray);

    // call align function
    my_output = align(seq1, seq2, my_settings);

    PyObject * retval = Py_BuildValue("ssi", my_output.aligned_seq1, my_output.aligned_seq2, my_output.alignment_score);
    return retval;
}

static PyMethodDef AlignmentMethods [] =
{
    {
        "align", align_wrapper, METH_VARARGS,
        "Pairwise alignment of nucleotide sequences."
    },
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initCgotoh2 (void) {
    (void) Py_InitModule("Cgotoh2", AlignmentMethods);
}

