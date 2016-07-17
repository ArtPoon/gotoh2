#include <Python.h>
#include "numpy/arrayobject.h"
#include <string.h>
#include <stdio.h>

struct align_settings {
    int is_global;
    int v;  // gap opening penalty
    int u;  // gap extension penalty
    int l;  // alphabet length
    const char * alphabet;
    int * d;  // weighting function
};

struct align_output {
    const char * aligned_seq1;
    const char * aligned_seq2;
    int alignment_score;
};


void map_ascii_to_alphabet (int * map, const char * alphabet) {
    for (int i = 0; i < 256; i++) {
        map[i] = -1;  // if left at -1, indicates non-alphabet character in sequence
    }
    for (int i = 0; i < (int) strlen(alphabet); i++) {
        map[(int)alphabet[i]] = i;  // e.g., "A" in sequence indexes into map[65] = 0
    }
}

int * encode_sequence (const char * seq, int * map) {
    int seqlen = (int) strlen(seq);
    int * encoded[seqlen];
    for (int i=0; i < seqlen; i++) {
        encoded[i] = map[(int)seq[i]];
    }
}



// main wrapper function
struct align_output align(const char * seq1, const char * seq2, struct align_settings m) {
    int map[256] = {0};
    struct align_output o;

    // 1. convert sequences into integer indices into alphabet
    map_ascii_to_alphabet(map, m.alphabet);

    int * sA = encode_sequence(seq1, map);
    int * sB = encode_sequence(seq2, map);

    // 2. generate D matrix


    // X. decode aligned integer sequences into alphabet sequences
    o.aligned_seq1 = seq1;
    o.aligned_seq2 = seq2;
    o.alignment_score = 0;
    return o;
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
    my_settings.v = gop;
    my_settings.u = gep;
    my_settings.is_global = (is_global > 0);
    my_settings.l = (int) strlen(alphabet);
    my_settings.alphabet = alphabet;

    fprintf (stdout, "my_settings.alphabet = %s\n", my_settings.alphabet);

    // parse NumPy array
    ndarray = PyArray_FROM_OTF(obj, NPY_INT, NPY_IN_ARRAY);
    if (ndarray == NULL) {
        return NULL;
    }
    my_settings.d = (int *) PyArray_DATA(ndarray);

    /*
     // display contents of matrix
    for (int i=0; i<my_settings.alphabet_length; i++) {
        for (int j=0; j<my_settings.alphabet_length; j++) {
            fprintf (stdout, "%1.1f ", my_settings.score_matrix[i*my_settings.alphabet_length + j]);
        }
        fprintf(stdout, "\n");
    }
    */

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
    import_array();  // required to avoid a segmentation fault
}

