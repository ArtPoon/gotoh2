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

void encode_sequence (const char * seq, int * map, int * encoded) {
    int seqlen = (int) strlen(seq);
    for (int i=0; i < seqlen; i++) {
        encoded[i] = map[(int)seq[i]];
    }
}


void populate_D_matrix(int * d, int nrows, int ncols, struct align_settings s) {
    int max_dim = nrows;
    if (nrows < ncols) {
        max_dim = ncols;
    }

    // D(0,0) = 0
    d[0] = 0;

    // linearize cache matrices so that (m,n) = (m * ncols + n)
    int p[nrows*ncols];
    int q[nrows*ncols];

    // initialize left column
    for (int m = 0; m < nrows; m++) {
        d[m*ncols] = p[m*ncols] = s.is_global ? (s.v+s.u*m) : 0;
    }

    // initialize the top row
    for (int n = 0; n < ncols; n++) {
        d[n] = q[n] = s.is_global ? (s.v+s.u*n) : 0;
    }

    // iterate through D matrix by diagonals
    // http://stackoverflow.com/questions/1779199/traverse-matrix-in-diagonal-strips
    for (int z, slice=0; slice<(2*max_dim-1); slice++) {
        z = (slice < max_dim) ? 0 : (slice-max_dim+1);  // which side are we on?
        for (int j=z; j<=slice-z; )
    }

    for (int m = 1; m < nrows; m++) {
        for (int n = 1; n < ncols; n++) {
            d[m*ncols+n] = 0;
        }
    }

}


// main wrapper function
struct align_output align(const char * seq1, const char * seq2, struct align_settings m) {
    int map[256] = {0};
    int l1 = (int) strlen(seq1);
    int l2 = (int) strlen(seq2);
    int sA[l1];
    int sB[l2];

    int d[(l1+1)*(l2+1)];
    struct align_output o;

    // 1. convert sequences into integer indices into alphabet
    map_ascii_to_alphabet(map, m.alphabet);

    encode_sequence(seq1, map, sA);
    encode_sequence(seq2, map, sB);

    // 2. generate D matrix
    populate_D_matrix(d, l1+1, l2+1, m);

    for (int i=0; i<l1+1; i++) {
        for (int j=0; j<l2+1; j++) {
            fprintf(stdout, "%d ", d[i*(l2+1)+j]);
        }
        fprintf(stdout, "\n");
    }

    // TODO: decode aligned integer sequences into alphabet sequences
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

