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


int min2(a, b) {
    if (a < b) {
        return (a);
    }
    return (b);
}

int min3(a, b, c) {
    if (a < b) {
        return min2(a, c);
    }
    return min2(b, c);
}



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


void populate_D_matrix(int * D, int nrows, int ncols, struct align_settings set,
                       int * a, int * b, int * bits) {
    int max_dim = nrows;
    int infty = 1000000;
    if (nrows < ncols) {
        max_dim = ncols;
    }

    // D(0,0) = 0
    D[0] = 0;

    // linearize cache matrices so that (m,n) = (m * ncols + n)
    int p[nrows*ncols];  // min within column
    int q[nrows*ncols];  // min within row

    // initialize left column (n = 0)
    for (int m = 0; m < nrows; m++) {
        q[m*ncols] = infty;
        D[m*ncols] = p[m*ncols] = set.is_global ? (set.v+set.u*m) : 0;
        // use a nested loop to zero out the bits matrix
        for (int n = 0; n < ncols; n++) {
            bits[m*ncols + n] = 0;
        }
    }
    bits[nrows*ncols-1] = 4;  // set c(l1+1, l2+1) = 1, i.e., 0000100

    // initialize the top row (m = 0)
    for (int n = 0; n < ncols; n++) {
        p[n] = infty;
        D[n] = q[n] = set.is_global ? (set.v+set.u*n) : 0;
    }

    // iterate through D matrix by diagonals
    for (int here, up, left, diag, n, m, offset=1; offset<ncols+nrows; offset++) {
        n = offset;
        for (m = 1; m < nrows; m++) {
            here = m*ncols + n;
            up = here - ncols;  // (m-1)*ncols + n
            left = here - 1;  // m*ncols + (n-1)
            diag = up - 1;  // (m-1)*ncols + (n-1)

            if (n < ncols) {
                p[here] = min2(D[up] + set.u+set.v, p[up] + set.u);
                if (p[here] == p[up] + set.u) {
                    bits[up] |= 8;  // set bit d(m-1,n) == 1
                } else {
                    bits[up] |= 16;  // set bit e(m-1,n) == 1
                }

                q[here] = min2(D[left] + set.u+set.v, q[left] + set.u);
                if (q[here] == q[left] + set.u) {
                    bits[left] |= 32;  // set bit f(m,n-1) == 1
                } else {
                    bits[left] |= 64;  // set bit g(m,n-1) == 1
                }

                D[here] = min3(D[diag] - set.d[a[m-1]*set.l+b[n-1]],
                               p[here],
                               q[here]);
                if (D[here] == p[here]) {
                    bits[here] |= 1;  // set bit a(m,n) to 1
                } else if (D[here] == q[here]) {
                    bits[here] |= 2;  // set bit b(m,n) to 1
                } else {
                    bits[here] |= 4;  // set bit c(m,n) to 1
                }
            }

            n -= 1;
            if (n == 0) {
                break;
            }
        }
    }
}


void traceback(int * D, const char * seq1, const char * seq2, char * aligned1, char * aligned2) {
    // try to implement Altschul-Erickson traceback


}

// main wrapper function
struct align_output align(const char * seq1, const char * seq2, struct align_settings m) {
    int map[256] = {0};
    int l1 = (int) strlen(seq1);  // sequence lengths
    int l2 = (int) strlen(seq2);
    int sA[l1];  // integer-encoded sequences
    int sB[l2];

    int D[(l1+1)*(l2+1)];
    int bits[(l1+1)*(l2+1)];  // stores Altschul and Erickson's seven traceback bits as integer
    struct align_output o;
    char aligned1[l1+l2];  // allocate enough space for completely non-overlapping sequences
    char aligned2[l1+l2];

    // 1. convert sequences into integer indices into alphabet
    map_ascii_to_alphabet(map, m.alphabet);

    encode_sequence(seq1, map, sA);
    encode_sequence(seq2, map, sB);

    // 2. generate D matrix
    populate_D_matrix(D, l1+1, l2+1, m, sA, sB, bits);

    // DEBUGGING - print D matrix to screen
    for (int i=0; i<l1+1; i++) {
        for (int j=0; j<l2+1; j++) {
            fprintf(stdout, "%d ", D[i*(l2+1)+j]);
        }
        fprintf(stdout, "\n");
    }

    // TODO: 3. traceback
    traceback(D, seq1, seq2, aligned1, aligned2);

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

