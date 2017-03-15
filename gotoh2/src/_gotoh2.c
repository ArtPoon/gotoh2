#include <Python.h>
#include "numpy/arrayobject.h"
#include <string.h>
#include <stdio.h>
#include <limits.h>


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

int min2(int a, int b) {
    if (a < b) {
        return (a);
    }
    return (b);
}

int min3(int a, int b, int c) {
    if (a < b) {
        return min2(a, c);
    }
    return min2(b, c);
}



void map_ascii_to_alphabet (int * map, const char * alphabet) {
    // map: integer vector that maps from ASCII into alphabet string index
    // alphabet: string to be indexed; e.g., "ACGT"
    for (int i = 0; i < 256; i++) {
        map[i] = -1;  // if left at -1, indicates non-alphabet character in sequence
    }
    for (int i = 0; i < (int) strlen(alphabet); i++) {
        map[(int)alphabet[i]] = i;  // e.g., "A" in sequence indexes into map[65] = 0
    }
}

void encode_sequence (const char * seq, int * map, int * encoded) {
    // seq: string to convert into integer indices
    // map: return value from map_ascii_to_alphabet()
    // encoded: container for return value
    int seqlen = (int) strlen(seq);
    for (int i=0; i < seqlen; i++) {
        encoded[i] = map[(int)seq[i]];
    }
}



void initialize(int * R, int * p, int * q, int nrows, int ncols, struct align_settings set, int * bits) {
    int i, j;  // row and column counters
    int infty = INT_MAX;  // maximum integer value from limits.h

    // initialize the top row (m = 0)
    for (j = 0; j < ncols; j++) {
        p[j] = infty;  // set P_0,j to +Inf
        R[j] = q[j] = set.is_global ? (set.v+set.u*j) : 0;
    }

    // initialize left column (n = 0)
    for (i = 0; i < nrows; i++) {
        q[i*ncols] = infty;  // set Q_i,0 to +Inf
        R[i*ncols] = set.is_global ? (set.v+set.u*i) : 0;
        // use a nested loop to zero out the bits matrix
        for (j = 0; j < ncols; j++) {
            bits[i*ncols + j] = 0;
        }
    }
    bits[nrows*ncols-1] = 4;  // set c(l1+1, l2+1) = 1, i.e., 0000100

    R[0] = 0;  // set R_0,0 to 0
}


// The D matrix has (m+1) rows and (n+1) columns for two sequences of lengths
// m and n, respectively.  Note I'm switching notation from D to R.
void cost_assignment(int * R, int * p, int * q, int nrows, int ncols, int * a, int * b, struct align_settings set, int * bits) {
    int i, j;  // row and column counters
    int here, up, left, diag;  // cache indices for linearized matrix
    int max_dim = nrows;
    if (nrows < ncols) {
        max_dim = ncols;
    }

    // iterate through cost matrix by diagonals
    for (int offset=1; offset<ncols+nrows; offset++) {
        j = offset;
        for (i = 1; i < nrows; i++) {
            here = i*ncols + j;
            up = here - ncols;  // (i-1)*ncols + j
            left = here - 1;  // i*ncols + (j-1)
            diag = up - 1;  // (i-1)*ncols + (j-1)

            if (j < ncols) {
                p[here] = set.u + min2(R[up]+set.v, p[up]);

                // can cost P_{i,j} be achieved by edge V_{i-1,j}?
                if (p[here] == p[up]+set.u) {
                    bits[up] |= 8;  // set bit d(i-1,j) == 1
                }

                // can P_{i,j} be achieved without using edge V_{i-1,j}?
                if (p[here] == R[up]+set.u+set.v) {
                    bits[up] |= 16;  // set bit e(i-1,j) == 1
                }

                // find minimum cost of path ending at i,j and using edge H_{i,j}
                q[here] = set.u + min2(R[left]+set.v, q[left]);

                // can cost Q_{i,j} be achieved using edge H_{i,j-1}?
                if (q[here] == q[left]+set.u) {
                    bits[left] |= 32;  // set bit f(i,j-1) == 1
                }
                if (q[here] == R[left]+set.u+set.v) {
                    bits[left] |= 64;  // set bit g(i,j-1) == 1
                }

                // find minimum cost of path ending at N_{i,j}
                // note we subtract (d) from (R) because we are minimizing
                R[here] = min3(R[diag] - set.d[a[i-1]*set.l+b[j-1]],
                               p[here],
                               q[here]);

                // can R_{i,j} be achieved using edges V_{i,j}, H_{i,j} or D_{i,j}?
                if (R[here] == p[here]) {
                    bits[here] |= 1;  // set bit a(m,n) to 1
                }
                if (R[here] == q[here]) {
                    bits[here] |= 2;  // set bit b(m,n) to 1
                }
                if (R[here] == R[diag] - set.d[a[i-1]*set.l+b[j-1]]){
                    bits[here] |= 4;  // set bit c(m,n) to 1
                }
            }

            j -= 1;
            if (j == 0) {
                break;
            }
        }
    }
}


// steps 8 through 11 inclusive of Altschul-Erickson algorithm
void edge_assignment(int * bits, int nrows, int ncols) {
    // bits: linearized matrix storing traceback bits
    //    e.g.,   0010010 = 18
    //            gfedcba

    // switching index integers from (m,n) to (i,j) because that's what
    // Altschul and Erickson use and it's easier to follow their pseudocode this way
    for (int here, down, right, diag, i = nrows-1; i >= 0; i--) {
        for (int j = ncols-1; j >= 0; j--) {
            here = i*ncols + j;
            down = here + ncols;  // (i+1)*ncols + j
            right = here + 1;  // i*ncols + (j+1)
            diag = down + 1;  // (i+1)*ncols + (j+1)

            if (
                (!(bits[down]&1) || !(bits[here]&(1<<4))) &&  // a[i+1,j] is 0  OR  e[i,j] is 0
                (!(bits[right]&(1<<1)) || !(bits[here]&(1<<6))) &&  // b[i,j+1] is 0  OR  g[i,j] is 0
                !(bits[diag]&(1<<2))  // c[i+i,j+1] is 0
            ) {
                bits[here] |= 1;  // set a[i,j] to 1
                bits[here] |= 2;  // set b[i,j] to 1
                bits[here] |= 4;  // set c[i,j] to 1
            }

            // if a[i+1,j] == b[i,j+1] == c[i+1,j+1] == 0, skip steps 10 and 11
            if (!(bits[down]&1) && !(bits[right]&(1<<1)) && !(bits[diag]&(1<<2))) {
                // step 10
                if (bits[down]&1 && bits[here]&(1<<3)) {  // if a[i+1,j] is 1 and d[i,j] is 1
                    // set d[i+1,j] to 1-e[i,j]
                    if (bits[here]&(1<<4)) {
                        bits[down] = bits[down] & ~(1<<3);
                    } else {
                        bits[down] |= 8;
                    }

                    // set e[i.j] to 1-a[i,j]
                    if (bits[here]&1) {
                        bits[here] = bits[here] & ~(1<<4);
                    } else {
                        bits[here] |= 16;
                    }

                    // set a[i,j] to 1
                    bits[here] |= 1;
                } else {
                    // otherwise, set d[i+1,j] and e[i,j] to 0
                    bits[down] = bits[down] & ~(1<<3);
                    bits[here] = bits[here] & ~(1<<4);
                }

                // step 11
                if ((bits[right]&2) && (bits[here]&32)) {  // if b[i,j+1] is 1  AND  f[i,j] is 1
                    // set f[i,j+1] to 1-g[i,j]
                    if (bits[here]&(1<<6)) {
                        bits[right] = bits[right] & ~(1<<5);
                    } else {
                        bits[right] |= 32;
                    }
                    // set g[i,j] to 1-b[i,j]
                    if (bits[here]&(1<<1)) {
                        bits[here] = bits[here] & ~(1<<6);
                    } else {
                        bits[here] |= 64;
                    }
                    // set b[i,j] to 1
                    bits[here] |= 2;
                } else {
                    // otherwise set f[i,j+1] and g[i,j] to 0
                    bits[right] = bits[right] & ~(1<<5);
                    bits[here] = bits[here] & ~(1<<6);
                }
            }
        }
    }
}

void traceback(int * bits, int nrows, int ncols,
               const char * seq1, const char * seq2,
               char * aligned1, char * aligned2) {
    // return all pairwise alignments given edge assignment
    // FIXME: for now, just return one path
    fprintf(stdout, "seq1: %s\nseq2: %s\n", seq1, seq2);

    // start from the lower right of our matrix
    int i = nrows-1,  // nrows is len(seq1)+1, adjust for zero-index
        j = ncols-1,
        here;
    int alen = 0;

    while (i>0 && j>0) {
        fprintf(stdout, "aligned1: '%s'\n", aligned1);
        here = i*ncols + j;
        if (bits[here]&1) {  // a[i,j]==1
            // an optimal path uses V(i,j)

            // append gap to aligned1
            aligned1[alen] = '-';
            // append base from seq2 to aligned2
            aligned2[alen] = seq2[j-1];

            fprintf(stdout, "%d %d V\n", i, j);
            i--;
        }
        else if (bits[here]&(1<<1)) {  // b[i,j]==1
            // an optimal path uses H(i,j)
            fprintf(stdout, "%d %d H\n", i, j);
            aligned1[alen] = seq1[i-1];
            aligned2[alen] = '-';
            j--;
        }
        else if (bits[here]&(1<<2)) {
            // an optimal path uses H(i,j)
            fprintf(stdout, "%d %d D\n", i, j);
            aligned1[alen] = seq1[i-1];
            aligned2[alen] = seq2[j-1];
            i--;
            j--;
        }
        else {
            fprintf(stdout, "uh oh, no optimal path?");
            exit(1);
        }
        alen++;
    }
    aligned1[alen] = '\0';
    fprintf(stdout, "aligned1: %s\n", aligned1);

    // reverse strings
    char temp[alen];
    for (i=0; i<alen; i++) {
        fprintf(stdout, "%d %c\n", i, aligned1[i]);
        temp[alen-i-1] = aligned1[i];
    }
}

// main wrapper function
struct align_output align(const char * seq1, const char * seq2, struct align_settings m) {
    int map[256] = {0};
    int l1 = (int) strlen(seq1);  // sequence lengths
    int l2 = (int) strlen(seq2);

    int sA[l1];  // integer-encoded sequences
    int sB[l2];

    int D[(l1+1)*(l2+1)];
    int P[(l1+1)*(l2+1)];
    int Q[(l1+1)*(l2+1)];

    int bits[(l1+1)*(l2+1)];  // stores Altschul and Erickson's seven traceback bits as integer
    struct align_output o;
    char aligned1[l1+l2];  // allocate enough space for completely non-overlapping sequences
    char aligned2[l1+l2];
    int alen;  // track lengths of aligned sequences

    // 1. convert sequences into integer indices into alphabet
    map_ascii_to_alphabet(map, m.alphabet);

    encode_sequence(seq1, map, sA);
    encode_sequence(seq2, map, sB);

    // 2. generate D matrix
    initialize(D, P, Q, l1+1, l2+1, m, bits);
    cost_assignment(D, P, Q, l1+1, l2+1, sA, sB, m, bits);

    // DEBUGGING - print cost matrix to screen
    for (int i=0; i<l1+1; i++) {
        for (int j=0; j<l2+1; j++) {
            fprintf(stdout, "%d ", D[i*(l2+1)+j]);
        }
        fprintf(stdout, "\n");
    }

    edge_assignment(bits, l1+1, l2+1);

    // TODO: 3. traceback
    traceback(bits, l1+1, l2+1, seq1, seq2, aligned1, aligned2);

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

