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

struct align_matrices {
    int nrows;  // all four matrices have the same dimensions
    int ncols;
    int * R;  // linearized cost matrix
    int * p;  //
    int * q;
    int * bits;
};


struct align_output {
    const char * aligned_seq1;
    const char * aligned_seq2;
    int alignment_score;
};

int min2(int a, int b) {
    if (a <= b) {
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



void initialize(struct align_matrices * mx, struct align_settings set) {
    int i, j;  // row and column counters
    int ix, ix2;  // indices into linearized matrix
    int infty = INT_MAX;  // maximum integer value from limits.h

    // initialize the top row (m = 0)
    for (j = 0; j < mx->ncols; j++) {
        mx->p[j] = infty;  // set P_0,j to +Inf
        mx->R[j] = mx->q[j] = set.is_global ? (set.v+set.u*j) : 0;
    }

    mx->R[0] = 0;  // set R_0,0 to 0

    // initialize left column (n = 0)
    for (i = 0; i < mx->nrows; i++) {
        ix = i*mx->ncols;
        mx->q[ix] = infty;  // set Q_i,0 to +Inf
        if (i > 0) {
            mx->R[ix] = set.is_global ? (set.v+set.u*i) : 0;
            ix2 = (i-1)*mx->ncols;
            mx->p[ix] = set.u + min2(mx->p[ix2], mx->R[ix2]+set.v);
        }

        // use a nested loop to zero out the bits matrix
        for (j = 0; j < mx->ncols; j++) {
            mx->bits[ix + j] = 0;
        }
    }

    // set c(l1+1, l2+1) = 1, i.e., 0000100
    mx->bits[mx->nrows * mx->ncols - 1] = 4;
}


// The D matrix has (m+1) rows and (n+1) columns for two sequences of lengths
// m and n, respectively.  Note I'm switching notation from D to R.
void cost_assignment(int * a, int * b, struct align_matrices * mx, struct align_settings set) {
    int i, j;  // row and column counters
    int here, up, left, diag;  // cache indices for linearized matrix

    // iterate through cost matrix by diagonals
    for (int offset=1; offset<mx->ncols+mx->nrows; offset++) {
        j = offset;
        for (i = 1; i < mx->nrows; i++) {
            here = i*mx->ncols + j;
            up = here - mx->ncols;  // (i-1)*ncols + j
            left = here - 1;  // i*ncols + (j-1)
            diag = up - 1;  // (i-1)*ncols + (j-1)

            if (j < mx->ncols) {
                mx->p[here] = set.u + min2(mx->R[up]+set.v, mx->p[up]);

                // can cost P_{i,j} be achieved by edge V_{i-1,j}?
                if (mx->p[here] == mx->p[up]+set.u) {
                    mx->bits[up] |= 8;  // set bit d(i-1,j) == 1
                }

                // can P_{i,j} be achieved without using edge V_{i-1,j}?
                if (mx->p[here] == mx->R[up]+set.u+set.v) {
                    mx->bits[up] |= 16;  // set bit e(i-1,j) == 1
                }

                // find minimum cost of path ending at i,j and using edge H_{i,j}
                mx->q[here] = set.u + min2(mx->R[left]+set.v, mx->q[left]);

                // can cost Q_{i,j} be achieved using edge H_{i,j-1}?
                if (mx->q[here] == mx->q[left]+set.u) {
                    mx->bits[left] |= 32;  // set bit f(i,j-1) == 1
                }
                if (mx->q[here] == mx->R[left]+set.u+set.v) {
                    mx->bits[left] |= 64;  // set bit g(i,j-1) == 1
                }

                // find minimum cost of path ending at N_{i,j}
                // note we subtract (d) from (R) because we are minimizing
                mx->R[here] = min3(mx->R[diag] - set.d[a[i-1]*set.l+b[j-1]],
                               mx->p[here],
                               mx->q[here]);

                // can R_{i,j} be achieved using edges V_{i,j}, H_{i,j} or D_{i,j}?
                if (mx->R[here] == mx->p[here]) {
                    mx->bits[here] |= 1;  // set bit a(m,n) to 1
                }
                if (mx->R[here] == mx->q[here]) {
                    mx->bits[here] |= 2;  // set bit b(m,n) to 1
                }
                if (mx->R[here] == mx->R[diag] - set.d[a[i-1]*set.l+b[j-1]]){
                    mx->bits[here] |= 4;  // set bit c(m,n) to 1
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
void edge_assignment(struct align_matrices * mx) {
    // bits: linearized matrix storing traceback bits
    //    e.g.,   0010010 = 18
    //            gfedcba
    int i, j;  // row, column indices
    int here, down, right, diag;  // indices into linearized arrays
    int nrows = mx->nrows,
        ncols = mx->ncols;  // unpack integers from struct
    int * bits = mx->bits;  // unpack arrays
    int cond1, cond2, cond3, cond4, cond5;

    for (i = nrows-1; i >= 0; i--) {
        for (j = ncols-1; j >= 0; j--) {
            here = i*ncols + j;
            down = here + ncols;  // (i+1)*ncols + j
            right = here + 1;  // i*ncols + (j+1)
            diag = down + 1;  // (i+1)*ncols + (j+1)

            // reset conditions
            cond1 = cond2 = cond3 = cond4 = cond5 = 0;

            // a[i+1,j] is 0 and not last row
            if (i < nrows-1) {
                if (!(bits[down] & 1) ) cond1 = 1;
            }

            // e[i,j] is 0
            if ( !(bits[here] & (1<<4)) ) {
                cond2 = 1;
            }

            // b[i,j+1] is 0 and not last column
            if (j < ncols-1) {
                if ( !(bits[right] & (1<<1)) ) cond3 =  1;
            }

            // g[i,j] is 0
            if ( !(bits[here] & (1<<6)) ) {
                cond4 = 1;
            }

            // c[i+1,j+1] is 0 and not last row nor last column
            if (i < nrows-1 && j < ncols-1) {
                if ( !(bits[diag] & (1<<2)) ) cond5 = 1;
            }

            // step 8
            if ( (cond1 || cond2) && (cond3 || cond4) && cond5 ) {
                bits[here] &= ~1;  // set a[i,j] to 0
                bits[here] &= ~(1<<1);  // set b[i,j] to 0
                bits[here] &= ~(1<<2);  // set c[i,j] to 0
            }

            // step 9. if a[i+1,j] == b[i,j+1] == c[i+1,j+1] == 0, skip steps 10 and 11
            if ( !cond1 && !cond3 && !cond5 ) {
                // step 10
                // if  a[i+1,j] is 1  AND  d[i,j] is 1
                if (i < nrows-1) {
                    if ( bits[down]&1 && bits[here]&(1<<3) ) {
                        // set d[i+1,j] to 1-e[i,j]
                        if (bits[here]&(1<<4)) {
                            bits[down] &= ~(1<<3); // d to 0
                        } else {
                            bits[down] |= (1<<3); // d to 1
                        }

                        // set e[i,j] to 1-a[i,j]
                        if (bits[here]&1) {
                            bits[here] &= ~(1<<4);  // e to 0
                        } else {
                            bits[here] |= (1<<4);  // e to 1
                        }

                        // set a[i,j] to 1
                        bits[here] |= 1;
                    } else {
                        // otherwise, set d[i+1,j] and e[i,j] to 0
                        bits[down] &= ~(1<<3);
                        bits[here] &= ~(1<<4);
                    }
                }

                // step 11
                // if b[i,j+1] is 1  AND  f[i,j] is 1
                if (j < ncols-1) {
                    if ( bits[right]&(1<<1) && bits[here]&(1<<5) ) {
                        // set f[i,j+1] to 1-g[i,j]
                        if (bits[here]&(1<<6)) {
                            bits[right] &= ~(1<<5);  // f to 0
                        } else {
                            bits[right] |= (1<<5);  // f to 1
                        }
                        // set g[i,j] to 1-b[i,j]
                        if (bits[here]&(1<<1)) {
                            bits[here] &= ~(1<<6);  // g to 0
                        } else {
                            bits[here] |= (1<<6);  // g to 1
                        }
                        // set b[i,j] to 1
                        bits[here] |= (1<<1);
                    } else {
                        // otherwise set f[i,j+1] and g[i,j] to 0
                        bits[right] &= ~(1<<5);
                        bits[here] &= ~(1<<6);
                    }
                }
            }
        }
    }
}


int traceback(struct align_matrices mx, struct align_settings set,
               const char * seq1, const char * seq2,
               char * aligned1, char * aligned2) {
    // return all pairwise alignments given edge assignment
    // FIXME: for now, just return one path
    int nrows = mx.nrows,
        ncols = mx.ncols;
    int i, j, here;  // index for linearized matrix
    int init_i, init_j;
    int score,  min_score = INT_MAX;
    int alen = 0;  // track length of pairwise alignment

    if (set.is_global) {
        // start at lower-right cell if global
        init_i = nrows-1;
        init_j = ncols-1;
        min_score = mx.R[init_i*ncols + init_j];
    } else {
        // search right-most column
        for (i=0; i < nrows; i++) {
            here = i*ncols + (ncols-1);
            score = mx.R[here];
            if (score < min_score) {
                init_i = i;
                init_j = ncols-1;
                min_score = score;
            }
        }
        // search bottom row
        for (j=0; j < ncols; j++) {
            here = (nrows-1)*ncols + j;
            score = mx.R[here];
            if (score < min_score) {
                init_i = nrows-1;
                init_j = j;
                min_score = score;
            }
        }
    }
    i = init_i;
    j = init_j;

    // pad right side of alignment if local
    if (i < nrows-1) {
        // add gaps in seq2
        for (alen=0; alen<nrows-1-i; alen++) {
            aligned1[alen] = seq1[nrows-1-alen-1];
            aligned2[alen] = '-';
        }
    }
    if (j < ncols-1) {
        for (alen=0; alen<ncols-1-j; alen++) {
            aligned1[alen] = '-';
            aligned2[alen] = seq2[ncols-1-alen-1];
        }
    }

    while (i>0 && j>0) {
        here = i*mx.ncols + j;
        if (mx.bits[here]&1) {  // a[i,j]==1
            // an optimal path uses V(i,j)
             //fprintf(stdout, "%d %d %d V\n", alen, i, j);
            aligned1[alen] = seq1[i-1];
            aligned2[alen] = '-';
            i--;
        }
        else if (mx.bits[here]&(1<<1)) {  // b[i,j]==1
            // an optimal path uses H(i,j)
            //fprintf(stdout, "%d %d %d H\n", alen, i, j);
            aligned1[alen] = '-';
            aligned2[alen] = seq2[j-1];
            j--;
        }
        else if (mx.bits[here]&(1<<2)) {
            // an optimal path uses D(i,j)
            //fprintf(stdout, "%d %d %d D\n", alen, i, j);
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

    // TODO: pad left side of alignment if local
    while (i>0) {
        aligned1[alen] = seq1[i-1];
        aligned2[alen] = '-';
        i--;
        alen++;
    }
    while (j>0) {
        aligned1[alen] = '-';
        aligned2[alen] = seq2[j-1];
        j--;
        alen++;
    }

    aligned1[alen] = '\0';  // null terminators
    aligned2[alen] = '\0';

    // reverse strings - FIXME: is there a more elegant way to do this?
    char temp[alen];
    for (i=0; i<alen; i++) temp[alen-i-1] = aligned1[i];
    for (i=0; i<alen; i++) aligned1[i] = temp[i];

    // do the same thing for aligned2
    for (i=0; i<alen; i++) temp[alen-i-1] = aligned2[i];
    for (i=0; i<alen; i++) aligned2[i] = temp[i];

    return -min_score;  // this is the alignment score
}


// main wrapper function
struct align_output align(const char * seq1, const char * seq2, struct align_settings set) {
    int map[256] = {0};
    int l1 = (int) strlen(seq1);  // sequence lengths
    int l2 = (int) strlen(seq2);

    int sA[l1];  // integer-encoded sequences
    int sB[l2];
    int align_score;

    struct align_matrices my_matrices;
    my_matrices.nrows = l1+1;
    my_matrices.ncols = l2+1;

    // dynamic allocation of arrays (strings)
    int array_length = (l1+1)*(l2+1);
    my_matrices.R = malloc(sizeof(int) * array_length);
    if (my_matrices.R==NULL) {
        fprintf(stdout, "malloc failed for R\n");
        exit(1);
    }
    my_matrices.p = malloc(sizeof(int) * array_length);
    if (my_matrices.p==NULL) {
        fprintf(stdout, "malloc failed for p\n");
        exit(1);
    }
    my_matrices.q = malloc(sizeof(int) * array_length);
    if (my_matrices.q==NULL) {
        fprintf(stdout, "malloc failed for q\n");
        exit(1);
    }
    my_matrices.bits = malloc(sizeof(int) * array_length);
    if (my_matrices.bits==NULL) {
        fprintf(stdout, "malloc failed for bits\n");
        exit(1);
    }

    struct align_output o;
    char aligned1[l1+l2];  // allocate enough space for completely non-overlapping sequences
    char aligned2[l1+l2];

    // 1. convert sequences into integer indices into alphabet
    map_ascii_to_alphabet(map, set.alphabet);

    encode_sequence(seq1, map, sA);
    encode_sequence(seq2, map, sB);

    // 2. generate D matrix
    // note we are passing align_matrices struct by reference to functions that will modify its contents
    initialize(&my_matrices, set);
    cost_assignment(sA, sB, &my_matrices, set);

    edge_assignment(&my_matrices);

    // DEBUGGING - print cost matrix to screen
    /*
    for (int i=0; i<l1+1; i++) {
        for (int j=0; j<l2+1; j++) {
            fprintf(stdout, "%d ", my_matrices.bits[i*(l2+1) + j]);
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
    */

    // TODO: 3. traceback
    align_score = traceback(my_matrices, set, seq1, seq2, aligned1, aligned2);

    // TODO: decode aligned integer sequences into alphabet sequences
    o.aligned_seq1 = aligned1;
    o.aligned_seq2 = aligned2;
    //o.alignment_score = -1 * my_matrices.R[(l1+1)*(l2+1)-1];  // FIXME: works for global only!
    o.alignment_score = align_score;

    free(my_matrices.R);
    free(my_matrices.p);
    free(my_matrices.q);
    free(my_matrices.bits);

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

    //fprintf (stdout, "my_settings.alphabet = %s\n", my_settings.alphabet);

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

