import _gotoh2
import pkg_resources as pkgres
from gotoh2_utils import *
import math
import argparse
import sys


class Aligner():
    def __init__(self, gop=10, gep=1, is_global=False, model='HYPHY_NUC'):
        """
        :param gop: Gap open penalty
        :param gep: Gap extension penalty
        :param is_global: if False, perform local alignment (no terminal gap penalties)
        :param model: a named substitution model - must be present in models/ directory as CSV file
        """

        # default settings
        self.gap_open_penalty = gop
        self.gap_extend_penalty = gep
        self.is_global = is_global

        # read models from files
        self.models = {}
        files = pkgres.resource_listdir('gotoh2', 'models')

        for f in files:
            model_name = f.replace('.csv', '')
            with pkgres.resource_stream('gotoh2', '/'.join(['models', f])) as handle:
                try:
                    mx, alpha = self.read_matrix_from_csv(handle)
                except:
                    print('Error importing matrix from file {}'.format(f))
                    raise
                self.models.update({model_name: (mx, alpha)})

        # set default model
        self.set_model(model)

    def __str__(self):
        # TODO: display useful information about alignment settings
        output = str(self.alphabet)
        output += '\n' + str(self.matrix)
        output += 'Gap open penalty: {}\nGap extend penalty: {}\n'.format(self.gap_open_penalty, self.gap_extend_penalty)
        return output

    def read_matrix_from_csv(self, handle):
        """
        CSV should contain column headers corresponding to the alphabet.  It
        should also be a square matrix (same number of row and column entries).
        :return: (NumPy matrix, str)
        """
        header = next(handle)
        if type(header) is bytes:
            header = header.decode('ascii')
        alphabet = ''.join(header.strip('\n').split(','))
        rows = []
        for line in handle:
            if type(line) is bytes:
                line = line.decode('ascii')
            values = map(int, line.strip('\n').split(','))
            rows.extend(list(values))

        return rows, alphabet

    def set_model(self, model):
        if model in self.models:
            self.matrix, self.alphabet = self.models[model]
        else:
            print('ERROR: Unrecognized model name {}'.format(model))

    def clean_sequence(self, seq):
        """
        Replace all non-alphabet characters with ambiguous symbol.
        :param seq:  str, sequence
        :return:  str
        """
        return re.sub(pattern='[^%s]' % (self.alphabet,), repl='?', string=seq.upper())

    def align(self, seq1, seq2):
        """
        Main wrapper function that passes data and parameters to C function.
        :param seq1: First sequence to align.
        :param seq2: Second sequence to align.
        :return:  (aligned seq1, aligned seq2, alignment score)
        """
        assert type(seq1) is str, 'seq1 must be a string'
        assert type(seq2) is str, 'seq2 must be a string'
        assert len(seq1) > 0, 'seq1 cannot be an empty string'
        assert len(seq2) > 0, 'seq2 cannot be an empty string'

        results = _gotoh2.align(
            self.clean_sequence(seq1),
            self.clean_sequence(seq2),
            self.gap_open_penalty,
            self.gap_extend_penalty,
            int(self.is_global),
            self.alphabet,
            self.matrix
        )
        return results


def procrust_align(ref, query, aligner):
    """
    Procrustean alignment of query against reference.  Any insertions in the query
    relative to the reference are excised and returned separately.  Filters
    sequences with low alignment scores (when scaled to sequence length, a good
    score should be around 5.0 for HYPHY_NUC - the minimum is 0).

    :param ref:  str, reference sequence
    :param query:  str, query sequence
    :param aligner:  gotoh2.Aligner object
    :return: str, list, float - aligned and trimmed query, list of insertions,
             and normalized alignment score
    """
    aref, aquery, ascore = aligner.align(ref, query)
    norm_score = ascore / len(aref)

    # if query is longer than reference, do not count terminal gaps as insertions
    left, right = len_terminal_gaps(aref)

    trim_seq = ''
    inserts = []
    for i in range(left, len(aref)-right):
        rn = aref[i]
        qn = aquery[i]
        if rn == '-':
            # insertion relative to reference
            inserts.append((i, qn))
            continue
        trim_seq += qn

    return trim_seq, inserts, norm_score


def map_coordinates(ref, query, aligner):
    """
    Generate a dictionary of query nucleotide coordinates to
    a reference coordinate system.
    :param ref:  str, reference sequence
    :param query:  str, query sequence
    :param aligner:  gotoh2.Aligner object
    :return:  dict
    """
    aref, aquery, _ = aligner.align(ref, query)
    qidx, ridx = 0, 0
    coords = {}
    for i, rn in enumerate(aref):
        qn = aquery[i]
        if rn == '-':
            qidx += 1  # insertion in query
        elif qn == '-':
            ridx += 1  # deletion in query
        else:
            coords.update({ridx: qidx})
            qidx += 1
            ridx += 1
    return coords



def codon_align(ref, query, paligner):
    """
    Codon-aware alignment of query to reference sequence.

    :param ref:  str, reference sequence - must be in reading frame
    :param query:  str, query sequence to align
    :param paligner:  gotoh2.Aligner() object configured for protein
                      sequences, i.e., Aligner(gop=5, model='EmpHIV25')
    :return:
    """
    refp = translate_nuc(ref, 0)
    max_score = -math.inf
    br, bq, bo = '', '', 0  # best ref, query and offset

    # determine reading frame of query
    for offset in range(3):
        p = translate_nuc(query, offset)
        aref, aquery, ascore = paligner.align(refp, p)
        if ascore > max_score:
            max_score = ascore
            br, bq, bo = aref, aquery, offset

    # apply AA alignment to nucleotide sequence
    r = apply_prot_to_nuc(br, ref)
    q = apply_prot_to_nuc(bq, '-'*bo + query)

    # excise overlapping region from query sequence
    left, right = get_boundaries(r)
    trimmed = q[left:right]

    return trimmed, max_score


def update_alignment(ref, src, dest, aligner, callback=None):
    """
    Stream sequences from <src> file, check if they are already
    present in <dest> file, and if not do pairwise alignment to
    <ref> and append to <dest>.
    It is assumed that <dest> is the product of Procrustean
    alignment!

    :param ref:  str, reference sequence
    :param src:  source file stream, read mode, FASTA format
    :param dest:  destination file stream, 'r+' mode
    :param aligner:  gotoh2.Aligner() object
    :param callback:  optional function for progress monitoring, assumed
                      to take str (message) argument
    :return:  int, number of sequences transferred to dest
    """
    # cache headers in current FASTA
    prev = {}
    for h, s in iter_fasta(dest):
        if len(s) != len(ref):
            print("Error in update_alignment: length of sequence {} "
                  "does not match reference.".format(h))
            sys.exit()
        prev.update({h: None})

    # iteration moves file pointer to end

    counter = 0
    for h, query in iter_fasta(src):
        if h in prev:
            continue
        if callback:
            callback("Aligning new sequence {}".format(h))
        aquery, _, _ = procrust_align(ref, query, aligner)
        dest.write('>{}\n{}\n'.format(h, aquery))
        counter += 1

    return counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Pairwise sequence alignment utilities for Python"
    )

    parser.add_argument('infile', type=argparse.FileType('r'),
                        help="input, FASTA-formatted file")
    parser.add_argument('ref', type=argparse.FileType('r'),
                        help="input, plain text file with reference sequence")

    parser.add_argument('--out', '-o', default=sys.stdout, type=argparse.FileType('w'),
                        help="output, destination file for alignment; "
                             "defaults to stdout.")
    
    parser.add_argument('--clean', '-c', action='store_true',
                        help="removes duplicates from input file.")

    parser.add_argument('--append', '-a', required=False,
                        type=argparse.FileType('r+'),
                        help="output, open this file in 'r+' (read and append) "
                             "mode to update with aligned sequences.")

    parser.add_argument('--aa', action='store_true',
                        help="Inputs and reference are amino acid sequences, "
                             "defaults to nucleotides.")
    parser.add_argument('--codon', action='store_true',
                        help="Codon-wise alignment of nucleotide sequences.")

    parser.add_argument('--quiet', '-q', action='store_true',
                        help="Suppress progress monitoring.")
    parser.add_argument('--threshold', '-t', type=float, default=0,
                        help="Normalized alignment score cutoff, defaults to 0"
                             " (no cutoff).")
    parser.add_argument('--insfile', type=argparse.FileType('w'),
                        required=False,
                        help="output, file to record insertions relative to "
                             "reference.")

    return parser.parse_args()


if __name__ == '__main__':
    """
    Command line interface
    """
    args = parse_args()

    # initialize Aligner object
    if args.aa:
        aligner = Aligner(model='BLOSUM62')
    else:
        # default nucleotide settings
        aligner = Aligner()

    # load reference sequence
    ref = read_seq(args.ref)

    callback = lambda x: None if args.quiet else print(x)

    if args.append:
        update_alignment(ref, src=args.infile, dest=args.append,
                         aligner=aligner, callback=callback)
    if args.clean:
        cleaned = {}
        
        for h, s in iter_fasta(args.infile):
            if h in cleaned:
                continue
            else:
                cleaned.update({h: None})
                args.out.write('>{}\n{}\n'.format(h, s))
        
    else:
        for h, s in iter_fasta(args.infile):
            trim_seq, inserts, norm_score = procrust_align(ref, s, aligner)
            if norm_score < args.threshold:
                continue

            args.out.write('>{}\n{}\n'.format(h, trim_seq))
            if args.insfile:
                for pos, nt in inserts:
                    args.insfile.write('{},{},{}\n'.format(h, pos, nt))
