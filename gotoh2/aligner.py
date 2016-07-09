import Cgotoh2
import numpy as np
from glob import glob
import os
import re

class Aligner():
    def __init__(self, model_path='models/'):
        # default (not very useful) settings
        self.gap_open_penalty = 1
        self.gap_extend_penalty = 1
        self.is_global = True

        # read models from files
        self.models = {}
        files = glob(model_path+'*.csv')
        for f in files:
            _, filename = os.path.split(f)
            model_name = filename.split('.')[0]
            with open(f, 'rU') as handle:
                mx, alpha = self.read_matrix_from_csv(handle)
                self.models.update({model_name: (mx, alpha)})

        # set default model
        self.matrix, self.alphabet = self.models['HYPHY_NUC']

    def __str__(self):
        # TODO: display useful information about alignment settings
        return self.__name__

    def read_matrix_from_csv(self, handle):
        """
        CSV should contain column headers corresponding to the alphabet.  It
        should also be a square matrix (same number of row and column entries).
        :return: (NumPy matrix, str)
        """
        alphabet = ''.join(handle.next().strip('\n').split(','))
        rows = []
        for line in handle:
            rows.append(map(int, line.strip('\n').split(',')))
        return np.matrix(rows), alphabet

    def clean_sequence(self, seq):
        # replace all non-alphabet characters with ambiguous symbol
        return re.sub(pattern='[^%s]' % (self.alphabet,), repl='?', string=seq.upper())

    def align(self, seq1, seq2):
        """
        Main wrapper function that passes data and parameters to C function.
        :param seq1: First sequence to align.
        :param seq2: Second sequence to align.
        :param gop: Gap open penalty.
        :param gep: Gap extension penalty.
        :param is_global: If False, perform local alignment.
        :param model: Key in self.models
        :return:  (aligned seq1, aligned seq2, alignment score)
        """
        assert type(seq1) is str, 'seq1 must be a string'
        assert type(seq2) is str, 'seq2 must be a string'

        results = Cgotoh2.align(
            self.clean_sequence(seq1),
            self.clean_sequence(seq2),
            self.gap_open_penalty,
            self.gap_extend_penalty,
            int(self.is_global),
            self.alphabet,
            self.matrix
        )
        return results
