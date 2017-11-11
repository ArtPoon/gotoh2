import Cgotoh2
import numpy as np
import os
import re
import pkg_resources as pkgres

class Aligner():
    def __init__(self):
        # default settings
        self.gap_open_penalty = 5
        self.gap_extend_penalty = 1
        self.is_global = True

        # read models from files
        self.models = {}
        files = pkgres.resource_listdir('gotoh2', 'models')

        for f in files:
            model_name = f.replace('.csv', '')
            with pkgres.resource_stream('gotoh2', '/'.join(['models', f])) as handle:
                mx, alpha = self.read_matrix_from_csv(handle)
                self.models.update({model_name: (mx, alpha)})

        # set default model
        self.set_model('HYPHY_NUC')

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
        alphabet = ''.join(handle.next().strip('\n').split(','))
        rows = []
        for line in handle:
            rows.append(map(int, line.strip('\n').split(',')))
        return np.array(rows, dtype=np.int32), alphabet

    def set_model(self, model):
        if model in self.models:
            self.matrix, self.alphabet = self.models[model]

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
        assert len(seq1) > 0, 'seq1 cannot be an empty string'
        assert len(seq2) > 0, 'seq2 cannot be an empty string'

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
