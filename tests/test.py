import unittest
from gotoh2.aligner import Aligner

class TestAligner(unittest.TestCase):
    def setUp(self):
        self.g2 = Aligner()
        self.g2.gap_open_penalty = 5

class TestAlignerSimpleGlobal(TestAligner):
    def runTest(self):
        aligned_ref, aligned_query, aligned_score = self.g2.align('ACGT', 'ACT')
        expected = 'ACGT'
        self.assertEqual(expected, aligned_ref)
        expected = 'AC-T'
        self.assertEqual(expected, aligned_query)
        expected = 5+5+(-5-1)+5
        self.assertEqual(expected, aligned_score)

class TestAlignerLongerGlobal(TestAligner):
    def runTest(self):
        aref, aquery, ascore = self.g2.align('ACGTACGTACGTACGT', 'ACGTACGTACTACGT')
        expected = 'ACGTACGTAC-TACGT'
        self.assertEqual(expected, aquery)
        # TODO: run progressively longer sequences

class TestAlignerSimpleLocal(TestAligner):
    def runTest(self):
        self.g2.is_global = False
        aligned_ref, aligned_query, aligned_score = self.g2.align('TACGTA', 'ACGT')
        expected = 'TACGTA'
        self.assertEqual(expected, aligned_ref)
        expected = '-ACGT-'
        self.assertEqual(expected, aligned_query)
        expected = 20
        self.assertEqual(expected, aligned_score)

class TestIssue5(TestAligner):
    def runTest(self):
        # this runs ok
        result = self.g2.align('ACGTT', 'ACGT')
        # this reproducibly crashes!
        result = self.g2.align('ACGT', 'ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTT')

class TestFlouri(TestAligner):
    """
    Evaluate test cases described in Flouri et al. bioRxiv 031500
    """
    def test_NWalign_example(self):
        self.g2.is_global = True
        self.g2.gap_open_penalty = 10
        self.g2.gap_extend_penalty = 1
        self.g2.set_model('NWALIGN')

        a1, a2, result = self.g2.align('GGTGTGA', 'TCGCGT')
        expected = -3
        self.assertEqual(expected, result)

    def test_Biopp_example(self):
        self.g2.is_global = True
        self.g2.gap_open_penalty = 4
        self.g2.gap_extend_penalty = 1
        self.g2.set_model('Biopp')

        a1, a2, score = self.g2.align('AAAGGG', 'TTAAAAGGGGTT')
        print '\n'+a1
        print a2
        print score

if __name__ == '__main__':
    unittest.main()
