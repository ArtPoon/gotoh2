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
        print result

if __name__ == '__main__':
    unittest.main()
