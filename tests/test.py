from gotoh2.aligner import Aligner
g2 = Aligner()
g2.gap_open_penalty = 2
print g2.align('ACGT', 'ACT')
