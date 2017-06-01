# gotoh2
A lightweight Python/C module for pairwise alignment of genetic sequences.

## Why did you bother writing yet another pairwise alignment program?
1. This is, in part, a hobby-coding project that I started to keep my C coding skills up (since I do everything in Python and R these days).  
2. There was a major issue where I used to work: the lab relied heavily on an old in-house C program for pairwise alignment that was integrated into most of the bioinformatic pipelines.  The code was difficult to read, poorly documented, and occasionally yielded inconsistent results (enough to become a problem).  
3. Over the years I've tried other pairwise alignment modules in Python and none<sup>1</sup> were as fast, customizable and easy to use as calling an alignment function from the [HyPhy shared library](https://github.com/veg/hyphy-python).  However, calling HyPhy -- an extensive software package for phylogenetic sequence analysis -- is kind of overkill.  
4. Even though there are plenty of powerful [multiple sequence alignment](https://en.wikipedia.org/wiki/Multiple_sequence_alignment) programs available, pairwise alignment can still be a useful tool -- particularly when processing a very large number of fairly long sequences.

<sup>1</sup> Other alternatives have come up since I started working on this in my spare time, so this is probably no longer true.

## Objectives
1. *It should be maintainable.*  I've tried to write accessible C code.
2. *It should be correct.*  I've been implementing a bunch of unit tests to validate the implementation of Altschul and Eriksson's modification of the Gotoh algorithm.
3. *It should be customizable.*  I've incorporated a basic interface at the Python level for importing different residue scoring matrices and for modifying gap penalties.
4. *It should be fast.*  The C should help here, but I still need to run some benchmark tests for comparison.  Ideally faster than calling out to another program.

## Usage example
Aligning two sequences under default settings:
```python
>>> from gotoh2.aligner import Aligner
>>> g2 = Aligner()
>>> g2.align('ACGT', 'ACT')
('ACGT', 'AC-T', 13)
```

## Benchmarking

Here are some times that were required for different programs to do pairwise alignment of the same 10 HIV-1 integrase sequences against the HXB2 reference sequence (867 nucleotides):

| Program | Computing time |
|---------|----------------|
| gotoh2  | 0.3706 seconds |
| Bio.pairwise2 | 89.5627 seconds |
| [alignment](https://github.com/eseraygun/python-alignment) | 18.6270 seconds |
| [alignment](alevchuk) | 7.7929 seconds |


## Requirements
* [Python](https://www.python.org/downloads/) - not sure about version support yet; this was developed in Python 2.7.
* [NumPy](http://www.numpy.org/) - my next objective is to eliminate this requirement
* [A build environment including gcc](https://en.wikipedia.org/wiki/GNU_Compiler_Collection), at least until I figure out how to distribute binaries
