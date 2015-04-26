Note
===

This software is quite old (2011). There are better techniques nowadays. A better de Bruijn graph construction algorithm is BCALM: https://github.com/Malfoy/bcalm (readable Python code here: https://github.com/rchikhi/python-bcalm).
Another is an implementation included in the book Genome-Scale Algorithm Design (Chapter 13.2, http://www.genome-scale.info/implementations.html). 

Debruijn
===

Software that constructs the de Bruijn graph of a set of reads (FASTA or FASTQ file).
Edges are the (k+1)-mers that appear in the reads, nodes are the k-mers.
No distintion is made between a k-mer (resp. (k+1)-mer) and its reverse-complement.
It returns a graph in the KisSplice format (ad-hoc) or DOT format (can be opened by 
most applications, including Zgrviewer and Gephi). 

The de Bruijn graph is useful for many next-generation sequencing applications,
including de novo genome assembly and variant detection.

For k <= 32, this implementation uses 

            max( 64*N/p, 64*G )

bits of memory, where:

- N is the number of k-mers present in the reads (N = number of reads * (read length - k + 1))  
- G is the number of distinct k-mers in the genome (G = roughly the size of the genome). 
- p is the number of passes (specified by option -p in the software). 

For 3 billion distinct k-mers, assuming sufficiently many passes, it should construct the graph in 24 Gb of memory. 

Note: for 32 <= k <= 64, the constant 64 becomes 128. K-mers larger than 64 nucleotides are not supported.

written by Rayan Chikhi (http://www.irisa.fr/symbiose/rayan_chikhi)
part of the KisSplice package (http://alcovna.genouest.org/kissplice/)
