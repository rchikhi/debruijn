all: debruijn3 merge_graphs

debruijn3: fasta.c kmers.c debruijn3.cpp kmers.h fasta.h
	g++ -Wall -g fasta.c kmers.c debruijn3.cpp -fopenmp -o $@

merge_graphs: fasta.c kmers.c merge_graphs.cpp kmers.h fasta.h
	g++ -Wall -g fasta.c kmers.c merge_graphs.cpp -fopenmp -o $@

clean:
	rm -f debruijn3 merge_graphs *.o
