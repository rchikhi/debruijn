all: debruijn2 merge_graphs

debruijn2: fasta.c kmers.c debruijn2.cpp kmers.h fasta.h
	g++ -Wall -g fasta.c kmers.c debruijn2.cpp -fopenmp -o debruijn2 

merge_graphs: fasta.c kmers.c merge_graphs.cpp kmers.h fasta.h
	g++ -Wall -g fasta.c kmers.c merge_graphs.cpp -fopenmp -o merge_graphs

clean:
	rm -f debruijn2 merge_graphs *.o
