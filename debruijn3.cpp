/**
 * Copyright INRIA, CNRS, ENS, main contributors: Chikhi, Lacroix, Peterlongo, Sacomoto
 * rayan.chikhi@irisa.fr
 * Vincent.Lacroix@univ-lyon1.fr
 * pierre.peterlongo@inria.fr
 * sacomoto@gmail.com
 *
 * The kisSplice software calls splicing events in a de bruijn bi-connected
 * graph created from one or two sets of HTS reads
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * http://www.cecill.info.

 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.

 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <omp.h>
#include <stdint.h>
#include <limits>

// to print uint64_t
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "fasta.h"
#include "kmers.h"

using namespace std;

#define MAX_THREADS 1

//#define KMERS_OVER_32 // disabled to compare with trinity
#if (defined KMERS_OVER_32 && defined _LP64)
typedef __uint128_t large_integer;
#else
typedef uint64_t large_integer;
#endif

typedef uint64_t large_index_integer;
typedef vector<large_integer> huge_vector;
typedef int abundance_type;
int max_abundance = numeric_limits<abundance_type>::max();
huge_vector* kmers_long[MAX_THREADS];
char ** reads_file_names;
int number_of_read_sets;



// default parameters
string out_prefix = "test";
unsigned int k=28;
unsigned int weight_threshold = 0;
short graph_format = 1;
bool dont_compute_nodes_abundance = false;
bool normalize_edges = true;
int partitions = 2;

// hard-coded stuff
string graph_file_suffix = ".graph";
string nodes_without_abundance_file_suffix = ".nodes_without_abundance";
string nodes_with_abundance_file_suffix = ".nodes";
string raw_edges_file_suffix = ".raw_edges";
string edges_file_suffix = ".edges";

/*
 *
 *
 */

int tables_init()
{
    //if (simple<3)
    //nb_threads=omp_get_max_threads()-1;
    int nb_threads=1;
    //omp_set_num_threads(nb_threads);

    for(int i=0; i<nb_threads; i++)
    {
        kmers_long[i]=new huge_vector;
    }
    return nb_threads;
}

void tables_free()
{   int nb_threads=1;
    for(int i=0; i<nb_threads; i++)
    {
        delete kmers_long[i];
    }
}

large_integer large_kmer_hash_from_tight(unsigned char *kmer_tight, unsigned int len_tightk_plus_bytesize)
{

    large_integer res;
    unsigned int i;
    res=0;
    for (i=0; i<len_tightk_plus_bytesize-TIGHT_KMER_BYTESIZE; i++)
    {
        res+=(((large_integer)(kmer_tight[TIGHT_KMER_BYTESIZE+i]))<<(large_integer)(i*8));
    }
    return res;
}

void large_kmer_unhash_to_tight(large_integer kmer_hashed, unsigned int k, unsigned char* result_kmer_tight)
{
    unsigned int i;
    unsigned int len_tightk = (k / 4) + (k % 4 ? 1 : 0);

    // set the size
    if (TIGHT_KMER_BYTESIZE==1)
        result_kmer_tight[0]=k;
    else if (TIGHT_KMER_BYTESIZE==2)
    {
        result_kmer_tight[0]=(k%0xFF);
        result_kmer_tight[1]=(k/0xFF);
    }

    // decode raw bits
    for (i=0; i<len_tightk; i++)
        result_kmer_tight[TIGHT_KMER_BYTESIZE+i]=((kmer_hashed>>(i*8))&0xFF);
}

bool normalize_kmer(unsigned char * kmer_tight)
{
    if (arbitrary_criterion(kmer_tight))
    {
        revcomp_tight(kmer_tight);
        return true;
    }
    return false;
}

// beware of polymorphism!
large_integer normalize_kmer(large_integer kmer_long, unsigned int k)
{
    unsigned int len_tightk_plus_bytesize = (k / 4) + (k % 4 ? 1 : 0) + TIGHT_KMER_BYTESIZE;
    unsigned char kmer_tight[MAX_TIGHT_SIZE_TIGHT];
    large_kmer_unhash_to_tight(kmer_long,k,kmer_tight);
    normalize_kmer(kmer_tight);
    return large_kmer_hash_from_tight(kmer_tight,len_tightk_plus_bytesize);
}

void insert_kmers_long(large_integer kmer_long, int thread_num)
{
    kmers_long[thread_num]->push_back(kmer_long);
    if (kmers_long[thread_num]->size()%20000000==0)
        printf("inserted %ul (possibly redundant) kmers\n",kmers_long[thread_num]->size());
}

large_index_integer get_kmer_number(large_integer kmer_long)
{
    int thread_num=0;
    large_index_integer index=0;
    index=upper_bound(kmers_long[thread_num]->begin(), kmers_long[thread_num]->end(), kmer_long) - kmers_long[thread_num]->begin() - 1;
    return index;
}

bool is_kmer_present(large_integer kmer_long)
{
    int thread_num=0;
    return binary_search(kmers_long[thread_num]->begin(), kmers_long[thread_num]->end(), kmer_long);
}

// output a single node to a file
FILE *graph_file,*nodes_file,*edges_file;
void print_node(large_index_integer index, char *ascii_kmer, char *revcomp_ascii_kmer)
{
    if (graph_format==0) // DOT format
        fprintf(graph_file,"%" PRIu64 " [label=\"%s / %s\"];\n",index,ascii_kmer,revcomp_ascii_kmer);
    else if (graph_format==1) // kissplice format
        fprintf(nodes_file,"%" PRIu64 "\t%s\n",index,ascii_kmer);

}

// output a single edge to a file
void print_edge(large_index_integer id, large_index_integer id2, string label, int weight)
{
    if (graph_format==0) // DOT format
        fprintf(graph_file,"%"PRIu64" -> %"PRIu64" [label=\"%s\" weight=%d];\n",id,id2,label.c_str(),weight);
    else if (graph_format==1) // kissplice format
        fprintf(edges_file,"%"PRIu64"\t%"PRIu64"\t%s\n",id,id2,label.c_str());

}

void sort_kmers(int thread_num)
{
    printf(" sorting\n");
    //omptl::
    sort(kmers_long[thread_num]->begin(),kmers_long[thread_num]->end());
}

void insert_pairs_of_kmers(large_integer kmer_long_1, large_integer kmer_long_2, unsigned short weight, unsigned int k)
{
    // we cannot recover the information of which (k+1)mer is at the very beginning or very end of a read, so we have to test-or-insert both kmers

    // algorithmic FIXME: there is a bug here: is_kmer_present works only if the list is sorted, which is not the case here. well, we will insert 2x too many kmers at worst. the nodes list will be sorted and redundancy-filtered later anyway, so this bug has no consequence on results.
	if (!is_kmer_present(kmer_long_1))
	    insert_kmers_long(normalize_kmer(kmer_long_1,k),0);
	if (!is_kmer_present(kmer_long_2))
	    insert_kmers_long(normalize_kmer(kmer_long_2,k),0);
}

/*
 *
 *
 *
 */

void write_node(large_index_integer index, large_integer kmer_long, int k)
{
    unsigned char kmer_tight[MAX_TIGHT_SIZE_TIGHT];
    char ascii_kmer[k+1];
    char revcomp_ascii_kmer[k+1];

    large_kmer_unhash_to_tight(kmer_long,k,kmer_tight);
    from_tight_DNA(kmer_tight,ascii_kmer);
    strcpy(revcomp_ascii_kmer,ascii_kmer);
    revcomp(revcomp_ascii_kmer,k);

    print_node(index,ascii_kmer,revcomp_ascii_kmer);
}


unsigned long nb_preliminary_edges = 0;
void write_preliminary_edge(large_integer edge_long, unsigned short weight, int weight_threshold, unsigned int k, FILE *dest)
{
    if (weight>=weight_threshold)
    {
        // decomposing a large_integer (k+1)mer into kmer_1->kmer_2
        large_integer kmer_long_1 = edge_long & ((large_integer(1)<<(k*2)) -1);
        large_integer kmer_long_2 = edge_long>> large_integer(2);
        fwrite(&kmer_long_1,sizeof(kmer_long_1),1,dest);
        fwrite(&kmer_long_2,sizeof(kmer_long_2),1,dest);
        fwrite(&weight,sizeof(weight),1,dest);
        nb_preliminary_edges+=1;
    }

}

unsigned long nb_edges = 0;
void write_edge(large_integer kmer_long_1, large_integer kmer_long_2, unsigned short weight, unsigned int k)
{
    unsigned char kmer_tight[MAX_TIGHT_SIZE_TIGHT];
    unsigned int len_tightk_plus_bytesize = (k / 4) + (k % 4 ? 1 : 0) + TIGHT_KMER_BYTESIZE;

    large_kmer_unhash_to_tight(kmer_long_1,k,kmer_tight);
    bool kmer1_reversed = normalize_kmer(kmer_tight);
    kmer_long_1 = large_kmer_hash_from_tight(kmer_tight,len_tightk_plus_bytesize);

    large_kmer_unhash_to_tight(kmer_long_2,k,kmer_tight);
    bool kmer2_reversed = normalize_kmer(kmer_tight);
    kmer_long_2 = large_kmer_hash_from_tight(kmer_tight,len_tightk_plus_bytesize);

    large_index_integer kmer_long_1_index = get_kmer_number(kmer_long_1);
    large_index_integer kmer_long_2_index = get_kmer_number(kmer_long_2);

	// kmer_long_1 and kmer_long_2 are raw,non-normalized kmers (normalize = pick reverse complement if it's lexically before)
	// the label associated to each normalized kmer is R if the kmer has been reversed, else F
    string label = (string)(kmer1_reversed?"R":"F") + (string)(kmer2_reversed?"R":"F");
    print_edge(kmer_long_1_index,kmer_long_2_index,label,weight);
    nb_edges++;

    // write revcomp normalized edge
    if (normalize_edges)
    {
        if (kmer_long_2_index == kmer_long_1_index) // except for self loops
            return;

		// since kmer1 -> kmer2 with label (l1,l2), the reverse edge is always: (invert_label(l2),invert_label(l1))
        string label = (string)(kmer2_reversed?"F":"R") + (string)(kmer1_reversed?"F":"R");
        print_edge(kmer_long_2_index,kmer_long_1_index,label,weight);
        nb_edges++;
    }
}


unsigned char node_tight[MAX_TIGHT_SIZE_TIGHT];
abundance_type *node_abundances;
void insert_node(unsigned long id, char node[MAX_TIGHT_SIZE], int k)
{
    unsigned int len_tightk_plus_bytesize = (k / 4) + (k % 4 ? 1 : 0) + TIGHT_KMER_BYTESIZE;
    int thread_num=0;
    large_integer current_kmer;
    // TODO: write kmer_hash_from_sequence directly
    to_tight_DNA(node,node_tight);
    current_kmer=large_kmer_hash_from_tight(node_tight,len_tightk_plus_bytesize);
    insert_kmers_long(current_kmer,thread_num);

}

/*
 *
 *
 *
 *
 *
 */

void process_sequence_with(string sequence, unsigned int k, void operation(large_integer kmer_long, int thread_num), int partition)
{
    unsigned int len_tightk_plus_bytesize = (k / 4) + (k % 4 ? 1 : 0) + TIGHT_KMER_BYTESIZE;

    //#pragma omp parallel default(shared), shared(stdout), firstprivate(scaffold,len_tightk,scaffold_id,k)
    {

        unsigned long i;
        int nb_threads = 1;
        int thread_num=0;//omp_get_thread_num();
        /*
        	if (thread_num==0)
        	printf("hashing scaffold %d using %d threads\n",scaffold_id,nb_threads);*/

        for (i=0; i<sequence.size()-k+1; i++)
        {
            unsigned char kmer_tight[MAX_TIGHT_SIZE];
            string kmer=sequence.substr(i,k);

            // discard kmers if they contain 'N'
            if (kmer.find('N')!=string::npos)
                continue;

            // TODO: write kmer_hash_from_sequence directly
            to_tight_DNA((char *)kmer.c_str(),kmer_tight);

            // normalize edges the k+1 mers (edges)
            if (normalize_edges)
                normalize_kmer(kmer_tight);

            // determine if the kmer should be processed by this partition
            if ((nb_threads!=1 || partitions!=1) && !is_my_kmer_tight_edition(kmer_tight,k+1,thread_num,nb_threads,partition,partitions))
                continue;

            large_integer kmer_long=large_kmer_hash_from_tight(kmer_tight,len_tightk_plus_bytesize);

            operation(kmer_long,thread_num);
        }
    }
}

void process_kmers_from_reads_by_doing(unsigned int k, void operation(large_integer kmer_long, int thread_num), int partition)
{
	FILE *f_reads;
	unsigned long nb_reads=0;
	char read[MAX_TIGHT_SIZE];


	int i;
	for(i=0;i<number_of_read_sets;i++){ // each set of reads
		printf("processing reads (reading %d-mers) from file %s\n",k, reads_file_names[i]);
		f_reads=fopen(reads_file_names[i],"r");
		if (f_reads==NULL)
		{
			printf("couldnt open reads file :%s \n",reads_file_names[i]);
			exit(1);
		}
		// while(!feof) is permitted here, because get_next_read returns feof when done
		while (!feof(f_reads))
		{
			get_next_read(f_reads,read);
			string read_str(read);
			if (read_str.size()>=k)
			{
				process_sequence_with(read_str, k, operation, partition);
				nb_reads++;
			}
		}
		if (nb_reads==0)
		{
			printf("no reads were useful for read file %s, k is too large?\n", reads_file_names[i]);
			exit(1);
		}
		fclose(f_reads);
	} // end each set of reads
}


void write_preliminary_edges(string out_prefix, int weight_threshold, int k, int partition)
{
	int thread_num=0;//omp_get_thread_num();

	sort_kmers(thread_num);

    FILE *dest;
    const char *mode = partition==0?"w":"a";
    dest = fopen((out_prefix+raw_edges_file_suffix).c_str(),mode);

    printf("outputting preliminary edges for partition %d/%d\n",partition+1,partitions);
    unsigned short weight=1;

    huge_vector::iterator it = kmers_long[thread_num]->begin();
    large_integer kmer_to_insert = *(it++);
    large_integer current_kmer = *it;

    huge_vector::iterator end_cached = kmers_long[thread_num]->end();
    for (; it!=end_cached; it++)
    {
        current_kmer = *it;
        if (current_kmer==kmer_to_insert)
        {
            weight++;
        }
        else
        {
            write_preliminary_edge(kmer_to_insert, weight, weight_threshold, k, dest);

            kmer_to_insert=current_kmer;
            weight=1;
        }
    }
    // write the last one
    write_preliminary_edge(current_kmer, weight, weight_threshold, k, dest);

    fclose(dest);
    printf(" wrote %lu edges sequences\n",nb_preliminary_edges);
}

void read_preliminary_edges_and_do(string out_prefix, int k, void operation(large_integer kmer_long_1, large_integer kmer_long_2, unsigned short weight, unsigned int k))
{
    FILE *file;
    file = fopen((out_prefix+raw_edges_file_suffix).c_str(),"r");
    unsigned long counter=0;
    while (1)
    {
        large_integer kmer_long_1,kmer_long_2;
        unsigned short weight;

        fread(&kmer_long_1,sizeof(kmer_long_1), 1, file);
        fread(&kmer_long_2,sizeof(kmer_long_2), 1, file);
        fread(&weight,sizeof(weight), 1, file);

        // why while(!feof()) doesnt work: http://www.gidnetwork.com/b-58.html
        if (feof(file))
            break;

        operation(kmer_long_1,kmer_long_2,weight,k);
        counter++;
    }
    fclose(file);
}

unsigned long nb_nodes = 0;
void write_nodes(string out_prefix, int k)
{
    int thread_num=0;

    sort_kmers(thread_num);

    large_index_integer index = 0;
    large_integer current_kmer ,last_kmer;

    huge_vector::iterator it=kmers_long[thread_num]->begin();
    current_kmer=*it;
    last_kmer = *(it++);
    unsigned short times_kmer_seen=1;

    huge_vector::iterator end_cached = kmers_long[thread_num]->end(); // dunno if it's really time-saving, but i'm caching end iterator..
    for (; it!=end_cached; ++it)
    {
        current_kmer=*it;
        if (current_kmer==last_kmer)
        {
            times_kmer_seen++;
            index++;
        }
        else
        {
            write_node(index,last_kmer,k);
            last_kmer=current_kmer;
            times_kmer_seen=1;
            index++;
            nb_nodes+=1;
        }
    }
    // update info for the last element
    write_node(index,current_kmer,k);
    nb_nodes+=1;
}

void read_nodes_and_do(string nodes_file_name, int k, void operation(unsigned long id, char node[MAX_TIGHT_SIZE], int k))
{
    FILE *file;
    file = fopen(nodes_file_name.c_str(),"r");
    char node[MAX_TIGHT_SIZE], revnode[MAX_TIGHT_SIZE];
    unsigned long id;
    int zero_abundance;

    while (fscanf(file, "%lu\t%s\t%s\t%d\n",&id,node,revnode,&zero_abundance) != EOF)
    {
        operation(id,node,k);
    }
    fclose(file);
}


void debruijn()
{

    int partition;
    for (partition = 0 ; partition < partitions; partition++)
    {
        tables_init();

        // process all the reads and output a list of all (k+1)mers seen in the reads, in the form: (kmer_1,kmer_2, number of times kmer is seen)
        process_kmers_from_reads_by_doing(k+1, insert_kmers_long, partition);
        write_preliminary_edges(out_prefix, weight_threshold, k, partition);

        // clear k+1 mers
        tables_free();
    }
    partitions = 1; // no more partitioning from here

    // start a graph
    string graph_file_name=(out_prefix+graph_file_suffix);
    string nodes_without_abundance_file_name=(out_prefix+nodes_without_abundance_file_suffix);
    string edges_file_name=(out_prefix+edges_file_suffix);
    if (graph_format==0)
    {
        graph_file = fopen(graph_file_name.c_str(),"w");
        fprintf(graph_file,"digraph debuijn {\n");
    }
    else
    {
        nodes_file = fopen(nodes_without_abundance_file_name.c_str(),"w");
        edges_file = fopen(edges_file_name.c_str(),"w");
    }

    // read the previous list and populate all the kmers in a vector, normalizing them on the fly
    tables_init();
    read_preliminary_edges_and_do(out_prefix, k, insert_pairs_of_kmers);

    // output the list of nodes, which are just the normalized kmers
    printf("extracting nodes\n");
    write_nodes(out_prefix,k);
    printf(" wrote %ld nodes (= distinct k-mers)\n",nb_nodes);

    // output the list of edges, which are (index(kmer_1),index(mer_2),label of both mers)
    printf("writing edges\n");
    read_preliminary_edges_and_do(out_prefix, k, write_edge);
    printf(" wrote %ld edges (= distinct (k+1)-mers)\n",nb_edges);

    if (graph_format==0)
    {
        fprintf(graph_file,"}\n");
        fclose(graph_file);
    }
    else
    {
        fclose(nodes_file);
        fclose(edges_file);
    }


    string cmd="mv "+nodes_without_abundance_file_name+" "+out_prefix+nodes_with_abundance_file_suffix;
    system(cmd.c_str());
    printf("finished building the de Bruijn graph in %s%s\n",out_prefix.c_str(),(graph_format==0)?graph_file_suffix.c_str():("["+nodes_with_abundance_file_suffix+"/"+edges_file_suffix+"]").c_str());
}

void print_usage_and_exit(char * name) {

    fprintf (stderr, "Builds the de bruijn graph for one or several set(s) of reads.\n");
    fprintf (stderr, "Nodes are the k-mers and edges correspond to the (k+1)-mers present in the reads.\n");
    fprintf (stderr, "Usage: %s <reads1.fasta/fastq> [<reads2.fasta/fastq> [<reads3.fasta/fastq> [...]]] [-k kmer_size] [-o out_graph] [-g graph_format] [-s solid]\n", name);
    fprintf (stderr, "\t -o output file: Default: test.fa\n");
    fprintf (stderr, "\t -k kmer_size: set the kmer size (max=%d). Default: %d\n",sizeof(large_integer)*4-1,k);
    fprintf (stderr, "\t -s edges_threshold: remove (k+1)mers (i.e. edges) seen < x times. Default: %d\n",weight_threshold);
    fprintf (stderr, "\t -g graph_format: 0 for DOT format, 1 for Kisssplice format. Default: %d\n",graph_format);
    fprintf (stderr, "\t -p nb_part: number of partitions. Default: %d\n",partitions);
    fprintf (stderr, "\t -n: do not automatically add the reverse edge (untested feature). Default: %s\n",normalize_edges?"reverse edges auto-added":"reverse edges not auto-added");
    fprintf (stderr, "\t -h: prints this message and exit\n");
    exit(0);
}

int main(int argc, char *argv[])
{
	char c;
	// get all the read files
	number_of_read_sets=0;
	while(number_of_read_sets+1<argc && argv[number_of_read_sets+1][0]!='-') number_of_read_sets++;
	reads_file_names = (char **) malloc(sizeof(char *)*number_of_read_sets);

	number_of_read_sets=0;
	while(number_of_read_sets+1<argc && argv[number_of_read_sets+1][0]!='-'){
		reads_file_names[number_of_read_sets]=strdup(argv[number_of_read_sets+1]);
		number_of_read_sets++;
	}

	printf("Creating de-bruijn graph for %d files: \n", number_of_read_sets);
	int i;
	for(i=0;i<number_of_read_sets;i++)
		printf("\t %s\n", reads_file_names[i]);

    // parsing the rest of the arguments
    while ((c = getopt (argc-number_of_read_sets, &argv[number_of_read_sets], "k:o:r:g:mt:s:vhnp:")) != -1)
        switch (c)
        {
        case 'k':
            k= atoi(optarg);
            break;
        case 'o':
            out_prefix=string(optarg);
            break;
        case 'g':
            graph_format= atoi(optarg);
            break;
        case 's':
            weight_threshold = atoi(optarg);
            break;
        case 'a':
            dont_compute_nodes_abundance = true;
            break;
        case 'n':
            normalize_edges = false;
            break;
        case 'p':
            partitions = atoi(optarg);
            break;
        case 'h':
            print_usage_and_exit(argv[0]);
            break;
        }

    if (number_of_read_sets==0)
        print_usage_and_exit(argv[0]);


    if (k+1>sizeof(large_integer)*4)
    {
        printf("max value of k on this compiled version: %d\n",sizeof(large_integer)*4-1);
        exit(1);
    }

    debruijn();
}

