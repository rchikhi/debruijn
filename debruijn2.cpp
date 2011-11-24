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
#ifdef _LP64
typedef __uint128_t large_integer;
#else
typedef uint64_t large_integer;
#endif
typedef uint64_t large_index_integer;
typedef vector<large_integer> huge_vector;
typedef int abundance_type;
int max_abundance = numeric_limits<abundance_type>::max();
huge_vector* kmers_long[MAX_THREADS];
string reads_filename;



// default parameters
string out_prefix = "test";
unsigned int k=28;
unsigned int weight_threshold = 0;
short graph_format = 1;
bool dont_compute_nodes_abundance = false;
bool normalize_edges = true;

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
    return kmer_hash_from_tight(kmer_tight,len_tightk_plus_bytesize);
}

void insert_kmers_long(large_integer kmer_long, int thread_num)
{
    // kmers_long cannot hold that many elements.. try another library or make many small vectors?
    kmers_long[thread_num]->push_back(kmer_long);
    if (kmers_long[thread_num]->size()%20000000==0)
        printf("inserted %d (possibly redundant) kmers\n",(int)kmers_long[thread_num]->size());
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
void print_node(large_index_integer index, char *ascii_kmer, char *revcomp_ascii_kmer, int abundance = 0)
{
    if (graph_format==0) // DOT format
        fprintf(graph_file,"%" PRIu64 " [label=\"%s / %s\"];\n",index,ascii_kmer,revcomp_ascii_kmer);
    else if (graph_format==1) // kissplice format
        fprintf(nodes_file,"%" PRIu64 "\t%s\t%s\t%d\n",index,ascii_kmer,revcomp_ascii_kmer,abundance);

}

// output a single edge to a file
void print_edge(large_index_integer id, large_index_integer id2, string label, int weight)
{
    if (graph_format==0) // DOT format
        fprintf(graph_file,"%"PRIu64" -> %"PRIu64" [label=\"%s\" weight=%d];\n",id,id2,label.c_str(),weight);
    else if (graph_format==1) // kissplice format
        fprintf(edges_file,"%"PRIu64"\t%"PRIu64"\t%s\t%d\n",id,id2,label.c_str(),weight);

}

void sort_kmers(int thread_num)
{
    printf(" sorting\n");
    //omptl::
    sort(kmers_long[thread_num]->begin(),kmers_long[thread_num]->end());
}

void insert_pairs_of_kmers(large_integer kmer_long_1, large_integer kmer_long_2, unsigned short weight, unsigned int k)
{
    // yes we are inserting ~2x too many kmers, but:
    // - we cannot recover the information of which (k+1)mer is at the very beginning or very end of a read
    // - since it's without their multiplicites, and coverage is likely to be >2, we are already upperbounded in memory by the first phase
    insert_kmers_long(normalize_kmer(kmer_long_1,k),0);
    insert_kmers_long(normalize_kmer(kmer_long_2,k),0);
}

/*
 *
 *
 *
 */

void write_node(large_index_integer index, large_integer kmer_long, int k, int abundance = 0)
{
    unsigned char kmer_tight[MAX_TIGHT_SIZE_TIGHT];
    char ascii_kmer[k+1];
    char revcomp_ascii_kmer[k+1];

    large_kmer_unhash_to_tight(kmer_long,k,kmer_tight);
    from_tight_DNA(kmer_tight,ascii_kmer);
    strcpy(revcomp_ascii_kmer,ascii_kmer);
    revcomp(revcomp_ascii_kmer,k);

    print_node(index,ascii_kmer,revcomp_ascii_kmer,abundance);
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
void create_edges(large_integer kmer_long_1, large_integer kmer_long_2, unsigned short weight, unsigned int k)
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

    string label = (string)(kmer1_reversed?"R":"F") + (string)(kmer2_reversed?"R":"F");
    print_edge(kmer_long_1_index,kmer_long_2_index,label,weight);
    nb_edges++;

	// write revcomp normalized edge
	if (normalize_edges)
	{
		string label = (string)(kmer1_reversed?"F":"R") + (string)(kmer2_reversed?"F":"R");
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

void write_node_with_abundance(unsigned long id, char node[MAX_TIGHT_SIZE], int k)
{
    unsigned int len_tightk_plus_bytesize = (k / 4) + (k % 4 ? 1 : 0) + TIGHT_KMER_BYTESIZE;
    large_integer current_kmer;
    // TODO: write kmer_hash_from_sequence directly
    to_tight_DNA(node,node_tight);
    current_kmer=large_kmer_hash_from_tight(node_tight,len_tightk_plus_bytesize);
    large_index_integer index=get_kmer_number(current_kmer);
    int abundance = node_abundances[index];
    write_node(id, current_kmer, k, abundance);
}


void count_kmer_abundance(large_integer kmer_long, int thread_num)
{
    large_integer normalized_kmer = normalize_kmer(kmer_long,k);
    if (!is_kmer_present(normalized_kmer))
        return;
    large_index_integer index=get_kmer_number(normalized_kmer);
    if (node_abundances[index]<max_abundance)
        node_abundances[index]++;
}

/*
 *
 *
 *
 *
 *
 */

void process_sequence_with(string sequence, unsigned int k, void operation(large_integer kmer_long, int thread_num))
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

            if (nb_threads!=1 && !is_my_kmer_tight_edition(kmer_tight,k,thread_num,nb_threads,0,1))
                continue;

            // if we normalize edges (because those k+1 mer are edges), we lose the distinction FF/RR and RF/FR 
			if (normalize_edges)
				normalize_kmer(kmer_tight);

            large_integer kmer_long=large_kmer_hash_from_tight(kmer_tight,len_tightk_plus_bytesize);

            operation(kmer_long,thread_num);
        }
    }
}

void process_kmers_from_reads_by_doing(const char * reads_file, unsigned int k, void operation(large_integer kmer_long, int thread_num))
{
    FILE *f_reads;
    unsigned long nb_reads=0;
    char read[MAX_TIGHT_SIZE];

    printf("processing reads (reading %d-mers)\n",k);
    f_reads=fopen(reads_file,"r");
    if (f_reads==NULL)
    {
        printf("couldnt open reads file :%s \n",reads_file);
        exit(1);
    }
    // feof is permitted here because of get_next_read. i should maybe avoid using fasta.c because it's ugly
    while (!feof(f_reads))
    {
        get_next_read(f_reads,read);
        string read_str(read);
        if (read_str.size()>=k)
        {
            process_sequence_with(read_str, k, operation);
            nb_reads++;
        }
    }
	if (nb_reads==0)
	{
		printf("no reads were useful, k is too large?\n");
		exit(1);
	}
    fclose(f_reads);
}


void write_preliminary_edges(string out_prefix, int weight_threshold, int k)
{
    int thread_num=0;//omp_get_thread_num();

    sort_kmers(thread_num);

    FILE *dest;
    dest = fopen((out_prefix+raw_edges_file_suffix).c_str(),"w");

    printf("outputting preliminary edges\n");
    unsigned short weight=1;

    huge_vector::iterator it = kmers_long[thread_num]->begin();
    large_integer kmer_to_insert = *(it++);
    large_integer current_kmer;

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

void compute_nodes_abundance(string nodes_without_abundance_file_name, string nodes_with_abundance_file_name, int k)
{
    if (graph_format!=1)
    {
        printf("cannot compute nodes abundance if graph_format!=1\n");
        exit(1);
    }

    read_nodes_and_do(nodes_without_abundance_file_name,k,insert_node);

    // allocate memory for abundance, it's a pseudo hash table: keys=kmers_long, values=node_abundances
    node_abundances=(abundance_type*)calloc((kmers_long[0]->size()+1),sizeof(abundance_type));

    // can't use the edge information to compute abundance. so we resort to the reads
    process_kmers_from_reads_by_doing(reads_filename.c_str(), k, count_kmer_abundance);

    // read the file again and this time write with abundance
    nodes_file = fopen(nodes_with_abundance_file_name.c_str(),"w");
    read_nodes_and_do(nodes_without_abundance_file_name,k,write_node_with_abundance);
    fclose(nodes_file);

    // remove the file containing nodes without abundance
    unlink(nodes_without_abundance_file_name.c_str());

    free(node_abundances);
}

void debruijn()
{
    tables_init();

    // process all the reads and output a list of all (k+1)mers seen in the reads, in the form: (kmer_1,kmer_2, number of times kmer is seen)
    process_kmers_from_reads_by_doing(reads_filename.c_str(), k+1, insert_kmers_long);
    write_preliminary_edges(out_prefix,weight_threshold,k);

    // clear k+1 mers
    tables_free();

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
    read_preliminary_edges_and_do(out_prefix, k, create_edges);
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

    if (!dont_compute_nodes_abundance && graph_format==1)
    {
        // clear kmers table
        tables_free();
        tables_init();

        printf("computing node abundances\n");
        string nodes_with_abundance_file_name=(out_prefix+nodes_with_abundance_file_suffix);
        compute_nodes_abundance(nodes_without_abundance_file_name,nodes_with_abundance_file_name,k);
    }
	else
		printf("node abundances are not computed\n");

    printf("finished building the de Bruijn graph in %s%s\n",out_prefix.c_str(),(graph_format==0)?graph_file_suffix.c_str():("["+nodes_with_abundance_file_suffix+"/"+edges_file_suffix+"]").c_str());
}

void print_usage_and_exit(char * name) {

    fprintf (stderr, "Builds the de bruijn graph for a set of reads.\n");
    fprintf (stderr, "Nodes are the k-mers and edges correspond to the (k+1)-mers present in the reads.\n");
    fprintf (stderr, "Usage: %s <reads.fasta> [-k kmer_size] [-o out_graph] [-g graph_format] [-s solid]\n", name);
    fprintf (stderr, "\t -o output file: Default: test.fa\n");
    fprintf (stderr, "\t -k kmer_size: set the kmer size (max=%d). Default: %d\n",sizeof(large_integer)*4-1,k);
    fprintf (stderr, "\t -s edges_threshold: remove (k+1)mers (i.e. edges) seen < x times. Default: %d\n",weight_threshold);
    fprintf (stderr, "\t -g graph_format: 0 for DOT format, 1 for Kisssplice format. Default: %d\n",graph_format);
    fprintf (stderr, "\t -a: don't compute abundance of nodes (faster execution). Default: %s\n",dont_compute_nodes_abundance?"abundance is not computed":"abundance is computed");
    fprintf (stderr, "\t -n: don't automatically add reverse edges. Default: %s\n",normalize_edges?"reverse edges are auto-added":"reverse edges are not auto-added");
    fprintf (stderr, "\t -h: prints this message and exit\n");
    exit(0);
}

int main(int argc, char *argv[])
{
    char c;

    while ((c = getopt (argc, argv, "k:o:r:g:mt:s:vhan")) != -1)
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
        case 'h':
            print_usage_and_exit(argv[0]);
            break;
        }

    if (argc - optind != 1)
        print_usage_and_exit(argv[0]);

    reads_filename=string(argv[optind]);

    if (k+1>sizeof(large_integer)*4)
    {
        printf("max value of k on this architecture: %d\n",sizeof(large_integer)*4-1);
        exit(1);
    }

    debruijn();
}
