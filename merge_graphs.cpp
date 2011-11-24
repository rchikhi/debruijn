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
#include <stdint.h>
#include <tr1/unordered_map>

// to print uint64_t
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "fasta.h"
#include "kmers.h"


using namespace std;
using namespace std::tr1;

#define MAX_THREADS 1
#ifdef _LP64
typedef __uint128_t large_integer;
#else
typedef uint64_t large_integer;
#endif

int k = 0;
int nb_abundances[2] = {0,0};

string nodes_file_suffix = ".nodes";
string edges_file_suffix = ".edges";
string correspondence_suffix = ".correspondence";

void check_graph_characteristics(string src_nodes_file_name, int graph_number)
{
    FILE *file;
    file = fopen(src_nodes_file_name.c_str(),"r");
    if (file == NULL)
        goto error;

    char node[MAX_TIGHT_SIZE], revnode[MAX_TIGHT_SIZE];
    unsigned long id;
    int abundance;

    if (fscanf(file, "%lu\t%s\t%s",&id,node,revnode) == 0)
        goto error;

    k=strlen(node);

    // can't use a single fscanf loop because it ignores whitespaces
    while (fgetc(file) == '\t')
    {
        fscanf(file, "%d",&abundance);
        nb_abundances[graph_number]++;
    }

    fclose(file);
    printf("detected that k=%d and nodes of graph %d have %d abundance values\n",k,graph_number,nb_abundances[graph_number]);

    return;

error:
    printf("prefix %s has bad graph format, make sure you specify exactly the [prefix] of ./debruijn2 -g 1 -o [prefix]\n",src_nodes_file_name.c_str());
    exit(1);

}

/*
 *
 * nodes
 *
 */

// store nodes in a hash table, with values = ((id1,(abundance1a,abundance1b,..)),(id2,(abundance2a,abundance2b,..)))
typedef pair<unsigned long,vector<int> > id_and_abundance_list ;
typedef unordered_map<string, pair<id_and_abundance_list,id_and_abundance_list> > graph_nodes_type;
graph_nodes_type graph_nodes;

void hash_nodes_with_abundances(string nodes_file_name, int graph_number)
{
    FILE *file;
    file = fopen(nodes_file_name.c_str(),"r");
    char node[MAX_TIGHT_SIZE], revnode[MAX_TIGHT_SIZE];
    unsigned long id;

    while (fscanf(file, "%lu\t%s\t%s",&id,node,revnode) != EOF)
    {
        vector<int> list_abundance;

        while (fgetc(file) == '\t')
        {
            int abundance;
            fscanf(file, "%d",&abundance);
            list_abundance.push_back(abundance);
        }

        id_and_abundance_list idl(id,list_abundance);

        if (graph_number==0)
            graph_nodes[node].first=idl;
        else
            graph_nodes[node].second=idl;
    }

    fclose(file);
}

void output_all_nodes_with_abundances(string dst_nodes_file_name)
{
    FILE *nodes_file;
    nodes_file=fopen(dst_nodes_file_name.c_str(),"w");

    graph_nodes_type::iterator it;
    unsigned long index = 0;

    for (it=graph_nodes.begin(); it!=graph_nodes.end(); it++)
    {
        char * node = strdup((it->first).c_str());
        fprintf(nodes_file,"%ld\t%s\t",index,node);
        revcomp(node,k);
        fprintf(nodes_file,"%s\t",node);
        free(node);

        vector<int>::iterator it_node;

        // test if the node is present in the reads set (i.e. has at least a non-empty abundance list)
        if (it->second.first.second.size()>0)
            for (it_node = it->second.first.second.begin(); it_node != it->second.first.second.end(); it_node++)
                fprintf(nodes_file,"%d\t",*it_node);
        else
            for (int i=0; i<nb_abundances[0]; i++)
                fprintf(nodes_file,"0\t");

        if (it->second.second.second.size()>0)
            for (it_node = it->second.second.second.begin(); it_node != it->second.second.second.end(); it_node++)
                fprintf(nodes_file,"%d\t",*it_node);
        else
            for (int i=0; i<nb_abundances[1]; i++)
                fprintf(nodes_file,"0\t");

        fprintf(nodes_file,"\n");
        index++;
    }
    fclose(nodes_file);
}

void output_correspondence(string dst_nodes_file_name, string correspondence_suffix)
{
    FILE *first_correspondences_file, *second_correspondences_file;
    first_correspondences_file=fopen((dst_nodes_file_name+".1"+correspondence_suffix).c_str(),"w");
    second_correspondences_file=fopen((dst_nodes_file_name+".2"+correspondence_suffix).c_str(),"w");

    graph_nodes_type::iterator it;
    unsigned long index = 0;

    for (it=graph_nodes.begin(); it!=graph_nodes.end(); it++)
    {
        if (it->second.first.second.size()>0)
        {
            unsigned long old_id_1= it->second.first.first;
            fprintf(first_correspondences_file,"%ld\t%ld\n",old_id_1,index);
        }

        if (it->second.second.second.size()>0)
        {
            unsigned long old_id_2= it->second.second.first;
            fprintf(second_correspondences_file,"%ld\t%ld\n",old_id_2,index);
        }

        index++;
    }
    fclose(first_correspondences_file);
    fclose(second_correspondences_file);
}

void merge_nodes(string src1_nodes_file_name, string src2_nodes_file_name, string dst_nodes_file_name)
{
    printf("reading %s\n",src1_nodes_file_name.c_str());
    hash_nodes_with_abundances(src1_nodes_file_name,0);

    printf("reading %s\n",src2_nodes_file_name.c_str());
    hash_nodes_with_abundances(src2_nodes_file_name,1);

    printf("writing %s\n",dst_nodes_file_name.c_str());
    output_all_nodes_with_abundances(dst_nodes_file_name);
}

/*
 *
 * edges
 *
 */

typedef unordered_map<unsigned long, unsigned long > correspondence_map;
correspondence_map graph_id_correspondence[2];

void load_correspondence(string correspondence_file_name, string correspondence_suffix, int graph_number)
{
    FILE *correspondences_file;
    // string(<int>) doesnt do what i want? i know about the stringstream but really? brb going back to python
    char midfix[5];
    sprintf(midfix,".%d",graph_number+1);
    correspondences_file=fopen((correspondence_file_name+string(midfix)+correspondence_suffix).c_str(),"r");

    unsigned long old_id, new_id;

    while (fscanf(correspondences_file, "%lu\t%lu\n",&old_id,&new_id) != EOF)
        graph_id_correspondence[graph_number][old_id]=new_id;

    fclose(correspondences_file);
}

// oh, just some specialization of tr1::hash for pairs<pairs<A,B>,C>
// follows a solution from there: http://stackoverflow.com/questions/7222143/unordered-map-hash-function-c
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
    std::tr1::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
namespace tr1 {
template<typename S, typename T, typename U> struct hash<pair<pair<S,U>, T> >
{
    inline size_t operator()(const pair<pair<S,U>, T> & v) const
    {
        size_t seed = 0;
        ::hash_combine(seed, v.first.first);
        ::hash_combine(seed, v.first.second);
        ::hash_combine(seed, v.second);
        return seed;
    }
};
}
}

// store edges in a hash table, with keys =(id1,id2,label) values = ((abundance1a,abundance1b,..),(abundance2a,abundance2b,..))
typedef pair< pair<unsigned long, unsigned long>, string> edge_type;
typedef unordered_map<edge_type, pair<vector<int>,vector<int> > > graph_edges_type;
graph_edges_type graph_edges;

void hash_edges_with_abundances(string nodes_file_name, int graph_number)
{
    FILE *file;
    file = fopen(nodes_file_name.c_str(),"r");
    char label[MAX_TIGHT_SIZE];
    unsigned long id1,id2;

    while (fscanf(file, "%lu\t%lu\t%s",&id1,&id2,label) != EOF)
    {
        vector<int> list_abundance;

        while (fgetc(file) == '\t')
        {
            int abundance;
            fscanf(file, "%d",&abundance);
            list_abundance.push_back(abundance);
        }
        unsigned long corresponding_id1 = graph_id_correspondence[graph_number][id1];
        unsigned long corresponding_id2 = graph_id_correspondence[graph_number][id2];

        edge_type edge(pair<unsigned long, unsigned long>(corresponding_id1,corresponding_id2),string(label));

        if (graph_number==0)
            graph_edges[edge].first=list_abundance;
        else
            graph_edges[edge].second=list_abundance;
    }

    fclose(file);
}

void output_all_edges_with_abundances(string dst_edges_file_name)
{
    FILE *edges_file;
    edges_file=fopen(dst_edges_file_name.c_str(),"w");

    graph_edges_type::iterator it;
    unsigned long index = 0;

    for (it=graph_edges.begin(); it!=graph_edges.end(); it++)
    {
        unsigned long id1 = it->first.first.first;
        unsigned long id2 = it->first.first.second;
        string label= it->first.second;
        fprintf(edges_file,"%lu\t%lu\t%s\t",id1,id2,label.c_str());

        vector<int>::iterator it_edge;

        // test if the node is present in the reads set (i.e. has at least a non-empty abundance list)
        if (it->second.first.size()>0)
            for (it_edge = it->second.first.begin(); it_edge != it->second.first.end(); it_edge++)
                fprintf(edges_file,"%d\t",*it_edge);
        else
            for (int i=0; i<nb_abundances[0]; i++)
                fprintf(edges_file,"0\t");

        if (it->second.second.size()>0)
            for (it_edge = it->second.second.begin(); it_edge != it->second.second.end(); it_edge++)
                fprintf(edges_file,"%d\t",*it_edge);
        else
            for (int i=0; i<nb_abundances[1]; i++)
                fprintf(edges_file,"0\t");


        fprintf(edges_file,"\n");
        index++;
    }
    fclose(edges_file);
}

void merge_edges(string src1_edges_file_name, string src2_edges_file_name, string dst_edges_file_name)
{
    printf("reading %s\n",src1_edges_file_name.c_str());
    hash_edges_with_abundances(src1_edges_file_name,0);

    printf("reading %s\n",src2_edges_file_name.c_str());
    hash_edges_with_abundances(src2_edges_file_name,1);

    printf("writing %s\n",dst_edges_file_name.c_str());
    output_all_edges_with_abundances(dst_edges_file_name);
}

/*
 *
 * main
 *
 */

void print_usage_and_exit(char * name) {

    fprintf (stderr, "Merges two de Bruijn graphs constructed using \"./debruijn2 -g 1\". \n");
    fprintf (stderr, "Usage: %s <prefix1> <prefix2> <destination>\n", name);
    fprintf (stderr, "\t -h: prints this message and exit\n");
    exit(0);
}

int main(int argc, char *argv[])
{
    char c;

    while ((c = getopt (argc, argv, "k:o:r:g:mt:s:vha")) != -1)
        switch (c)
        {
        case 'h':
            print_usage_and_exit(argv[0]);
            break;
        }

    if (argc - optind != 3)
        print_usage_and_exit(argv[0]);

    string graph1_filename, graph2_filename, dest_graph_filename;
    graph1_filename=string(argv[optind]);
    graph2_filename=string(argv[optind+1]);
    dest_graph_filename=string(argv[optind+2]);

    check_graph_characteristics(graph1_filename+nodes_file_suffix,0);
    check_graph_characteristics(graph2_filename+nodes_file_suffix,1);

    printf("-NODES-\n");
    merge_nodes(graph1_filename+nodes_file_suffix,graph2_filename+nodes_file_suffix,dest_graph_filename+nodes_file_suffix);

    printf("writing nodes correspondence\n");
    output_correspondence(dest_graph_filename,correspondence_suffix);

    graph_nodes.clear();

    printf("loading nodes correspondence\n");
    load_correspondence(dest_graph_filename,correspondence_suffix,0);
    load_correspondence(dest_graph_filename,correspondence_suffix,1);

    printf("-EDGES-\n");
    merge_edges(graph1_filename+edges_file_suffix,graph2_filename+edges_file_suffix,dest_graph_filename+edges_file_suffix);

    printf("finished merging graphs\n");
}
