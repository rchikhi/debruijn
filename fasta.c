#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "kmers.h"
#include "fasta.h"

// very basic!
unsigned long get_next_read(FILE *file, char *read)
{
	unsigned char line[MAX_TIGHT_SIZE];
	char *rv;
	char *p;
	int nextchar=0;
	rv=fgets((char *)line,MAX_TIGHT_SIZE,file); // read comment ('>read00xxxx...\n')
	if (rv==NULL)
	{
		printf("(reading read header) error or end of file reached!\n");
		return 0;
	}
	rv=fgets((char *)read,MAX_TIGHT_SIZE,file); //
	if (rv==NULL)
	{
		printf("(reading read content) error or end of file reached!\n");
		return 0;
	}
	p = (char *)strchr((char*)read, '\n');
        if (p) *p = '\0';
	nextchar=fgetc(file); // cheat, reads the next '>' character in order to induce EOF
	if (nextchar!='>' && !feof(file))
	{
		printf("error: bad input (line='%s', read='%s').\n btw this program cannot process genomes. the input file must be reads, of length <%d, one line per read\n",line,read,MAX_TIGHT_SIZE);
		exit(1);
	}
	if ((unsigned long)(p-read)>MAX_TIGHT_SIZE)
	{
		printf("error loading read %s (%p, end of read %p) deduced len %ld\n",read,read,p,(unsigned long)(p-read));
	}
	return (unsigned long)(p-read); // readlen
}

// get the fasta comment also
// WARNING : suppose the sequence on a single texte line

unsigned long get_next_read_and_comment(FILE *file, char *read, char *comment)
{
	char *rv;
	char *p;
	rv=fgets((char *)comment,MAX_TIGHT_SIZE,file);
	if (rv==NULL) return 0;
	rv=fgets((char *)read,MAX_TIGHT_SIZE,file);
	if (rv==NULL) return 0;
	p = (char *)strchr((char*)read, '\n');
        if (p) *p = '\0';
	return (unsigned long)(p-read); // readlen
}


// select readlen such that most of the reads are longer than this length
// btw nb_reads is the number of reads added, not total number of reads
void determine_min_readlen(unsigned long *readlen_histogram, unsigned long nb_reads, uint32_t *min_readlen, uint32_t *max_readlen )
{
	unsigned long accumulated=0;
	// uncomment for another way of detemrining, less precise.
	//unsigned long threshold=(nb_reads*99)/100;
	unsigned long threshold=(10000);
	int i;
	if (min_readlen)
		*min_readlen=0;
	if (max_readlen)
		*max_readlen=0;

	//for (i=max_readlen;i>0;i--)
	for (i=0;i<MAX_TIGHT_SIZE;i++)
	{
		accumulated+=readlen_histogram[i];
		if (accumulated>threshold && min_readlen!=NULL && *min_readlen==0)
			*min_readlen=i;
		if (accumulated>(nb_reads*0.9) && max_readlen!=NULL && *max_readlen==0)
			*max_readlen=i;
	}
	if (*min_readlen==0)
		*min_readlen=*max_readlen;
	if ( (min_readlen!=NULL && *min_readlen==0) || (max_readlen!=NULL && *max_readlen==0))
	{
		printf("couldn't find either min_readlen (%d) or max_readlen (%d)!",*min_readlen,*max_readlen);
		exit(1);
	}
}
void get_min_max_readlens(char *filename1, char *filename2, uint32_t *min_readlen, uint32_t *max_readlen, unsigned long *nb_reads)
{
	unsigned long *readlen_histogram;
	*nb_reads=0L;
	uint32_t cur_readlen;
	static char read[MAX_TIGHT_SIZE];

	readlen_histogram=(unsigned long *)calloc(MAX_TIGHT_SIZE,sizeof(unsigned long));

	printf("getting min/max readlen from reads %s %s\n",filename1,filename2);

	FILE *files[2];
	int nb_files=0;
	files[0]=fopen(filename1,"r");
	if (files[0]!=NULL)
		nb_files++;
	else
		printf("error opening file %s\n",filename1);
	files[1]=fopen(filename2,"r");
	if (files[1]!=NULL)
		nb_files++;
	else
		if (strlen(filename2)!=0)
			printf("error opening file %s\n",filename2);
	int file_index;
	for (file_index=0;file_index<nb_files;file_index++)
	{
		FILE *file=files[file_index];
		if (file==NULL)
		{
			printf("null file!\n");
			return;
		}

		while (!feof(file))
		{
			cur_readlen=get_next_read(file, read);
			readlen_histogram[cur_readlen]++;
			(*nb_reads)++;
		}
		fclose(file);
	}

	determine_min_readlen(readlen_histogram,*nb_reads,min_readlen,max_readlen);

	print_histogram(readlen_histogram, MAX_TIGHT_SIZE, 1, (char *)"Readlen histogram");
	printf("\nDetermined minimal acceptable readlen: %d, maximal acceptable readlen: %d\n",*min_readlen, *max_readlen);
	free(readlen_histogram);
}


void print_histogram(unsigned long * histogram, int size, int print_total, char *legend)
{
	// parameters
	int nb_columns=20;

	int min_nonzero_x=size;
	int max_nonzero_x=0;
	unsigned long min_nonzero_y=LONG_MAX;
	unsigned long max_nonzero_y=0;
	unsigned long i,j;

	// get bounding box
	for (i=0;i<(unsigned long)size;i++)
	{
		if (min_nonzero_x==size && histogram[i]!=0)
			min_nonzero_x=i;

		if (histogram[i]!=0)
			max_nonzero_x=i;

		if (histogram[i]>max_nonzero_y)
			max_nonzero_y=histogram[i];

		if (histogram[i] != 0 && histogram[i]<min_nonzero_y)
			min_nonzero_y=histogram[i];
	}
	printf("histogram min_x=%d max_x=%d min_y=%ld max_y=%ld\n",min_nonzero_x,max_nonzero_x,min_nonzero_y,max_nonzero_y);
	
	// compute mean/stdev for better printing
	int mean,stdev;
	compute_mean_stdev(histogram, max_nonzero_x, &mean, &stdev);
	min_nonzero_x=mean-4*stdev;
	if (min_nonzero_x<0)
		min_nonzero_x=0;
	if (mean+4*stdev<max_nonzero_x)
		max_nonzero_x=mean+4*stdev;

	// now print it
	int delta_x=stdev/2;
	if (delta_x<=0)
		delta_x=1;
	int delta_y=(max_nonzero_y-min_nonzero_y)/nb_columns;
	if (delta_y<=0)
		delta_y=1;

	printf("%s (mean=%d stdev=%d)\n",legend,mean,stdev);
	for (i=(unsigned long)min_nonzero_x;i<=(unsigned long)max_nonzero_x;i+=(unsigned long)delta_x)
	{
		unsigned long mean_for_delta=0;
		for (j=i;j<i+delta_x;j++)
		{
			if (j>(unsigned long)max_nonzero_x)
				continue;
			mean_for_delta+=histogram[j];
		}
		if (print_total==0)
			mean_for_delta/=delta_x;
		printf("[%2ld-%2ld]:\t%ld\t",i,i+delta_x-1,mean_for_delta);
		for (j=min_nonzero_y;j<mean_for_delta;j+=delta_y*delta_x)
			printf("*");
		printf("\n");
	}
}

void compute_mean_stdev(unsigned long *histogram, int max_x, int *mean, int *stdev)
{
	uint64_t total_kmers=0;
	uint64_t weighted_total=0;
	int i;
        *mean=0;
        *stdev=0;
	for (i=1;i<=max_x;i++)
	{
		total_kmers+=histogram[i];
		weighted_total+=i*histogram[i];
        }
        if (total_kmers==0)
                return;
        *mean=(int)(weighted_total/total_kmers);
        
	uint64_t tmp_stdev=0;
	for (i=1;i<max_x;i++)
		tmp_stdev+=histogram[i]*(i-*mean)*(i-*mean);
	tmp_stdev/=total_kmers;
	*stdev=(int)sqrt(tmp_stdev)+1;
}

// very basic!
void get_next_scaffold(FILE *file, char *scaffold, char *buffer, unsigned int buffer_length)
{
	char *rv;
	char *p;
	int nextchar=0;
	rv=fgets((char *)buffer,buffer_length,file); // read comment ('>read00xxxx...\n')
	if (rv==NULL)
	{
		printf("(reading seq header) error or end of file reached!\n");
		exit(1);
	}
	rv=fgets((char *)scaffold,buffer_length,file); // read beginning of scaffold
	if (rv==NULL)
	{
		printf("(reading seq content) error or end of file reached!\n");
		exit(1);
	}
	p = (char *)strchr((char*)scaffold, '\n');
        if (p) *p = '\0';
	nextchar=fgetc(file); // cheat, reads the next '>' character in order to induce EOF
				// for the next scaffold; it'll just skip the name
	// TODO: that code is very slow (strcat's), should use onelinefasta.pl before
	while (nextchar!='>' && !feof(file))
	{
		buffer[0]=nextchar;
		buffer[1]=0;
		strcat(scaffold,buffer);
		rv=fgets((char *)buffer,buffer_length,file); // get the rest of the line
		p = (char *)strchr((char*)buffer, '\n');
	        if (p) *p = '\0';
		strcat(scaffold,buffer);
		if (strlen(scaffold)>buffer_length-10)
			printf("(reading seq content) scaffold too long!\n");
		nextchar=fgetc(file);
		if (nextchar=='\n') // empty line reached, read again
			nextchar=fgetc(file);
			
	}
	return;
}
