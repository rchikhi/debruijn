#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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


