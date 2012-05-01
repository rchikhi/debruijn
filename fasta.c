#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fasta.h"
#include "kmers.h"

// parse a set of fasta or fastq reads
unsigned long get_next_read(FILE *file, char *read)
{
	unsigned char line[MAX_TIGHT_SIZE];
	char *rv;
	char *p;
    char *read_buffer = read;
	int nextchar=0;
	rv=fgets((char *)line,MAX_TIGHT_SIZE,file); // read comment ('>read00xxxx...\n') and discard it
	if (rv==NULL)
	{
		printf("(reading read header) error or end of file reached!\n");
        exit(1);
	}
    int reading_quality = 0;
    int readlen = 0;
    *read='\0';
    // read the read content (can be multi-line)
	do
    {
        if (reading_quality)
            rv=fgets((char *)line,MAX_TIGHT_SIZE,file); // read quality and discard it
        else
            rv=fgets((char *)read_buffer,MAX_TIGHT_SIZE,file); // read the read content to read_buffer (prevents a strcat)

        if (rv==NULL)
        {
            printf("(reading read content) error, or unexpected end of file reached (current read: %s)!\n",read);
            exit(1);
        }

        // remove trailing \n char; also update readlen
        p = (char *)strchr((char*)read_buffer, '\n');
        if (p) 
        {
            *p = '\0';
            readlen += (int)(p-read_buffer);
        }
        else
            readlen += strlen(read_buffer);

        // remove trailing \r char
        p = (char *)strchr((char*)read_buffer, '\r');
        if (p)
        {
            *p = '\0';
            readlen -= 1;
        }

        // advance the read buffer:
        read_buffer=read+readlen;

        // reads the next '>|@|+' character to: (i) know whether the read continues, (ii) whether there is a quality score, and (iii) in order to induce EOF
        do
            nextchar=fgetc(file);
        while (nextchar == '\n' || nextchar == '\r' ); // also, keep reading newlines.

        if (feof(file)) // end of file (no newline at end of read)
            break;

        if (nextchar == '+')
            reading_quality = 1;
        else
        {
            if ((nextchar != '>') && (nextchar != '@'))
            {
                // the read content continues at the next line
                if (!reading_quality)
                {
                    // if it's not quality, append that read nucleotide
                    readlen++;
                    *read_buffer=nextchar;
                    read_buffer++;
                }
            }
            else // the next line is a new read
                break;
        }

        if (readlen >= MAX_TIGHT_SIZE)
        {
            p = (char *)strchr((char*)line, '\n');
            if (p) *p = '\0';

            printf("error: bad input (line='%s', read='%s').\nthe input file must be reads, of length <%d\n",line,read,MAX_TIGHT_SIZE);
            exit(1);
        }
    }
    while (1);

    //printf("reading read %s of length %d (== %d)\n",read,readlen,strlen(read));

	return readlen;
}


