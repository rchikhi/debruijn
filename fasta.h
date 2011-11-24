#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void determine_min_readlen(unsigned long *readlen_histogram, unsigned long nb_reads, uint32_t *min_readlen, uint32_t *max_readlen );
unsigned long get_next_read(FILE *file, char *read);
unsigned long get_next_read_and_comment(FILE *file, char *read, char *comment);
void print_histogram(unsigned long * histogram, int size, int print_total, char *legend);
void get_min_max_readlens(char *filename1, char *filename2, uint32_t *min_readlen, uint32_t *max_readlen, unsigned long *nb_reads);
void get_next_scaffold(FILE *file, char *scaffold, char *buffer, unsigned int buffer_length);
void compute_mean_stdev(unsigned long *histogram, int max_x, int *mean, int *stdev);

#ifdef __cplusplus
}
#endif

