#define TIGHT_KMER_BYTESIZE 1
#define MAX_TIGHT_SIZE (1<<(8*TIGHT_KMER_BYTESIZE))
#define MAX_TIGHT_SIZE_TIGHT (1+(1<<(7*TIGHT_KMER_BYTESIZE)))
#define READ_GROUP_SIZE 50000000L

#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
struct read_t{
	short len;
	char read[MAX_TIGHT_SIZE];
};
typedef struct read_t read_t;

int arbitrary_criterion(unsigned char *s);
unsigned long kmer_hash(char *kmer, unsigned int k);
uint64_t kmer_hash_from_tight(unsigned char *kmer_tight, unsigned int len_tightk_plus_bytesize);
void kmer_unhash_to_tight(uint64_t kmer_hashed, unsigned int k, unsigned char* result_kmer_tight);
void test_tight_kmers(void);
void extract_tight_kmer(unsigned char *src, unsigned int k, unsigned int j, unsigned char *kmer);
void extract_kmer(char *src, int k, int j, char *kmer);
int from_tight_DNA(unsigned char *s, char *decoded);
int to_tight_DNA(char *s,unsigned char *encoded);
void revcomp_tight(unsigned char s[]);
void revcomp(char s[],int len);
int is_palindromic(char *s);
int is_low_complexity(unsigned char *s);
int is_my_kmer_tight_edition(unsigned char * kmer, int k, int thread_num, int nb_threads, int node_num, int nb_nodes);
int get_thread_num_for_my_kmer_tight_edition(unsigned char * kmer, int k, int nb_threads);
int kmer_has_Ns(char * read, int start, int k);
void right_modify_tight_kmer(unsigned char *src, int ii, unsigned char *kmer); // DL
void left_modify_tight_kmer(unsigned char *src, int ii, unsigned char *kmer); // DL
void correct_tight_read(unsigned char *read_tight, int pos, int val); // DL

#ifdef __cplusplus
}
#endif

