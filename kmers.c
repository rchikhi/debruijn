#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "kmers.h"
#include "fasta.h"
#include <math.h>

#define using_precomputed_tables 1
#ifdef using_precomputed_tables
#include "tight_tables.h"
#endif


void revcomp(char s[], int len)
{
#define CHAR_REVCOMP(a,b) {switch(a){\
	case 'A': b='T';break;case 'C': b='G';break;case 'G': b='C';break;case 'T': b='A';break;default: b=a;break;}}
		  int i;
		  unsigned char t;
		  for (i=0;i<len/2;i++)
		  {
			  t=s[i];
			  CHAR_REVCOMP(s[len-i-1],s[i]);
			  CHAR_REVCOMP(t,s[len-i-1]);
		  }
		  if (len%2==1)
			  CHAR_REVCOMP(s[len/2],s[len/2]);

}

void revcomp_tight(unsigned char s[])
{

	unsigned int i;
	unsigned int len;
	if (TIGHT_KMER_BYTESIZE==1)
		len=s[0];
	if (TIGHT_KMER_BYTESIZE==2)
		len=s[0]+(0xFF*s[1]);
	unsigned char cur1,cur2,next1,next2;
	next1=0;next2=0;
	for (i=0;i<len/2;i++)
	{
		cur1=(s[TIGHT_KMER_BYTESIZE+(i/4)]>>(2*(i%4)))&3;
		cur2=(s[TIGHT_KMER_BYTESIZE+((len-i-1)/4)]>>(2*((len-i-1)%4)))&3;
		next1+=(3-cur2)<<(2*(i%4));
		next2+=(3-cur1)<<(2*((len-i-1)%4));
		if (i%4==3)
		{
			s[TIGHT_KMER_BYTESIZE+(i/4)]=next1;
			next1=0;
		}
		if ((len-i-1)%4==0)
		{
			s[TIGHT_KMER_BYTESIZE+((len-i-1)/4)]=next2;
			next2=0;
		}
	}
	if (len%2==1)
		next1+=((3-((s[TIGHT_KMER_BYTESIZE+((len/2)/4)])>>(2*((len/2)%4))))&3)<<(2*((len/2)%4));
	if (next1+next2 || (len%2==0 && len%8!=0) || len%2==1) // fixme: horrible hack
		s[TIGHT_KMER_BYTESIZE+(len/2)/4]=next1+next2;
}

int to_tight_DNA(char *s, unsigned char *encoded)
{
	int i,resi,idx;
	unsigned char cur;
	cur=0;
	idx=0;
	resi=TIGHT_KMER_BYTESIZE;
	i=0;
	if (s==NULL || encoded==NULL)
	{
		printf("error: one tightDNA parameter is null\n");
		exit(1);
	}
	while (s[i]!=0)
	{
		switch (s[i])
		{
			case 'C': cur+=1<<idx; break;
			case 'G': cur+=2<<idx; break;
			case 'T': cur+=3<<idx; break;
			case 'A': break;
			case 'N': cur+=(rand()%4)<<idx; break; 
			// creative way to treat N's.. can i think of a better one?
		}

		idx+=2;
		if (idx>=8)
		{
			idx=0;
			encoded[resi++]=cur;
			cur=0;
		}
		if (resi>MAX_TIGHT_SIZE_TIGHT)
		{
			printf("overflow in tightDNA! %d\n",resi);
			exit(1);
		}
		i++;
	}
	unsigned int dnalen=i;
	//unsigned int lenres=(dnalen/4)+(dnalen%4?1:0);
	if (dnalen>MAX_TIGHT_SIZE)
	{
		printf("cannot call totightDNA with this string length (%d)",dnalen);
		return -1;
	}
	if (TIGHT_KMER_BYTESIZE==1)
		encoded[0]=dnalen;
	else
		if (TIGHT_KMER_BYTESIZE==2)
		{
			encoded[0]=(dnalen%0xFF);
			encoded[1]=(dnalen/0xFF);
		}

	if (idx!=0)
		encoded[resi++]=cur;
	encoded[resi]=0;
	return 0;
}

void construct_tight_DNA_tables(void)
{
	unsigned int s;
	FILE *tables;
	unsigned char tight_test[3];
	char decoded_test[5];
	tables=fopen("tight_tables.h","w");
	if (tables == NULL)
	{
		printf("can't open tight_tables.h\n");
		exit(1);
	}
	tight_test[0]=4;
	tight_test[2]='\0';
	fprintf(tables,"const char *tight_DNA_table[256]={\\\n");
	for (s=0;s<=0xFF;s++)
	{
		tight_test[1]=s;
		from_tight_DNA(tight_test,decoded_test);
		fprintf(tables,"%s\"%s\"",s!=0?",\n":"",decoded_test);
	}
	fprintf(tables,"};\n");
	fclose(tables);
}

int from_tight_DNA(unsigned char *s, char *decoded)
{
	unsigned int dnalen;
	if (TIGHT_KMER_BYTESIZE==1)
		dnalen=s[0];
	if (TIGHT_KMER_BYTESIZE==2)
		dnalen=s[0]+(s[1]*0xFF);

	//unsigned int lenres=(dnalen/4)+(dnalen%4?1:0);
	int i;
	if (s==NULL || decoded==NULL)
	{
		printf("error: one tightDNA parameter is null\n");
		exit(1);
	}

	if (dnalen>MAX_TIGHT_SIZE)
	{
		printf("cannot call fromtightDNA with this string length (%d)",dnalen);
		return -1;
	}
	decoded[dnalen]=0;

	#if using_precomputed_tables
	{
		int j=0;
		const char *decoded_bytes;
		int nb_regular_copy=dnalen/4;
		int nb_remaining=dnalen%4;
		/*
		 * full bytes in the tight sequence
		 */
		for (i=0;i<nb_regular_copy;i++)
		{
			decoded_bytes=tight_DNA_table[s[TIGHT_KMER_BYTESIZE+i]];
			decoded[j++]=decoded_bytes[0];
			decoded[j++]=decoded_bytes[1];
			decoded[j++]=decoded_bytes[2];
			decoded[j++]=decoded_bytes[3];
		}
		/*
		 * last bits of the tight sequence
		 */
		decoded_bytes=tight_DNA_table[s[TIGHT_KMER_BYTESIZE+nb_regular_copy]];
		for (i=0;i<nb_remaining;i++)
		{
			decoded[j++]=decoded_bytes[i];
		}
	}
	#else
	unsigned char cur;
	{
		for (i=0;i<dnalen;i++)
		{
			cur=(s[TIGHT_KMER_BYTESIZE+(i/4)]>>(2*(i%4)))&3;
			switch (cur)
			{
				case 0:	decoded[i]='A';break;
				case 1:	decoded[i]='C';break;
				case 2:	decoded[i]='G';break;
				case 3:	decoded[i]='T';break;
			}
		}
	}
	#endif
	return 0;
}

void extract_kmer(char *src, int k, int j, char *kmer)
{
	int l;
	for (l=0;l<k;l++)
		kmer[l]=src[j+l];
	kmer[k]=0;
}

void extract_tight_kmer(unsigned char *src, unsigned int k, unsigned int j, unsigned char *kmer)
{
	unsigned int resi,i,idx;
	unsigned char cur;
	resi=TIGHT_KMER_BYTESIZE;
	if (TIGHT_KMER_BYTESIZE==1)
		kmer[0]=k;
	else
	{
		kmer[0]=k%0xFF;
		kmer[1]=k/0xFF;
	}
	idx=0;
	cur=0;

	for (i=j;i<j+k;i++)
	{
		cur+=(src[TIGHT_KMER_BYTESIZE+(i/4)]>>(2*(i%4))&3)<<idx;
		idx+=2;
		if (idx>=8)
		{
			idx=0;
			kmer[resi++]=cur;
			cur=0;
		}
	}
	if (idx!=0)
		kmer[resi++]=cur;
}

/*
 * well it's a perfect hash
 */
unsigned long kmer_hash(char *kmer, unsigned int k)
{
	unsigned long res;
	unsigned long idx,i;
	idx=0;
	res=0;
	for (i=0;i<k;i++)
	{
		if (kmer[i]=='C') res+=1L<<idx;
		if (kmer[i]=='G') res+=2L<<idx;
		if (kmer[i]=='T') res+=3L<<idx;
		idx+=2;
	}
	return res;
}

uint64_t kmer_hash_from_tight(unsigned char *kmer_tight, unsigned int len_tightk_plus_bytesize)
{

	uint64_t res;
	uint64_t i;
	res=0;
	for (i=0;i<len_tightk_plus_bytesize-TIGHT_KMER_BYTESIZE;i++)
	{
		res+=(((uint64_t)(kmer_tight[TIGHT_KMER_BYTESIZE+i]))<<(i*8));
	}
	return res;
}

void kmer_unhash_to_tight(uint64_t kmer_hashed, unsigned int k, unsigned char* result_kmer_tight)
{
	uint64_t i;
	unsigned int len_tightk = (k / 4) + (k % 4 ? 1 : 0);

	// set the size
	if (TIGHT_KMER_BYTESIZE==1)
		result_kmer_tight[0]=k;
	else
		if (TIGHT_KMER_BYTESIZE==2)
		{
			result_kmer_tight[0]=(k%0xFF);
			result_kmer_tight[1]=(k/0xFF);
		}

	// decode raw bits
	for (i=0;i<len_tightk;i++)
		result_kmer_tight[TIGHT_KMER_BYTESIZE+i]=(kmer_hashed>>(i*8))&0xFF;
}

int arbitrary_criterion(unsigned char *s)
{
	// whenever a kmer k1 is such that revcomp(k1) < k1, hash revcomp(k1), else (k1 < revcomp(k1)) hash k1
	// the order "<" is given here: (it's lexicographical)
	// idea: (original) A..A <-> (revcomp) T..T so let's choose A..A
	// A..T <-> A..T can't decide, let's see what's next inside
	// A..C <-> G..T let's pick (A,C)<(G,T)
	// A..G <-> C..T well let's say A<C<G<T, so A..G wins
	// G..A <-> T..C so G..A wins
	// here A=0,C=1,G=2,T=3
	unsigned int i;
	unsigned int len;
	if (TIGHT_KMER_BYTESIZE==1)
		len=s[0];
	if (TIGHT_KMER_BYTESIZE==2)
		len=s[0]+(s[1]*0xFF);

	unsigned char cur1,cur2;
	i=0;
	while (i<len/2)
	{
		cur1=(s[TIGHT_KMER_BYTESIZE+(i/4)]>>(2*(i%4)))&3; // i-th nucleotide of s
		cur2=(s[TIGHT_KMER_BYTESIZE+((len-i-1)/4)]>>(2*((len-i-1)%4)))&3; // len(s)-i-th nucleotide of s
		//printf("%d %d\n",cur1,cur2);
		if (cur1==cur2) return (cur1>1); // A..A T..T G..G C..C cases
		if (cur1!=(3-cur2)) // every other case:
			return !(cur1==0 || cur2==0); // the 'A' always wins =)
		i++;
	}
	// only one last middle nucleotide to disambiguate
	if (len%2==1)
	{
		cur1=(((s[TIGHT_KMER_BYTESIZE+((len/2)/4)])>>(2*((len/2)%4))))&3;
		return cur1>1;
	}

	// if we land here, kmer and its revcomp are equal..
	return 2;
}

int is_palindromic(char *s)
{
	unsigned char encoded[MAX_TIGHT_SIZE];
	to_tight_DNA(s,encoded);
	return arbitrary_criterion(encoded)==2;
}

// detects if the tight kmer contains AAAAAAA,CCCCCCCC,GGGGGGG or TTTTTTT at a 4-nucleotides aligned position
// remark: kmers are tight!
// (deprecated) not using it anymore, hurts dbg-based assembly
int low_complexity_too_stringent(unsigned char *s)
{
	unsigned int len=s[0];
	int j;
	unsigned char sameletter[4]={0,0x55,0xAA,0xFF};
	unsigned int i=0;
	while (i<len/4-1)
	{
		for (j=0;j<4;j++)
			if (s[TIGHT_KMER_BYTESIZE+i]==sameletter[j])
				if (s[TIGHT_KMER_BYTESIZE+i+1]==sameletter[j])
					return 1;
		i++;
	}
	return 0;
}

// returns 1 if the sequence contains only <=2 different nucleotides (e.g. AAAAAAAAAAA, ATATATATAT, AAAAGGAAAAA)
int is_low_complexity(unsigned char *s)
{
	unsigned int dnalen;
	if (TIGHT_KMER_BYTESIZE==1)
		dnalen=s[0];
	if (TIGHT_KMER_BYTESIZE==2)
		dnalen=s[0]+(s[1]*0xFF);

	int count[4]={0};
	unsigned int i;
	unsigned char cur;
	for (i=0;i<dnalen;i++)
	{
		cur=(s[TIGHT_KMER_BYTESIZE+(i/4)]>>(2*(i%4)))&3;
		count[cur]++;
	}
	int nb_zero_count=0;
	for (i=0;i<4;i++)
		if (count[i]==0)
			nb_zero_count++;
	return (nb_zero_count>=2);
}

// test procedure to verify that kmers-related procedures work well
void test_tight_kmers()
{
	// setting a random read
	int i=0,max_i;
	char src[500];
	srand ( time(NULL) );
	//max_i=rand()%10+100;
	max_i=100;
	for (i=0;i<max_i;i++)
	{
		switch (rand()%4){case 0: src[i]='A';break;case 1: src[i]='C';break;
		case 2: src[i]='T';break;case 3: src[i]='G';break;}
	}
	src[max_i]=0;

// sample problematic cases
//	strcpy(src,"AAAAAAGACGGGGGGGGGGGGGGGGGCAGTAAAAAAA");
//	strcpy(src,"AAAAAATCTTATTCAGCAGTTTTTTGATGAGGTCGTAAA");
//	strcpy(src,"CTTGCGCTAATTTTTTGTCATCAAACCTGT");
//	strcpy(src,"ATTGATGCATTTTAACCTTAACCGTTTGGTTAGGGTA");
	unsigned char dest[500];
	unsigned char dest2[500];
	unsigned char tmp[500];
	int kmer_size=20;
	to_tight_DNA(src,tmp);
	revcomp_tight(tmp);
	from_tight_DNA(tmp,(char *)dest);
	printf("original seq: %s (len:%d) revcomp seq: %s\n",src,tmp[0],dest);
	if (tmp[0]!=strlen((const char*)dest) || strlen((const char*)dest)!=strlen((const char*)src))
		printf("*** LENGTH ERROR!****\n");

	extract_kmer(src,kmer_size,6,(char *)dest);
	printf("\noriginl: %s\n",dest);
	revcomp((char *)dest,kmer_size);
	printf("\nrevcomp: %s\n",dest);
	to_tight_DNA(src,tmp);
	extract_tight_kmer(tmp,kmer_size,6,dest2);
	revcomp_tight(dest2);
	from_tight_DNA(dest2,(char *)tmp);
	printf("\ndecoded: %s\n",tmp);

	if (strcmp((char *)tmp,(char *)dest)!=0)
		printf("*** error: not agreeing!\n");
	printf("arbitrary_criterion: %d (revcomp) should be equal to..\n",arbitrary_criterion(dest));
	revcomp_tight(dest);
	printf("arbitrary_criterion: %d (!original)\n",!arbitrary_criterion(dest));

	// testing kmer_hash and kmer_unhash
	extract_kmer(src,kmer_size,6,(char *)dest);
	to_tight_DNA((char *)dest,tmp);
	int len_tightk_p_b = (kmer_size/ 4) + (kmer_size% 4 ? 1 : 0)+TIGHT_KMER_BYTESIZE;
	uint64_t hash=kmer_hash_from_tight(tmp,len_tightk_p_b);
	printf("kmer_hash_from_tight of %s: %llX\n",dest,hash);
	kmer_unhash_to_tight(hash,kmer_size,tmp);
	from_tight_DNA(tmp,(char *)dest2);
	printf("kmer_unhash_from_tight %s\n",dest2);
	if (strcmp((char *)dest,(char *)dest2)!=0)
		printf("*** error: not agreeing!\n");
}

int is_my_kmer_tight_edition(unsigned char * kmer, int k, int thread_num, int nb_threads, int node_num, int nb_nodes)
{
	/*
	 * kmer -> long -> % nb_nodes -> % nb_threads, perfect partition of the kmer space into nodes/threads
	 */
	unsigned long kmer_hash=0;
	int len_tightk = (k / 4) + (k % 4 ? 1 : 0);
	int i;

	// compute a kmer hash
	for (i=0;i<len_tightk;i++)
		kmer_hash+=(kmer[TIGHT_KMER_BYTESIZE+i]);

	// testing if this node is concerned (note that if nb_nodes=1, all kmers should pass this test)
	int is_my_node=((int)(kmer_hash%nb_nodes)==node_num);
	if (!is_my_node)
		return 0;

	// now testing if this thread is concerned
	kmer_hash/=nb_nodes;
	return ((int)(kmer_hash%nb_threads)==thread_num);
}

int get_thread_num_for_my_kmer_tight_edition(unsigned char * kmer, int k, int nb_threads)
{
	/*
	 * kmer -> long -> % nb_threads
	 */
	unsigned long kmer_hash=0;
	int len_tightk = (k / 4) + (k % 4 ? 1 : 0);
	int i;

	// compute a kmer hash
	for (i=0;i<len_tightk;i++)
		kmer_hash+=(kmer[TIGHT_KMER_BYTESIZE+i]);

	return (int)(kmer_hash%nb_threads);
}

int kmer_has_Ns(char * read, int start, int k)
{
	int i;
	for (i=0;i<k;i++)
		if (read[start+i]=='N')
			return 1;
	return 0;
}


//// kmer modification functions used in read corrector

void right_modify_tight_kmer(unsigned char *src, int ii, unsigned char *kmer)
{
  int i,j,y;
  int x[4];
  int dnalen;
  dnalen = src[0]; // give the length of the kmer
  for (i=0; i<TIGHT_KMER_BYTESIZE+1+(dnalen/4); i++) kmer[i]=src[i];
  i = dnalen-1;
  y = kmer[TIGHT_KMER_BYTESIZE+(i/4)];
  for (j=0; j<4; j++) x[j] = (y>>(2*j))&3;
  j =(i%4)&3;  
  x[j] = (x[j]+ii)%4;
  y=0;
  for (j=0; j<4; j++) y = y + (x[j]<<(2*j));
  kmer[TIGHT_KMER_BYTESIZE+(i/4)] = y;
}

void left_modify_tight_kmer(unsigned char *src, int ii, unsigned char *kmer)
{
  int dnalen, i,y;
  dnalen = src[0]; // give the length of the kmer
  for (i=0; i<TIGHT_KMER_BYTESIZE+1+(dnalen/4); i++) kmer[i]=src[i];
  y = kmer[TIGHT_KMER_BYTESIZE];
  y = (((y&3)+ii)%4)|(y&0xFC);
  kmer[TIGHT_KMER_BYTESIZE] = y;
}

void correct_tight_read(unsigned char* read_tight, int pos, int val)
{
  int j,y;
  int x[4];

  y = read_tight[TIGHT_KMER_BYTESIZE+(pos/4)];
  for (j=0; j<4; j++) x[j] = (y>>(2*j))&3;
  j =(pos%4)&3;  
  x[j] = (x[j]+val)%4;
  y=0;
  for (j=0; j<4; j++) y = y + (x[j]<<(2*j));
  read_tight[TIGHT_KMER_BYTESIZE+(pos/4)] = y;
}
