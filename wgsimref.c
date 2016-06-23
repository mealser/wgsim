#include <zlib.h>  
#include <stdio.h>  
#include "kseq.h"  

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>

// STEP 1: declare the type of file handler and the read() function  
KSEQ_INIT(gzFile, gzread)  



char* substring(const char* str, int start, int end) 
{ 
  if (str == 0 || strlen(str) == 0 || strlen(str) < start || strlen(str) < (start+end)) 
    return 0; 

  return strndup(str + start, end-start+1); 
} 

void prepare_genome(const char *fn, const char *fn_new)
{
	FILE *fpnew;
	kseq_t *ks;
	gzFile fp_fa;
	uint64_t tot_len;
	int l, n_ref;

	fpnew = fopen(fn_new, "w");
	fp_fa = gzopen(fn, "r");
	ks = kseq_init(fp_fa);
	tot_len = n_ref = 0;
	fprintf(stderr,"[%s] calculating the total length of the reference sequence...\n", __func__);
	while ((l = kseq_read(ks)) >= 0) {
		fprintf(fpnew,"%s", ks->seq.s); 
		tot_len += l;
		++n_ref;
	}
	
	fprintf(stderr,"[%s] %d sequences, total length: %llu\n", __func__, n_ref, (long long)tot_len);
	kseq_destroy(ks);
	gzclose(fp_fa);
	fclose(fpnew);
}


char* read_genome(FILE *genomenew,FILE *fpnew, int start, int end)
{
	char* S="hi";
	char mystring [end-start+1];
	fseek ( genomenew , start-1 , SEEK_SET );
	fgets (mystring , end-start+1 , genomenew);
	fprintf(stderr,"[%s] Extracting the reference sequence...%d\n%d\n%d\n %s\n", __func__,start, end, end-start+1, mystring);
	return S;
}


/*
char* read_genome(const char *fn,FILE *fpnew, int start, int end)
{
	kseq_t *ks;
	gzFile fp_fa;
	uint64_t tot_len;
	int l, n_ref;
	char* S="hi";
	char *RefSeq;
//char *dest= malloc(end-start+1);
	
	l = end-start+1;


	fp_fa = gzopen(fn, "r");
	ks = kseq_init(fp_fa);
	tot_len = n_ref = 0;
	fprintf(stderr,"[%s] calculating the total length of the reference sequence...\n", __func__);
	while ((l = kseq_read(ks)) >= 0) {
		fprintf(fpnew,"%s", ks->seq.s); 
		tot_len += l;
		++n_ref;
	}
	
	RefSeq=ks->seq.s;
        //strncpy(dest, RefSeq + start, end - start+1);

  //char*       substr = substring(RefSeq, start, end); 

fprintf(fpnew,"seq: %d\n %d\n %s\n", start,end,RefSeq); 
	fprintf(stderr,"[%s] %d sequences, total length: %llu\n", __func__, n_ref, (long long)tot_len);
	kseq_destroy(ks);
	gzclose(fp_fa);
	fclose(fpnew);
return S;
}
*/



int main(int argc, char *argv[])  
{  
    gzFile fp;  
    kseq_t *seq;  
    int l;  
    int start, end;
    char *token,*token2,*token3, *token4;
    FILE *fpnew,*genomenew;

    if (argc - optind < 4) {  
        //prepare the genome file
	prepare_genome(argv[2],argv[3]);
	fprintf(stderr, "Usage: %s <in.fastq> <human_g1k_v37.fasta> <human_g1k_v37_modified.fasta> <output.fastq>\n", argv[0]);
	fprintf(stderr, "Example: ./wgsimref read1.fq ../human_g1k_v37.fasta human_new.fasta read3.fq\n");  
        return 1;  
    }  
    
    fpnew = fopen(argv[4], "w");
    genomenew = fopen(argv[3], "r");
    fp = gzopen(argv[1], "r"); // STEP 2: open the file handler  
    seq = kseq_init(fp); // STEP 3: initialize seq  
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence  
        fprintf(fpnew,"name: %s\n", seq->name.s);  
        if (seq->comment.l) 
		fprintf(fpnew,"comment: %s\n",seq->comment.s); 
	//split the name 
	token = strtok(seq->name.s, "_");
	token2 = strtok(NULL, "_");
	token3 = strtok(NULL, "_");
	start = atoi(token2);
	end = atoi(token3);
	token4 = read_genome(genomenew, fpnew, start, end);

	fprintf(fpnew,"%s\n%s\n%s\n%s\n", token,token2,token3,token4);
		 
	fprintf(fpnew,"seq: %s\n", seq->seq.s);
	fprintf(fpnew,"ref: %s\n", seq->seq.s);  
        if (seq->qual.l) 
		fprintf(fpnew,"qual: %s\n", seq->qual.s);  
	
    }  
    printf("return value: %d\n", l);  
    kseq_destroy(seq); // STEP 5: destroy seq  
    gzclose(fp); // STEP 6: close the file handler 
    fclose(fpnew);
    fclose(genomenew);
    return 0;  
}  
