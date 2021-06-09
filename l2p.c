
/*
 * todo  :

1) make sure custom works.
2) updategenes()
3) ensembl2gene()
4) pretty print  , pvals  if want raw data frame DOUBLE CHECK AGAINST EMAIL
5) Another column will be called “%Gene Hits per Pathway”, explained as  (A/(A+B))
6) bench marks permute vs. non . plot results.

vi l2p2.c ; make
gcc --Wall -Os -o l2p2 pwgenes.c small.c

example command line debug: cat deglist | valgrind -v --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes  ./l2p -permute

NUMBER_HITS … becomes Significant genes IN Pathway (A)
NUMBER_MISSES … becomes Non-Significant genes IN Pathway (B)
NUMBER_USER_GENES … becomes Significant genes NOT IN Pathway (C)
TOTAL_GENES_MINUS_INPUT … becomes Non-Significant Genes NOT IN Pathway (D)
     
      new_name                       old_name                   description
   1  pval                           pval                       fisher's 2x2 exact test
   2  fdr                            fdr                        false discovery rate: default is benjamini hochberg, GPCC method if permute=1 
   3  enrichment_score               ratio                      same as old but multiplied by 100 : ((number_hits /(number_hits+number_misses)) - (number_user_genes/(number_user_genes+total_gens_minus_input))) * 100
** 4  percent_gene_hits_per_pathway  NEW FIELD                  (number_hits/(number_hits+number_misses))
   5  number_hits                    pwhitcount                 number of genes hit in pathway
   6  number_misses                  pwnohitcount               pathway number of genes in the pathway
   7  number_user_genes              inputcount                 total count of user genes (user input)
   8  total_genes_minus_input        pwuniverseminuslist        total number of unique genes in all pathways
   9  pathway_id                     pathwayaccessionidentifier canonical accession ( if available, otherwise assigned by us )
  10  category                       category                   KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  *was "source"*
  11  pathway_name                   pathwayname                Name of pathway
  12  pathway_type                   pathwaytype                functional_set,pathway,structural_complex,custom
  13  genes                          genes_space_separated      HUGO genes from user that hit the pathway


The Fisher’s exact test on the permutation table is basically stand-alone, it doesn’t need a permutation test initially.  But it looks like I’m testing something different from what is usually tested; for my test the total number of genes in the pathway isn’t relevant.  Good, I want to do something different.  

The contingency table: So I can think clearly I’ll set this up the way I’m used to (from epidemiology, first column is no disease, second is disease; first row in no exposure, second is exposure).  For us “disease” means occurring in the pathway in question; “exposure” means being on the investigator’s list of genes. 

First without the correction:
A:  Count of genes in the universe but not on the list and not in the pathway
B:  Count of genes in the universe but not on the list, in the pathway
C:  Count of genes on the list but not in the pathway
D:  Count of genes on the list and in the pathway                                            
c d b a
a     4  pwhitcount           number of genes hit in pathway
b     5  pwnohitcount         pathway number of genes in the pathway
c     6  inputcount           total count of user genes (user input)
d     7  pwuniverseminuslist  unique genes in universe minus inputcount

Now with the correction:
B and D are unchanged.
For A and C:  Sum over the genes, the number of pathways (in the universe of pathways) each gene occurs in.
For example suppose there are 3 genes in C (must be a very short list of genes!) so in the uncorrected table the entry in C is 3.
The first gene occurs in 5 pathways, the second in 2, the third in 10.  The entry in C is 17.  
Same for A, but this will always be a sum over a large number of genes.

The sum of the entries in the table will now be larger, so the table will need to be corrected, but ignore this for now. 

original                               GPCC
a=hits                                  d
b=pathway-hits                          b
c=list-hits                             c
d=universe-pathway-list+hits            a

so ...
   UUUUUUUUUUUUUUUUUUUUUUUUUUUUU 
   U                           U 
   U universe   PPPPPPPPPPPP   U 
   U            P          P   U 
   U            P pathway  P   U 
   U            P          P   U 
   U       LLLLLLLLLLLLLL  P   U 
   U       L    P       L  P   U 
   U       L    P hits  L  P   U 
   U       L    P       L  P   U 
   U       L    PPPPPPPPLPPP   U 
   U       L            L      U 
   U       L  list      L      U 
   U       L            L      U 
   U       LLLLLLLLLLLLLL      U 
   U                           U 
   UUUUUUUUUUUUUUUUUUUUUUUUUUUUU 

   UUUUUUUUUUUUUUUUUUUUUUUUUUUUU    UUUUUUUUUUUUUUUUUUUUUUUUUUUUU    UUUUUUUUUUUUUUUUUUUUUUUUUUUUU 
   U                           U    U                           U    U                           U 
   U universe                  U    U            PPPPPPPPPPPP   U    U                           U 
   U                           U    U            P          P   U    U                           U 
   U                           U    U            P pathway  P   U    U                           U 
   U                           U    U            P          P   U    U                           U 
   U                           U    U            P          P   U    U       LLLLLLLLLLLLLL      U 
   U                           U    U            P          P   U    U       L            L      U 
   U                           U    U            P          P   U    U       L            L      U 
   U                           U    U            P          P   U    U       L            L      U 
   U                           U    U            PPPPPPPPPPPP   U    U       L            L      U 
   U                           U    U                           U    U       L            L      U 
   U                           U    U                           U    U       L  list      L      U 
   U                           U    U                           U    U       L            L      U 
   U                           U    U                           U    U       LLLLLLLLLLLLLL      U 
   U                           U    U                           U    U                           U 
   UUUUUUUUUUUUUUUUUUUUUUUUUUUUU    UUUUUUUUUUUUUUUUUUUUUUUUUUUUU    UUUUUUUUUUUUUUUUUUUUUUUUUUUUU 

A = u - p - h
B = u - l 
C = l - h 
D = h

dcba

see https://www.biostars.org/p/252937/ on how DAVID does it.

list=300 hits=3  pw=40 genome=30001
         usergenes     genome
ipw         3-1         40
notinpw     297       29960

a=hits
b=list-hits
c=list
d=u-list

problem #1

problem #2
   One sided test - count

               diseased healthy
exposed           DE       HE
notexposed        DN       HN
OR=(DE/HE)/(DN/HN)

a/b / c/d

*/
#define NELSON_C 1

// RICH:  for me running C program:
#if NELSON_C
#define DEGLISTNAME "brca_up_in_meta500.txt"
#define OUTPUTFILE "up_in_meta BRCA gpct.txt"
//#define DEGLISTNAME "deglistupinreactive.txt"
// RICH: will need to be passed paramter:
#endif
#define SCALE_UNIVERSE 1
#define NUM_PERMUTES 10000

// QVAL_PERMUTE still being tested, keep 0
#define NUM_QVAL_PERMUTES 0
#define SIG_P 0.05

#define NELSON 0
   // George Nelson's method

#ifdef L2P_USING_R
#include <R.h>
#include <Rdefines.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

#if 1
 // support fisher.c GPL code 
#include <inttypes.h>
#endif

#include "pathworks.h"
#ifdef WEBASSEMBLY
#include "small.h"
#else

#if L2PUSETHREADS
#include <pthread.h>
#endif

#include "big.h"
#endif

#define MAX_INGENES 40000

extern unsigned short int pwgenes[];
extern struct smallgenetype genes[];
extern struct pwtype pws[];
extern unsigned int numpws;
extern unsigned int numgenes;
extern unsigned int numpwgenes;

struct smallgenetype *by_egids; // a copy of "genes" ,but sorted by egid (entrez gene id )

// George's globals
struct tree_with_count *head = (struct tree_with_count *)0;
struct tree_with_count **aug_treepointers;
struct tree_with_count **treepointers;
struct used_path_type **augreverse = (struct used_path_type **)0;
struct used_path_type **sig_path_ptr = (struct used_path_type **)0;
unsigned int deg_count = 0;
unsigned int ingenes[MAX_INGENES];
unsigned int ingenecnt;
unsigned int *gene_path_cts;
unsigned int *deg_path_cts;
unsigned int real_genect = 0;
unsigned int aug_gene_ct = 0;
unsigned int *auguniverse;
unsigned int n_sig_paths = 0;
double mean_u_pathcounts;
#if NELSON_C
FILE *path_gene_countsfile;
FILE *gene_path_countsfile;
#endif



#if 0
 // debug routine
void dump_used_paths(struct used_path_type *u, int num, char *msg)
{
    struct used_path_type *uptr; // used path pointer
    int i,j;
    char *z;

fprintf(stderr,"in dump_used_paths() num=%d msg=[%s]\n",num,msg); fflush(stderr); 
    for (i=0 ; i<100 ; i++)   // for each line of custom file
    {
        uptr = (u+i);
        fprintf(stderr,"acc %d %s\n",i,uptr->acc);  fflush(stderr); 
        fprintf(stderr,"name %d %s\n",i,uptr->name); 
        fprintf(stderr,"numgenes %d %u\n",i,uptr->numgenes); 
        fprintf(stderr,"numfixed %d %u\n",i,uptr->numfixedgenes); 
        fprintf(stderr,"genehits %d %p\n",i,uptr->genehits); 
        fprintf(stderr,"egids %d %p\n",i,uptr->egids); 
        for (j=0;j<10;j++)
       // for (j=0;j<uptr->numfixedgenes;j++)
        {
            z = egid2hugo(*(uptr->egids+j));
            fprintf(stderr,"gene %d %d = %u = %s [%s]\n",i,j,*(uptr->egids+j),z,msg);
        }
        fflush(stderr);
    }
    return;
}
#endif

char *type2string(int type)
{
   static char functional_set[] = "functional_set";
   static char pathway  [] = "pathway";
   static char structural_complex  [] = "structural_complex";
   static char custom_string  [] = "custom";
   static char null_info  [] = "NA";

   if (type == type_functional_set) return &functional_set[0];
   else if (type == type_pathway) return &pathway[0];
   else if (type == type_structural_complex) return &structural_complex[0];
   else if (type == type_custom) return &custom_string[0];
   return &null_info[0];
}

int cmp_hugo(const void *a, const void *b)
{
    return strcmp(((struct smallgenetype *)a)->hugo, ((struct smallgenetype *)b)->hugo);
}

int cmp_by_egid(const void *a, const void *b) // compare entrez gene id
{
    if      ( ((struct smallgenetype *)a)->egid < ((struct smallgenetype *)b)->egid ) return -1;
    else if ( ((struct smallgenetype *)a)->egid > ((struct smallgenetype *)b)->egid ) return 1;
    return 0;
}

char *egid2hugo(int egid)
{
    struct smallgenetype glocal;
    struct smallgenetype *gptr;

    glocal.egid = egid;
    gptr = (struct smallgenetype *)bsearch(&glocal,by_egids,numgenes,sizeof(struct smallgenetype),cmp_by_egid);
    if (gptr) return gptr->hugo;
    else return (void *)0;
}

unsigned int hugo2egid(char *h)
{
    struct smallgenetype glocal;
    struct smallgenetype *gptr;

    glocal.hugo = h;
    gptr = (struct smallgenetype *)bsearch(&glocal,genes,numgenes,sizeof(struct smallgenetype),cmp_hugo);
    if (gptr) return gptr->egid;
    else return (unsigned int)UINT_MAX;
}

struct smallgenetype *hugo2geneptr(char *h)
{
    struct smallgenetype glocal;
    struct smallgenetype *gptr;

    glocal.hugo = h;
    gptr = (struct smallgenetype *)bsearch(&glocal,genes,numgenes,sizeof(struct smallgenetype),cmp_hugo);
    return gptr;
}

unsigned short int hugo2usedindex(char *h)
{
    struct smallgenetype glocal;
    struct smallgenetype *gptr;
    unsigned short int ret;
    int idx;

    glocal.hugo = h;
    gptr = (struct smallgenetype *)bsearch(&glocal,genes,numgenes,sizeof(struct smallgenetype),cmp_hugo);
    if (gptr)
    {
        idx = gptr - &genes[0];
        if (idx > 20000) {fprintf(stderr,"ERROR in hub2usedindex %s gptr=%p idx=%d\n",h,gptr,idx); fflush(NULL); exit(0); }
        ret = (unsigned short int)idx;
    }
    else
    {
        ret = (unsigned short int)USHRT_MAX;
    }
    return ret;
}

int cmp_ui(const void *a, const void *b)
{
    if      ( *(unsigned int *)a < *(unsigned int *)b ) return -1;
    else if ( *(unsigned int *)a > *(unsigned int *)b ) return 1;
    return 0;
}

#if 0
int cmp_usi(const void *a, const void *b)
{
    if      ( *(unsigned short int *)a < *(unsigned short int *)b ) return -1;
    else if ( *(unsigned short int *)a > *(unsigned short int *)b ) return 1;
    return 0;
}

int cmp_double(const void *a, const void *b)
{
    if      ( *(double *)a < *(double *)b ) return -1;
    else if ( *(double *)a > *(double *)b ) return 1;
    return 0;
}

 not used
int cmp_float(const void *a, const void *b)
{
    if      ( *(float *)a < *(float *)b ) return -1;
    else if ( *(float *)a > *(float *)b ) return 1;
    return 0;
}
#endif

#if 0
/*
ReservoirSample(S[1..n], R[1..k])
  R[1] := S[1]
  for i from 2 to k do
      j := randomInteger(1, i)  // inclusive range
      R[i] := R[j]
      R[j] := S[i]
  for i from k + 1 to n do
      j := randomInteger(1, i)  // inclusive range
      if (j <= k)
          R[j] := S[i]
*/

static inline void reservoir(unsigned int s[], int n, unsigned int r[],int k)
{
    int i,j;

    r[0] = s[0];
    for (i=1;i<k;i++)
    {
        j = rand() % i;
        r[i] = r[j];
        r[j] = s[i];
    }
    for (i=k+1;i<n;i++)
    {
        j = rand() % i;
        if (j <= k)
            r[j] = s[i];
    }
    return;
}
#endif

// -- for benjamini hochberg FDR ...
struct ordertype
{
    double val;
    int order;
};

static int cmp_ordertype_by_order(const void *a, const void *b)
{
    struct ordertype *aa;
    struct ordertype *bb;
    aa = (void *)a;
    bb = (void *)b;
    if      (aa->order >  bb->order) return 1;
    else if (aa->order <  bb->order) return -1;
    return 0;
}

int cmp_ordertype_by_val_REV(const void *a, const void *b)
{
    struct ordertype *aa;
    struct ordertype *bb;
    aa = (void *)a;
    bb = (void *)b;
    if      (aa->val <  bb->val) return 1;
    else if (aa->val >  bb->val) return -1;
    // if (aa>bb) return 1; else if (aa<bb) return -1;
    return 0;
}

static void benjaminihochberg(int n,double pvals[], double returnpvals[])
{
/*
 here's the code from R that I re-imagined
                   i <- lp:1L
                   o <- order(p, decreasing = TRUE)
                   ro <- order(o)
                   pmin(1, cummin( n / i * p[o] ))[ro]
*/
    int j,k;
    struct ordertype *i;
    struct ordertype *o;
    struct ordertype *po;
    struct ordertype *cummin;
//    struct ordertype *ro;
//    struct ordertype *intermed;

    i = (struct ordertype *)malloc((sizeof(struct ordertype))*n);
    for (k=n,j=0;j<n;j++,k--) (i+j)->order=k;

#define RDEBUG 0
#if RDEBUG
FILE *fp;
fp = fopen("test.pvals","w");
#endif
    o = (struct ordertype *)malloc((sizeof(struct ordertype))*n);
    for (j=0 ; j<n ; j++)
    {
#if RDEBUG
fprintf(fp,"%20.18f\n",pvals[j]);
#endif
        (o+j)->val=pvals[j];
        (o+j)->order=j+1;
    }
#if RDEBUG
fclose(fp);
#endif
    qsort(o,n,sizeof(struct ordertype),cmp_ordertype_by_val_REV);

#if 0
    ro = (struct ordertype *)malloc((sizeof(struct ordertype))*n);
    for (j=0;j<n;j++)
    {
        (ro+j)->val = (double)(o+j)->order;
        (ro+j)->order = j+1;
    }
    qsort(ro,n,sizeof(struct ordertype),cmp_ordertype_by_val);
#endif
    po = (struct ordertype *)malloc((sizeof(struct ordertype))*n);
    memset(po,0,(sizeof(struct ordertype))*n);
    for (j=0;j<n;j++)
    {
        (po+j)->val = (double)pvals[j];
        (po+j)->order = (o->order); // why the hell isn't this ro? what the what?
    }
    qsort(po,n,sizeof(struct ordertype),cmp_ordertype_by_val_REV); // == p[o]

    cummin = (struct ordertype *)malloc((sizeof(struct ordertype))*n); // holds n / i * po
    for (j=0;j<n;j++)
    {
        (cummin+j)->val = (double)n / (double)(i+j)->order * ((po+j)->val) ;
    }
                   // Rcode: pmin(1, cummin( n / i * p[o] ))[ro]           ******************
    for (j=1;j<n;j++)
    {
        if ((cummin+j)->val > (cummin+j-1)->val)
            (cummin+j)->val = (cummin+j-1)->val;
    }
    for (j=0;j<n;j++)
    {
        if ((cummin+j)->val > 1)
            (cummin+j)->val = 1;
        (cummin+j)->order = (o+j)->order ;
    }
    qsort(cummin,n,sizeof(struct ordertype),cmp_ordertype_by_order);
#if RDEBUG
FILE *fp2;
fp2 = fopen("test.fdrs","w");
#endif
    for (j=0;j<n;j++)
    {
        returnpvals[j] = (cummin+j)->val;
#if RDEBUG
fprintf(fp2,"%20.18f\n",returnpvals[j]);
#endif
    }
#if RDEBUG
fclose(fp2);
#endif
    if (i) free(i);
    if (o) free(o);
    if (po) free(po);
    if (cummin) free(cummin);

    return;
}



#if 0

// modified from  Darel Rex Finley https://alienryderflex.com/quicksort/ - public domain
static inline int quickSort(unsigned int *arr, int elements) 
{
#define  MAX_LEVELS  1000

  // int  piv, beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R ;
  unsigned int  piv;
  int  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R ;

  beg[0]=0; end[0]=elements;
  while (i>=0) 
  {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=arr[L]; if (i==MAX_LEVELS-1) return -1;
      while (L<R) 
      {
        while (arr[R]>=piv && L<R) 
            R--; 
        if (L<R) arr[L++]=arr[R];
        while (arr[L]<=piv && L<R) 
            L++; 
        if (L<R) arr[R--]=arr[L]; 
      }
      arr[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L; }
    else {
      i--; }
  }
  return 0; 
}
static inline void subsamp(unsigned int s[], int n, int k)
{
    int i,j,dncnt;
    unsigned int  ui;

    for (dncnt=n,i=0;i<k;i++,dncnt--)
    {
        j = (rand() % dncnt) + i;
        ui = s[i];
        s[i] = s[j];
        s[j] = ui;
    }
    return;
}
int permute(int idx,struct used_path_type *uptr,unsigned int real_universe_cnt,unsigned int *real_universe, int ingenecnt)
{
    unsigned int r[40002];
    double   p[NUM_PERMUTES];
    int i,j,k;
    double d;
    int tmphitcnt;
    unsigned int ui_uk,ui_ej;

    memcpy(r,real_universe,real_universe_cnt*sizeof(unsigned int));
    for (i=0;i<NUM_PERMUTES;i++)
    {
        subsamp(r, real_universe_cnt,ingenecnt); // inlined if optimization on 
        quickSort(r,ingenecnt); // inlined if optimization on 
        j = k = tmphitcnt = 0;
        while ((j<uptr->numfixedgenes) && (k < ingenecnt))
        {
            ui_ej = *(uptr->egids+j);
            ui_uk = *(r+k);
            if (ui_ej == ui_uk)
            {
                tmphitcnt++;
                k++;
                j++;
                continue;
            }
            else if (ui_ej < ui_uk) j++;
            else                    k++;
        }
        d = exact22((int)tmphitcnt,uptr->numfixedgenes-tmphitcnt,ingenecnt,real_universe_cnt-ingenecnt);
        p[i] = d;
    }
    qsort(p,NUM_PERMUTES,sizeof(double),cmp_double);
    d = 1.0;
    for (j=0;j<NUM_PERMUTES;j++) 
    {
        if (p[j] > uptr->pval) { d = p[j]; break; }
    }
    d = (double)j/(double)NUM_PERMUTES;
    uptr->fdr2 = d;
    return 0;
}

#endif



void test22(int a,int b,int c,int d,double pv)
{
    FILE *fp;
    char cmd[512];
    char s[512];
    char junks[512];
    char junk2[512];
    char p[512];
    double diff;
    
    fp = fopen("tmp.R","w");
    
    fprintf(fp,"data = matrix(c(%d,%d,%d,%d), nrow = 2)\n",a,b,c,d);
//    fprintf(fp,"fisher.test(data,alternative=\"two.sided\")\n");
    fprintf(fp,"fisher.test(data,alternative=\"greater\")\n");
//    fprintf(fp,"fisher.test(data,alternative=\"less\")\n");
//    fprintf(fp,"library(\"exact2x2\")\n");
//    fprintf(fp,"exact2x2(data,alternative=\"two.sided\")\n");
//    fprintf(fp,"exact2x2(data,alternative=\"greater\")\n");
//    fprintf(fp,"exact2x2(data,alternative=\"less\")\n");
    fclose(fp);

    strcpy(cmd,"R -q --vanilla < tmp.R");
    fp = popen(cmd,"r");
    while (fgets(s,sizeof(s),fp))
    {
        if (strstr(s,"p-value"))
        {
            sscanf(s,"%s %s %s\n",junks,junk2,p);
            diff = pv - atof(p);
            if (diff < 0) diff = 0 - diff;
// fprintf(stderr,"%d %d %d %d %f %s %f\n",a,b,c,d,pv,p,diff);
        }
    }
    fclose(fp);
    return;
}


int do_pvals_and_bh(unsigned int ingenecnt, struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt,int oneside)
{
    double pv;
    double pv2;
    double enrichment_score;
    unsigned int i;
    double *pvals;
    double *fdrs;
    unsigned int localhitcnt ;
    struct used_path_type *uptr; // used path pointer 
    int a2 , b2 , c2 , d2;


// fprintf(stderr,"debug in do_pvals_and_bh() ingenecnt=%d num_used_paths=%d real_universe_cnt=%d oneside=%d\n",ingenecnt,num_used_paths,real_universe_cnt,oneside); fflush(stderr); 

    pv = 1.0;
// fprintf(stderr,"in do_pvals_and_bh() ingenecnt=%d num_used_paths=%d real_universe_cnt=%d 2\n",ingenecnt,num_used_paths,real_universe_cnt); fflush(stderr); 
    pvals = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));
// fprintf(stderr,"in do_pvals_and_bh() 3 pvals=%p\n",pvals); fflush(stderr); 
    if (!pvals) { fprintf(stderr,"ERROR: no memory\n"); fflush(stderr); return -1; }
// fprintf(stderr,"in do_pvals_and_bh() before loop to %d (num_used_paths)\n",num_used_paths); fflush(stderr); 
    for (i=0 ; i<num_used_paths;i++)
    {
        pv = pv2 = 1.0; // initialize 
// fprintf(stderr,"in do_pvals_and_bh() path=%d of %d\n",i,num_used_paths); fflush(stderr); 
        uptr = (usedpaths+i);
        localhitcnt = uptr->hitcnt;
        if (localhitcnt)
        {
            if (localhitcnt > (unsigned int)uptr->numfixedgenes) 
            {
                 fprintf(stderr,"ERROR: more hits than genes (localhitcnt %d > %d for %s)\n", 
                         localhitcnt , uptr->numfixedgenes,uptr->acc);
                 fprintf(stderr,"details: index=%d , hitcnt=%u\n", i,localhitcnt); 
                 fflush(stderr);
                 return -2;
            }
        }
#if NELSON
/*
A = u - p - h
B = u - l 
C = l - h 
D = h
*/
int u,p,h,l,a,b,c,d;
u = real_universe_cnt;
p = uptr->numfixedgenes;
h = localhitcnt;
l = ingenecnt;
a = u - p - h;
b = u - l;
c = l - h;
d = h;
        if (oneside)
        {
            pv = fisher22_1sided((uint32_t)a, (uint32_t)b, (uint32_t)c, (uint32_t)d, (uint32_t)1, (uint32_t)1);
#if 0
 if ((pv < 0.0) || (pv > 1.0))
 {
 fprintf(stderr,"error pv=%f\n"pv);
 exit(0);
 }
#endif
        }
        else
        {
            pv = exact22(a,b,c,d);
        }
        enrichment_score = ( ((double)localhitcnt / (double)uptr->numfixedgenes) - ((double)ingenecnt/(double)real_universe_cnt) );
	enrichment_score = 100.0 * enrichment_score;
	// FIX fix 
	// 6.	Another column will be called “%Gene Hits per Pathway”, explained as  (A/(A+B))
        *(pvals+i) = uptr->pval = pv;
        uptr->a = a;
        uptr->b = b;
        uptr->c = c;
        uptr->d = d;
        uptr->enrichment_score = enrichment_score;
#endif // nelson 
/*
     4  pwhitcount            number of genes hit in pathway
     5  pwnohitcount          pathway    number of genes in the pathway
     6  inputlist                  total count of user genes (user input)
     7  pwuniverseminuslist    tinlistpathwaysuniverse    total number of unique genes in all pathways
a=hits
b=thispw-hits
c=list
d=u-list
*/
        a2 = localhitcnt;
        b2 = uptr->numfixedgenes - localhitcnt;
        c2 = ingenecnt;
        d2 = real_universe_cnt - ingenecnt;

        if (oneside)
        {
// yyy
#if 0
            pv = exact22_oneside((uint32_t)a, (uint32_t)b, (uint32_t)c, (uint32_t)d,0); // last arg is debug flag
 // testing 
double pv2 = exact22(a,b,c,d);
fprintf(stderr,"after exact22_oneside pv1 = %20.17f , pv2 = %20.17f %d %d %d %d\n",pv,pv2,a,b,c,d); fflush(NULL);
// test22(a,b,c,d,pv);
#endif
#if 1
/*
int u,p,h,l;
u = real_universe_cnt;
p = uptr->numfixedgenes;
h = localhitcnt;
l = ingenecnt;
int a2=hits
int b2=pathway - hits
int c2=list - hits
int d2=universe - pathway - list + hits
*/
a2 = localhitcnt;
b2 = uptr->numfixedgenes - localhitcnt;
c2 = ingenecnt-localhitcnt;
d2 = real_universe_cnt - ( uptr->numfixedgenes ) - ingenecnt + localhitcnt;
            pv2 = exact22_oneside((uint32_t)a2, (uint32_t)b2, (uint32_t)c2, (uint32_t)d2,0); // last arg is debug flag
// fprintf(stderr,"dbg pv = %f pv2=%f %d %d %d %d %d %d %d %d %s\n",pv,pv2,u,p,h,l,a2,b2,c2,d2,uptr->acc); 
            pv = pv2;
#endif
        }
        else
        {
            pv = exact22(a2,b2,c2,d2);
// fprintf(stderr,"called exact22(), pv=%f\n",pv);fflush(stderr); 
        }
        uptr->pval2 = *(pvals+i) = pv2;
        uptr->pval = *(pvals+i) = pv;
// fprintf(stderr,"rpf setting *(pvals+i) %d  %f\n",i,*(pvals+i) ); 
        uptr->a = a2;
        uptr->b = b2;
        uptr->c = c2;
        uptr->d = d2;
        // enrichment_score = ( ((double)localhitcnt / (double)uptr->numfixedgenes) - ((double)ingenecnt/(double)real_universe_cnt) );
        enrichment_score = ( ((double)a2 / (double)(a2+b2)) - ((double)c2/(double)(c2+d2)) );

        uptr->enrichment_score = enrichment_score;
    }
#if 0
for (i=0 ; i<num_used_paths ; i++)
{
 fprintf(stderr,"rpf checking 1 *(pvals+i) %d  %f\n",i,*(pvals+i) ); 
}
#endif

    fdrs = (double *)malloc((size_t)(sizeof (double)*num_used_paths)); 
    if (!fdrs) { free(pvals); /* clean up */ fprintf(stderr,"ERROR: no memory in do_pvals_and_bh() 2\n"); return -3; }

// fprintf(stderr,"rpf before benjaminihochberg\n"); fflush(stderr); 

#if 0
for (i=0 ; i<num_used_paths ; i++)
{
fprintf(stderr,"rpf checking 2 *(pvals+i) %d  %f\n",i,*(pvals+i) ); 
}
#endif
    benjaminihochberg(num_used_paths,pvals,fdrs);
// fprintf(stderr,"rpf after benjaminihochberg\n"); fflush(stderr); 
    for (i=0 ; i<num_used_paths ; i++)
    {
         usedpaths[i].fdr = *(fdrs+i);
// fprintf(stderr,"%d %f\n", i,usedpaths[i].fdr);
    }
    free(pvals);
    free(fdrs);
// fprintf(stderr,"in do_pvals_and_bh() end\n"); fflush(stderr); 
    return 0;
}


int do_just_bh(unsigned int ingenecnt, struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt)
{
    double enrichment_score;
    unsigned int i;
    double *pvals;
    double *fdrs;
    unsigned int localhitcnt ;
    struct used_path_type *uptr; // used path pointer 
    int a,b,c,d;


// fprintf(stderr,"debug in do_just_bh() ingenecnt=%d num_used_paths=%d real_universe_cnt=%d \n",ingenecnt,num_used_paths,real_universe_cnt); fflush(stderr); 

    pvals = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));
    if (!pvals) { fprintf(stderr,"ERROR: no memory\n"); fflush(stderr); return -1; }
    for (i=0 ; i<num_used_paths;i++)
    {
        uptr = (usedpaths+i);
        localhitcnt = uptr->hitcnt;
        if (localhitcnt)
        {
            if (localhitcnt > (unsigned int)uptr->numfixedgenes) 
            {
                 fprintf(stderr,"ERROR: more hits than genes (localhitcnt %d > %d for %s)\n", 
                         localhitcnt , uptr->numfixedgenes,uptr->acc);
                 fprintf(stderr,"details: index=%d , hitcnt=%u\n", i,localhitcnt); 
                 fflush(stderr);
                 return -2;
            }
        }
/*
     4  pwhitcount            number of genes hit in pathway
     5  pwnohitcount          pathway    number of genes in the pathway
     6  inputlist                  total count of user genes (user input)
     7  pwuniverseminuslist    tinlistpathwaysuniverse    total number of unique genes in all pathways
a=hits
b=thispw-hits
c=list
d=u-list
*/
        a = localhitcnt;
        b = uptr->numfixedgenes - localhitcnt;
        c = ingenecnt;
        d = real_universe_cnt - ingenecnt;

  /* new way : */
       *(pvals+i) = uptr->pval;
        uptr->a = a;
        uptr->b = b;
        uptr->c = c;
        uptr->d = d;
        // enrichment_score = ( ((double)localhitcnt / (double)uptr->numfixedgenes) - ((double)ingenecnt/(double)real_universe_cnt) );
        enrichment_score = ( ((double)a / (double)(a+b)) - ((double)c/(double)(c+d)) );
        uptr->enrichment_score = enrichment_score;
    }
    fdrs = (double *)malloc((size_t)(sizeof (double)*num_used_paths)); 
    if (!fdrs) { free(pvals); /* clean up */ fprintf(stderr,"ERROR: no memory in do_just_bh() 2\n"); return -3; }
#if 0
 // debug
fprintf(stderr,"before benjaminihochberg %d %p %p\n",num_used_paths,pvals,fdrs); fflush(NULL); 
for (i=0;i<num_used_paths;i++)
{
fprintf(stderr," pv %d %f\n",i,*(pvals+i)); 
}
#endif
    benjaminihochberg(num_used_paths,pvals,fdrs);
    for (i=0 ; i<num_used_paths ; i++)
    {
         usedpaths[i].fdr = *(fdrs+i);
    }
    free(pvals);
    free(fdrs);
    return 0;
}


#if L2PUSETHREADS
#define NUM_THREADS 4
int permute_range( unsigned int num_used_paths,struct used_path_type usedpaths[],unsigned int real_universe_cnt,unsigned int *real_universe, int ingenecnt, int lo, int hi)
{
    int i;

    for (i=lo ; i<=hi ; i++)
    {
// fprintf(stderr, "%d from %d to %d of %d %s \n",i,lo,hi,num_used_paths,usedpaths[i].name); fflush(stderr); 
        permute(i,&usedpaths[i],real_universe_cnt,real_universe,ingenecnt);
    }
    return 0;
}

struct thread_data_t {
    int tid;
    int myid;
    int lo;
    int hi;
    unsigned int ingenecnt;
    struct used_path_type *usedpaths;
    unsigned int num_used_paths;
    unsigned int real_universe_cnt;
    unsigned int *real_universe;
    unsigned int *sampling;
    float *pspace;
};

#if 0
/* thread function */
void *thr_func(void *arg) 
{
    struct thread_data_t *data = (struct thread_data_t *)arg;
    unsigned int *r;
   // double   p[NUM_PERMUTES];
    int i,j,k;
    double d;
    int tmphitcnt;
    unsigned int ui_uk,ui_ej;
    // float *pspace = data->pspace;
    unsigned int ingenecnt = data->ingenecnt;
    struct used_path_type *uptr = data->usedpaths;

    r = data->sampling;
    for (i=data->lo;i<data->hi;i++)
    {
        j = k = tmphitcnt = 0;
        while ((j<uptr->numfixedgenes) && (k < ingenecnt))
        {
            ui_ej = *(uptr->egids+j);
            ui_uk = *(r+k);
            if (ui_ej == ui_uk)
            {
                tmphitcnt++;
                k++;
                j++;
                continue;
            }
            else if (ui_ej < ui_uk) j++;
            else                    k++;
        }

        d = exact22((int)tmphitcnt,uptr->numfixedgenes,ingenecnt,       data->real_universe_cnt);
	//fix here
        // *(pspace + (i*NUM_PERMUTES)+data->used_path_index; = d;
    }
#if 0
    qsort(p,NUM_PERMUTES,sizeof(double),cmp_double);
    d = 1.0;
    for (j=0;j<NUM_PERMUTES;j++) 
    {
        if (p[j] > uptr->pval) { d = p[j]; break; }
    }
    d = (double)j/(double)NUM_PERMUTES;
    uptr->fdr2 = d;
#endif
  pthread_exit(NULL);
}
#endif


#if 0
int oldthreaded_permutes(unsigned int ingenecnt, struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt,unsigned int *real_universe)
{
    pthread_t thr[NUM_THREADS];
    unsigned int r[40002];
    struct thread_data_t thr_data[NUM_THREADS];
    int i,j,k;
    int rc;
    unsigned int *pspace;

    pspace = malloc(sizeof(float)*NUM_PERMUTES*num_used_paths); // no dire need to memset
    memcpy(r,real_universe,real_universe_cnt*sizeof(unsigned int));
fprintf(stderr,"in oldthreaded_permutes()\n"); fflush(stderr); 
    for (i=0;i<NUM_PERMUTES;i++)
    {
        subsamp(r, real_universe_cnt,ingenecnt); // inlined if optimization on 
        quickSort(r,ingenecnt); // inlined if optimization on 
        k = NUM_PERMUTES/NUM_THREADS;
        for (j = 0; j < NUM_THREADS; ++j) 
        {     /* create threads */
            thr_data[j].lo = j*k;
            thr_data[j].hi = (j*k)+k-1;
            if (j == (NUM_THREADS-1))
                thr_data[j].hi = num_used_paths-1; // note "-1", 
fprintf(stderr,"%d of %d [ %d .. %d ] of %d size= %d\n",j,NUM_THREADS,thr_data[j].lo,thr_data[j].hi,num_used_paths,thr_data[j].hi-thr_data[j].lo); 
            thr_data[j].myid = j;
            thr_data[j].ingenecnt = ingenecnt;
            thr_data[j].sampling = r;
            // FIX thr_data[j].usedpathsindex = usedpaths;
            thr_data[j].usedpaths = usedpaths;
            thr_data[j].num_used_paths = num_used_paths;
fprintf(stderr,"creating thread %d\n",j); fflush(stderr); 
            if ((rc = pthread_create(&thr[j], NULL, thr_func, &thr_data[j]))) 
            {
              fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
              return EXIT_FAILURE;
            }
        }
      /* block until all threads complete */
        for (j = 0; j < NUM_THREADS; ++j) 
        {
fprintf(stderr,"blocking %d \n",j);
            pthread_join(thr[j], NULL);
        }
fprintf(stderr,"after blocking  \n");
    }
  return EXIT_SUCCESS;
}
#endif

#endif

inline void shuffle(unsigned int s[], int n)
{
    int i,j;
    unsigned int ui;

// for i from 0 to n−2 do
//      j ← random integer such that i ≤ j < n
//      exchange a[i] and a[j]
// or
// for i from n−1 downto 1 do
//      j ← random integer such that 0 ≤ j ≤ i
//      exchange a[j] and a[i]

    for (i=n-1;i>=1;i--)
    {
#ifdef L2P_USING_R
// must used GetRNGstate() before calls to unif_rand AND ut PutRNGstate after calls 
        j = (int)(unif_rand() * (double)i); // unif_rand() appears to return between 0.0 and 1.0
#else
        j = rand() % i;
#endif
        ui = s[i];
        s[i] = s[j];
        s[j] = ui;
    }
    return;
}


struct tree_type
{
    unsigned int val; // entrez gene id
    struct tree_type *left;
    struct tree_type *right;
};

void free_tree(struct tree_type *node)
{
    if (node==(void *)0) return;
    if (node->left != (void *)0)  free_tree(node->left);
    if (node->right != (void *)0) free_tree(node->right);
    free(node); 
    return;
}

void put_tree_to_array(struct tree_type *node, unsigned int puthere[], int *index_ptr)
{
    if (node != (void *)0)
    {
        put_tree_to_array(node->left,puthere,index_ptr);
        puthere[*index_ptr] = node->val;
        *index_ptr = *index_ptr + 1;
        put_tree_to_array(node->right,puthere,index_ptr);
    }
    return;
}
void free_tree_with_count(struct tree_with_count *node)
{
    if (node)
    {
        if (node->all_gene_paths) free(node->all_gene_paths);
        if (node->left) free_tree_with_count(node->left); // recursive calls
        if (node->right)free_tree_with_count(node->right);
        free(node);
    }
    return;
}
int put_to_tree(unsigned int this_gene, struct tree_with_count *head) // tree below head keeps growing; we will keep using it
{
    // has to:
    // 1. create a node (after the first, head) for each new gene (fill in gene number to node)
    //      --in left-right order
    // 2. cound the number of unique genes (one for each new node, n.b. first node created above)
    // 3. add how many instances there are of each gene to ->count
     struct tree_with_count *z, *n;
     z = head;
     //printf("this_gene = %d\n", this_gene);
    while (1) // follow the tree down to the bottom switching right or left by gene number. New node if needed. Then outer loop takes to next gene
     {
          if (z->val == this_gene)
          {
              z->count += 1;
              //if (z->count == 1) ++real_genect; // shouldn't happen, we've seen this node before // this wastes time...
              //head->all_gene_paths = malloc(sizeof(struct used_path_type));// not yet
              //head->all_gene_paths = *(augreverse+ui);
              return(1); //gene has been placed, return
          }
          else if (this_gene < z->val)
          {
              if (z->left == (void *)0)
              {
                  n = z->left = malloc(sizeof(struct tree_with_count));
                  n->val=this_gene;
                  //n->all_gene_paths = &this_path; ? have to run genes and malloc to # repeats?
                  ++real_genect; // OK, it's a new node
                  n->count = 1;
                  n->deg = 0;
                  n->left=(void *)0;
                  n->right=(void *)0;
                  n->all_gene_paths=(void *)0;
                  n->pathindex = 0;
                  return(1);//gene has been placed, return
              }
              z = z->left; // and back to top of while
          }
          else  /// i.e.  (z->val > this_gene)
          {
              if (z->right == (void *)0)
              {
                  n = z->right = malloc(sizeof(struct tree_with_count));
                  n->val=this_gene;
                  ++real_genect;
                  n->count = 1;
                  n->deg = 0;
                  n->left=(void *)0;
                  n->right=(void *)0;
                  n->all_gene_paths=(void *)0;
                  n->pathindex = 0;
                  return(1);//gene has been placed, return
              }
              z = z->right; // and back to top of while
          }
     }
}
// RICH: Beginning of my code
void malloc_pathpointers(struct tree_with_count *node) // counts aligned with universe (real_universe)
{
    // int i;
// rpf     static int degindex = 0;
    if (node)
    {
        malloc_pathpointers(node->left);
        node->all_gene_paths = malloc(sizeof(struct used_path_type *)*(node->count));
        if (!node->all_gene_paths)
        {
            fprintf(stderr,"ERROR: not enough memory in malloc_pathpointers()\n"); fflush(stderr);
        }
        malloc_pathpointers(node->right);
    }
    return;
}
int path_pointers_to_tree(unsigned int this_gene, struct used_path_type *this_path, struct tree_with_count *head)
// following put_to_tree, but tree is already there
// I think right, have to find the unsorted genes here, this is guide for filling the pointer
// could have been done while the tree was created, but don't have space for the pointers malloc'd yet
{
    // all nodes exist,
    // has to:
     struct tree_with_count *z;
     // struct tree_with_count *n;
     z = head;
     while (1) // follow the tree down to the bottom switching right or left by gene number. Tree created already, fails if new node needed.
     {
          if (z->val == this_gene)
          {
              z->all_gene_paths[z->pathindex] = this_path;
              ++(z->pathindex);
              if (z->pathindex > z->count) // note this catches trouble after it's started (if z->pathindex > z->count, overwrote z->all_gene_paths above).  Could put above but wastes time.
              {
                  printf("bad pathway count for gene %ds, exiting\n", this_gene);
                  exit(0);
              }
              return(1);
          }
          else if (this_gene < z->val)
          {
              if (z->left == (void *)0)
              {
                  printf("gene %d not found while sorting pathpointers, exiting\n", this_gene);
                  exit(0);
              }
              z = z->left; // and back to top of while
          }
          else  // i.e.  (z->val > this_gene)
          {
            if (z->right == (void *)0)
            {
                  printf("gene %d not found while sorting pathpointers, exiting\n", this_gene);
            }
            z = z->right; // and back to top of while
          }
     }
}

static inline void put_count_tree_to_array(struct tree_with_count *node, unsigned int universe[], unsigned int *index_ptr, unsigned int *aug_ptr) // counts aligned with universe (real_universe)
{
    // no trap for going off the tree
    unsigned int ui;
    if (node != (void *)0) //  trap against going off the tree, but are we catching genes with no match?
    {
        //printf("left\t");
        put_count_tree_to_array(node->left,universe,index_ptr,aug_ptr);
        universe[*index_ptr] = node->val;
        // printf("tst tree %d\n",node->val);
        gene_path_cts[*index_ptr] = node->count;
        ++*index_ptr;
#if NELSON_C
        // debug
        fprintf(gene_path_countsfile, "%d\t%d\n", node->val, node->count);
#endif
        for (ui = 0; ui < node->count; ui++)
        {
            auguniverse[*aug_ptr] = node->val; // to randomly select genes for permutation test, weighted by pathway count
            aug_treepointers[*aug_ptr] = treepointers[*index_ptr] = node; // multiple pointers to same node,
            //printf("%d %d\n", *aug_ptr, node->val);
           //augreverse[*aug_ptr] = node; // replaced by
            *aug_ptr = *aug_ptr + 1;
        }
        //printf("right\t");
        put_count_tree_to_array(node->right,universe,index_ptr,aug_ptr);
    }
    //printf("up\t");
    return;
}
static inline void put_deg_tree_to_array(struct tree_with_count *node, unsigned int *index_ptr) // counts aligned with universe (real_universe)
{
    // no trap for going off the tree
    if (node != (void *)0) //  trap against going off the tree, but are we catching genes with no match?
    {
        put_deg_tree_to_array(node->left,index_ptr);
        if (node->deg)
        {
            deg_path_cts[*index_ptr] = node->count;
            ingenes[*index_ptr] = node->val;
            ++*index_ptr;
        }
        put_deg_tree_to_array(node->right,index_ptr);
    }
    return;
}
int degs_to_tree(unsigned int this_gene, struct tree_with_count *head) // tree below head keeps growing; we will keep using it
{
    int i;
    struct tree_with_count *z, *node;
    struct used_path_type *thispath;
     z = head;
     while (1) // follow the tree down to the bottom switching right or left by gene number. If hit void, gene not in universe. Then outer loop takes to next gene
     {
          if (z->val == this_gene)
          {
              z->deg = 1;
              node = z;
              ++deg_count;
              for (i=0; i<node->count; i++)
              {
                  thispath = node->all_gene_paths[i];
                  thispath->gpcountsum += node->count; // this sums the gene gp cts for all hits in the pathway
              }
              return(1);
          }
          else if (this_gene < z->val)
          {
              if (z->left == (void *)0)
              {
                  printf("deg %d not in universe\n", this_gene);
                  return(0); // deg not in universe
              }
              z = z->left; // and back to top of while
          }
          else  /// i.e.  (z->val > this_gene)
          {
              if (z->right == (void *)0)
              {
                  printf("deg %d not in universe\n", this_gene);
                  return(0); // deg not in universe
              }
              z = z->right; // and back to top of while
          } // but if coded wrong could go on forever, no trap; but ok, keeps going down till a null
     }
}

static void chooseAugRandList(unsigned int n, unsigned int m, unsigned int seed, unsigned int *randlist) // n is size of augmented gene universe, m number to choose
{
    unsigned int i,j;
    unsigned int randint;
    unsigned int randgene;
    double randrslt;
    double scale = (double) (n-1)/(RAND_MAX);
#ifdef L2P_USING_R
       // method to shut up compiler warnings
    #define UNUSED(x) (void)(x)
    UNUSED(seed);
        // user can set the seed, by calling a function like this  set.seed(12345) , before calling
#else
    static int inited_random_flag = 0;
    if (inited_random_flag == 0)
    {      // should go elsewhere; where do we set for test run (fixed seed)?
           // for non R, the seed gets passed in as a paramter.  0 means let this code decide (randome based on time and process id "PID" )
        if (seed != 0) srand(seed);
        else srand(time(NULL) + getpid());
        inited_random_flag = 1;
    }
 // PutRNGstate(); // restore random state
#endif
    for (i=0;i<m;)
    {
#ifdef L2P_USING_R
        randrslt = (int)(unif_rand() * (double)i)+1;   // unif_rand() appears to return between 0.0 and 1.0
#else
        randrslt = rand();
#endif
        randint = randrslt*scale+1;
        randgene = auguniverse[randint];
        j = 0;
        for (j=0;j<i;j++)
        {
            if (randgene == randlist[j]) break;
        }
        if (j == i)
        {
            randlist[i] = randint;
            //fprintf(stdout, "i, j, chosen index, gene %f %d %d %d %d \n",randrslt, i, j, randint, randgene);
            ++i;
        }
    }
}


void chooseRandList(unsigned int n, unsigned int m, int seed, unsigned int *randlist) // n is size of augmented gene universe, m number to choose
{
    unsigned int i;
    double randrslt;
    double scale = (double) (n-1)/(RAND_MAX);
    // unsigned int j;
    // unsigned int randint;
    // unsigned int randgene;
    //double scale = (double) n/RAND_MAX;
#ifdef L2P_USING_R
UNUSED(seed);
    //      GetRNGstate();
#else
    static int inited_random_flag = 0;
    if (inited_random_flag == 0)  // should go elsewhere; where do we set for test run (fixed seed)?
    {
        if (seed != 0) srand(seed);
        else srand(time(NULL) + getpid());
        inited_random_flag = 1;
    }
#endif
    for (i=0 ; i<m ; i++)
    {
#ifdef L2P_USING_R
        randrslt = (int)(unif_rand());
#else
        randrslt = rand();
#endif
        randlist[i] = randrslt*scale+1;
    }
}/**/


// rpf int tablecalc(struct used_path_type usedpaths[], unsigned int num_used_paths, unsigned int *real_universe, unsigned int real_universe_cnt,int oneside)
int gpcount_correction(struct used_path_type usedpaths[], unsigned int num_used_paths) // ,int oneside)
{
    int i;
    double mean_gpcts;
    struct used_path_type *uptr; // used path pointer
    for (i=0 ; i<num_used_paths;i++)
    {
        uptr = (usedpaths+i);
        mean_gpcts =uptr->gpcountsum/uptr->hitcnt;
        //if (uptr->hitcnt == 0) uptr->pval4 = 99;
        //else
            uptr->pval4 = pow(mean_gpcts/mean_u_pathcounts, log(uptr->hitcnt))*uptr->pval;
    }
    return 0;
}

int tablecalc(struct used_path_type usedpaths[], unsigned int num_used_paths, unsigned int real_universe_cnt) // ,int oneside)
{
    double enrichment_score, OR, pv;
    unsigned int i, I=0;
    double *pvals;
    double *fdrs;
    unsigned int localhitcnt;
    struct used_path_type *uptr; // used path pointer
    int a,b,c,d;
    
    //fprintf(stderr,"debug in tablecalc() ingenecnt=%d num_used_paths=%d real_universe_cnt=%d \n",ingenecnt,num_used_paths,real_universe_cnt); fflush(stderr);
    
    pvals = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));  // all in pathway struct, don't need
    if (!pvals) { fprintf(stderr,"ERROR: no memory\n"); fflush(stderr); return -1; }
    memset(pvals,0,sizeof (double)*(num_used_paths));
    n_sig_paths = 0;
    for (i=0 ; i<num_used_paths;i++)
    {
        uptr = (usedpaths+i);
        localhitcnt = uptr->hitcnt;
        //if (!localhitcnt) continue; //  Let's only look at overexpression
        if (localhitcnt > (unsigned int)uptr->numfixedgenes)
        {
             fprintf(stderr,"ERROR: more hits than genes (localhitcnt %d > %d for %s)\n",
                     localhitcnt , uptr->numfixedgenes,uptr->acc);
             fprintf(stderr,"details: index=%d , hitcnt=%u\n", i,localhitcnt);
             fflush(stderr);
             return -2;
        }
/*
     4  pwhitcount            number of genes hit in pathway
     5  pwnohitcount          pathway    number of genes in the pathway
     6  inputlist                  total count of user genes (user input)
     7  pwuniverseminuslist    tinlistpathwaysuniverse    total number of unique genes in all pathways
    a=hits
    b=thispw-hits
    c=list
    d=u-list
*/
        /*a = localhitcnt;
        b = uptr->numfixedgenes - localhitcnt;
        c = ingenecnt;
        d = real_universe_cnt - ingenecnt;*/
        
        d = localhitcnt;
        c = ingenecnt - localhitcnt;
        b = uptr->numfixedgenes - localhitcnt;
        a = real_universe_cnt - uptr->numfixedgenes - ingenecnt + localhitcnt;

        /* new way : */
        //j *(pvals+i) = uptr->pval;
        uptr->a = a;
        uptr->b = b;
        uptr->c = c;
        uptr->d = d;
        // enrichment_score = ( ((double)localhitcnt / (double)uptr->numfixedgenes) - ((double)ingenecnt/(double)real_universe_cnt) );
        enrichment_score = ( ((double)a/(double)b) - ((double)c/(double)d) );
        uptr->enrichment_score = enrichment_score;
        if(c*b != 0)
        {
            OR = ((double)a*(double)d)/((double)c*(double)b);
            uptr->OR = OR;
#if NELSON_C
            if (OR > 1000)
            {
                fprintf(stderr,"Note:big OR %f, path %s\n", OR, uptr->acc);
            }
#endif
        }
        else uptr->OR = OR = 99999;
        pv = exact22(a,b,c,d);
        //pv = exact22_oneside((uint32_t)a, (uint32_t)b, (uint32_t)c, (uint32_t)d,0); // last arg is debug flag
        /*if (!strncmp("basal lamina", uptr->name, 30))
        {
            printf("%s", uptr->name);
        }*/
        uptr->pval = pv;
        if ((pv < SIG_P) & (OR > 1)) ++n_sig_paths;

    }
    sig_path_ptr = malloc((sizeof(struct used_path_type *)*n_sig_paths));
    // HERE WE WOULD WANT TO SORT THE SIGNIFICANT P VALUES, FOR SPEED LATER
//    struct used_path_type *tstptr;
    I = 0;
    for (i=0; i<num_used_paths; i++)
    {
        uptr = (usedpaths+i);
        if((uptr->pval < SIG_P) & (uptr->OR > 1))
        {
            *(sig_path_ptr+I) = uptr; // since unlabelled, essential that this array of pointers has the same order as auguniverse
//           tstptr = *(sig_path_ptr+I);
            //printf("%11.2e\n", uptr->pval);
            uptr->rand_paths_more_sig = 0;
            I++;
            //printf("%11.2e %11.2e %s\n", tstptr->pval, tstptr->pval3, tstptr->name);
            //printf("%11.2e %11.2e %s\n", (*(sig_path_ptr+I))->pval, (*(sig_path_ptr+I))->pval3, (*(sig_path_ptr+I))->name);
        }
        if(I>n_sig_paths)
        {
            printf("trouble in setting up sig paths\n");  //CLEAN UP
            exit(1);
        }
    }
    //  I HAVEN'T LOOKED AT THE FDR CODE
    fdrs = (double *)malloc((size_t)(sizeof (double)*num_used_paths));
    if (!fdrs) { free(pvals); /* clean up */ fprintf(stderr,"ERROR: no memory in do_just_bh() 2\n"); return -3; }
    #if 0
     // debug
    fprintf(stderr,"before benjaminihochberg %d %p %p\n",num_used_paths,pvals,fdrs); fflush(NULL);
    for (i=0;i<num_used_paths;i++)
    {
    fprintf(stderr," pv %d %f\n",i,*(pvals+i));
    }
    #endif
    benjaminihochberg(num_used_paths,pvals,fdrs);
    for (i=0 ; i<num_used_paths ; i++)
    {
         usedpaths[i].fdr = *(fdrs+i);
    }
    free(pvals);
    free(fdrs);
    return 0;
}

#if 0
unsigned int GPCC(struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt, unsigned int *real_universe)
#endif
unsigned int GPCC(struct used_path_type usedpaths[], unsigned int num_used_paths, unsigned int real_universe_cnt, unsigned int *real_universe, int seed)
{
    /*
     Create a vector of gene names, with each name appearing as often as the number of pathways the gene appears in.
     From this list randomly select as many genes as on your DEG list.
     Replace repeats, need a list of unique genes
     ii.     Create the contingency table for each pathway of interest (list of top pathways)
     iii.     Obtain and record its p value
     iv.     For each pathway of interest keep count of how many of the 200 permutations produced more significant p values than the real data
     
     Repeat a) etc. 200 times
     The counts/200 are approximate p values for the pathways.
     */
    unsigned int i,I,j,k,K;
    int utilctr;
    unsigned int this_gene;
    unsigned int ui = 0;
    unsigned int d_up = 0, d_down = 0, d_equal = 0;
    unsigned int a, rand_d, rand_c, rand_b;
    int var_d;
    double randOR;
    unsigned int sig_rand_p_ct = 0;
    unsigned int aug_gene_ct = 0;
    unsigned int *randgenesindex = (unsigned int *)0;
    unsigned int *pool = (unsigned int *)0;
    double pv;
    struct tree_with_count *head = (struct tree_with_count *)0;
    struct used_path_type **augreverse = (struct used_path_type **)0;
    struct used_path_type  *this_path = (struct used_path_type *)0;
    unsigned int *index_ptr = (unsigned int *)0;
    unsigned int *aug_ptr = (unsigned int *)0;
    struct used_path_type *uptr = (struct used_path_type *)0;    // used path pointer
    struct used_path_type *thispath;
    struct tree_with_count *thisnode = (struct tree_with_count *)0;
    double *pspace = (double *)0;   // array of pvalues[num_used_paths*NUM_PERMUTES]
    FILE *deg_path_count_file = (FILE *)0;
    FILE *predict_fp = (FILE *)0;
    FILE *pathgenes = (FILE *)0;


    //int *pool = (unsigned int *)0;
    // Set up augmented universe, count_tree, and pointers to pathways
#if NELSON_C
    gene_path_countsfile = fopen("gene_pathway counts.txt", "w");
    deg_path_count_file = fopen("deg pathway counts.txt", "w");
    predict_fp = fopen("tst predictor of p vals.txt", "w");
#endif
    //path_gene_countsfile = fopen("pathway_gene counts.txt", "w");
    if (!usedpaths)
    {
        fprintf(stderr,"ERROR: null used_path_type passed to get_used_universe()\n"); fflush(stderr);
        return (unsigned int)0 ;
    }
    // get count of all genes with repititions in pathways
    #if 0
        pathgenes = fopen("paths_genes_names.txt", "w");
    #endif
    for (i=aug_gene_ct=0 ; i<num_used_paths ; i++)
    {         // get number of genes in all pathways, count duplicated gene names
        uptr = (usedpaths+i);
        aug_gene_ct += uptr->numfixedgenes;
        //fprintf(path_gene_countsfile, "%s \t %d\n", uptr->name, uptr->numfixedgenes);
        /*for (j=0;j<uptr->numfixedgenes;j++)
        {
            fprintf(pathgenes, "%s\t %s\t %d\n", uptr->name, uptr->acc, uptr->egids[j]);
        }*/
    }
        if (pathgenes) fflush(pathgenes);
    #if 1
    // fix fix fix  rpf
        auguniverse = malloc((aug_gene_ct+1) * sizeof(unsigned int)); // seems to be short by one uint , so just add it
    #else
        auguniverse = malloc(aug_gene_ct * sizeof(unsigned int));
    #endif
    mean_u_pathcounts = aug_gene_ct/(float) real_universe_cnt;
    if (!auguniverse)
    {
        fprintf(stderr,"ERROR: not enough memory in get_used_universe()\n"); fflush(stderr);
        return (unsigned int)0 ;
    }
    // for fast permute:  I need pointers from each gene (repeated) to all the pathways containing the gene
    // need to be sorted in gene number order
    // XXXXXXXXXX
    // NOTE MANY VARIANTS ON PUTTING TO TREE AND GETTING FROM TREE FOLLOW HERE
    augreverse = malloc(sizeof(struct used_path_type *)*(aug_gene_ct)); // temp, sorted contents go to ordrd_pathpointer
    for (i=k=0;i<num_used_paths;i++)
    {
        uptr = (usedpaths+i);
        for (j=0 ; j < uptr->numfixedgenes ; j++)
        {
            *(auguniverse+k) = *(uptr->egids + j);
            *(augreverse+k) = uptr; // since unlabelled, essential that this array of pointers has the same order as auguniverse
            k++;
            if (k > aug_gene_ct)
            {
                fprintf(stderr,"in GPCC ERROR i=%d k=%d\n",i,k);  fflush(stderr);
                exit(0);
            }
        }
    }
    head = malloc(sizeof(struct tree_with_count));
    head->val = *(auguniverse); // entrez gene id
    real_genect = 1;
    head->count = 1; // each node starts with count 1
    head->deg = 0;
    head->left=(void *)0;
    head->right=(void *)0;
    // run first in get used universe, or separate function
    // this gets counts
    // malloc to counts
    head->all_gene_paths = (void *)0;
    head->pathindex = 0;
    for (ui=1 ; ui<aug_gene_ct ; ui++)
     {
         this_gene = *(auguniverse+ui);
         this_path = *(augreverse+ui);
         put_to_tree(this_gene, head);
      }
    
    index_ptr = malloc(sizeof(unsigned int));
    memset(index_ptr,0,sizeof(unsigned int));
    aug_ptr = malloc(sizeof(unsigned int)*2); // wtf?
    memset(aug_ptr,0,sizeof(unsigned int)*2);
    aug_treepointers = malloc(sizeof(struct tree_with_count **)*(aug_gene_ct+1)); // valgrind complains
    memset(aug_treepointers ,0,sizeof(struct tree_with_count **)*(aug_gene_ct+1)); // valgrind complains
    // rpf treepointers = malloc(sizeof(struct tree_with_count **)*real_universe_cnt);
    treepointers = malloc(sizeof(struct tree_with_count **)*(real_universe_cnt+1));
    memset(treepointers ,0,sizeof(struct tree_with_count **)*(real_universe_cnt+1));
    gene_path_cts = malloc(sizeof(unsigned int)*(real_genect+1)); // valgrind complains so add 1
    memset(gene_path_cts ,0,sizeof(unsigned int)*(real_genect+1));

     malloc_pathpointers(head);
     // fprintf(stderr,"test after malloc_pathpointers\n"); fflush(NULL);
    utilctr = 0;
    for (ui=0 ; ui<aug_gene_ct ; ui++)
    {
         this_gene = *(auguniverse+ui);
         this_path = *(augreverse+ui);
         path_pointers_to_tree(this_gene, this_path, head); // recall path_pointers are needed for permutation test
        if(this_gene == 7040) ++utilctr;
    }
    *index_ptr = 0; // fcn passes down pointer
    *aug_ptr = 0; // fcn passes down pointer
    put_count_tree_to_array(head,real_universe,index_ptr, aug_ptr); // n.b. auguniverse is now a global variable--sorted by this recursive function
    deg_count = 0;
    for(i=0;i<ingenecnt;i++)
    {
        this_gene = ingenes[i];
        degs_to_tree(this_gene, head);
    }
// xxx
    ingenecnt = deg_count; // should just use deg_count everywhere
    deg_path_cts = malloc(sizeof(unsigned int)*deg_count);
    *index_ptr = 0;
// fprintf(stderr,"in GPCC 2.0 head=%p index_ptr=%p\n",head,index_ptr); fflush(NULL);
    put_deg_tree_to_array(head,index_ptr); // not getting deg count here, should I?
// fprintf(stderr,"in GPCC 2.1 deg_path_count_file=%p \n",deg_path_count_file); fflush(NULL);
    for (i=0;i<deg_count;i++) 
    {
         if (deg_path_count_file) fprintf(deg_path_count_file, "%d\n  ", deg_path_cts[i]);
    }
    
    // now I will select randomly from auguniverse
    // this produces a set of gene numbers
    // For each gene number, in I have an array of pointers to the pathways--
    // --in the pathway count list
    // rows will correspond to pathways in ptr_from_unique. Need a consistent list of pathways (1-n); the row refers to this list
    randgenesindex = malloc(sizeof(unsigned int *)*deg_count);
    //thispath = malloc(sizeof(struct used_path_type *));
//     int oneside = 1;
    // tablecalc(usedpaths, num_used_paths, real_universe, real_universe_cnt, oneside);
    tablecalc(usedpaths, num_used_paths,  real_universe_cnt);
    // FET has been done in tablecalc? Now do correction
    gpcount_correction(usedpaths, num_used_paths);
    //*setseed = 12345; // for testing, but get an error; never malloc'd integer storage
    for (i = 0; i < NUM_PERMUTES; i++)
    {
        // test: at i = 0, put in the actual genes
        // Choose ingenecnt unique genes
//        if ((i % 500) == 0) { fprintf(stderr,"permute i = %d\n", i); fflush(NULL); }
        for (I=0 ; I<num_used_paths ; I++) // test, remove
        {
            uptr = (usedpaths+I);
            if(uptr->randhits > 0)
                fprintf(stderr,"randhits not initialized\n");
        }
        // following returns the index of the chosen genes; this goes to pointer to all of paths for that gene
        // NEED R SWITCH HERE, FOR SEED
        if (SCALE_UNIVERSE) chooseAugRandList(aug_gene_ct, ingenecnt, 0, randgenesindex);
        else chooseRandList(real_universe_cnt, ingenecnt, seed, randgenesindex); // FIX SEED ISSUE
        //for(j=0;j<ingenecnt;j++) printf("geneindex, gene, %d %d\n", randgenesindex[j], auguniverse[randgenesindex[j]]);
       /* if(i == 1)
        for(j=0;j<ingenecnt;j++)
        {
            printf("%d\n", aug_treepointers[randgenesindex[j]]->val);
        }
        printf("\n\n");*/
        for(j=0;j<ingenecnt;j++)
        {
            if (SCALE_UNIVERSE) thisnode = aug_treepointers[randgenesindex[j]];
            else thisnode = treepointers[randgenesindex[j]];

            //thisnode = treepointers[randgenesindex[j]]; //no needs to point to tree;
            /*if (j < ingenecnt-1 && thisnode->val == treepointers[randgenesindex[j+1]]->val)
            {
                printf("bad\n");
            }*/
            for(k=0;k<thisnode->count;k++)
            {
                thispath = thisnode->all_gene_paths[k];
                (thispath+k)->randhits++;
                // printf("%s %d %d %d %d %d %d\n",thispath->acc, thispath->randhits, thisnode->count, thisnode->val, j,k, randgenesindex[j]);
                // if(!strcmp(thispath->acc, "GO:0034612") && i == 14) printf("%d\n", thisnode->val);
            }
        }
        // now compare randhits with hits
        // could save time looking only at paths with randhits, hits
        for (I=0 ; I<num_used_paths ; I++)
        {         // get number of genes in all pathways, count duplicated gene names
            uptr = (usedpaths+I);
            /*if (!strcmp(uptr->acc, "P00001"))
            {
                printf("comparing counts   \n");
                printf("%d %d %s\n", uptr->hitcnt, uptr->randhits, uptr->acc);
            }*/
            //printf("hits and randhits %d %d\n", uptr->d, uptr->randhits);
            if (uptr->d < uptr->randhits) ++uptr->countover; // more random hits than real hits in this simulation instance; if never, p = 0
            //if(!strcmp(uptr->name, "basal lamina"))
                //printf("%11.2e %d %d %d\n", uptr->pval, uptr->hitcnt, uptr->randhits, uptr->countover);
            else if (uptr->d == uptr->randhits) ++uptr->countequal;
            else ++uptr->countunder;
            // quasi fdr; do for fewer permutes
            if (0) // if (i < NUM_QVAL_PERMUTES)
            {
                rand_d = uptr->randhits;
                var_d = rand_d - uptr->d;
                rand_c = uptr->c - var_d;
                rand_b = uptr->b - var_d;
                a = uptr->a;
                if (rand_b*rand_c == 0) randOR = 99;
                else randOR = a*rand_d/(rand_b*rand_c);
                if (randOR < 1)
                {
                    uptr->randhits = 0;
                    continue;
                }
                //orscore = (randOR - .9)*uptr->numfixedgenes;
                if (var_d > 0) ++d_up;
                else if(var_d < 0) ++d_down;
                else ++d_equal;
                pv = exact22(a,rand_b,rand_c,rand_d);
                //fprintf(predict_fp, "%11.2e \t %11.2e \n", orscore, pv);

                if (pv < SIG_P)
                {
                    ++sig_rand_p_ct;
                    //printf("\n\n %s %d %d %d %d %d %11.2e \n",uptr->name, a, rand_b, rand_c, rand_d, i, pv);
                    for (K=0; K<n_sig_paths; K++)
                    {
                        if(sig_path_ptr[K]->pval >= pv) ++sig_path_ptr[K]->rand_paths_more_sig;
#if 0
                    uptrx = sig_path_ptr[K];
                        unsigned int l;
                        if(!strcmp(sig_path_ptr[K]->name, "basal lamina") && pv < 3.7e-6)
                        {
                            printf("%d\t%d\t%11.2e\t%d\t%d\t%d\t%d\t%d\t%d\t%11.9f\t%11.9f\t%11.7f\t%s\t%s",
                                        i, sig_path_ptr[K]->rand_paths_more_sig, pv, uptr->a, uptr->b, uptr->c, uptr->d, uptr->randhits,
                                       uptr->numfixedgenes,
                                        uptr->pval,  uptr->p_permute_over, uptr->p_permute_under,uptr->acc, uptr->name);
                            //for (l=0;l<uptr->numfixedgenes; l++) printf("%d\n", uptr->egids[l]);
                            printf("\n");
                        }
#endif
                                       
                        

#if 0
gn rpf not used
                        uptrx = sig_path_ptr[K];
#endif
                        //if(!strcmp(sig_path_ptr[K]->name, "basal lamina"))
                        {
                            //printf("%s %11.2e %d\n", (&sig_path_ptr[K])->name, (&sig_path_ptr[K])->pval,(&sig_path_ptr[K])->rand_paths_more_sig);
                        }
                    }
                    //printf("%d %f11.3 %s %d\n", i, pv, (&sig_path_ptr[K])->name, (&sig_path_ptr[K])->rand_paths_more_sig);
                }/**/
            }
            //printf("n, nsig %d %d\n", i, sig_rand_p_ct);
            uptr->randhits = 0;
        }
    }
//fprintf(stderr,"in GPCC 4 predict_fp=%p\n",predict_fp); fflush(NULL);
    if (predict_fp) { fclose(predict_fp); predict_fp = (FILE *)0; }

//fprintf(stderr,"in GPCC 4.0.1\n"); fflush(NULL);
    for (I=0 ; I<num_used_paths ; I++)
    {
//fprintf(stderr,"in GPCC 4.1 %d to %d\n",I,num_used_paths); fflush(NULL);
        uptr = (usedpaths+I);
        uptr->p_permute_over  = (uptr->countover  + 0.5*uptr->countequal + 0.5)/(float) NUM_PERMUTES; // low p value if few random hits > real hits
        uptr->p_permute_under = (uptr->countunder + 0.5*uptr->countequal + 0.5)/(float) NUM_PERMUTES; // casts don't seem to be needed; shouldn't be for large ints
    }
//fprintf(stderr,"in GPCC 5\n"); fflush(NULL);
    for (K=0; K<n_sig_paths; K++)
    {
        sig_path_ptr[K]->pval3  = ((double) sig_path_ptr[K]->rand_paths_more_sig)/((double) NUM_QVAL_PERMUTES);
        //if (sig_path_ptr[K]->pval3 > 0) printf("%s, %d, %11.2e, %11.2e\n", sig_path_ptr[K]->name, sig_path_ptr[K]->rand_paths_more_sig, sig_path_ptr[K]->pval, sig_path_ptr[K]->pval3);
    }
//fprintf(stderr,"in GPCC 7\n"); fflush(NULL);

#if 0
    for (K=0; K<num_used_paths; K++) // DOES THIS DO ANYTHING?
    {
        uptr = (usedpaths+K);
        // if(uptr->pval3 > 0 | uptr->pval < 0.001)  printf("%s\t %d\t %11.2e\t %11.2e\n", uptr->name, uptr->rand_paths_more_sig, uptr->pval, uptr->pval3);
        //if ((uptr->pval3 >= 0) || (uptr->pval < 0.001))  printf("%s\t %d\t %11.2e\t %11.2e\n", uptr->name, uptr->rand_paths_more_sig, uptr->pval, uptr->pval3);
    }
#endif
#ifdef L2P_USING_R
    //PutRNGstate(); // restore random state
#endif
        if (randgenesindex) free(randgenesindex);
        if (head)  free_tree_with_count(head);
        if (pool)     free(pool);
        if (pspace)   free(pspace);
        if (index_ptr) free(index_ptr);
        if (aug_ptr) free(aug_ptr);
        if (aug_treepointers) free (aug_treepointers);
        if (treepointers) free(treepointers);
        if (gene_path_cts) free(gene_path_cts);
        if (auguniverse) free(auguniverse);
        if (augreverse) free(augreverse);
        if (sig_path_ptr) free(sig_path_ptr);
        if (deg_path_cts) free (deg_path_cts);
    // fprintf(stderr,"end GPCC\n"); fflush(NULL);

    return(0);
}

#if 0
old
int fast_permutes(unsigned int ingenecnt, struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt,unsigned int *real_universe)
{
    // pthread_t thr[NUM_THREADS];
    // struct thread_data_t thr_data[NUM_THREADS];
    static unsigned int r[80002];
    double d;
    double oddsratio;
    unsigned int i,j,k;
    double *pspace;
    struct used_path_type *uptr;
    unsigned int usedpathindex;
    double *column;
    double *bestrow;
    unsigned int tmphitcnt;
    unsigned int ui_ej;
    unsigned int ui_uk;

    memcpy(r,real_universe,real_universe_cnt*sizeof(unsigned int));

    pspace = malloc(sizeof(double)*NUM_PERMUTES*num_used_paths); // no dire need to memset
    memset(pspace,0,sizeof(double)*NUM_PERMUTES*num_used_paths);

// fprintf(stderr,"debug in fast_permutes()\n"); fflush(stderr); 
    for (i=0;i<NUM_PERMUTES;i++)
    {
        subsamp(r, real_universe_cnt,ingenecnt); // inlined if optimization on 
        quickSort(r,ingenecnt);                  // inlined if optimization on 
// for (j=0 ; j<ingenecnt ; j++) fprintf(stderr,"ss: %d %d egid= %u\n",i,j,*(r+j));
        for (usedpathindex=0;usedpathindex<num_used_paths;usedpathindex++)
        {
            uptr = (usedpaths+usedpathindex);
            j = k = tmphitcnt = 0;
            while ((j<uptr->numfixedgenes) && (k < ingenecnt))
            {              // get the number of hits 
                ui_ej = *(uptr->egids+j); // entrez gene ids in the jth used pathway
                ui_uk = *(r+k);           // entrez gene from random sampling 
                if (ui_ej == ui_uk)
                {
                    tmphitcnt++;
                    k++;
                    j++;
                    continue;
                }
                else if (ui_ej < ui_uk) j++;
                else                    k++;
            }
	        /* d = exact22((int)tmphitcnt,uptr->numfixedgenes-tmphitcnt,ingenecnt,real_universe_cnt-ingenecnt);
                   *(pspace + (usedpathindex*NUM_PERMUTES) + i) = (double)d; */
	    d = ((double)tmphitcnt / (double)ingenecnt) / 
		 ((double)(uptr->numfixedgenes-tmphitcnt) / (double)(real_universe_cnt-ingenecnt));
// fprintf(stderr,"OR=%20.18f from (%u/%u)/(%u/%u)\n",oddsratio,tmphitcnt,ingenecnt,uptr->numfixedgenes-tmphitcnt,real_universe_cnt-ingenecnt);
            *(pspace + (usedpathindex*NUM_PERMUTES) + i) = (double)d;

/*
               diseased healthy
exposed           DE       HE            A   C
notexposed        DN       HN            B   D
OR=(DE/HE)/(DN/HN) = A/C / B/D

list=300 hits=3  pw=40 genome=30001
         usergenes     genome
ipw         3-1         40               tmphitcnt                        ingenecnt
notinpw     297       29960              uptr->numfixedgenes-tmphitcnt    real_universe_cnt-ingenecnt 
*/


#if 0
fprintf(stderr,"exact22 %f %d %d %d %d upi=%d i=%d %s\n",d,
		tmphitcnt,uptr->numfixedgenes-tmphitcnt,ingenecnt,real_universe_cnt-ingenecnt,
		usedpathindex,i,
		uptr->name);
#endif
        }

    }

#if 0
 for (i=0 ; i<NUM_PERMUTES ; i++) 
 { 
         fprintf(stderr,"perm:%d ",i);
         for (usedpathindex=0 ; usedpathindex<num_used_paths ; usedpathindex++)
         {
             d = *(pspace + (usedpathindex*NUM_PERMUTES) + i);
             fprintf(stderr,"%9.7f ",d);
         }
         fprintf(stderr,"\n");
 }
#endif

    column = malloc(sizeof(double)*num_used_paths); // no dire need to memset
    bestrow = malloc(sizeof(double)*NUM_PERMUTES); // no dire need to memset
    for (i=0;i<NUM_PERMUTES;i++)
    {
        for (usedpathindex=0;usedpathindex<num_used_paths;usedpathindex++)
        {
            d = *(pspace + (usedpathindex*NUM_PERMUTES) + i);
            *(column+usedpathindex) = d;
        }
        qsort(column,num_used_paths,sizeof(double),cmp_double);
        *(bestrow+i) = *(column+0);
    }
    qsort(bestrow,NUM_PERMUTES,sizeof(double),cmp_double);
#if 0
for (j=0 ; j<NUM_PERMUTES ; j++) { fprintf(stderr,"%d %f\n",j,*(bestrow+j)); }
#endif
    d = 1.0;
    for (usedpathindex=0 ; usedpathindex<num_used_paths ; usedpathindex++)
    {
        uptr = (usedpaths+usedpathindex);
        oddsratio = ((double)uptr->hitcnt / (double)ingenecnt) /
		 ((double)(uptr->numfixedgenes-(uptr->hitcnt)) / (double)(real_universe_cnt-ingenecnt));
        for (j=0 ; j<NUM_PERMUTES ; j++) 
        {
            if ( *(bestrow+j) < oddsratio )      // if (*(bestrow+j) > uptr->pval) 
            {
                d = *(bestrow+j); 
                break;
            }
        }
        d = (double)j/(double)NUM_PERMUTES;
        uptr->fdr2 = d;
    }

    free(pspace);
    free(column);
    free(bestrow);

    return 0;

#if 0
        k = NUM_PERMUTES/NUM_THREADS;
        for (j = 0; j < NUM_THREADS; ++j) 
        {     /* create threads */
            thr_data[j].lo = j*k;
            thr_data[j].hi = (j*k)+k-1;
            if (j == (NUM_THREADS-1))
                thr_data[j].hi = num_used_paths-1; // note "-1", 
fprintf(stderr,"%d of %d [ %d .. %d ] of %d size= %d\n",j,NUM_THREADS,thr_data[j].lo,thr_data[j].hi,num_used_paths,thr_data[j].hi-thr_data[j].lo); 
            thr_data[j].myid = j;
            thr_data[j].ingenecnt = ingenecnt;
            thr_data[j].sampling = r;
            thr_data[j].usedpathsindex = usedpaths;
            thr_data[j].usedpaths = usedpaths;
            thr_data[j].num_used_paths = num_used_paths;
fprintf(stderr,"creating thread %d\n",j); fflush(stderr); 
            if ((rc = pthread_create(&thr[j], NULL, thr_func, &thr_data[j]))) 
            {
              fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
              return EXIT_FAILURE;
            }
        }
      /* block until all threads complete */
        for (j = 0; j < NUM_THREADS; ++j) 
        {
fprintf(stderr,"blocking %d \n",j);
            pthread_join(thr[j], NULL);
        }
fprintf(stderr,"after blocking  \n");
    return EXIT_SUCCESS;
#endif
}
#endif
 

          /* Be sure to free egids if you call setup_by_egids() !!! */
int setup_by_egids(void)
{
    size_t sz = sizeof(struct smallgenetype) * numgenes;

    by_egids = (struct smallgenetype *)malloc(sz);
    if (!by_egids) return -1;
    memcpy(by_egids,genes,sz);
    qsort(by_egids,numgenes,sizeof(struct smallgenetype),cmp_by_egid);
    return 0;
}


unsigned int  tmpgenes[MAX_INGENES];
 
#if 0
int l2pfunc(struct used_path_type *usedpaths,unsigned int num_used_paths,unsigned int real_universe_cnt,unsigned int *real_universe, int permute_flag, int *user_incnt_ptr, int oneside)
#endif
int l2pfunc(struct used_path_type *usedpaths,unsigned int num_used_paths,unsigned int real_universe_cnt,unsigned int *real_universe, int permute_flag, int *user_incnt_ptr, int oneside,unsigned int seed)
{
    char s[512];
    struct used_path_type *uptr; // used path pointer 
    unsigned int prev;
    unsigned int ui;
    unsigned int ui_ej;
    unsigned int user_incnt = 0;
    unsigned int i,j,k,ll;
    int ret = 0;
    int overflowflag = 0;
    unsigned int *uiptr = (void *)0;
// struct timespec time1, time2;
    // char *z;

    *user_incnt_ptr = j = 0;
    // ADDING INPUT FILE FOR C RUN, INPUT GENE LIST
    // struct timespec time1, time2;
    // char *z;
    // printf("input DEG list\n");
#if NELSON_C
    FILE *deglist;
    if ((deglist = fopen (DEGLISTNAME, "r")) == NULL){
        printf ("can't open deglist, exiting");
        exit (1);
    }
    // READING FROM DEGLIST
    while ( fgets(s, sizeof(s), deglist) ) // gets() function is deprecated
#else
    while ( fgets(s, sizeof(s), stdin) ) // gets() function is deprecated
#endif
    {
        for (i=0 ; s[i] ; i++) 
        { 
             if ((s[i] == '\n')||(s[i] == '\r')) s[i] = (char)0; 
             else if ((s[i] < ' ') || (s[i] >= 127)) s[i] = ' ';
        }
        ui = hugo2egid(s);
        if (ui == (unsigned int)UINT_MAX)
        {
             fprintf(stderr,"Note: invalid gene \"%s\" in user input\n",s);  
             continue;
        }
        uiptr = (unsigned int *)bsearch(&ui,real_universe,real_universe_cnt,sizeof(unsigned int),cmp_ui);
        if (!uiptr)
        {
             fprintf(stderr,"Note: gene not in universe : \"%s\" \n",s);  
             continue;
        }
        if (j < MAX_INGENES)
        {
            tmpgenes[j++] = ui;
        }
        else
        {
            if (overflowflag == 1) { fprintf(stderr,"NOTE: too many genes, input gene number limited to %d, rest are ignored, incnt=%d\n",MAX_INGENES,j); fflush(NULL); }
            overflowflag = 1;
        }
    }
    qsort(tmpgenes,j,sizeof(unsigned int),cmp_ui); // sort user input egids (entrez gene ids)
    prev = 0;
    for (i=user_incnt=0;i<j;i++)
    {
        ui = tmpgenes[i];
        if (ui != prev)
        {
            ingenes[user_incnt++] = ui;
        }
        prev = ui;
    }
    ingenecnt = user_incnt;
// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    for (i=0 ; i<num_used_paths ; i++)
    {
        uptr = (usedpaths+i);
        if (!uptr->egids) continue;
        j = k = ll = 0;
        while ((j<uptr->numfixedgenes) && (k < user_incnt))
        {
            ui_ej = *(uptr->egids+j);
            ui = ingenes[k];
            if (ui_ej == ui)
            {
                *((uptr->genehits) + (ll++)) = ui; // remember, because need to print out later
                k++;
                j++;
                continue;
            }
            else if (ui_ej < ui) j++;
            else                 k++;
        }
        uptr->hitcnt = ll;
    }
#if 0
clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
#define BILLION  1000000000.0;
double time_spent = (time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec) / BILLION;
fprintf(stderr,"Time elapsed is %f seconds\n", time_spent);  fflush(NULL);

fprintf(stderr,"here , before do_pvals %d %d %p %d \n",real_universe_cnt , user_incnt, usedpaths,num_used_paths); 
fflush(NULL);
#endif

    *user_incnt_ptr = user_incnt;
    if (permute_flag == 0)
    {
           do_pvals_and_bh(user_incnt, usedpaths,num_used_paths, real_universe_cnt,oneside);
    }
    else
    {
#if 0
        GPCC(usedpaths,num_used_paths,real_universe_cnt, real_universe);
#endif
        //GPCC(user_incnt,usedpaths,num_used_paths,real_universe_cnt,seed);
//fprintf(stderr,"in lp2func(), before GPCC\n"); fflush(NULL);
        GPCC(usedpaths,num_used_paths,real_universe_cnt, real_universe, seed);
//fprintf(stderr,"in lp2func(), after  GPCC\n"); fflush(NULL);

    // RICH: this code does all comparisons out of GPCC
    //do_just_bh(user_incnt,usedpaths,num_used_paths,real_universe_cnt);
    }
    return ret;
}




#ifdef WEBASSEMBLY
int testwasm(void)
{
int i;
// need to bring some symbols for functions to test out wasm compilation
for (i=0;i<10;i++) printf("%s",genes[i].hugo); // wasm 
for (i=0;i<10;i++) printf("%s",pws[i].name); // wasm 
for (i=0;i<10;i++) printf("%d",pwgenes[i]); // wasm 
return 0;
}
#else
void usage(void)
{
fprintf(stderr,"l2p : \"list to pathways\" program.\n");
fprintf(stderr,"Usage: cat listofHUGOgenes_one_per_line.txt | l2p [optional args]\n");
fprintf(stderr,"possible optional args are ...\n");
fprintf(stderr," -help\n");
fprintf(stderr," -precise\n");
fprintf(stderr," -justheader\n");
fprintf(stderr," -noheader\n");
fprintf(stderr," -universe=Universefile_one_gene_per_line.txt\n");
fprintf(stderr," -categories=listofcatgories  (comma separated)\n");
fprintf(stderr,"    example: -categories=KEGG,REACTOME,BIOCYC,PANTH  (i.e. only use genes in those 4 categories\n");
fprintf(stderr,"    another -categories example: \"-categories=H,C6\"  (i.e. only use msigdb's Hallmark and C6 (cancer) category pathways\n");
fprintf(stderr,"    available categories are :\n");
fprintf(stderr,"    BIOCYC  - organism specific Pathway/ Genome Databases (PGDBs)  - https://biocyc.org/\n");
fprintf(stderr,"    GO  - initiative to unify representation of gene and gene product attributes -  http://geneontology.org\n");
fprintf(stderr,"    KEGG - databases dealing with genomes, biological pathways, - https://www.kegg.jp/\n");
fprintf(stderr,"    PANTH - databases for protein analysis through evolutionary relationships - http://www.pantherdb.org/\n");
fprintf(stderr,"    PID  - Pathway interaction database: legacy database from Carl Schaefer & buddies at NCI\n");
fprintf(stderr,"    REACTOME - curated database of biological pathways - https://reactome.org/\n");
fprintf(stderr,"    WikiPathways - community resource for biological pathways - https://www.wikipathways.org\n");
fprintf(stderr,"    C1 - MSigDB only, positional gene sets for each human chromosome and cytogenetic band.\n");
fprintf(stderr,"    C2 - MSigDB only, curated gene sets from online pathway databases, publications in PubMed, and experts.\n");
fprintf(stderr,"    C3 - MSigDB only, motif gene sets based on conserved cis-regulatory motifs from comparative analysis\n");
fprintf(stderr,"    C4 - MSigDB only, computational gene sets defined by mining large collections of cancer-oriented microarray data.\n");
fprintf(stderr,"    C5 - MSigDB only, gene sets  consist of genes annotated by the same GO terms.\n");
fprintf(stderr,"    C6 - MSigDB only, oncogenic gene sets defined directly from microarray data from cancer gene perturbations.\n");
fprintf(stderr,"    C7 - MSigDB only, immunologic gene sets  from microarray data from immunologic studies.\n");
fprintf(stderr,"    H - MSigDB only, hallmark gene sets: signatures from MSigDB gene sets to represent biological processes.\n");
fprintf(stderr,"Example run : printf \"TP53\\nPTEN\\nAPC\\nKRAS\\nNRAS\\n\" | ./l2p -categories=PID | sort -k2,2n -k1,1n -k3,3nr | head\n");
            fflush(NULL);
}

int bitCount(int n)
{
    int cnt = 0;
    while (n)
    {
        cnt += n % 2;
        n >>= 1;
    }
    return cnt;
}

static int parsecats(char *z, unsigned int *catspat)
{
    char ts[16][16];
    int bit;
    int j,k;
    int toks;

    for (j=0 ; *(z+j) ; j++)
    {
        if (*(z+j) == ',') *(z+j) = ' ';
    }
    memset(ts,0,sizeof(ts));
    toks = sscanf(z,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s ",
          ts[0], ts[1], ts[2], ts[3], ts[4], ts[5], ts[6], ts[7], ts[8], ts[9],
          ts[10], ts[11], ts[12], ts[13], ts[14], ts[15]);
    for (k=0;k<toks;k++)
    {
        bit=string_to_category_code(ts[k]);
        if (bit)
            *catspat |= bit;
        else
        {
            fprintf(stderr,"ERROR: invalid category = %s\n",ts[k]);
            usage();
            exit(0);
        }
    }
    j = bitCount(*catspat) ;
    if (j == 0)
    {
#ifdef L2P_USING_R
        return 0;
#else
        fprintf(stderr,"ERROR: no categories specified\n");
        usage();
        return -1;
#endif
    }
    return 0;
}


void print_header(void)
{
 /* 1  */ printf("pval\t"); 
 /* 2  */ printf("fdr\t");
 /* 3  */ printf("enrichment_score\t");        // if positive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
 /* 4  */ printf("pwhitcount\t");              // number of genes hit in pathway
 /* 5  */ printf("pwnohitcount\t");            // pathway number of genes in the pathway
 /* 6  */ printf("inputcount\t");              // total count of user genes (user input)
 /* 7  */ printf("pwuniverseminuslist\t");     // total number of unique genes in all pathways
 /* 8  */ printf("pathwayaccessionidentifier\t");// canonical accession ( if available, otherwise assigned by us )
 /* 9  */ printf("category\t");                // KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  *was "source"*
 /* 10 */ printf("pathwayname\t");             // Name of pathway
 /* 11 */ printf("genes\t");               // genes_space_separated   HUGO genes from user that hit the pathway
          printf("\n");
}

int l2p_init_C(int argc,char *argv[],unsigned int *catspatptr,int *precise_flag, int *permute_flag,int *no_header_flag, char universe_file[], char custom_file[], int *get_universe_flag, int *oneside, unsigned int *seed)
{
    char *z;
    int i;
    int oneerr = 0;

    *get_universe_flag = 0;
    *precise_flag = 0;
    *permute_flag = 0;
    *no_header_flag = 0;
    *oneside = 0; // default is two-sided fe
    *seed = 0;
    universe_file[0] = (char)0;
    custom_file[0] = (char)0;
    for (i=1 ; i<argc ; i++)
    {
        if ((strcmp(argv[i],"-help") == 0) || (strcmp(argv[i],"--help") == 0) || ((strcmp(argv[i],"help") == 0) ) )
        {
             usage();
             exit(1);
        }
        else if (strncmp(argv[i],"-categories=",12) == 0)
        {
            z = argv[i],
            z += 12;
            parsecats(z,catspatptr);
        }
        else if (strcmp(argv[i],"-precise") == 0)
              *precise_flag = 1; // print more digits out some user doesn't complain about "real pvals"
        else if (strcmp(argv[i],"-permute") == 0)
              *permute_flag = 1; // print more digits out some user doesn't complain about "real pvals"
        else if (strcmp(argv[i],"-noheader") == 0)
              *no_header_flag = 1;
        else if (strcmp(argv[i],"-justheader") == 0)
        {
            print_header();
            exit(0);
        }
        else if (strncmp(argv[i],"-oneside=",9) == 0)
        {
            oneerr = 0;
            if (argv[i] + 9)
            {
                *oneside = atoi(argv[i] + 9);
                if ( (*oneside < 0) || (*oneside > 2) )
                    oneerr = 1;
            }
            else 
            {
                oneerr = 1;
            }
            if (oneerr)
            {
                oneside = 0;
                fprintf(stderr,"ERROR: invalid one sided flag , using two-sided\n"); fflush(stderr);
            }
        }
        else if (strncmp(argv[i],"-universe=",10) == 0)
        {
            strcpy(universe_file,argv[i] + 10);
            // user_universe_flag = 1;
fprintf(stderr,"note: using \"%s\" as universe file\n",universe_file);
        }
        else if (strncmp(argv[i],"-seed=",6) == 0)
        {
            *seed = (unsigned int)atoi(argv[i]+6);
fprintf(stderr,"seed is %u\n",*seed); 
        }
        else if (strcmp(argv[i],"-getuniverse") == 0)
        {
            *get_universe_flag = 1;
        }
        else if (strncmp(argv[i],"-customfile=",12) == 0)
        {
            strcpy(custom_file,argv[i] + 12);
//             customflag++;
fprintf(stderr,"note: using \"%s\" as custom pathway gene listfile\n",custom_file);
        }
        else 
        {
fprintf(stderr,"Note: Invalid argument \"%s\" : ignored.\n",argv[i]);
        }
    }
    if (*catspatptr == 0) category_set_all(catspatptr);
    return 0;
}


void print_universe(unsigned int real_universe_cnt, unsigned int *ru)
{              // not a debug routine, used in command line
    unsigned int i;
    char *z;
// fprintf(stderr,"in print_universe real_universe_cnt = %d \n",real_universe_cnt);  fflush(stderr); 
    for (i=0;i<real_universe_cnt;i++)
    {
        z = egid2hugo(*(ru+i));
        if (z) printf("%s\n",z);
    }
    return;
}


unsigned int *get_used_universe(struct used_path_type *u, unsigned int num_used, unsigned int *real_universe_cnt_ptr)
{
    unsigned int i,j,k;
    size_t sz = 0;
    unsigned int prev = 0;
    unsigned int ui = 0;
    unsigned int maxelems = 0;
    unsigned int *tempverse  = (unsigned int *)0;
    unsigned int *real_universe  = (unsigned int *)0;
    struct used_path_type *uptr = (void *)0;    // used path pointer


    if (!u)
    {
        fprintf(stderr,"ERROR: null used_path_type passed to get_used_universe()\n"); fflush(stderr);
        return (unsigned int *)0 ;
    }

    maxelems = numpwgenes * 2;
    sz = (maxelems * sizeof(unsigned int)); //  maximum possible number of genes
    tempverse = malloc(sz);
    if (!tempverse)
    {
        fprintf(stderr,"ERROR: not enough memory in get_used_universe()\n"); fflush(stderr);
        return (unsigned int *)0 ;
    }
    memset(tempverse,0,sz);
    for (i=k=0;i<num_used;i++)
    {
        uptr = (u+i);
        for (j=0 ; j < uptr->numfixedgenes ; j++)
        {
            *(tempverse+k) = *(uptr->egids + j);
            k++;
            if (k >= maxelems)
            {
fprintf(stderr,"in get_used_universe ERROR i=%d k=%d\n",i,k);  fflush(stderr);
exit(0);
            }
        }
    }
    qsort(tempverse,k,sizeof(unsigned int),cmp_ui);
    real_universe = (unsigned int *)malloc(sz);
    prev = 0;
    for (i=j=0;i<k;i++)
    {
       ui = *(tempverse+i);
       if (ui != prev)
       {
           *(real_universe+j) = ui;
           j++;
       }
       prev = ui;
    }
    *real_universe_cnt_ptr = j;
    free(tempverse);
    return (unsigned int *)real_universe;
}

unsigned int *load_user_universe_file(char universe_file[], unsigned int *universe_cnt_ret)
{
    char s[130000];   // 127570 is number of bytes in universe file, so should handle anything
    unsigned int ui;
    FILE *universe_fp;
    unsigned int *universe_genes = (unsigned int *)0;
    unsigned int universe_cnt;
    int i,j;

    *universe_cnt_ret = j = 0;
    if (universe_file[0])
    {                 // we have a universe file
        universe_fp = fopen(universe_file,"r");
        if (!universe_fp) 
        {
fprintf(stderr,"ERROR: Can't open universe_file named \"%s\"\n",universe_file); 
            return (void *)0;
        }
        universe_cnt = 0;
        while ( fgets(s, sizeof(s), universe_fp) )
            universe_cnt++;
        universe_genes = malloc(sizeof(unsigned int *)*universe_cnt);
        rewind( universe_fp);
        j = 0;
        while ( fgets(s, sizeof(s), universe_fp) )
        {
            for (i=0 ; s[i] ; i++) { if ((s[i] == '\n')||(s[i] == '\r')) s[i] = (char)0; }
// if (ui==29974) { fprintf(stderr,"got A1CF 29974\n"); exit(0); }
            ui = hugo2egid(s);
            if (ui == (unsigned int)UINT_MAX)
            {
                fprintf(stderr,"Note: invalid gene \"%s\" in universe file.\n",s);  
                continue;
            }
            *(universe_genes+j) = ui;
            j++;
        }
        fclose (universe_fp);
    }
    qsort(universe_genes,j,sizeof(unsigned int),cmp_ui);
    *universe_cnt_ret = j; // repair this , some genes may have been "invalid". initial count may have been wrong.
    return universe_genes;
}

 
int count_lines_in_file(char *fn)
{
    FILE *fp;
    int lines = 0;
    int ch = 0;

    fp = fopen(fn,"r");
    if (!fp) return 0;
    while(!feof(fp))
    {
        ch = fgetc(fp);
        if (ch == '\n') lines++;
    }
    fclose(fp);
    return lines;
}


struct used_path_type *setup_used_paths(unsigned int *num_used_paths, unsigned int catspat, char universe_file[], unsigned int in_universe_cnt,unsigned int *in_universe, char custom_file[], unsigned int *real_universe_cnt_ptr,unsigned int **real_universe,unsigned int lencust,struct custom_type *mycustompw)
{
    char s[100000];
    char tmps[5120];
    char custname[5120];
    unsigned int i,j,k,ll,used_index;
    int tokcnt;
    struct custom_type *kust = (struct custom_type *)0;
    char *z = (void *)0;
    int this_cust_pw_cnt = 0;
    unsigned int ui = 0;
    unsigned int *uiptr = (void *)0;
    unsigned short int usi = (unsigned short int)0;
    struct used_path_type *u = (struct used_path_type *)0;
    struct used_path_type *uptr; // used path pointer
    unsigned int custom_cnt = 0;
    unsigned int universe_cnt = 0;
    unsigned int *user_universe_genes = (unsigned int *)0;
    unsigned int *this_egids = (unsigned int *)0;
    unsigned int *tmpugenes = (unsigned int *)0;
    unsigned int prev = 0;
    unsigned int newcnt = 0;
    int keepgoing;
    FILE *fp;
        // char **custom_data = (void *)0;
        // struct timespec time1, time2;

// fprintf(stderr,"rpf debug mycustompw=%p lencust=%d\n",mycustompw,lencust); fflush(stderr);  
// fprintf(stderr,"in setup_used_paths() start catspat=%d = 0x%x\n",catspat,catspat); fflush(stderr); 
        // clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    *real_universe = (unsigned int *)0;

              // just need to know how many lines in gmt file
    if (custom_file[0]) 
    {
        custom_cnt = count_lines_in_file(custom_file);
// fprintf(stderr,"rpf debug custom_cnt=%d from %s\n",custom_cnt,custom_file); fflush(stderr);  
    }
    else
    {
        if (mycustompw) 
        {
// fprintf(stderr,"rpf debug mycustompw=%p\n",mycustompw); fflush(stderr);  
        }
    }

    if (universe_file[0])
    {
        user_universe_genes = load_user_universe_file(universe_file,&universe_cnt); // may return null 
    }
    else if (in_universe_cnt != 0)
    {
        user_universe_genes = in_universe;
	universe_cnt = in_universe_cnt;
    }
    u = (struct used_path_type *)malloc(sizeof(struct used_path_type) * (numpws+custom_cnt+lencust));
    if (!u) { fprintf(stderr,"ERROR: not enough memory.\n"); return (void *)0; }
    memset(u,0,sizeof(struct used_path_type) * (numpws+custom_cnt+lencust));
    used_index = 0;
    for (i=0 ; i<numpws ; i++)
    {
        if ((pws[i].category & catspat) == 0) continue;
        this_egids = (unsigned int *)malloc(sizeof(unsigned int)*pws[i].numgenes);
        if (this_egids == (void *)0)
        {
                 exit(0);
        }
        if (universe_cnt)   // get rid of unused genes - use only universe genes
        {
            for (newcnt = j = 0 ; j < pws[i].numgenes ; j++)
            {
                usi = pwgenes[pws[i].pwgenesindex+j];
                ui = genes[usi].egid;
                uiptr = (unsigned int *)bsearch(&ui,user_universe_genes,universe_cnt,sizeof(unsigned int),cmp_ui);
                if (uiptr)
                {
                    *(this_egids+newcnt) = ui; // put this in the right place
                    newcnt++;
                }
                else
                {
                }
            }
            if (newcnt < 3)  // not enough genes
            {
                free(this_egids);
                this_egids = (unsigned int *)0;
                continue;
            }
        }
        else
        {
            for (newcnt = j = 0 ; j<pws[i].numgenes ; j++)
            {
                ll = pwgenes[pws[i].pwgenesindex + j];
                if (ll == USHRT_MAX) // masked out as not in universe
                {
                    fprintf(stderr,"ERROR: %u \n",ll); exit(0);
                    continue;
                }
                *(this_egids+newcnt) = genes[ll].egid;
                newcnt++;
            }
        }
        uptr = (u+used_index);

        qsort(this_egids,newcnt,sizeof(unsigned int),cmp_ui);
        uptr->egids = this_egids;
        uptr->numfixedgenes = newcnt;
        uptr->genehits = (unsigned int *)malloc(sizeof(unsigned int)*newcnt);
        uptr->category = pws[i].category;
        uptr->acc  = strdup(pws[i].acc);       // can't just assign (?)  . because need to free if custom
        uptr->name = strdup(pws[i].name);      // can't just assign
        uptr->numgenes = pws[i].numgenes;
        uptr->pwgenesindex = pws[i].pwgenesindex; // not going to care about this probably. info is now in egids
        uptr->hitcnt = 0;
        uptr->gpcountsum = 0;
        used_index++;
    }
// fprintf(stderr,"rpf used_index = %d\n",used_index); fflush(stderr); 

#define MAXUGENES 40000
    if (custom_cnt)
    {
        if (tmpugenes) { free(tmpugenes); tmpugenes = (unsigned int *)0; }
        tmpugenes=malloc(sizeof(unsigned int)*MAXUGENES);
        if (!tmpugenes) { fprintf(stderr,"can't malloc %u bytes \n",(unsigned int)(sizeof(unsigned int)*MAXUGENES));  fflush(stderr); }
    fp = fopen(custom_file,"r");
    if (!fp) { fprintf(stderr,"NOTE: can not open \"%s\" - ignoring\n",custom_file); fflush(stderr); }
    while ( fgets(s, sizeof(s)-2, fp) ) // for each line of custom file
    {
        for (k=0 ; s[k] ; k++) { if ((s[k] == '\n')||(s[k] == '\r')) s[k] = (char)0; } // strip newline
        uptr = (u+used_index); // used_index points to next available slot 
        z = &s[0];
        j = tokcnt = k = 0;
        tmps[0] = (char)0;
        keepgoing = 1;
        newcnt = 0;
        while (keepgoing == 1)
        {
            if ( (*(z+j) == '\t') || (*(z+j) == (char)0) )
            {   // first token is name, 2nd is optional and rest are gene names separated by tabs
               if (tokcnt == 0) strcpy(custname,tmps); 
               else if (tokcnt == 1) { /* ignore */ } 
               else
               {
                   if (tmps[0])
                   {
                       ui = hugo2egid(tmps);
                       if (ui == (unsigned int)UINT_MAX)
                       {
                           fprintf(stderr,"Note: invalid gene \"%s\" in custom pathway file.\n",tmps);  
                       }
                       else
                       {
                           if (universe_cnt) // must be in user universe
                           {
                               uiptr = (unsigned int *)bsearch(&ui,user_universe_genes,universe_cnt,sizeof(unsigned int),cmp_ui);
                               if (uiptr)
                               {
                                   *(tmpugenes+newcnt) = ui; // put this in the right place
                                   newcnt++;
                                   if (newcnt == MAXUGENES) { fprintf(stderr,"ERROR: too many genes on custom pathway gmt file.\n"); fflush(stderr); }
                               }
                           }
                           else 
                           {
                               *(tmpugenes+newcnt) = ui; // put this in the right place
                               newcnt++;
                           }
                       }
                   }
               }
               tokcnt++;
               k = 0;
               tmps[0] = (char)0;
               if (*(z+j) == (char)0) keepgoing = 0;
            }
            else
            {
                tmps[k++] = *(z+j); 
                tmps[k] = (char)0;
            }
            j++;
        }
        if (newcnt == 0) continue;
        qsort(tmpugenes,newcnt,sizeof(unsigned int),cmp_ui);
        this_cust_pw_cnt = 0;
        prev = 0;
        this_egids = (unsigned int *)malloc(sizeof(unsigned int)*newcnt);
        for (k=0;k<newcnt;k++) // de-duplicate, put only UNIQUE genes in
        {
            ui = *(tmpugenes+k); 
            if (ui != prev)
            {
                *(this_egids+this_cust_pw_cnt) = ui;
                this_cust_pw_cnt++;
            }
            prev = ui;
        }
        uptr->acc = strdup(custname); // remember to free this
        uptr->name =strdup(custname);  // remember to free this
        uptr->egids = this_egids; // remember to free this
        uptr->numgenes = uptr->numfixedgenes = this_cust_pw_cnt;
        uptr->pwgenesindex = USHRT_MAX;
        uptr->genehits = (unsigned int *)malloc(sizeof(unsigned int)*this_cust_pw_cnt); // remember to free this
        uptr->hitcnt = 0;
        uptr->category = CAT_CUSTOM | (string2type("custsom") << 28);
        used_index++;
    }
    if (tmpugenes) { free(tmpugenes); tmpugenes = (unsigned int *)0; }
    }

    for (i=0;i<lencust;i++) // add R user custom lists passed as a list of vectors (with "gmt" lines)
    {
        uptr = (u+used_index); // used_index points to next available slot 
        kust = (mycustompw + i);
        this_egids = (unsigned int *)malloc(sizeof(unsigned int) * kust->numgenes);
        for (j=0 ; j<kust->numgenes ; j++)   // for each line of custom file
        {
            *(this_egids+j) = *(kust->genes+j);
// fprintf(stderr,"rpf debug i=%d cust=%d egid=%d\n",i,j,*(this_egids+j)); fflush(stderr);  
        }
        uptr->acc = strdup(kust->name);
        uptr->name = strdup(kust->name);
        uptr->egids = this_egids; // remember to free this
        uptr->numgenes = uptr->numfixedgenes = kust->numgenes;
        uptr->pwgenesindex = USHRT_MAX;
        uptr->genehits = (unsigned int *)malloc(sizeof(unsigned int)*kust->numgenes);
        uptr->hitcnt = 0;
        uptr->category = CAT_CUSTOM | (string2type("custsom") << 28);
        used_index++;
    }
    *num_used_paths = used_index;
    if (user_universe_genes) free(user_universe_genes);

    *real_universe = get_used_universe(u, used_index, real_universe_cnt_ptr);
#if 0
clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
#define BILLION  1000000000.0;
double time_spent = (time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec) / BILLION;
fprintf(stderr,"Time elapsed is %f seconds\n", time_spent); 
#endif
    return u;
}


void print_permute_header(void)
{
    printf("a\t"); // uptr->a 
    printf("c\t"); // uptr->c 
    printf("b\t"); // uptr->b 
    printf("d\t"); // uptr->d 
    printf("countunder\t"); // uptr->countunder 
    printf("countequal\t"); // uptr->countequal 
    printf("countover\t"); // uptr->countover 
    printf("numfixedgenes\t"); // uptr->numfixedgenes 
    printf("OR\t"); // uptr->OR 
    printf("pval\t"); // uptr->pval 
    printf("p_permute_over\t"); // uptr->p_permute_over 
    printf("p_permute_under\t"); // uptr->p_permute_under 
    printf("pval3\t"); // uptr->pval3 
    printf("acc\t"); // uptr->acc 
    printf("tmps\t"); // tmps 
    printf("name"); // uptr->name); 
    printf("\n");
} 

int main(int argc,char *argv[])
{
    char universe_file[PATH_MAX];
    char custom_file[PATH_MAX];
    char tmps[PATH_MAX];
    int user_incnt = 0;
    int no_header_flag = 0;
    int precise_flag = 0;
    int permute_flag = 0;
    int oneside = 0;
    unsigned int seed = 0;
    unsigned int catspat = 0;
    unsigned int num_used_paths;
    struct used_path_type *u;
    struct used_path_type *uptr;
    int status;
    unsigned int *real_universe;
    unsigned int real_universe_cnt;
    int get_universe_flag = 0;
    unsigned int i = 0;
    unsigned int j = 0;
    char *z;
    double fdr_for_output;
#if NELSON_C
    double corr_tst, mean_gpcts;
    FILE *pathcalls;
#endif


    (void)setup_by_egids();
    l2p_init_C(argc,argv,&catspat,&precise_flag,&permute_flag,&no_header_flag,universe_file,custom_file,&get_universe_flag,&oneside,&seed);
    u = setup_used_paths(&num_used_paths, catspat,universe_file, 0,(void *)0,custom_file,&real_universe_cnt,&real_universe,0,(struct custom_type *)0);
    if (get_universe_flag == 1)
    {
        print_universe(real_universe_cnt,real_universe);
        return 0;
    }
#if 0
    status = l2pfunc(u,num_used_paths,real_universe_cnt,real_universe,permute_flag,&user_incnt,oneside);
#endif
    status = l2pfunc(u,num_used_paths,real_universe_cnt,real_universe,permute_flag,&user_incnt,oneside,seed);

#if NELSON_C
    pathcalls = fopen(OUTPUTFILE, "w");
    fprintf(pathcalls, "a\t c\t v\t d\t countunder \t countequar \t countover \t pwgenect \t OR \t FET p \t permute p \t mean_gpcts \t gp corr p\n");

#endif
    if (permute_flag)
    {
        if (no_header_flag == 0) print_permute_header();
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = (u+i);
            categories_pattern_to_strings(uptr->category,tmps);
            fdr_for_output = uptr->fdr;
            printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%11.2f\t%11.2e\t%11.2e\t%11.2e\t%11.2e\t%s\t%s\t%s\t",
                          uptr->a, uptr->c, uptr->b, uptr->d,
                          uptr->countunder, uptr->countequal, uptr->countover, uptr->numfixedgenes, uptr->OR,
                          uptr->pval,  uptr->p_permute_over, uptr->p_permute_under, uptr->pval3,
                          uptr->acc, tmps, uptr->name);
#if NELSON_C
            mean_gpcts = uptr->gpcountsum/((float) uptr->hitcnt);
            //corr_tst = pow(uptr->gpcountsum/uptr->hitcnt/mean_u_pathcounts, uptr->hitcnt/5);
                fprintf(pathcalls, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%11.2f\t%11.2e\t%11.2e\t%11.2e\t%11.2e\t%s\t%s\t%s\t\n",
                 uptr->a, uptr->c, uptr->b, uptr->d,
                uptr->countunder, uptr->countequal, uptr->countover, uptr->numfixedgenes, uptr->OR,
                  uptr->pval,  uptr->p_permute_over,  mean_gpcts,
                        uptr->pval4,  uptr->acc, tmps, uptr->name);
#endif
fprintf(stderr,"hitcnt = %d\n",uptr->hitcnt); 
            for (j=0 ; j<uptr->hitcnt ; j++)
            {
                z = egid2hugo(*((uptr->genehits)+j));
                printf("%s ",z);
            }
            printf("\n");
        }
    }
    else
    {
        if (no_header_flag == 0) print_header();
        for (i=0 ; i<num_used_paths ; i++)
        {
           uptr = (u+i);
           categories_pattern_to_strings(uptr->category,tmps);
           fdr_for_output = uptr->fdr;
           if (precise_flag)
           {
               printf("%20.18f\t%20.18f\t%11.9f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t",
                   uptr->pval,     fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, type2string(uptr->category >> 28));
           }
           else
           {
               printf("%11.9f\t%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t",
                   uptr->pval,     
                   uptr->pval2,     
                   fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, type2string(uptr->category >> 28));
#if 0
               printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
                   uptr->pval,     fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, type2string(uptr->category >> 28));
           printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
               uptr->pval,     fdr_for_output, uptr->enrichment_score,
               uptr->a, uptr->b, uptr->c, uptr->d,
               uptr->acc, tmps,
               uptr->name, type2string(uptr->category >> 28));
#endif
            }
            for (j=0 ; j<uptr->hitcnt ; j++)
            {
                z = egid2hugo(*((uptr->genehits)+j));
                printf("%s ",z);
            }
            printf("\n");
        }
    }
    
    fflush(NULL); 

    if (u)
    {
        for (i=0 ; i<num_used_paths ; i++)
        {
           uptr = (u+i);
           if (uptr->genehits) { free(uptr->genehits); uptr->genehits = (unsigned int *)0; }
           if (uptr->egids)    { free(uptr->egids);    uptr->egids = (unsigned int *)0; }
           if (uptr->acc)      { free(uptr->acc);      uptr->acc = (char *)0; }
           if (uptr->name)     { free(uptr->name);     uptr->name = (char *)0; }
        }
        free(u);
    }
    if (by_egids) free(by_egids);
    if (real_universe) free(real_universe);
    return status;
}
#endif


/*

Hi Rich,

Create a vector of gene names, with each name appearing as often as the number of pathways the gene appears in.  
From this list randomly select as many genes as on your DEG list.
Replace repeats, need a list of unique genes
                                                             ii.     Create the contingency table for each pathway of interest (list of top pathways)

                                                            iii.     Obtain and record its p value

                                                            iv.     For each pathway of interest keep count of how many of the 200 permutations produced more significant p values than the real data

Repeat a) etc. 200 times
The counts/200 are approximate p values for the pathways.

REFERENCE NOTES on L2P:


1.Fisher’s Exact Test is as below
 
Odds ratio = (A*D)/(B*C)

A and D are high where there’s true significance in the pathway (a good pathway).

2. Adjusted p-value (GPCC): gene pathway count correction): another column for L2P


Calculate odds ratio for each pathway….    ->    Get all Odds ratios > 1

Select from a list of “augmented set of random genes” – represent same proportions as in the total pathways  

Recreate Fisher based on 500 “random genes” -  calculate the odds ratio 1000x or more

Is the Odds Ratio the same or more extreme in the same direction.

                Count these instances where (Odds Ratio)rand  > =  (Odds Ratio)real and divide by the number of permutations (ie. random tests).

The p value is the permutation p value

3. Multiple Comparison Correction (BH pvalue): another column for L2P

Rich’s own implementation of BH in C

4. Zeroes in the L2P output


?formatC
?apply
<0.001  -> 1x10-4
0.05 <- 5x10-2

colNumber new_name (l2p output) current_MAAPster_name old_name (l2p) description

1 Pathway_Name Description pathwayname Name of pathway 
2 Category Source category KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database) *was soruce* 
3 P_Value P_Value pval fisher's 2x2 exact test 
4 FDR FDR fdr false discovery rate: default is benjamini hochberg, 
**5 Permutation_p-value Permute_pval See detailed methods 
6 Enrichment_Score Ratio ratio (A/B)/(C/D) * 100 
**7 Percent_Gene_Hits_Per_Pathway NEW FIELD (A/(A+B)) * 100 
8 Significant_Genes_IN_Pathway Number_Hits pwhitcount (A) 
9 Non-Significant_genes_IN_Pathway Number_Misses pwnohitcount (B) 
10 Significant_genes_NOT_IN_Pathway Number_User_Genes inputcount (C) 
11 Non-Significant_Genes_NOT_IN_Pathway Total_Genes_Minus_Input pwuniverseminuslist (D) 
12 Pathway_ID Pathway_ID pathwayaccessionidentifier canonical accession ( if available, otherwise assigned by us ) 
13 Pathway_Type pathwaytype functional_set, pathway, structural_complex , custom 
14 Gene_List Gene_List genes_space_separated HUGO genes from user list that hit the pathway

Table referenced from attached email:
Significant Genes Non-Significant Genes
IN Pathway A B
NOT IN Pathway C D
 
*/


