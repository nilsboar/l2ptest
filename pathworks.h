
#ifndef   PATHWORKS_H 
#define   PATHWORKS_H 1
#endif

// #define MAXBSID 25000
     // last count: 24376  pathworks.txt
     //
#define MAXBSIDPOSSIBLE 851568   
     // for pwharvest

#define MAXGENE 70000
     // last count 61141 pathworksgenes.txt 

  // types
#define type_functional_set 1
#define type_pathway 2
#define type_structural_complex 3
#define type_custom 4
#define type_unknown 9

   // scope
#define conserved_biosystem_scope 10
#define organism_specific_biosystem_scope 11



#define CAT_NCBI_BIOCYC    (1<<0)     
                                           // count = 294
#define CAT_NCBI_GO        (1<<1)     
                                           // count = 13515 
#define CAT_NCBI_KEGG      (1<<2)   
                                           // count = 485 
#define CAT_NCBI_PANTH     (1<<3)     
                                           // count = 129 
#define CAT_NCBI_PID       (1<<4)                 
                                           // count = 183 
#define CAT_NCBI_REACTOME      (1<<5)  
                                           // count = 1548 
#define CAT_NCBI_WikiPathways  (1<<6)  
                                           // count = 345 
#define CAT_MSIG_C1            (1<<7)  
                                           // count = 325 
#define CAT_MSIG_C2            (1<<8)  
                                           // count = 3777 
#define CAT_MSIG_C3            (1<<9)  
                                           // count = 836 
#define CAT_MSIG_C4            (1<<10)  
                                           // count = 858 
#define CAT_MSIG_C5            (1<<11)  
                                           // count = 5871 
#define CAT_MSIG_C6            (1<<12)  
                                           // count = 144 
#define CAT_MSIG_C7            (1<<13)  
                                           // count = 1888 
#define CAT_MSIG_H             (1<<14)  
                                           // count = 50 
#define CAT_MSIG_C8            (1<<14)  
#define CAT_CUSTOM             (1<<16)     
#if 0
#define CAT_MSIG_ARCHIVED      (1<<15)  
                                           // count = 858 
#endif
 
#define NCBI_PAT = (CAT_NCBI_BIOCYC|CAT_NCBI_GO|CAT_NCBI_KEGG|CAT_NCBI_PANTH|CAT_NCBI_PID|CAT_NCBI_REACTOME|CAT_NCBI_WikiPathways);
#define MSIG_PAT = (CAT_MSIG_C1|CAT_MSIG_C2|CAT_MSIG_C3|CAT_MSIG_C4|CAT_MSIG_C5|CAT_MSIG_C6|CAT_MSIG_C7|CAT_MSIG_C8|CAT_MSIG_H);

struct binpathouttype     // the "binary" pathway information file  
{ 
    int bsid;         // 32 bit integer. note: originally "bs" was for "biosystems"
    int category;     // bit patern for each category (particular may be bit set to turn on) examples:CAT_NCBI_GO
    int accession;    // spill from char * this is an "offset from spill space start" points to 
    int name;  // spill for char * . 
    int type; // 1 byte, from char *
    int scope; // 32 bit int  from char *
    int taxid; // 32 bit int , taxonomyid
    int desc;  // spill from char *    "points" (really offset) into spill space
    unsigned int numgenes;

// little tricky here: "hits" is not set at record creation (it is set to null), then, later (i.e when running l2p) ,  
// it is used in processing when the binpath[] data is read in. use this "hits" field for the count of 
// genes that hit this pathway
    int offset2geneids;   // pointer to "numgenes" geneids 

        // ***** NOTE: Different C compilers produce different sized records for this structure (binpathouttype).
        // ***** The output file for this will only contain the important to save fields
        // ***** we only need to writeout the above 10 fields. So output record size is 10*4=40 bytes.

#if 0
//  hits is a ptr to an array with 
//    firstelemen=[0]=numhits, then the rest of the array is [1...n] ptrs to struct of generecs (bingentype? right?)
    void *hits;   // used in  l2p for user "hits" to this pathway
#endif
};

#define MAXGENENAME 26
// maximum length ARHGAP27P1-BPTFP1-KPNA2P3 = 25

struct bingenetype
{
    int geneid; // entrez gene id
    char hugo[MAXGENENAME];
    char ensembl[MAXGENENAME];
    int pathcount;   // count of paths, ids are in int array famous at "pathplace"
    int pathplace;   // index to path (to a struct binpathouttype record,see above)
    int categories;  // bit patterns 
};

struct genelisttype // used by harvest programs 
{
    int geneid;
    struct genelisttype *n;
};

struct bstype // biosystems id and info - input into this array  -- used by harvest programs 
{
    int bsid;
    int  category; // use CAT_ bitpattern defines (above) 
    char *accession;
    char *name;
    char *type;
    char *scope;
    int taxid;
    char *desc;
    int redundant;       // flag for checking to see if this pathway is duplicated by another pathway 
// next two fields get values from other file
    int numgenes;       // "count of" in next line of code line (i.e. number of genes)
    struct genelisttype *geneslinkedlist; // a linked list
};


struct hugo_type
{
    char *hugo;
    struct bingenetype *generec_ptr; 
    int status;   // this can be used for various purposes, initial reason is to use for "universe" masking 
};

struct genetype // from ncbi
{                              /// this is only used in pwharvest, l2p uses bingenetype
    int geneid;
    char *hugo;
    char *ensembl;
    int categories;
};

struct hit_type
{
    unsigned int hitcnt;
    unsigned int maxhits;
    unsigned int *hitsindexes;
};


// pathway commons
#define chemical_affects                             1
#define in_complex_with                              2
#define catalysis_precedes                           4
#define controls_expression_of                       8
#define controls_state_change_of                     16
#define controls_production_of                       32
#define consumption_controlled_by                    64
#define controls_phosphorylation_of                  128
#define used_to_produce                              256
#define transport                                    512
#define reacts_with                                  1024
#define interacts_with                               2048
#define reference                                    4096
#define other                                        8192


#define MAXPC 1915769
        // latest 1915769 PathwayCommons12.All.hgnc.txt

struct pctype
{
    int ID_Interactor_A;
    int ID_Interactor_B;
    int interaction_type;
};

#define MAXBIOGRID 303568

// bits for "interaction_type" field ...                         // count name  
#define association                                           1  // 8931 psi-mi:"MI:0914(association)"
#define colocalization                                        2  // 44101 psi-mi:"MI:0403(colocalization)"
#define synthetic_genetic_interaction_defined_by_inequality   4  // 50045 psi-mi:"MI:0794(synthetic genetic interaction defined by inequality)"
#define suppressive_genetic_interaction_defined_by_inequality 8  // 197811 psi-mi:"MI:0796(suppressive genetic interaction defined by inequality)"
#define direct_interaction                                    16 // 206875 psi-mi:"MI:0407(direct interaction)"
#define physical_association                                  32 // 329721 psi-mi:"MI:0915(physical association)"
#define additive_genetic_interaction_defined_by_inequality    64 // 535593 psi-mi:"MI:0799(additive genetic interaction defined by inequality)"
struct biogridtype
{
    int ID_Interactor_A;
    int ID_Interactor_B;
    int interaction_type;
};

struct smallgenetype
{
    char *hugo;            // hugo = human gene name nomenclature authority  ("official gene name")
    unsigned int egid;              // entrez gene id
};
struct tree_with_count
{
    unsigned int val; // entrez gene id
    unsigned int count; // times
    unsigned int deg; // 1 on deglist, 0 not on
    struct tree_with_count *left;
    struct tree_with_count *right;
    struct used_path_type **all_gene_paths; // all gene paths is an array of pointers
    unsigned int pathindex; // which array member gets the pointer to pathway?
};
struct used_path_type
{
    unsigned int category;
    char *acc;
    char *name;
    unsigned int numgenes;
    unsigned int numfixedgenes;
    unsigned int hitcnt;
    double gpcountsum;
    unsigned int *genehits; // put hits here, reason: need to print them out
    double OR;
    double pval;
    double pval2; // alt
    double pval3; // permute
    double pval4; // gp count corrected FET
    double fdr;
    int rand_paths_more_sig;
    double enrichment_score;                  // ratio
    unsigned int pwgenesindex;
    unsigned int *egids;
    // orginal george int a,b,c,d;  // a=universe-userinput-pwgenes-list b=pw-hits, c=degs-hits , d = number of hits
    unsigned int a,b,c,d;  // a=universe-userinput-pwgenes-list b=pw-hits, c=degs-hits , d = number of hits
    unsigned int randhits;
    unsigned int countover; // data hits value > permutation p hits
    unsigned int countequal;
    unsigned int countunder; // redundant
    double p_permute_over;
    double p_permute_under; // redundant
};

struct custom_type
{
    char *name;
    char *optional;
    unsigned int numgenes;
    unsigned int *genes;
};

struct ens2gene_type
{
    char *ens;
    char *symbol;
};

struct a2a_type {
   int taxid1; 
   int taxid2; 
   int ensidx1; 
   int ensidx2; 
};

struct synonym_type {
     char *Synonym;
     int GeneID;
     char *Symbol;
     int status;
};



void category_set_all(unsigned int *pat);
void category_code_to_string(unsigned int cat,char puthere[]);
int string_to_category_code(char cats[]);
void categories_pattern_to_strings(unsigned int cat,char puthere[]);
double exact22(int n11_,int n12_,int n21_,int n22_);    // fishers exact 
double exact22_oneside(int n11_,int n12_,int n21_,int n22_, int dbg);
unsigned int string2type(char *s);
int bitCount(int n);
int setup_by_egids(void);
char *egid2hugo(int egid);
unsigned int hugo2egid(char *h);
char *type2string(int type);
int cmp_ui(const void *a, const void *b);
unsigned int *get_used_universe(struct used_path_type *u, unsigned int num_used, unsigned int *real_universe_cnt);
int cmp_ordertype_by_val_REV(const void *a, const void *b);
int cmp_usi(const void *a, const void *b);
struct used_path_type *setup_used_paths(unsigned int *num_used_paths, unsigned int catspat, char universe_file[], unsigned int in_universe_cnt,unsigned int *in_universe, char custom_file[], unsigned int *real_universe_cnt_ptr,unsigned int **real_universe,unsigned int lencust,struct custom_type *mycustompw);
int do_pvals_and_bh(unsigned int ingenecnt, struct used_path_type usedpaths[], unsigned int numusedpaths,unsigned int real_universe_cnt, int oneside);
unsigned int GPCC(struct used_path_type usedpaths[], unsigned int num_used_paths, unsigned int real_universe_cnt, unsigned int *real_universe, int seed);
int do_just_bh(unsigned int ingenecnt, struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt);
// void malloc_pathpointers(struct tree_with_count *node); // counts aligned with universe (real_universe)
