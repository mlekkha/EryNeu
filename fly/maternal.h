/*****************************************************************
 *                                                               *
 *   maternal.h                                                  *
 *                                                               *
 *****************************************************************
 *                                                               *
 *   written by JR, modified by Yoginho                          *
 *                                                               *
 *****************************************************************
 *                                                               *
 * This header contains several kinds of things:                 *
 *                                                               *
 * 1. structs and constants for things used throughout the fly   *
 *    model specific part of the code (since this file will al-  *
 *    ways get included when we run the fly model); therefore    *
 *    all fly-specific things that used to be in global.h are    *
 *    now in here                                                *
 * 2. stuff that deals with that part of the blastoderm which is *
 *    fixed by the maternal genotype, i.e. division schedules    *
 *    (including the rules for zygotic.c, the number of nuclei   *
 *    at each cleavage cycle and the lineage number of the most  *
 *    anterior nucleus for each cleavage cycle), bicoid gra-     *
 *    dients and diffusion schedules                             *
 * 3. bias-related stuff is also in here since it's needed to    *
 *    run the model                                              *
 *                                                               *
 *****************************************************************/

#ifndef MATERNAL_INCLUDED
#define MATERNAL_INCLUDED 

/* this def needed for func. defs that refer to (* FILE) *******************/
#ifndef _STDIO_INCLUDED
#include <stdio.h>    
#endif

/* following for structures & consts used thruout **************************/

#ifndef GLOBAL_INCLUDED
#include <global.h>
#endif



/*** CONSTANTS *************************************************************/

/* The following is the maximum stepsize for solvers (now set to the whole */
/* time it takes to model to run (100.8 minutes))                          */

#define MAX_STEPSIZE        100.8

/* masks for lineage numbers */

#define     CYCLE1            1
#define     CYCLE2            2
#define     CYCLE3            4
#define     CYCLE4            8
#define     CYCLE5           16
#define     CYCLE6           32
#define     CYCLE7           64
#define     CYCLE8          128
#define     CYCLE9          256
#define     CYCLE10         512
#define     CYCLE11        1024
#define     CYCLE12        2048
#define     CYCLE13        4096
#define     CYCLE14        8192

/* define total number of divisions including the ones before cycle 10
*/

#define TOTAL_DIVS 6

/* Data points that have a concentration of -1 will be ignored! */
/* this needs to be here rather than in score.h, since it's used to */
/* calculate the number of data points in ReadData */

#define     IGNORE          -1.  



/*** STRUCTS ***************************************************************/

/* general struct used for sized array of doubles */

typedef struct DArrPtr {
  int      size;     
  double   *array;  
} DArrPtr;

/* this is the problem at hand */

typedef struct TheProblem {
  int    ngenes;                            /* number of genes to consider */
  int    egenes;                            /* number of external inputs */
  char   *gene_ids;                       /* pointer to char with gene IDs */
  char   *egene_ids;                       /* pointer to char with external gene IDs */
  int    ndivs;                                /* number of cell divisions */
  int    nnucs;                        /* number of nuclei at gastrulation */
  char   diff_schedule;                              /* diffusion schedule */
} TheProblem;

/* two structs for bicoid gradients */

typedef struct BcdGrad {
  int      ccycle;                                   /* the cleavage cycle */
  DArrPtr  gradient;                  /* the array of concs for a gradient */
} BcdGrad;

typedef struct BArrPtr {      /* points to an array of BcdStates of 'size' */
  int      size;                   
  BcdGrad  *array;           
} BArrPtr;

/* two structures for bias data and Blastoderm() output (solution) */

typedef struct NucState {
  double   time;                                               /* the time */
  DArrPtr  state;                /* the array of v's, and how many of them */
} NucState; 

typedef struct NArrPtr {
  int      size;                   /* How many NucState elements in array? */
  NucState *array;    /* points to 1st element of array of NucState elems. */
} NArrPtr; 

/* The following three structs are used to store facts data. ***************
 * The idea is to encode sparse data against non-sparse v[index] state     *
 * arrays. This will in most of the present cases use more memory, but     *
 * that won't always be true. Also, it is more convenient for those        *
 *'don't care' points.                                                     *
 * NOTE: needs to be in this generic header for union definition below     *
 ***************************************************************************/

typedef struct DataPoint {
  int        index;
  double     conc;
} DataPoint;

typedef struct DataRecord {
  int        size;
  double     time;
  DataPoint  *array;
} DataRecord;

typedef struct DataTable {
  int        size;
  DataRecord *record;
} DataTable;

/* struct and union for genotype number and pointer to data ****************/
/* used for bicoid, bias, facts and time tables                            */

typedef union DataPtr {
  DArrPtr          times;                   /* used for bias and tab times */
  BArrPtr          bicoid;                  /* for bicoid stuff            */
  NArrPtr          bias;                    /* for the bias                */
  DataTable        *facts;                  /* for facts                   */
} DataPtr;

typedef struct GenoType {      
  char             *genotype;         
  DataPtr          ptr;
} GenoType;  
 
/* linked lists for reading bicoid (Blist) and bias/facts (Dlist) */

typedef struct Blist {                  /* linked list used to read bicoid */
  unsigned int lineage;
  double       conc;
  struct Blist *next;
} Blist;


typedef struct Dlist {            /* linked list used to read bias & facts */
  unsigned int lineage;
  double       *d;
  struct Dlist *next;
} Dlist;


/* linked list used to read in genotype specific sections from datafile */

typedef struct Slist {
  char         *genotype;           /* genotype string (see dataformatX.X) */
  char         *bcd_section;                          /* bcd section title */
  char         *bias_section;                        /* bias section title */
  char         *fact_section;                        /* fact section title */
  char         *hist_section;                        /* fact section title */
  char         *ext_section;                        /* fact section title */
  struct Slist *next;             
} Slist;





/*** GLOBALS ** ************************************************************/

int    debug;                                            /* debugging flag */
int    olddivstyle;                    /* flag: old or new division times? */
int    nalleles;                  /* number of alleles (genotypes) in data */
                               
double maxconc;                /* max prot conc: 12 (old-), 255 (newstyle) */
double custom_gast;                  /* custom gastrulation time set by -S */ 


/***************************************************************************
 * - defs:     is the global problem struct; must be visible to scoring    *
 *             functions, to Blastoderm in integrate.c, to stuff in mater- *
 *             nal.c and to functions in zygotic.c (need any more justifi- *
 *             cation for this to be global?)                              *
 ***************************************************************************/

TheProblem  defs;





/*** FUNCTION PROTOTYPES ***************************************************/

/* Initialization Functions */

/*** InitBicoid: Copies the Blist read by ReadBicoid into the DArrPtr ******
 *               structure. The bicoid DArrPtr contains pointers to bcd    *
 *               arrays for each cell division cycle. Also initializes the *
 *               lin_start array for keeping track of lineage numbers.     *
 ***************************************************************************/

void         InitBicoid(FILE *fp);   

/*** InitBias:  puts bias records in a form where get_bias can use them. ***
 *              It expects times in increasing order. It expects a non-    *
 *              sparse entry, with no genes or nuclei missing.             *
 ***************************************************************************/

void         InitBias(FILE *fp);       

/*** InitBTs: initializes the static BT struct that holds all times for ****
 *            which we have bias.                                          *
 ***************************************************************************/

void         InitBTs(void);

/*** InitNNucs: takes the global defs.nnucs and calculates number of nucs **
 *              for each cleavage cycle which are then stored in reverse   *
 *              order in the static nnucs[] array                          *
 *   CAUTION:   defs struct needs to be initialized before this!           *
 ***************************************************************************/

void         InitNNucs(void);





/* Functions that return info about embryo */

/*** GetBicoid: returns a bicoid gradients (in form of a DArrPtr) for a ****
 *              specific time and genotype.                                *
 ***************************************************************************/

DArrPtr      GetBicoid(double time, int genindex);

/*** GetBias: This function returns bias values for a given time and *******
 *            genotype.                                                    *
 ***************************************************************************/

DArrPtr      GetBias(double time, int genindex);
 
/*** GetBTimes: returns times for which there's bias ***********************
 ***************************************************************************/

DArrPtr      GetBTimes(char *genotype);

/*** GetNNucs: returns number of nuclei for a given time *******************
 ***************************************************************************/

int          GetNNucs(double t);

/*** GetStartLin: returns the lineage number of the most anterior nucleus **
 *                for a given time                                         *
 ***************************************************************************/

int          GetStartLin(double t);

/*** GetCCycle: returns cleavage cycle for a given time ********************
 ***************************************************************************/

unsigned int GetCCycle(double time);
                        
/*** ParseLineage: takes lineage number as input and returns the cleavage **
 *                 cycle the nucleus belongs to.                           *
 ***************************************************************************/

unsigned int ParseLineage(unsigned int lin);
 
/*** GetDivtable: returns times of cell divisions depending on ndivs and ***
 *                olddivstyle; returns NULL in case of an error            *
 ***************************************************************************/

double *GetDivtable(void);

/*** GetDurations: returns pointer to durations of cell divisions de- ******
 *                 pending on ndivs and olddivstyle; returns NULL in case  *
 *                 of an error                                             *
 ***************************************************************************/

double *GetDurations(void);

/*** GetGastTime: returns time of gastrulation depending on ndivs and ******
 *                olddivstyle; returns 0 in case of an error               *
 ***************************************************************************/

double GetGastTime(void);

/*** GetD: returns diffusion parameters D according to the diff. params. ***
 *         in the data file and the diffusion schedule used.               *
 *   NOTE: Caller must allocate D_tab                                      *
 ***************************************************************************/

void         GetD(double t, double *d, char diff_schedule, double *D_tab);

/*** GetRule: returns the appropriate rule for a given time; used by the ***
 *            derivative function                                          *
 ***************************************************************************/

int          GetRule(double time);

/*** Theta: the value of theta(t) used by DvdtOrig ***
 *            for autonomous eqns.                                          *
 ***************************************************************************/

int          Theta(double time);

/*** GetIndex: this functions returns the genotype index for a given *******
 *             genotype number for reading the GenoType struct.            *
 ***************************************************************************/

int          GetIndex(char *genotype);
                               
/*** MakeTable: Returns the appropriate time table from uftimes.c **********
 ***************************************************************************/

DArrPtr      MakeTable(double p_stepsize);





/* Following functions read data from file into structs for bias and bcd */

/*** ReadTheProblem: reads the problem section of a data file into the *****
 *                   TheProblem struct.                                    *
 ***************************************************************************/

TheProblem   ReadTheProblem(FILE *fp);

/*** ReadGenotypes: This function reads all the genotypes in a datafile & **
 *                  returns an SList with genotype number and pointers to  *
 *                  the corresponding section titles for bias, bcd & facts *
 ***************************************************************************/

Slist        *ReadGenotypes(FILE *fp);

/*** ReadBicoid: Reads the bcd section of a data file into a linked list, **
 ***************************************************************************/

Blist        *ReadBicoid(FILE *fp, char *section); 

/*** ReadData: reads in a data or bias section and puts it in a linked *****
 *             list of arrays, one line per array; ndp is used to count    *
 *             the number of data points in a data file (ndp), which is    *
 *             used to calculate the root mean square (RMS) if required    *
 *                                                                         * 
 *             ReadData allows comments that start with a letter or punc-  *
 *             tuation mark, and treats lines starting with a number or    *
 *             .num or -.num as data. It expects to read an int and        *
 *             ngenes + 1 data points: lineage, time and ngenes protein    *
 *             concentrations. It expects data to be in increasing spatial *
 *             and temporal order.                                         *
 ***************************************************************************/

Dlist        *ReadData(FILE *fp, char *section, int *ndp);       
                
/*** ReadTimes: reads a time table from a file and returns a DArrPtr *******
 * FILE FORMAT: one time per line separated with newlines                  *
 *        NOTE: max. times tab size is 1000                                *
 ***************************************************************************/

DArrPtr      ReadTimes(char *timefile);

/*** List2Bicoid: takes a Blist and returns the corresponding BArrPtr ******
 *                structure.                                               *
 ***************************************************************************/

BArrPtr      List2Bicoid(Blist *inlist);

/*** ReadGuts: reads the $gutsdefs section in a data file into an array ****
 *             of strings which then need to get parsed                    *
 ***************************************************************************/

char         **ReadGuts(FILE *fp);

/*** List2Bias: takes a Dlist and returns the corresponding DArrPtr ********
 *              structure.                                                 *
 ***************************************************************************/
       
NArrPtr      List2Bias(Dlist *inlist);





/* Following are functions that are needed by many other reading funcs */

/*** FindSection: This function finds a given section of the input file & **
 *                returns a pointer positioned to the first record of that *
 *                section. Section titles should be passed without the pre-*
 *                ceding '$'. If it can't find the right section, the func-*
 *                tion returns NULL.                                       *
 ***************************************************************************/

FILE        *FindSection(FILE *fp, char *section);

/*** KillSection: erases the section with 'title' from the 'filename' ******
 ***************************************************************************/

void        KillSection(char *filename, char *title);






/* Following functions are utility functions for different linked lists ****
 * which are used to read in data of unknown size. All utility functions   *
 * follow the same scheme (X stands for the different letters below):      *
 *                                                                         *
 * - init_Xlist:         allocates first element and returns pointer to it *
 * - addto_Xlist:        adds the adduct to an already existing linkd list *
 * - free_Xlist:         frees memory of linked list                       *
 * - count_Xlist:        counts elements in the linked list                *
 *                                                                         *
 ***************************************************************************/

/* Utility functions for Blist (for reading bicoid) */

Blist        *init_Blist    (void);
Blist        *addto_Blist   (Blist *start, Blist *adduct);
void         free_Blist     (Blist *start);
int          count_Blist    (Blist *start);


/* Utility functions for Dlist (for reading bias and facts) */

Dlist 		 *init_Dlist	(int size);  
Dlist        *addto_Dlist   (Dlist *start, Dlist *adduct);
void         free_Dlist     (Dlist *start);
int          count_Dlist    (Dlist *start);
 

/* Utility functions for Slist (for reading genotypes ) */

Slist        *init_Slist    (void);
Slist        *addto_Slist   (Slist *start, Slist *adduct);
void         free_Slist     (Slist *start);
int          count_Slist    (Slist *start);

void InitFullNNucs(void);
int Index2StartLin(int index);
int Index2NNuc(int index);
int GetStartLinIndex(double t);
double *Get_Theta_Discons(int *theta_discon_size);
DataTable *List2Interp(Dlist *inlist, int num_genes);
int descend(int *x, int *y);
Dlist *ReadInterpData(FILE *fp , char *section, int num_genes, 
												int *ndp);

#endif
