/*****************************************************************
 *                                                               *
 *   score.h                                                     *
 *                                                               *
 *****************************************************************
 *                                                               *
 *   written by JR, modified by Yoginho                          *
 *                                                               *
 *****************************************************************
 *                                                               *
 * This header file contains stuff that is needed for reading    *
 * and defining facts and limits and for scoring functions. The  *
 * functions declared here initialize or manipulate facts or     *
 * data time tables, read and initialize limits and penalty (if  *
 * needed and do the actual scoring.                             *
 *                                                               *
 * Should be included in code that does scoring (e.g. printscore *
 * or annealing code).                                           *
 *                                                               *
 *****************************************************************/


#ifndef SCORE_INCLUDED
#define SCORE_INCLUDED

/* this def needed for old func. defs that refer to (* FILE) */
#ifndef _STDIO_INCLUDED
#include <stdio.h>    
#endif

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include <global.h>
#endif

/* maternal.h should ALWAYS be included when score.h is used */
#ifndef MATERNAL_INCLUDED
#include <maternal.h>
#endif



/*** STRUCTURES ************************************************************/

/* range struct for limits and such */

typedef struct Range {       
  double   lower;    
  double   upper;    
} Range;

/***************************************************************************
 * The following describes the search space to the scoring function.       *
 * There are two ways to specify limits. Lambda, R, and d always get a     *
 * range for each variable---an upper & lower limit. Elements of T, h      *
 * and m that contribute to u in g(u) can be given a range.                *
 * This is probably the way to go as it definitely results in an           *
 * ergodic search space. However, as I write this code (10/95) we have     *
 * only one set of runs using this method. The alternative, which has      *
 * been mostly used, is to treat T, h and m with a penalty function on     *
 * u of the form:                                                          *
 *                                                                         *
 *           | 0 if exp(Lambda*\sum((T^ij)v_max_j)^2 + h^2 +               *
 *           |                                     (m_i mmax)^2 - 1 < 0    *
 * penalty = |                                                             *
 *           | exp(Lambda*\sum((T^ij)v_max_j)^2 + h^2 +(m_i mmax)^2        *
 *           |                                                otherwise    *
 *                                                                         *
 * where vmax_i is the highest level of gene i in the data, similarly      *
 * for mmax and bcd. This method can be non-ergodic if h and T or m are    *
 * being altered and variables change one at a time. In any case, one      *
 * scheme or the other must be used and the programmer must insure that    *
 * one or another set of pointers (see below) are NULL.                    *
 *                                                                         *
 * NOTE: - Lambda is NOT the same stuff as lambda in equation params or    *
 *         the lambda of the Lam schedule (don't get confused!)            *
 ***************************************************************************/

typedef struct {
  double    *pen_vec;    /* pointer to array that defines penalty function */
  Range     **Rlim;                      /* limits fore promoter strengths */
  Range     **Tlim;        /* limits for T matrix, NULL if pen_vec != NULL */
  Range     **Elim;        /* limits for E matrix, NULL if pen_vec != NULL */
  Range     **mlim;               /* limits for m, NULL if pen_vec != NULL */
  Range     **hlim;               /* limits for h, NULL if pen_vec != NULL */
  Range     **dlim;                 /* limit(s) for diffusion parameter(s) */
  Range     **lambdalim;           /* limits for lambda (prot. decay rate) */
  Range     **taulim;           /* limits for tau (delays) */
} SearchSpace;

/* Yousong's GutInfo struct for writing square diff guts */

typedef struct GutInfo{
  int      flag;                    /* for setting the gut flag in score.c */
  int      ndigits;                                /* gut output precision */
} GutInfo;



/* FUNCTION PROTOTYPES *****************************************************/

/* Initialization Functions */
 
/*** InitScoring: intializes a) facts-related structs and TTimes and *******
 *                           b) parameter ranges for the Score function.   *
 ***************************************************************************/ 

void        InitScoring   (FILE *fp);

/*** InitFacts: puts facts records into the appropriate DataTable. *********
 *              Returns a pointer to a DataTable, which in turn points to  *
 *              a sized array of DataRecords, one for each time step.      *
 ***************************************************************************/

void        InitFacts     (FILE *fp);

/*** InitLimits: reads limits section from the data file into the struct ***
 *               limits, which is static to score.c. Then, it initiali-    *
 *               zes the penalty function if necessary.                    *
 *   NOTE:       lambda limits are stored as protein half lives in the     *
 *               data file and therefore get converted upon reading        *
 *               (compare to ReadParameters in zygotic.c for EqParms)      *
 ***************************************************************************/

void        InitLimits    (FILE *fp);      

/*** InitPenalty: initializes vmax[] and mmax static to score.c; these *****
 *                variables are used to calculate the penalty function     *
 *                for scoring.                                             *
 *         NOTE: Should only get called, if penalty function is used       *
 ***************************************************************************/

void        InitPenalty   (FILE *fp);

/*** InitTTs: Initializes the time points for which we need model output ***
 *            i.e. the time points for which we have data.                 *
 *      NOTE: we do not allow data at t=0 but it still gets included into  *
 *            TTs                                                          *
 ***************************************************************************/

void        InitTTs       (void);

/*** InitStepsize: the only thing this function does is making step *******
 *                 static to score.c so the Score function can use it.    *
 **************************************************************************/

void InitStepsize(double step, double accuracy, FILE *slog, char* infile, char* cost_function);





/* Actual Scoring Functions */

/*** Score: as the name says, score runs the simulation, gets a solution   *
 *          and then compares it to the data using the Eval least squares  *
 *          function.                                                      *
 *   NOTE:  both InitZygote and InitScoring have to be called first!       *
 ***************************************************************************/

double      Score         (void);     

/*** Eval: scores the summed squared differences between equation solution *
 *         and data. Because the times for states written to the Solution  *
 *         structure are read out of the data file itself, we do not check *
 *         for consistency of times in this function---all times with data *
 *         will be in the table, but the table may also contain additional *
 *         times.                                                          *
 ***************************************************************************/

double      Eval          (NArrPtr Solution, int gindex);
 




/*** Scoregut functions */

/* SetGuts: sets the gut info in score.c for printing out guts *************
 ***************************************************************************/

void SetGuts (int srsflag, int ndigits);

/*** GutEval: this is the same as Eval, i.e it calculates the summed squa- *
 *            red differences between equation solution and data, with the *
 *            addition that individual squared differences between data-   *
 *            points are written to STDOUT in the unfold output format     *
 ***************************************************************************/

double GutEval(NArrPtr Solution, int gindex);





/* Functions that return info about score.c-specific stuff */

/*** GetTTimes: this function returns the times for which there's data *****
 *              for a given genotype                                       *
 *     CAUTION: InitTTs has to be called first!                            *
 ***************************************************************************/

DArrPtr     GetTTimes     (char *genotype);     

/*** GetLimits:  returns a pointer to the static limits struct in score.c **
 *     CAUTION:  InitScoring must be called first!                         *
 ***************************************************************************/

SearchSpace *GetLimits    (void);

/*** GetPenalty: calculates penalty from static limits, vmax and mmax ******
 *      CAUTION: InitPenalty must be called first!                         *
 ***************************************************************************/

double      GetPenalty    (void);

/*** GetLimits:  returns the number of data points (static to score.c) *****
 *     CAUTION:  InitScoring must be called first!                         *
 ***************************************************************************/

int         GetNDatapoints(void);




/* A function for converting penalty to explicit limits */

/*** Penalty2Limits: uses the inverse function of g(u) to calculate upper **
 *                   and lower limits for T, m and h using the penalty     *
 *                   lambda parameter; these limits can only be an appro-  *
 *                   ximation, since all T, m and h are added up to yield  *
 *                   u for g(u); we try to compensate for this summation   *
 *                   dividing the limits by sqrt(n); this function then    *
 *                   sets the penalty vector to NULL and supplies explicit *
 *                   limits for T, m and h, which can be used for scram-   *
 *                   bling parameters and such                             *
 ***************************************************************************/

void Penalty2Limits(SearchSpace *limits);




/* Functions used to read stuff into structs */

/*** ReadLimits: reads the limits section of a data file and returns the  **
 *               approriate SearchSpace struct to the calling function     *
 ***************************************************************************/

SearchSpace *ReadLimits   (FILE *fp);

/*** List2Facts: takes a Dlist and returns the corresponding DataTable *****
 *               structure we use for facts data.                          *
 ***************************************************************************/

DataTable   *List2Facts   (Dlist *inlist);

void FreeFacts(DataTable *D);
void InitHistory(FILE *fp);
void InitExternalInputs(FILE *fp);
void GetInterp(FILE *fp, char *title, int num_genes, 
										DataTable **interp_tables);
void FreeHistory(void);
void FreeExternalInputs(void);

#endif






