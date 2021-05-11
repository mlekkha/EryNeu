/*****************************************************************
 *                                                               *
 *   zygotic.h                                                   *
 *                                                               *
 *****************************************************************
 *                                                               *
 *   written by JR, modified by Yoginho and Manu                 *
 *                                                               *
 *****************************************************************
 *                                                               *
 * This file is for functions that deal with the right hand side *
 * of the ODEs. We have an initializing function called          *
 * InitZygote that reads parameters, defs data and search space  *
 * limits from a data file and installs a couple of things like  *
 * the solver and temporary arrays for the dvdt function.        *
 * Then we have dvdt_orig, which is the inner loop of the model. *
 * It propagates the equations from one time point to another    *
 * and gets called whenever Blastoderm is in PROPAGATE mode.     *
 * Last but not least, we have a few mutator functions in this   *
 * file. They change the local lparm EqParm struct.              *
 *                                                               *
 *****************************************************************/
                                                          
#ifndef ZYGOTIC_INCLUDED
#define ZYGOTIC_INCLUDED

/* this def needed for func. defs that refer to (* FILE) */
#ifndef _STDIO_INCLUDED
#include <stdio.h>    
#endif

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include <global.h>
#endif

/* following for InterpObject */
#ifndef INTEGRATE_INCLUDED
#include <integrate.h>
#endif

#include <time.h>
#include <sys/resource.h>


/*** CONSTANTS *************************************************************/

/* these are the propagation rules for dvdt_orig */

#define INTERPHASE              0
#define MITOSIS                 1



/*** AN ENUM ***************************************************************/

/* This is the g(u)-function enum which describes the different types of   *
 * g(u) functions we can use in derivative functions                       */

typedef enum GFunc {
  Sqrt,
  Tanh,
  Exp,
  Hvs
} GFunc;



/*** A GLOBAL **************************************************************/

GFunc    gofu;                            /* the g(u) function we're using */



/*** STRUCTS ***************************************************************
 *   The following two structs are newstyle structs that define the pro-   *
 *   blem at hand (TheProblem) and hold the equation parameters (EqParms). *
 ***************************************************************************/

typedef struct {
  double *R;                    /* strength of each promoter--always >= 0. */
  double *T;                            /* the genetic interconnect matrix */
  double *E;                            /* the external input
  regulatory matrix */
  double *m;                      /* regulatory coefficient of bcd on gene */
  double *h;           /* reg. coeff. for generic TFs on synthesis of gene */
  double *d;          /* spatial interaction at gastrulation--always >= 0. */
  double *lambda;                       /*protein half lives--always >= 0. */
  double *tau;			/* delay times for the proteins */
} EqParms;



/*** FUNCTION PROTOTYPES ***************************************************/

/* Initialization Functions */

/*** InitZygote: makes pm and pd visible to all functions in zygotic.c and *
 *               reads EqParms and TheProblem. It then initializes bicoid  *
 *               and bias (including BTimes) in maternal.c. Lastly, it     *
 *               allocates memory for structures used by the derivative    *
 *               function DvdtOrig that are static to zygotic.c.           *
 ***************************************************************************/

void InitZygote(FILE *fp, void (*pd)(), void (*pj)(), char* parm_section);

/*** InitGenotype: function used to make genotype static to zygotic.c ******
 *                 since DvdtOrig later needs to read it for getting bcd   *
 ***************************************************************************/

void InitGenotype(int g_index); 





/* Cleanup functions */

/*** FreeZygote: frees memory for D, vinput, bot2 and bot arrays ***********
 ***************************************************************************/

void FreeZygote(void);

/*** FreeMutant: frees mutated parameter struct ****************************
 ***************************************************************************/

void FreeMutant(void);





/* Derivative Function(s) */

/* Derivative functions calculate the derivatives for the solver. **********/

/*** DvdtOrig: the original derivative function; implements the equations **
 *             as published in Reinitz & Sharp (1995), Mech Dev 49, 133-58 *
 *             plus different g(u) functions as used by Yousong Wang in    *
 *             spring 2002.                                                *
 ***************************************************************************/

void DvdtOrig(double *v, double t, double *vdot, int n);


void DvdtDelay(double *v, double **vd, double t, double *vdot, int n);



/* Jacobian Function(s) */

/* Calculate the Jacobian for a given model at a given time; these funcs ***
 * are used by certain implicit solvers                                    */

/*** JacobnOrig: Jacobian function for the DvdtOrig model; calculates the **
 *               Jacobian matrix (matrix of partial derivatives) for the   *
 *               equations at a give time t; input concentrations come in  *
 *               v (of size n), the Jacobian is returned in jac; note that *
 *               all dfdt's are zero in our case since our equations are   *
 *               autonomous (i.e. have no explicit t in them)              *
 ***************************************************************************/

void JacobnOrig(double t, double *v, double *dfdt, double **jac, int n);





/*** GUTS FUNCTIONS ********************************************************/

/*** CalcGuts: calculates guts for genotpye 'gtype' using unfold output in *
 *             'table' and the guts string 'gutsdefs'; it returns the num- *
 *             ber of columns we'll need to print and the guts table in    *
 *             'gtable'                                                    *
 ***************************************************************************/

int CalcGuts(int gindex, char *gtype, 
				InterpObject hist_interrp, InterpObject extinp_interrp, 
				NArrPtr table, NArrPtr *gtable, char *gutsdefs);

/*** CalcRhs: calculates the components of the right-hand-side (RHS) of ****
 *            the equatiion which make up guts (e.g. regulatory contribu-  *
 *            tions of specific genes, diffusion or protein decay); this   *
 *            function makes use of the derivative function and calculates *
 *            everything outside g(u) in reverse, i.e. starting from the   *
 *            derivative and calculating the desired gut properties back-  *
 *            wards from there; everything within g(u) is simply recon-    *
 *            structed from concentrations and parameters                  *
 ***************************************************************************/

void CalcRhs(double *v, double t, double *guts, int n, int gn, 
	     int numguts, int which, unsigned long *gutcomps);

/*** ParseString: parses a line of the $gutsdefs section (dataformatX.X ****
 *                has the details); this function then returns two things: *
 *                1) the number of 'words' separated by whitespace in that *
 *                   string is the function's return value                 *
 *                2) the 'words' themselves are returned in arginp         *
 ***************************************************************************/

int ParseString(char *v, char **arginp);

/*** GetGutsComps: takes a geneID string and an array of ID strings; then **
 *                 returns an array of long ints with bit flags set that   *
 *                 tell CalcRhs() what to calculate for each column of     *
 *                 guts output; this function also returns a pointer to    *
 *                 the current gene in the gene ID string                  *
 ***************************************************************************/

char *GetGutsComps(char *geneidstring, char *egeneidstring, 
					char **specsin, unsigned long *specsout);


/* Mutator functions */

/*** Mutate: calls mutator functions according to genotype string **********
 ***************************************************************************/

void Mutate(char *g_type);

/*** T_Mutate: mutates genes by setting all their T matrix entries *********
 *             to zero. Used to simulate mutants that express a            *
 *             non-functional protein.                                     *
 ***************************************************************************/

void T_Mutate(int g_type);

/*** R_Mutate: mutates genes by setting their promotor strength R **********
 *             to zero, so there will be no transcription at all           *
 *             anymore. Used to simulate mutants that don't pro-           *
 *             duce any protein anymore.                                   *
 ***************************************************************************/

void R_Mutate(int g_type);

/*** RT_Mutate: mutates gene by setting both promoter strength R ***********
 *              and T matrix entries to zero. Useful, if you want          *
 *              to be really sure that there is NO protein pro-            *
 *              duction and that maternal contribution also don't          *
 *              contribute anything to gene interactions.                  *
 ***************************************************************************/

void RT_Mutate(int g_type);

/*** CopyParm: copies all the parameters into the lparm struct *************
 ***************************************************************************/

EqParms CopyParm(EqParms orig_parm);






/* A function that sets static stuff in zygotic.c */

/*** SetRule: sets the static variable rule to MITOSIS or INTERPHASE *******
 ***************************************************************************/

void SetRule(int r);


/* A function that return static stuff from zygotic.c */

/*** GetParameters: returns the parm struct to the caller; note that this **
 *                  function returns the ORIGINAL PARAMETERS as they are   *
 *                  in the data file and NOT THE MUTATED ONES; this is im- *
 *                  portant to prevent limit violations in Score()         *
 ***************************************************************************/

EqParms *GetParameters(void);

/*** GetMutParameters: same as above but returns the mutated copy of the ***
 *                     parameter struct; important for writing guts        *
 ***************************************************************************/

EqParms *GetMutParameters(void);




/* A function that reads EqParms from the data file */

/*** ReadParamters: reads the parameters for a simulation run from the *****
 *                  eqparms or input section of the data file as indicated *
 *                  by the section_title argument and does the conversion  *
 *                  of protein half lives into lambda parameters.          *
 ***************************************************************************/

EqParms ReadParameters(FILE *fp, char *section_title);





/* Functions that write or print EqParms */

/*** WriteParameters: writes the out_parm struct into a new section in the *
 *                    file specified by filename; the new 'eqparms' sec-   *
 *                    tion is inserted right after the 'input' section;    *
 *                    to achieve this, we need to write to a temporary     *
 *                    file which is then renamed to the output file name   *
 *              NOTE: lambdas are converted back into protein half lives!! *
 ***************************************************************************/

void WriteParameters(char *filename, EqParms *p, char *title, int ndigits);

/*** PrintParameters: prints an eqparms section with 'title' to the stream *
 *                    indicated by fp                                      *
 ***************************************************************************/

void PrintParameters(FILE *fp, EqParms *p, char *title, int ndigits);



/***** WriteDerivLog: write to solver log file *****************************/

void WriteDerivLog(char *deriv, int rule, int num_nuc);

/**** tvsub: Calculate the user cpu time in seconds ************/

double tvsub (struct rusage a, struct rusage b);

#endif
