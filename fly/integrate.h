/*****************************************************************
 *                                                               *
 *   integrate.h                                                 *
 *                                                               *
 *****************************************************************
 *                                                               *
 *   written by JR, modified by Yoginho                          *
 *                                                               *
 *****************************************************************
 *                                                               *
 * This file has problem-specific stuff not needed for right     *
 * hand of ODE. The function Blastoderm runs the fly model and   *
 * calls the appropriate solver for propagating the equations.   *
 * PrintBlastoderm formats and prints the output of the Blasto-  *
 * derm function and FreeSolution frees the solution allocated   *
 * by Blastoderm. integrate.h also comes with a few utility      *
 * functions for TLists which are linked lists used to initia-   *
 * lize the time/mode table for Blastoderm.                      *
 *                                                               *
 *****************************************************************
 *                                                               *
 * NOTE: all right-hand-of-ODE-specific stuff is in zygotic.h    *
 *                                                               * 
 *****************************************************************/


#ifndef INTEGRATE_INCLUDED
#define INTEGRATE_INCLUDED

/* this def needed for func. defs that refer to (* FILE) *******************/
#ifndef _STDIO_INCLUDED
#include <stdio.h>    
#endif


/* following for structures & consts used thruout **************************/
#ifndef GLOBAL_INCLUDED
#include <global.h>
#endif


/* following for NArrPtr */
#ifndef MATERNAL_INCLUDED
#include <maternal.h>
#endif


/* CONSTANTS FOR DBL_EPSILON ***********************************************/
/* EPSILON is used as minimal time step and corresponds to the minimal     */
/* double (or float) that can still be handled by C                        */

/* maybe the following declarations 
 * can be restricted, but for now they go here */

#ifndef _FLOAT_INCLUDED
#include <float.h>
#endif


#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/* DBL_EPSILON is about 2 x 10^-16. so BIG_EPSILON is ~10^-11 min. */ 

#define      EPSILON     (DBL_EPSILON * 1000.)  
                                  /* This assumes all times < 100. min. !! */

#ifdef       FLOAT_EVERYTHING
#undef       EPSILON
#define      EPSILON     (FLT_EPSILON * 100.)             /* see above */
#endif

#define      BIG_EPSILON (EPSILON * 1000.)
#define      HALF_EPSILON (EPSILON * .5)

/* MORE CONSTANTS: OPS FOR BLASTODERM -- WITH PRIORITIES *******************/

#define      ADD_BIAS     1    /* can happen any step */

#define      NO_OP        2    /* ops below: first found set executed      */
#define      DIVIDE       4    /* NO_OP does nothing, DIVIDE divides nucs  */
#define      PROPAGATE    8    /* and PROPAGATE propagates the equations   */
#define      MITOTATE     16    /* and MITOTATE carries forward the
								system by epsilon (to handle rounding 
								errors from the solver   */





/* A STRUCT ****************************************************************/

typedef struct TList {          /* Tlist is a linked list we use to set up */
  double       time;            /* the solution struct in integrate.c; the */
  int          n;               /* TList has an element for each time con- */
  int          op;              /* taining the number of elements for the  */
  struct TList *next;           /* solution struct at that time (n) and a  */
} TList;                        /* rule (op) to tell the solver what to do */



typedef struct InterpObject {

	
	double 				*fact_discons;
	int 				fact_discons_size;
	NArrPtr				func;
	NArrPtr				slope;
	int 				maxsize; 
	double 				maxtime;


} InterpObject;


/*** NEW GLOBAL ************************************************************/

void         (*ps)(double *, double *, double, double, 
		   double, double, int, FILE *);
                                                     /* this is the solver */

double         (*eval)(NArrPtr, int); /* evaluation function */

double SquareEval (NArrPtr Solution, int gindex); /* square evaluation function 
                                                   */ 
double AbsoluteEval (NArrPtr Solution, int gindex); /* absolute evaluation
                                                     * function */

/* FUNCTION PROTOTYPES *****************************************************/

/* Blastoderm Functions */

/*** Blastoderm: runs embryo model and returns an array of concentration ***
 *               arrays for each requested time (given by TabTimes) using  *
 *               stephint as a suggested stepsize and accuracy as the re-  *
 *               quired numerical accuracy (in case of adaptive stepsize   *
 *               solvers or global stepsize control) for the solver.       *
 *         NOTE: TabTimes *must* start from 0 and have increasing times.   *
 *               It includes times when bias is added, cell division times *
 *               and times for which we have data and ends with the time   *
 *               of gastrulation.                                          *
 ***************************************************************************/


NArrPtr Blastoderm(int genindex, char *genotype, 
						InterpObject hist_interrp,
						InterpObject extinp_interrp, DArrPtr tabtimes,
		   				double stephint, double accuracy, FILE *slog);
/*** PrintBlastoderm: writes the output of the model to a stream specified *
 *                    by the fp file pointer. The Table is a solution of   *
 *                    the model as returned by Blastoderm, the id speci-   *
 *                    fies the title of the output and ndigits specifies   *
 *                    the floating point precision to be printed.          *
 *                    PrintBlastoderm adjusts its format automatically to  *
 *                    the appropriate number of genes.                     *
 ***************************************************************************/

void PrintBlastoderm(FILE *fp, NArrPtr table, char *id, 
		     int ndigits, int columns);

/*** ConvertAnswer: little function that gets rid of bias times, division **
 *                  times and such and only returns the times in the tab-  *
 *                  times struct as its output. This is used to produce    *
 *                  unfold which is nothing but the requested times.       *
 ***************************************************************************/

NArrPtr ConvertAnswer(NArrPtr answer, DArrPtr tabtimes);

/*** FreeSolution: frees memory of the solution structure created by *******
 *                 Blastoderm() or gut functions                           *
 ***************************************************************************/

void FreeSolution(NArrPtr *solution);





/* TList Utility Functions */

/*** InitTList: initializes TList and adds first (t=0) and last ************
 *              (t=gastrulation) element of Tlist.                         *
 ***************************************************************************/

TList *InitTList(void);

/*** InsertTList: takes pointer to first element of TList plus time and ****
 *                desired op for new TList element and inserts a new TList *
 *                element for time at the appropriate place within the     *
 *                linked list. This function returns a pointer to the      *
 *                first element of the TList if insertion worked out fine. *
 ***************************************************************************/

TList *InsertTList(TList *first, double time, int op);

/*** CountEntries: counts how many entries we have in a TList **************
 ***************************************************************************/

int CountEntries(TList *first);

/*** FreeTList: frees the memory of a TList ********************************
 ***************************************************************************/

void FreeTList(TList *first);

NArrPtr Dat2NArrPtr(DataTable *table, int *maxind);
void Go_Forward(double *output, double *input, int output_ind, int
input_ind, int num_genes);
void Go_Backward(double *output, double *input, int output_ind, int
input_ind, int num_genes);
void History(double t, double t_size, double *yd, int n);
double *GetFactDiscons(int *sss);
void FreeInterpObject(InterpObject *interp_obj);
void DoInterp(DataTable *interp_dat, InterpObject *interp_res, 
												int num_genes);
void SetHistoryInterp(InterpObject interp_info);
void SetExternalInputInterp(InterpObject interp_info);
void ExternalInputs(double t, double t_size, double *yd, int n);
void SetFactDiscons(void);
void FreeFactDiscons(void);
void TestInterp(int num_genes, int type);
void FreeExternalInputTemp(void);
void  FreeHistoryTemp(void);


#endif
