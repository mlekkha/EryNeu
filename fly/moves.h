/*****************************************************************
 *                                                               *
 *   moves.h                                                     *
 *                                                               *
 *****************************************************************
 *                                                               *
 *   written by JR, modified by Yoginho                          *
 *                                                               *
 *****************************************************************
 *                                                               *
 * moves.h contains problem-specific stuff for the Lam annealer. *
 * There are functions for move generation and translation of    *
 * model data structures into annealing parameter arrays. More-  *
 * over there are a couple of functions that read annealing and  *
 * tune parameters and the tweak table from the data file.       * 
 *                                                               *
 *****************************************************************/

#ifndef MOVES_INCLUDED
#define MOVES_INCLUDED

/* this def needed for func. defs that refer to (* FILE) */
#ifndef _STDIO_INCLUDED
#include <stdio.h>    
#endif

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include <global.h>
#endif

/* following for StopStyle enum style */
#ifndef SA_INCLUDED
#include <sa.h>
#endif

/* ... and this for the Range struct */
#ifndef SCORE_INCLUDED
#include <score.h>
#endif



/*** CONSTANTS *************************************************************/

#define THETA_MIN  0.            /* minimum value of theta_bar (move size) */
#define THETA_INIT 5.0      /* initial value for all theta_bar (move size) */

#define LOWBITS    0x330E                /* following two are for drand to */
#define BYTESIZE   8                                   /* erand conversion */

 



/*** TONS OF STRUCTS *******************************************************/

/* The following are annealing parameters that are not specific to the Lam *
 * algorithm. In general they should be used in moves.c or fly_sa.c but    *
 * *not* in lsa.c. In the data file, the members of the struct labeled RO  *
 * are read from the $annealing_input section. They are used for initial   *
 * conditions of the annealer and do not change during a run. Members      *
 * labeled OUT are written to the $annealing_output section upon comple-   *
 * tion of a run.                                                          */

typedef struct {
  long   seed;                      /* seed for random number generator RO */
  double start_tempr;          /* the initial equilibration temperature RO */
  double gain;            /* gain for proportional control of move size RO */
  double stop_energy;                /* the final energy of the answer OUT */
  int    max_count;                      /* total number of iterations OUT */
  int    interval;       /* number of sweeps between updating theta_bar RO */
/*int    distribution;    1 - uniform; 2 - exp; 3 - normal; 4 - lorentz RO */
} AParms;

/* Struct for acceptance statistics, which are used for keeping the accep- *
 * tance ratio for each parameter as close as possible to .44; there's one *
 * such struct for each parameter to be tweaked                            */

typedef struct {
  double acc_ratio;                      /* acceptance ratio for parameter */
  double theta_bar;              /* theta bar is proportional to move size */
  int    hits;       /* number of moves since last call to UpdateControl() */
  int    success;              /* number of these moves that were accepted */
} AccStats;


/* Some annealers need to know about the search space to make the right    *
 * moves. For penalty type ranges, the penalty needs to be converted into  *
 * explicit range limits using Penalty2Limits (in score.c). Probably, this *
 * should be called from UpdateControl and give the ranges corresponding   *
 * to a defined penalty level                                              *
 * BUT: PASSING PARAM LIMITS TO THE ANNEALER IS CURRENTLY NOT SUPPORTED    */

typedef struct {
  double    *param;                /* pointers to parameters to be tweaked */
  Range     *param_range;        /* pointers to corresponding range limits */
} ParamList;


typedef struct {
  int       size;                           /* size of the ParamList array */
  ParamList *array;            /* points to 1st element of ParamList array */
  double    *pen_vec;   /* penalty vector: see score.h, struct SearchSpace */
} PArrPtr;



/* following struct contains copies of the static variables of moves.c     *
 * together with the values of parameters undergoing annealing             */

typedef struct { 
  ParamList *pt;  /* Used during a save to point to annealed-on parameters */
  AccStats  *acc_tab_ptr;            /* points to current acceptance stats */
  double    *newval; /* points to array of annealed-on doubles for restore */
  double    old_energy;                     /* energy before the last move */
  int       nparams;                      /* # of parameters to be tweaked */
  int       index;      /* index of parameter to be tweaked during a sweep */
  int       nhits;                         /* number of moves already made */
  int       nsweeps;                         /* number of completed sweeps */
} MoveState;



/* Tweak struct is for tweaking individual parameters or not; each pointer *
 * below points to an array of ints which represent each paramter          *
 * 1 means tweak, 0 means leave it alone                                   *
 * see ../doc/dataformatX.X for further details on the tweak section       */
 
typedef struct {
  int       *Rtweak;                             /* which Rs to be tweaked */
  int       *Ttweak;                             /* which Ts to be tweaked */
  int       *Etweak;                             /* which Es to be tweaked */
  int       *mtweak;                             /* which ms to be tweaked */
  int       *htweak;                             /* which hs to be tweaked */
  int       *dtweak;                             /* which ds to be tweaked */
  int       *lambdatweak;                   /* which lambdas to be tweaked */
  int       *tautweak;                   /* which taus to be tweaked */
} Tweak;


/* Opts struct is for saving command line options in savestate.c */

typedef struct {
  char      *inname;                             /* filename of input file */
  char      *outname;                           /* filename of output file */
  char      *argv;                                /* original command line */
  char      *derivfunc;                             /* derivative function */
  char      *solver;                                    /* solver function */
  StopStyle stop_flag;                                   /* stop criterion */
  int       prolix_flag;                                    /* prolix flag */
  int       landscape_flag;                              /* landscape flag */
  int       time_flag;                             /* flag for timing code */
  int       log_flag;                                  /* log display flag */
  long      state_write;              /* frequency for writing state files */
  long      print_freq;         /* frequency for printing status to stdout */
  long      captions;                         /* opt for printing captions */
  int       olddivstyle;                            /* division style flag */
  int       precision;                                 /* output precision */
  double    stepsize;                                  /* solver step size */
  int       quenchit;          /* flag for quenchit mode (T=0 immediately) */
  int       equil;              /* flag for equilibration mode (T = const) */
#ifdef MPI
  int       tuning;                                /* flag for tuning mode */
  int       covar_index; /* index for sample interval (=covar_index * tau) */
  int       write_tune_stat;    /* how many times to write tune statistics */
  int       auto_stop_tune;                      /* auto-stop tuning runs? */
#endif
} Opts;




/*** FUNCTION PROTOTYPES ***************************************************/

/* fly_sa.c: I/O functions for miscellaneous stuff */

/*** InitEquilibrate: reads the equilibrate section of the data file, ******
 *                    which is needed for equilibration runs               *
 ***************************************************************************/

void InitEquilibrate(FILE *fp);

/*** ReadTune: reads the tune_parameters section in a data file and ********
 *             turns a SAType structure to the caller                      *
 ***************************************************************************/

SAType ReadTune(FILE *fp);

/*** ReadAParameters: reads the AParm struct from an annealing_input sec- **
 *                    tion; these are the annealing parameters that are    *
 *                    not Lam-specific (and should NOT go into lsa.c       *
 ***************************************************************************/

AParms ReadAParameters(FILE *fp);

/*** WriteAParameters: writes the aparm struct into a new section in the ***
 *                     file specified by filename; the new 'annealing_out- *
 *                     put section is inserted right after the 'tune_para- *
 *                     meters' section; to achieve this, we need to write  *
 *                     to a temporary file which is then renamed to the    *
 *                     output file name                                    *
 ***************************************************************************/

void WriteAParameters(char *filename, AParms aparm);

/*** PrintAParameters: writes an 'annnealing_output' section with 'title' **
 *                     to the stream specified by fp                       *
 ***************************************************************************/

void PrintAParameters(FILE *fp, AParms aparm, char *title);

/*** WriteEquil: writes the equilibrate_variance section to the data file **
 *               right after the $equilibrate section                      *
 ***************************************************************************/

void WriteEquil(char *filename, double *equil_var);

/*** PrintEquil: writes an 'equilibrate_variance' section with 'title' *****
 *               to the stream specified by fp                             *
 ***************************************************************************/

void PrintEquil(FILE *fp, double *equil_var, char *title);

/*** PrintTimes: writes two (parallel: three) times sections ***************
 ***************************************************************************/

void PrintTimes(FILE *fp, double *delta);

/*** WriteVersion: prints the version and the complete command line used ***
 *                 to run fly_sa into the $version section of the data     *
 *                 file                                                    *
 ***************************************************************************/

void WriteVersion(char *filename);



/* functions that communicate with savestate.c */

/*** GetOptions: returns command line options to savestate.c ***************
 ***************************************************************************/

Opts *GetOptions(void);

/*** RestoreOptions: restores the values of the command line opt variables *
 ***************************************************************************/

void RestoreOptions(Opts *options);





/* translate.c: establishes communication between model and annealer */

/*** InitTweak: installs tweak as a static variable in translate.c; tweak **
 *              is read from the $tweak section in newstyle data files     *
 ***************************************************************************/

void InitTweak(FILE *fp);

/*** Translate: creates an array of pointers that point to the parameters **
 *              to be tweaked; the PArrPtr that is returned also includes  *
 *              pointers to the corresponding parameter ranges for each    *
 *              rameters, although this feature is not yet used anywhere   *
 *              in the annealing code                                      *
 *     CAUTION: InitZygote and InitScoring have to be called first!        *
 ***************************************************************************/

PArrPtr Translate(void);

/*** ReadTweak: reads the tweak array from the $tweak section in the data **
 *              file; this array has a value of 1 or 0 for each parameter  *
 *              in the model and is used by Translate to create the array  *
 *              of pointers to the parameters-to-be-tweaked                *
 ***************************************************************************/

Tweak ReadTweak(FILE *fp);





/* moves.c: functions for move generation */

/* Initializing and restoring functions */

/*** InitMoves: initializes the following moves.c-specific stuff: **********
 *              - static annealing parameter struct (ap)                   *
 *              - static tweak struct (tweak) in translate.c               *
 *              - initializes random number generator in lsa.c             *
 *              - receives parameter list from Translate, stores nparams   *
 *              - initializes acc_tab for acceptance statistics            *
 *              - set mixing interval in lsa.c (parallel code only)        *
 *                                                                         *
 *              it then returns the initial temperature to the caller      *
 ***************************************************************************/

double InitMoves(FILE *fp);

/*** RestoreMoves: restores move generator from state file *****************
 *           NOTE: InitMoves will be called before this function during    *
 *                 a restore                                               *
 ***************************************************************************/

void RestoreMoves(MoveState *MovePtr);

/*** RestoreProlix: restores prolix file after a run has been interrupted **
 ***************************************************************************/

void RestoreProlix(void);



/* a function for finalizing a run */

/*** GetFinalInfo: collects stop energy and final count for output to the **
 *                 data file                                               *
 ***************************************************************************/

AParms GetFinalInfo(void);



/* move generation functions that are used in moves.c, but not lsa.c */

/*** Move: tweaks the parameter-to-be-tweaked according to the current *****
 *         value of index; also calls UpdateControl if necessary           *
 ***************************************************************************/

void Move(void);

/*** UpdateControl: each interval number of steps, acceptance stats are ****
 *                  dated here; this function also prints prolix stuff, if *
 *                  required                                               *
 ***************************************************************************/

void UpdateControl(void);



/* functions that communicate with other source files */

/*** SetProlix: sets flag for printing prolix output on acceptance stats ***
 *              and initializes the prolix file if required                *
 ***************************************************************************/

void SetProlix(int value, char *file, int init_flag);

/*** MoveSave: returns a MoveState struct in which the current state of ****
 *             moves is saved; use for writing state file                  *
 ***************************************************************************/

MoveState *MoveSave(void);





/* savestate.c */

/* funcs that read/write/erase the state file */

/*** StateRead: reads Lam statistics, move state and erand state from a ****
 *              state file and restores the annealer's state to the same   *
 *              state it was in before it got interrupted                  *
 *     CAUTION: InitMoves should be called before StateRead!               *
 ***************************************************************************/

void StateRead(char *statefile, Opts *options, MoveState *move_ptr,
	       double *stats, unsigned short *rand, double *times);

/*** StateRm: removes the state file after the run has been completed ******
 ***************************************************************************/

void StateRm(void);

#endif
