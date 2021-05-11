/*******************************************************************
 *                                                                 *
 *   fly_sa.c                                                      *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 *   written by JR, modified by Yoginho                            *
 *   -D landscape option by Lorraine Greenwald Oct 2002            *
 *   -g option by Yousong Wang, Feb 2002                           *
 *   -a option by Marcel Wolf, Apr 2002                            *
 *                                                                 *
 ******************************************************************* 
 *                                                                 *
 * Although main() is in lsa.c, this is the file that 'defines'    *
 * the fly_sa program, since it contains most of its problem-      *
 * specific code (except for move generation -> moves.c, saving    *
 * of intermediate state files -> savestate.c and communication    *
 * with the specific cost function that is used -> translate.c).   *
 *                                                                 *
 * After I've told you all that's NOT in this file, here's what    *
 * the funcs below actually do: parsing fly_sa command line opts   *
 * is one of its jobs; there are funcs that make the first and     *
 * last moves and funcs that read and write Lam and Lam-indepen-   *
 * dent annealing parameters to the problem-specific data file.    *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 * Copyright (C) 1989-2003 John Reinitz                            *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License as  *
 * published by the Free Software Foundation; either version 2 of  *
 * the License, or (at your option) any later version.             *
 *                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 * GNU General Public License for more details.                    *
 *                                                                 *
 * You should have received a copy of the GNU General Public       *
 * License along with this program; if not, write to the           *
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,     *
 * Boston, MA  02111-1307, U.S.A.                                  *
 *                                                                 *
 *******************************************************************/

#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>                       /* for command line option stuff */

#include <error.h>                                 /* error handling funcs */
#include <distributions.h>               /* DistP.variables and prototypes */
#include <integrate.h>
#include <maternal.h>                          /* for olddivstyle and such */
#include <moves.h>                     /* problem-specific annealing funcs */
#include <random.h>                                     /* for InitERand() */
#include <sa.h>                     /* problem-independent annealing funcs */
#include <score.h>                             /* for init and Score funcs */
#include <solvers.h>                           /* for name of solver funcs */
#include <zygotic.h>            /* for init, mutators and derivative funcs */

#ifdef MPI                 /* this inludes parallel-specific stuff for MPI */
#include <mpi.h>                     /* this is the official MPI interface */
#include <MPI.h>  /* our own structs and such only needed by parallel code */
#endif

/*** Constants *************************************************************/

/* command line option string */

/* #define  OPTS       ":a:b:Bc:C:d:De:Ef:g:hi:lLnopQr:s:StTvw:W:y:" */
                                            
#define  OPTS       ":a:b:Bc:C:De:Ef:g:hi:k:lLnNopQr:s:StTvw:W:y:"
                                             /* command line option string */
                                             /* D will be debug, like scramble, score */
                     /* must start with :, option with argument must have a : following */


/*** STATIC VARIABLES ******************************************************/

/* Help, usage and version messages */

#ifdef MPI
static const char usage[]    =
"Usage: fly_sa.mpi [-b <bkup_freq>] [-B] [-C <covar_ind>] \n"
"                  [-D] [-e <freeze_crit>][-E] [-f <param_prec>] [-g <g(u)>]\n"
"                  [-h] [-i <stepsize>] [-k <cost_function>] [-l] [-L] [-n] [-N] [-p] [-s <solver>]\n"
"                  [-S] [-t] [-T] [-v] [-w <out_file>]\n"
"                  [-W <tune_stat>] [-y <log_freq>]\n"
"                  <datafile>\n";
#else
static const char usage[]    =
"Usage: fly_sa [-a <accuracy>] [-b <bkup_freq>] [-B] [-e <freeze_crit>] [-E]\n"
"              [-f <param_prec>] [-g <g(u)>] [-h] [-i <stepsize>]\n"
"              [-l] [-L] [-n] [-N] [-p] [-Q] [-s <solver>] [-t] [-v]\n"
"              [-w <out_file>] [-y <log_freq>]\n"
"              <datafile>\n";
#endif

static const char help[]     =

#ifdef MPI
"Usage: fly_sa.mpi [options] <datafile>\n\n"
#else
"Usage: fly_sa [options] <datafile>\n\n"
#endif

"Argument:\n"
"  <datafile>          input data file\n\n"

"Options:\n"
"  -a <accuracy>       solver accuracy for adaptive stepsize ODE solvers\n"
"  -b <bkup_freq>      write state file every <bkup_freq> * tau moves\n"
"  -B                  run in benchmark mode (only do fixed initial steps)\n"
#ifdef MPI
"  -C <covar_ind>      set covar sample interval to <covar_ind> * tau\n"
#endif
"  -D                  debugging mode, prints all kinds of debugging info\n"
"  -e <freeze_crit>    set annealing freeze criterion to <freeze_crit>\n"
"  -E                  run in equilibration mode\n"
"  -f <param_prec>     float precision of parameters is <param_prec>\n"
"  -g <g(u)>           chooses g(u): e = exp, h = hvs, s = sqrt, t = tanh\n"
"  -h                  prints this help message\n"
"  -i <stepsize>       sets ODE solver stepsize (in minutes)\n"
"  -k <cost_function> select error function for calculation the differences between\n"
"                      equation and data: sq=squared differences, abs=absolute differences\n"
"                      default=sq\n"
"  -l                  echo log to the terminal\n"

#ifdef MPI
"  -L                  write local logs (llog files)\n"
#endif
"  -n                  nofile: don't print .log or .state files\n"
"  -N                  generates landscape to .landscape file in equilibrate mode \n"
"  -o                  use oldstyle cell division times (3 div only)\n"
"  -p                  prints move acceptance stats to .prolix file\n"
#ifndef MPI
"  -Q                  quenchit mode, T is lowered immediately to zero\n"
#endif
"  -s <solver>         choose ODE solver\n"
#ifdef MPI
"  -S                  disable tuning stop flag\n"
#endif
"  -t                  write timing information to .times file\n"
#ifdef MPI
"  -T                  run in tuning mode\n"
#endif
"  -v                  print version and compilation date\n"
"  -w <out_file>       write output to <out_file> instead of <datafile>\n"
#ifdef MPI
"  -W <tune_stat>      tuning stats written <tune_stat> times per interval\n"
#endif
"  -y <log_freq>       write log every <log_freq> * tau moves\n\n"

"Please report bugs to <yoginho@usa.net>. Thank you!\n";

static char version[MAX_RECORD];                 /* version gets set below */

static const char verstring[] = 

"compiled by:      %s\n"
"         on:      %s\n"
"      using:      %s\n"
"      flags:      %s\n"
"       date:      %s at %s\n";

/* other static variables */

static char   *inname;                           /* filename of input file */
static char   *outname;                         /* filename of output file */

static char   *argvsave;          /* static string for saving command line */
 
static double stepsize    = 1.;                     /* stepsize for solver */
static double accuracy    = 0.001;   /* accuracy for solver (not used yet) */
static int    precision   = 8;                    /* precision for eqparms */
static int    prolix_flag = 0;               /* to prolix or not to prolix */
static int    landscape_flag = 0;        /* generate energy landscape data */
static char   *cost_function = "sq"; /* error function utilized in error
                                       * evaluation in Eval method      */
           /* set the landscape flag (and the landscape filename) in lsa.c */

/* the following lines define a pointers to:                               */
/*            - pd:    dvdt function, currently only DvdtOrig in zygotic.c */
/*            - pj:    Jacobian function, in zygotic.c                     */
/* This pointer needs to be static since both InitialMove and Restore-     */
/* State need it; it gets set in ParseCommandLine()                        */
/*                                                                         */
/* NOTE: ps (solver) is declared as global in integrate.h                  */

static void (*pd)(double *, double, double *, int);
static void (*pj)(double, double *, double *, double **, int);





/*** FUNCTIONS *************************************************************/

/*** COMMAND LINE OPTS ARE ALL PARSED HERE *********************************/

/*** ParseCommandLine: well, parses the command line and returns an index **
 *                     to the 1st argument after the command line options  *
 ***************************************************************************/
 
int ParseCommandLine(int argc, char **argv)
{
  int         c;                     /* used to parse command line options */
  
/* external declarations for command line option parsing (unistd.h) */

  extern char *optarg;                     /* command line option argument */
  extern int  optind;                /* pointer to current element of argv */
  extern int  optopt;               /* contain option character upon error */

/* set the version string */

#ifdef MPI    
      sprintf(version, "fly_sa version %s parallel", VERS);
#else
#ifdef ALPHA_DU
      sprintf(version, "fly_sa version %s serial dec-kcc-dxml", VERS);
#else
      sprintf(version, "fly_sa version %s serial", VERS);
#endif
#endif

/* following part sets default values for command line options */

  pd              = DvdtOrig;               /* default derivative function */
  pj              = JacobnOrig;               /* default Jacobian function */
  ps              = Rk4;                                 /* default solver */
  eval            = SquareEval;             /* default evaluation function */

  captions        = 100000000;  /* default freq for writing captions (off) */
  print_freq      = 100;           /* default freq for writing to log file */
  state_write     = 1000;        /* default freq for writing to state file */

  stop_flag       = absolute_freeze;             /* type of stop criterion */
  time_flag       = 0;                                  /* flag for timing */
  log_flag        = 0;              /* flag for writing logs to the screen */
  nofile_flag     = 0;       /* flog for not writing .log and .state files */

#ifdef MPI
  tuning          = 0;                    /* tuning mode is off by default */
  covar_index     = 1;      /* covariance sample will be covar_index * tau */
  write_tune_stat = 1;         /* how many times do we write tuning stats? */
  auto_stop_tune  = 1;               /* auto stop tuning runs? default: on */
  write_llog      = 0; /* write local llog files when tuning; default: off */
#endif

/* following part parses command line for options and their arguments      */
  
  optarg = NULL;
  while( (c = getopt(argc, argv, OPTS)) != -1 ) {
    switch(c) {
    case 'a':
      accuracy = atof(optarg);
      if ( accuracy <= 0 )
	error("fly_sa: accuracy (%g) is too small", accuracy);
      break;
    case 'b':            /* -b sets backup frequency (to write state file) */
      state_write = strtol(optarg, NULL, 0);
      if ( state_write < 1 )
	error("fly_sa: max. backup frequency is every tau steps i.e. -b 1");
      if ( state_write == LONG_MAX )
	error("fly_sa: argument for -b too large");
      break;
    case 'B':      /* -B runs in benchmark mode -> quit after intial steps */
      bench = 1;
      time_flag = 1;
      break;
    case 'c':               /* -c sets the frequency for printing captions */
      error("fly_sa: -c is not supported anymore, captions off for good");
/* if you want to be able to insert captions into .log files, uncomment    *
 * the following lines; CAUTION: make sure that RestoreLog() works proper- *
 * ly with captions before you do this!                                    */
/*      captions = strtol(optarg, NULL, 0);
      if ( captions < 1 )
        error("fly_sa: can't print captions more than every line (-c 1)");
      if ( captions == LONG_MAX )
        error("fly_sa: argument for -c too large");                        */
      break;
    case 'C':             /* parallel code: -C sets the covar sample index */
#ifdef MPI
/* -C is currently disabled since I (Yogi) could not get it to work and I  *
 * got tired of debugging and I didn't really need it at that time and so  *
 * on and so forth; try running on a smaller number of nodes or get it to  *
 * work yourself; that's what I say. Grmbl.                                */
/*    covar_index = atoi(optarg); 
      if ( covar_index < 1 )
        error("fly_sa: covariation sample index must be >= 1");            */
      error("fly_sa: -C does not work yet, try to run on fewer nodes");
#else
      error("fly_sa: can't use -C in serial, tuning only in parallel");
#endif
      break;
    case 'D':      /* -D is debug mode */
      debug = 1;
      break;
    case 'e':                            /* -e sets the stopping criterion */
      if( !(strcmp(optarg,"pfreeze")) ) 
	stop_flag = proportional_freeze;  
      else if( !(strcmp(optarg,"afreeze")) )
	stop_flag = absolute_freeze;
      else if( !(strcmp(optarg,"abs")) )
	stop_flag = absolute_energy;
      else
	error("fly_sa: valid stopping criteria are pfreeze, afreeze, abs");
      break;
    case 'E':                                /* -E does equilibration runs */
      equil = 1;
      break;
    case 'f':
      precision = atoi(optarg);           /* -f determines float precision */
      if ( precision < 0 )
	error("fly_sa: what exactly would a negative precision be???");
      if ( precision > MAX_PRECISION )
	error("fly_sa: max. float precision is %d!", MAX_PRECISION);
      break;
   case 'g':                                    /* -g choose g(u) function */
      if( !(strcmp(optarg, "s")) )
	gofu = Sqrt;                          
      else if ( !(strcmp(optarg, "t"))) 
	gofu = Tanh;
      else if ( !(strcmp(optarg, "e"))) 
	gofu = Exp;
      else if ( !(strcmp(optarg, "h"))) 
	gofu = Hvs;
      else 
	error("fly_sa: %s is an invalid g(u), should be e, h, s or t",
	      optarg); 
     break;
    case 'h':                                            /* -h help option */
      PrintMsg(help, 0);
      break;
    case 'i':                                      /* -i sets the stepsize */
      stepsize = atof(optarg);
      if ( stepsize < 0 )
	error("fly_sa: going backwards? (hint: check your -i)");
      if ( stepsize == 0 ) 
	error("fly_sa: going nowhere? (hint: check your -i)");
      if ( stepsize > MAX_STEPSIZE )
	error("fly_sa: stepsize %g too large (max. is %g)", 
	      stepsize, MAX_STEPSIZE);
      break;
    
    case 'k':       
      if ( !(strcmp(optarg, "sq")) )
          eval = SquareEval;
      if ( !(strcmp(optarg, "abs")) )
          eval = AbsoluteEval;
      break;
      
    case 'l':                         /* -l displays the log to the screen */
      log_flag = 1;
      break;
    case 'L':                   /* -L writes local .llog files when tuning */
#ifdef MPI
      write_llog = 1;
#else
      error("fly_sa: can't use -L in serial, tuning only in parallel");
#endif
      break;
    case 'n':                       /* -n suppresses .state and .log files */
      nofile_flag = 1; 
      break;
    case 'N':                  /* -N sets laNdscape flag and Equilibrate mode */
       equil = 1;
       landscape_flag = 1;
      break;        
    case 'o':             /* -o sets old division style (ndivs = 3 only! ) */
      olddivstyle = 1;
      break;
    case 'p':                                   /* -p sets prolix printing */
      prolix_flag = 1;
      break;
    case 'Q':                  /* -Q sets quenchit mode (serial code only) */
#ifdef MPI
      error("fly_sa: can't use -Q in parallel, quenchit only in serial");
#else
      quenchit = 1;
#endif
      break;
    case 'r':
      error("fly_sa: -r is currently not supported, use -g instead");
      break;
    case 's':                                 /* -s sets solver to be used */
      if ( !(strcmp(optarg, "a")) )
	ps = Adams;
      else if ( !(strcmp(optarg, "bd")) )
	ps = BaDe;
      else if ( !(strcmp(optarg, "bs")) )
	ps = BuSt;
      else if ( !(strcmp(optarg, "e")) )
	ps = Euler;
      else if ( !(strcmp(optarg, "h")) )
	ps = Heun;
      else if ( !(strcmp(optarg, "mi")) || !(strcmp(optarg, "m")) )
	ps = Milne;
      else if ( !(strcmp(optarg, "me")) )
	ps = Meuler;
      else if ( !(strcmp(optarg, "r4")) || !(strcmp(optarg, "r")) )
	ps = Rk4;
      else if ( !(strcmp(optarg, "r2")) )
	ps = Rk2;
      else if ( !(strcmp(optarg, "rck")) )
	ps = Rkck;
      else if ( !(strcmp(optarg, "rf")) )
	ps = Rkf;
      else if ( !(strcmp(optarg, "sd")) )
	ps = SoDe;
      else 
	error("fly_sa: invalid solver (%s), use: a,bs,e,h,mi,me,r{2,4,ck,f}",
	      optarg);
      break;
    case 'S':                         /* -S unsets the auto_stop_tune flag */
#ifdef MPI
      auto_stop_tune = 0;
#else
      error("fly_sa: can't use -S in serial, tuning only in parallel");
#endif   
      break;
    case 't':                   /* -t saves times in data and .times files */
      time_flag = 1;
      break;
    case 'T':                       /* -T does tuning (parallel code only) */
#ifdef MPI
      tuning = 1;
#else
      error("fly_sa: can't use -T in serial, tuning only in parallel");
#endif
      break;
 
    case 'v':                                  /* -v prints version message */
      fprintf(stderr, "%s\n", version);
      fprintf(stderr, verstring, 
	      USR, MACHINE, COMPILER, FLAGS, __DATE__, __TIME__);
      exit(0);
    case 'w':                                       /* -w sets output file */
      outname = (char *)calloc(MAX_RECORD, sizeof(char)); 
      outname = strcpy(outname, optarg);
      SetOutname(outname);                /* communicates outname to lsa.c */
      break;
    case 'W':       /* -W determines the frequency of writing tuning stats */
#ifdef MPI
      write_tune_stat = atoi(optarg);
      if ( write_tune_stat < 1 )
	error("fly_sa: frequency of writing tune stats must be >= 1");
#else
      error("fly_sa: can't use -W in serial, tuning only in parallel");
#endif
      break;
    case 'y':                             /* -y set frequency to print log */
      print_freq = strtol(optarg, NULL, 0);
      if ( print_freq < 1 )
	error("fly_sa: can't print status more than every tau steps (-y 1)");
      if ( print_freq == LONG_MAX )
	error("fly_sa: argument for -y too large");
      break;
    case ':':
      error("fly_sa: need an argument for option -%c", optopt);
      break;
    case '?':
    default:
      error("fly_sa: unrecognized option -%c", optopt);
    }
  }

/* error checking here */

#ifdef MPI
  if ( (tuning == 1) && (equil == 1) )
    error("fly_sa: can't combine -E with -T");
  if ( write_llog && !tuning )
    error("fly_sa: -L only makes sense when tuning");
#else
  if ( (quenchit == 1) && (equil == 1) )
    error("fly_sa: can't combine -E with -Q");
#endif

  if ( ((argc - (optind - 1)) != 2) ) 
    PrintMsg(usage, 1);

  return optind;
}




/*** FUNCTIONS THAT INITIALIZE/RESTORE THINGS ******************************/

/**********************************************************************
 * This routine reads the distribution parameters  from the input file*
 * and stores them in DistP.xx from distributions.h                   *
 * LG 03-02: need q for gen visiting distribution input file          *
 * LG 05-02: set factors only dependent on q for general visiting     *
 * distribution by calling qgt2_init or qlt2_init from distributions.c*
 **********************************************************************/

void InitDistribution(FILE *fp)
{ /* begin init distribution */

   /* read in the data from file */
   /* read in header info then output it*/
   /* At the end of file exit from the loop */


  fp = FindSection(fp,"distribution_parameters"); /* find distribution section */

  if( !fp ) /* to be backward compatible set default to exp */
    {
    printf("ReadTune: could not locate distribution parameters:using exponential.\n");
    DistP.distribution = 1;
    DistP.q == 1.0;
    } /* end of default to exp */

  else 
    { /* input user selected distribution parameters */
      fscanf(fp,"%*s\n");    /* advance past title line no blanks allowed! */
      /* read distribution stuff */
      if ( 2!= (fscanf(fp,"%d %lf\n",  &(DistP.distribution), &(DistP.q) )) )
	error ("ReadTune: error reading distribution stuff"); 
 
      if ((DistP.distribution > 11) || (DistP.distribution < 1))
	{error ("fly_sa: distribution must be int from 1 to 11 \n");
	}
      else if  ((DistP.distribution == 4)||(DistP.distribution == 3))
	{error ("fly_sa: PLEASE use 5 for Lorentz or 10 for normal distribution \n");
	}
      else if ((DistP.distribution == 6)||(DistP.distribution == 9))
	{error ("fly_sa: 6=poisson or 9=pareto distribution returns positive values--do not use for fly \n");
	}
      else if (DistP.distribution == 7 )
	{  /* general distribution */
	  if ((DistP.q >= 3.0) || (DistP.q <= 1.0))
	    {error ("tsp_sa: q must be between 1 and 3 \n");}
	  else if (DistP.q == 2.0)
	    {DistP.distribution = 5;
	    /* fly needs lorentz, tsp use abs lorentz(4)*/
	    printf ("fly_sa: q=2 is lorentz--setting distribution to 5\n");}
	  else if (DistP.q > 2.0)
	    {qgt2_init(); }
	  else 
	    {qlt2_init(); }
	} /* end of general distribution */
	    /***************LG 05-02***************/ 
	    /* calculate q dependent factors that */
	    /* do not change for entire run.      */
	    /* these live in distributions.h      */
	    /***************LG 05-02***************/ 

    } /* end of input user selected distribution parameters */


}  /* end InitDistribution */

/*** InitialMove: initializes the following stuff: *************************
 *                - various static things for filenames and such           *
 *                - Lam parameter struct for use in lsa.c (from tune sect) *
 *                - model and scoring funcs & solver stepsize              *
 *                - move generation in moves.c                             *
 *                - sets initial energy and checks validity of initial     *
 *                  parameters according to limit ranges                   *
 *                then it returns the initial temperature to the caller    *
 ***************************************************************************/

double InitialMove(int argc, char **argv, int opt_index, 
		   NucStatePtr state_ptr, double *p_chisq)
{
  FILE    *infile;

  int     i;

  char    *p;

  double  i_temp;
  double  energy;

  SAType  in_tune;          /* temporary struct to read in tune parameters */



/* save some things in static storage for FinalMove() and StateWrite() */

  inname  = (char *)calloc(MAX_RECORD + 1, sizeof(char));    /* input file */
  inname  = strcpy(inname, argv[opt_index]);

  if ( !outname ) {             /* in case there's no -w, outname = inname */
    outname = (char *)calloc(MAX_RECORD + 1, sizeof(char));
    outname = strcpy(outname, inname);
  }                                      

  argvsave = (char *)calloc(MAX_RECORD, sizeof(char));
  for(i=0; i < argc; i++){
    if ( i>0 )
      argvsave = strcat(argvsave, " ");
    argvsave = strcat(argvsave, argv[i]);
  }
  argvsave = strcat(argvsave, "\n");

/* set the prolix flag (and the .prolix filename) in moves.c, if necessary */

  if ( prolix_flag )
    SetProlix(prolix_flag, outname, 1);       /* 1 means: new .prolix file */

  if ( landscape_flag )                      /* 1 means gen landscape files */
    InitLandscape(landscape_flag, outname);    /*lives in lsa.c */

/* initialize some Lam/Greening structures */

  p = state_ptr->tune.progname;     /* tune.progname contains program name */
  p = strcpy(p, version);

  state_ptr->tune.debuglevel = 0;          /* following stuff not used now */
  p = state_ptr->tune.tunename;
  p = strcpy(p, "The Other One");        /* Grateful Dead tune, what else? */
  state_ptr->tune.tunefile = NULL;

/* read data file and initialize different things for scoring and move     */
/* generation (in zygotic.c, score.c and moves.c                           */

  infile = fopen(argv[opt_index], "r");
  if ( !infile )
    file_error("fly_sa");

  InitZygote(infile, pd, pj, "input");     
     /* init problem, parms, bias, bcd, static structs for deriv and deriv */
  InitScoring(infile);                     /* initializes facts and limits */
  InitHistory(infile);

  InitExternalInputs(infile);

  InitStepsize(stepsize, accuracy, NULL, argv[opt_index], cost_function); 

  in_tune = ReadTune(infile);               /* read tune_parameter section */
  i_temp  = InitMoves(infile);   /* set initial temperature and initialize */
                       /* random number generator and annealing parameters */
  InitDistribution(infile);   /* initialize distribution stuff */
  fclose(infile);

/* initialize Lam parameters (see sa.h for further detail) */

  state_ptr->tune.lambda              = in_tune.lambda;
  state_ptr->tune.lambda_mem_length_u = in_tune.lambda_mem_length_u;
  state_ptr->tune.lambda_mem_length_v = in_tune.lambda_mem_length_v;
  state_ptr->tune.control             = in_tune.control;
  state_ptr->tune.initial_moves       = in_tune.initial_moves;
  state_ptr->tune.tau                 = in_tune.tau;
  state_ptr->tune.freeze_count        = in_tune.freeze_count;
  state_ptr->tune.update_S_skip       = in_tune.update_S_skip;
  state_ptr->tune.criterion           = in_tune.criterion;
#ifdef MPI
  state_ptr->tune.mix_interval        = in_tune.mix_interval;
#endif

  energy = Score();                      /* score first time to check if:  */
                                         /* - parameters within range?     */
  if ( energy == FORBIDDEN_MOVE )        /* - set initial energy (p_chisq) */
    error("fly_sa: the initial state was forbidden, cannot proceed");

  *p_chisq = energy;                                  /* set initial score */

  return(i_temp);                                   /* initial temperature */
}



/*** RestoreState: called when an interrupted run is restored; does the ****
 *                 following:                                              *
 *                 - stores various static things for filenames and such   *
 *                 - initializes Lam parameters for lsa.c                  *
 *                 - initializes model and scoring funcs & solver stepsize *
 *                 - initializes move generation in moves.c                *
 *                 - restores state at which previous run was interrupted  *
 *                   according to state file                               *
 *                                                                         *
 * Comment by JR:  RestoreState was originally going to be implemented with*
 * branches in InitialMove. That won't work because when when this func.   *
 * returns, control should go right to Loop(), skipping all the initiali-  *
 * zation stuff in Initialize. Hence most of the code in InitialMove is    *
 * just repeated here.                                                     *
 ***************************************************************************/

void RestoreState(char *statefile, NucStatePtr state_ptr, double *p_chisq)
{
  FILE           *infile;                /* pointer to original input file */
  char           *p;                                   /* temporary string */
  double         i_temp;
  SAType         in_tune;          /* temporary struct for tune parameters */

  Opts           *options;         /* used to restore command line options */
  MoveState      *move_ptr;                       /* used to restore moves */
  double         *stats;                      /* used to restore Lam stats */
  unsigned short *rand;                         /* used to restore ERand48 */
  double         delta[2];                        /* used to restore times */

/* allocate memory for structures that will be returned (stats gets allo-  *
 * cated in StateRead(), since we need to know if we're tuning or not)     */

  options = (Opts *)malloc(sizeof(Opts));
  options->inname    = (char *)calloc(MAX_RECORD, sizeof(char));
  options->outname   = (char *)calloc(MAX_RECORD, sizeof(char));
  options->argv      = (char *)calloc(MAX_RECORD, sizeof(char));
  options->derivfunc = (char *)calloc(MAX_RECORD, sizeof(char));
  options->solver    = (char *)calloc(MAX_RECORD, sizeof(char));

  stats    = (double *)calloc(31, sizeof(double));
  move_ptr = (MoveState *)malloc(sizeof(MoveState));
  rand     = (unsigned short *)calloc(3, sizeof(unsigned short));

  StateRead(statefile, options, move_ptr, stats, rand, delta);

/* restore options in fly_sa.c (and some in lsa.c) */

  RestoreOptions(options);

/* initialize some Lam/Greening structures */

  p = state_ptr->tune.progname;     /* tune.progname contains program name */
  p = strcpy(p, version);

  state_ptr->tune.debuglevel = 0;          /* following stuff not used now */
  p = state_ptr->tune.tunename;
  p = strcpy(p,"The Other One");         /* Grateful Dead tune, what else? */
  state_ptr->tune.tunefile = NULL;

/* read data file and initialize different things for scoring and move     */
/* generation (in zygotic.c, score.c and moves.c                           */

  infile = fopen(inname, "r");
  if ( !infile )
    file_error("fly_sa");
  
  InitZygote(infile, pd, pj, "input");      
                                  /* init initial cond., mutator and deriv */
  InitScoring(infile);                            /* init facts and limits */
  InitHistory(infile);

  InitExternalInputs(infile);
  
  InitStepsize(stepsize, accuracy, NULL, inname, cost_function);

  in_tune = ReadTune(infile);               /* read tune_parameter section */
  i_temp  = InitMoves(infile);   /* set initial temperature and initialize */
                                              /* various things in moves.c */
  InitDistribution(infile);   /* initialize distribution stuff */
  fclose(infile);

/* initialize Lam parameters (see sa.h for further detail) */

  state_ptr->tune.lambda              = in_tune.lambda;
  state_ptr->tune.lambda_mem_length_u = in_tune.lambda_mem_length_u;
  state_ptr->tune.lambda_mem_length_v = in_tune.lambda_mem_length_v;
  state_ptr->tune.control             = in_tune.control;
  state_ptr->tune.initial_moves       = in_tune.initial_moves;
  state_ptr->tune.tau                 = in_tune.tau;
  state_ptr->tune.freeze_count        = in_tune.freeze_count;
  state_ptr->tune.update_S_skip       = in_tune.update_S_skip;
  state_ptr->tune.criterion           = in_tune.criterion;
#ifdef MPI
  state_ptr->tune.mix_interval        = in_tune.mix_interval;
#endif

  RestoreMoves(move_ptr);
  RestoreLamstats(stats);
  if ( time_flag )
    RestoreTimes(delta);
  InitERand(rand);
  if ( prolix_flag )
    RestoreProlix();
}





/*** THE FINAL MOVE FUNCTION ***********************************************/

/*** FinalMove: reads final energy and move count, then prints $version, ***
 *              $annealing_output and $eqparms sections and removes the    *
 *              state file                                                 *
 ***************************************************************************/

void FinalMove(void)
{
#ifdef MPI
  int      i;                                              /* loop counter */
#endif

  AParms   ap;
  EqParms  *parm;

  double   equil_var[2];         /* array for results of equilibration run */

  char     *shell_cmd;           /* shell command to be executed by system */

#ifdef MPI
  int      winner = 0;                           /* id of the winning node */
  double   minyet = DBL_MAX;         /* minimum score, used to find winner */

  double   *final_e;               /* array of final energies of all nodes */

  final_e = (double *)calloc(nnodes, sizeof(double));
#endif 

/* get the answer and some additional info */

  ap   = GetFinalInfo();              /* reads final energy and move count */
  parm = GetParameters();
  
  if ( equil )
    GetEquil(equil_var);            /* get the final equilibration results */

#ifdef MPI

/* parallel code: find the node with the lowest score */

  for(i=0; i<nnodes; i++)                   /* initialize the energy array */
    final_e[i] = 0;  
                                /* collect the final scores from all nodes */
  MPI_Allgather(&ap.stop_energy, 1, MPI_DOUBLE, final_e, 1, MPI_DOUBLE, 
		MPI_COMM_WORLD);

  for(i=0; i<nnodes; i++) {     /* evaluate the node with the lowest score */
    if ( final_e[i] <= minyet ) {
      minyet = final_e[i];
      winner = i;
    }
  }

/* write the answer */

  if ( myid == winner ) {
#endif                                   /* outfile different from infile? */
    if ( strcmp(inname, outname) ) {           /* -> create copy of infile */
      shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));
      sprintf(shell_cmd, "cp -f %s %s", inname, outname);
      if ( -1 == system(shell_cmd) )
	error("FinalMove: error creating output file %s", outname);
      free(shell_cmd);
    }

/* all the funcs below write a section at its appropriate position in the  */
/* data file; to achieve this, they create a temporary file which is then  */
/* renamed to the final output file name                                   */

    WriteVersion(outname);
    if ( equil ) 
      WriteEquil(outname, equil_var); 
    else {
      WriteParameters(outname, parm, "eqparms", precision);
      WriteAParameters(outname, ap);
    }

#ifdef MPI
  }
#endif

/* clean up the state file and free memory */

  if ( !equil && !nofile_flag )
#ifdef MPI
    if ( ! tuning )
#endif
      StateRm();

  FreeZygote();

  FreeHistory();

  FreeExternalInputs();

}






/*** FILE FUNCTIONS FOR READING AND WRITING MISCELLANEOUS STUFF ************/

/*** InitEquilibrate: reads the equilibrate section of the data file, ******
 *                    which is needed for equilibration runs; then puts    *
 *                    the parameters into a static struct in lsa.c         *
 ***************************************************************************/

void InitEquilibrate(FILE *fp)
{
  ChuParam   l_equil_param;                    /* local equil_param struct */

  fp = FindSection(fp, "equilibrate");                /* find tune section */
  if( !fp )
    error("ReadEquilibrate: could not locate equilibrate section");
  
  fscanf(fp,"%*s\n");                         /* advance past title line 1 */
  
  if ( 1 != fscanf(fp, "%lf\n", &(l_equil_param.end_T)) )
    error("ReadEquilibrate: error reading equilibration temperature");

  fscanf(fp,"%*s\n");                         /* advance past title line 2 */

  if ( 1 != fscanf(fp, "%d\n", &(l_equil_param.fix_T_skip)) )
    error("ReadEquilibrate: error reading fixed temperature skip");

  fscanf(fp,"%*s\n");                         /* advance past title line 3 */

  if ( 1 != fscanf(fp, "%d\n", &(l_equil_param.fix_T_step)) )
    error("ReadEquilibrate: error reading fixed temperature step");

  SetEquilibrate(l_equil_param);
}



/*** ReadTune: reads the tune_parameters section in a data file and ********
 *             turns a SAType structure to the caller                      *
 ***************************************************************************/

SAType ReadTune(FILE *fp)
{
  SAType in_tune;              /* this is the SAType struct that we return */

  int intbuf1;         /* following four are used as temp. buffer for ints */
  int intbuf2;
  int intbuf3;
  int intbuf4;

  double dbuf1;     /* following four are used as temp. buffer for doubles */
  double dbuf2;
  double dbuf3;
  double dbuf4;



  fp = FindSection(fp, "tune_parameters");            /* find tune section */
  if( !fp )
    error("ReadTune: could not locate tune_parameters");

  fscanf(fp,"%*s\n");                         /* advance past title line 1 */

                                                  /* read Lam lambda stuff */
  if ( 4!= (fscanf(fp,"%lg %lg %lg %lg\n", &dbuf1, &dbuf2, &dbuf3, &dbuf4)) )
    error ("ReadTune: error reading Lam lambda stuff");

  in_tune.lambda              = dbuf1;
  in_tune.lambda_mem_length_u = dbuf2;
  in_tune.lambda_mem_length_v = dbuf3;
  in_tune.control             = dbuf4;

  fscanf(fp,"%*s\n");                         /* advance past title line 2 */

                        /* read initial moves, tau, freeze count and Sskip */
  if ( 4 !=(fscanf(fp,"%d %d %d %d\n",&intbuf1,&intbuf2,&intbuf3,&intbuf4)))
    error ("ReadTune: error reading line 2 of tune parameters");

  in_tune.initial_moves = intbuf1;
  in_tune.tau           = intbuf2;
  in_tune.freeze_count  = intbuf3;
  in_tune.update_S_skip = intbuf4;

  fscanf(fp,"%*s\n");                         /* advance past title line 3 */

                                                    /* read stop criterion */
  if ( 1 != (fscanf(fp,"%lg\n", &dbuf1)) )
    error ("ReadTune: error reading stop criterion");

  in_tune.criterion = dbuf1;

#ifdef MPI
  fscanf(fp,"%*s\n");                         /* advance past title line 4 */

  if ( 1 != fscanf(fp, "%d\n", &intbuf1) )
    error("ReadTune: error reading mixing interval");

  in_tune.mix_interval = intbuf1;
#endif
  
  in_tune.tunefile  = NULL;                  /* We're not using this stuff */
  in_tune.debuglevel = 0;
  
  return (in_tune);
}



/*** ReadAParameters: reads the AParm struct from an annealing_input sec- **
 *                    tion; these are the annealing parameters that are    *
 *                    not Lam-specific (and should NOT go into lsa.c       *
 ***************************************************************************/

AParms ReadAParameters(FILE *fp)
{
  AParms     l_aparms;                              /* local Aparms struct */

  fp = FindSection(fp, "annealing_input");            /* find tune section */
  if( !fp )
    error("ReadAParameters: could not locate tune_parameters");
  
  fscanf(fp,"%*s\n");                         /* advance past title line 1 */
  
  if ( 1 != fscanf(fp, "%li\n", &(l_aparms.seed)) )
    error("ReadAParameters: error reading random number generator seed");

  fscanf(fp,"%*s\n");                         /* advance past title line 2 */

  if ( 1 != fscanf(fp, "%lf\n", &(l_aparms.start_tempr)) )
    error("ReadAParameters: error reading start temperature");

  fscanf(fp,"%*s\n");                         /* advance past title line 3 */

  if ( 1 != fscanf(fp, "%lf\n", &(l_aparms.gain)) )
    error("ReadAParameters: error reading gain");

  fscanf(fp,"%*s\n");                         /* advance past title line 4 */

  if ( 1 != fscanf(fp, "%d\n", &(l_aparms.interval)) )
    error("ReadAParameters: error reading interval");

  return l_aparms;
}




/*** WriteAParameters: writes the aparm struct into a new section in the ***
 *                     file specified by filename; the new 'annealing_out- *
 *                     put section is inserted right after the 'tune_para- *
 *                     meters' section; to achieve this, we need to write  *
 *                     to a temporary file which is then renamed to the    *
 *                     output file name                                    *
 ***************************************************************************/

void WriteAParameters(char *filename, AParms aparm)
{
  char   *temp;                                     /* temporary file name */
  char   *record;                         /* record to be read and written */
  char   *record_ptr;        /* pointer used to remember record for 'free' */
  char   *saverec;                 /* used to save following section title */
  char   *shell_cmd;                             /* used by 'system' below */

  FILE   *outfile;                                  /* name of output file */
  FILE   *tmpfile;                               /* name of temporary file */
    

  temp      = (char *)calloc(MAX_RECORD, sizeof(char));
  record    = (char *)calloc(MAX_RECORD, sizeof(char));
  saverec   = (char *)calloc(MAX_RECORD, sizeof(char));
  shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));

  record_ptr = record;            /* this is to remember record for 'free' */

/* open output and temporary file */

  outfile = fopen(filename, "r");              /* open outfile for reading */
  if ( !outfile )                              /* sorry for the confusion! */
    error("WriteAParameters: error opening %s", filename);

  temp = strcpy(temp,"aparmXXXXXX");              /* required by mkstemp() */
  if ( mkstemp(temp) == -1 )              /* get unique name for temp file */
    error("WriteAParameters: error creating temporary file name");

  tmpfile = fopen(temp, "w");               /* ... and open it for writing */
  if ( !tmpfile )
    error("WriteAParameters: error opening temporary file %s", temp);

  if ( FindSection(outfile, "annealing_output") ) { 
    fclose(outfile);                   /* erase section if already present */
    KillSection(filename, "annealing_output");          
    outfile = fopen(filename, "r");        
  }
  rewind(outfile);

/* the follwoing three loops look for the appropriate file position to     */
/* write the annealing_output section                                      */

  while ( strncmp(record=fgets(record, MAX_RECORD, outfile), 
		  "$tune_parameters", 16) )
      fputs(record, tmpfile);
  fputs(record, tmpfile);

  while ( strncmp(record=fgets(record, MAX_RECORD, outfile), "$$", 2) )
    fputs(record, tmpfile);
  fputs(record, tmpfile);

  do {
    record = fgets(record, MAX_RECORD, outfile);
    if ( !record ) break;
  } while ( strncmp(record, "$", 1) );

  fputs("\n", tmpfile);

  if ( record )
    saverec = strcpy(saverec, record);

/* now we write the eqparms section into the tmpfile */

  PrintAParameters(tmpfile, aparm, "annealing_output");

  fprintf(tmpfile, "\n");

/* ... and then write all the rest */

  if ( record )
    fputs(saverec, tmpfile);

  while ( (record=fgets(record, MAX_RECORD, outfile)) )
    fputs(record, tmpfile);          
  
  fclose(outfile);
  fclose(tmpfile);
  
/* rename tmpfile into new file */

  sprintf(shell_cmd, "cp -f %s %s", temp, filename);

  if ( -1 == system(shell_cmd) )
    error("WriteAParameters: error renaming temp file %s", temp);

  if ( remove(temp) )
    warning("WriteAParameters: temp file %s could not be deleted", 
	    temp);

/* clean up */

  free(record_ptr);
  free(saverec);
  free(temp);
  free(shell_cmd);
}
 


/*** PrintAParameters: writes an 'annnealing_output' section with 'title' **
 *                     to the stream specified by fp                       *
 ***************************************************************************/

void PrintAParameters(FILE *fp, AParms aparm, char *title)
{

  fprintf(fp, "$%s\n", title);
  fprintf(fp, "final_energy:\n");
  fprintf(fp, "%.3f\n", aparm.stop_energy);
  fprintf(fp, "max_count:\n");
  fprintf(fp, "%d\n", aparm.max_count);
  fprintf(fp, "$$\n");

}
  
  

/*** WriteEquil: writes the equilibrate_variance section to the data file **
 *               right after the $equilibrate section                      *
 ***************************************************************************/

void WriteEquil(char *filename, double *equil_var)
{
  char   *temp;                                     /* temporary file name */
  char   *record;                         /* record to be read and written */
  char   *record_ptr;        /* pointer used to remember record for 'free' */
  char   *saverec;                 /* used to save following section title */
  char   *shell_cmd;                             /* used by 'system' below */

  FILE   *outfile;                                  /* name of output file */
  FILE   *tmpfile;                               /* name of temporary file */
    

  temp      = (char *)calloc(MAX_RECORD, sizeof(char));
  record    = (char *)calloc(MAX_RECORD, sizeof(char));
  saverec   = (char *)calloc(MAX_RECORD, sizeof(char));
  shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));

  record_ptr = record;            /* this is to remember record for 'free' */

/* open output and temporary file */

  outfile = fopen(filename, "r");              /* open outfile for reading */
  if ( !outfile )                              /* sorry for the confusion! */
    error("WriteEquil: error opening %s", filename);

  temp = strcpy(temp,"equilXXXXXX");              /* required by mkstemp() */
  if ( mkstemp(temp) == -1 )              /* get unique name for temp file */
    error("WriteEquil: error creating temporary file name");

  tmpfile = fopen(temp, "w");               /* ... and open it for writing */
  if ( !tmpfile )
    error("WriteEquil: error opening temporary file %s", temp);

  if ( FindSection(outfile, "equilibrate_variance") ) { 
    fclose(outfile);                   /* erase section if already present */
    KillSection(filename, "equilibrate_variance");          
    outfile = fopen(filename, "r");        
  }
  rewind(outfile);

/* the follwoing three loops look for the appropriate file position to     */
/* write the equilibrate_variance section                                  */

  while ( strncmp(record=fgets(record, MAX_RECORD, outfile), 
		  "$equilibrate", 12) )
      fputs(record, tmpfile);
  fputs(record, tmpfile);

  while ( strncmp(record=fgets(record, MAX_RECORD, outfile), "$$", 2) )
    fputs(record, tmpfile);
  fputs(record, tmpfile);

  do {
    record = fgets(record, MAX_RECORD, outfile);
    if ( !record ) break;
  } while ( strncmp(record, "$", 1) );

  fputs("\n", tmpfile);

  if ( record )
    saverec = strcpy(saverec, record);

/* now we write the eqparms section into the tmpfile */

  PrintEquil(tmpfile, equil_var, "equilibrate_variance");

  fprintf(tmpfile, "\n");

/* ... and then write all the rest, if there is any */

  if ( record )
    fputs(saverec, tmpfile);
  
  while ( (record=fgets(record, MAX_RECORD, outfile)) )
    fputs(record, tmpfile);          
  
  fclose(outfile);
  fclose(tmpfile);
  
/* rename tmpfile into new file */

  sprintf(shell_cmd, "cp -f %s %s", temp, filename);

  if ( -1 == system(shell_cmd) )
    error("WriteEquil: error renaming temp file %s", temp);

  if ( remove(temp) )
    warning("WriteEquil: temp file %s could not be deleted", 
	    temp);

/* clean up */

  free(record_ptr);
  free(saverec);
  free(temp);
  free(shell_cmd);
}



/*** PrintEquil: writes an 'equilibrate_variance' section with 'title' *****
 *               to the stream specified by fp                             *
 ***************************************************************************/

void PrintEquil(FILE *fp, double *equil_var, char *title)
{
  fprintf(fp, "$%s\n", title);
  fprintf(fp, "variance:\n");
  fprintf(fp, "%g\n", equil_var[0]);
  fprintf(fp, "equilibrate_final_energy:\n");
  fprintf(fp, "%g\n", equil_var[1]);
  fprintf(fp, "$$\n");
}



/*** WriteVersion: prints the version and the complete command line used ***
 *                 to run fly_sa into the $version section of the data     *
 *                 file                                                    *
 ***************************************************************************/

void WriteVersion(char *filename)
{
  char   *temp;                                     /* temporary file name */
  char   *record;                         /* record to be read and written */
  char   *record_ptr;        /* pointer used to remember record for 'free' */
  char   *convline;                 /* temporarily saves conv version line */
  char   *shell_cmd;                             /* used by 'system' below */

  FILE   *outfile;                                  /* name of output file */
  FILE   *tmpfile;                               /* name of temporary file */
    
  int    versionflag = 0;       /* version section already present or not? */



  temp      = (char *)calloc(MAX_RECORD, sizeof(char));
  record    = (char *)calloc(MAX_RECORD, sizeof(char));
  convline  = (char *)calloc(MAX_RECORD, sizeof(char));
  shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));

  record_ptr = record;            /* this is to remember record for 'free' */

  outfile = fopen(filename, "r");              /* open outfile for reading */
  if ( !outfile )                              /* sorry for the confusion! */
    error("WriteVersion: error opening %s", filename);

  if (FindSection(outfile, "version"))  /* version section already there? */
    versionflag = 1;
  rewind(outfile);                               
  
  temp = strcpy(temp,"versionXXXXXX");            /* required by mkstemp() */
  if ( mkstemp(temp) == -1 )              /* get unique name for temp file */
    error("WriteVersion: error creating temporary file");

  tmpfile = fopen(temp, "w");               /* ... and open it for writing */
  if ( !tmpfile )
    error("WriteVersion: error opening temporary file %s", temp);

  if ( !versionflag ) {                  /* CASE 1: no version section yet */
    fprintf(tmpfile, "$version\n");       /* -> write new one at beginning */
    fprintf(tmpfile, "%s\n", version);                 /* of the data file */
    fprintf(tmpfile, "%s", argvsave);
    fprintf(tmpfile, "$$\n\n");
    while ( (record = fgets(record, MAX_RECORD, outfile)) )
      fputs(record, tmpfile);
                                /* CASE 2: version section already present */
  } else {                      /* -> append new lines to existing section */
    while ( strncmp(record=fgets(record, MAX_RECORD, outfile), 
		    "$version", 8) )      /* first write everything before */
      fputs(record, tmpfile);                           /* version section */
    fputs(record, tmpfile);

    while ( strncmp(record=fgets(record, MAX_RECORD, outfile),
		    "$$", 2) )             /* save existing converter line */
      if ( !strncmp(record, "converted", 9) )         /* skip all the rest */
	convline = strcpy(convline, record);
    
    fprintf(tmpfile, "%s\n", version);       /* new lines are written here */
    fprintf(tmpfile, "%s", argvsave);

    if ( strlen(convline) > 0 )      /* conversion line is appended at end */
      fprintf(tmpfile, "%s", convline);

    fprintf(tmpfile, "$$\n");

    while ( (record=fgets(record, MAX_RECORD, outfile)) )
      fputs(record, tmpfile);           /* ... and then write all the rest */
  }
  
  fclose(outfile);
  fclose(tmpfile);

/* rename tmpfile into new file */

  sprintf(shell_cmd, "cp -f %s %s", temp, filename);

  if ( -1 == system(shell_cmd) )
    error("WriteVersion: error renaming temp file %s", temp);

  if ( remove(temp) )
    warning("WriteVersion: temp file %s could not be deleted", temp);

/* clean up */
  free(temp);
  free(record_ptr);
  free(convline);
  free(shell_cmd);
}



/*** WriteTimes: writes the time-structure to a .times file ****************
 ***************************************************************************/

void WriteTimes(double *times)
{
  char   *timefile;                             /* name of the .times file */
  FILE   *timeptr;                         /* file pointer for .times file */
  
           /* create time file name by appending .times to input file name */
  timefile = (char *)calloc(MAX_RECORD, sizeof(char));
  timefile = strcpy(timefile, outname);
  timefile = strcat(timefile, ".times");

  timeptr = fopen(timefile, "w");
  if ( !timeptr )
    file_error("main");

  PrintTimes(timeptr, times);                /* write times to .times file */

  fclose(timeptr);                                             /* clean up */
  free(timefile);
}



/*** PrintTimes: writes two (parallel: three) times sections ***************
 ***************************************************************************/

void PrintTimes(FILE *fp, double *times)
{
  fprintf(fp, "wallclock: %.3f\n", times[0]);
  fprintf(fp, "user:      %.3f\n", times[1]);
}





/*** FUNCTIONS THAT COMMUNICATE WITH SAVESTATE.C ***************************/

/*** GetOptions: returns command line options to savestate.c ***************
 *               for the detailed meaning of all these options see Parse-  *
 *               CommandLine() above); Opts struct defined in moves.h      *
 ***************************************************************************/

Opts *GetOptions(void)
{  
  Opts       *options;
  
  options = (Opts *)malloc(sizeof(Opts));

  options->inname  = inname;
  options->outname = outname;
  options->argv    = argvsave;

  options->derivfunc = (char *)calloc(MAX_RECORD, sizeof(char));

  if ( pd == DvdtOrig )
    options->derivfunc = strcpy(options->derivfunc, "DvdtOrig");
  else
    error("GetOptions: unknown derivative function");

  options->solver    = (char *)calloc(MAX_RECORD, sizeof(char));
  
  if ( ps == Adams )
    options->solver = strcpy(options->solver, "Adams");
  else if ( ps == BaDe )
    options->solver = strcpy(options->solver, "BaDe");
  else if ( ps == BuSt )
    options->solver = strcpy(options->solver, "BuSt");
  else if ( ps == Euler )
    options->solver = strcpy(options->solver, "Euler");
  else if ( ps == Heun )
    options->solver = strcpy(options->solver, "Heun");
  else if ( ps == Milne )
    options->solver = strcpy(options->solver, "Milne");
  else if ( ps == Meuler )
    options->solver = strcpy(options->solver, "Meuler");
  else if ( ps == Rk2 )
    options->solver = strcpy(options->solver, "Rk2");
  else if ( ps == Rk4 )
    options->solver = strcpy(options->solver, "Rk4");
  else if ( ps == Rkck )
    options->solver = strcpy(options->solver, "Rkck");
  else if ( ps == Rkf )
    options->solver = strcpy(options->solver, "Rkf");
  else if ( ps == SoDe )
    options->solver = strcpy(options->solver, "SoDe");
  else
    error("GetOptions: unknown solver function");

  options->stop_flag   = stop_flag;
  options->prolix_flag = prolix_flag;
  options->landscape_flag = landscape_flag;
  options->log_flag    = log_flag;
  options->time_flag   = time_flag;
  options->state_write = state_write;
  options->print_freq  = print_freq;
  options->captions    = captions;
  options->olddivstyle = olddivstyle;
  options->precision   = precision;
  options->stepsize    = stepsize;
  options->quenchit    = quenchit;
  options->equil       = equil;

#ifdef MPI
  options->tuning          = tuning;
  options->covar_index     = covar_index;
  options->write_tune_stat = write_tune_stat;
  options->auto_stop_tune  = auto_stop_tune;
#endif

  return options;
}



/*** RestoreOptions: restores the values of the command line opt variables *
 *                   from the Opts struct (used for restoring a run)       *
 ***************************************************************************/

void RestoreOptions(Opts *options)
{

/* restore input/output file names and the full command line string; note  *
 * that the output file name needs to be communicated to lsa.c, since we   * 
 * need it there for setting up log and tuning file names                  */

  inname  = options->inname;
  outname = options->outname;
  SetOutname(outname); 
  argvsave = options->argv;

/* all the other options */

  if ( !strcmp(options->derivfunc, "DvdtOrig") )
    pd = DvdtOrig;
  else
    error("RestoreOptions: unknown deriv function %s", options->derivfunc);

  if ( !strcmp(options->solver, "Adams") )
    ps = Adams;
  else if ( !strcmp(options->solver, "BaDe") )
    ps = BuSt;
  else if ( !strcmp(options->solver, "BuSt") )
    ps = BuSt;
  else if ( !strcmp(options->solver, "Euler") )
    ps = Euler;
  else if ( !strcmp(options->solver, "Heun") )
    ps = Heun;
  else if ( !strcmp(options->solver, "Milne") )
    ps = Milne;
  else if ( !strcmp(options->solver, "Meuler") )
    ps = Meuler;
  else if ( !strcmp(options->solver, "Rk2") )
    ps = Rk2;
  else if ( !strcmp(options->solver, "Rk4") )
    ps = Rk4;
  else if ( !strcmp(options->solver, "Rkck") )
    ps = Rkck;
  else if ( !strcmp(options->solver, "Rkf") )
    ps = Rkf;
  else if ( !strcmp(options->solver, "SoDe") )
    ps = SoDe;
  else
    error("RestoreOptions: unknown solver %s", options->solver);

  stop_flag   = options->stop_flag;

  prolix_flag = options->prolix_flag;
  if ( prolix_flag )
    SetProlix(prolix_flag, outname, 0);

  landscape_flag = options->landscape_flag;
  if ( landscape_flag )
    error("RestoreOptions: cannot restore an equilibration (lanDscape) run");

  log_flag    = options->log_flag;
  time_flag   = options->time_flag;

  state_write = options->state_write;
  print_freq  = options->print_freq;
  captions    = options->captions;

  olddivstyle = options->olddivstyle;
  precision   = options->precision;
  stepsize    = options->stepsize;

  quenchit    = options->quenchit;
  equil       = options->equil;
  if ( equil )
    error("RestoreOptions: cannot restore an equilibration run");

#ifdef MPI
  tuning          = options->tuning;
  if ( tuning )
    error("RestoreOptions: cannot restore a tuning run");
  covar_index     = options->covar_index;
  write_tune_stat = options->write_tune_stat;
  auto_stop_tune  = options->auto_stop_tune;
#endif

  free(options->derivfunc);
  free(options->solver);
  free(options);
}
