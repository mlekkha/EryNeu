/*******************************************************************
 *                                                                 *
 *   printscore.c                                                  *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 *   written by JR, modified by Yoginho                            *
 *   -g option by Yousong Wang, Feb 2002                           *
 *   -a option by Marcel Wolf, Apr 2002                            *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 * printscore.c contains main() for the printscore utility.        *
 *                                                                 *
 * printscore prints score and penalty, as well as the root mean   *
 * square (RMS) of a gene circuit to STDOUT                        *
 *                                                                 *
 * See the dataformatX.X file for further details on the current   *
 * data file format (X.X stands for the current code version).     *
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>                                          /* for getopt */

#include <error.h>
#include <integrate.h>
#include <maternal.h>
#include <score.h>
#include <solvers.h> 
#include <zygotic.h>


/*** Constants *************************************************************/

#define  OPTS      ":a:Df:g:Ghi:j:k:opqr:s:vx:"  /* command line option string */


/*** Help, usage and version messages **************************************/

static const char usage[]  =

"Usage: printscore [-a <accuracy>] [-D] [-f <float_prec>] [-g <g(u)>] [-G]\n"
"                  [-h] [-i <stepsize>] [-k <cost_function>] [-o] [-p] [-s <solver>] [-v]\n"
"                  [-x <sect_title>]\n"
"                  <datafile>\n";

static const char help[]   =

"Usage: printscore [options] <datafile>\n\n"

"Arguments:\n"
"  <datafile>          data file for which we evaluate score and RMS\n\n"

"Options:\n"
"  -a <accuracy>       solver accuracy for adaptive stepsize ODE solvers\n"
"  -D                  debugging mode, prints all kinds of debugging info\n"
"  -f <float_prec>     float precision of output is <float_prec>\n"
"  -g <g(u)>           chooses g(u): e = exp, h = hvs, s = sqrt, t = tanh\n"
"  -G                  gut mode: prints squared diffs for all datapoints\n"
"  -h                  prints this help message\n"
"  -i <stepsize>       sets ODE solver stepsize (in minutes)\n"
"  -k <cost_function> select error function for calculation the differences between\n"
"                      equation and data: sq=squared differences, abs=absolute differences\n"
"                      default=sq\n"
"  -o                  use oldstyle cell division times (3 div only)\n"
"  -p                  prints penalty in addition to score and RMS\n"
"  -s <solver>         choose ODE solver\n"
"  -v                  print version and compilation date\n"
"  -x <sect_title>     uses equation paramters from section <sect_title>\n\n"

"Please report bugs to <yoginho@usa.net>. Thank you!\n";

static const char verstring[] = 

"%s version %s\n"
"compiled by:      %s\n"
"         on:      %s\n"
"      using:      %s\n"
"      flags:      %s\n"
"       date:      %s at %s\n";


static char   *cost_function = "sq"; /* default cost function

/*** Main program **********************************************************/

int main(int argc, char **argv)
{
  int           c;                   /* used to parse command line options */
  FILE          *infile;                     /* pointer to input data file */

  char          *dumpfile;               /* name of debugging model output */
  FILE          *dumpptr;          /* pointer for dumping raw model output */

  char          *slogfile;                      /* name of solver log file */
  FILE          *slog;                          /* solver log file pointer */

/* the follwoing few variables are read as arguments to command line opts  */

  int           ndigits     = 12;     /* precision of output, default = 12 */
  int           gutndigits  = 6;      /* precision for gut output, def = 6 */
  int           penaltyflag = 0;              /* flag for printing penalty */
  int           rmsflag     = 1;     /* flag for printing root mean square */
  int           gutflag     = 0;         /* flag for root square diff guts */

  double        stepsize    = 1.;                   /* stepsize for solver */
  double        accuracy    = 0.001;                /* accuracy for solver */
  char          *cost_function = "sq"; /* error function utilized in error
                                         * evaluation in Eval method      */

  char          *section_title;                  /* parameter section name */

/* two format strings */

  char          *format;                           /* output format string */
  char          *precision;   /* precision string for output format string */

/* output values */

  double        chisq       = 0.; /* sum of squared differences for output */

  double        penalty     = 0.;                  /* variable for penalty */

  double        rms         = 0.;                /* root mean square (RMS) */
  double        ms          = 0.;             /* mean square (= rms * rms) */
  int           ndp         = 0;                   /* number of datapoints */
 
/* the following lines define a pointers to:                               */
/*            - pd:    dvdt function, currently only DvdtOrig in zygotic.c */
/*            - pj:    Jacobian function, in zygotic.c                     */
/*                                                                         */
/* NOTE: ps (solver) is declared as global in integrate.h                  */

  void (*pd)(double *, double, double *, int);
  void (*pj)(double, double *, double *, double **, int);

/* external declarations for command line option parsing (unistd.h) */

  extern char   *optarg;                   /* command line option argument */
  extern int    optind;              /* pointer to current element of argv */
  extern int    optopt;             /* contain option character upon error */


/* following part sets default values for deriv, Jacobian and solver funcs */

  pd = DvdtOrig;
  pj = JacobnOrig;
  ps = Rk4;
  
  eval = SquareEval; /* default evaluation function */

  section_title = (char *)calloc(MAX_RECORD, sizeof(char));
  section_title = strcpy(section_title, "eqparms");  /* default is eqparms */

/* following part parses command line for options and their arguments      */
  
  optarg = NULL;
  while ( (c = getopt(argc, argv, OPTS)) != -1 ) 
    switch (c) {
    case 'a':
      accuracy = atof(optarg);
      if ( accuracy <= 0 )
	error("fly_sa: accuracy (%g) is too small", accuracy);
      break;
    case 'D':                                 /* -D runs in debugging mode */
      debug = 1;
      break;
    case 'f':
      ndigits    = atoi(optarg);          /* -f determines float precision */
      gutndigits = atoi(optarg);
      if ( ndigits < 0 )
	error("printscore: what exactly would a negative precision be???");
      if ( ndigits > MAX_PRECISION )
	error("printscore: max. float precision is %d!", MAX_PRECISION);
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
	error("printscore: %s is an invalid g(u), should be e, h, s or t",
	      optarg); 
      break;
    case 'G':                                              /* -G guts mode */
      gutflag = 1;
      break;
    case 'h':                                            /* -h help option */
      PrintMsg(help, 0);
      break;
    case 'i':                                      /* -i sets the stepsize */
      stepsize = atof(optarg);
      if ( stepsize < 0 )
	error("printscore: going backwards? (hint: check your -i)");
      if ( stepsize == 0 ) 
	error("printscore: going nowhere? (hint: check your -i)");
      if ( stepsize > MAX_STEPSIZE )
	error("printscore: stepsize %g too large (max. is %g)", 
	      stepsize, MAX_STEPSIZE);
      break;
    case 'k':       
      cost_function = (char *)calloc(MAX_RECORD, sizeof(char)); 
      cost_function = strcpy(cost_function, optarg);
      if ( !(strcmp(optarg, "sq")) )
          eval = SquareEval;
      if ( !(strcmp(optarg, "abs")) )
          eval = AbsoluteEval;
      break;
    case 'o':             /* -o sets old division style (ndivs = 3 only! ) */
      olddivstyle = 1;
      break;
    case 'p':
      penaltyflag = 1;
      break;
    case 'q':
      warning("printscore: RMS now printed by default, no need for -q");
      break;
    case 'r':
      error("printscore: -r is not supported anymore, use -g instead");
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
	error("printscore: bad solver (%s), use: a,bd,bs,e,h,mi,me,r{2,4,ck,f}",
	      optarg);
      break;
    case 'v':                                  /* -v prints version number */
      fprintf(stderr, verstring, *argv, 
	      VERS, USR, MACHINE, COMPILER, FLAGS, __DATE__, __TIME__);
      exit(0);
    case 'x':
      if ( (strcmp(optarg, "input")) && (strcmp(optarg, "eqparms")) &&
	   (strcmp(optarg, "parameters")) )
	error("unfold: invalid section title (%s)", optarg);
      section_title = strcpy(section_title, optarg);
      break;
    case ':':
      error("printscore: need an argument for option -%c", optopt);
      break;
    case '?':
    default:
      error("printscore: unrecognized option -%c", optopt);
    }
  
/* error check */

  if ( (argc - (optind - 1 )) != 2 )
    PrintMsg(usage, 1);
  
/* dynamic allocation of output format strings */

  precision = (char *)calloc(MAX_RECORD, sizeof(char));
  format    = (char *)calloc(MAX_RECORD, sizeof(char));

/* initialize guts */

  if ( gutflag ) 
    SetGuts(gutflag, gutndigits);
 
/* let's get started and open data file here */

  infile = fopen(argv[optind],"r");
  if ( !infile )
    file_error("printscore");
  
/* if debugging: initialize solver log file */

  if ( debug ) {

    slogfile = (char *)calloc(MAX_RECORD, sizeof(char));
    sprintf(slogfile, "%s.slog", argv[optind]);

    slog = fopen(slogfile, "w");              /* delete existing slog file */
    fclose(slog);

    slog = fopen(slogfile, "a");            /* now keep open for appending */
 } 

/* Initialization code here: InitZygote() initializes everything needed    *
 * for running the model, InitScoring() everything for scoring; InitStep-  *
 * size sets solver stepsize and accuracy in score.c to pass it on to      *
 * Blastoderm                                                              */

  InitZygote(infile, pd, pj, section_title);      
  InitScoring(infile);          

  InitHistory(infile);

  InitExternalInputs(infile);

  fclose(infile);       
  InitStepsize(stepsize, accuracy, slog, argv[optind], cost_function); 

/* Scoring happens here (and also printing of guts if -G) */

  chisq = Score();  

  if ( gutflag ) {              /* exit here if we've printed guts already */
    FreeZygote();

    FreeHistory();

	FreeExternalInputs();

    free(precision);
    free(format);
    free(section_title);
    return 0;
  }

/* next few lines format & print the output to the appropriate precision   */

  sprintf(precision, "%d", ndigits);     
                                                       
  format = strcpy(format, " chisq = %.");   
  format = strcat(format, precision);        
  format = strcat(format, "f");     
                                             
  printf(format, chisq);           

  if ( penaltyflag ) {                 /* in case of -p, print penalty too */

    penalty = GetPenalty();

    if ( penalty == -1 ) {
      printf("    explicit limits !");
    } else {
      format = strcpy(format, "     penalty = %.");
      format = strcat(format, precision);
      format = strcat(format, "f");
      
      printf(format, penalty);
    }
  }
  
  if ( rmsflag ) {                                            /* print rms */

    ndp = GetNDatapoints();
    penalty = GetPenalty();

    ms  = (chisq - penalty) / (double)ndp;

    if ( !(strcmp(cost_function, "sq")) ){
        rms = sqrt(ms);
    }else if ( !(strcmp(cost_function, "abs")) ){
        rms = ms;
    }
    
    format = strcpy(format, "     rms = %.");
    format = strcat(format, precision);
    format = strcat(format, "f");

    printf(format, rms);
  }

  printf("\n");

/* clean up before you go home... */

  FreeZygote();

  FreeHistory();

  FreeExternalInputs();

  free(precision);
  free(format);
  free(section_title);

  return 0;
}
