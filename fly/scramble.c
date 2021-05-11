/*******************************************************************
 *                                                                 *
 *   scramble.c                                                    *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 *   written by JR, modified by Yoginho                            *
 *                                                                 * 
 *******************************************************************
 *                                                                 *
 * scramble.c contains main for the scramble program. This is a    *
 * standalone program that reads the 'input' or 'eqparms' sec-     * 
 * tion of a data file and returns random parameter values that    * 
 * are within the limits specified by the $limits section of the   *
 * data file. Note that only those parameters will get scrambled   *
 * which are indicated as to-be-tweaked in the $tweak section of   *
 * the data file.                                                  *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 * A short comment on penalty limits (JJ, Aug 7, 2001):            *
 *                                                                 *
 * In order to produce random values of T, m and h that are        *
 * within reasonable limits when no explicit limits are used for   *
 * them, we do the following approximation:                        *
 *                                                                 *
 * First, we determine the limits of u in g(u) within which        *
 * there is no penalty assigned. This is true for all g(u) be-     *
 * tween Lambda and (1-Lambda). Therefore, we calculate the in-    *
 * verse function of g(u) for those to values to get upper and     * 
 * lower limits for u (based on ideas by Eric Mjolsness).          * 
 * Since we need to sum up all T * vmax, m * mmax and h to get     *
 * u, we'll compensate by dividing the limits by the sqrt of the   *
 * number of genes in the problem. This way we think, we'll get    *
 * reasonable limits for single parameters (idea by JR).           *
 * All the above happens in the Penalty2Limits function in sco-    *
 * re.c. When comparing parameters to limits, don't forget to      *
 * multiply Ts with vmax and ms with mmax. This happens in main    *
 * of scramble below. hs are compared as they are.                 *
 *                                                                 *
 * See JJs lab notes for further detail on g(u)-inverse and such.  *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include <moves.h>                                            /* for tweak */
#include <error.h>
#include <maternal.h>                                          /* for defs */
#include <score.h>                                           /* for limits */
#include <zygotic.h>                                        /* for EqParms */


#define  OPTS       ":f:hvw:x:"              /* command line option string */



/*** Help, usage and version messages **************************************/

static const char usage[]  =

"Usage: scramble [-f <float_prec>] [-h] [-v]\n"
"                [-w <out_file>] [-x <sect_title>]\n"
"                <datafile>\n";

static const char help[]   =

"Usage: scramble [options] <datafile>\n\n"

"Argument:\n"
"  <datafile>          data file to be scambled\n\n"

"Options:\n"
"  -f <float_prec>     float precision of output is <float_prec>\n"
"  -h                  prints this help message\n"
"  -v                  print version and compilation date\n"
"  -w <out_file>       write output to <out_file> instead of <datafile>\n"
"  -x <sect_title>     scrambles eq params from section <sect_title>\n\n"

"Please report bugs to <yoginho@usa.net>. Thank you!\n";

static const char verstring[] = 

"%s version %s\n"
"compiled by:      %s\n"
"         on:      %s\n"
"      using:      %s\n"
"      flags:      %s\n"
"       date:      %s at %s\n";




/*** Main program **********************************************************/

int main(int argc, char **argv)
{
  int         c;                     /* used to parse command line options */
  FILE        *infile;                   /* parameter file to be scrambled */

  char        *section;        /* section title of section to be scrambled */
  char        *outname = NULL;                  /* output file name for -o */
  int         ndigits = 8;             /* precision of output, default = 8 */

  int         i,j;                                  /* local loop counters */

  EqParms     param;                             /* local parameter struct */
  SearchSpace *limits;                              /* local limits struct */
  Tweak       tweak;                                 /* local tweak struct */

  long        pid;         /* process ID of scramble: used to seed drand48 */
  
  int         penaltyflag = 0;

  double      out;                 /* temporary storage for random numbers */
  double      T_pen;      /* temp var for T * vmax in case penalty is used */
  double      m_pen;      /* temp var for m * mmax in case penalty is used */

  char        *shell_cmd;                        /* used by 'system' below */

/* external declarations for command line option parsing (unistd.h) */

  extern char *optarg;                     /* command line option argument */
  extern int  optind;                /* pointer to current element of argv */
  extern int  optopt;               /* contain option character upon error */
  

/* initialize section to 'input' as default */

  section = (char *)calloc(MAX_RECORD, sizeof(char));
  section = strcpy(section, "input");          /* input section is default */

/* following part parses command line for options and their arguments      */
  
  optarg = NULL;
  while ( (c = getopt(argc, argv, OPTS)) != -1 ) 
    switch (c) {
    case 'f':
      ndigits = atoi(optarg);             /* -f determines float precision */
      if ( ndigits > MAX_PRECISION )
	error("scramble: max. float precision is %d!", MAX_PRECISION);
      break;
    case 'h':                                            /* -h help option */
      PrintMsg(help, 0);
      break;
    case 'v':                                  /* -v prints version number */
      fprintf(stderr, verstring,
	      *argv, VERS, USR, MACHINE, COMPILER, FLAGS, __DATE__, __TIME__);
      exit(0);
    case 'w':                                       /* -o sets output file */
      outname = calloc(MAX_RECORD + 1, sizeof(char)); 
      outname = strcpy(outname, optarg);
      break;
    case 'x':
      if ( (strcmp(optarg, "input")) && (strcmp(optarg, "eqparms")) )
	error("scramble: invalid section title (%s), must be input or eqparms",
	      optarg);
      section = strcpy(section, optarg);
      break;
    case ':':
      error("scramble: need an argument for option -%c", optopt);
      break;
    case '?':
    default:
      error("scramble: unrecognized option -%c", optopt);
    }

/* error check */  

  if(argc - (optind - 1) != 2)
    PrintMsg(usage, 1);

/* now, let's get started */

  infile = fopen(argv[optind], "r");
  if( !infile )
    file_error("scramble");
  
  defs     = ReadTheProblem(infile);          /* read problem (for ngenes) */
  param    = ReadParameters(infile, section); /* read the input parameters */
  tweak    = ReadTweak(infile);                  /* read the tweak section */

  InitScoring(infile);          /* initializes Limits & penalty in score.c */

  fclose(infile);

  limits   = GetLimits();         /* get the limits struct (incl. penalty) */

  if ( limits->pen_vec != 0 ) {                          /* using penalty? */
    Penalty2Limits(limits);                  /* convert to explicit limits */
    penaltyflag = 1;                             /* see also comment above */
  }

  pid = getpid();                  /* get the process ID for seeding erand */
  srand48(pid);                                            /* seed drand48 */
  

/* the following part of the code assigns each parameter-to be tweaked a   */
/* random value within the limits specified by the limits struct using the */
/* formula:                                                                */
/*             new_value = lowerlim + rand_real * (upperlim - lowerlim)    */
/*                                                                         */
/* where rand_real is a random number between 0 and 1 (drand48())          */
/* (for the penalty case, see comment above)                               */


  for ( i=0; i<defs.ngenes; i++ ) {

    if ( tweak.Rtweak[i] == 1 ) {                      /* scramble Rs here */
      out        = drand48();
      param.R[i] = limits->Rlim[i]->lower + 
	out * (limits->Rlim[i]->upper - limits->Rlim[i]->lower);
    }

    for ( j=0; j<defs.ngenes; j++ ) {                  /* scramble Ts here */
      if ( tweak.Ttweak[(i*defs.ngenes)+j] == 1 ) {
	out   = drand48();
	if ( !penaltyflag ) {
	  param.T[(i*defs.ngenes)+j] = 
	    limits->Tlim[(i*defs.ngenes)+j]->lower +
	    out * (limits->Tlim[(i*defs.ngenes)+j]->upper - 
		   limits->Tlim[(i*defs.ngenes)+j]->lower);
	} else {
	  T_pen = limits->Tlim[(i*defs.ngenes)+j]->lower +
	    out * (limits->Tlim[(i*defs.ngenes)+j]->upper -
		   limits->Tlim[(i*defs.ngenes)+j]->lower);
	  param.T[(i*defs.ngenes)+j] = T_pen / limits->pen_vec[j+2];
	}                                      /* above is T_pen / vmax[i] */
      }
    }

    for ( j=0; j<defs.egenes; j++ ) {                  /* scramble Es here */
      if ( tweak.Etweak[(i*defs.egenes)+j] == 1 ) {
	out   = drand48();
	if ( !penaltyflag ) {
	  param.E[(i*defs.egenes)+j] = 
	    limits->Elim[(i*defs.egenes)+j]->lower +
	    out * (limits->Elim[(i*defs.egenes)+j]->upper - 
		   limits->Elim[(i*defs.egenes)+j]->lower);
	} else {
	  T_pen = limits->Elim[(i*defs.egenes)+j]->lower +
	    out * (limits->Elim[(i*defs.egenes)+j]->upper -
		   limits->Elim[(i*defs.egenes)+j]->lower);
	  param.E[(i*defs.egenes)+j] = T_pen /
	  								limits->pen_vec[defs.ngenes+j+2];
	}                                      /* above is T_pen / vmax[i] */
      }
    }

    if ( tweak.mtweak[i] == 1 ) {                      /* scramble ms here */
      out = drand48();
      if ( !penaltyflag ) {
	param.m[i] = limits->mlim[i]->lower + 
	  out * (limits->mlim[i]->upper - limits->mlim[i]->lower);
      } else {
	m_pen = limits->mlim[i]->lower + 
	  out * (limits->mlim[i]->upper - limits->mlim[i]->lower);
	param.m[i] = m_pen / limits->pen_vec[1];           /* m_pen / mmax */
      }
    }

    if ( tweak.htweak[i] == 1 ) {                      /* scramble hs here */
      out        = drand48();
      param.h[i] = limits->hlim[i]->lower +
	out * (limits->hlim[i]->upper - limits->hlim[i]->lower);
    }
    
    if ( tweak.lambdatweak[i] == 1 ) {            /* scramble lambdas here */
      out             = drand48();
      param.lambda[i] = limits->lambdalim[i]->lower + 
	out * (limits->lambdalim[i]->upper - limits->lambdalim[i]->lower);
    }

    if ( tweak.tautweak[i] == 1 ) {                      /* scramble
	delays here */
      out        = drand48();
      param.tau[i] = limits->taulim[i]->lower +
	out * (limits->taulim[i]->upper - limits->taulim[i]->lower);
    }
  }

/* ds need to be srambled separately, since for diffusion schedules A & C  */
/* there's just one single d                                               */

  if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) {
    if ( tweak.dtweak[0] == 1 ) {
      out        = drand48();
      param.d[0] = limits->dlim[0]->lower +
	out * (limits->dlim[0]->upper - limits->dlim[0]->lower);
    }
  } else {
    for ( i=0; i<defs.ngenes; i++ ) {
      if ( tweak.dtweak[i] == 1 ) {
	out        = drand48();
	param.d[i] = limits->dlim[i]->lower + 
	  out * (limits->dlim[i]->upper - limits->dlim[i]->lower);
      }
    }
  }

  if ( outname != NULL ) {                                     /* -o used? */
    shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));
    sprintf(shell_cmd, "cp -f %s %s", argv[optind], outname);
    if ( -1 == system(shell_cmd) )
      error("scramble: error creating the output file %s", outname);
    free(shell_cmd);
    WriteParameters(outname, &param, section, ndigits);    
  } else {
    WriteParameters(argv[optind], &param, section, ndigits);
  }

/* clean up */

  free(section);

  free(param.R);
  free(param.T);

  if (defs.egenes > 0)
	  free(param.E);

  free(param.m);
  free(param.h);
  free(param.d);
  free(param.lambda);
  free(param.tau);

  free(tweak.Rtweak);
  free(tweak.Ttweak);

  if (defs.egenes > 0)
	  free(tweak.Etweak);

  free(tweak.mtweak);
  free(tweak.htweak);
  free(tweak.dtweak);
  free(tweak.lambdatweak);
  free(tweak.tautweak);

  for (i=0;i<defs.ngenes;i++) {
    free(limits->Rlim[i]);
    for (j=0;j<defs.ngenes;j++) 
      free(limits->Tlim[(i*defs.ngenes)+j]);
    for (j=0;j<defs.egenes;j++) 
      free(limits->Elim[(i*defs.egenes)+j]);
    free(limits->mlim[i]);
    free(limits->hlim[i]);
    free(limits->lambdalim[i]);
    free(limits->taulim[i]);
  }
  if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) 
    free(limits->dlim[0]);
  else
    for (i=0;i<defs.ngenes;i++)
      free(limits->dlim[i]);

  if ( limits->pen_vec )
    free(limits->pen_vec);
  free(limits->Rlim);
  free(limits->Tlim);

  if (defs.egenes > 0)
	  free(limits->Elim);

  free(limits->mlim);
  free(limits->hlim);
  free(limits->dlim);
  free(limits->lambdalim);
  free(limits->taulim);
  
  free(limits);

  if ( outname )
    free(outname);

  return 0;
}
