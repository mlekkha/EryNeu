/*******************************************************************
 *                                                                 *
 *   unfold.c                                                      *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 *   Version 9.3.2                                                 *
 *   written by JR, modified by Yoginho                            *
 *   -g option by Yousong Wang, Feb 2002                           *
 *   -a option by Marcel Wolf, Apr 2002                            *
 *   -G (guts) option by Manu, June 2002                           *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 * unfold.c contains main for the unfold utility. unfold runs      *
 * the fly model using the parameters and the maternal contribu-   *
 * tions stored for the specified genotype in the data file        *
 * specified by filename.                                          *
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
#include <unistd.h>                                          /* for getopt */
#include <time.h>							/* for time calculation */
#include <sys/resource.h>			/*		for time calculation */

#include <error.h>   
#include <integrate.h>
#include <maternal.h>
#include <solvers.h>
#include <zygotic.h>
#include <score.h>
 

/*** Constants *************************************************************/

#define  OPTS      ":a:Df:g:Ghi:j:op:r:s:t:vx:z:"   /* cmd line opt string */
 
/*** Help, usage and version messages **************************************/

static const char usage[]  =

"Usage: unfold [-a <accuracy>] [-D] [-f <float_prec>] [-g <g(u)>] [-G]\n"
"              [-h] [-i <stepsize>] [-j <timefile>] [-o] [-p <pstep>]\n"
"              [-s <solver>] [-t <time>] [-v] [-x <sect_title>]\n"
"              [-z <gast_time>]\n"
"              <datafile> [<genotype>]\n";

static const char help[]   =

"Usage: unfold [options] <datafile> [<genotype>]\n\n"

"Arguments:\n"
"  <datafile>          data file which contains equation parameters\n"
"  <genotype>          genotype string (W, R, S or T for each gene)\n\n"

"Options:\n"
"  -a <accuracy>       solver accuracy for adaptive stepsize ODE solvers\n"
"  -D                  debugging mode, prints all kinds of debugging info\n"
"  -f <float_prec>     float precision of output is <float_prec>\n"
"  -g <g(u)>           chooses g(u): e = exp, h = hvs, s = sqrt, t = tanh\n"
"  -G                  prints guts instead of model output\n"
"  -h                  prints this help message\n"
"  -i <stepsize>       sets ODE solver stepsize (in minutes)\n"
"  -j <timefile>       reads output fimes from <timefile>\n"
"  -o                  use oldstyle cell division times (3 div only)\n"
"  -p <pstep>          prints output every <pstep> minutes\n"
"  -s <solver>         choose ODE solver\n"
"  -t <time>           prints output for <time>\n"
"  -v                  print version and compilation date\n"
"  -x <sect_title>     uses equation paramters from section <sect_title>\n\n"
"  -z <gast_time>      set custom gastrulation time (max. 10'000'000 (min));\n"
"                      custom gastrulation MUST be later than the normal one\n"

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
  char            *infile;                   /* pointer to input file name */
  FILE            *fp;                       /* pointer to input data file */
 
  int             c;                 /* used to parse command line options */
  int             i,ii;                                       /* loop counter */
  int		  numguts;           /* number of guts elements calculated */

  char            *dumpfile;             /* name of debugging model output */
  FILE            *dumpptr;        /* pointer for dumping raw model output */

  char            *slogfile;                    /* name of solver log file */
  FILE            *slog;                        /* solver log file pointer */

/* the follwoing few variables are read as arguments to command line opts */

  int             ndigits  = 6;        /* precision of output, default = 6 */
  int             guts     = 0;                  /* flag for printing guts */

  double          stepsize   = 1.;                  /* stepsize for solver */
  double          accuracy   = 0.001;              /* accuracy for solvers */
  double          p_stepsize = 0.;                      /* output stepsize */

  double          time = -999999999.;                /* time for -t option */
  char            *timefile = NULL;                  /* file for -j option */

  char            *section_title;                /* parameter section name */

/* stuff used as input/output for blastoderm */

  NArrPtr         answer;                /* model output is stored in this */
  NArrPtr         outtab;                            /* times to be output */
  NArrPtr	  goutput;	       /* output of guts is stored in this */
  DArrPtr         tt;           /* array with times for which there's data */

  int			   genindex = 0;			/* genotype index */
  char            *genotype  = NULL;                    /* genotype string */
  Slist			  *geno, *curr; 		/* We will read in the genotypes to
  set facts in score.c */
  char 		  **gutsdefs = NULL; /* array of strings to hold guts defs */
  char            *gutstitle = NULL;
  int             length;                     /* length of genotype string */
 struct rusage 	begin,end;		/*	structs for measuring time */

 DataTable 		**interp_dat; /* for providing
										 history to blastoderm */
 double 			  *theta_discons;
 int 			  theta_discons_size;
 double		   	  *temp_divtable;
 double		      *temp_durations;
 InterpObject 	  *polation;
 InterpObject 	  *extinp_polation;

/* the following lines define a pointers to:                               */
/*            - pd:    dvdt function, currently only DvdtOrig in zygotic.c */
/*            - pj:    Jacobian function, in zygotic.c                     */
/*                                                                         */
/* NOTE: ps (solver) is declared as global in integrate.h                  */

  void (*pd)(double *, double, double *, int);
  void (*pdd)(double *, double **, double, double *, int);
  void (*pj)(double, double *, double *, double **, int);

/* external declarations for command line option parsing (unistd.h) */

  extern char *optarg;                     /* command line option argument */
  extern int  optind;                /* pointer to current element of argv */
  extern int  optopt;               /* contain option character upon error */

/* following part sets default values for deriv, Jacobian and solver funcs */
/* and section title                                                       */

  getrusage (RUSAGE_SELF, &begin);	/*		get start time */
  
  pd = DvdtOrig;
  pdd = DvdtDelay;
  pj = JacobnOrig;
  ps = Rk4;

  section_title = (char *)calloc(MAX_RECORD, sizeof(char));
  section_title = strcpy(section_title, "eqparms");  /* default is eqparms */

/* following part parses command line for options and their arguments      */
/* modified from original getopt manpage                                   */
 
  optarg = NULL;
  while ( (c = getopt(argc, argv, OPTS)) != -1 )
    switch (c) {
    case 'a':
      accuracy = atof(optarg);
      if ( accuracy <= 0 )
	error("unfold: accuracy (%g) is too small", accuracy);
      break;
    case 'D':                                 /* -D runs in debugging mode */
      debug = 1;
      break;
    case 'f':
      ndigits = atoi(optarg);             /* -f determines float precision */
      if ( ndigits < 0 )
	error("unfold: what exactly would a negative precision be???");
      if ( ndigits > MAX_PRECISION )
	error("unfold: max. float precision is %d!", MAX_PRECISION);
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
	error("unfold: %s is an invalid g(u), should be e, h, s or t",
	      optarg); 
     break;
    case 'G':                                           /* -G guts options */
      guts = 1;
      break;
    case 'h':                                            /* -h help option */
      PrintMsg(help, 0);
      break;
    case 'i':                           /* -i sets stepsize for the solver */
      stepsize = atof(optarg);
      if ( stepsize < 0 )
	error("unfold: going backwards? (hint: check your -i)");
      if ( stepsize == 0 ) 
	error("unfold: going nowhere? (hint: check your -i)");
      if ( stepsize > MAX_STEPSIZE )
	error("unfold: stepsize %g too large (max. is %g)", 
	      stepsize, MAX_STEPSIZE);
      break;
    case 'j':
      if ( p_stepsize > 0. )
	error("unfold: can't use -p and -j together!");
      if ( time != -999999999. )
	error("unfold: can't use -t and -j together!");
      timefile = (char *)calloc(MAX_RECORD, sizeof(char));
      timefile = strcpy(timefile, optarg);
      break;
    case 'o':             /* -o sets old division style (ndivs = 3 only! ) */
      olddivstyle = 1;
      break;
    case 'p':                         /* -p sets print stepsize for output */
      if ( timefile )
	error("unfold: can't use -j and -p together!");
      if ( time != -999999999. )
	error("unfold: can't use -t and -p together!");
      p_stepsize = atof(optarg);
      if ( p_stepsize < 0.001 ) 
	error("unfold: output stepsize (%g) too small (min 0.001)",
	      p_stepsize);
      break;
    case 'r':
      error("unfold: -r is not supported anymore, use -g instead");
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
	error("unfold: invalid solver (%s), use: a,bd,bs,e,h,mi,me,r{2,4,ck,f}",
	      optarg);
      break;
    case 't':
      if ( timefile )
	error("unfold: can't use -j and -t together!");
      if ( p_stepsize > 0. )
	error("unfold: can't use -p and -t together!");
      time = atof(optarg);
      if ( (time < 0) && (time != -999999999))
	error("unfold: the time (%g) doesn't make sense", time);
      break;
    case 'v':                                  /* -v prints version message */
      fprintf(stderr, verstring, *argv, 
	      VERS, USR, MACHINE, COMPILER, FLAGS, __DATE__, __TIME__);
      exit(0);
    case 'x':
      if ( (strcmp(optarg, "input")) && (strcmp(optarg, "eqparms")) &&
	   (strcmp(optarg, "parameters")) )
	error("unfold: invalid section title (%s)", optarg);
      section_title = strcpy(section_title, optarg);
      break;
    case 'z':
      custom_gast = atof(optarg);
      if ( custom_gast < 0. )
	error("unfold: gastrulation time must be positive");
      if ( custom_gast > 10000000. )
	error("unfold: gastrulation time must be smaller than 10'000'000");
      break;
    case ':':
      error("unfold: need an argument for option -%c", optopt);
    case '?':
    default:
      error("unfold: unrecognized option -%c", optopt);
    }

/* error check */
  
  if ( stepsize > p_stepsize && p_stepsize != 0 )
    error("unfold: print-stepsize (%g) smaller than stepsize (%g)!",
	  p_stepsize, stepsize);

  if ( (argc - (optind - 1)) < 2 || (argc - (optind - 1)) > 4 )
     PrintMsg(usage, 1);

/* let's get started and open the data file here */

  infile = argv[optind];
  fp = fopen(infile, "r");      
  if( !fp )
    file_error("unfold");
 
/* Initialization code here */
  
  if ( (optind+1) < argc ) 
  	genindex = atoi(argv[optind+1]);

  if ( (optind+2) < argc)  {
    genotype = (char *)calloc(MAX_RECORD, sizeof(char));
    genotype = strcpy(genotype, argv[optind+2]);        /* genotype string */
  }

  InitZygote(fp, pd, pj, section_title); 
                                     /* initializes parameters, defs, bias */
                                     /* and bicoid and defines static vars */
                                     /* for the dvdt function              */

/* initialize genotype if necessary, otherwise check for errors */

  if ( !(genotype) ) {
    genotype = (char *)calloc(MAX_RECORD, sizeof(char));
    for (i=0; i<defs.ngenes; i++)
      genotype = strcat(genotype, "W");    /* construct wt genotype string */
  } else {
    length = (int)strlen(genotype);
    if ( length != defs.ngenes )
      error("unfold: genotype length doesn't match number of genes (%d)", 
	    defs.ngenes);
    for (i=0; i<defs.ngenes; i++) {
      c = (int)*(genotype+i);
      if ( ( c != 87 ) && (( c < 82 ) || ( c > 84 )) )
	error("unfold: genotype string can only contain R, S, T or W");
    }
  }

/* Now we will read and traverse the genotypes list and find the
facts section tilte for the genotype and set the facts list in
integrate.c */

	geno = ReadGenotypes(fp);
	curr = geno;

	if ( !(genindex < count_Slist(geno)))
		error("Unfold: genindex more than last allele!\n");

	i = 0;
	while (i++ < genindex) curr = curr -> next;

    if ( !(interp_dat=(DataTable **)calloc(1, 
                                            sizeof(DataTable *))) )
        error("Unfold: could not allocate interp_dat struct");

    if ( !(polation=(InterpObject *)calloc(1, sizeof(InterpObject))) )
        error("Unfold: could not allocate polation struct");
    

    /* Initialize the history is we are using the delay solver */
    if (ps == SoDe)
    {

        
        GetInterp(fp, curr->hist_section, defs.ngenes, interp_dat);
            
        DoInterp(*interp_dat, polation, defs.ngenes);
        FreeFacts(*interp_dat);

        theta_discons = Get_Theta_Discons(&theta_discons_size);


        polation->fact_discons = (double *) realloc(polation->fact_discons,
            (polation->fact_discons_size+theta_discons_size)*sizeof(double));

        for (ii=0; ii < theta_discons_size; ii++)
            polation->fact_discons[polation->fact_discons_size + ii] =
                                                    theta_discons[ii];

        polation->fact_discons_size += theta_discons_size;
        free(theta_discons);

        if ( defs.ndivs > 0 ) {
            if ( !(temp_divtable = GetDivtable()) )
                error("Unfold: error getting temp_divtable");
            if ( !(temp_durations = GetDurations()) )
                error("Unfold: error getting division temp_durations");
          
            for (ii=0; ii<defs.ndivs; ii++) { 
                polation->fact_discons = (double *) realloc(polation->fact_discons,
                                (polation->fact_discons_size+4)*sizeof(double));
        
                polation->fact_discons[polation->fact_discons_size] =
                                    temp_divtable[ii];
                polation->fact_discons[polation->fact_discons_size+1] =
                                    temp_divtable[ii] + EPSILON;
                polation->fact_discons[polation->fact_discons_size+2] =
                                    temp_divtable[ii]-temp_durations[ii];
                polation->fact_discons[polation->fact_discons_size+3] =
                                    temp_divtable[ii]-temp_durations[ii] + EPSILON;

                polation->fact_discons_size += 4;
            }
        }


        qsort((void *) polation->fact_discons, 
            polation->fact_discons_size, 
            sizeof(double),(int (*) (void*,void*)) compare);
    }
/* Below we prepare the interpolant for external inputs */		

	if ( !(extinp_polation=(InterpObject *)calloc(1, sizeof(InterpObject))) )
		error("Unfold: could not allocate extinp_polation struct");

    if (defs.egenes > 0) {
		
		GetInterp(fp, curr->ext_section, defs.egenes, interp_dat);
			
		DoInterp(*interp_dat, extinp_polation, defs.egenes);

		FreeFacts(*interp_dat);
	}


	free_Slist(geno);

	
  fclose(fp);
 

/* initialize debugging file names */

  if ( debug ) {

    slogfile = (char *)calloc(MAX_RECORD, sizeof(char));
    dumpfile = (char *)calloc(MAX_RECORD, sizeof(char));

    sprintf(slogfile, "%s.slog", infile);
    sprintf(dumpfile, "%s.%s.%d.uout", infile, genotype, genindex);

    slog = fopen(slogfile, "w");              /* delete existing slog file */
    fclose(slog);

    slog = fopen(slogfile, "a");            /* now keep open for appending */
 } 

/* the following establishes the tabulated times according to either:      */
/* -t (output single time step)                                            */
/* -j (read times from file)                                               */
/* -p (stepsize set by argument)                                           */
/* default behavior is producing output for gastrulation time only         */

  if ( time != -999999999. ) {
    if ( time > GetGastTime() )
      error("unfold: time (%g) larger than gastrulation time!", time);
    tt.size  = 1;
    tt.array = (double *)calloc(1,sizeof(double));
    tt.array[0] = time;
  } else if ( timefile != NULL ) {
    tt = ReadTimes(timefile);
    free(timefile);
  } else if ( p_stepsize != 0. ) {
    tt = MakeTable(p_stepsize);   
  } else {
    tt.size  = 1;
    tt.array = (double *)calloc(1,sizeof(double));
    tt.array[0] = GetGastTime();
  }

/* Run the model... */

  answer = Blastoderm(genindex, genotype, *polation, *extinp_polation, tt, stepsize, accuracy, slog);

/* if debugging: print out the innards of the model to unfold.out */

  if ( debug ) {
    dumpptr = fopen(dumpfile, "w");      
    if( !dumpptr ) {
      perror("unfold");
      exit(1);
    }
    PrintBlastoderm(dumpptr, answer, "debug_output", MAX_PRECISION, 
		    defs.ngenes);
    fclose(dumpptr);
  }
                                              
/* strip output of anything that's not in tt */
   
  outtab = ConvertAnswer(answer, tt);

/* code for printing guts... first read gutsdefs section */

  if ( guts ) { 
 		
    fp = fopen(infile, "r");      
    if( !fp ) {
      perror("unfold");
      exit(1);
    }
    gutsdefs = ReadGuts(fp);	
    fclose(fp);

/* Make sure that interp_info in integrate.c is pointing at the
 * external inputs InterpObject */
 
	if (defs.egenes > 0)
		SetExternalInputInterp(*extinp_polation);

    gutstitle = (char *)calloc(MAX_RECORD, sizeof(char));

/* then calculate guts and print them */

    for (i=0; *(gutsdefs+i); i++) {
      if (numguts = CalcGuts(genindex, genotype, *polation, *extinp_polation, outtab, &goutput, *(gutsdefs+i))) { 
	gutstitle = strcpy(gutstitle, "guts_for_");
	PrintBlastoderm(stdout, goutput, strcat(gutstitle, *(gutsdefs+i)), 
		       ndigits, numguts); 
	FreeSolution(&goutput);
      }
      free(*(gutsdefs+i));
    }
    free(gutsdefs);		
    free(gutstitle);

/* code for printing model output */

  } else 
    PrintBlastoderm(stdout, outtab, "output\n", ndigits, defs.ngenes);
                                                 
/* ... and then clean up before you go home. */

  getrusage (RUSAGE_SELF, &end);	 /*	get end time */
  printf("Unfold ran for %.13f seconds\n", tvsub(end, begin));

  FreeSolution(&answer);

  FreeZygote();

  /* Free stuff allocated for delay solver */
  if (ps == SoDe)
  {    
      FreeInterpObject(polation);
  }  

  free(polation);
  free(interp_dat);


  if (defs.egenes > 0) {
	  FreeInterpObject(extinp_polation);
  }	  
  free(extinp_polation);

  free(outtab.array);
  free(tt.array);

  free(section_title);
  free(genotype);

  if ( debug ) {
    fclose(slog);
    free(dumpfile);
    free(slogfile);
  }

  return 0;
}
