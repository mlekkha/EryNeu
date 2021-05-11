/*******************************************************************
 *                                                                 *
 *   newtonraphson.c                                                      *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 *   Version 0.1	                                               *
 *   written by Manu, adapted from newtonraphson.c				   *	
 *                                                                 *
 *******************************************************************
 *                                                                 *
 * It calculates the 1D unstable manifolds of the saddle points    *
 * of one-nucleus circuits of the fly equations. 				   *	
 * It takes the datafile and the stationary points file (output    *
 * of newtonraphson) and outputs the trajectories of the one-      *
 * dimensional manifolds to standard output. It reads in the saddle*
 * point (x_0) specified on the command line and its associated    *
 * eigenvectors and calculates trajectories starting from x_0      *
 * perturbed by +-epsilon in the direction of the eigenvector of   *
 * the positive eigenvalue.									 		   *
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

#include <float.h> 								/* for DBL_EPSILON */
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <unistd.h>                                          /* for getopt */
#include <math.h>
#include <time.h>							/* for time calculation */
#include <sys/resource.h>			/*		for time calculation */

#include <error.h>   
#include <integrate.h>
#include <maternal.h>
#include <solvers.h>
#include <zygotic.h>
#include <score.h>


/*** Constants *************************************************************/

#define  OPTS      ":a:e:f:hi:n:p:s:t:T:vx:"   /* cmd line opt string */
 
/*** Help, usage and version messages **************************************/

static const char usage[]  =

"Usage: newtonraphson  [-f <float_prec>] [-g <grid_start>] \n"
"              [-G <grid_end>] [-h] [-I <max_iterations>] \n"
"              [-s <grid size>] [-t <external input time>] \n"
"			   [-T <tolerance>] [-v] [-x <sect_title>] <datafile> \n";

static const char help[]   =

"Usage: newtonraphson [options] <datafile> \n\n"

"Arguments:\n"
"  <datafile>          data file which contains equation parameters\n"
"Options:\n"
"  -f <float_prec>     float precision of output is <float_prec>\n"
"  -g <grid_start>     Starting value for uniform grid\n"
"  -G <grid_end>       Ending value for uniform grid\n"
"  -h                  prints this help message\n"
"  -I <max_iterations> Maximum number of newton raphson iterations\n"
"  -s <grid_size>      Number of points on the grid including "
"                      start and end points\n"
"  -t <extinput time>  Choose which time to pick external input"
"                      values from\n"
"  -T <tolerance>      Tolerance for newton raphson\n"
"  -v                  print version and compilation date\n"
"  -x <sect_title>     uses equation paramters from section <sect_title>\n\n"

"Please report bugs to <manu@odd.bio.sunysb.edu>. Thank you!\n";

static const char verstring[] = 

"%s version %s\n"
"compiled by:      %s\n"
"         on:      %s\n"
"      using:      %s\n"
"      flags:      %s\n"
"       date:      %s at %s\n";

static const char settings[] = 

"$settings.%d\n"
"Input file: %s\n"
"Stationary points file %s\n"
"Time point from which external input values taken %f\n"
"Time interval [%f,%f]\n"
"Increment %f\n"
"Solver stephint/stepsize %f\n"
"Accuracy %f\n"
"Stationary point #%d\n"
"$$\n\n";

/* Static global variables in newtonraphson.c, for rhs() and
 * Jacobian() */  

  static double *R;                    /* strength of each promoter--always >= 0. */
  static double *T;                            /* the genetic interconnect matrix */
  static double *E;                            /* the external input
  regulatory matrix */
  static double *m;                      /* regulatory coefficient of bcd on gene */
  static double *h;           /* reg. coeff. for generic TFs on synthesis of gene */
  static double *lambda;                       /*protein half lives--always >= 0. */
  static int	ngenes;			/* number of genes */
  static int	egenes;			/* number of external inputs */
  static double	*v_ext;		/* to hold external input concentrations */
  static double	v_bcd;		/* to hold bcd concentration */
  static double *bot, *bot2, *vinput; 	/* intermediates in rhs() */


  static double	*xin;			/* starting value of newtonraphson */
  static double	*xout;			/* result from newtonraphson 	*/

/* function declarations */

void rhs(double *x, double t, double *y, int n);
void ReadStatpoints(FILE *fp, char *section_title, double *point);
void ReadEvs(FILE *fp, char *section_title, double *vectoresr, 
												double *vectoresi);



/*** Main program **********************************************************/
 
int main(int argc, char **argv)
{
  char            *infile;                   /* pointer to input file name */
  char            *statptfile;       /* pointer to stationary points file name */
  FILE            *fp;                       /* pointer to input data file */
 
  int             c;                 /* used to parse command line options */
  int             num_nucs;                 /* number of nuclei */
  int             i,ii,jj,kk;                      /* loop counter */
  int 			  length;			/* to check genotype's validity */	
  double		  grid_increment;

/* the follwoing few variables are read as arguments to command line opts */

  int			  statpt_no=-1;		/* stationary point no. of the
									  saddle */
  int             ndigits  = 14;        /* precision of output, default = 6 */
  double          stepsize   = 4.;                  /* stepsize for solver */
  double          accuracy   = 0.001;              /* accuracy for solvers */
  double		  tinitial = 0.0;			/* always 0 */
  double		  tfinal = 1000.0;			/* the time until which
								  the unstable manifold is calculated */
  double		  tin, tout;		/*starting and ending times for
										  call to the solver */
  double          p_stepsize = 1.;              /* time increment */
  double		  epsilon = 0.1;				/* perturbation to
												  stattionary point */

  double		  extinput_time = 67.975;		/* time for which the
  external inputs value should be taken */

  char            *section_title;                /* parameter section name */

  int			   genindex = 0;			/* genotype index */
  char            *genotype  = NULL;                    /* genotype string */
  Slist			  *geno, *curr; 		/* We will read in the genotypes to
  set facts in score.c */
 struct rusage 	begin,end;		/*	structs for measuring time */

 DataTable 		**interp_dat; /* for providing
										 history to blastoderm */
 InterpObject 	  *extinp_polation;
 DArrPtr		bcd_concs;		/* to get the bicoid from maternal.c */
 EqParms		  parms;

 double 		*stationary_point=NULL;	/*stat pt from file */
 double			*evr,*evi;					/* eigenvectors		*/
 double			*vector;				/* point to the right ev */



/* external declarations for command line option parsing (unistd.h) */

  extern char *optarg;                     /* command line option argument */
  extern int  optind;                /* pointer to current element of argv */
  extern int  optopt;               /* contain option character upon error */

/* following part sets default values for deriv, Jacobian and solver funcs */
/* and section title                                                       */


  ps = BuSt;
  p_deriv = rhs;

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
	error("unstable-manifold: accuracy (%g) is too small", accuracy);
      break;
  case 'e':
      epsilon = atof(optarg);
      if ( epsilon <= 0 )
	error("unstable-manifold: epsilon (%g) is too small", epsilon);
      break;
  case 'f':
      ndigits = atoi(optarg);    
      if ( ndigits < 0 )
	error("newtonraphson: what exactly would a negative precision be???");
      if ( ndigits > MAX_PRECISION )
	error("newtonraphson: max. float precision is %d!", MAX_PRECISION);
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
     case 'n':                         /* which statpoint to pick
					 from 0 to N-1 where the number of points is N */
      statpt_no = atoi(optarg);
      if ( statpt_no < 0 ) 
	error("unstable-manifold: invalid stattionary point no.%d",
														statpt_no);
      break;
    case 'p':                         /* -p sets print stepsize for output */
      p_stepsize = atof(optarg);
      if ( p_stepsize < 0.001 ) 
	error("unstable-manifold: increments (%g) too small (min 0.001)",
	      p_stepsize);
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
      extinput_time = atof(optarg);  
      if ( extinput_time < 0 )
	error("newtonraphson: the external input time (%g) doesn't make sense", extinput_time);
     break;
    case 'T':
      tfinal = atof(optarg);             
      if ( tfinal < DBL_EPSILON )
	error("newtonraphson: the tfinal (%g) is too small", tfinal);
     break;
     case 'v':                         /* -v prints version message */
      fprintf(stderr, verstring, *argv, 
	      VERS, USR, MACHINE, COMPILER, FLAGS, __DATE__, __TIME__);
      exit(0);
   case 'x':
      if ( (strcmp(optarg, "input")) && (strcmp(optarg, "eqparms")) &&
	   (strcmp(optarg, "parameters")) )
	error("newtonraphson: invalid section title (%s)", optarg);
      section_title = strcpy(section_title, optarg);
      break;
    case ':':
      error("newtonraphson: need an argument for option -%c", optopt);
      PrintMsg(usage, 1);
	  break;
    case '?':
    default:
      error("newtonraphson: unrecognized option -%c", optopt);
      PrintMsg(usage, 1);
	  break;
    }


  if ( optind != argc - 2 )
     PrintMsg(usage, 1);

  if (statpt_no < 0)
  	 error("No stationary point specified via -n\n");

  if (tfinal - tinitial <= p_stepsize)
  	 error("Increment (%g) too big for total duration (%g)",
	 							p_stepsize, tfinal-tinitial);

/* let's get started and open the data file here */

  infile = argv[optind];
  statptfile = argv[optind + 1];

  fp = fopen(infile, "r");      
  if( !fp )
    file_error("unstable-manifold");
 
/* Initialization code here */

/* Below are the things that normally (in unfold for eg.) are done by
 * InitZygote(). Params and defs are needed for newton-raphson, while
 * the others are needed to read the external input section */

	defs = ReadTheProblem(fp);     /* initializes defs */
	parms = ReadParameters(fp, section_title);  /* initialize parameters */	
	InitBicoid(fp);				
	InitBias(fp);
	InitNNucs();
	

/* initialize genotype if necessary, otherwise check for errors */

  if ( !(genotype) ) {
    genotype = (char *)calloc(MAX_RECORD, sizeof(char));
    for (i=0; i<defs.ngenes; i++)
      genotype = strcat(genotype, "W");    /* construct wt genotype string */
  } else {
    length = (int)strlen(genotype);
    if ( length != defs.ngenes )
      error("newtonraphson: genotype length doesn't match number of genes (%d)", 
	    defs.ngenes);
    for (i=0; i<defs.ngenes; i++) {
      c = (int)*(genotype+i);
      if ( ( c != 87 ) && (( c < 82 ) || ( c > 84 )) )
	error("newtonraphson: genotype string can only contain R, S, T or W");
    }
  }


/* Now we will read and traverse the genotypes list and find the
external inputs title and set the interpolation structure in
integrate.c. The genindex is set to its default value 0, so we are
only doing the first genotype here */

	geno = ReadGenotypes(fp);
	curr = geno;

	if ( !(genindex < count_Slist(geno)))
		error("newtonraphson: genindex more than last allele!\n");

	i = 0;
	while (i++ < genindex) curr = curr -> next;

/* Below we prepare the interpolant for external inputs */		

  	if ( !(interp_dat=(DataTable **)calloc(1, 
  											sizeof(DataTable *))) )
    	error("newtonraphson: could not allocate interp_dat struct");
	

  	if ( !(extinp_polation=(InterpObject *)calloc(1, sizeof(InterpObject))) )
    	error("newtonraphson: could not allocate extinp_polation struct");
	
	GetInterp(fp, curr->ext_section, defs.egenes, interp_dat);
		
  	DoInterp(*interp_dat, extinp_polation, defs.egenes);

	FreeFacts(*interp_dat);

	free_Slist(geno);
	
  fclose(fp);
 
  SetExternalInputInterp(*extinp_polation);

  num_nucs = GetNNucs(extinput_time);
  
  v_ext = (double *) calloc(num_nucs*defs.egenes, sizeof(double));
  ExternalInputs(extinput_time, extinput_time, v_ext, num_nucs*defs.egenes);

  bcd_concs = GetBicoid(extinput_time, genindex);
  v_bcd = bcd_concs.array[0];	

  /* assign params to static variables so rhs() can see them.
   * v_bcd and v_ext already done.
   * Make static various variables and allocate arrays 
   * newtonraphson(), traverse_grid() can see them.  */
   
/* for rhs() and Jac() */   

  R = parms.R;
  T = parms.T;
  E = parms.E;
  lambda = parms.lambda;
  h = parms.h;
  m = parms.m;

  ngenes = defs.ngenes;
  egenes = defs.egenes;

  bot = (double *) calloc(ngenes, sizeof(double));
  bot2 = (double *) calloc(ngenes, sizeof(double));
  vinput = (double *) calloc(ngenes, sizeof(double));

/* for traverse_grid() and newtonraphson () */

  xin = (double *) calloc(ngenes, sizeof(double));
  xout = (double *) calloc(ngenes, sizeof(double));

  stationary_point = (double *) calloc(ngenes, sizeof(double));
  evr = (double *) calloc(ngenes*ngenes, sizeof(double));
  evi = (double *) calloc(ngenes*ngenes, sizeof(double));

  /* lets read in the stationary point and associated eigenvectors */

  fp = fopen(statptfile, "r");      
  if( !fp )
    file_error("unstable-manifold");
	
  sprintf(section_title, "zero.%d", statpt_no);
  ReadStatpoints(fp, section_title, stationary_point);

  sprintf(section_title, "eigenvectors.%d", statpt_no);
  ReadEvs(fp, section_title, evr, evi);

  fclose(fp);

/* print the settings default and from commandline */

  printf(settings,statpt_no,infile,statptfile,extinput_time,
			  tinitial,tfinal,p_stepsize,stepsize,accuracy,statpt_no);


 /* pick the first eigenvector - associated with the largest
  * eigenvalue, and hence the only positive one. This is because the
  * newtonraphson orders the eigenvalues and eigenvectors in this way.
  * To be added: read eigenvalues also and find the position of the
  * positive one.*/

  vector = evr;

/* The first stating point, xin = stationary points + 
 *						epsilon*(eigenvector of positive eigenvalue) */

 for (ii=0; ii < ngenes; ii++)
 	xin[ii] = vector[ii]*epsilon + stationary_point[ii];

 tin = tinitial;
 tout = tinitial+p_stepsize;
 memcpy(xout, xin, sizeof(*xin) * ngenes);
 
 printf("$um_1.%d\n",statpt_no);
 while (tout <= tfinal) {

	printf(" %9.3f",tin); 
	for (ii = 0; ii < ngenes; ii++)
		printf(" %*.*f",ndigits+5,ndigits,xout[ii]); 
	printf("\n");

	(*ps)(xin, xout, tin, tout, stepsize, accuracy, ngenes, stderr);

	tin = tout;
	memcpy(xin, xout, sizeof(*xout) * ngenes);
	tout += p_stepsize;
 }
 printf("$$\n\n");

/* The second stating point, xin = stationary points - 
 *						epsilon*(eigenvector of positive eigenvalue) */

 for (ii=0; ii < ngenes; ii++)
 	xin[ii] = stationary_point[ii] - vector[ii]*epsilon;

 tin = tinitial;
 tout = tinitial+p_stepsize;
 memcpy(xout, xin, sizeof(*xin) * ngenes);
 
 printf("$um_2.%d\n",statpt_no);
 while (tout <= tfinal) {

	printf(" %9.3f",tin); 
	for (ii = 0; ii < ngenes; ii++)
		printf(" %*.*f",ndigits+5,ndigits,xout[ii]); 
	printf("\n");

	(*ps)(xin, xout, tin, tout, stepsize, accuracy, ngenes, stderr);

	tin = tout;
	memcpy(xin, xout, sizeof(*xout) * ngenes);
	tout += p_stepsize;
 }
 printf("$$\n\n");

  FreeInterpObject(extinp_polation);
  free(extinp_polation);
  free(interp_dat);
  free(v_ext);
  free(bot);
  free(bot2);
  free(vinput);
  free(xin);
  free(xout);
  free(stationary_point);
  free(evr);
  free(evi);

  free(section_title);
  free(genotype);


  return 0;
}

/* Right hand side of a one-nucleus circuit */

void rhs(double *x, double t, double *y, int n)
{

  int        j;                                   /* local loop counter */
  int        k;                   /* index of gene k */
  double 	 vinput1 = 0;     		/* accumulators */
  double 	 vdot1, g1;
  



		for (k=0; k < ngenes; k++) {

		  
		  vinput1  = h[k];
		  vinput1 += m[k] * v_bcd;   

		  for(j=0; j < egenes; j++)
			vinput1 += E[(k*egenes)+j] * v_ext[j];

		  for(j=0; j < ngenes; j++)
			vinput1 += T[(k*ngenes)+j] * x[j];

		  bot2[k]   = 1 + vinput1 * vinput1;
		  vinput[k] = vinput1;
		}

    for(k=0; k < ngenes; k++)                  /* slow traditional style sqrt */
		bot[k] = sqrt(bot2[k]);

            /* next loop does the rest of the equation (R, Ds and lambdas) */
                                                 /* store result in y[] */



	  for( k=0; k<ngenes; k++ ) { 
	  
	    vdot1   = -lambda[k] * x[k];
	    g1      = 1 + vinput[k]/bot[k];         
	    vdot1  += R[k] * 0.5 * g1;     
	    y[k] = vdot1; 
	  }
	
	  
	  
	return;
}

void ReadStatpoints(FILE *fp, char *section_title, double *point)
{

  int               i;                               /* local loop counter */
  int               c;                         /* used to parse text lines */
  int               lead_punct;            /* flag for leading punctuation */

  char              *base;          /* pointer to beginning of line string */
  char              *record;    /* string for reading whole line of params */

  char              **fmt;   /* array of format strings for reading params */
  char              *skip;               /* string of values to be skipped */

  const char        read_fmt[] = "%lg";                   /* read a double */
  const char        skip_fmt[] = "%*lg ";               /* ignore a double */

  base = (char *)calloc(MAX_RECORD, sizeof(char));

  skip = (char *)calloc(MAX_RECORD, sizeof(char));
  fmt  = (char **)calloc(ngenes, sizeof(char *));

/* create format strings according to the number of genes */

  for ( i=0; i<ngenes; i++ ) {  
    fmt[i] = (char *)calloc(MAX_RECORD, sizeof(char));   
    fmt[i] = strcpy(fmt[i], skip);
    fmt[i] = strcat(fmt[i], read_fmt);
    skip   = strcat(skip, skip_fmt);
 } 

  fp = FindSection(fp, section_title);             /* find eqparms section */
  if( !fp )
    error("ReadStatpoints: cannot locate %s section", section_title);

  while ( strncmp(( base=fgets(base, MAX_RECORD, fp)), "$$", 2)) {
    
    record = base;
    lead_punct = 0;                                           /* of string */

    c = (int)*record;
    while ( c != '\0' ) {
      
      if ( isdigit(c) ) {                            /* line contains data */
	record = base;

/* usually read ngenes parameters, but for diff. schedule A or C only read *
 * one d parameter                                                         */

	  for ( i=0; i < ngenes; i++ ) {
		if ( 1 != sscanf(record, fmt[i], point+i) )
		  error("ReadStatpoints: error reading stationary points");
	  }

	break;                           /* don't do rest of loop anymore! */
      } 
      else if ( isalpha(c) ) {                     /* letter means comment */
	break;
      }
      else if ( c == '-' ) {                  /* next two elsifs for punct */
	if ( ((int) *(record+1)) == '.' )
	  record++;
	lead_punct = 1;
	c=(int)*(++record);
      }
      else if ( c == '.' ) {
	lead_punct = 1;
	c=(int)*(++record);
      }
      else if ( ispunct(c) ) {                /* other punct means comment */
	break;
      }
      else if ( isspace(c) ) {               /* ignore leading white space */
	if ( lead_punct )                 /* white space after punct means */
	  break;                                              /* comment */
	else {
	  c=(int)*(++record);              /* get next character in record */
	}
      }	
      else {
	error("ReadStatpoints: illegal character in %s");
      }
    }
  }
  
  free(base);
  free(skip);
  
  for (i=0; i<ngenes; i++)
    free(fmt[i]);
  free(fmt);

  return;
}

void ReadEvs(FILE *fp, char *section_title, double *vectoresr, 
												double *vectoresi)
{

  int               i;                               /* local loop counter */
  int               c;                         /* used to parse text lines */
  int				rowno = 0;			/*the row # you are on */
  int               lead_punct;            /* flag for leading punctuation */

  char              *base;          /* pointer to beginning of line string */
  char              *record;    /* string for reading whole line of params */

  char              **fmt;   /* array of format strings for reading params */
  char              *skip;               /* string of values to be skipped */

  const char        read_fmt[] = "%lg%lg%*c";                   /* read a double */
  const char        skip_fmt[] = "%*lg%*lg%*c ";               /* ignore a double */

  base = (char *)calloc(MAX_RECORD, sizeof(char));

  skip = (char *)calloc(MAX_RECORD, sizeof(char));
  fmt  = (char **)calloc(ngenes, sizeof(char *));

/* create format strings according to the number of genes */

  for ( i=0; i<ngenes; i++ ) {  
    fmt[i] = (char *)calloc(MAX_RECORD, sizeof(char));   
    fmt[i] = strcpy(fmt[i], skip);
    fmt[i] = strcat(fmt[i], read_fmt);
    skip   = strcat(skip, skip_fmt);
 } 

  fp = FindSection(fp, section_title);             /* find eqparms section */
  if( !fp )
    error("ReadEvs: cannot locate %s section", section_title);

  while ( strncmp(( base=fgets(base, MAX_RECORD, fp)), "$$", 2)) {
    
    record = base;
    lead_punct = 0;                                           /* of string */

    c = (int)*record;
    while ( c != '\0' ) {
      
      if ( isdigit(c) ) {                            /* line contains data */
	record = base;

/* Read in the eigenvectors. The imaginary parts will not be used as
 * this program concerns itself only with stationary points that have
 * one positive real eigenvalue and the rest real negative ones. Also,
 * since in the output of newtonraphson the eigenvectors are printed
 * column-wise, we transpose them while reading them in.*/

	  if (rowno < ngenes) {
		  for ( i=0; i < ngenes; i++ ) {
			if ( 2 != 
				sscanf(record, fmt[i], &vectoresr[i*ngenes+rowno],
										&vectoresi[i*ngenes+rowno]) 
			)
			  error("ReadEvs: error reading eigenvectors");
		  }
	  
		  rowno++;
	  }	else 
	  	error("Too many rows in the eigenvector matrix\n");

	break;                         /* don't do rest of loop anymore! */
      } 
      else if ( isalpha(c) ) {                     /* letter means comment */
	break;
      }
      else if ( c == '-' ) {                  /* next two elsifs for punct */
	if ( ((int) *(record+1)) == '.' )
	  record++;
	lead_punct = 1;
	c=(int)*(++record);
      }
      else if ( c == '.' ) {
	lead_punct = 1;
	c=(int)*(++record);
      }
      else if ( ispunct(c) ) {                /* other punct means comment */
	break;
      }
      else if ( isspace(c) ) {               /* ignore leading white space */
	if ( lead_punct )                 /* white space after punct means */
	  break;                                              /* comment */
	else {
	  c=(int)*(++record);              /* get next character in record */
	}
      }	
      else {
	error("ReadStatpoints: illegal character in %s");
      }
    }
  }
  
  free(base);
  free(skip);
  
  for (i=0; i<ngenes; i++)
    free(fmt[i]);
  free(fmt);

  if (rowno < ngenes - 1)
  	error("Not enough rows in the eigenvector matrix\n");

  return;
}
