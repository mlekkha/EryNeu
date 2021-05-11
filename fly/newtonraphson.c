/*******************************************************************
 *                                                                 *
 *   newtonraphson.c                                                      *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 *   Version 0.1	                                               *
 *   written by Manu, adapted from unfold.c						   *	
 *                                                                 *
 *******************************************************************
 *                                                                 *
 * It calculates the stationary points of the fly equations  	   *
 * by brute force newton raphson. It samples uniformly spaced      *
 * starting points on a m^N grid where m is the number of points   *
 * in one axis and N is the number of genes. For time-varying      *
 * external inputs (non-autonomous) circuits the values for the    *
 * inputs to be used is specified by giving the time on the        *
 * command line. It will only for single-nucleus circuits 		   *
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

#define  OPTS      ":f:g:G:hI:s:t:T:vx:"   /* cmd line opt string */
#define  MAX_COND  1e3			  /* to decide if the jacobian is
							badly conditioned during newton-raphson */	
 
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

"$settings\n"
"Input file: %s\n"
"Time point from which external input values taken %f\n"
"grid size %d^%d\n"
"Grid: [%f,%f]\n"
"Precision %d\n"
"Tolerance %f\n"
"Max iterations %d\n"
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
  static double	*jacobian;		/* jacobian for eigenvalues */

  static double **stationary_points=NULL;	/*stat pts found */
  static double **didnotconverge=NULL;		/* starting pts that did not*/
  static double **badlyconditioned=NULL;  /* points where the jacobian
										  was badly conditioned and 
										  newton raphson was 
										  abandoned */
  static int	num_stat_pts=0;
  static int	num_nonconvergent=0;
  static int	num_badcond=0; 

  static double *grid;					  /* grid */
  static int	grid_size; 

  static double tolerance;				  /* tolerance for NR */
  static int	max_iterations;			  /* for NR */ 

/* function declarations */

void rhs(double *x, double *y);
void Jac(double *x, double *jac);
void sum_vectors(double *y, double *x, double a, int size);
double twonorm(double *y, double *x, int size);
int newtonraphson(double *result, double *x0, double tol, int max_iter);
void traverse_grid(int index);
void DoEwEv(int index, double *point, int ndigits);

/* for Lapack */

void dgeco_(double *a, int *lda, int *n, int *ipvt, double *rcond, double *z);
void dgesl_(double *a, int *lda, int *n, int *ipvt, double *b,  int *job);
void dgeev_(char *jobvl, char *jobvr, int *n, double *jacobian, 
			int *lda, double *wr, double *wi, double *vl, 
			int *ldvl, double *vr, int *ldvr, double *work, 
			int *lwork, int *info);	



/*** Main program **********************************************************/
 
int main(int argc, char **argv)
{
  char            *infile;                   /* pointer to input file name */
  FILE            *fp;                       /* pointer to input data file */
 
  int             c;                 /* used to parse command line options */
  int             num_nucs;                 /* number of nuclei */
  int             i,ii,jj,kk;                      /* loop counter */
  int 			  length;			/* to check genotype's validity */	
  double		  grid_increment;

/* the follwoing few variables are read as arguments to command line opts */

  int             ndigits  = 14;        /* precision of output, default = 6 */

  double		  extinput_time = 67.975;		/* time for which the
  external inputs value should be taken */

  char            *section_title;                /* parameter section name */
  int			  gsize = 5;	/* grid size - will be made global in
														  grid_size */
  double		  grid_start=0.0;	
  double		  grid_end=250.0;

  double 		  tol=1e-6;				  /* tolerance for NR */
  int 	 		  max_iter=100;			  /* for NR */ 

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




/* external declarations for command line option parsing (unistd.h) */

  extern char *optarg;                     /* command line option argument */
  extern int  optind;                /* pointer to current element of argv */
  extern int  optopt;               /* contain option character upon error */

/* following part sets default values for deriv, Jacobian and solver funcs */
/* and section title                                                       */


  section_title = (char *)calloc(MAX_RECORD, sizeof(char));
  section_title = strcpy(section_title, "eqparms");  /* default is eqparms */

/* following part parses command line for options and their arguments      */
/* modified from original getopt manpage                                   */
 
  optarg = NULL;
  while ( (c = getopt(argc, argv, OPTS)) != -1 )
    switch (c) {
    case 'f':
      ndigits = atoi(optarg);    
      if ( ndigits < 0 )
	error("newtonraphson: what exactly would a negative precision be???");
      if ( ndigits > MAX_PRECISION )
	error("newtonraphson: max. float precision is %d!", MAX_PRECISION);
      break;
    case 'g':
      grid_start = atof(optarg);  
      break;
    case 'G':
      grid_end = atof(optarg);     
      break;
    case 'h':                                            /* -h help option */
		PrintMsg(help, 0);
	  	break;
    case 'I':
      max_iter = atoi(optarg);    
      if ( max_iter < 1 )
	error("newtonraphson: Atleast one iteration required");
	  break;	
    case 's':
   		gsize = atoi(optarg);       
		if ( gsize < 2 )
			error("newtonraphson: Need atleast two grid points!");
     	break;
    case 't':
      extinput_time = atof(optarg);  
      if ( extinput_time < 0 )
	error("newtonraphson: the external input time (%g) doesn't make sense", extinput_time);
     break;
    case 'T':
      tol = atof(optarg);             
      if ( tol < EPSILON )
	error("newtonraphson: the tolerance (%g) is too small", tol);
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


  if ( optind != argc - 1 )
     PrintMsg(usage, 1);

/* let's get started and open the data file here */

  infile = argv[optind];
  fp = fopen(infile, "r");      
  if( !fp )
    file_error("newtonraphson");
 
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

  if (grid_start + EPSILON > grid_end)
  	error("newtonraphson: grid start %f and grid end %f too close",
												grid_start, grid_end);

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

  /* now lets do newton-raphson */

  /* assign params to static variables so rhs(), Jac() can see them.
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
  jacobian = (double *) calloc(ngenes*ngenes, sizeof(double));

  grid_size = gsize;
  grid = (double *) calloc(grid_size, sizeof(double));
  grid_increment = (grid_end - grid_start)/((double) (grid_size - 1));

  tolerance = tol;
  max_iterations = max_iter;

  for (ii=0; ii<grid_size; ii++) {
  	grid[ii] = (double) grid_start+ii*grid_increment;
  }

  getrusage (RUSAGE_SELF, &begin);  /*      get start time */
  
  traverse_grid(0);

  getrusage (RUSAGE_SELF, &end);     /* get end time */

/* print the settings default and from commandline */

  printf(settings,infile,extinput_time, gsize, ngenes, grid_start, 
  								grid_end, ndigits, tol, max_iter);

  printf("$cpu_time\n%.14f\n$$\n\n",tvsub(end,begin));

  for (ii=0; ii<num_stat_pts; ii++) {
  	printf("$zero.%d\n",ii);
	for (jj=0; jj<defs.ngenes; jj++)	
		printf(" %*.*f",ndigits+5,ndigits,stationary_points[ii][jj]); 
  	printf("\n$$\n\n",ii);
	DoEwEv(ii,stationary_points[ii],ndigits);
  }	

  printf("$did_not_converge\n"); 
  for (ii=0; ii<num_nonconvergent; ii++) {
	for (jj=0; jj<defs.ngenes; jj++)	
		printf(" %*.*f",ndigits+5,ndigits,didnotconverge[ii][jj]); 
  	printf("\n");
  }	
  printf("$$\n\n");

  printf("$badly_conditioned\n"); 
  for (ii=0; ii<num_badcond; ii++) {
	for (jj=0; jj<defs.ngenes; jj++)	
		printf(" %*.*f",ndigits+5,ndigits,badlyconditioned[ii][jj]); 
  	printf("\n");
  }	
  printf("$$\n\n");
 

  
/* end of doing newton-raphson */

  FreeInterpObject(extinp_polation);
  free(extinp_polation);
  free(interp_dat);
  free(v_ext);
  free(bot);
  free(bot2);
  free(vinput);
  free(xin);
  free(xout);
  free(jacobian);	
  free(grid);

/* the following three were allocates in traverse_grid() */

  if (num_stat_pts)	
	  for (ii=0; ii<num_stat_pts; ii++)
		free(stationary_points[ii]);
	  free(stationary_points);	

  if (num_nonconvergent)	
	  for (ii=0; ii<num_nonconvergent; ii++)
		free(didnotconverge[ii]);
	  free(didnotconverge);	

  if (num_badcond)	
	  for (ii=0; ii<num_badcond; ii++)
		free(badlyconditioned[ii]);
	  free(badlyconditioned);	

  free(section_title);
  free(genotype);


  return 0;
}

/* Function to calculate and print eigenvalues and eigenvectors */

void DoEwEv(int index, double *point, int ndigits)
{

/* eigenvalue and eigenvector arrays */

 double			*ewr,*ewi,*ev;		/* eigenvalues (ewr+i*ewi) and
													 eigenvectors ev */

 /* for dgeev() the lapack eigenvalue/vector driver routine */

 double			*dummy_ev;			/*dummy array for the right
									 eigenvectors that dgeev does 
									 not reference */
 char			jobvl='V';			/* to tell it to calculate the
									 left eigenvectors */
 char			jobvr='N';			/* tell it not to calculates the right
									 eigenvectors */
 double			*work_vector;		/* space for it to do its thing */
 int			lwork;				/* size of the work array */
 int			exit_status;		/* exit status */

 int 			ii,jj,kk;			/* counters */

/* for calculating eigenvalues and eigenvectors */

  lwork = 10*ngenes;	
  ewr = (double *) calloc(ngenes, sizeof(double));
  ewi = (double *) calloc(ngenes, sizeof(double));
  ev = (double *) calloc(ngenes*ngenes, sizeof(double));
  work_vector = (double *) calloc(lwork, sizeof(double));

  Jac(point, jacobian);

 
 /* below we calculate the eigenvalues and eigenvectors. With dgeev_
  * you have the option of calculating the "left" or "right"
  * eigenvectors, or both (see dgeev.f in the linpack directory). It
  * seems the "right" eigenvectors are what we want AX=LX. But it
  * actually turns out that it is the "left" ones that satisfy that
  * relationship. So we are calculating those. Also the matrix needs
  * to be transposed */

  dgeev_(&jobvl, &jobvr, &ngenes, jacobian, &ngenes, ewr, ewi,
  			ev, &ngenes, dummy_ev, &ngenes, work_vector, &lwork,
													&exit_status);	

  if (exit_status == 0) {
	  printf("$eigenvalues.%d\n",index); 
	  for (jj=0; jj<ngenes; jj++) {
		if (ewi[jj] == 0.0) 
			printf(" %*.*f +%1.1fi",ndigits+5,ndigits,
									ewr[jj],ewi[jj]); 
		else if (ewi[jj] > 0.0) 
			printf(" %*.*f +%.*fi",ndigits+5,ndigits,
									ewr[jj],ndigits,ewi[jj]); 
		else 
			printf(" %*.*f %.*fi",ndigits+5,ndigits,
									ewr[jj],ndigits,ewi[jj]); 
	  }	
	  printf("\n$$\n\n");
	  printf("$eigenvectors.%d\n",index); 
	  for (jj=0; jj<ngenes; jj++) {
		kk = 0;
		while (kk<ngenes) {	
			if (fabs(ewr[kk] - ewr[kk+1]) > DBL_EPSILON) {
				printf(" %*.*f +%1.1fi",ndigits+5,ndigits,
									ev[(kk*ngenes)+jj],0.0); 
				kk+=1;
			} else {
				if (ev[((kk+1)*ngenes)+jj] > 0.0) {
					printf(" %*.*f +%.*fi",ndigits+5,ndigits,
							ev[(kk*ngenes)+jj],ndigits,
							fabs(ev[((kk+1)*ngenes)+jj])); 
					printf(" %*.*f -%.*fi",ndigits+5,ndigits,
							ev[(kk*ngenes)+jj],ndigits,
							fabs(ev[((kk+1)*ngenes)+jj])); 
				} else {
					printf(" %*.*f -%.*fi",ndigits+5,ndigits,
							ev[(kk*ngenes)+jj],ndigits,
							fabs(ev[((kk+1)*ngenes)+jj])); 
					printf(" %*.*f +%.*fi",ndigits+5,ndigits,
							ev[(kk*ngenes)+jj],ndigits,
							fabs(ev[((kk+1)*ngenes)+jj])); 
				}
				kk+=2;
			}	
		}	
		printf("\n");
	  }	
	  printf("$$\n\n");
  } else if (exit_status > 0) {
		printf("$comment.%d\nQR did not converge\n$$\n\n",index);
  }	else {
  		printf("$comment.%d\nIllegal argument #%d to dgeev_\n$$\n\n",
											index,-1*exit_status);
  }		



  free(ewr);
  free(ewi);
  free(ev);
  free(work_vector);


}

/* Right hand side of a one-nucleus circuit */

void rhs(double *x, double *y)
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


/*** Jac: Jacobian function for the DvdtOrig model single nucleus	*
 * model; calculates the Jacobian matrix for the equations at a 	*
 * given time t; input concentrations come in x (of size ngenes), 	*
 * the Jacobian is returned in jac; 								*
 *                                                                  *
 * The Equation: dv/dx = Rg(u)       + lambda()                     *
 *                                                                  *
 * The Jacobian: (in this example: nnucs = 3, ngenes = 3            *
 *                                                                  *
 *         b =          1   2   3   								*
 *                                                                  *
 *         a = 1      X11 Y12 Y13 									*
 *         a = 2      Y21 X22 Y23									*
 *         a = 3      Y31 Y32 X33									*
 *                                                                  *
 * Where:  Xab = Rg'(u) - lambda{a} 		                        *
 *         Yab = Rg'(u)                                             *
 *                                                                  *
 ********************************************************************/

void Jac(double *x, double *jac)
{
  int        j;                                  /* local loop counters */
  int        k, kk;               /* index of gene k in a specific nucleus */
  double 	 vinput1 = 0;                    /* used to calculate u */
  double 	 gdot1, vdot1;     /* used to calculate Xab's and Yab's */


/*********************************************************************
 *                                                                   *
 *  g(u)  = 1/2 * ( u / sqrt(1 + u^2) + 1)                           *
 *  g'(u) = 1/2 * ( 1 / ((1 + u^2)^3/2)) * T{ab}                     *
 *                                                                   *
 *********************************************************************/

	for (k=0; k < ngenes; k++) {
	  
	  vinput1  = h[k];
	  vinput1 += m[k] * v_bcd;

	  for(j=0; j < egenes; j++)
		vinput1 += E[(k*egenes)+j] * v_ext[j];

	  for(j=0; j < ngenes; j++)
	    vinput1 += T[(k*ngenes)+j] * x[j];

	  bot2[k] = 1 + vinput1 * vinput1;
	}

/* now calculate sqrt(1+u^2); store it in bot[] */

      for(k=0; k < ngenes; k++)                  /* slow traditional style sqrt */
	bot[k] = sqrt(bot2[k]);

/* resume loop after vector sqrt above; we finish calculating g'(u) and    *
 * place D's in the diagonal of the off-diagonal blocks (cf diagram above) */

	for (k=0; k < ngenes; k++) {
	  
	  gdot1  = 1 / (bot[k] * bot2[k]);
	  gdot1 *= R[k] * 0.5;
	  
	  for (kk=0; kk < ngenes; kk++) {
	    
	    vdot1 = T[(k*ngenes)+kk] * gdot1;

	    if ( k == kk ) 
		  vdot1 -= lambda[k];

	    jac[(k*ngenes)+kk] = vdot1;
	    
	  }
	}
   
  return;

}

/* sum_vectors(): adds x*a to y where x and y are arrays and a is a
 * scalar */ 

void sum_vectors(double *y, double *x, double a, int size)
{

	int ii;								/* counter */

	for (ii=0; ii < size; ii++)
		y[ii] = y[ii] + a*x[ii];

	return;
}	

double twonorm(double *y, double *x, int size)
{

	int ii;								/* counter */
	double norm=0.0;						

	for (ii = 0; ii < size; ii++)
		norm += (y[ii]-x[ii])*(y[ii]-x[ii]);
	
	return sqrt(norm);

}	

int newtonraphson(double *result, double *x0, double tol, int max_iter)
{

	/* for Linpack */

	 double		 	*xi, *xiplus1, *eff;
 	 double			*J;		/* stores the jacobian */
	 int			num_iter = 1;
	 int 			status;
	 int			ii;		/* counter */

	 int			*pivot;			/* pivots for GAussian elimination
	 returned by dgeco_, for use in dgesl_ */
	 double			*work_vector;	/* work vector for dgeco_, usually
	 unimportant */
	 double			cond_number;	/* condition number returned by dgeco_
	 */
	 int			dgesl_flag=1;	/* 0 for solving A*x=b, non-zero for
	 A'*x= b */

	 /* vectors needed for the algorithm */

	 xi = (double *) calloc(ngenes, sizeof(double));
	 xiplus1 = (double *) calloc(ngenes, sizeof(double));
	 eff = (double *) calloc(ngenes, sizeof(double));
	 J = (double *) calloc(ngenes*ngenes, sizeof(double));


	/* allocating the arrays that linpack requires	*/ 
		
	 pivot = (int *) calloc(ngenes, sizeof(int));
	 work_vector = (double *) calloc(ngenes, sizeof(double));
	  
	 memcpy(xiplus1, x0, sizeof(*xiplus1) * ngenes);
	 
	 for (ii=0; ii < ngenes; ii++)
	 	xi[ii] = xiplus1[ii] + 1.0;

	 while ( (twonorm(xi, xiplus1, ngenes) > tol) && (num_iter <=
														 max_iter)) 	
	 {

		 memcpy(xi, xiplus1, sizeof(*xi) * ngenes);
				 
		 rhs(xi, eff);
		 Jac(xi, J);
	 	
		 dgeco_(J, &ngenes, &ngenes, pivot, &cond_number, work_vector);	

		 if (cond_number > MAX_COND)
		 	break;
		 
		 dgesl_(J, &ngenes, &ngenes, pivot, eff, &dgesl_flag);

		 sum_vectors(xiplus1, eff, -1.0, ngenes);

//		 printf("Iteration #%d:%f\n",num_iter,cond_number);
		 num_iter += 1;
		 cond_number = 0;

	 }	 

	 if (num_iter > max_iter) {
		 memcpy(result, x0, sizeof(*x0) * ngenes);
		 status = 1;
	 } else if (cond_number > MAX_COND) {
		 memcpy(result, xiplus1, sizeof(*x0) * ngenes);
		 status = 2;
	 } else {
		 memcpy(result, xiplus1, sizeof(*x0) * ngenes);
		 status = 0;
	 }

	 free(xi);
	 free(xiplus1);
	 free(eff);
	 free(J);

	/* free linpack stuff */
	 free(pivot);
	 free(work_vector);

	 return status;

}	 

/* Recursive function to traverse the R^n grid by generating n-tuples
 * (rather ngenes-tuples) sampled from the grid array. 
 * The first call to it should be with index = 0. It will do a 
 * for loop for each index over the grid points and will call itself
 * with the next higher index. Once it reaches the deepest level
 * (ngenes - 1), it will call newton-raphson with the constructed
 * n-tuple (stored in xin) as the starting point. The returned
 * stationary point is checked in the compiled list of stationary
 * points and is added to the list if unique. It similarly builds
 * lists for non-converging starting points and points that have badly
 * conditioned jacobians */


void traverse_grid(int index)
{

	int ii, jj, kk;	
	int	exit_status;

	if (index < ngenes - 1) {
		for (ii=0; ii<grid_size; ii++) {
			xin[index] = grid[ii];
			traverse_grid(index+1);
		}	
	} else {
		for (ii=0; ii<grid_size; ii++) {
			xin[index] = grid[ii];

/*		  printf("index:%d, ii:%d\n",index,ii);
		  for (kk=0; kk<ngenes; kk++)
			printf("Input %d:%f\t",kk,xin[kk]); 
		  printf("\n");*/

		  exit_status = newtonraphson(xout, xin, tolerance,
													  max_iterations);
		  switch (exit_status) {
			case 0:
				if (!num_stat_pts) {
				
					stationary_points = (double **)
											malloc(sizeof(double *));
					stationary_points[0] = (double *) calloc(ngenes,
														sizeof(double));
														
					 memcpy(stationary_points[0], xout, sizeof(*xout) * ngenes);
					 num_stat_pts = 1;
				} else {

					 jj=0;
					 while ( 
						(jj < num_stat_pts) && 
					 	(twonorm(xout, stationary_points[jj], ngenes) 
									> tolerance)
							)
						jj++;

					 if (jj == num_stat_pts) {
					 	
						num_stat_pts++;
						stationary_points = (double **)
							realloc(stationary_points, sizeof(double *)
														*num_stat_pts);
					 	stationary_points[num_stat_pts-1] = (double *)
									calloc(ngenes, sizeof(double));

						memcpy(stationary_points[num_stat_pts-1], 
										xout, sizeof(*xout) * ngenes);
					 }	
				}
				break;
			case 1:
				if (!num_nonconvergent) {
					
					didnotconverge = (double **) 
											malloc(sizeof(double *));
					didnotconverge[0] = (double *) calloc(ngenes,
													sizeof(double));
													
					 memcpy(didnotconverge[0], xout, sizeof(*xout) * ngenes);
					 num_nonconvergent = 1;

				} else {
					num_nonconvergent++;
					didnotconverge = (double **)
						realloc(didnotconverge, sizeof(double *)
													*num_nonconvergent);
					didnotconverge[num_nonconvergent-1] = (double *)
								calloc(ngenes, sizeof(double));

					memcpy(didnotconverge[num_nonconvergent-1], 
									xout, sizeof(*xout) * ngenes);

				}
				break;
			case 2:
				if (!num_badcond) {
					
					badlyconditioned = (double **) 
											malloc(sizeof(double *));
					badlyconditioned[0] = (double *) calloc(ngenes,
													sizeof(double));
													
					 memcpy(badlyconditioned[0], xout, sizeof(*xout) * ngenes);
					 num_badcond = 1;

				} else {
					num_badcond++;
					badlyconditioned = (double **)
						realloc(badlyconditioned, sizeof(double *)
													*num_badcond);
					badlyconditioned[num_badcond-1] = (double *)
								calloc(ngenes, sizeof(double));

					memcpy(badlyconditioned[num_badcond-1], 
									xout, sizeof(*xout) * ngenes);

				}
				break;
		  	}		
		}	  
	}
}	
