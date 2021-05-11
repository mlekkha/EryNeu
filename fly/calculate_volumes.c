/*******************************************************************
 *                                                                 *
 *   calculate_volumes.c                                           *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 *   Version 0.1	                                               *
 *   written by Manu, adapted from unfold.c						   *	
 *                                                                 *
 *******************************************************************
 * Calculates the volume occupied by solutions that start within   *
 * an IC box. It does so by discretizing the ODEs into a map and   *
 * then treating the repeated application of the map as repeated   *
 * curvilinear coordinate transforms. The using the determinant of *
 * the Jacobian of the map as the scale factor, calculating the    *
 * volume in the transformed coordinates as the integral over the  *
 * IC box.
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
#include <integrate.h>				/*provides BIG_EPSILON among other
											things */
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

"Usage: calculate_volumes  [-f <float_prec>] [-g <grid_start>] \n"
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
  static int    num_nucs;                 /* number of nuclei */
  static double	extinput_time = 67.975;		/* time for which the
							  external inputs value should be taken */
  static double	*v_ext;		/* to hold external input concentrations */
  static double	v_bcd;		/* to hold bcd concentration */
  static double *bot, *bot2, *vinput; 	/* intermediates in rhs() */


  static double	*xi;			/* starting value of newtonraphson */
  static double	*eff;			/* result from newtonraphson 	*/
  static double	*jacobian;		/* jacobian for eigenvalues */
	 
/* for dgeco and dgedi */	 
  static int	*pivot;			/* pivots for GAussian elimination
							 returned by dgeco_, for use in dgesl_ */
  static double	*work_vector;	/* work vector for dgeco_, usually
													 unimportant */

  static double tstep;				  /* time step for map */
  static double timepoints[9] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, \
  								 60.0, 70.0, 80.0};
								 /* time points at which volumes
								  * desired. Fix to be read from file
								  * in future */

/* function declarations */

void rhs(double *x, double *y);
void rhsmit(double *x, double *y);
void Jac(double *x, double *jac);
void sum_vectors(double *y, double *x, double a, int size);
double twonorm(double *y, double *x, int size);
double trapzd(double (*func)(double *, int), double *point, int dim, \
					double *left, double *right, int index, int n);
double tfunc(double *point, int dim);
void trapzdv(void (*func)(double *, int, double *, int), \
				double *answer, int ansdim, double *point, int dim, \
				double *left, double *right, int index, int n);
				
void tfuncv(double *answer, int ansdim, double *point, int dim);
void volelement(double *answer, int ansdim, double *point, int dim);

/* for Lapack */

void dgeco_(double *a, int *lda, int *n, int *ipvt, double *rcond, double *z);
void dgedi_(double *a, int *lda, int *n, int *ipvt, double *det, \
												double *z, int *job);



/*** Main program **********************************************************/
 
int main(int argc, char **argv)
{
  char            *infile;                   /* pointer to input file name */
  FILE            *fp;                       /* pointer to input data file */
 
  int             c;                 /* used to parse command line options */
  int             i,ii,jj,kk;                      /* loop counter */
  int 			  length;			/* to check genotype's validity */	

/* the follwoing few variables are read as arguments to command line opts */

  int             ndigits  = 14;        /* precision of output, default = 6 */


  char            *section_title;                /* parameter section name */
  double 		  timestep=0.1;				  /* time step for map */

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

 double			  points[4] = {0.0,0.0,0.0,0.0};
 double			  leftlim[4] = {0.0,0.0,0.0,0.0};
 double			  rightlim[4] = {20.0,80.0,80.0,80.0};
 double			  answer[2] = {0.0, 0.0};
 double			  volumes[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
 								0.0, 0.0, 0.0};





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
	error("calculate_volumes: what exactly would a negative precision be???");
      if ( ndigits > MAX_PRECISION )
	error("calculate_volumes: max. float precision is %d!", MAX_PRECISION);
      break;
    case 'h':                                            /* -h help option */
		PrintMsg(help, 0);
	  	break;
    case 't':
      extinput_time = atof(optarg);  
      if ( extinput_time < 0 )
	error("calculate_volumes: the external input time (%g) doesn't make sense", extinput_time);
     break;
    case 'T':
      timestep = atof(optarg);             
      if ( timestep < BIG_EPSILON )
	error("calculate_volumes: the time step (%g) is too small",\
															timestep);
     break;
     case 'v':                         /* -v prints version message */
      fprintf(stderr, verstring, *argv, 
	      VERS, USR, MACHINE, COMPILER, FLAGS, __DATE__, __TIME__);
      exit(0);
   case 'x':
      if ( (strcmp(optarg, "input")) && (strcmp(optarg, "eqparms")) &&
	   (strcmp(optarg, "parameters")) )
	error("calculate_volumes: invalid section title (%s)", optarg);
      section_title = strcpy(section_title, optarg);
      break;
    case ':':
      error("calculate_volumes: need an argument for option -%c", optopt);
      PrintMsg(usage, 1);
	  break;
    case '?':
    default:
      error("calculate_volumes: unrecognized option -%c", optopt);
      PrintMsg(usage, 1);
	  break;
    }


  if ( optind != argc - 1 )
     PrintMsg(usage, 1);

/* let's get started and open the data file here */

  infile = argv[optind];
  fp = fopen(infile, "r");      
  if( !fp )
    file_error("calculate_volumes");
 
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
      error("calculate_volumes: genotype length doesn't match number of genes (%d)", 
	    defs.ngenes);
    for (i=0; i<defs.ngenes; i++) {
      c = (int)*(genotype+i);
      if ( ( c != 87 ) && (( c < 82 ) || ( c > 84 )) )
	error("calculate_volumes: genotype string can only contain R, S, T or W");
    }
  }

/* Now we will read and traverse the genotypes list and find the
external inputs title and set the interpolation structure in
integrate.c. The genindex is set to its default value 0, so we are
only doing the first genotype here */

	geno = ReadGenotypes(fp);
	curr = geno;

	if ( !(genindex < count_Slist(geno)))
		error("calculate_volumes: genindex more than last allele!\n");

	i = 0;
	while (i++ < genindex) curr = curr -> next;

/* Below we prepare the interpolant for external inputs */		

  	if ( !(interp_dat=(DataTable **)calloc(1, 
  											sizeof(DataTable *))) )
    	error("calculate_volumes: could not allocate interp_dat struct");
	

  	if ( !(extinp_polation=(InterpObject *)calloc(1, sizeof(InterpObject))) )
    	error("calculate_volumes: could not allocate extinp_polation struct");
	
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

  xi = (double *) calloc(ngenes, sizeof(double));
  eff = (double *) calloc(ngenes, sizeof(double));
  jacobian = (double *) calloc(ngenes*ngenes, sizeof(double));
  pivot = (int *) calloc(ngenes, sizeof(int));
  work_vector = (double *) calloc(ngenes, sizeof(double));
	  

  tstep = timestep;

  getrusage (RUSAGE_SELF, &begin);  /*      get start time */
  
  volelement(volumes, 9, points, ngenes);

  trapzdv(volelement, volumes, 9, points, 4, leftlim, rightlim, 0, 81);
  for (ii=0; ii < 9; ii++)
  	printf("%f\t%f\n",timepoints[ii],volumes[ii]);

  getrusage (RUSAGE_SELF, &end);     /* get end time */

/* print the settings default and from commandline */

/*  printf(settings,infile,extinput_time, gsize, ngenes, grid_start, 
  								grid_end, ndigits, timestep, );*/

  printf("$cpu_time\n%.14f\n$$\n\n",tvsub(end,begin));

  
/* end of doing newton-raphson */

  FreeInterpObject(extinp_polation);
  free(extinp_polation);
  free(interp_dat);
  free(v_ext);
  free(bot);
  free(bot2);
  free(vinput);
  free(xi);
  free(eff);
  free(jacobian);	
  free(pivot);
  free(work_vector);

  free(section_title);
  free(genotype);


  return 0;
}

/* evaluate the kernel in the volume integral at a particular initial
 * condition */

void volelement(double *answer, int ansdim, double *point, int dim)
{
	
	int ii,jj;					/* counters */
	int ansindex = 0;			/* to parse the asnwers array */
	double detjac = 1.0;		/* determinant of the jacobian at time
								   t */
	double detprod = 1.0;		/* products of the determinant until
							       time t-tstep	*/ 
	double detarr[2];			/*determinanyt array for dgedi */							   

	double t;					/* time */
    double cond_number;			/* condition number returned by dgeco_
	 */
	int job=10;					/* tell dgedi to make determinant */ 

	answer[0] = detprod;
	ansindex++;

	for (ii=0; ii<dim; ii++)
		xi[ii] = point[ii];

	for (t=0.0; t <= timepoints[ansdim - 1]; t += tstep)
	{

		detprod *= detjac;
		
/*		if (fabs(t-floor(t)) < BIG_EPSILON)
		{

			printf(" 8226\t%.3f\t",t);
			for (jj=0; jj< dim; jj++)
				printf("%.6f\t", xi[jj]);
			printf("%.6f\n",detprod);
		}*/

		if (fabs(t - timepoints[ansindex]) < BIG_EPSILON)
		{
			answer[ansindex] = detprod;
			ansindex++;

		}

		ExternalInputs(t, extinput_time, v_ext, num_nucs*egenes);

		/* implement mitosis as a kludge */
		
		if ((t >= 16.0) && (t < 21.1))
		{
			rhsmit(xi, eff);
			
			detjac = 1.0;

			for (jj=0;jj<ngenes;jj++)
				detjac *= (1.0 - lambda[jj]*tstep);
			
		} else {	
			rhs(xi, eff);
			Jac(xi, jacobian);

/*			if (fabs(t-floor(t)) < BIG_EPSILON)
			{
				for (jj=0; jj< dim; jj++)
					printf("%.6f\t", xi[jj]);
				printf("\n");	
				for (ii=0; ii < ngenes; ii++) {
					for (jj=0; jj < ngenes; jj++)
						printf(" %.8f", jacobian[ii*dim+jj]);
					printf("\n");		
				}
			printf("\n");		
			}*/

			for (ii=0; ii < ngenes; ii++)
				for (jj=0; jj < ngenes; jj++)
					if (ii == jj)
					{
					
						jacobian[ii*dim+jj] *= tstep;
						jacobian[ii*dim+jj] += 1.0;
					
					} else {
					
						jacobian[ii*dim+jj] *= tstep;
						
					}


			dgeco_(jacobian, &ngenes, &ngenes, pivot, \
										&cond_number, work_vector);	

    	    if (cond_number > MAX_COND)
				error("calculate_volumes: Badly conditioned jacobian");
		 
			dgedi_(jacobian, &ngenes, &ngenes, pivot, \
										detarr, work_vector, &job);	

			detjac = detarr[0]*exp(log(10.0)*detarr[1]);

		}

		sum_vectors(xi, eff, tstep, dim);

	}	
	
	return;	

}


/* trapezoidal rule adapted from Numerical Recipes, 2nd Edition, Page
   137. Returns vector as answer (do multiple functions together.*/

void trapzdv(void (*func)(double *, int, double *, int), \
				double *answer, int ansdim, double *point, int dim, \
				double *left, double *right, int index, int n)
{

	double x,del;
	double *sum, *lans, *temps;	/* To hold the partial sums (sum), and 
							   array to pass to func/trapzdv for
							   answer (lans). */	
	int j;
	double a,b;

	a = left[index];
	b = right[index];
	del = (b-a)/(n-1);

	sum = (double *) calloc(ansdim, sizeof(double));
	lans = (double *) calloc(ansdim, sizeof(double));
	temps = (double *) calloc(ansdim, sizeof(double));

	if (dim - index <= 1) {

		point[index] = a;
		(*func)(temps, ansdim, point, dim);
		
		sum_vectors(sum, temps, 0.5*del, ansdim);

		point[index] = b;
		(*func)(temps, ansdim, point, dim);

		sum_vectors(sum, temps, 0.5*del, ansdim);

	} else {

		point[index] = a;
		trapzdv(func, temps, ansdim, point, dim, left, right, \
														index+1, n);
		
		sum_vectors(sum, temps, 0.5*del, ansdim);
		
		point[index] = b;
		trapzdv(func, temps, ansdim, point, dim, left, right, \
														index+1, n);
		
		sum_vectors(sum, temps, 0.5*del, ansdim);
	}

	x = a + del;

	for (j=1;j<n-1;j++,x+=del)
	{

		if (dim - index <= 1) {

			point[index] = x;
			(*func)(temps, ansdim, point, dim);
			
			sum_vectors(sum, temps, del, ansdim);
		
		} else {

			point[index] = x;
			trapzdv(func, temps, ansdim, point, dim, left, right, \
														index+1, n);

			sum_vectors(sum, temps, del, ansdim);
		}

	}
	
	for (j=0;j<ansdim;j++)
		answer[j] = sum[j];
	
	free(sum);
	free(lans);
	free(temps);

	return;

}	

void tfuncv(double *answer, int ansdim, double *point, int dim)
{

	double *accum;
	int j;

	accum = (double *) calloc(ansdim, sizeof(double));

	for (j=0;j<dim;j++)
		accum[0] = accum[0]+point[j];

	accum[1] = sin(3.14159265358*point[3])*(point[0]*point[1]*point[2]);

	for (j=0;j<ansdim;j++)
		answer[j] = accum[j];

	free(accum);
	return;	

}


/* trapezoidal rule adapted from Numerical Recipes, 2nd Edition, Page
   137. 	*/

double trapzd(double (*func)(double *, int), double *point, int dim, \
					double *left, double *right, int index, int n)
{

	double x,sum,del;
	int j;
	double a,b, tempa, tempb;

	a = left[index];
	b = right[index];
	sum = 0.0;	
	del = (b-a)/(n-1);

	if (dim - index <= 1) {

		point[index] = a;
		tempa = (*func)(point, dim);

		point[index] = b;
		tempb = (*func)(point, dim);

		sum = 0.5*del*(tempa+tempb);

	} else {

		point[index] = a;
		tempa = trapzd(func, point, dim, left, right, index+1, n);
		
		point[index] = b;
		tempb = trapzd(func, point, dim, left, right, index+1, n);
		
		sum = 0.5*del*(tempa+tempb);
	}

	x = a + del;

	for (j=1;j<n-1;j++,x+=del)
	{

		if (dim - index <= 1) {

			point[index] = x;
			sum += del*(*func)(point, dim);
		
		} else {

			point[index] = x;
			sum += del*trapzd(func, point, dim, left, right, \
												index+1, n);
		}

	}
	
	return sum;

}	

double tfunc(double *point, int dim)
{

	double accum=0.0;
	int j;

/*	for (j=0;j<dim;j++)
		accum = accum+point[j];*/

	accum = cos(3.14159265358*point[3])*(point[0]*point[1]*point[2]);

	return accum;	

}

/* Right hand side of a one-nucleus circuit during mitosis*/

void rhsmit(double *x, double *y)
{

  int        k;                   /* index of gene k */
  double 	 vdot1;

	  for( k=0; k<ngenes; k++ ) { 
	  
	    vdot1   = -lambda[k] * x[k];
	    y[k] = vdot1; 
	  }
	
	return;
}


/* Right hand side of a one-nucleus circuit during interphase*/

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
