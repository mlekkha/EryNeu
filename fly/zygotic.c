/*****************************************************************
 *                                                               *
 *   zygotic.c                                                   *
 *                                                               *
 *****************************************************************
 *                                                               *
 *   written by JR, modified by Yoginho                          *
 *   additional g(u)'s by Yousong Wang, Feb 2002                 *
 *   guts code by Manu, June/July 2002                           *
 *                                                               *
 *****************************************************************
 *                                                               *
 * This file is for functions that deal with the right hand side *
 * of the ODEs. We have an initializing function called          *
 * InitZygote that reads parameters, defs data and search space  *
 * limits from a data file and installs a couple of things like  *
 * the solver and temporary arrays for the dvdt function.        *
 * Then we have the DvdtOrig's, which represent the inner loop   *
 * of the model. Each of these derivative functions calculates   *
 * the derivatives of a particular right-hand-side at a given    *
 * time point. They are called by the solver and need to know if *
 * we're in mitosis or interphase, since the equations differ    *
 * in each case.                                                 *
 * Last but not least, we have a few mutator functions in this   *
 * file. They change the local lparm EqParm struct.              *
 *                                                               *
 *****************************************************************   
 *                                                               *
 * Copyright (C) 1989-2003 John Reinitz                          *
 * the full GPL copyright notice can be found in lsa.c           *
 *                                                               *
 *****************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>   

#include <error.h>                                   /* for error handling */
#include <maternal.h>  /* need this for BArrPtr (for the Bicoid structure) */
#include <solvers.h>                                        /* for p_deriv */
#include <integrate.h>                                        /* for
ExternalInputs */
#include <zygotic.h>                                          /* obviously */



/* following is a superfast vector square root function available for DEC  */

#ifdef ALPHA_DU
extern void vsqrt_();
extern void vexp_();
#endif




/* STATIC VARIABLES PRIVATE TO THIS FILE ***********************************
 * most of these static variables are used to speed up the inner loop (the *
 * dvdt function of the model for annealing purposes                       *
 ***************************************************************************/

/* following two structs are used for equation paramters; parm holds       */
/* parameters BEFORE mutation and is returned by GetParameters (i.e. to    */
/* Score(), which would produce an error on mutated parameters); lparm is  */
/* a copy of the struct that gets mutated (either Rs or Ts get set to zero */
/* and therefore violate limits if sent to Score()) and is then used by    */
/* the derivative function DvdtOrig                                        */

static EqParms parm;              /* static struct for equation parameters */
static EqParms lparm;         /* local copy of parameter that gets mutated */
                                                             
/* following arrays are used by DvdtOrig */

static double  *D;                 /* contains info about diffusion sched. */
static double  *vinput;               /* vinput, bot2 and bot are used for */
static double  *bot2, *bot;       /* storing intermediate stuff for vector */
                                    /* functions used by the dvdt function */

/* two local copies for propagation rule and genotype number */

static int     rule;                      /* propagation rule for DvdtOrig */
static int 	   genindex;        /* genotype index, needed by DvdtOrig */
                                /* for getting bicoid from maternal.c */



/*** INITIALIZATION FUNCTIONS **********************************************/

/*** InitZygote: makes pm and pd visible to all functions in zygotic.c and *
 *               reads EqParms and TheProblem; it then initializes bicoid  *
 *               and bias (including BTimes) in maternal.c; lastly, it     *
 *               allocates memory for structures used by the derivative    *
 *               functions                                                 *
 ***************************************************************************/

void InitZygote(FILE *fp, 
		void (*pd)(double *, double, double *, int),
		void (*pj)(double, double *, double *, double **, int),
		char *parm_section)
{



/* install the dvdt function: makes pd global (solvers need to access it) */

  p_deriv   = pd;
  p_jacobn  = pj;
  d_deriv   = DvdtDelay;

/* read equation parameters and the problem */
 
  defs = ReadTheProblem(fp);                
  parm = ReadParameters(fp, parm_section);  


/* install bicoid and bias and nnucs in maternal.c */

  InitBicoid(fp);
  InitBias(fp);  
  InitNNucs(); 


                          
/*** The following things are static to zygotic.c: D contains diffusion ****
 *   parameters converted according to the diffusion schedule; vinput,     *
 *   bot2 and bot are arrays for intermediate results used by DvdtOrig,    *
 *   which are collected for fast DEC vector functions (see DvdtOrig for   *
 *   details)                                                              *
 ***************************************************************************/

  D      = (double *)calloc(defs.ngenes, sizeof(double));
  vinput = (double *)calloc(defs.ngenes * defs.nnucs, sizeof(double));
  bot2   = (double *)calloc(defs.ngenes * defs.nnucs, sizeof(double));
  bot    = (double *)calloc(defs.ngenes * defs.nnucs, sizeof(double));
}



/*** InitGenotype: function used to make genotype static to zygotic.c ******
 *                 since derivative functions later need to read it for    *
 *                 getting the right bicoid gradients from maternal.c      *
 ***************************************************************************/

void InitGenotype(int g_index)
{
	genindex = g_index;
}





/*** CLEANUP FUNCTIONS *****************************************************/

/*** FreeZygote: frees memory for D, vinput, bot2 and bot arrays ***********
 ***************************************************************************/

void FreeZygote(void)
{
  free(D);
  free(vinput);
  free(bot2);
  free(bot);
}



/*** FreeMutant: frees mutated parameter struct ****************************
 ***************************************************************************/

void FreeMutant(void)
{
  free(lparm.R);
  free(lparm.T);
  if (defs.egenes > 0)
	  free(lparm.E);
  free(lparm.m);
  free(lparm.h);
  free(lparm.d);
  free(lparm.lambda);
  free(lparm.tau);
}






/*** DERIVATIVE FUNCTIONS **************************************************/

/*** Derivative functions: *************************************************
 *                                                                         *
 *   These functions calculate derivatives for the solver function based   *
 *   on the state variables in array v for the time t; the output is       *
 *   written into array vdot; n describes the size of both v and vdot      *
 *                                                                         *
 *   There are different g(u) functions that can be used in each deriva-   *
 *   tive function (see comments in the code below). They are set by the   *
 *   command line option -g.                                               * 
 *                                                                         *
 *   There are two rules for the derivative function, one for INTERPHASE   *
 *   and one for MITOSIS. The difference between the two is that only de-  *
 *   cay and diffusion happen during mitosis and the regulation term is    *
 *   only included during INTERPHASE. The rule gets set in Blastoderm      *
 *   (using SetRule) to ensure that the derivative function is in the      *
 *   right rule.                                                           *
 *                                                                         *
 *   JR: We get rid of 3 of 4 divisions in the inner loop by nesting. We   *
 *   get rid of an if by adding one iteration of the loop before and one   *
 *   after for the diffusion (wall) boundary conditions.                   *
 *                                                                         *
 ***************************************************************************/


/*** DvdtOrig: the original derivative function; implements the equations **
 *             as published in Reinitz & Sharp (1995), Mech Dev 49, 133-58 *
 *             plus different g(u) functions as used by Yousong Wang in    *
 *             spring 2002.                                                *
 ***************************************************************************/

void DvdtOrig(double *v, double t, double *vdot, int n)
{
    double vinput1 = 0;     
  int        m;                                        /* number of nuclei */
  int        ap;                 /* nuclear index on AP axis [0,1,...,m-1] */
  int        i,j;                                   /* local loop counters */
  int        k;                   /* index of gene k in a specific nucleus */
  int        base, base1;            /* index of first gene in a specific nucleus */
  int 		 *l_rule;		/* for autonomous implementation */
#ifdef ALPHA_DU
  int        incx = 1;        /* increment step size for vsqrt input array */
  int        incy = 1;       /* increment step size for vsqrt output array */
#endif

  static int num_nucs  = 0;      /* store the number of nucs for next step */
  static int bcd_index = 0;   /* the *next* array in bicoid struct for bcd */
  static DArrPtr bcd;              /* pointer to appropriate bicoid struct */
  double *v_ext;			/* array to hold the external input
									  concentrations at time t */



/* get D parameters and bicoid gradient according to cleavage cycle */

  m = n / defs.ngenes;                        /* m is the number of nuclei */
  if (m != num_nucs) {      /* time-varying quantities only vary by ccycle */
    GetD(t, lparm.d, defs.diff_schedule, D);  
                      /* get diff coefficients, according to diff schedule */
    if (num_nucs > m)                  /* started a new iteration in score */
      bcd_index = 0;                           /* -> start all over again! */
    num_nucs = m;                         /* store # of nucs for next step */
    bcd = GetBicoid(t, genindex);                   /* get bicoid gradient */
    if( bcd.size != num_nucs) 
     error("DvdtOrig: %d nuclei don't match Bicoid!", num_nucs);
    bcd_index++;                   /* store index for next bicoid gradient */
  }

  l_rule = (int *) calloc(defs.ngenes, sizeof(int));
  for (i = 0; i < defs.ngenes; i++)
 	l_rule[i] = !(Theta(t));	

/* Here we retrieve the external input concentrations into v_ext */
/* 01/13/10 Manu: If there are no external inputs, don't
 * allocate v_ext or populate it with external input
 * concentrations */  

    if (defs.egenes > 0) {
		v_ext = (double *) calloc(m*defs.egenes, sizeof(double));
		ExternalInputs(t, t, v_ext, m*defs.egenes);
    }

/* This is how it works (by JR): 
         
    ap      nucleus position on ap axis
    base    index of first gene (0) in current nucleus
    k       index of gene k in current nucleus 

            Protein synthesis terms are calculated according to g(u)

            First we do loop for vinput contributions; vinput contains
	    the u that goes into g(u)

	    Then we do a separate loop or vector func for sqrt or exp

	    Then we do vdot contributions (R and lambda part)

	    Then we do the spatial part (Ds); we do a special case for 
            each end 

            Note, however, that before the real loop, we have to do a 
            special case for when rhere is one nuc, hence no diffusion

	    These loops look a little funky 'cause we don't want any 
            divides'                                                       */

	/* 01/13/10 Manu: If there are no external inputs, the for
	 * loops for updating vinput1 with external input terms
	 * are never entered */

/***************************************************************************
 *                                                                         *
 *        g(u) = 1/2 * ( u / sqrt(1 + u^2) + 1)                            *
 *                                                                         *
 ***************************************************************************/
    if ( gofu == Sqrt ) {

      for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
	for (i=base; i < base + defs.ngenes; i++) {

	  k = i - base;
	  
	  vinput1  = lparm.h[k];
	  vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

	  for(j=0; j < defs.egenes; j++)
	    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[base1 + j];

	  for(j=0; j < defs.ngenes; j++)
	    vinput1 += lparm.T[(k*defs.ngenes)+j] * v[base + j];

	  bot2[i]   = 1 + vinput1 * vinput1;
	  vinput[i] = vinput1;
	}
      }

                            /* now calculate sqrt(1+u2); store it in bot[] */
#ifdef ALPHA_DU
      vsqrt_(bot2,&incx,bot,&incy,&n);    /* superfast DEC vector function */
#else
      for(i=0; i < n; i++)                  /* slow traditional style sqrt */
	bot[i] = sqrt(bot2[i]);
#endif
            /* next loop does the rest of the equation (R, Ds and lambdas) */
                                                 /* store result in vdot[] */

      if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

	register double vdot1, g1;

	for( base=0; base<n; base+=defs.ngenes ) {
	  for( i=base; i<base+defs.ngenes; i++ ) { 
	  
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      = 1 + vinput[i]/bot[i];         
	    vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;     
	    vdot[i] = vdot1; 
	  }
	}
	
      } else {                    /* then for multiple nuclei -> diffusion */

	register double vdot1,g1;

	for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
	  k       = i;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      = 1 + vinput[i]/bot[i];
	  vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;
	  vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}    
                                                     /* then middle nuclei */
	for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
	  for(i=base; i < base + defs.ngenes; i++){
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      =  1 + vinput[i]/bot[i];
	    vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;
	    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) + 
			      (v[i + defs.ngenes] - v[i]) );
	    vdot[i] = vdot1;
	  }
	}                                    
	                                   /* last: posterior-most nucleus */
	for(i=base; i < base + defs.ngenes; i++){
	  k       = i - base;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      =  1 + vinput[i]/bot[i];
	  vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;
	  vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}
      }

/***************************************************************************
 *                                                                         *
 *        g(u) = 1/2 * (tanh(u) + 1) )                                     *
 *                                                                         *
 ***************************************************************************/

    } else if ( gofu == Tanh ) {

      for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
	for(i=base; i < base + defs.ngenes; i++){
	  
	  k = i - base;
	  
	  vinput1  = lparm.h[k];
	  vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

	  for(j=0; j < defs.egenes; j++)
	    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[base1 + j];

	  for(j=0; j < defs.ngenes; j++)
	    vinput1  += lparm.T[(k*defs.ngenes)+j] * v[base + j];

	  vinput[i] = vinput1;
	}
      }

            /* next loop does the rest of the equation (R, Ds and lambdas) */
                                                 /* store result in vdot[] */

      if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

	register double vdot1, g1;

	for( base=0; base<n; base+=defs.ngenes ) {
	  for( i=base; i<base+defs.ngenes; i++ ) { 
	  
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      = tanh(vinput[i]) + 1;         
	    vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;     
	    vdot[i] = vdot1; 
	  }
	}

      } else {                    /* then for multiple nuclei -> diffusion */

	register double vdot1,g1;
	
	for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
	  k       = i;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      = tanh(vinput[i]) + 1;
	  vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;
	  vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}     
	                                             /* then middle nuclei */
	for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
	  for(i=base; i < base + defs.ngenes; i++){
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      = tanh(vinput[i]) + 1;
	    vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;
	    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) + 
			      (v[i + defs.ngenes] - v[i]) );
	    vdot[i] = vdot1;
	  }
	}                                    
                                           /* last: posterior-most nucleus */
	for(i=base; i < base + defs.ngenes; i++){
	  k       = i - base;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      = tanh(vinput[i]) + 1;
	  vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;
	  vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}
      }

/***************************************************************************
 *                                                                         *
 *        g(u) = 1 / (1 + exp(-2u))                                        *
 *                                                                         *
 ***************************************************************************/

    } else if ( gofu == Exp ) {

      for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
	for(i=base; i < base + defs.ngenes; i++){

	  k = i - base;

	  vinput1  = lparm.h[k];
	  vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

	  for(j=0; j < defs.egenes; j++)
	    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[base1 + j];

	  for(j=0;  j < defs.ngenes; j++)
	    vinput1  += lparm.T[(k*defs.ngenes)+j] * v[base + j];

	  vinput[i] = -2.0 * vinput1;	
	}
      }

                               /* now calculate exp(-u); store it in bot[] */
#ifdef ALPHA_DU
      vexp_(vinput,&incx,bot,&incy,&n);   /* superfast DEC vector function */
#else
      for (i=0; i<n; i++)                    /* slow traditional style exp */
	bot[i] = exp(vinput[i]);
#endif

            /* next loop does the rest of the equation (R, Ds and lambdas) */
                                                 /* store result in vdot[] */

      if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

	register double vdot1, g1;

	for( base=0; base<n; base+=defs.ngenes ) {
	  for( i=base; i<base+defs.ngenes; i++ ) { 
	  
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      = 1 / (1 + bot[i]);         
	    vdot1  += l_rule[k]*lparm.R[k] * g1;     
	    vdot[i] = vdot1; 
	  }
	}

      } else {                    /* then for multiple nuclei -> diffusion */

	register double vdot1,g1;

	for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
	  k       = i;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      = 1 / (1 + bot[i]);
	  vdot1  += l_rule[k]*lparm.R[k] * g1;
	  vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}     
                                                     /* then middle nuclei */
	for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
	  for(i=base; i < base + defs.ngenes; i++){
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      = 1 / (1 + bot[i]);
	    vdot1  += l_rule[k]*lparm.R[k] * g1;
	    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) + 
			      (v[i + defs.ngenes] - v[i]) );
	    vdot[i] = vdot1;
	  }
	}                                    
                                           /* last: posterior-most nucleus */
	for(i=base; i < base + defs.ngenes; i++){
	  k       = i - base;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      = 1 / (1 + bot[i]);
	  vdot1  += l_rule[k]*lparm.R[k]  * g1;
	  vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}
      }

/***************************************************************************
 *                                                                         *
 *        g(u) = 0 if u<0, 1 if u>=0 (Heaviside function)                  *
 *                                                                         *
 *        this makes the model quasi boolean and the equations locally     *
 *        linear for both u>0 and u<0                                      *
 ***************************************************************************/

    } else if ( gofu == Hvs ) {

      for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
	for(i=base; i < base + defs.ngenes; i++){

	  k = i - base;

	  vinput1  = lparm.h[k];
	  vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

	  for(j=0; j < defs.egenes; j++)
	    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[base1 + j];

	  for(j=0; j < defs.ngenes; j++)
	    vinput1  += lparm.T[(k*defs.ngenes)+j] * v[base + j];
	  vinput[i] = vinput1;
	}
      }

     /* next loop does the rest of the equation (R, Ds and lambdas) */
                                                 /* store result in vdot[] */

      if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

	register double vdot1, g1;

	for( base=0; base<n; base+=defs.ngenes ) {
	  for( i=base; i<base+defs.ngenes; i++ ) { 
	    
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    if (vinput[i] >= 0.) 
	      g1 = 1.;
	    else 
	      g1 = 0.;         
	    vdot1  += l_rule[k]*lparm.R[k] * g1;     
	    vdot[i] = vdot1; 
	  }
	}

      } else {                    /* then for multiple nuclei -> diffusion */

	register double vdot1,g1;

	for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
	  k       = i;
	  vdot1   = -lparm.lambda[k] * v[i];
	  if (vinput[i] >= 0.) 
	    g1 = 1.;
	  else 
	    g1 = 0.;
	  vdot1  += l_rule[k]*lparm.R[k] * g1;
	  vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}     
	                                             /* then middle nuclei */
	for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
	  for(i=base; i < base + defs.ngenes; i++){
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    if (vinput[i] >= 0.) 
	      g1 = 1.;
	    else 
	      g1 = 0.;
	    vdot1  += l_rule[k]*lparm.R[k] * g1;
	    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) + 
			      (v[i + defs.ngenes] - v[i]) );
	    vdot[i] = vdot1;
	  }
	}                                    
	                                   /* last: posterior-most nucleus */
	for(i=base; i < base + defs.ngenes; i++){
	  k       = i - base;
	  vdot1   = -lparm.lambda[k] * v[i];
	  if (vinput[i] >= 0.) 
	    g1 = 1.;
	  else 
	    g1 = 0.;
	  vdot1  += l_rule[k]*lparm.R[k] * g1;
	  vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}
      }

    } else
      error("DvdtOrig: unknown g(u)");
    
/* during mitosis only diffusion and decay happen */
    
	free(l_rule);

	if (defs.egenes > 0)
		free(v_ext);
	
  return;
}



/*** JACOBIAN FUNCTION(S) **************************************************/

/*** JacobnOrig: Jacobian function for the DvdtOrig model; calculates the **
 *               Jacobian matrix (matrix of partial derivatives) for the   *
 *               equations at a give time t; input concentrations come in  *
 *               v (of size n), the Jacobian is returned in jac; note that *
 *               all dfdt's are zero in our case since our equations are   *
 *               autonomous (i.e. have no explicit t in them)              *
 ***************************************************************************
 *                                                                         *
 * The Equation: dv/dx = Rg(u) + D() + lambda()                            *
 *                                                                         *
 * The Jacobian: (in this example: nnucs = 3, ngenes = 3                   *
 *                                                                         *
 *                        i = 1       i = 2        i = 3                   *
 *         b =          1   2   3   1   2   3    1   2   3                 *
 *                                                                         *
 *         a = 1      X11 Y12 Y13  D1   0   0    0   0   0                 *
 * i = 1   a = 2      Y21 X22 Y23   0  D2   0    0   0   0                 *
 *         a = 3      Y31 Y32 X33   0   0  D3    0   0   0                 *
 *                                                                         *
 *         a = 1       D1   0   0 X11 Y12 Y13   D1   0   0                 *
 * i = 2   a = 2        0  D2   0 Y21 X22 Y23    0  D2   0                 *
 *         a = 3        0   0  D3 Y31 Y32 X33    0   0  D3                 *
 *                                                                         *
 *         a = 1        0   0   0  D1   0   0  X11 Y12 Y13                 *
 * i = 3   a = 2        0   0   0   0  D2   0  Y21 X22 Y23                 *
 *         a = 3        0   0   0   0   0   0  Y31 Y32 Y33                 *
 *                                                                         *
 * Where:  Xab = Rg'(u) - 2D{a} - lambda{a}                                *
 *         Yab = Rg'(u)                                                    *
 *         Da  = D{a}                                                      *
 *                                                                         *
 ***************************************************************************/

void JacobnOrig(double t, double *v, double *dfdt, double **jac, int n)
{
  int        m;                                        /* number of nuclei */
  int        ap;                 /* nuclear index on AP axis [0,1,...,m-1] */
  int        i, j;                                  /* local loop counters */
  int        k, kk;               /* index of gene k in a specific nucleus */
  int        base;            /* indices of 1st gene in a specific nucleus */
#ifdef ALPHA_DU
  int        incx = 1;        /* increment step size for vsqrt input array */
  int        incy = 1;       /* increment step size for vsqrt output array */
#endif

  static int num_nucs  = 0;      /* store the number of nucs for next step */
  static int bcd_index = 0;   /* the *next* array in bicoid struct for bcd */
  static DArrPtr bcd;              /* pointer to appropriate bicoid struct */


/* get D parameters and bicoid gradient according to cleavage cycle */

  m = n / defs.ngenes;                        /* m is the number of nuclei */
  if (m != num_nucs) {      /* time-varying quantities only vary by ccycle */
    GetD(t, lparm.d, defs.diff_schedule, D);  
                      /* get diff coefficients, according to diff schedule */
    if (num_nucs > m)                  /* started a new iteration in score */
      bcd_index = 0;                           /* -> start all over again! */
    num_nucs = m;                         /* store # of nucs for next step */
    bcd = GetBicoid(t, genindex);                   /* get bicoid gradient */
    if( bcd.size != num_nucs) 
     error("JacobnOrig: %d nuclei don't match Bicoid!", num_nucs);
    bcd_index++;                   /* store index for next bicoid gradient */
  }

/*** INTERPHASE rule *******************************************************/

  if (rule == INTERPHASE) {          

    register double vinput1 = 0;                    /* used to calculate u */
    register double gdot1, vdot1;     /* used to calculate Xab's and Yab's */


/***************************************************************************
 *                                                                         *
 *  g(u)  = 1/2 * ( u / sqrt(1 + u^2) + 1)                                 *
 *  g'(u) = 1/2 * ( 1 / ((1 + u^2)^3/2)) * T{ab}                           *
 *                                                                         *
 ***************************************************************************/

    if ( gofu == Sqrt ) {

/* this looks confusing, but it's actually quite simple: the two loops be- *
 * low are in reality just one that loops over the first dimension of the  *
 * Jacobian (size n) 'nuclear block' by 'nuclear block' (the i's in the    *
 * long introductory comment above); ap keeps track of the nucleus number, * 
 * k keeps track of which gene we're dealing with; bot2 saves 1+u^2        */

      for (base=0, ap=0; base<n ; base+=defs.ngenes, ap++) {
	for (i=base; i < base+defs.ngenes; i++) {

	  k = i - base;
	  
	  vinput1  = lparm.h[k];
	  vinput1 += lparm.m[k] * bcd.array[ap];

	  for(j=0; j < defs.ngenes; j++)
	    vinput1 += lparm.T[(k*defs.ngenes)+j] * v[base + j];

	  bot2[i] = 1 + vinput1 * vinput1;
	}
      }

/* now calculate sqrt(1+u^2); store it in bot[] */

#ifdef ALPHA_DU
      vsqrt_(bot2,&incx,bot,&incy,&n);    /* superfast DEC vector function */
#else
      for(i=0; i < n; i++)                  /* slow traditional style sqrt */
	bot[i] = sqrt(bot2[i]);
#endif

/* resume loop after vector sqrt above; we finish calculating g'(u) and    *
 * place D's in the diagonal of the off-diagonal blocks (cf diagram above) */

      for (base=0; base < n; base+=defs.ngenes) {
	for (i=base; i < base+defs.ngenes; i++) {
	  
	  k = i - base;
	  
	  gdot1  = 1 / (bot[i] * bot2[i]);
	  gdot1 *= lparm.R[k] * 0.5;
	  
	  for (j=base; j < base+defs.ngenes; j++) {
	    
	    kk = j - base;
	    
	    vdot1 = lparm.T[(k*defs.ngenes)+kk] * gdot1;
	    if ( k == kk ) {
	      if ( n > defs.ngenes )
		if ( base > 0 && base < n-defs.ngenes )
		  vdot1 -= 2. * D[k];
		else
		  vdot1 -= D[k];
	      vdot1 -= lparm.lambda[k];
	    }
	    jac[i][j] = vdot1;
	    
	  }
	  
	  if ( base > 0 ) 
	    jac[i][i-defs.ngenes] = D[k];
	  
	  if ( base < n-defs.ngenes ) 
	    jac[i][i+defs.ngenes] = D[k];
	  
	}
      }
      
/***************************************************************************
 *                                                                         *
 * g(u)  = 1/2 * (tanh(u) + 1) )  or  g(u)  = 1 / ( 1 + e^(-2u))           *
 * g'(u) = Sech(u)^2 /2           or  g'(u) = 2e^(-2u) / (1 + e^(-2u))^2   *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 * These are actually the same function in different guises; we implement  *
 * Exp below since it's faster                                             *
 *                                                                         *
 ***************************************************************************/

    } else if ( (gofu == Tanh) || (gofu == Exp) ) {

/* this looks confusing, but it's actually quite simple: the two loops be- *
 * low are in reality just one that loops over the first dimension of the  *
 * Jacobian (size n) 'nuclear block' by 'nuclear block' (the i's in the    *
 * long introductory comment above); ap keeps track of the nucleus number, * 
 * k keeps track of which gene we're dealing with; bot2 saves 1+u^2        */

      for (base=0, ap=0; base<n ; base+=defs.ngenes, ap++) {
	for (i=base; i < base+defs.ngenes; i++) {

	  k = i - base;
	  
	  vinput1  = lparm.h[k];
	  vinput1 += lparm.m[k] * bcd.array[ap];

	  for(j=0; j < defs.ngenes; j++)
	    vinput1 += lparm.T[(k*defs.ngenes)+j] * v[base + j];

	  bot[i] = -2.0 * vinput1;
	}
      }

/* now calculate exp(-u); store it in bot[] */
#ifdef ALPHA_DU
      vexp_(bot,&incx,bot2,&incy,&n);   /* superfast DEC vector function */
#else
      for (i=0; i<n; i++)                    /* slow traditional style exp */
	bot2[i] = exp(bot[i]);
#endif

/* resume loop after vector exp above; we finish calculating g'(u) and     *
 * place D's in the diagonal of the off-diagonal blocks (cf diagram above) */

      for (base=0; base < n; base+=defs.ngenes) {
	for (i=base; i < base+defs.ngenes; i++) {
	  
	  k = i - base;
	  
	  gdot1  = 2. * bot2[i];
	  gdot1 /= (1. + bot2[i]) * (1. + bot2[i]);
	  gdot1 *= lparm.R[k];
	  
	  for (j=base; j < base+defs.ngenes; j++) {
	    
	    kk = j - base;
	    
	    vdot1 = lparm.T[(k*defs.ngenes)+kk] * gdot1;
	    if ( k == kk ) {
	      if ( n > defs.ngenes )
		if ( base > 0 && base < n-defs.ngenes )
		  vdot1 -= 2. * D[k];
		else
		  vdot1 -= D[k];
	      vdot1 -= lparm.lambda[k];
	    }
	    jac[i][j] = vdot1;
	    
	  }
	  
	  if ( base > 0 ) 
	    jac[i][i-defs.ngenes] = D[k];
	  
	  if ( base < n-defs.ngenes ) 
	    jac[i][i+defs.ngenes] = D[k];
	  
	}
      }
      
/*** semi-implicit solvers are NOT allowed with heaviside g(u) *************/

    } else if ( gofu == Hvs ) {
      error("JacobnOrig: can't use semi-implicit solver on heaviside g(u)!");

    } else {
      error("JacobnOrig: unknown g(u) function!\n");
    }
    
/* during mitosis only diffusion and decay happen */
    
  } else if (rule == MITOSIS) {

    register double vdot1;
      
       for (base=0; base < n; base+=defs.ngenes) {
	for (i=base; i < base+defs.ngenes; i++) {
	  k = i - base;
	  
	  for (j=base; j < base+defs.ngenes; j++) {
	    kk = j - base;
	    
	    if ( k == kk ) {
	      vdot1  = -lparm.lambda[k];
	      if ( n < defs.ngenes )
		if ( base > 0 && base < n-defs.ngenes )
		  vdot1 -= 2. * D[k];
		else
		  vdot1 -= D[k];
	    } else 
	      vdot1 = 0.;

	    jac[i][j] = vdot1;
	    
	  }
	  
	  if ( base > 0 ) 
	    jac[i][i-defs.ngenes] = D[k];
	  
	  if ( base < n-defs.ngenes ) 
	    jac[i][i+defs.ngenes] = D[k];
	  
	}
      }
 
  } else 
    error("JacobnOrig: Bad rule %i sent to JacobnOrig", rule);

  return;

}



/*** GUTS FUNCTIONS ********************************************************/

/*** CalcGuts: calculates guts for genotpye 'gtype' using unfold output in *
 *             'table' and the guts string 'gutsdefs'; it returns the num- *
 *             ber of columns we'll need to print and the guts table in    *
 *             'gtable'                                                    *
 ***************************************************************************/

int CalcGuts(int gindex, char *gtype, 
				InterpObject hist_interrp, InterpObject extinp_interrp, 
				NArrPtr table, NArrPtr *gtable, char *gutsdefs)
{
  int           m;			  /* number of nuclei at each time */
  int           i, j, k;			     /* just for the loops */
  int           numguts;               /* number of guts columns requested */
  int           which;	               /* for which gene to calculate guts */
  int           l_rule;                             /* local copy for rule */
  unsigned long *gutcomps = NULL;  /* bit flags for each word of ID string */

  NArrPtr       gutsy;	       /* temporary place for guts, to be returned */

  char          **w     = NULL;       /* w(ord) string array after parsing */
  char          **dummy = NULL;        /* beginning of w(ord) string array */
  char          *targetptr;     /* for parsing the gutsdefs, and a counter */


	
/* allocate memory for parsed gut strings, then parse the string */

  dummy = (char **)calloc(MAX_RECORD, sizeof(char *));
  w = dummy;
  numguts = (ParseString(gutsdefs, w) - 1);   /* w contains array of words */  

/* if string is empty -> return zero */

  if ( numguts == -1 ) {
    while (*++w)
      free(*w);
    free(*dummy);
    free(dummy);
    return numguts+1;
  }

/* if no guts specified but gene name in string -> error */

  if ( numguts == 0 ) 
    error("CalcGuts: no guts specified for target gene %s", w[0]);
	
/* allocate the gutcomps array and fill it with the bit flags correspon-   *
 * ding to the things that need to get calculated below; each word in an   *
 * ID string corresponds to one long int in the gutcomps array; targetptr  *
 * points to the name of the gene in the geneID string for which we cal-   *
 * culate guts                                                             */

  gutcomps = (unsigned long *)calloc(numguts, sizeof(unsigned long));
  if ( !(targetptr = GetGutsComps(defs.gene_ids, defs.egene_ids, w, gutcomps)) ) 
    error("CalcGuts: target gene (%s) is not defined in the geneID string",
	  w[0]);

/* allocate gutsy array and set the right genotype in score.c & zygotic.c */

  gutsy.size  = table.size;	
  gutsy.array = (NucState *)calloc(table.size, sizeof(NucState));
  InitGenotype(gindex);
  InitDelaySolver();
  SetHistoryInterp(hist_interrp);

  if (defs.egenes > 0)
	  SetExternalInputInterp(extinp_interrp);

  SetFactDiscons();
  Mutate(gtype);
 	
/* which tells us which gene we calculate guts for */

  which = targetptr - defs.gene_ids;
	
/* free gut ID strings, we don't need them anymore since we have gutcomps */

  while (*++w)
    free(*w);
  free(*dummy);
  free(dummy);
	
/* then start calculating guts */

  for (i=0; i < table.size; i++) {
    gutsy.array[i].time = table.array[i].time;
	
/* calculate size of 2D array that holds guts and allocate memory */

    m = table.array[i].state.size/defs.ngenes; 
    gutsy.array[i].state.size  = numguts * m;
    gutsy.array[i].state.array = 
      (double *)calloc(gutsy.array[i].state.size, sizeof(double));
	
/* set whether it is MITOSIS or INTERPHASE when calculating guts */   

    l_rule = GetRule(gutsy.array[i].time);
    SetRule(l_rule);

/* call the function that calculates the internals of the RHS */

    CalcRhs(table.array[i].state.array, 
	    gutsy.array[i].time, gutsy.array[i].state.array, 
	    table.array[i].state.size, gutsy.array[i].state.size, 
	    numguts, which, gutcomps);
  }
	
/* clean up */

  FreeMutant();
  FreeFactDiscons();
  FreeDelaySolver();
  free(gutcomps);
		
/* return the guts and number of columns for PrintGuts() */

  *gtable = gutsy;
  return numguts;
}



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
	     int numguts, int which, unsigned long *gutcomps)
{
  int    m;                                            /* number of nuclei */
  int    ap;                     /* nuclear index on AP axis [0,1,...,m-1] */

  int    i, j;                                            /* loop counters */
  int    base, base1, base2;     /* index counters for data, guts, and external inputs respectively */
  unsigned long *tempctr;     /* temporary counter for traversing gutcomps */

  double *deriv;           /* array for the derivatives at a specific time */

  static int num_nucs  = 0;      /* store the number of nucs for next step */
  static int bcd_index = 0;   /* the *next* array in bicoid struct for bcd */

  static     DArrPtr   bcd;        /* pointer to appropriate bicoid struct */
  double *v_ext;			/* array to hold the external input
									  concentrations at time t */
	


/* the following calcululates appropriate D's according to the diffusion   *
 * schedule, gets the bicoid gradient and complains if anything is wrong;  *
 * this all only happens after a cell division or at the beginning         */

  m = n / defs.ngenes;                /* m is the current number of nuclei */

  if (m != num_nucs) {     
    GetD(t, lparm.d, defs.diff_schedule, D);  
    num_nucs = m;         
    bcd = GetBicoid(t, genindex); 
    if ( bcd.size != num_nucs) 
      error("CalcRhs: %d nuclei don't match Bicoid!", num_nucs);
    bcd_index++;    
  }

/* call the derivative function to get the total derivative */

  deriv = (double *)calloc(n, sizeof(double));
  p_deriv(v, t, deriv, n);

/* Here we retrieve the external input concentrations into v_ext */
/* 01/13/10 Manu: If there are no external inputs, don't
 * allocate v_ext or populate it with external input
 * concentrations */  


    if (defs.egenes > 0) {
		v_ext = (double *) calloc(m*defs.egenes, sizeof(double));
		ExternalInputs(t, t, v_ext, m*defs.egenes);
	}	

/* all the code below calculates the requested guts (by checking the bits  *
 * set in the gutcomps array); it does so by forward reconstructing all    *
 * the guts that are part of u and then reverse calculating (starting from *
 * the derivative) the guts outside g(u) such that guts work for any g(u); *
 * most of this is directly derived from the derivative function(s) above  */

/* the next block of code calculates guts which are within u, i.e. all the *
 * regulatory contributions or the threshold of a certain gene             */

	/* 01/13/10 Manu: If there are no external inputs, the for
	 * loops for updating guts with external input terms
	 * are never entered */


  for (base = 0, base1 = 0, base2 = 0, ap=0; 
       base < n, base1 < gn; 
       base += defs.ngenes, base1 += numguts, base2 += defs.egenes, ap++) {

    tempctr = gutcomps;
    for (i=base1; i < base1+numguts; i++) {
      for(j=0; j < defs.ngenes; j++) 
	if ( (*tempctr >> j) % 2) 
	  guts[i] += lparm.T[(which*defs.ngenes)+j]*v[base + j];
								
	  for(j=0; j < defs.egenes; j++)
	if ( (*tempctr >> defs.ngenes+j) % 2) 
	  guts[i] += lparm.E[(which*defs.egenes)+j] * v_ext[base2 + j];

      if ( (*tempctr >> defs.ngenes+defs.egenes) % 2) 
	guts[i]   += lparm.m[which] * bcd.array[ap];
								
      if ( (*tempctr >> defs.ngenes+defs.egenes+1) % 2) 
	guts[i]   += lparm.h[which];
      tempctr++;
    }	
  }

/* the code below 'reverse calculates' everything outside g(u) */

  if ( n == defs.ngenes ) {           /* first part: one nuc, no diffusion */

    for( base=0, base1=0; 
	 base < n, base1 < gn; 
	 base += defs.ngenes, base1 += numguts) {

      tempctr = gutcomps;
      for (i=base1; i < base1+numguts; i++) {
								
	if ( (*tempctr >> defs.ngenes+defs.egenes+2) % 2) 
	  guts[i] += -lparm.lambda[which] * v[base+which];
								
	if ( (*tempctr >> defs.ngenes+defs.egenes+5) % 2) 
	  guts[i] += deriv[base+which];
									
	if ( (*tempctr >> defs.ngenes+defs.egenes+6) % 2) 
	  guts[i] += deriv[base+which] 
                  + lparm.lambda[which]*v[base+which];
								
	if ( (*tempctr >> defs.ngenes+defs.egenes+7) % 2) 
	  guts[i] += (deriv[base+which] 
                  + lparm.lambda[which] * v[base+which]) 
                  / lparm.R[which];
								
	tempctr++;
      }	
    }

  } else {                        /* then for multiple nuclei -> diffusion */

    tempctr = gutcomps;                     /* first anterior-most nucleus */
    for (i=0; i < numguts; i++) {
		       					
      if ( (*tempctr >> defs.ngenes+defs.egenes+2) % 2) 
	guts[i] += -lparm.lambda[which] * v[which];
								
      if ( (*tempctr >> defs.ngenes+defs.egenes+4) % 2) 
	guts[i] += 
	  D[which] * (v[which+defs.ngenes] - v[which]);
									
      if ( (*tempctr >> defs.ngenes+defs.egenes+5) % 2) 
	guts[i] += deriv[which];
									
      if ( (*tempctr >> defs.ngenes+defs.egenes+6) % 2) 
	guts[i] += deriv[which] 
                + lparm.lambda[which] * v[which] 
	        - D[which] * (v[which+defs.ngenes] - v[which]);
								
      if ( (*tempctr >> defs.ngenes+defs.egenes+7) % 2) 
	guts[i] += (deriv[which] 
		+ lparm.lambda[which] * v[which] 
		- D[which] * (v[which+defs.ngenes] - v[which])) 
                / lparm.R[which];
								
      tempctr++;
    }	
        
                                                     /* then middle nuclei */
    for(base = defs.ngenes, base1 = numguts; 
	base < n - defs.ngenes, base1 < gn - numguts; 
	base += defs.ngenes, base1 += numguts) {

      tempctr = gutcomps;
      for (i=base1; i < base1+numguts; i++) {
								
	if ( (*tempctr >> defs.ngenes+defs.egenes+2) % 2) 
	  guts[i] += -lparm.lambda[which] * v[base+which];
								
	if ( (*tempctr >> defs.ngenes+defs.egenes+3) % 2) 
	  guts[i] += D[which] * (v[base+which-defs.ngenes] - v[base+which]); 
								
	if ( (*tempctr >> defs.ngenes+defs.egenes+4) % 2) 
	  guts[i] += D[which] * (v[base+which+defs.ngenes] - v[base+which]);
									
	if ( (*tempctr >> defs.ngenes+defs.egenes+5) % 2) 
	  guts[i] += deriv[base+which];
									
	if ( (*tempctr >> defs.ngenes+defs.egenes+6) % 2) 
	  guts[i] += deriv[base+which] 
	          + lparm.lambda[which] * v[base+which] 
	          - D[which] * ((v[base+which-defs.ngenes] - v[base+which]) 
		              + (v[base+which+defs.ngenes] - v[base+which]));
								
	if ( (*tempctr >> defs.ngenes+defs.egenes+7) % 2) 
	  guts[i] += (deriv[base+which] 
		  + lparm.lambda[which] * v[base+which] 
		  - D[which] * ((v[base+which-defs.ngenes] - v[base+which]) 
			      + (v[base+which+defs.ngenes] - v[base+which])))
	          / lparm.R[which];

	tempctr++;
      }	
    }                                    
                                           /* last: posterior-most nucleus */
    tempctr = gutcomps;
    for (i=base1; i < base1 + numguts; i++) {
								
      if ( (*tempctr >> defs.ngenes+defs.egenes+2) % 2) 
	guts[i] += -lparm.lambda[which] * v[base+which];
								
      if ( (*tempctr >> defs.ngenes+defs.egenes+3) % 2) 
	guts[i] += D[which] * (v[base+which-defs.ngenes] - v[base+which]);
									
      if ( (*tempctr >> defs.ngenes+defs.egenes+5) % 2) 
	guts[i] += deriv[base+which];
									
      if ( (*tempctr >> defs.ngenes+defs.egenes+6) % 2) 
	guts[i] += deriv[base+which] 
	        + lparm.lambda[which] * v[base+which] 
	        - D[which] * (v[base+which-defs.ngenes] - v[base+which]);
								
      if ( (*tempctr >> defs.ngenes+defs.egenes+7) % 2) 
	guts[i] += (deriv[base+which] 
                + lparm.lambda[which] * v[base+which] 
	        - D[which] * (v[base+which-defs.ngenes] - v[base+which])) 
                / lparm.R[which];
								
      tempctr++;
    }	
  }
  
  free(deriv);

  if (defs.egenes > 0)
	  free(v_ext);

  return;
}



/*** ParseString: parses a line of the $gutsdefs section (dataformatX.X ****
 *                has the details); this function then returns two things: *
 *                1) the number of 'words' separated by whitespace in that *
 *                   string is the function's return value                 *
 *                2) the 'words' themselves are returned in arginp         *
 ***************************************************************************/

int ParseString(char *v, char **arginp)
{
  const int IN  = 1;                   /* are we inside or outside a word? */
  const int OUT = 0;

  int   state   = OUT;          /* state can be either 'IN' or 'OUT' above */
  int   n       = 0;                                  /* counter for words */
  char *mover   = v;        /* pointer-counter for pos within whole string */
  char **myargs = arginp;           /* pointer-counter for different words */
  char *tempo   = NULL;     /* pointer-counter for pos within current word */

/* parsing starts here */

  n = 0;
  state = OUT;
  while ((*mover != '\0') && (*mover != '\n')) {   /* loop thru whole line */
    if (*mover == ' ' || *mover == '\t') { /* whitespace -> done with word */
      if (state == IN) {                            /* word just finished: */
	*tempo = '\0';                       /* finish string by adding \0 */
	myargs++;                                       /* go to next word */
      }
      state = OUT;                           /* whitespace -> outside word */
    } else { 
      if (state == OUT) {                                    /* a new word */
	state   = IN;                           /* we're inside a word now */
	*myargs = (char *)calloc(MAX_RECORD, sizeof(char));
	tempo   = *myargs;      /* allocate and assign tempo to new string */
	++n;
      } 
      *tempo = *mover;	                         /* characters copied here */
      tempo++;                     /* tempo keeps track of pos within word */
    }		
    mover++;               /* mover keeps track of pos within whole string */
  }  

  if (state == IN) {       /* if we finished within a word, add a final \0 */
    *tempo = '\0';
    myargs++;
  }  

  return n;                    /* return the number of words in the string */
}
	



/*** GetGutsComps: takes a geneID string and an array of ID strings; then **
 *                 returns an array of long ints with bit flags set that   *
 *                 tell CalcGuts() what to calculate for each column of    *
 *                 guts output; this function also returns a pointer to    *
 *                 the current gene in the gene ID string                  *
 ***************************************************************************/

char *GetGutsComps(char *geneidstring, char *egeneidstring, 
					char **specsin, unsigned long *specsout)
{
  const char    *gutsids = "BALXYDJZU";                  /* legal guts IDs */

  char          *strcntr     = NULL;       /* ptr-counter within ID string */
  char          *ids         = NULL;          /* string with all legal IDs */
  char          **strarrcntr = NULL;     /* ptr-counter for gut ID strings */
  unsigned long *arrcntr     = NULL;           /* ptr-counter for gutcomps */
  char          *i           = NULL;  /* ptr for finding IDs in ids string */
	
	
/* ids = geneIDs plus gutsIDs */

  ids = (char *)calloc(MAX_RECORD, sizeof(char));


  ids = strcat(strcat(strcat(ids, geneidstring), egeneidstring), gutsids);


/* set pointer-counters to both ID strings and gutcomps */

  arrcntr = specsout;
  strarrcntr = specsin;
  strarrcntr++;                            /* skip first entry (gene name) */

/* the following loop parses the words of the ID string and sets according *
 * bits in the gutscomp array                                              */

  while (strcntr = *strarrcntr) {
    while ( *strcntr != '\0' ) {
      if ( !(i = strchr(ids, *strcntr)) ) 
	error("GetGutsComps: unknown gene in gutsdefs string");
      else if (*strcntr == 'U') {
	if ((*(strcntr+1) != '\0') || (strcntr != *strarrcntr) )
	  error("GetGutsComps: U not alone in gutsdefs string");
       	*arrcntr = (1 << (strlen(geneidstring)+strlen(egeneidstring)+2)) - 1;
      } else
	*arrcntr |= (1 << (i-ids));
      strcntr++;
    }
    arrcntr++;   	
    strarrcntr++;   	
  }

/* clean up and return a pointer to the gut gene in the gene ID string */
  free(ids);
  return strchr(geneidstring, **specsin);
}




/*** MUTATOR FUNCTIONS *****************************************************/

/*** Mutate: calls mutator functions according to genotype string **********
 ***************************************************************************/

void Mutate(char *g_type)
{
  int  i;
  int  c;

  char *record;

  lparm = CopyParm(parm);   /* make local copy of parameters to be mutated */

  record = g_type;
  c=(int)*record;

  for ( i=0; c != '\0'; i++, c=(int)*(++record) ) {
    if ( c == 'W' )
      continue;
    else if ( c == 'R' ) 
      R_Mutate(i);
    else if ( c == 'S' )
      RT_Mutate(i);
    else if ( c == 'T' )
      T_Mutate(i);
    else 
      error("Mutate: unrecognized letter in genotype string!");
  }
}



/*** T_Mutate: mutates genes by setting all their T matrix entries *********
 *             to zero. Used to simulate mutants that express a            *
 *             non-functional protein.                                     *
 ***************************************************************************/

void T_Mutate(int gene)
{
  int      i;

  for( i=0; i<defs.ngenes; i++)
    lparm.T[(i*defs.ngenes)+gene]=0;
}



/*** R_Mutate: mutates genes by setting their promotor strength R **********
 *             to zero, so there will be no transcription at all           *
 *             anymore. Used to simulate mutants that don't pro-           *
 *             duce any protein anymore.                                   *
 ***************************************************************************/

void R_Mutate(int gene)
{
    lparm.R[gene]=0;
}



/*** RT_Mutate: mutates gene by setting both promoter strength R ***********
 *              and T matrix entries to zero. Useful, if you want          *
 *              to be really sure that there is NO protein pro-            *
 *              duction and that maternal contribution also don't          *
 *              contribute anything to gene interactions.                  *
 ***************************************************************************/

void RT_Mutate(int gene)
{
  int      i;

  for( i=0; i<defs.ngenes; i++)
    lparm.T[(i*defs.ngenes)+gene]=0;
  lparm.R[gene]=0;
  /* Don't need to zero param.thresh in (trans acting) mutants */
}


/*** CopyParm: copies all the parameters into the lparm struct *************
 ***************************************************************************/

EqParms CopyParm(EqParms orig_parm)
{
  int           i,j;                                /* local loop counters */

  EqParms       l_parm;              /* copy of parm struct to be returned */

  l_parm.R = (double *)calloc(defs.ngenes, sizeof(double));
  l_parm.T = (double *)calloc(defs.ngenes * defs.ngenes, sizeof(double));

  if (defs.egenes > 0)
	  l_parm.E = (double *)calloc(defs.ngenes * defs.egenes, sizeof(double));

  l_parm.m = (double *)calloc(defs.ngenes, sizeof(double));
  l_parm.h = (double *)calloc(defs.ngenes, sizeof(double));
  if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) {
    l_parm.d = (double *)malloc(sizeof(double));
  } else {
    l_parm.d = (double *)calloc(defs.ngenes, sizeof(double));
  }
  l_parm.lambda = (double *)calloc(defs.ngenes, sizeof(double));
  l_parm.tau = (double *)calloc(defs.ngenes, sizeof(double));
  
  for (i=0; i<defs.ngenes; i++ ) {
    l_parm.R[i] = orig_parm.R[i];
    for (j=0; j<defs.ngenes; j++ )
      l_parm.T[(i*defs.ngenes)+j] = orig_parm.T[(i*defs.ngenes)+j];
    for (j=0; j<defs.egenes; j++ )
      l_parm.E[(i*defs.egenes)+j] = orig_parm.E[(i*defs.egenes)+j];
    l_parm.m[i] = orig_parm.m[i];
    l_parm.h[i] = orig_parm.h[i];
    l_parm.lambda[i] = orig_parm.lambda[i];
    l_parm.tau[i] = orig_parm.tau[i];
  }

  if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) {
    l_parm.d[0] = orig_parm.d[0];
  } else {
    for (i=0; i<defs.ngenes; i++) 
      l_parm.d[i] = orig_parm.d[i];
  }

  return l_parm;
}




 

/*** A FUNCTION THAT SETS STATIC STUFF IN ZYGOTIC.C ************************/

/*** SetRule: sets the static variable rule to MITOSIS or INTERPHASE *******
 ***************************************************************************/

void SetRule(int r) {
  rule = r;
}





/*** A FUNCTION THAT RETURN STATIC STUFF FROM ZYGOTIC.C ********************/

/*** GetParameters: returns the parm struct to the caller; note that this **
 *                  function returns the ORIGINAL PARAMETERS as they are   *
 *                  in the data file and NOT THE MUTATED ONES; this is im- *
 *                  portant to prevent limit violations in Score()         *
 ***************************************************************************/

EqParms *GetParameters(void)
{
  return &parm;
}



/*** GetMutParameters: same as above but returns the mutated copy of the ***
 *                     parameter struct; important for writing guts        *
 ***************************************************************************/

EqParms *GetMutParameters(void)
{
  return &lparm;
}




/*** A FUNCTION THAT READS PARAMETERS FROM THE DATA FILE INTO STRUCTS ******/

/*** ReadParamters: reads the parameters for a simulation run from the *****
 *                  eqparms or input section of the data file as indicated *
 *                  by the section_title argument and does the conversion  *
 *                  of protein half lives into lambda parameters.          *
 ***************************************************************************/

EqParms ReadParameters(FILE *fp, char *section_title)
{
  EqParms           l_parm;                 /* local copy of EqParm struct */
  double            *tempparm,*tempparm1;            /* temporary array to read parms */

  int               i;                               /* local loop counter */
  int               c;                         /* used to parse text lines */
  int               lead_punct;            /* flag for leading punctuation */
  int               linecount = 0;        /* keep track of # of lines read */
  int               Tcount    = 0;           /* keep track of T lines read */
  int               Ecount    = 0;           /* keep track of E lines read */

  char              *base;          /* pointer to beginning of line string */
  char              *record;    /* string for reading whole line of params */

  char              **fmt,**fmt1;   /* array of format strings for reading params */
  char              *skip,*skip1;               /* string of values to be skipped */

  const char        read_fmt[] = "%lg";                   /* read a double */
  const char        skip_fmt[] = "%*lg ";               /* ignore a double */

  base = (char *)calloc(MAX_RECORD, sizeof(char));

  skip = (char *)calloc(MAX_RECORD, sizeof(char));
  
  skip1 = (char *)calloc(MAX_RECORD, sizeof(char));
  
  fmt  = (char **)calloc(defs.ngenes, sizeof(char *));
  
  if (defs.egenes > 0)
	  fmt1  = (char **)calloc(defs.egenes, sizeof(char *));

  tempparm = (double *)calloc(defs.ngenes, sizeof(double));
  
  tempparm1 = (double *)calloc(defs.egenes, sizeof(double));

/* create format strings according to the number of genes */

  for ( i=0; i<defs.ngenes; i++ ) {  
    fmt[i] = (char *)calloc(MAX_RECORD, sizeof(char));   
    fmt[i] = strcpy(fmt[i], skip);
    fmt[i] = strcat(fmt[i], read_fmt);
    skip   = strcat(skip, skip_fmt);
 } 

/* create format strings according to the number of external inputs */

  for ( i=0; i<defs.egenes; i++ ) {  
    fmt1[i] = (char *)calloc(MAX_RECORD, sizeof(char));   
    fmt1[i] = strcpy(fmt1[i], skip1);
    fmt1[i] = strcat(fmt1[i], read_fmt);
    skip1   = strcat(skip1, skip_fmt);
 } 

/* initialize the EqParm struct */

  l_parm.R = (double *)calloc(defs.ngenes, sizeof(double));
  l_parm.T = (double *)calloc(defs.ngenes * defs.ngenes, sizeof(double));

  if (defs.egenes > 0)
	  l_parm.E = (double *)calloc(defs.ngenes * defs.egenes, sizeof(double));

  l_parm.m = (double *)calloc(defs.ngenes, sizeof(double));
  l_parm.h = (double *)calloc(defs.ngenes, sizeof(double));
  if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) {
    l_parm.d = (double *)malloc(sizeof(double));
  } else {
    l_parm.d = (double *)calloc(defs.ngenes, sizeof(double));
  }
  l_parm.lambda = (double *)calloc(defs.ngenes, sizeof(double));
  l_parm.tau = (double *)calloc(defs.ngenes, sizeof(double));


  fp = FindSection(fp, section_title);             /* find eqparms section */
  if( !fp )
    error("ReadParameters: cannot locate %s section", section_title);

  while ( strncmp(( base=fgets(base, MAX_RECORD, fp)), "$$", 2)) {
    
    record = base;
    lead_punct = 0;                                           /* of string */

    c = (int)*record;
    while ( c != '\0' ) {
      
      if ( isdigit(c) ) {                            /* line contains data */
	record = base;

/* usually read ngenes parameters, but for diff. schedule A or C only read *
 * one d parameter                                                         */

	  /* If no external inputs, we must be on the next set of parameters, therefore increase linecount by one. */
	if ((linecount == 2) && (defs.egenes == 0)) { 
			linecount++;
	}	  

	if ((linecount == 5) &&          
	     ((defs.diff_schedule=='A') || (defs.diff_schedule=='C'))) {
	  if ( 1 != sscanf(record, fmt[0], &tempparm[0]) )
	    error("ReadParameters: error reading parms");
	} else if (linecount == 2) {     
			  for ( i=0; i < defs.egenes; i++ ) { 
				if ( 1 != sscanf(record, fmt1[i], &tempparm1[i]) )
				  error("ReadParameters: error reading parms");
			  }
	} else {     
	  for ( i=0; i < defs.ngenes; i++ ) {
	    if ( 1 != sscanf(record, fmt[i], &tempparm[i]) )
	      error("ReadParameters: error reading parms");
	  }
	}

	switch (linecount) {  /* copy read parameters into the right array */
	case 0:
	  for ( i=0; i < defs.ngenes; i++ )                           /* R */
	    l_parm.R[i] = tempparm[i];
	  linecount++;
	  break;
	case 1:                 /* T: keep track of read lines with Tcount */
	  for ( i=0; i < defs.ngenes; i++ )
	    l_parm.T[i+Tcount*defs.ngenes] 
	      = tempparm[i];
	  Tcount++;
	  if ( Tcount == defs.ngenes ) 
	    linecount++;	  
	  break;
	case 2:                 /* E: keep track of read lines with Ecount */
	  if (defs.egenes > 0) {
		  for ( i=0; i < defs.egenes; i++ )
			l_parm.E[i+Ecount*defs.egenes] 
			  = tempparm1[i];
		  Ecount++;
		  if ( Ecount == defs.ngenes ) 
			linecount++;	  
	  }	  
	  break;
	case 3:                                                       /* m */
	  for ( i=0; i < defs.ngenes; i++ )
	    l_parm.m[i] = tempparm[i];
	  linecount++;
	  break;
	case 4:
	  for ( i=0; i < defs.ngenes; i++ )                           /* h */
	    l_parm.h[i] = tempparm[i];
	  linecount++;
	  break;
	case 5:                              /* d: consider diff. schedule */
	  if ((defs.diff_schedule == 'A') || (defs.diff_schedule == 'C' )) {
	    l_parm.d[0] = tempparm[0];
	  } else {
	    for ( i=0; i < defs.ngenes; i++ )
	      l_parm.d[i] = tempparm[i];
	  }
	  linecount++;
	  break;
	case 6:                                                  /* lambda */
	  for ( i=0; i < defs.ngenes; i++ ) {
	    l_parm.lambda[i] = tempparm[i];
	    l_parm.lambda[i] = log(2.) / l_parm.lambda[i]; 
          }                                        /* conversion done here */
	  linecount++;
	  break;
	case 7:                                                  /* tau */
	  for ( i=0; i < defs.ngenes; i++ ) 
	    l_parm.tau[i] = tempparm[i];
	  linecount++;
	  break;
	default:
	  error("ReadParameters: too many lines in parameter section");
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
	error("ReadParameters: illegal character in %s");
      }
    }
  }
  
  free(tempparm);
  free(tempparm1);
  free(base);
  free(skip);
  free(skip1);
  
  for (i=0; i<defs.ngenes; i++)
    free(fmt[i]);
  free(fmt);

  for (i=0; i<defs.egenes; i++)
    free(fmt1[i]);

  if (defs.egenes > 0)
	  free(fmt1);

  return l_parm;
}



/*** A FUNCTION THAT WRITES PARAMETERS INTO THE DATA FILE ******************/

/*** WriteParameters: writes the out_parm struct into a new section in the *
 *                    file specified by filename; the new 'eqparms' sec-   *
 *                    tion is inserted right after the 'input' section;    *
 *                    to achieve this, we need to write to a temporary     *
 *                    file which is then renamed to the output file name   *
 *              NOTE: lambdas are converted back into protein half lives!! *
 ***************************************************************************/

void WriteParameters(char *filename, EqParms *p, char *title, int ndigits)
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
    error("WriteParameters: error opening output file");

  temp = strcpy(temp,"parmXXXXXX");               /* required by mkstemp() */
  if ( mkstemp(temp) == -1 )              /* get unique name for temp file */
    error("WriteParameters: error creating temporary file");

  tmpfile = fopen(temp, "w");               /* ... and open it for writing */
  if ( !tmpfile )
    error("WriteParameters: error opening temporary file");

  if ( FindSection(outfile, title) ) { /* erase section if already present */
    fclose(outfile);                       /* this is a little kludgey but */
    KillSection(filename, title);      /* since KillSection needs to be in */
    outfile = fopen(filename, "r");           /* total control of the file */
  }
  rewind(outfile);

/* the follwoing two loops look for the appropriate file position to write */
/* the eqparms section (alternatives are input and eqparms)                */

  if ( !strcmp(title, "input") ) {
    while ( strncmp(record=fgets(record, MAX_RECORD, outfile), 
		    "$genotypes", 10) )
      fputs(record, tmpfile);
  } else if ( !strcmp(title, "eqparms") ) {
    while ( strncmp(record=fgets(record, MAX_RECORD, outfile), 
		    "$input", 6) )
      fputs(record, tmpfile);
  }
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

  PrintParameters(tmpfile, p, title, ndigits);

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
    error("WriteParameters: error renaming temp file %s");

  if ( remove(temp) )
    warning("WriteParameters: temp file %s could not be deleted", temp);

/* clean up */

  free(temp);
  free(record_ptr);
  free(saverec);
  free(shell_cmd);
}




/*** PrintParameters: prints an eqparms section with 'title' to the stream *
 *                    indicated by fp                                      *
 ***************************************************************************/

void PrintParameters(FILE *fp, EqParms *p, char *title, int ndigits)
{
  int    i, j;                                      /* local loop counters */
  double lambda_tmp;                           /* temporary var for lambda */


  fprintf(fp, "$%s\n", title);
  fprintf(fp, "promoter_strengths:\n");             /* Rs are written here */
  
  for ( i=0; i<defs.ngenes; i++ )  
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->R[i]);

  fprintf(fp, "\n");
  fprintf(fp, "genetic_interconnect_matrix:\n");        /* Ts written here */
  
  for ( i=0; i<defs.ngenes; i++ ) {
    for ( j=0; j<defs.ngenes; j++ ) 
      fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->T[(i*defs.ngenes)+j]);
    fprintf(fp, "\n");
  }

  fprintf(fp, "external_input_strengths:\n");        /* Es written here */
  
  for ( i=0; i<defs.ngenes; i++ ) {
    for ( j=0; j<defs.egenes; j++ ) 
      fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->E[(i*defs.egenes)+j]);
    fprintf(fp, "\n");
  }

  fprintf(fp, "maternal_connection_strengths:\n");      /* ms written here */

  for ( i=0; i<defs.ngenes; i++ ) 
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->m[i]);

  fprintf(fp, "\n");
  fprintf(fp, "promoter_thresholds:\n");            /* hs are written here */

  for ( i=0; i<defs.ngenes; i++ )
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->h[i]);

  fprintf(fp, "\n");
  fprintf(fp, "diffusion_parameter(s):\n");         /* ds are written here */
  
  if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') )
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->d[0]);
  else
    for ( i=0; i<defs.ngenes; i++ )
      fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->d[i]);
      
  fprintf(fp, "\n");
  fprintf(fp, "protein_half_lives:\n");        /* lambdas are written here */
  
  for ( i=0; i<defs.ngenes; i++ ) {
    lambda_tmp = log(2.) / p->lambda[i];           /* conversion done here */
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, lambda_tmp);
  }
  
  fprintf(fp, "\n");
  fprintf(fp, "translational_transcriptional_delays:\n");        /* taus are written here */
  
  for ( i=0; i<defs.ngenes; i++ ) 
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->tau[i]);
  
  fprintf(fp, "\n$$\n");
  fflush(fp);

}


/* Function given the starting (fb) and ending (fa) time structures,
calculates the use cpu time in seconds */

double tvsub(struct rusage a, struct rusage b)
{
  double fa, fb ;

    fa = a.ru_utime.tv_sec + a.ru_stime.tv_sec +
	    (a.ru_utime.tv_usec + a.ru_stime.tv_usec) / 1e6 ;
	fb = b.ru_utime.tv_sec + b.ru_stime.tv_sec +
		(b.ru_utime.tv_usec + b.ru_stime.tv_usec) / 1e6 ;

	return fa - fb ;
}  

/*** DvdtDelay: the delay derivative function; implements the equations **
 *             as discussed in the delay report,                    *
 *             plus different g(u) functions.					    *
 ***************************************************************************/

void DvdtDelay(double *v, double **vd, double t, double *vdot, int n)
{
    double vinput1 = 0;     
  int        m;                                        /* number of nuclei */
  int        ap;                 /* nuclear index on AP axis [0,1,...,m-1] */
  int        i,j;                                   /* local loop counters */
  int        k;                   /* index of gene k in a specific nucleus */
  int        base, base1;            /* index of first gene in a specific nucleus */
  int 		 *l_rule;		/* for autonomous implementation */
#ifdef ALPHA_DU
  int        incx = 1;        /* increment step size for vsqrt input array */
  int        incy = 1;       /* increment step size for vsqrt output array */
#endif

  static int num_nucs  = 0;      /* store the number of nucs for next step */
  static int bcd_index = 0;   /* the *next* array in bicoid struct for bcd */
  static DArrPtr bcd;              /* pointer to appropriate bicoid struct */
  double **v_ext;			/* array to hold the external input
									  concentrations at time t */


/* get D parameters and bicoid gradient according to cleavage cycle */

  m = n / defs.ngenes;                        /* m is the number of nuclei */
  if (m != num_nucs) {      /* time-varying quantities only vary by ccycle */
    GetD(t, lparm.d, defs.diff_schedule, D);  
                      /* get diff coefficients, according to diff schedule */
    if (num_nucs > m)                  /* started a new iteration in score */
      bcd_index = 0;                           /* -> start all over again! */
    num_nucs = m;                         /* store # of nucs for next step */
    bcd = GetBicoid(t, genindex);                   /* get bicoid gradient */
    if( bcd.size != num_nucs) 
     error("DvdtDelay: %d nuclei don't match Bicoid!", num_nucs);
    bcd_index++;                   /* store index for next bicoid gradient */
  }

  l_rule = (int *) calloc(defs.ngenes, sizeof(int));
  
  for (i = 0; i < defs.ngenes; i++)
 	l_rule[i] = !Theta(t - lparm.tau[i]);/*for autonomous equations */

/* Here we retrieve the external input concentrations into v_ext */
/* 01/13/10 Manu: If there are no external inputs, don't
 * allocate v_ext or populate it with external input
 * concentrations */  


    if (defs.egenes > 0) {
		v_ext = (double **) calloc(defs.ngenes, sizeof(double *));
		
		for (i = 0; i < defs.ngenes; i++) {
			
			v_ext[i] = (double *) calloc(m*defs.egenes, sizeof(double));
			ExternalInputs(t - lparm.tau[i], t, v_ext[i], m*defs.egenes);
		}
	}

/* This is how it works (by JR): 
         
    ap      nucleus position on ap axis
    base    index of first gene (0) in current nucleus
    k       index of gene k in current nucleus 

            Protein synthesis terms are calculated according to g(u)

            First we do loop for vinput contributions; vinput contains
	    the u that goes into g(u)

	    Then we do a separate loop or vector func for sqrt or exp

	    Then we do vdot contributions (R and lambda part)

	    Then we do the spatial part (Ds); we do a special case for 
            each end 

            Note, however, that before the real loop, we have to do a 
            special case for when rhere is one nuc, hence no diffusion

	    These loops look a little funky 'cause we don't want any 
            divides'                                                       */

	/* 01/13/10 Manu: If there are no external inputs, the for
	 * loops for updating vinput1 with external input terms
	 * are never entered */


/***************************************************************************
 *                                                                         *
 *        g(u) = 1/2 * ( u / sqrt(1 + u^2) + 1)                            *
 *                                                                         *
 ***************************************************************************/
    if ( gofu == Sqrt ) {

      for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
	for (i=base; i < base + defs.ngenes; i++) {

	  k = i - base;
	  
	  vinput1  = lparm.h[k];
	  vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

	  for(j=0; j < defs.egenes; j++)
	    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[k][base1 + j];

	  for(j=0; j < defs.ngenes; j++)
	    vinput1 += lparm.T[(k*defs.ngenes)+j] * vd[k][base + j];

	  bot2[i]   = 1 + vinput1 * vinput1;
	  vinput[i] = vinput1;
	}
      }

                            /* now calculate sqrt(1+u2); store it in bot[] */
#ifdef ALPHA_DU
      vsqrt_(bot2,&incx,bot,&incy,&n);    /* superfast DEC vector function */
#else
      for(i=0; i < n; i++)                  /* slow traditional style sqrt */
	bot[i] = sqrt(bot2[i]);
#endif
            /* next loop does the rest of the equation (R, Ds and lambdas) */
                                                 /* store result in vdot[] */

      if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

	register double vdot1, g1;

	for( base=0; base<n; base+=defs.ngenes ) {
	  for( i=base; i<base+defs.ngenes; i++ ) { 
	  
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      = 1 + vinput[i]/bot[i];         
	    vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;     
	    vdot[i] = vdot1; 
	  }
	}
	
      } else {                    /* then for multiple nuclei -> diffusion */

	register double vdot1,g1;

	for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
	  k       = i;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      = 1 + vinput[i]/bot[i];
	  vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;
	  vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}    
                                                     /* then middle nuclei */
	for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
	  for(i=base; i < base + defs.ngenes; i++){
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      =  1 + vinput[i]/bot[i];
	    vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;
	    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) + 
			      (v[i + defs.ngenes] - v[i]) );
	    vdot[i] = vdot1;
	  }
	}                                    
	                                   /* last: posterior-most nucleus */
	for(i=base; i < base + defs.ngenes; i++){
	  k       = i - base;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      =  1 + vinput[i]/bot[i];
	  vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;
	  vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}
      }

/***************************************************************************
 *                                                                         *
 *        g(u) = 1/2 * (tanh(u) + 1) )                                     *
 *                                                                         *
 ***************************************************************************/

    } else if ( gofu == Tanh ) {

      for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
	for(i=base; i < base + defs.ngenes; i++){
	  
	  k = i - base;
	  
	  vinput1  = lparm.h[k];
	  vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

	  for(j=0; j < defs.egenes; j++)
	    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[k][base1 + j];

	  for(j=0; j < defs.ngenes; j++)
	    vinput1  += lparm.T[(k*defs.ngenes)+j] * vd[k][base + j];

	  vinput[i] = vinput1;
	}
      }

            /* next loop does the rest of the equation (R, Ds and lambdas) */
                                                 /* store result in vdot[] */

      if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

	register double vdot1, g1;

	for( base=0; base<n; base+=defs.ngenes ) {
	  for( i=base; i<base+defs.ngenes; i++ ) { 
	  
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      = tanh(vinput[i]) + 1;         
	    vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;     
	    vdot[i] = vdot1; 
	  }
	}

      } else {                    /* then for multiple nuclei -> diffusion */

	register double vdot1,g1;
	
	for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
	  k       = i;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      = tanh(vinput[i]) + 1;
	  vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;
	  vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}     
	                                             /* then middle nuclei */
	for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
	  for(i=base; i < base + defs.ngenes; i++){
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      = tanh(vinput[i]) + 1;
	    vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;
	    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) + 
			      (v[i + defs.ngenes] - v[i]) );
	    vdot[i] = vdot1;
	  }
	}                                    
                                           /* last: posterior-most nucleus */
	for(i=base; i < base + defs.ngenes; i++){
	  k       = i - base;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      = tanh(vinput[i]) + 1;
	  vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;
	  vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}
      }

/***************************************************************************
 *                                                                         *
 *        g(u) = 1 / (1 + exp(-2u))                                        *
 *                                                                         *
 ***************************************************************************/

    } else if ( gofu == Exp ) {

      for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
	for(i=base; i < base + defs.ngenes; i++){

	  k = i - base;

	  vinput1  = lparm.h[k];
	  vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

	  for(j=0; j < defs.egenes; j++)
	    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[k][base1 + j];

	  for(j=0;  j < defs.ngenes; j++)
	    vinput1  += lparm.T[(k*defs.ngenes)+j] * vd[k][base + j];

	  vinput[i] = -2.0 * vinput1;	
	}
      }

                               /* now calculate exp(-u); store it in bot[] */
#ifdef ALPHA_DU
      vexp_(vinput,&incx,bot,&incy,&n);   /* superfast DEC vector function */
#else
      for (i=0; i<n; i++)                    /* slow traditional style exp */
	bot[i] = exp(vinput[i]);
#endif

            /* next loop does the rest of the equation (R, Ds and lambdas) */
                                                 /* store result in vdot[] */

      if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

	register double vdot1, g1;

	for( base=0; base<n; base+=defs.ngenes ) {
	  for( i=base; i<base+defs.ngenes; i++ ) { 
	  
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      = 1 / (1 + bot[i]);         
	    vdot1  += l_rule[k]*lparm.R[k] * g1;     
	    vdot[i] = vdot1; 
	  }
	}

      } else {                    /* then for multiple nuclei -> diffusion */

	register double vdot1,g1;

	for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
	  k       = i;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      = 1 / (1 + bot[i]);
	  vdot1  += l_rule[k]*lparm.R[k] * g1;
	  vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}     
                                                     /* then middle nuclei */
	for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
	  for(i=base; i < base + defs.ngenes; i++){
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    g1      = 1 / (1 + bot[i]);
	    vdot1  += l_rule[k]*lparm.R[k] * g1;
	    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) + 
			      (v[i + defs.ngenes] - v[i]) );
	    vdot[i] = vdot1;
	  }
	}                                    
                                           /* last: posterior-most nucleus */
	for(i=base; i < base + defs.ngenes; i++){
	  k       = i - base;
	  vdot1   = -lparm.lambda[k] * v[i];
	  g1      = 1 / (1 + bot[i]);
	  vdot1  += l_rule[k]*lparm.R[k]  * g1;
	  vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}
      }

/***************************************************************************
 *                                                                         *
 *        g(u) = 0 if u<0, 1 if u>=0 (Heaviside function)                  *
 *                                                                         *
 *        this makes the model quasi boolean and the equations locally     *
 *        linear for both u>0 and u<0                                      *
 ***************************************************************************/

    } else if ( gofu == Hvs ) {

      for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
	for(i=base; i < base + defs.ngenes; i++){

	  k = i - base;

	  vinput1  = lparm.h[k];
	  vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

	  for(j=0; j < defs.egenes; j++)
	    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[k][base1 + j];

	  for(j=0; j < defs.ngenes; j++)
	    vinput1  += lparm.T[(k*defs.ngenes)+j] * vd[k][base + j];
	  vinput[i] = vinput1;
	}
      }

     /* next loop does the rest of the equation (R, Ds and lambdas) */
                                                 /* store result in vdot[] */

      if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

	register double vdot1, g1;

	for( base=0; base<n; base+=defs.ngenes ) {
	  for( i=base; i<base+defs.ngenes; i++ ) { 
	    
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    if (vinput[i] >= 0.) 
	      g1 = 1.;
	    else 
	      g1 = 0.;         
	    vdot1  += l_rule[k]*lparm.R[k] * g1;     
	    vdot[i] = vdot1; 
	  }
	}

      } else {                    /* then for multiple nuclei -> diffusion */

	register double vdot1,g1;

	for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
	  k       = i;
	  vdot1   = -lparm.lambda[k] * v[i];
	  if (vinput[i] >= 0.) 
	    g1 = 1.;
	  else 
	    g1 = 0.;
	  vdot1  += l_rule[k]*lparm.R[k] * g1;
	  vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}     
	                                             /* then middle nuclei */
	for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
	  for(i=base; i < base + defs.ngenes; i++){
	    k       = i - base;
	    vdot1   = -lparm.lambda[k] * v[i];
	    if (vinput[i] >= 0.) 
	      g1 = 1.;
	    else 
	      g1 = 0.;
	    vdot1  += l_rule[k]*lparm.R[k] * g1;
	    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) + 
			      (v[i + defs.ngenes] - v[i]) );
	    vdot[i] = vdot1;
	  }
	}                                    
	                                   /* last: posterior-most nucleus */
	for(i=base; i < base + defs.ngenes; i++){
	  k       = i - base;
	  vdot1   = -lparm.lambda[k] * v[i];
	  if (vinput[i] >= 0.) 
	    g1 = 1.;
	  else 
	    g1 = 0.;
	  vdot1  += l_rule[k]*lparm.R[k] * g1;
	  vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
	  vdot[i] = vdot1;
	}
      }

    } else
      error("DvdtDelay: unknown g(u)");
    
/* during mitosis only diffusion and decay happen */
    
	free(l_rule);

    if (defs.egenes > 0) {
		for (i = 0; i < defs.ngenes; i++)
			free(v_ext[i]);
			
		free(v_ext);
	}	

  return;
}

