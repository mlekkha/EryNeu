/*****************************************************************
 *                                                               *
 *   maternal.c                                                  *
 *                                                               *
 *****************************************************************
 *                                                               *
 *   written by JR, modified by Yoginho                          *
 *                                                               *
 *****************************************************************
 *                                                               *
 * This file contains several kinds of things:                   *
 *                                                               *
 * 1. some small I/O functions (FindSection(), KillSection())    *
 *    that are used throughout in the fly code                   *
 * 2. stuff that deals with that part of the blastoderm which is *
 *    fixed by the maternal genotype, i.e. functions dealing     *
 *    with division schedules (including the rules for zygotic.c,*
 *    the number of nuclei at each cleavage cycle and the li-    *
 *    neage number of the most anterior nucleus for each clea-   * 
 *    vage cycle), bicoid gradients and diffusion schedules      *
 * 3. bias-related stuff is also in here since it's needed to    *
 *    run the model                                              *
 *                                                               *
 *****************************************************************
 *                                                               *
 * Copyright (C) 1989-2003 John Reinitz                          *
 * the full GPL copyright notice can be found in lsa.c           *
 *                                                               *
 *****************************************************************/


#include <float.h>                                      /* for DBL_EPSILON */
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>   

#include <error.h>                 /* for error and linked list functions  */
#include <maternal.h>
#include <integrate.h>			/* for HALF_EPSILON */
#include <zygotic.h>                      /* for derivative function rules */



/*** MITOSIS SCHEDULE: hard-wired cell division tables *********************
 *                                                                         *
 * From Foe & Alberts, '83 (Table 1, TOTAL ELAPSED added by JR and JJ):    *
 *                                                                         *
 * The authors observed anterior tips of growing embryos during cycles 10  *
 * to 14 and measured the time of somatic bud cycles and recorded the      *
 * visibility of the nuclei, which indicate the duration of mitoses.       * 
 *                                                                         *
 * TOTAL ELAPSED is the total time which elapsed since the start of a sim- *
 * ulation. I.e. we assume                                                 *
 *                                                                         *
 * t=0 at end of cell cycle 10 for NDIVS=3                                 *
 * t=0 1 min into interphase of cell cycle 10 for NDIVS=4 (see (*) below)  *
 *                                                                         *
 * Times are given in minutes, standard devs are between parentheses)      *
 *                                                                         *
 * CYCLE  TOTAL DURATION      NO NUCLEUS            TOTAL ELAPSED          *
 *                                              NDIVS=3       NDIVS=4      *
 *                                                                         *
 *  10 (*)  7.8 (0.6)          3.3 (0.9)                        0.0        *
 *  11      9.5 (0.7)          3.0 (0.9)          0.0           7.8        *
 *  12     12.4 (0.9)          3.3 (0.9)          9.5          17.3        * 
 *  13     21.1 (1.5)          5.1 (0.9)         21.9          29.7        *
 *     	             	       				                   *
 *  14        50+          mitosis post gast     43.0          50.8        *
 *                                                                         *
 *  gastrulation time:                           93.0         100.8        *
 *                                                                         *
 * (*) NOTE: cell cycle 10 actually takes 8.8 minutes, but the migration   *
 *           of nuclei to the cell surface is only completed 1 min into    *
 *           interphase of cell cycle 10. This is our simulation starting  *
 *           point for NDIVS=4.                                            * 
 *                                                                         *
 * The following arrays are the hard-wired stuff about times and durations *
 * of cell divisions and time of gastrulation; each function in maternal.c *
 * that returns information about timing in the blastoderm needs to choose *
 * the appropriate tables depending on the problem, i.e. the mitosis sche- *
 * dule used.                                                              *
 *                                                                         *
 * Oldstyle division times were a coarse approximation based on the 4 min  *
 * time units intially used by JR and DS. They were refined into 4/3 min   *
 * units in pre 9.2 code.                                                  *
 *                                                                         *
 * IMPORTANT: Post 9.2 code does NOT require rounded division times any-   *
 * more. These old division times are only included for backward compati-  *
 * bility and should not be used for annealing to new data!                *
 *                                                                         *
 * NOTE: The times must be in reverse order!!!                             *
 * ALSO NOTE: divtimes are at the END of the cell division!!!              *
 * LAST NOTE: there are 3 arrays each, which are (in order of appearance): *
 *            - oldstyle (see above)                                       *
 *            - 3-celldivision schedule (starting at cycle 11)             *
 *            - 4-celldivision schedule (starting at cycle 10)             *
 *                                                                         *
 * JJ (04/04/02):                                                          *
 *            I've added a 'zero'-division schedule which we will use to   *
 *            check if we can get patterns without cell divisions; since   *
 *            this division schedule doesn't have any divisions, it also   *
 *            doesn't need any division tables below                       *
 *                                                                         *
 ***************************************************************************/


/* division times: at end of mitosis! */

static const double old_divtimes[3]     = {52.0, 28.0, 12.0};

static const double divtimes1[1]        = {21.1};
static const double divtimes2[2]        = {33.5, 12.4};
static const double divtimes3[3]        = {43.0, 21.9,  9.5};
static const double divtimes4[4]        = {50.8, 29.7, 17.3,  7.8};

/* division durations */

static const double old_div_duration[3] = { 4.0,  4.0,  4.0};

static const double div_duration1[1]    = { 5.1};
static const double div_duration2[2]    = { 5.1,  3.3};
static const double div_duration3[3]    = { 5.1,  3.3,  3.0};
static const double div_duration4[4]    = { 5.1,  3.3,  3.0,  3.3};

/* gastrulation times */

static const double old_gast_time       = 88.;

/*** FIX gast_time0**/
//static const double gast_time0          = 50.;   
static const double gast_time0          = 170.;
static const double gast_time1          = 71.1;
static const double gast_time2          = 83.5;     
static const double gast_time3          = 93.;     
static const double gast_time4          = 100.8;

/* full division times: including t<0 */

static const double full_divtimes0[TOTAL_DIVS]        = {0.0, -21.1,
-33.5, -43.0, -51.8, -57.8};
static const double full_divtimes1[TOTAL_DIVS]        = {21.1, 0.0, -12.4, -21.9, -30.7, -36.7};
static const double full_divtimes2[TOTAL_DIVS]        = {33.5, 12.4, 0.0, -9.5, -18.3, -24.3};
static const double full_divtimes3[TOTAL_DIVS]        = {43.0, 21.9,  9.5, 0.0, -8.8, -14.8};
static const double full_divtimes4[TOTAL_DIVS]        = {50.8, 29.7, 17.3, 7.8, -1.0, -7.0};

/* division durations */

static const double full_div_durations[TOTAL_DIVS]    = { 5.1,  3.3,  3.0, 3.3, 3.0, 3.0};





/*** STATIC VARIABLES ******************************************************/

/* following is number of nucs in each cleavage cycle, reverse order */

static int         *nnucs;
static int		   *full_nnucs=NULL;

/* following contains lineage numbers at which nuclei at each ccycle start */

static int         *lin_start; 
static int		   *full_lin_start=NULL;
static int			full_ccycles;

/* The following two store bicoid gradients and bias static to maternal.c  */
/* bias is found here because it contains the maternal contributions to    */
/* some zygotic genes (such as hunchback), but it can also be used to add  */
/* heatshocks during the simulation                                        */

static GenoType    *bcdtype;   
static GenoType    *biastype;

/* these two static structs are used for storing the times for which       */
/* there is bias; these times can be retrieved by using GetBTimes          */

static GenoType    *bt;                    /* bias times for each genotype */

static int         bt_init_flag = 0;                   /* flag for BTtable */
static int         d_flag       = 0;        /* flag for first call to GetD */
static int         rule_flag    = 0;     /* flag for first call to GetRule */
static int         theta_flag    = 0;     /* flag for first call to
Theta */






/*** INITIALIZATION FUNCTIONS **********************************************/

/*** InitBicoid: Copies the Blist read by ReadBicoid into the DArrPtr ******
 *               structure; the bicoid DArrPtr contains pointers to bcd    *
 *               arrays for each cell division cycle                       *
 ***************************************************************************/

void InitBicoid(FILE *fp)
{
  int               i;                               /* local loop counter */

  Blist             *inlist;              /* temporary linked list for bcd */

  Slist             *genotypes;         /* temporary linked list for geno- */
  Slist             *current;                      /* types from data file */


  genotypes = ReadGenotypes(fp);
  if ( nalleles == 0 )
    nalleles  = count_Slist(genotypes);

  if ( !(bcdtype=(GenoType *)calloc(nalleles, sizeof(GenoType))) )
    error("InitBicoid: Could not allocate bcdtype struct");

/*** for loop: read bicoid for each genotype *******************************/

  for (current=genotypes, i=0; current; current=current->next, i++) {

    if ( !(inlist=ReadBicoid(fp, current->bcd_section)) )   /* read bicoid */
      error("InitBicoid: error reading %s", current->bcd_section);
    else {

      if ( !(bcdtype[i].genotype = (char *)calloc(MAX_RECORD, sizeof(char))) )
	error("InitBicoid: could not allocate bcd genotype string");
      bcdtype[i].genotype   = strcpy(bcdtype[i].genotype, current->genotype);
      bcdtype[i].ptr.bicoid = List2Bicoid(inlist);

      free_Blist(inlist);
    }
  }
  free_Slist(genotypes);
}



/*** InitBias:  puts bias records in a form where get_bias can use them; ***
 *              it expects times in increasing order; it expects a non-    *
 *              sparse entry, with no genes or nuclei missing              *
 ***************************************************************************/

void InitBias(FILE *fp)
{
  int       i;                                            /* loop counters */

  Dlist     *inlist;                     /* temporary linked list for bias */
                                   
  Slist     *genotypes;                 /* temporary linked list for geno- */
  Slist     *current;                              /* types from data file */

  int       ndp = 0;  /* dummy for ReadData, no need to count datapts here */
  
  genotypes = ReadGenotypes(fp);
  if ( nalleles == 0 )
    nalleles  = count_Slist(genotypes);
  
  if ( !(biastype=(GenoType *)calloc(nalleles, sizeof(GenoType))) )
    error("InitBias: Could not allocate biastype struct");

/*** for loop: read bicoid for each genotype *******************************/

  for (current=genotypes, i=0; current; current=current->next, i++) {
 
    if( !(inlist=ReadData(fp, current->bias_section, &ndp)) ) /* read bias */
      error("InitBias: error reading %s", current->bias_section);
    else {
      
      if (!(biastype[i].genotype=(char *)calloc(MAX_RECORD, sizeof(char))))
	error("InitBias: could not allocate bias genotype string");      
      biastype[i].genotype = strcpy(biastype[i].genotype, current->genotype);
      biastype[i].ptr.bias = List2Bias(inlist);
      
      free_Dlist(inlist);
    }
  }
    
  InitBTs();                     /* initialize the static bias time struct */
  free_Slist(genotypes);
}

	

/*** InitBTs: initializes the static BT struct that holds all times for ****
 *            which we have bias.                                          *
 ***************************************************************************/

void InitBTs(void)
{
  int       i, j;
 
  if ( !bt_init_flag ) {
    if ( !(bt=(GenoType *)calloc(nalleles, sizeof(GenoType))) )
      error("InitBTs: could not allocate bt struct");
    
    for (i=0; i<nalleles; i++) {
      
      if ( !(bt[i].genotype = (char *)calloc(MAX_RECORD, sizeof(char))) )
	error("InitBTs: could not allocate BT genotype string");
      bt[i].genotype       = strcpy(bt[i].genotype, biastype[i].genotype);
      bt[i].ptr.times.size = biastype[i].ptr.bias.size;
      if ( !(bt[i].ptr.times.array = 
	   (double *)calloc(bt[i].ptr.times.size, sizeof(double))) )
	error("InitBTs: could not allocate bt array");
      
      for (j=0; j<bt[i].ptr.times.size; j++)
	bt[i].ptr.times.array[j] = biastype[i].ptr.bias.array[j].time;
    } 
    
    bt_init_flag = 1;
  }
}
  


/*** InitNNucs: takes the global defs.nnucs and calculates number of nucs **
 *              for each cleavage cycle which are then stored in reverse   *
 *              order in the static nnucs[] array                          *
 *   CAUTION:   defs struct and lin_start need to be initialized before!   *
 ***************************************************************************/

void InitNNucs(void)
{
  int           i;                                         /* loop counter */
  int           n; /* used to calculate number of nucs for each cell cycle */

  if ( !(nnucs = (int *)calloc(defs.ndivs+1, sizeof(int))) )
    error("InitNNucs: could not allocate nnucs array");

  n = defs.nnucs;

/* below we have to take into account two cases: a) most anterior lineage  *
 * number is odd-numbered -> always add a nucleus to the earlier cycle or  *
 * b) most anterior lineage number is even-numbered: just add an additio-  *
 * nal nucleus if last nucleus is odd-numbered (see also exhaustive com-   *
 * ments about this at the DIVIDE rule in Blastoderm() in integrate.c)     */

  for(i=0; i<=defs.ndivs; i++) {
    nnucs[i] = n;
    if ( lin_start[i] % 2 )
      n = n/2 + 1;
    else
      n = (n % 2) ? n/2 + 1 : n/2;
  }

}

void InitFullNNucs(void)
{
  int           i;                                         /* loop counter */
  int           n; /* used to calculate number of nucs for each cell cycle */

  if (full_nnucs)
	return;
  
  if ( !(full_nnucs = (int *)calloc(full_ccycles, sizeof(int))) )
    error("InitFullNNucs: could not allocate full_nnucs array");

  n = defs.nnucs;

/* below we have to take into account two cases: a) most anterior lineage  *
 * number is odd-numbered -> always add a nucleus to the earlier cycle or  *
 * b) most anterior lineage number is even-numbered: just add an additio-  *
 * nal nucleus if last nucleus is odd-numbered (see also exhaustive com-   *
 * ments about this at the DIVIDE rule in Blastoderm() in integrate.c)     */

  for(i=0; i < full_ccycles; i++) {
    full_nnucs[i] = n;
    if ( full_lin_start[i] % 2 )
      n = n/2 + 1;
    else
      n = (n % 2) ? n/2 + 1 : n/2;
  }

/*  for (i=0; i < full_ccycles; i++)
  	printf("History lineages %d, nnucs %d\n", full_lin_start[i],full_nnucs[i]);*/
  

}






/*** FUNCTIONS THAT RETURN INFO ABOUT THE EMBRYO **************************/

/*** GetBicoid: returns a bicoid gradients (in form of a DArrPtr) for a ****
 *              specific time and genotype.                                *
 ***************************************************************************/

DArrPtr GetBicoid(double time, int genindex)
{
  int             i;
  unsigned int    ccycle;
  
  if ( genindex < 0 || genindex >= nalleles )   
    error("GetBicoid: invalid genotype index %d", genindex);

  ccycle = GetCCycle(time);
  
  for (i=0; i <= bcdtype[genindex].ptr.bicoid.size; i++) 
    if ( ccycle == bcdtype[genindex].ptr.bicoid.array[i].ccycle )
      return bcdtype[genindex].ptr.bicoid.array[i].gradient;
  
  error("GetBicoid: no bicoid gradient for ccycle %d", ccycle);
  
  return bcdtype[genindex].ptr.bicoid.array[i].gradient;   
                                       /* just to make the compiler happy! */
}



/*** GetBias: This function returns bias values for a given time and *******
 *            genotype.                                                    *
 ***************************************************************************/


DArrPtr GetBias(double time, int genindex)
{
  int      i;

  if ( genindex < 0 || genindex >= nalleles )   
    error("GetBias: invalid genotype index %d", genindex);

  for(i=0; i < biastype[genindex].ptr.bias.size; i++) {
    if (biastype[genindex].ptr.bias.array[i].time == time)
      break;
  }

  return biastype[genindex].ptr.bias.array[i].state;
}



/*** GetBTimes: returns a sized array of times for which there is bias *****
 ***************************************************************************/

DArrPtr GetBTimes(char *genotype)
{
  int index;                                               /* loop counter */

  for(index=0; index<nalleles; index++) 
    if ( !(strcmp(biastype[index].genotype, genotype)) )
      break;

/* if no explicit bias times for this genotype -> use wt bias times */

  if ( index == nalleles ) 
    index = 0;

/* check if we actually have biastimes at all or otherwise -> error! */

  if ( bt_init_flag ) 
     return bt[index].ptr.times;
   else
     error("GetBTimes: called without initialized BTs");

  return bt[index].ptr.times;                          
                                       /* just to make the compiler happy */ 
}
  


/*** GetNNucs: reads number of nuclei for a given time *********************
 ***************************************************************************/

int GetNNucs(double t)
{
  int      i;                                              /* loop counter */
  double   *table;                         /* local copy of divtimes table */
  
/* assign 'table' to the appropriate division schedule */

  if ( olddivstyle ) {
    if ( defs.ndivs != 3 ) 
      error("GetNNucs: only 3 cell divisions allowed for oldstyle (-o)");
    table = (double *)old_divtimes;
  } else 
    if ( defs.ndivs == 0 )
      return nnucs[0];
    else if ( defs.ndivs == 1 )
      table = (double *)divtimes1;
    else if ( defs.ndivs == 2 )
      table = (double *)divtimes2;
    else if ( defs.ndivs == 3 )
      table = (double *)divtimes3;
    else if ( defs.ndivs == 4 ) 
      table = (double *)divtimes4;
    else
      error("GetNNucs: can't handle %d cell divisions!", defs.ndivs);

/* evaluate nnucs for current time; note that for the *exact* time of cell *
 * division, we'll return the number of nuclei before the division has ac- *
 * tually occurred                                                         */

  for (i=0; i<defs.ndivs; i++)              
    if ( t>table[i] )
      return nnucs[i];

  return nnucs[i];
}



/*** GetStartLin: returns the lineage number of the most anterior nucleus **
 *                for a given time                                         *
 ***************************************************************************/

int GetStartLin(double t)
{
  int      i;                                              /* loop counter */
  double   *table;                         /* local copy of divtimes table */
  
/* assign 'table' to the appropriate division schedule */

  if ( olddivstyle ) {
    if ( defs.ndivs != 3 ) 
      error("GetStartLin: only 3 cell divisions allowed for oldstyle (-o)");
    table = (double *)old_divtimes;
  } else 
    if ( defs.ndivs == 0 )
      return lin_start[0];
    else if ( defs.ndivs == 1 )
      table = (double *)divtimes1;
    else if ( defs.ndivs == 2 )
      table = (double *)divtimes2; 
    else if ( defs.ndivs == 3 )
      table = (double *)divtimes3;
    else if ( defs.ndivs == 4 ) 
      table = (double *)divtimes4;
    else
      error("GetStartLin: can't handle %d cell divisions!", defs.ndivs);
  
/* evaluate lineage number of most anterior nucleus for current time; note *
 * that for the *exact* time of cell division, we'll return the lineage    *
 * number of the most anterior nucleus of the previous cell cycle          */

  for (i=0; i<defs.ndivs; i++)              
    if ( t>table[i] )
      return lin_start[i];
  
  return lin_start[i];
}

int Index2StartLin(int index)
{

	return full_lin_start[index];

}

int Index2NNuc(int index)
{

	return full_nnucs[index];

}

/*** GetStartLinIndex: returns the index of the lineage array *                for a given time                                         *
 ***************************************************************************/

int GetStartLinIndex(double t)
{
  int      i;                                              /* loop counter */
  double   *table;                         /* local copy of divtimes table */
  
/* assign 'table' to the appropriate division schedule */

  if ( olddivstyle ) {
    if ( defs.ndivs != 3 ) 
      error("GetStartLin: only 3 cell divisions allowed for oldstyle (-o)");
    table = (double *)old_divtimes;
  } else 
    if ( defs.ndivs == 0 )
      table = (double *)full_divtimes0;
    else if ( defs.ndivs == 1 )
      table = (double *)full_divtimes1;
    else if ( defs.ndivs == 2 )
      table = (double *)full_divtimes2; 
    else if ( defs.ndivs == 3 )
      table = (double *)full_divtimes3;
    else if ( defs.ndivs == 4 ) 
      table = (double *)full_divtimes4;
    else
      error("GetStartLin: can't handle %d cell divisions!", defs.ndivs);
  
/* evaluate lineage number of most anterior nucleus for current time; note *
 * that for the *exact* time of cell division, we'll return the lineage    *
 * number of the most anterior nucleus of the previous cell cycle          */

  for (i=0; i<full_ccycles-1; i++)              
    if ( t > table[i]+HALF_EPSILON )
      return i;
  
  return i;
}




/*** GetCCycle: returns cleavage cycle number for a given time *************
 ***************************************************************************/

unsigned int GetCCycle(double time)
{
  int      i;                                              /* loop counter */
  double   *table;                         /* local copy of divtimes table */
  
/* assign 'table' to the appropriate division schedule */

  if ( olddivstyle ) {
    if ( defs.ndivs != 3 ) 
      error("GetCCycle: only 3 cell divisions allowed for oldstyle (-o)");
    table = (double *)old_divtimes;
  } else 
    if ( defs.ndivs == 0 )
      return 14;           /* if defs.ndivs == 0: we're always in cycle 14 */
    else if ( defs.ndivs == 1 ) 
      table = (double *)divtimes1;
    else if ( defs.ndivs == 2 ) 
      table = (double *)divtimes2;
    else if ( defs.ndivs == 3 ) 
      table = (double *)divtimes3;
    else if ( defs.ndivs == 4 ) 
      table = (double *)divtimes4;
    else
      error("GetCCycle: can't handle %d cell divisions!", defs.ndivs);
    
/* evaluate number of cell cycle for current time; note that for the exact *
 * time of cell division, we'll return the number of the previous cell cy- *
 * cyle                                                                    */

  for (i=0; i<defs.ndivs; i++)
    if ( time > table[i] )
      return 14-i;
  
  return 14-(i++);
}



/*** ParseLineage: takes lineage number as input and returns the cleavage **
 *                 cycle the nucleus belongs to.                           *
 ***************************************************************************/

unsigned int ParseLineage(unsigned int lin)
{
  if      ( lin & CYCLE14 )
    return 14;
  else if ( lin & CYCLE13 )
    return 13;
  else if ( lin & CYCLE12 )
    return 12;
  else if ( lin & CYCLE11 )
    return 11;
  else if ( lin & CYCLE10 )
    return 10;
  else
    error("ParseLineage: illegal lineage number %d", lin);
  
  return 0;                            /* just to make the compiler happy! */
} 



/*** GetDivtable: returns times of cell divisions depending on ndivs and ***
 *                olddivstyle; returns NULL in case of an error            *
 ***************************************************************************/

double *GetDivtable(void)
{
  if ( olddivstyle ) {
    if ( defs.ndivs != 3 ) 
      error("GetDivtable: only 3 cell divisions allowed for oldstyle (-o)");
    return (double *)old_divtimes;
  } else 
    if ( defs.ndivs == 0 )
      return NULL;                           /* return error if ndivs == 0 */
    else if ( defs.ndivs == 1 ) 
      return (double *)divtimes1;
    else if ( defs.ndivs == 2 ) 
      return (double *)divtimes2;
    else if ( defs.ndivs == 3 ) 
      return (double *)divtimes3;
    else if ( defs.ndivs == 4 ) 
      return (double *)divtimes4;
    else 
      error
	("GetDivtable: can't handle %d cell divisions!", defs.ndivs);

  return NULL;
}



/*** GetDurations: returns pointer to durations of cell divisions de- ******
 *                 pending on ndivs and olddivstyle; returns NULL in case  *
 *                 of an error                                             *
 ***************************************************************************/

double *GetDurations(void)
{
  if ( olddivstyle ) {
    if ( defs.ndivs != 3 ) 
      error("GetDurations: only 3 cell divisions allowed for oldstyle (-o)");
    return (double *)old_div_duration;
  } else 
    if ( defs.ndivs == 0 )
      return NULL;                           /* return error if ndivs == 0 */
    else if ( defs.ndivs == 1 ) 
      return (double *)div_duration1;
    else if ( defs.ndivs == 2 ) 
      return (double *)div_duration2;
    else if ( defs.ndivs == 3 ) 
      return (double *)div_duration3;
    else if ( defs.ndivs == 4 ) 
      return (double *)div_duration4;
    else
      error("GetDurations: can't handle %d cell divisions!", defs.ndivs);
  
  return NULL;
}



/*** GetGastTime: returns time of gastrulation depending on ndivs and ******
 *                olddivstyle; returns 0 in case of an error; if a custom  *
 *                gastrulation time is chosen with -S, it will be returned *
 *                only if it's bigger than the normal gastrulation time    *
 ***************************************************************************/

double GetGastTime(void)
{
  if ( olddivstyle ) { 
    if ( defs.ndivs != 3 ) 
      error("GetDurations: only 3 cell divisions allowed for oldstyle (-o)");

    if ( custom_gast > old_gast_time )
      return custom_gast;
    else
      return old_gast_time;

  } else 

    if ( defs.ndivs == 0 )

      if ( custom_gast > gast_time0 )
	return custom_gast;
      else
	return gast_time0;

    else if ( defs.ndivs == 1 ) 

      if ( custom_gast > gast_time1 )
	return custom_gast;
      else
	return gast_time1;

    else if ( defs.ndivs == 2 ) 

      if ( custom_gast > gast_time2 )
	return custom_gast;
      else
	return gast_time2;

    else if ( defs.ndivs == 3 ) 

      if ( custom_gast > gast_time3 )
	return custom_gast;
      else
	return gast_time3;

    else if ( defs.ndivs == 4 )

      if ( custom_gast > gast_time4 )
	return custom_gast;
      else
	return gast_time4;

    else 
      error("GetGastTime: can't handle %d cell divisions!", defs.ndivs);
  
  return 0;
}



/*** GetD: returns diffusion parameters D according to the diff. params. ***
 *         in the data file and the diffusion schedule used                *
 *   NOTE: Caller must allocate D_tab                                      *
 ***************************************************************************/

void GetD(double t, double *d, char diff_schedule, double *D_tab)
{
  static double *table;                          /* local copy of divtimes */

  int           i;                                         /* loop counter */
  double        gast;                                 /* gastrulation time */
  double        cutoff;                                     /* cutoff time */
  double        lscale = 1;                              /* scaling factor */


/* first time GetD is called: set pointer to the right division table */

  if ( !d_flag ) {
    if ( olddivstyle ) {
      if ( defs.ndivs != 3 ) 
	error("GetD: only 3 cell divisions allowed for oldstyle (-o)");
      table = (double *)old_divtimes;
    } else 
      if ( defs.ndivs == 0 )
	;                  /* no need for table, if there are no cell divs */
      else if ( defs.ndivs == 1 ) 
	table = (double *)divtimes1;
      else if ( defs.ndivs == 2 ) 
	table = (double *)divtimes2;
      else if ( defs.ndivs == 3 ) 
	table = (double *)divtimes3;
      else if ( defs.ndivs == 4 ) 
	table = (double *)divtimes4;
      else 
	error("GetD: can't handle %d cell divisions!", defs.ndivs);
    d_flag=1;
  }

/* this loop takes lscale square for each cell division, i.e. the earlier  *
 * we are the bigger lscale (and the smaller the D's that we return        */

  for(i=0; i<defs.ndivs; i++)  
    if (t<table[i])         
      lscale *= 2; 

/* diffusion schedule A: all Ds always the same */

  if(diff_schedule == 'A') 
    for(i=0; i<defs.ngenes; i++)
      D_tab[i] = d[0];

/* diffusion schedule B: all genes have different D's that depend on in-   *
 * verse l-square                                                          */

  else if (diff_schedule == 'B') 
    for(i=0; i<defs.ngenes; i++ )
      D_tab[i] = d[i] / (lscale * lscale);
  
/* diffusion schedule C: all genes have the same D that depends on inverse *
 * l-square                                                                */

  else if (diff_schedule == 'C') 
    for(i=0; i<defs.ngenes; i++)
      D_tab[i] = d[0] / (lscale * lscale);
  
/* diffusion schedule D: all genes have different D's which don't change   *
 * over time                                                               */

  else if (diff_schedule == 'D')  
    for(i=0; i<defs.ngenes; i++ ) 
      D_tab[i] = d[i];
                                    
/* diffusion schedule E: used cutoff at gast-12 otherwise just like B */

  else if (diff_schedule == 'E') {            
    if ( olddivstyle ) {
      if ( defs.ndivs != 3 ) 
	error("GetD: only 3 cell divisions allowed for oldstyle (-o)");
      gast=old_gast_time;
    } else 
      if ( defs.ndivs == 0 )
	gast=gast_time0;
      else if ( defs.ndivs == 1 ) 
	gast=gast_time1;
      else if ( defs.ndivs == 2 ) 
	gast=gast_time2;
      else if ( defs.ndivs == 3 ) 
	gast=gast_time3;
      else if ( defs.ndivs == 4 ) 
	gast=gast_time4;
      else
	error("GetD: can't handle %d cell divisions!", defs.ndivs);
    
    cutoff = gast - 12.0;      /* This value probably wrong; see Merrill88 */
    for(i=0; i<defs.ngenes; i++ )
      D_tab[i] = (t < cutoff) ?  d[i] / (lscale * lscale) : 0.;

    
/* any other diffusion schedule: error! */

  } else 
    error("GetD: no code for schedule %c!", diff_schedule); 

}



/*** Theta: Returns the value of theta(t) in the autonomous ***
 *            version of the equations                              *
 ***************************************************************************/

int Theta(double time)
{
  int i;

  static double *dt;                     /* pointer to division time table */
  static double *dd;                 /* pointer to division duration table */



  if ( !theta_flag ) {                                 /* only do this once */
    if ( olddivstyle ) {            /* get pointers to division time table */
      dt=(double *)old_divtimes;            /* and division duration table */
      dd=(double *)old_div_duration;
    } else {
      if ( defs.ndivs == 0 ) {
	dt=(double *)full_divtimes0;
	dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 1 ) {
	dt=(double *)full_divtimes1;
	dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 2 ) {
	dt=(double *)full_divtimes2;
	dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 3 ) {
	dt=(double *)full_divtimes3;
	dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 4 ) {
	dt=(double *)full_divtimes4;
	dd=(double *)full_div_durations;
      } else 
	error("Theta: can't handle %d cell divisions!", defs.ndivs);
    }
    theta_flag=1;
  }

/* checks if we're in a mitosis; we need the 10*DBL_EPSILON kludge for gcc *
 * on Linux which can't handle truncation errors very well                 */

  for(i=0; i<TOTAL_DIVS; i++)            
  {
/*  	printf("Thetatime=%.16f,[%.16f,%.16f]\n",time,(*(dt+i) - *(dd+i) + HALF_EPSILON), (*(dt+i) + HALF_EPSILON)); */
    if ((time <= (*(dt+i) + HALF_EPSILON)) 
	&& (time >= (*(dt+i) - *(dd+i) + HALF_EPSILON))) 
      return MITOSIS;                          
  }

  return INTERPHASE;
}

/*** GetRule: returns the appropriate rule for a given time; used by the ***
 *            derivative function                                          *
 ***************************************************************************/

int GetRule(double time)
{
  int i;

  static double *dt;                     /* pointer to division time table */
  static double *dd;                 /* pointer to division duration table */



/* the following block of code are implemented like this (instead of       */
/* calling GetDivtable or GetDurations) to save function calls and increa- */
/* se performance (GetRule is called from within the inner loop)           */


  if ( !rule_flag ) {                                 /* only do this once */
    if ( olddivstyle ) {            /* get pointers to division time table */
      dt=(double *)old_divtimes;            /* and division duration table */
      dd=(double *)old_div_duration;
    } else {
      if ( defs.ndivs == 0 ) {
	return INTERPHASE;                /* no cell division? no MITOSIS! */
      } else if ( defs.ndivs == 1 ) {
	dt=(double *)divtimes1;
	dd=(double *)div_duration1;
      } else if ( defs.ndivs == 2 ) {
	dt=(double *)divtimes2;
	dd=(double *)div_duration2;
      } else if ( defs.ndivs == 3 ) {
	dt=(double *)divtimes3;
	dd=(double *)div_duration3;
      } else if ( defs.ndivs == 4 ) {
	dt=(double *)divtimes4;
	dd=(double *)div_duration4;
      } else 
	error("GetRule: can't handle %d cell divisions!", defs.ndivs);
    }
    rule_flag=1;
  }

/* checks if we're in a mitosis; we need the 10*DBL_EPSILON kludge for gcc *
 * on Linux which can't handle truncation errors very well                 */

  for(i=0; i<defs.ndivs; i++)            
  {
/*  	printf("GetRuletime=%.16f,[%.16f,%.16f]\n",time,(*(dt+i) - *(dd+i) + HALF_EPSILON), (*(dt+i) + HALF_EPSILON)); */
    if ((time <= (*(dt+i) + HALF_EPSILON)) 
	&& (time >= (*(dt+i) - *(dd+i) + HALF_EPSILON))) 
      return MITOSIS;                          
  }

  return INTERPHASE;
}



/*** GetIndex: this functions returns the genotype index for a given *******
 *             genotype number for reading the GenoType struct.            *
 ***************************************************************************/

int GetIndex(char *genotype)
{
  int     i;                                         /* local loop counter */

  for (i=0; i < nalleles; i++)               /* nalleles static to score.c */
    if( !(strcmp(bcdtype[i].genotype, genotype)) )
      return i;

  error("GetIndex: could not find index for genotype %s", genotype);
  return -1000;
}



/*** MakeTable: this function constructs a time table for Blastoderm *******
 *             based on a comand line option (print_stepsize) and the      *
 *             cell division tables here in maternal.c                     *
 * DISCLAIMER: this thing is written in a very bad way; but hey it does    *
 *             its job!!!!!!!                                              *
 ***************************************************************************/

DArrPtr MakeTable(double p_stepsize)
{
  int        i;                                      /* local loop counter */
  int        t;                                            /* time counter */

  double     time;                                       /* double counter */
  double     gast;                                    /* gastrulation time */

  DArrPtr    table;                           /* time table to be returned */


  
  gast = GetGastTime();                           /* get gastrulation time */

  if ( p_stepsize > gast )
    error("GetTable: output stepsize can't be larger than gast time");

  t=0;
  for ( time=0; time<gast; time+=p_stepsize ) {
    t++;
  }

  table.size = t+1;                /* add one for gastrulation time itself */
  if ( !(table.array=(double *)calloc(table.size, sizeof(double))) )
    error("GetTable: error allocating times array");   
  
  time=0;
  for ( i=0; i<table.size-1; i++ ) {
    table.array[i] = time;
    time+=p_stepsize;
  }
  table.array[i] = gast;
  
  return table; 
}






/*** FUNCTIONS THAT READ DATA FROM FILE INTO STRUCTS (BIAS & BCD) *********/

/*** ReadTheProblem: reads the problem section of a data file into the *****
 *                   TheProblem struct.                                    *
 ***************************************************************************/

TheProblem ReadTheProblem(FILE *fp)
{
  TheProblem          l_defs;           /* local copy of TheProblem struct */

  fp = FindSection(fp, "problem");                 /* find problem section */
  if( !fp )
    error("ReadTheProblem: cannot locate problem section");

  fscanf(fp, "%*s\n");             /* advance pointer past first text line */

  if ( 1 != (fscanf(fp, "%d\n", &l_defs.ngenes)) )      /* read # of genes */
    error("ReadTheProblem: error reading problem section (ngenes)");
	     
  fscanf(fp, "%*s\n");    /* advance the pointer past the second text line */

  if ( 1 != (fscanf(fp, "%d\n", &l_defs.egenes)) )      /* read # of
  external inputs */
    error("ReadTheProblem: error reading problem section (egenes)");
	     
  fscanf(fp, "%*s\n");    /* advance the pointer past the second text line */

  l_defs.gene_ids = (char *)calloc(l_defs.ngenes+1, sizeof(char));
  if ( 1 != (fscanf(fp, "%s\n", l_defs.gene_ids)) ) /* read geneID string */
    error("ReadTheProblem: error reading problem section (gene_ids)");

  fscanf(fp, "%*s\n");    /* advance the pointer past the second text line */

  if (l_defs.egenes > 0 ) {
	  
	  l_defs.egene_ids = (char *)calloc(l_defs.egenes+1, sizeof(char));
	  if ( 1 != (fscanf(fp, "%s\n", l_defs.egene_ids)) ) /* read geneID string */
		error("ReadTheProblem: error reading problem section (egene_ids)");
  } else {
	  l_defs.egene_ids = (char *)calloc(1, sizeof(char));
	  l_defs.egene_ids = "\0";
	  fscanf(fp, "%*s\n");                       /* next line (ignore empty external gene ID line) */
  }	  
  
  fscanf(fp, "%*s\n");                       /* next line (ignore comment) */

  if ( 1 != (fscanf(fp, "%d\n", &l_defs.ndivs)) )   /* read # of cell divs */
    error("ReadTheProblem: error reading problem section (ndivs)");

  fscanf(fp, "%*s\n");     /* advance the pointer past the third text line */
  if ( 1 != (fscanf(fp, "%d\n", &l_defs.nnucs)) )/* read the max # of nucs */
    error("ReadTheProblem: error reading problem section (nnucs)");

  fscanf(fp, "%*s\n");                    /* advance the pointer once more */

  if ( 1 != (fscanf(fp, "%c\n", &l_defs.diff_schedule)))/* read diff sched */
    error("ReadTheProblem: error reading problem section (diff. schedule)");

  return l_defs;
}



/*** ReadGenotypes: This function reads all the genotypes in a datafile & **
 *                  returns an SList with genotype number and pointers to  *
 *                  the corresponding section titles for bias, bcd & facts *
 ***************************************************************************/

Slist *ReadGenotypes(FILE *fp)
{
                                      /* Buffer strings for reading:       */
  char      biasbuf[MAX_RECORD];      /* section title of bias section     */
  char      factsbuf[MAX_RECORD];     /* section title of data section     */
  char      matbuf[MAX_RECORD];       /* section title of bcd section      */
  char      histbuf[MAX_RECORD];       /* section title of history section      */
  char      extbuf[MAX_RECORD];       /* section title of external
  input section      */
  char      gtbuf[MAX_RECORD];        /* genotype string                   */

  char      *record;                  /* pointer to current data record    */

  Slist     *current;                 /* holds current element of Slist    */
  Slist     *first;                   /* pointer to first element of Slist */
  Slist     *last = NULL;             /* pointer to last element of Slist  */


/*** open the data file and locate genotype section ************************/
  
  fp = FindSection(fp, "genotypes");
  if( !fp )
    error("ReadGenotypes: cannot locate genotypes");
  
  if ( !(record=(char *)calloc(MAX_RECORD, sizeof(char))) )
     error("ReadGenotypes: error allocating record");
  first   = init_Slist();
  current = first;
  
/*** read titles of bias, data and bcd sections and genotype number for    *
 *   each genotype *********************************************************/
  
  while ( strncmp((record=fgets(record, MAX_RECORD, fp)), "$$", 2)) { 
    
    if ( 6 != sscanf(record, "%s %s %s %s %s %s", 
		     biasbuf, factsbuf, matbuf, histbuf, extbuf, gtbuf) )
      error("ReadGenotypes: error reading %s", record);

/* we eventually want to get rid of the hard wired genotype identifiers    *
 * and replace them by some kind of mechanism by which we can include any  *
 * genes we want in a scoring or annealing run. Think about this!          */

    if ( strlen(gtbuf) != defs.ngenes ) 
      error("ReadGenotypes: bad genotype string %s (does not match ngenes)",
	    gtbuf);

    if ( !current )
      current = init_Slist();
 
    if (!(current->bias_section = (char *)calloc(MAX_RECORD, sizeof(char))))
      error("ReadGenotypes: error allocating bias_section");
    current->bias_section = strcpy(current->bias_section, biasbuf);
  
    if (!(current->fact_section = (char *)calloc(MAX_RECORD, sizeof(char))))
      error("ReadGenotypes: error allocating fact_section");
    current->fact_section = strcpy(current->fact_section, factsbuf);
  
    if (!(current->bcd_section  = (char *)calloc(MAX_RECORD, sizeof(char))))
      error("ReadGenotypes: error allocating bcd_section");
    current->bcd_section  = strcpy(current->bcd_section, matbuf);
  
    if (!(current->hist_section  = (char *)calloc(MAX_RECORD, sizeof(char))))
      error("ReadGenotypes: error allocating hist_section");
    current->hist_section  = strcpy(current->hist_section, histbuf);
  
    if (!(current->ext_section  = (char *)calloc(MAX_RECORD, sizeof(char))))
      error("ReadGenotypes: error allocating ext_section");
    current->ext_section  = strcpy(current->ext_section, extbuf);
  
    if (!(current->genotype     = (char *)calloc(MAX_RECORD, sizeof(char))))
      error("ReadGenotypes: error allocating genotype string");
    current->genotype     = strcpy(current->genotype, gtbuf);

    addto_Slist(last, current);

    last = current;
    current = NULL;
  }

  free(record);

  return first;
}



/*** ReadBicoid: reads the bcd section of a data file into a linked list; **
 *               also determines maxconc from the bicoid gradient          *
 ***************************************************************************/

Blist *ReadBicoid(FILE *fp, char *section)
{
  int       c;                           /* holds char for parser          */
  int       lead_punct;                  /* flag for leading punctuation   */
  
  double    maxv = -1.;                  /* maximum v (protein conc.)      */

  char      *base;                       /* pointer to beginning of string */
  char      *record;                     /* pointer to string, used as     */
                                         /* counter                        */
  Blist     *current;                    /* holds current element of Blist */
  Blist     *inlist;                     /* holds whole read Blist         */
   
  if ( (fp=FindSection(fp, section)) ) {                /* position the fp */
   
    base    = (char *)calloc(MAX_RECORD, sizeof(char));
    current = NULL;
    inlist  = NULL;

/* while loop: reads and processes lines from file until sections ends *****/

    while ( strncmp(( base=fgets(base, MAX_RECORD, fp)), "$$", 2 )) {

      record     = base;                    /* base always points to start */
      lead_punct = 0;                                         /* of string */

/* for loop: parses and processes each line from the data file *************/

      c=(int)*record;
      while ( c != '\0' ) {
      
	if( isdigit(c) ) {                            /* number means data */
	
	  record = base;                  /* reset pointer to start of str */
	  current = init_Blist();	  /* allocate memory for lnkd list */

	  if ( 2 != sscanf(record, "%d %lg", 
			   &(current->lineage), &(current->conc)) )
	    error("ReadBicoid: error reading %s", base);
	  
	  if ( current->conc > maxv )
	    maxv = current->conc;

	  inlist = addto_Blist(inlist, current);
	  break;

	}
	else if ( isalpha(c) ) {                   /* letter means comment */
	  break;
	}
	else if ( c == '-' ) {                /* next two elsifs for punct */
	  if ( ((int) *(record+1)) == '.' )
	    record++;
	  lead_punct = 1;
	  c=(int)*(++record);
	}
	else if ( c == '.' ) {
	  lead_punct = 1;
	  c=(int)*(++record);
	}
	else if ( ispunct(c) ) {              /* other punct means comment */
	  break;
	}
	else if ( isspace(c) ) {             /* ignore leading white space */
	  if ( lead_punct )               /* white space after punct means */
	    break;                                              /* comment */
	  else {
	    c=(int)*(++record);            /* get next character in record */
	  }
	}	
	else {
	  error("ReadBicoid: illegal character in %s", base);
	}
      }
    }
    
    if ( maxv > 12. )                        /* oldstyle or newstyle data? */
      maxconc = 255.;
    else
      maxconc = 12.;

    free(base);
    return inlist;

  } else {
    
    return NULL;
  }
}



/*** ReadData: reads in a data or bias section and puts it in a linked *****
 *             list of arrays, one line per array; ndp is used to count    *
 *             the number of data points in a data file (ndp), which is    *
 *             used to calculate the root mean square (RMS) if required    *
 *                                                                         * 
 *             ReadData allows comments that start with a letter or punc-  *
 *             tuation mark, and treats lines starting with a number or    *
 *             .num or -.num as data. It expects to read an int and        *
 *             ngenes + 1 data points: lineage, time and ngenes protein    *
 *             concentrations. It expects data to be in increasing spatial *
 *             and temporal order.                                         *
 ***************************************************************************/

Dlist *ReadData(FILE *fp , char *section, int *ndp)
{
  int     c;                           /* holds char for parser            */
  int     lead_punct;                  /* flag for leading punctuation     */
  int     i;                           /* loop counter                     */
  
  char    *base;                       /* pointer to beginning of string   */
  char    *record;                     /* pointer to string, used as       */
                                       /* counter                          */
  Dlist   *current;                    /* holds current element of Dlist   */
  Dlist   *inlist;                     /* holds whole read Dlist           */
   
/* the following chars are used for parsing lines of data (record)         */
/* using sscanf                                                            */

  char       *fmt  = NULL;             /* used as format for sscanf        */
  char       *skip = NULL;             /* skip this stuff!                 */

  const char init_fmt[] = "%*d ";      /* initial format string            */
  const char skip_fmt[] = "%*lg ";     /* skip one more float              */
  const char read_fmt[] = "%lg ";      /* skip one more float              */
  
  if ( (fp=FindSection(fp,section)) ) {                 /* position the fp */
  
    base    = (char *)calloc(MAX_RECORD, sizeof(char));
    fmt     = (char *)calloc(MAX_RECORD, sizeof(char));
    skip    = (char *)calloc(MAX_RECORD, sizeof(char));
    
    current = NULL;
    inlist  = NULL;
    
/* while loop: reads and processes lines from file until sections ends *****/
    
    while ( strncmp(( base=fgets(base, MAX_RECORD, fp)), "$$", 2)) {
      
      record     = base;                  /* base always points to start   */
      lead_punct = 0;                     /* of string                     */
      
/* while loop: parses and processes each line from the data file *************/
      
      c=(int)*record;
      
      while ( c != '\0' ) {
	
	if( isdigit(c) ) {                            /* number means data */
	  
	  record = base;                  /* reset pointer to start of str */
	  current = init_Dlist(defs.ngenes+1);
	  
	  if ( 1 != sscanf(record, "%d ", &(current->lineage)) )
	    error("ReadData: error reading %s", base);
	  
/* the following loop reads a line of data of variable length into a d-    */
/* array, skipping lineage number and all previously read data values      */
	  
	  skip = strcpy(skip, init_fmt); /* format to be skipped by sscanf */
	  
	  for(i=0; i < (defs.ngenes + 1); i++) {
	    
	    fmt  = strcpy(fmt, skip);     /* all this stuff here is to get */
	    fmt  = strcat(fmt, read_fmt);  /* the fmt string for sscanf to */
	    skip = strcat(skip, skip_fmt);           /* the correct format */
	    
	    if ( 1 != sscanf(record, (const char *)fmt, &(current->d[i])) )
	      error("ReadData: error reading %s", base);

                                           /* update number of data points */
	    if ( ( i != 0) && (current->d[i] != IGNORE) )
	      (*ndp)++;
	    
	  }
	                                  /* now add this to the lnkd list */
	  inlist = addto_Dlist(inlist, current);
	  break;
	}                        
	
	else if ( isalpha(c) ){                    /* letter means comment */
	  break;
	}
	else if ( c == '-' ) {                /* next two elsifs for punct */
	  if ( ((int)*(record+1)) == '.')
	    record++;
	  lead_punct = 1;
	  c=(int)*(++record);
	}
	else if ( c == '.' ) {
	  lead_punct = 1;
	  c=(int)*(++record);
	}
	else if ( ispunct(c) ){               /* other punct means comment */
	  break;
	}
	else if ( isspace(c) ) {             /* ignore leading white space */
	  if ( lead_punct )               /* white space after punct means */
	    break;                                              /* comment */
	  else {                          
	    c=(int)*(++record);            /* get next character in record */
	  }      
	}                                 
	else {
	  error ("ReadData: Illegal character in %s", base);
	}
      }
    }
    
    free(base);
    free(fmt);
    free(skip);
    
    return inlist;

  } else {
    
    return NULL;

  }
}



/*** ReadTimes: reads a time table from a file and returns a DArrPtr *******
 * FILE FORMAT: one time per line separated with newlines                  *
 *        NOTE: max. times tab size is 10000                               *
 ***************************************************************************/

DArrPtr ReadTimes(char *timefile)
{
  FILE         *fp;

  int          count = 0;
  double       gast;

  DArrPtr      tabtimes;
  double       *tabptr;

  gast = GetGastTime();

  tabtimes.size  = 10000;
  tabtimes.array = (double *)calloc(10000, sizeof(double));
  tabptr = tabtimes.array;

  fp = fopen(timefile, "r");
  if ( !fp )
    file_error("ReadTimes");

  if ( 1 != (fscanf(fp, "%lg\n", tabptr)) )
    error("ReadTimes: time file %s empty!", timefile);
  tabptr++;
  count++;

  while ( (fscanf(fp, "%lg\n", tabptr)) != EOF ) {
    if ( (*tabptr < 0) || (*tabptr<=*(tabptr-1)) || (*tabptr > gast) )
      error("ReadTimes: invalid time(s) in %s!", timefile);
    tabptr++;
    count++;
  }
      
  fclose(fp);

  tabtimes.size = count;
  tabtimes.array = (double *)realloc(tabtimes.array, count*sizeof(double));

  return tabtimes;
}



/*** ReadGuts: reads the $gutsdefs section in a data file into an array ****
 *             of strings which then need to get parsed                    *
 ***************************************************************************/

char **ReadGuts(FILE *fp)
{
  char      **gutsbuf = NULL;                /* buffer strings for reading */
  char       *record  = NULL;  
  int	      i,j;                                        /* loop counters */

  if ( !(gutsbuf = (char **)calloc(MAX_RECORD, sizeof(char *)) ) ) 
    error("ReadGuts: error allocating memory for gutsdefs");

  if ( !(record  = (char *)calloc(MAX_RECORD, sizeof(char))) )
    error("ReadGuts: error allocating record"); 
  
/*** locate gutsdefs section ***********************************************/
  
  fp = FindSection(fp, "gutsdefs");
  if( !fp )
    error("ReadGuts: cannot locate gutsdefs");
			
/* Remember the location of the first line and count the number of lines **/

  i = 0; 
  while ( strncmp( (record=fgets(record, MAX_RECORD, fp)), "$$", 2) ) { 
    if (! (*(gutsbuf+i) = (char *)calloc(MAX_RECORD, sizeof(char))))
      error("ReadGuts: error allocating memory for gutsdefs strings");   
    *(gutsbuf+i) = strcpy(*(gutsbuf+i), record);
    i++;
  }

  if (!(*gutsbuf)) 
    error("ReadGuts: gutsdefs section is empty!");

  free(record);

  return gutsbuf;
}



/*** List2Bicoid: takes a Blist and returns the corresponding BArrPtr ******
 *                structure; also initializes the lin_start array which    *
 *                contains the lineage number of the most anterior nucleus *
 *                for each cell cycle (we need this for cell division and  *
 *                printing model output)                                   *
 ***************************************************************************/

BArrPtr List2Bicoid(Blist *inlist)
{
  int          i = 0;                                /* local loop counter */
  int          n;                          /* used to evaluate # of nuclei */
  int          lin_count = 0;        /* counter for initializing lin_start */

  unsigned int ccycle;                   /* what cleavage cycle are we in? */
  unsigned int samecycle = 0;            /* same cleavage cycle as before? */

  Blist        *current;                       /* current element of Blist */
  Blist        *start;                     /* used to evaluate # of nuclei */

  BArrPtr      bicoid;                           /* BArrPtr to be returned */

  bicoid.size  = 0;
  bicoid.array = NULL;

  if ( !(lin_start = (int *)calloc(defs.ndivs+1, sizeof(int))) )
    error("List2Bicoid: could not allocate lin_start array");
  
/*** for loop: step through linked list and copies values into an array    *
 *             of BcdGrads; there's one BcdGrad for each cleavage cycle;   *
 *             each BcdGrad has: - ccycle (the cleavage cycle number)      *
 *                               - gradient.size (# of nuclei for grad.)   *
 *                               - gradient.array (pointer to array)       *
 ***************************************************************************/

  for ( current=inlist; current; current=current->next ) {
	
    ccycle = ParseLineage( current->lineage );       /* which cycle is it? */
    
/* allocate new gradient struct for new cell cycle */

    if ( ccycle != samecycle ) { 

      samecycle = ccycle;
      
      bicoid.size++;
      bicoid.array = 
	(BcdGrad *)realloc(bicoid.array, bicoid.size*sizeof(BcdGrad));
      
/* this loop determines the number of nuclei in each cycle */

      n = 0;             
      for ( start=current; ParseLineage(current->lineage) == samecycle;
	    current=current->next ) {
	n++;
	if ( !(current->next) )    /* don't count garbage -> core dumps... */
	  break;                                        
      }
      current = start;                   /* reset pointer to where we were */
	
/* initialize lin_start: this array is later used by Blastoderm and such   */

      lin_start[defs.ndivs-lin_count] = current->lineage;
      lin_count++;

/* allocate array for gradient here and reset the counter */

      bicoid.array[bicoid.size-1].ccycle = samecycle;        /* next three */
                           /* lines define BcdGrad for each cleavage cycle */
      bicoid.array[bicoid.size-1].gradient.array = 
	(double *)calloc(n, sizeof(double));      
      bicoid.array[bicoid.size-1].gradient.size  = n;
      i = 0;

    }

/* in any case: read concentration into gradient array */

    bicoid.array[bicoid.size-1].gradient.array[i] = current->conc;
    i++;
  }

  return bicoid;
}  
  
    

/*** List2Bias: takes a Dlist and returns the corresponding DArrPtr ********
 *              structure.                                                 *
 ***************************************************************************/
       
NArrPtr List2Bias(Dlist *inlist)
{
  int         i = 0;
  int         j;                                    /* local loop counters */

  int         n;                           /* used to evaluate # of nuclei */
  double      now = -999999999.0;                     /* variable for time */

  Dlist       *current;                  /* holds current element of Dlist */
  Dlist       *start;                  /* pointer used for evaluating # of */

  NArrPtr     bias;                              /* DArrPtr to be returned */
  
  bias.size  = 0;
  bias.array = NULL;

/*** for loop: steps through linked list and copies values into an array   *
 *             of NucStates. There's one NucState for each time step       *
 *             Each NucState has: - time (the time)                        *
 *                                - state.size (# of genes * # of nucs)    *
 *                                - state.array (pointer to array)         *
 ***************************************************************************/
    
  for ( current=inlist; current; current=current->next) {
       
    if ( current->d[0] != now ) {            /* a new time point: allocate */
      now = current->d[0];                             /* the time is now! */
      bias.size++;                      /* one NucState for each time step */
      bias.array =                       /* allocate for one more NucState */
	(NucState *)realloc(bias.array, bias.size * sizeof(NucState));

/* determine number of nuclei per time step */

      n = 0;               
      for ( start=current; current->d[0] == now; current=current->next ) {
	n++;
	if ( !(current->next) )      /* don't count garbage and cause dis- */
	  break;                                  /* may and core dumps... */
      }
      current = start;                  /* reset list ptr to where we were */
	
/* allocate a bias array for each biastime */

      bias.array[bias.size-1].time = now;                                          
      bias.array[bias.size-1].state.array =  
	(double *)calloc(n*defs.ngenes, sizeof(double)); 
      bias.array[bias.size-1].state.size = n*defs.ngenes;  

      i=0;
    }
      
/* always: read concs into array */

    for( j=1; j <= defs.ngenes; j++ ) {
      bias.array[bias.size - 1].state.array[i]
	= current->d[j];
      i++;
    }
  }

  return bias;
}


DataTable *List2Interp(Dlist *inlist, int num_genes)
{ 


  int       i = 0;
  int       j;                                      /* local loop counters */

  double    now = -999999999.;            /* assigns data to specific time */
  
  Dlist     *current;                    /* holds current element of Dlist */
  
  DataTable *D;                                 /* local copy of DataTable */

  if (!full_lin_start) {
  	
	if ( !(full_lin_start = (int *)calloc(defs.ndivs+1, 
												sizeof(int))) )
    	error("List2Interp: could not allocate full_lin_start array");

 	full_ccycles =  defs.ndivs + 1;
	
	for (i=0; i < full_ccycles; i++)
		full_lin_start[i] = lin_start[i];
  }		

  D = (DataTable *)malloc(sizeof(DataTable));
                                         /* Initialize DataTable structure */
  D->size = 0;
  D->record = NULL;
      
/*** for loop: steps through linked list and transfers facts into Data-    *
 *             Records, one for each time step                             *
 ***************************************************************************/

  for (current=inlist; current; current=current->next) {              

    if ( current->d[0] != now ) {             /* a new time point: allocate */
      now = current->d[0];                             /* the time is now! */
      D->size++;                           /* one DataRecord for each time */
      D->record =                                   /* allocate DataRecord */
	(DataRecord *)realloc(D->record,D->size*sizeof(DataRecord));
      
      D->record[D->size-1].time = now;          /* next three lines define */
      D->record[D->size-1].size = 0;                /* DataRecord for each */
      D->record[D->size-1].array = NULL;                      /* time step */
      i = 0;
    }
	
    for(j=1; j <= num_genes; j++) {     /* always: read concs into array */
	D->record[D->size-1].size++;            /* one more in this record */
	D->record[D->size-1].array =       /* reallocate memory for array! */
	  realloc(D->record[D->size-1].array, 
		  D->record[D->size-1].size * sizeof(DataPoint));
	    
/* the following two lines assign concentration value and index ************/
	    
	D->record[D->size-1].array[D->record[D->size-1].size-1].conc = 
	  current->d[j];
	D->record[D->size-1].array[D->record[D->size-1].size-1].index = i;

      i++;  

    }

/* initialize lin_start: this array is later used by Blastoderm and such   */

	if (ParseLineage(full_lin_start[full_ccycles - 1]) !=
			ParseLineage(current->lineage)  ) {
	  	full_ccycles++;

  		if ( !(full_lin_start = (int *)realloc(full_lin_start, 
											full_ccycles*sizeof(int))) )
		    error("List2Interp: could not allocate full_lin_start array");

		full_lin_start[full_ccycles-1] = current->lineage;

	}

  }

	qsort((void *) full_lin_start, full_ccycles, sizeof(int),(int (*)
	(void*,void*)) descend);

/*  for (i=0; i < full_ccycles; i++)
  	printf("History lineages before removing dups %d\n", full_lin_start[i]);*/
  
/* Now lets remove duplicates */

	i = 0;

	while (i < full_ccycles-1 )
	{
		if (full_lin_start[i] == full_lin_start[i+1])
		{
			memmove((full_lin_start+i),(full_lin_start+i+1),
								(full_ccycles-i-1)*sizeof(int));
/*			printf("Shifted %d elements to %d\n",full_ccycles-i-1,i);*/
			full_ccycles--;
			i--;
		} 

	i++;
	full_lin_start = (int *) realloc(full_lin_start,
									full_ccycles*sizeof(int));
	}

  return D;
}

int descend(int *x, int *y)
{

	if (*x < *y) return 1;
	else if (*x > *y) return -1;
	else return 0;

}



/* A FUNCTION WHICH IS NEEDED BY ALL OTHER READING FUNCS *******************/

/*** FindSection: This function finds a given section of the input file & **
 *                returns a pointer positioned to the first record of that *
 *                section. Section titles should be passed without the pre-*
 *                ceding '$'. If it can't find the right section, the func-*
 *                tion returns NULL.                                       *
 ***************************************************************************/

FILE *FindSection(FILE *fp, char *input_section)
{
  int      c;                      /* input happens character by character */
  int      nsought;                       /* holds length of section title */
  char     *base;                              /* string for section title */
  long     looksite;                            /* file position indicator */
  
  rewind(fp);                /* start looking at the beginning of the file */
  
  nsought = strlen(input_section);
  base = (char *)calloc(MAX_RECORD, sizeof(char));
  
/*** while loop goes through the file character by character... ************/
   while ( (c=getc(fp)) != EOF) {               /* ...until if finds a '$' */
    
    if ( c == '$') {                  /* found a sectioning control string */
      looksite = ftell(fp);                               /* where are we? */
      base = fgets(base, MAX_RECORD, fp);   /* get sect title (without $)*/
      if ( !(strncmp(base, input_section, nsought)) ) {
        fseek(fp, looksite, 0);        /* found the sought string: reposi- */
        fscanf(fp, "%*s\n");              /* tion the pointer to '$', then */
	free(base);
        return(fp);                    /* advance the pointer to after the */
      }                                     /* section title and return it */
      else {                          /* didn't find it: skip this control */
        fseek(fp, looksite, 0);                    /* record, keep looking */
        fscanf(fp, "%*s");                 /* NOTE: "%*s" advances pointer */
      }                                              /* without assignment */
    }
  }
   
  free(base);
 
  return(NULL);                         /* couldn't find the right section */
}



/*** KillSection: erases the section with 'title' from the file 'fp' *******
 ***************************************************************************/

void KillSection(char *filename, char *title)
{
  size_t length;                                 /* length of title string */

  char   *fulltitle;                                      /* title incl. $ */
  char   *temp;                                     /* temporary file name */
  char   *record;                         /* record to be read and written */
  char   *record_ptr;        /* pointer used to remember record for 'free' */

  char   *shell_cmd;                               /* used by system below */

  FILE   *fp;             /* name of file where section needs to be killed */
  FILE   *tmpfile;                               /* name of temporary file */
    

  fulltitle = (char *)calloc(MAX_RECORD, sizeof(char));
  temp      = (char *)calloc(MAX_RECORD, sizeof(char));
  record    = (char *)calloc(MAX_RECORD, sizeof(char));
  shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));

  record_ptr = record;            /* this is to remember record for 'free' */

  fp = fopen(filename, "r");                      /* open file for reading */
  if ( !fp )                           
    error("KillSection: error opening file %s", filename);

  temp = strcpy(temp,"killXXXXXX");               /* required by mkstemp() */
  if ( mkstemp(temp) == -1 )              /* get unique name for temp file */
    error("KillSection: error creating temporary file");

  tmpfile = fopen(temp, "w");               /* ... and open it for writing */
  if ( !tmpfile )
    error("KillSection: error opening temporary file");

  if ( !FindSection(fp, title) )
    error("KillSection: section to be killed not found");
  rewind(fp);
  
/* remove section by simply ignoring it while copying the file */

  fulltitle = strcpy(fulltitle, "$");
  fulltitle = strcat(fulltitle, title);
  length    = strlen(fulltitle);

  while (strncmp((record=fgets(record, MAX_RECORD, fp)), 
		 fulltitle, length))
    fputs(record, tmpfile);

  while (strncmp((record=fgets(record, MAX_RECORD, fp)), "$$", 2))
    ;

  do {
    record = fgets(record, MAX_RECORD, fp);
    if ( !record ) break;
  } while ( strncmp(record, "$", 1) );

  if ( record )
    fputs(record, tmpfile);

  while ( (record=fgets(record, MAX_RECORD, fp)) )
    fputs(record, tmpfile);

  fclose(fp);
  fclose(tmpfile);
  
/* rename tmpfile into new file */

  sprintf(shell_cmd, "cp -f %s %s", temp, filename);

  if ( -1 == system(shell_cmd) )
    error("KillSection: error renaming temp file %s", temp);

  if ( remove(temp) )
    warning("KillSection: temp file %s could not be deleted", temp);


/* clean up */

  free(record_ptr);
  free(temp);
  free(shell_cmd);
  free(fulltitle);
}





/* FOLLOWING FUNCTIONS ARE UTILITY FUNCTIONS FOR DIFFERENT LINKED LISTS ****
 * which are used to read in data of unknown size. All utility functions   *
 * follow the same scheme (X stands for the different letters below):      *
 *                                                                         *
 * - init_Xlist:         allocates first element and returns pointer to it *
 * - addto_Xlist:        adds the adduct to an already existing linkd list *
 * - free_Xlist:         frees memory of linked list                       *
 * - count_Xlist:        counts elements in the linked list                *
 *                                                                         *
 ***************************************************************************/

/*** Utility functions for Blist *******************************************/

Blist *init_Blist(void)
{
  Blist     *p;
  
  if( (p =(Blist *)malloc(sizeof(Blist))) ) {

    p->lineage = 0;
    p->conc    = 0;
    p->next    = NULL;

  }
  else
    error("init_Blist: couldn't allocate!");
	  
  return p;
}



Blist *addto_Blist(Blist *start, Blist *adduct)
{
  Blist     *current;

  if( !start )
    return adduct;
  
  current = start;

  while (current->next) {
    current = current->next;
  }

  current->next = adduct;
  return start;
}



void free_Blist(Blist *start)
{
  if(start->next)
    free_Blist(start->next);
  
  free(start);
}



int count_Blist(Blist *start)
{
  int        n = 0;

  while( start != NULL ) {
    n++;
    start = start->next;
  }

  return n;
}



/*** Utility functions for Dlist *******************************************/

Dlist *init_Dlist(int size)  
{
  Dlist    *p;                               

  if( (p=(Dlist *)malloc(sizeof(Dlist))) ) {
    
    if ( !(p->d=(double *)calloc(size, sizeof(double))) )
      error("initDlist: couldn't allocate!");    

    p->lineage = 0;
    p->next    = NULL;
  
}
  else 
    error("initDlist: couldn't allocate!");

  return p;
}



Dlist *addto_Dlist(Dlist *start, Dlist *adduct)
{
  Dlist         *current;

  if(!start)
    return adduct;

  current = start;

  while (current->next) {
    current = current->next;
  }

  current->next = adduct;
  return start;
}


void free_Dlist(Dlist *start)
{
  if(start->next)
    free_Dlist(start->next);
  
  free(start->d);
  free(start);
}



int count_Dlist(Dlist *start)
{
  int n = 0;
  while(start != NULL) {
    n++;
    start = start->next;
  }

  return n;
}



/*** Utility functions for Slist *******************************************/

Slist *init_Slist(void)
{
  Slist     *p;

  if( (p=(Slist *)malloc(sizeof(Slist))) ) {

    p->bias_section =       NULL;
    p->fact_section =       NULL;
    p->bcd_section  =       NULL;
    p->hist_section  =       NULL;
    p->ext_section  =       NULL;
    p->genotype     =       NULL;
    p->next         =       NULL;

  }
  else 
    error("initSlist: couldn't allocate");

  return p;
}



Slist *addto_Slist(Slist *start, Slist *adduct)
{
  Slist      *current;

  if(!start)
    return adduct;

  current = start;

  while (current->next) {
    current = current->next;
  }

  current->next = adduct;
  return start;
}



void free_Slist(Slist *start)
{
  if(start->next)
    free_Slist(start->next);

  free(start->bias_section);
  free(start->fact_section);
  free(start->bcd_section);
  free(start->hist_section);
  free(start->ext_section);
  free(start->genotype);
  free(start);
}



int count_Slist(Slist *start)
{
  int n = 0;
  while(start != NULL) {
    n++;
    start = start->next;
  }

  return n;
}

double *Get_Theta_Discons(int *theta_discon_size)
{

	int s,i,j;
	double *dt,*dd,*disc_array;
	
	s = (TOTAL_DIVS - defs.ndivs);
	disc_array = (double *) calloc(2*s, sizeof(double));

    if ( defs.ndivs == 0 ) {
		dt=(double *)full_divtimes0;
		dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 1 ) {
		dt=(double *)full_divtimes1;
		dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 2 ) {
		dt=(double *)full_divtimes2;
		dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 3 ) {
		dt=(double *)full_divtimes3;
		dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 4 ) {
		dt=(double *)full_divtimes4;
		dd=(double *)full_div_durations;
      } else 
		error("Get_Theta_Discons: can't handle %d cell divisions!", defs.ndivs);
	
	for (i=0,j=0;i<s;i++,j+=2) {	
		*(disc_array+j) = dt[defs.ndivs+i]; 
		*(disc_array+j+1) = dt[defs.ndivs+i] -
		dd[defs.ndivs+i];
	} 

	*theta_discon_size = 2*s;
	return disc_array;
}

Dlist *ReadInterpData(FILE *fp , char *section, int num_genes, int *ndp)
{
  int     c;                           /* holds char for parser            */
  int     lead_punct;                  /* flag for leading punctuation     */
  int     i;                           /* loop counter                     */
  
  char    *base;                       /* pointer to beginning of string   */
  char    *record;                     /* pointer to string, used as       */
                                       /* counter                          */
  Dlist   *current;                    /* holds current element of Dlist   */
  Dlist   *inlist;                     /* holds whole read Dlist           */
   
/* the following chars are used for parsing lines of data (record)         */
/* using sscanf                                                            */

  char       *fmt  = NULL;             /* used as format for sscanf        */
  char       *skip = NULL;             /* skip this stuff!                 */

  const char init_fmt[] = "%*d ";      /* initial format string            */
  const char skip_fmt[] = "%*lg ";     /* skip one more float              */
  const char read_fmt[] = "%lg ";      /* skip one more float              */
  
  if ( (fp=FindSection(fp,section)) ) {                 /* position the fp */
  
    base    = (char *)calloc(MAX_RECORD, sizeof(char));
    fmt     = (char *)calloc(MAX_RECORD, sizeof(char));
    skip    = (char *)calloc(MAX_RECORD, sizeof(char));
    
    current = NULL;
    inlist  = NULL;
    
/* while loop: reads and processes lines from file until sections ends *****/
    
    while ( strncmp(( base=fgets(base, MAX_RECORD, fp)), "$$", 2)) {
      
      record     = base;                  /* base always points to start   */
      lead_punct = 0;                     /* of string                     */
      
/* while loop: parses and processes each line from the data file *************/
      
      c=(int)*record;
      
      while ( c != '\0' ) {
	
	if( isdigit(c) ) {                            /* number means data */
	  
	  record = base;                  /* reset pointer to start of str */
	  current = init_Dlist(num_genes+1);
	  
	  if ( 1 != sscanf(record, "%d ", &(current->lineage)) )
	    error("ReadInterpData: error reading %s", base);
	  
/* the following loop reads a line of data of variable length into a d-    */
/* array, skipping lineage number and all previously read data values      */
	  
	  skip = strcpy(skip, init_fmt); /* format to be skipped by sscanf */
	  
	  for(i=0; i < (num_genes + 1); i++) {
	    
	    fmt  = strcpy(fmt, skip);     /* all this stuff here is to get */
	    fmt  = strcat(fmt, read_fmt);  /* the fmt string for sscanf to */
	    skip = strcat(skip, skip_fmt);           /* the correct format */
	    
	    if ( 1 != sscanf(record, (const char *)fmt, &(current->d[i])) )
	      error("ReadInterpData: error reading %s", base);

                                           /* update number of data points */
	    if ( ( i != 0) && (current->d[i] != IGNORE) )
	      (*ndp)++;
	    
	  }
	                                  /* now add this to the lnkd list */
	  inlist = addto_Dlist(inlist, current);
	  break;
	}                        
	
	else if ( isalpha(c) ){                    /* letter means comment */
	  break;
	}
	else if ( c == '-' ) {                /* next two elsifs for punct */
	  if ( ((int)*(record+1)) == '.')
	    record++;
	  lead_punct = 1;
	  c=(int)*(++record);
	}
	else if ( c == '.' ) {
	  lead_punct = 1;
	  c=(int)*(++record);
	}
	else if ( ispunct(c) ){               /* other punct means comment */
	  break;
	}
	else if ( isspace(c) ) {             /* ignore leading white space */
	  if ( lead_punct )               /* white space after punct means */
	    break;                                              /* comment */
	  else {                          
	    c=(int)*(++record);            /* get next character in record */
	  }      
	}                                 
	else {
	  error ("ReadInterpData: Illegal character in %s", base);
	}
      }
    }
    
    free(base);
    free(fmt);
    free(skip);
    
    return inlist;

  } else {
    
    return NULL;

  }
}
