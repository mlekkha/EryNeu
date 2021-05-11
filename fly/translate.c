/*****************************************************************
 *                                                               *
 *   translate.c                                                 *
 *                                                               *
 *****************************************************************
 *                                                               *
 *   written by JR, modified by Yoginho                          *
 *                                                               *
 *****************************************************************
 *                                                               *
 * SENDING PARAMATERS TO THE ANNEALER                            *
 *                                                               *
 * This file is for making the list of pointers that the         *
 * annealer uses to make moves. The idea is to conceal as much   *
 * problem specific details as possible from the annealer. Some  *
 * of the more general ways of doing this must wait for future   *
 * code. e.g., when there are multiple possible right hand-side- *
 * of-the-ODE derivative functions, each of these must have a    *
 * translation function, and possibly a small family of other    *
 * functions as well. These would need to be handled as arrays   *
 * of pointers to functions for each choice, etc. Also, the      *
 * ParamList structure now has ranges, because it may be impor-  *
 * tant in some contexts for the annealer to know the limits of  *
 * the search space. If penalty funcs are being used (pen_vec != *
 * NULL), a little more work is required, and the coder might,   *
 * say, write a function that goes in this file which calculates *
 * the effective range (say for penalty <= 10000) of one param   *
 * while asuming the values of the others to be fixed.           *
 *                                                               *
 * SPECIFYING WHICH PARAMETERS ARE ANNEALED, WHICH HELD CONSTANT *
 *                                                               *
 * Parameters are specified as fixed one by one in the parameter *
 * file, by setting their entries in the $tweak section to zero. *
 * this entry is then used to construct the array of pointers to *
 * the parameters that are to be tweaked.                        *
 *                                                               *
 *****************************************************************
 *                                                               *
 * Copyright (C) 1989-2003 John Reinitz                          *
 * the full GPL copyright notice can be found in lsa.c           *
 *                                                               *
 *****************************************************************/
                                                          
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include <error.h>
#include <maternal.h>                               /* for the defs global */
#include <moves.h>
#include <score.h>                                      /* for SearchSpace */
#include <zygotic.h>                                        /* for EqParms */





/*** A STATIC VARIABLE *****************************************************/

Tweak         tweak;       /* tells the annealer which parameters to tweak */




/*** FUNCTION DEFINITIONS **************************************************/

/* Initialization */

/*** InitTweak: installs tweak as a static variable in translate.c; tweak **
 *              is read from the $tweak section in newstyle data files     *
 ***************************************************************************/

void InitTweak(FILE *fp)
{
  tweak = ReadTweak(fp);
}





/* Translation: installing an array of pointer to params-to-be-tweaked */

/*** Translate: creates an array of pointers that point to the parameters **
 *              to be tweaked; the PArrPtr that is returned also includes  *
 *              pointers to the corresponding parameter ranges for each    *
 *              rameters, although this feature is not yet used anywhere   *
 *              in the annealing code                                      *
 *     CAUTION: InitZygote and InitScoring have to be called first!        *
 ***************************************************************************/

PArrPtr Translate(void)
{
  int         i, j;
  int         index    = 0;

  PArrPtr     plist;                          /* local copy of the PArrPtr */
  ParamList   *p;                          /* local copy of parameter list */

  EqParms     *parm;                          /* local copy of the EqParms */
  SearchSpace *limits;                    /* local copy of the SearchSpace */

  int         max_size = 0;             /* maximum size of parameter array */
  int         max_elem = 0;    /* maximum # of elements of parameter array */
 
/* initially we allocate max_elem members of the PArrPtr array in max_size */
/* bytes. The actual space needed will be less than or equal to this.      */

  max_elem = ((defs.ngenes*defs.ngenes) 
  					+ (defs.ngenes*defs.egenes) 
							+ 6 * defs.ngenes); 

  max_size = max_elem * sizeof(ParamList);

/* initialize the ParamList */

  p = (ParamList *)malloc(max_size); 

  for (i=0; i<max_elem; i++) {
    p[i].param       = NULL;
    p[i].param_range = NULL;
  }

  plist.size = 0;

/* Get limits and parameters */

  parm   = GetParameters();   /* get pointer to EqParm struct in zygotic.c */
  limits = GetLimits();           /* get pointer to SearchSpace in score.c */

/* if we are using a penalty, not all of Limits will have been allocated   */
/* and we must take care not to take indices of unallocated arrays!!!      */

  if ( limits->pen_vec != NULL )  /* we are using a penalty, hence pen_vec */
    plist.pen_vec = limits->pen_vec;  
  else                               /* or: explicit ranges for everything */
    plist.pen_vec = NULL;                

  for ( i=0; i < defs.ngenes; i++ )                 /* pointers to R stuff */
    if ( tweak.Rtweak[i] == 1 ) {    
      p[index].param = &parm->R[i];
      p[index].param_range = limits->Rlim[i];
      index++;
    }

  for ( i=0; i < defs.ngenes; i++ )                 /* pointers to T stuff */
    for ( j=0; j < defs.ngenes; j++ )
      if ( tweak.Ttweak[i*defs.ngenes+j] == 1 ) {
	p[index].param = &parm->T[j+i*defs.ngenes];
	if ( plist.pen_vec == NULL )
	  p[index].param_range = limits->Tlim[j+i*defs.ngenes];
	index++;
      }

  /* 01/13/10 Manu: If there are no external inputs, the for
   * loop below should never be entered */

  for ( i=0; i < defs.ngenes; i++ )                 /* pointers to E stuff */
    for ( j=0; j < defs.egenes; j++ )
      if ( tweak.Etweak[i*defs.egenes+j] == 1 ) {
	p[index].param = &parm->E[j+i*defs.egenes];
	if ( plist.pen_vec == NULL )
	  p[index].param_range = limits->Elim[j+i*defs.egenes];
	index++;
      }

  for ( i=0; i < defs.ngenes; i++ )                 /* pointers to m stuff */
    if ( tweak.mtweak[i] == 1 ) {
      p[index].param = &parm->m[i];
      if ( plist.pen_vec == NULL )
	p[index].param_range = limits->mlim[i];
      index++;
    }

  for ( i=0; i < defs.ngenes; i++ )                 /* pointers to h stuff */
    if ( tweak.htweak[i] == 1 ) {
      p[index].param = &parm->h[i];
      if ( plist.pen_vec == NULL )
	p[index].param_range = limits->hlim[i];
      index++;
    }

  for ( i=0; i < defs.ngenes; i++ )                 /* pointers to d stuff */
    if ( tweak.dtweak[i] == 1 ) {
      p[index].param = &parm->d[i];
      p[index].param_range = limits->dlim[i];
      index++;
      if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') )
	break;
    }

  for ( i=0; i < defs.ngenes; i++ )            /* pointers to lambda stuff */
    if ( tweak.lambdatweak[i] == 1 ) {
      p[index].param = &parm->lambda[i];
      p[index].param_range = limits->lambdalim[i];
      index++;
    }

  for ( i=0; i < defs.ngenes; i++ )            /* pointers to tau stuff */
    if ( tweak.tautweak[i] == 1 ) {
      p[index].param = &parm->tau[i];
      p[index].param_range = limits->taulim[i];
      index++;
    }


/* reallocate memory for smaller-than-or-equal-to-maxsize array */

  p = (ParamList *)realloc((void *)p, (index * sizeof(ParamList))); 

/* finish up initializing PArrPtr */

  plist.size  = index;
  plist.array = p;

  return plist;
}





/* Reading: following reads tweak struct from file ($tweak section) */

/*** ReadTweak: reads the tweak array from the $tweak section in the data **
 *              file; this array has a value of 1 or 0 for each parameter  *
 *              in the model and is used by Translate to create the array  *
 *              of pointers to the parameters-to-be-tweaked                *
 ***************************************************************************/

Tweak ReadTweak(FILE *fp)
{
  Tweak             l_tweak;                         /* local Tweak struct */
  int               *temptweak,*temptweak1;          /* temporary arrays to read tweaks */

  int               i;                               /* local loop counter */
  int               c;                         /* used to parse text lines */
  int               linecount = 0;        /* keep track of # of lines read */
  int               Tcount    = 0;           /* keep track of T lines read */
  int               Ecount    = 0;           /* keep track of E lines read */

  char              *base;          /* pointer to beginning of line string */
  char              *record;    /* string for reading whole line of params */

  char              **fmt;   /* array of format strings for reading params */
  char              **fmt1;   /* array of format strings for reading
  E tweaks */
  char              *skip,*skip1;     /* string of values to be skipped */

  const char        read_fmt[] = "%d";                      /* read an int */
  const char        skip_fmt[] = "%*d ";                  /* ignore an int */

  base = (char *)calloc(MAX_RECORD, sizeof(char *));

  skip = (char *)calloc(MAX_RECORD, sizeof(char *));

  skip1 = (char *)calloc(MAX_RECORD, sizeof(char *));

  fmt  = (char **)calloc(defs.ngenes, sizeof(char *));
  if (defs.egenes > 0)
	  fmt1  = (char **)calloc(defs.egenes, sizeof(char *));

  temptweak = (int *)calloc(defs.ngenes, sizeof(int *));
  temptweak1 = (int *)calloc(defs.egenes, sizeof(int *));

/* create format strings according to the number of genes */

  for ( i=0; i<defs.ngenes; i++ ) {  
    fmt[i] = (char *)calloc(MAX_RECORD, sizeof(char));   
    fmt[i] = strcpy(fmt[i], skip);
    fmt[i] = strcat(fmt[i], read_fmt);
    skip   = strcat(skip, skip_fmt);
  }

/* create format strings according to the number of external inputs */


    if (defs.egenes > 0) {

	  for ( i=0; i<defs.egenes; i++ ) {  
		fmt1[i] = (char *)calloc(MAX_RECORD, sizeof(char));   
		fmt1[i] = strcpy(fmt1[i], skip1);
		fmt1[i] = strcat(fmt1[i], read_fmt);
		skip1   = strcat(skip1, skip_fmt);
	  }

	}  

/* initialize the Tweak struct */

  l_tweak.Rtweak = (int *)calloc(defs.ngenes, sizeof(int));
  l_tweak.Ttweak = (int *)calloc(defs.ngenes * defs.ngenes, sizeof(int));

  if (defs.egenes > 0)
	  l_tweak.Etweak = (int *)calloc(defs.ngenes * defs.egenes, sizeof(int));

  l_tweak.mtweak = (int *)calloc(defs.ngenes, sizeof(int));
  l_tweak.htweak = (int *)calloc(defs.ngenes, sizeof(int));
  if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) {
    l_tweak.dtweak = (int *)malloc(sizeof(int));
  } else {
    l_tweak.dtweak = (int *)calloc(defs.ngenes, sizeof(int));
  }
  l_tweak.lambdatweak = (int *)calloc(defs.ngenes, sizeof(int));
  l_tweak.tautweak = (int *)calloc(defs.ngenes, sizeof(int));


  fp = FindSection(fp, "tweak");                     /* find tweak section */
  if( !fp )
    error("ReadTweak: could not locate tweak\n");
  
  while ( strncmp(( base=fgets(base, MAX_RECORD, fp)), "$$", 2)) {
    
    record = base;

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
	  if ( 1 != sscanf(record, fmt[0], &temptweak[0]) )
	    error("ReadTweak: error reading tweaks");
	} else if (linecount == 2) {     
	  for ( i=0; i < defs.egenes; i++ ) {
	    if ( 1 != sscanf(record, fmt1[i], &temptweak1[i]) )
	      error("ReadTweak: error reading tweak variables");
	  }
	} else {     
	  for ( i=0; i < defs.ngenes; i++ ) {
	    if ( 1 != sscanf(record, fmt[i], &temptweak[i]) )
	      error("ReadTweak: error reading tweak variables");
	  }
	}

	switch (linecount) {  /* copy read parameters into the right array */
	case 0:
	  for ( i=0; i < defs.ngenes; i++ )                    /* R tweaks */
	    l_tweak.Rtweak[i] = temptweak[i];
	  linecount++;
	  break;
	case 1:          /* T tweaks: keep track of read lines with Tcount */
	  for ( i=0; i < defs.ngenes; i++ )
	    l_tweak.Ttweak[i+Tcount*defs.ngenes] 
	      = temptweak[i];
	  Tcount++;
	  if ( Tcount == defs.ngenes ) 
	    linecount++;	  
	  break;
	case 2:          /* E tweaks: keep track of read lines with Ecount */
		if (defs.egenes > 0) {
		  for ( i=0; i < defs.egenes; i++ )
			l_tweak.Etweak[i+Ecount*defs.egenes] 
			  = temptweak1[i];
		  Ecount++;
		  if ( Ecount == defs.ngenes ) 
			linecount++;	  
		}

	  break;
	case 3:                                                /* m tweaks */
	  for ( i=0; i < defs.ngenes; i++ )
	    l_tweak.mtweak[i] = temptweak[i];
	  linecount++;
	  break;
	case 4:
	  for ( i=0; i < defs.ngenes; i++ )                    /* h tweaks */
	    l_tweak.htweak[i] = temptweak[i];
	  linecount++;
	  break;
	case 5:                       /* d tweaks: consider diff. schedule */
	  if ((defs.diff_schedule == 'A') || (defs.diff_schedule == 'C' )) {
	    l_tweak.dtweak[0] = temptweak[0];
	  } else {
	    for ( i=0; i < defs.ngenes; i++ )
	      l_tweak.dtweak[i] = temptweak[i];
	  }
	  linecount++;
	  break;
	case 6:                                           /* lambda tweaks */
	  for ( i=0; i < defs.ngenes; i++ )
	    l_tweak.lambdatweak[i] = temptweak[i];
	  linecount++;
	  break;
	case 7:                                           /* lambda tweaks */
	  for ( i=0; i < defs.ngenes; i++ )
	    l_tweak.tautweak[i] = temptweak[i];
	  linecount++;
	  break;
	default:
	  error("ReadTweak: too many data lines in tweak section");
	}
	break;                           /* don't do rest of loop anymore! */
      } 

      else if ( isspace(c) ) {               /* ignore leading white space */
	c = (int)*(++record);
      }

      else {                  /* anything but space or digit means comment */
	break;
      } 
    }
  }
    
  free(temptweak);
  free(temptweak1);
  free(base);
  free(skip);
  free(skip1);
  
  for (i=0; i<defs.ngenes; i++)
    free(fmt[i]);
  free(fmt);

  if (defs.egenes > 0) {
	  for (i=0; i<defs.egenes; i++)
		free(fmt1[i]);
	  free(fmt1);
  }	  

  return l_tweak;
}
