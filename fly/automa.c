#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <unistd.h>                                   /* for getopt */
#include <ctype.h>

#include "error.h"   
#include "integrate.h"
#include "maternal.h"
#include "solvers.h"
#include "zygotic.h"
 

#define UCOL 2	/* location (column) of u */


Dlist *ReadGutsData(FILE *fp , char *section, int *ndp, int numguts);

NArrPtr List2Guts(Dlist *inlist, int numguts);

Dlist *init_Dlistflex(int numguts);

void PrintGutsDumb(FILE *fp, NArrPtr table, char *id, int ndigits, 
														int numguts);
void FindUDrivers(double *v, int arrsize, int arrwid, char
					**components, char **iden, float ratio);

void UDrSrch(double *v, int arrsize, int arrwid, char **iden, 
int framemax, int framemin, int start, int stop, int numcomps, 
char **components, float ratio);

void DiffsMethod(double *v, int arrsize, int arrwid, char **iden, 
int numcomps, char **components, float ratio);

double TrapInt(double *v, int column, int start, int end, int
															arrwid); 
double norm_2(double *v, int column, int start, int end, int
															arrwid); 

char *CalcNorms(double *lc, int stripenum, int segstart, int segend, 
int arrsize, int arrwid, char **iden, float ratio);




int main(int argc, char **argv)
{


char            *infile;               /* pointer to input file name */
FILE            *fp;                   /* pointer to input data file */

int           	i,j,k, ndp;          /* loop counter and dummy variable
													for ReadGutsData*/
int				numguts, numnucs; /* Number of Guts Elements calculated
									 and number of nuclei */

NArrPtr	  		goutput;	      /* Output of Gutz is stored in this */
Dlist			*inlist; /* ReadGutsData reads guts into this */
char 			**componenres=NULL; /*Array of Strings to hold
									components */

/* Array of Strings to hold Guts defs guts title and parsed strings*/

char 			**gutsdefs = NULL, *gutstitle = NULL, 
				**dummy = NULL, **w = NULL; 



/* let's get started and open the data file here */

if (argc != 5) 
	error("Automa: Usage: automa <datafile> <gutsoutputfile> <epsilon>\
	<no. of guts being searched>");

infile = argv[1];
fp = fopen(infile, "r");      
if( !fp ) {
  perror("unfold");
  exit(1);
}

/* use zygotic, maternal routines to read data */

gutsdefs = ReadGuts(fp);	
fclose(fp);

infile = argv[2];
fp = fopen(infile, "r");      
if( !fp ) {
	perror("unfold");
   	exit(1);
}

gutstitle = (char *) calloc(MAX_RECORD, sizeof(char));
dummy = (char **) calloc(MAX_RECORD, sizeof(char *));

for (i=0; *(gutsdefs+i); i++) {
	
	w = dummy;
	numguts = ParseString(*(gutsdefs+i),w)-1;

	if ((numguts == -1)) {
		while (*++w) {
			free(*w);
		}
		free(*dummy);
		break;
	}

	if ((numguts == 0)) {
		printf("Ignoring %s as no guts requirement specified", 
		*(gutsdefs+i));
		while (*++w) {
			free(*w);
		}
		free(*dummy);
		break;
	}
	
	gutstitle = strcat(strcpy(gutstitle, "guts_for_"),*(gutsdefs+i));
	
/*	printf("Reading the guts section, %s, %d\n", gutstitle,numguts);*/



    if( !(inlist=ReadGutsData(fp, gutstitle, &ndp, numguts)) )  
      error("Automa: error reading ");
    

    goutput = List2Guts(inlist, numguts);
      
	for (j = 0; j < goutput.size; j++) {

/*		printf("Doing Analysis for t = %f\n",goutput.array[j].time);*/

		/* Allocating space for the array of strings to hold components
		for each nucleus */

		numnucs = goutput.array[j].state.size / (numguts + 2); 
/*		printf("%d \n", numnucs);*/
		componenres = (char **) calloc(numnucs, sizeof(char *));

		/* Allocating space for the string that would hold the results
		for a particular nucleus - the strlen cannot be greater than the
		number of guts! */
		
/*		for (k=0; k < numnucs; k++) 
			*(componenres+k) = (char *) calloc(numguts, sizeof(char));
		printf("Abdee Kabdee!\n");	*/

/* The first algorithm tried - based on detecting eve stripes, by sign
changes in u, now commented out */

/*		FindUDrivers(goutput.array[j].state.array,
						goutput.array[j].state.size,
						numguts+2,componenres,w, atof(argv[3]));
	


		for (k=0; k < numnucs; k++)
			printf("%d %s\n",
			(int) goutput.array[j].state.array[(numguts+2)*k], 
			*(componenres+k));

		for (k=0; k < numnucs; k++) 
			free(*(componenres+k));

		free(componenres);
		componenres = (char **) calloc(numnucs, sizeof(char *));

		
		for (k=0; k < numnucs; k++) 
			*(componenres+k) = (char *) calloc(numguts,
			sizeof(char));*/

/* The second approach tried - to search through the different frames
possible, leads to non-unique results, now commented out */

/*		UDrSrch(goutput.array[j].state.array,
						goutput.array[j].state.size,
						numguts+2,w,1,1,0,numnucs - 1,
						atoi(argv[4]), componenres, atof(argv[3]));


		for (k=0; k < numnucs; k++)
			printf("%d %s\n",
			(int) goutput.array[j].state.array[(numguts+2)*k], 
			*(componenres+k));

		for (k=0; k < numnucs; k++) 
			free(*(componenres+k));

		free(componenres); */

		componenres = (char **) calloc(numnucs, sizeof(char *));

		for (k=0; k < numnucs; k++) 
			*(componenres+k) = (char *) calloc(numguts, sizeof(char));

/* The current favourite, calculates which components contribute more
than a certain ratio of the change component to change in u from nucleus to nucleus,
ration supplied as "epsilon squared" in command line arguments. 
Further, it
takes a second argument - the maximum number of components to be
reported (n). If the number of components above epsilon squared is
greater than n, they are not reported, anything less than or equal to
it is. */

		DiffsMethod(goutput.array[j].state.array, 
							goutput.array[j].state.size, 
							numguts+2,w, atoi(argv[4]), 
							componenres, atof(argv[3]));

		/* print the output */

		for (k=1; k < numnucs; k++)
			printf("%d %s\n",
			(int) goutput.array[j].state.array[(numguts+2)*k], 
			*(componenres+k));

		for (k=0; k < numnucs; k++) 
			free(*(componenres+k));

		free(componenres);
 	}

	/*	FreeGuts(goutput); */ 
	free(*(gutsdefs+i));
	while (*++w) {
		free(*w);
	}
	free(*dummy);
    free_Dlist(inlist);
}
free(gutsdefs);		
free(gutstitle);
free(dummy);
fclose(fp);

return 0;
}



Dlist *ReadGutsData(FILE *fp , char *section, int *ndp, int numguts)
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
	  current = init_Dlistflex(numguts);
	  
	  if ( 1 != sscanf(record, "%d ", &(current->lineage)) )
	    error("ReadGutsData: error reading %s", base);
	  
/* the following loop reads a line of data of variable length into a d-    */
/* array, skipping lineage number and all previously read data values      */
	  
	  skip = strcpy(skip, init_fmt); /* format to be skipped by sscanf */
	  
	  for(i=0; i < (numguts + 1); i++) {
	    
	    fmt  = strcpy(fmt, skip);     /* all this stuff here is to get */
	    fmt  = strcat(fmt, read_fmt);  /* the fmt string for sscanf to */
	    skip = strcat(skip, skip_fmt);           /* the correct format */
	    
	    if ( 1 != sscanf(record, (const char *)fmt, &(current->d[i])) )
	      error("ReadGutsData: error reading %s", base);

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
	  error ("ReadGutsData: Illegal character in %s", base);
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





NArrPtr List2Guts(Dlist *inlist, int numguts)
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
	(double *)calloc(n*(numguts+2), sizeof(double)); 
      bias.array[bias.size-1].state.size = n*(numguts+2);  

      i=0;
    }
      
/* always: read lin numbers, times and concs into array */

bias.array[bias.size - 1].state.array[i++] = (double)
current->lineage;

bias.array[bias.size - 1].state.array[i++] = current->d[0];

    for( j=1; j <= numguts; j++ ) {
      bias.array[bias.size - 1].state.array[i]
	= current->d[j];
      i++;
    }
  }

  return bias;
}


Dlist *init_Dlistflex(int numguts)  
{
  Dlist    *p;                               

  if( (p=(Dlist *)malloc(sizeof(Dlist))) ) {
    
    if ( !(p->d=(double *)calloc(numguts+1, sizeof(double))) )
      error("initDlistflex: couldn't allocate!");    

    p->lineage = 0;
    p->next    = NULL;
  
}
  else 
    error("initDlist: couldn't allocate!");

  return p;
}


void PrintGutsDumb(FILE *fp, NArrPtr table, char *id, int ndigits, int numguts)
{
  int                i, j, k;                       /* local loop counters */

  int                lineage;                /* lineage number for nucleus */
  int                nnucs  = 0;         /* how many nuclei at given time? */

  fprintf(fp, "$%s\n\n", id); 
  
  for (i=0; i<table.size; i++) {
    for(j=0; j<(table.array[i].state.size/(numguts+2)); j++) {     
      fprintf(fp, "%5d %7.3f", (int)table.array[i].state.array[(j*(numguts+2))], table.array[i].state.array[1+(j*(numguts+2))]);
      for (k=2; k < (numguts+2); k++)
	fprintf(fp, " %*.*f", ndigits+4, ndigits,
		table.array[i].state.array[k+(j*(numguts+2))]);
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }                             
  fprintf(fp,"$$\n");
  fflush(fp);
}

void UDrSrch(double *v, int arrsize, int arrwid, char **iden, 
int framemax, int framemin, int start, int stop, int numcomps, 
char **components, float ratio)
{

	int i,j,k; 					/* loop counters */
	int len;					/* Length of segment */
	int framesize; 				/* The current framesize being looked at
																	*/
	int flag = 0;				/* Flag to decide about exiting loop 
																	*/
	double *lc;					/* local copy of v */
	char compsans[MAX_RECORD]; 	/* Temporary place to store the
								components identified for a frame */




	lc = (double *) calloc(arrsize*arrwid, sizeof(double));

	for (framesize = framemax; framesize >= framemin; framesize--)
		for (i=start; i <= stop - framesize + 1; i++) {

			for (j = 0; j < arrwid; j++)
				for (k = 0; k < arrsize / arrwid; k++)
					lc[j + k*arrwid] = v[j + k*arrwid];

			strcpy(compsans, CalcNorms(lc, 0, i,
					(i + framesize - 1), arrsize, arrwid, iden, ratio));
			
			if (strlen(compsans) == numcomps) {

				printf("Found %s in %d - %d\n",compsans,i,i +
				framesize - 1);

				for (k = i; k < i + framesize; k++)
					strcpy(*(components + k), compsans);

				/* Call UDrSrch recursively below, first left of
				determined segment */

				len = i;
				if (len >= framemin)
					if (len <= framesize - 1) 
					
						UDrSrch(v, arrsize, arrwid, iden, len, framemin,
						start, i - 1, numcomps, components, ratio);
					
					else UDrSrch(v, arrsize, arrwid, iden, framesize
						 - 1, framemin, start, i - 1, numcomps, 
						 components, ratio);
				
				/* Now do the same for the region right of the
				determined segment */
				
				len = stop  - (i + framesize) + 1; 
				
				if (len >= framemin)
					if (len <= framesize)

						UDrSrch(v, arrsize, arrwid, iden, len, framemin,
						i + framesize, stop, numcomps, components,
						ratio);
					
					else UDrSrch(v, arrsize, arrwid, iden, framesize,
						 framemin, i + framesize, stop, numcomps, 
						 components, ratio);
				
				return; /* Your job is done, return to the calling
								routine */
			}
		}

	free(lc);
}

void DiffsMethod(double *v, int arrsize, int arrwid, char **iden, 
int numcomps, char **components, float ratio)
{

	int i,j,k; 					/* Counters */
	char compsans[MAX_RECORD]; 	/* To store the components that
										Calcnorms returns */
	
	double *lc;		/* local copy of v */


	lc = (double *) calloc(arrsize*arrwid, sizeof(double));


	for (j = 2; j < arrwid; j++)
		for (k = 1; k < arrsize / arrwid; k++)
			lc[j + k*arrwid] = v[j + k*arrwid] - v[j + (k-1)*arrwid];

	for (i=1; i < arrsize / arrwid; i++) {

		strcpy(compsans, CalcNorms(lc, 0, i, i, arrsize, arrwid, iden, 
														ratio));
			
		if (strlen(compsans) <= numcomps) {

/*			printf("Found %s in %d - %d\n",compsans,i-1,i);*/
			strcpy(*(components + i), compsans);
		}
	}
	
	free(lc);
}

void FindUDrivers(double *v, int arrsize, int arrwid, char
**components, char **iden, float ratio)
{

	int in = 0; 	/* whether you are in a +ve region or not */
	int i,j,k; 		/* loop counters */
	int segstart, segend=0; 			/* start and end of the most 
										recent +ve region */
	int stripenum = 0; /* to count which stripe we have discovered
																*/

	double *lc;		/* local copy of v */
	char compsans[MAX_RECORD];
	NArrPtr temp;	/* to print the normalized guts */

	lc = (double *) calloc(arrsize*arrwid, sizeof(double));

	for (i=0; i < arrsize / arrwid; i++) {

	/* if you encounter a +ve element and you are outside a +ve region,
	start the counter and toggle in, set end to 0 so that you don't
	process this region until the end is reached (and defined) */

		if (v[UCOL + arrwid*i] >= 0.) {
			if (!in) {				
				
				segstart = i;
				segend = 0;
				in = !in;

			} else if (i+1 == arrsize / arrwid) {
					
					stripenum++;
					printf("Found Stripe No. %d, It is\n", stripenum);
					for (k = segstart; k <= i; k++)
						printf("k: %d, u(k)=%f\n",k,v[UCOL +
														k*arrwid]);

					for (j = 0; j < arrwid; j++)
						for (k = 0; k < arrsize / arrwid; k++)
							lc[j + k*arrwid] = v[j + k*arrwid];

					strcpy(compsans, CalcNorms(lc, stripenum, segstart,
					i, arrsize, arrwid, iden, ratio));
					
					printf("The drivers for stripe %d, are "
					"%s\n", stripenum, compsans);
					
					for (k=segstart; k<=i; k++)
						strcpy(*(components + k), compsans);
			}
		}
	/* if it is negative and you are in a +ve region, define the
	end, and toggle in */

		else {
			
			if (in) {
				
				in = !in;
				segend = i-1;
				
				if (segend - segstart) { 
					
					stripenum++;
					printf("Found Stripe No. %d, It is\n", stripenum);
					
					for (k = segstart; k <= segend; k++)
						printf("k: %d, u(k)=%f\n",k,v[UCOL +
														k*arrwid]);

					for (j = 0; j < arrwid; j++)
						for (k = 0; k < arrsize / arrwid; k++)
							lc[j + k*arrwid] = v[j + k*arrwid];

					strcpy(compsans, CalcNorms(lc, stripenum, segstart,
					segend, arrsize, arrwid, iden, ratio));
					
					printf("The drivers for stripe %d, are "
					"%s, %d\n", stripenum, compsans,segstart + (segend
					- segstart)/2);
				
					for (k=segstart; k<=segend; k++)
						strcpy(*(components + k), compsans);

				} else segend = 0;

			} else segend = 0;
		}

	}

	free(lc);
}

/* Given the start of a segment and its end, calculates the components
of u that were above ratio of u. It first reduces the average to zero,
to leave constant terms. Then it just simple compares the areas under
the curves of the different components with u, determining the
dominant ones, and returns a string of symbols, of the dominant ones
based on the gene id string. */

char *CalcNorms(double *lc, int stripenum, int segstart, int segend, 
int arrsize, int arrwid, char **iden, float ratio)
{

	int colcount;	/* counter to traverse the different u
										contributions */
	int i,j,k; 		/* loop counters */

	double average;			/* average to be found by TrapInt to 
							remove the constant fourier component */
	double unorm, sqnorm;

	char segcomps[MAX_RECORD]="";


/*			for (colcount = UCOL; colcount < arrwid; colcount++) {

				sqnorm = norm_2(lc, colcount, segstart, segend, arrwid);
				
				printf("Before Normalization For %s, the norm for
stripe %d is %f\n", *(iden+colcount-1),stripenum,sqnorm);  
			}*/

			/* Remove the average if segment has more than one nuclei
			*/

			if (segend - segstart) {
			
				average = TrapInt(lc, UCOL, segstart, segend, arrwid)/
												(segend - segstart);	
				for (k = 0; k < arrsize / arrwid; k++)
					lc[UCOL + k*arrwid] = lc[UCOL + k*arrwid] - average;
			}

			/* Calculate the norm of u expression in stripe */
			unorm = norm_2(lc, UCOL, segstart, segend, arrwid);
				
/*			printf("For %s, the norm for stripe %d is %f\n",
			*(iden+1),stripenum,unorm);*/ 

			/* Calculate the norms of the others, and see if they
			are expressed */


			for (colcount = UCOL+1; colcount < arrwid; colcount++) {

				/* Like above remove average if it exists */
				
				if (segend - segstart) {
				
					average = TrapInt(lc, colcount, segstart, segend, 
										arrwid)/(segend - segstart);	

					for (k = 0; k < arrsize / arrwid; k++)
						lc[colcount + k*arrwid] = lc[colcount + 
													k*arrwid] - average;
				}
				
				sqnorm = norm_2(lc, colcount, segstart, segend, arrwid);
				
				if (sqnorm > ratio*unorm) {
					
					strcat(segcomps,*(iden+colcount-1));
/*					printf("Found...%s\n",*(iden+colcount-1));*/

				}

/*				printf("After Normalization For %s, the norm for
stripe %d is %f\n", *(iden+colcount-1),stripenum,sqnorm); */
			} 
	
		return segcomps;
}

/* Calculates the area under a curve using the trapezoid rule */

double TrapInt(double *v, int column, int start, int end, int
arrwid) 
{

	int i;			/* loop counter */
	double area = 0; /* where you store the integral */

	if (end - start) {

		area = 0.5 * (v[column + start*arrwid] + v[column +
													end*arrwid]);

		if (end - start > 1)
			for (i = start+1; i < end; i++)
				area += v[column + i*arrwid];
	}

	return area;
}

/* Calculates the 2-norm of a vector */

double norm_2(double *v, int column, int start, int end, int
arrwid) 
{

	int i;		/* loop counter */
	double normsq = 0;/* where you hold the square of the norm */

	for (i = start; i <= end; i++) 
		normsq += v[column + i*arrwid]*v[column + i*arrwid];

	return normsq;
}		
