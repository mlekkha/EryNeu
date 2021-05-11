/*****************************************************************
 *                                                               *
 *   solvers.h                                                   *
 *                                                               *
 *****************************************************************
 *                                                               *
 * solvers.h contains the interface to the solver functions, i.e *
 * solver function prototypes and the p_deriv global which ser-  *
 * ves as a definition of the interface to the derivative func.  *
 *                                                               *
 *****************************************************************
 *                                                               *
 * NOTE: *ONLY* general solvers allowed here; they *MUST* comply *
 *       to the generic solver interface (see solvers.c for de-  *
 *       tails)                                                  *
 *                                                               *
 *****************************************************************/

#ifndef SOLVERS_INCLUDED
#define SOLVERS_INCLUDED



/*** GLOBAL VARIABLES ******************************************************/

void (*p_deriv)(double *, double, double *, int);    
void (*d_deriv)(double *, double **, double, double *, int);
void (*p_jacobn)(double, double *, double *, double **, int);




/*** FUNCTION PROTOTYPES ***************************************************/

/*** Euler: propagates vin (of size n) from tin to tout by the Euler *******
 *          method; the result is returned by vout                         *
 ***************************************************************************/

void Euler(double *vin, double *vout, double tin, double tout, 
	   double stephint, double accuracy, int n, FILE *slog);



/*** Meuler: propagates vin (of size n) from tin to tout by the Modified ***
 *           Euler method (this is NOT the midpoint method, see Rk2());    *
 *           the result is returned by vout                                *
 ***************************************************************************/

void Meuler(double *vin, double *vout, double tin, double tout, 
	    double stephint, double accuracy, int n, FILE *slog);



/*** Heun: propagates vin (of size n) from tin to tout by Heun's method ****
 *         the result is returned by vout                                  *
 ***************************************************************************/

void Heun(double *vin, double *vout, double tin, double tout, 
	  double stephint, double accuracy, int n, FILE *slog);



/*** Rk2: propagates vin (of size n) from tin to tout by the Midpoint or ***
 *        Second-Order Runge-Kutta method; the result is returned by vout  *
 ***************************************************************************/

void Rk2(double *vin, double *vout, double tin, double tout, 
	 double stephint, double accuracy, int n, FILE *slog);



/*** Rk4: propagates vin (of size n) from tin to tout by the Fourth-Order **
 *        Runge-Kutta method; the result is returned by vout               *
 ***************************************************************************
 *                                                                         *
 * written by Joel Linton (somewhere around 1998)                          *
 * fixed and modified by Yoginho (somewhere around 2001)                   *
 *                                                                         *
 ***************************************************************************/

void Rk4(double *vin, double *vout, double tin, double tout, 
	 double stephint, double accuracy, int n, FILE *slog);



/*** Rkck: propagates vin (of size n) from tin to tout by the Runge-Kutta **
 *         Cash-Karp method, which is an adaptive-stepsize Rk method; it   *
 *         uses a fifth-order Rk formula with an embedded forth-oder for-  *
 *         mula for calucalting the error; its result is returned by vout  *
 ***************************************************************************
 *                                                                         *
 * This solver was written by Marcel Wolf, Spring 2002.                    *
 *                                                                         *
 ***************************************************************************/

void Rkck(double *vin, double *vout, double tin, double tout, 
	  double stephint, double accuracy, int n, FILE *slog);



/*** Rkf: propagates vin (of size n) from tin to tout by the Runge-Kutta ***
 *        Fehlberg method, which is a the original adaptive-stepsize Rk    *
 *        method (Cash-Karp is an improved version of this); it uses a     *
 *        fifth-order Rk formula with an embedded forth-oder formula for   *
 *        calucalting the error; its result is returned by vout            *
 ***************************************************************************
 *                                                                         *
 * This solver was written by Marcel Wolf, Spring 2002.                    *
 *                                                                         *
 ***************************************************************************/

void Rkf(double *vin, double *vout, double tin, double tout, 
	 double stephint, double accuracy, int n, FILE *slog);



/***** Milne: propagates vin (of size n) from tin to tout by Milne-Simpson *
 *            which is a predictor-corrector method; the result is retur-  *
 *            ned by vout                                                  *
 ***************************************************************************
 *                                                                         *
 * This solver was implemented by Konstantin Koslov, Dec 2001/Jan 2002     *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 * THIS SOLVER SEEMS TO BE BUGGY FOR SOME REASON, DO NOT USE IT!!!!!       *
 *                                                                         *
 ***************************************************************************/

void Milne(double *vin, double *vout, double tin, double tout, 
	   double stephint, double accuracy, int n, FILE *slog);



/***** Adams: propagates vin (of size n) from tin to tout by Adams-Moulton *
 *            which is an implicit predictor-corrector method of second    *
 *            order; the result is returned by vout                        *
 ***************************************************************************
 *                                                                         *
 * This solver was implemented by Konstantin Koslov, Spring 2002           *
 * Slightly modified by Manu, July 2002                                    *
 *                                                                         *
 ***************************************************************************/

void Adams(double *vin, double *vout, double tin, double tout, 
	   double stephint, double accuracy, int n, FILE *slog); 



/***** BuST: propagates v(t) from t1 to t2 by Bulirsch-Stoer; this method **
 *           uses Richardson extrapolation to estimate v's at a hypothe-   *
 *           tical stepsize of 0; the extrapolation also yields an error   *
 *           estimate, which is used to adapt stepsize and change the or-  *
 *           der of the method as required; the result is returned by vout *
 ***************************************************************************
 *                                                                         *
 * This solver was implemented by Manu, July 2002                          *
 *                                                                         *
 ***************************************************************************/

void BuSt(double *vin, double *vout, double tin, double tout, 
	  double stephint, double accuracy, int n, FILE *slog);



/***** BaDe: propagates v(t) from t1 to t2 by Bader-Deuflhard; this method *
 *           uses Richardson extrapolation to estimate v's at a hypothe-   *
 *           tical stepsize of 0; the extrapolation also yields an error   *
 *           estimate, which is used to adapt stepsize and change the or-  *
 *           der of the method as required; the result is returned by vout *
 ***************************************************************************
 *                                                                         *
 * This solver was implemented by Yogi, based on BuSt, Aug 2002            *
 *                                                                         *
 ***************************************************************************/

void BaDe(double *vin, double *vout, double tin, double tout, 
	  double stephint, double accuracy, int n, FILE *slog);



/***** WriteSolvLog: write to solver log file ******************************/ 

void WriteSolvLog(char *solver, double tin, double tout, double h, int n, 
		  int nderivs, FILE *slog);

double *Construct_Discont_Array(double range, double *taus, int n, double *starts, int sn, int *disc_size);

int compare(double *x, double *y);

void SoDe(double *vin, double *vout, double tin, double tout, 
	 double stephint, double accuracy, int n, FILE *slog);


int y_delayed(double ***vd, int n, double *rktimes, double *tau,
double *grid, double **vdone, double **deriv1, double **deriv2, double
**deriv3, double **deriv4, int gridsize, double accu);


void DCERk32(double **vatt, int n, double *tarray, int tpoints, double
			*darray, int dpoints, double stephint, double accuracy);


void CE(double t, double *vans, double tbegin, double *v_at_tbegin,
double ech, double *d1, double *d2, double *d3, double *d4, int n);

void DivideHistory(double t1, double t2);
void FreeDelaySolver(void);
void InitDelaySolver(void);

#endif
