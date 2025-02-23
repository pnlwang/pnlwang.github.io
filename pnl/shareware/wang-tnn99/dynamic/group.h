       /*****************************************************************/
       /* This is the header file for the program that models           */
       /* speech segregation. By DeLiang Wang, 8/97, at OSU             */
       /*****************************************************************/
        
#include <stdio.h>
#include <math.h>
#define STEPS   200	/* number of steps to get output */
#define delta   0.12	/* inteval for obtaining an output */
#define baseline -0.02  /* baseline oscillator input */
#define N1	127	/* number of rows */
#define N2      150     /* number of columns */
#define SEG     95     /* number of segments + 1 */
#define error   0.0001  /* tolerable error to save computing */

#define M	64 	/* number of displayed oscillators */
#define Ms	258 	/* number of displayed oscillators for speech */
#define Mp	216 	/* number of displayed oscillators for phone */
#define COL     44       /* the relative column to be displayed */
