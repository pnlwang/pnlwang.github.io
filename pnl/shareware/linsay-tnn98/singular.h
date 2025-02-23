       /*****************************************************************/
       /* This is the header file for the C program "singular.c"        */
       /* Written by DeLiang Wang, 7/96, at OSU.                        */
       /*****************************************************************/
        
#include <stdio.h>
#include <math.h>
#define STEPS   300	/* number of steps to get output */
#define delta   0.12	/* inteval for obtaining an output */
#define baseline -0.02  /* baseline oscillator input */
#define N1	50	/* number of rows */
#define N2      50       /* number of columns */
#define error   0.0001    /* tolerable error to save computing */

#define M	2501   	/* number of displayed oscillators */
