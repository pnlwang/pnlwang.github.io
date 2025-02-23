       /*****************************************************************/
       /* This is a header file for the program that implements the     */
       /* selective gating mechanism. By DeLiang Wang, 8/93, at OSU     */
       /*****************************************************************/
       
#include <stdio.h>
#include <math.h>
#define WIDTH   8000	/* integration time */
#define LEFT    0	/* starting display time */
#define WIDTH2  320	/* display time steps, 320, or 500 previously */
#define N1	20	/* number of rows */
#define N2      20       /* number of columns */
#define M	401  	/* number of displayed oscillators */

