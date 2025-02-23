       /*****************************************************************/
       /* This is the header file for the LEGION program that does      */
       /* 2D image segmentation. By DeLiang Wang, 7/94, at OSU          */
       /*****************************************************************/
       
#include <stdio.h>
#include <math.h>
#define WIDTH   1200	/* integration time */
#define N1	160	/* number of rows */
#define N2      160       /* number of columns */
#define KS      5      /* the diameter of a kernel */
#define LP      0      /* value of the left point in the limit cycle */
#define RP      4      /* value of the right point in the limit cycle */
#define LK      1      /* value of the left knee in the limit cycle */
#define RK      3      /* value of the right knee in the limit cycle */
#define M	25600   	/* number of displayed oscillators = N1xN2 */
int intv[WIDTH];         /* step sizes */

