       /*****************************************************************/
       /* This is the first header file for the program that models     */
       /* 2D segmentation. By DeLiang Wang, 8/93, at OSU                */
       /*****************************************************************/
       
#include <stdio.h>
#include <math.h>
#define WIDTH   1600	/* integration time */
#define N1	127	/* number of rows */
#define N2      150       /* number of columns */
#define KS      5      /* the diameter of a kernel */
#define LP      0      /* value of the left point in the limit cycle */
#define RP      4      /* value of the right point in the limit cycle */
#define LK      1      /* value of the left knee in the limit cycle */
#define RK      3      /* value of the right knee in the limit cycle */
#define M	19050   	/* number of displayed oscillators = N1xN2 */
int intv[WIDTH];         /* step sizes */

