       /*****************************************************************/
       /* This program is for 2D segmentation using the LEGION network  */
       /* The equations were modeled as a PDP algorithm for efficiency  */
       /* by DeLiang Wang in 7/94. For a detailed description of the    */
       /* algorithm see D.L. Wang and D. Terman: "Image segmentation    */
       /* based on oscillatory correlation", Neural Computation, Vol. 9 */
       /* pp. 805-836, 1997 (For errata see Neural Computation 9,       */
       /* 1623-1626, 1997), on page 822.                                */
       /*****************************************************************/

#include "segment.h"

#define MAXVAL 32768.0
               /* 2147483648.0 or 32768.0,machine (HP or SUN) dependable*/
#define MAX 32767 /* 2147483647 or 32767, machine dependable */
#define GRAY 255  /* Maximum value of a pixel */

int bi();
int jump[2], z, Num = 1;             /* so we can change step size at will */
                               /* Num: number of the segments */
FILE *outfile, *fopen();
float x[N1][N2][2], lateral[N1][N2];
                      /* lateral: the sum of input from lateral connections */
float W[N1][N2][8];   /* strength of the interactions for 4 directions */
int stim[N1][N2], matrix[N1][N2];
int y[N1][N2], cycle[N1][N2], trigger[N1][N2];
     /* y: integer of x; cycle: test for 1st cycle; trigger: jumping signal */


inputs()
{
   int i, j;
   int width, height, c;
   char magic[25];
   float stimulus = 0.2;
   FILE *IN;

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++)
       matrix[i][j] = 0;

   IN = fopen("lake1.pgm", "r");
   fscanf(IN,"%s\n%d %d\n%d\n",magic,&width,&height,&c);
   fprintf(stderr,"[Type %s %dx%d colours %d]\n",magic,width,height,c);

   for(i=0;i<160;i++)  
      for(j=0;j<160;j++)
         stim[i][j] = getc(IN);

   fclose(IN);
}

connections()
{
   int i, j, k;
   float baseline = 0.0, Wtotal = 6.0, sum;

   for (i = 0; i < N1; i++) /* for eight nearest neighbors */
     for (j = 0; j < N2; j++) {
       W[i][j][0] = W[i][j][1] = W[i][j][2] = W[i][j][3] = 0.0;
       W[i][j][4] = W[i][j][5] = W[i][j][6] = W[i][j][7] = 0.0;
     }

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++) {
	 sum = 0.0;

         if (j-1 >= 0)
	    W[i][j][0] = (float) GRAY/(1+abs(stim[i][j]-stim[i][j-1]));
         if (i-1 >= 0 && j-1 >= 0)
	    W[i][j][1] = (float) GRAY/(1+abs(stim[i][j]-stim[i-1][j-1]));
         if (i-1 >= 0)
	    W[i][j][2] = (float) GRAY/(1+abs(stim[i][j]-stim[i-1][j]));
         if (i-1 >= 0 && j+1 < N2)
	    W[i][j][3] = (float) GRAY/(1+abs(stim[i][j]-stim[i-1][j+1]));
         if (j+1 < N2)
	    W[i][j][4] = (float) GRAY/(1+abs(stim[i][j]-stim[i][j+1]));
         if (i+1 < N1 && j+1 < N2) 
	    W[i][j][5] = (float) GRAY/(1+abs(stim[i][j]-stim[i+1][j+1]));
         if (i+1 < N1)
	    W[i][j][6] = (float) GRAY/(1+abs(stim[i][j]-stim[i+1][j]));
         if (i+1 < N1 && j-1 >= 0)
	    W[i][j][7] = (float) GRAY/(1+abs(stim[i][j]-stim[i+1][j-1]));

         for (k = 0; k < 8; k++) lateral[i][j] += W[i][j][k];
    }
}

initiate()
{
   int i,j,t, temp;  /* labels for setting up display intervals */

   for (i = 0; i < N1; i++)   /* randomize the initial phases */
     for (j = 0; j < N2; j++) {
        y[i][j] = rand();
        x[i][j][0] = y[i][j]/(MAXVAL + 1.0);
	y[i][j] -= 2;
	if (y[i][j] < 0)
	  y[i][j] = 0;
	trigger[i][j] = 0;
	cycle[i][j] = 0;
     }
   z = 0; 

 /* initial display values for 2D */
   for (i = 0; i < N1; i++)
      for (j = 0; j < N2; j++) {
        temp = x[i][j][0]*RP;       
        fprintf(outfile, "%d ", temp);
      }
}

select()
{
   int i, j, max, imax=0, jmax=0;
   float baseline = 1200.0, temp;
/* baseline 1200 -> 17 regions; 1000 -> 23 regions. Both reasonable */

   max = -1;
   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++) {
       trigger[i][j] = 0;
       if (lateral[i][j] > baseline && y[i][j] != MAX && y[i][j] > max) {
	 max = y[i][j];
	 imax = i;
	 jmax = j;
       }
     }

   if (cycle[imax][jmax]) {     /* find out jumpers */
     for (i = 0; i < N1; i++)
       for (j = 0; j < N2; j++)
	 if (y[i][j] == max && cycle[i][j])
	   trigger[i][j] = 1;
   } else
     trigger[imax][jmax] = 1;

   for (i = 0; i < N1; i++)   /* advance oscillators on the left branch */
     for (j = 0; j < N2; j++) 
       if (y[i][j] != MAX && (cycle[i][j] || lateral[i][j] > baseline)) {
	 y[i][j] = y[i][j] + (MAX - max);
         x[i][j][0] = y[i][j]/(MAXVAL + 1.0);
       }
}

evaluate()
{
   float temp, thetae = 2.0, temp2, temp3;
   int i,j,t, cat, count;  /* labels for setting up display intervals */

   jump[0] = 0;

   select();
   for ( t = 0; t < WIDTH-1; t++) {
      count = 0;
      for (i = 0; i < N1; i++)
	for (j = 0; j < N2; j++) { /* compute net input */
	  if (x[i][j][0] > thetae && jump[0] == 0) {  /* jumping down */
             x[i][j][1] = LP;
	     y[i][j] = -1;
	     cycle[i][j] = 1;
	     count++;
	  } else if (x[i][j][0] < thetae) {  /* summing up overall input */
             temp =  -20.0 * bi((float) z,0.5);  /* used to be 50 */

/* the following is based on traversing through strongest link */
             temp2 = 0.0;
             if (j-1 >= 0) {
                temp3 = W[i][j][0] * bi(x[i][j-1][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             if (i-1 >= 0 && j-1 >= 0) {
                temp3 = W[i][j][1] * bi(x[i-1][j-1][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             if (i-1 >= 0) {
                temp3 = W[i][j][2] * bi(x[i-1][j][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             if (i-1 >= 0 && j+1 < N2) {
                temp3 = W[i][j][3] * bi(x[i-1][j+1][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             if (j+1 < N2) {
                temp3 = W[i][j][4] * bi(x[i][j+1][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             if (i+1 < N1 && j+1 < N2) {
                temp3 = W[i][j][5] * bi(x[i+1][j+1][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             if (i+1 < N1) {
                temp3 = W[i][j][6] * bi(x[i+1][j][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             if (i+1 < N1 && j-1 >= 0) {
                temp3 = W[i][j][7] * bi(x[i+1][j-1][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             temp += temp2;

	     if (trigger[i][j] || temp > 0.25) {        /* jumping up */
	        x[i][j][1] = RP; 
		y[i][j] = MAX;
	        jump[1] = 1;
	        z++;
	     } else
	        x[i][j][1] = x[i][j][0];
	  } else 
	     x[i][j][1] = x[i][j][0];
	}
      z -= count;

      if (jump[1] == 0 && z == 0) 
        select();

    /* store for display in 2d case*/
      if (jump[1] == 0 && jump[0] == 1) {
printf("we are here \n");
	Num++;
        for (i = 0; i < N1; i++)
          for (j = 0; j < N2; j++) {
             if (z > 0) {
	          if (x[i][j][1] < thetae)
                     cat = 0;       
		  else cat = x[i][j][1];       
	     } else  
		  cat = x[i][j][1];       
             fprintf(outfile, "%d ", cat);

             if (Num < 25 && cat > 0)
                matrix[N1-1-i][j] = 26-Num;  /* 17 regions, or 23 */
	  }
      }


      for (i = 0; i < N1; i++)
         for (j = 0; j < N2; j++)  /* replace OLD by NEW */
           x[i][j][0] = x[i][j][1];
      jump[0] = jump[1];
      jump[1] = 0;
   }
}

store()
{
   int i, j;
   FILE *out, *out1, *fopen();

   printf("The number of the segments/cycles is %d\n", Num);
   
   out = fopen("NUMBER", "w");
   out1 = fopen("Result", "w");
   fprintf(out, "%d ", Num);
   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++)
       fprintf(out1, "%d ", matrix[i][j]);
   fclose(out);
   fclose(out1);
}

main()
{
    int seed;

    printf("Input a seed now\n");
    scanf("%d",&seed); 
    srand(seed);

    outfile = fopen("Out", "w");

    inputs();
    connections();
    initiate();
    evaluate();
    store();
    fclose(outfile);
}
   
bi(x, threshold)
   float x, threshold;
{
    if (x > threshold) return (1);
    else return (0);
}
