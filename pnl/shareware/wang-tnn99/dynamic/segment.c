       /*****************************************************************/
       /* This program is for 2D segmentation using oscillators         */
       /* based on Terman-Wang equations.  These equations were modeled */
       /* as a PDP algorithm for efficiency by DeLiang Wang, 7/94       */
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
float W[N1][N2][4];   /* strength of the interactions for 4 directions */
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

   IN = fopen("v8n6.pitchmask.pgm", "r");
   fscanf(IN,"%s\n%d %d\n%d\n",magic,&width,&height,&c);
   fprintf(stderr,"[Type %s %dx%d colours %d]\n",magic,width,height,c);

   for(j=0;j<N2;j++)
     for(i=0;i<N1;i++) {
        fscanf(IN, "%d",&stim[i][j]); 
        if (stim[i][j] > 0)
          stim[i][j] = 1;
     }

   fclose(IN);
}

connections()
{
   int i, j, k;
   float baseline = 0.0, Wtotal = 6.0, sum;

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++) {
       W[i][j][0] = W[i][j][1] = W[i][j][2] = W[i][j][3] = 0.0;
     }

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++) {
	 sum = 0.0;

/* the following is for absolute difference */
         if (j-1 >= 0 && stim[i][j] == stim[i][j-1])
	    W[i][j][0] = (float) GRAY;
         if (i-1 >= 0 && stim[i][j] == stim[i-1][j])
	    W[i][j][1] = (float) GRAY;
         if (j+1 < N2 && stim[i][j] == stim[i][j+1])
	    W[i][j][2] = (float) GRAY;
         if (i+1 < N1 && stim[i][j] == stim[i+1][j])
	    W[i][j][3] = (float) GRAY;

	 lateral[i][j] = 0.0;
	 if (stim[i][j] > 0) {
	   if (j-1 >= 0 && stim[i][j-1] > 0)
	     lateral[i][j] += 1.0;
	   if (j+1 < N2 && stim[i][j+1] > 0)
	     lateral[i][j] += 1.0;
	 }
    }

  printf("weights are \n");
  printf("W[0][0] are %.2f %.2f %.2f %.2f\n", 
	 W[0][0][0], W[0][0][1], W[0][0][2],W[0][0][3]);
  printf("W[1][1] are %.2f %.2f %.2f %.2f\n", 
	 W[1][1][0], W[1][1][1], W[1][1][2],W[1][1][3]);
  printf("W[14][85] are %.2f %.2f %.2f %.2f\n", 
	 W[14][85][0], W[14][85][1], W[14][85][2],W[14][85][3]);
  printf("\n");
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
   float baseline = 1.5, temp;

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
     printf("repeating now and Num = %d\n", Num);
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

/* -20 is good for sqrt(pixel) */
/* -20 is good for absolute distance */
/* the following is based on traversing through strongest link */

             temp2 = 0.0;
             if (j-1 >= 0) {
                temp3 = W[i][j][0] * bi(x[i][j-1][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             if (i-1 >= 0) {
                temp3 = W[i][j][1] * bi(x[i-1][j][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             if (j+1 < N2) {
                temp3 = W[i][j][2] * bi(x[i][j+1][0], thetae);
                if (temp2 < temp3) temp2 = temp3; }
             if (i+1 < N1) {
                temp3 = W[i][j][3] * bi(x[i+1][j][0], thetae);
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
printf("we are here. Num = %d\n", Num);
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

             if (Num < 95 && cat > 0)
                matrix[N1-1-i][j] = 96-Num;  /* 90 segments for phone; 94 for v8n6 */
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
