       /*****************************************************************/
       /* This program is for 2D image segmentation using a LEGION      */
       /* algorithm (no dynamics) described in IEEE Trans. on Neural    */
       /* Networks, vol. 11, No. 5, pp. 1106-1123, 2000, in an article  */
       /* entitled "Weight adaptation and oscillatory correlation       */
       /* for image segmentation". Re-implemented by DeLiang Wang       */
       /* after the paper was published. May 2001.                      */
       /*****************************************************************/

#include <stdio.h>
#include <math.h>
#define N1      230     /* number of rows: height */
#define N2      240     /* number of columns: width */
#define GRAY 255  /* Maximum value of a pixel */

   /* Adjustble parameters: 1st 3 for leader selection; Wz for segmentation */
#define Rp 9            /* reasonable for runway extraction */
#define Tmu 2.0         /* 0.1 in paper */
#define Tsigma 10.0     /* 2.0 in paper */
#define Wz 30.0         /* 40.0 for "noise_40" */

int result[N1][N2], Num;            
 /* result: segmentation results; Num: number of the segments */
float x[N1][N2], W[N1][N2][8]; /* W: weights for 8 neighbors */
int I[N1][N2]; /* Input image */
int leader[N1][N2], z; /* leader: selected seeds; z: global inhibitor */

inputs() /* input an image to segment */
{
   int i, j;
   int width, height, c;
   char magic[25];
   FILE *IN, *fopen();

   IN = fopen("result/noise_1000.pgm", "r"); /* input smoothed version of original */
   fscanf(IN,"%s\n%d %d\n%d\n",magic,&width,&height,&c);
   fprintf(stderr,"[Type %s %dx%d colours %d]\n",magic,width,height,c);

   for(i=0;i<N1;i++)  
     for(j=0;j<N2;j++) /* read a text file */
         fscanf(IN,"%d",&I[i][j]);
   fclose(IN);
}

connections() /* set up weights */
{
   int i, j, Imin = GRAY, Imax = 0, Wmax;

   for (i = 0; i < N1; i++) /* for eight nearest neighbors */
     for (j = 0; j < N2; j++) {
       if (I[i][j] < Imin)
	 Imin = I[i][j];
       if (I[i][j] > Imax)
	 Imax = I[i][j];
       W[i][j][0] = W[i][j][1] = W[i][j][2] = W[i][j][3] = 0.0;
       W[i][j][4] = W[i][j][5] = W[i][j][6] = W[i][j][7] = 0.0;
     }

   Wmax = Imax - Imin;
   fprintf(stderr,"Wmax (Imax - Imin) = %d\n", Wmax);

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++) {
         if (j-1 >= 0)
	    W[i][j][0] = (float) Wmax/(1+abs(I[i][j]-I[i][j-1]));
         if (i-1 >= 0 && j-1 >= 0)
	    W[i][j][1] = (float) Wmax/(1+abs(I[i][j]-I[i-1][j-1]));
         if (i-1 >= 0)
	    W[i][j][2] = (float) Wmax/(1+abs(I[i][j]-I[i-1][j]));
         if (i-1 >= 0 && j+1 < N2)
	    W[i][j][3] = (float) Wmax/(1+abs(I[i][j]-I[i-1][j+1]));
         if (j+1 < N2)
	    W[i][j][4] = (float) Wmax/(1+abs(I[i][j]-I[i][j+1]));
         if (i+1 < N1 && j+1 < N2) 
	    W[i][j][5] = (float) Wmax/(1+abs(I[i][j]-I[i+1][j+1]));
         if (i+1 < N1)
	    W[i][j][6] = (float) Wmax/(1+abs(I[i][j]-I[i+1][j]));
         if (i+1 < N1 && j-1 >= 0)
	    W[i][j][7] = (float) Wmax/(1+abs(I[i][j]-I[i+1][j-1]));
    }
}

initiate()
{
   int i,j; 

   Num = 0;
   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++)
       result[i][j] = 0;
}

select() /* leader or seed selection */
{
   int i, j, m, n, winsize;
   float temp, mean1, meanR, var1, varR;
   FILE *out, *fopen();

   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++) {
       leader[i][j] = 0;

       winsize = 0; mean1 = 0.0; var1 = 0.0;
       for (m = -1; m < 2; m++)
         for (n = -1; n < 2; n++)
           if ((i+m >=0) && (j+n>=0) && (i+m < N1) && (j+n < N2)) {
               mean1 += I[i+m][j+n];
               var1 += I[i+m][j+n] * I[i+m][j+n];
               winsize++;
           }
       mean1 = mean1/winsize;
       var1 = var1/winsize - mean1*mean1;

       winsize = 0; meanR = 0.0; varR = 0.0;
       for (m = -Rp; m < Rp + 1; m++)
         for (n = -Rp; n < Rp + 1; n++)
           if ((i+m >=0) && (j+n>=0) && (i+m < N1) && (j+n < N2)) {
               meanR += I[i+m][j+n];
               varR += I[i+m][j+n] * I[i+m][j+n];
               winsize++;
           }
       meanR = meanR/winsize;
       varR = varR/winsize - meanR*meanR;

       if (abs(meanR-mean1) <= Tmu && abs(varR-var1) <= Tsigma)
	 leader[i][j] = 1;
     }

   /* print the leaders image */
   out = fopen("Leader.pgm", "w");
   fprintf(out,"%s\n%d %d\n%d\n","P2",N2,N1,1);
   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++)
       fprintf(out, "%d ", leader[i][j]);
   fclose(out);
}

evaluate()
{
   float Sij;
   int i, j, m, n, nonstop = 1, ilead, jlead, jump, count;  

   while (nonstop) {
     nonstop = 0;
     z = 0;
     ilead = -1; jlead = -1; 
     for (i = N1-1; i >= 0; i--) /* leaders are activated systematically */
       for (j = N2-1; j >= 0; j--) {
	 x[i][j] = 0;
	 if (leader[i][j]) {
	   ilead = i;
	   jlead = j;
	 }
       }

     if (ilead>=0 && jlead>=0) {
       x[ilead][jlead] = 1;
       nonstop = 1;
       Num++;
       z = 1;
       jump = 1;

       while (jump) {
         jump = 0;
         for (i = 0; i < N1; i++) /* jump more oscillators */
           for (j = 0; j < N2; j++) 
	     if (result[i][j] == 0 && x[i][j] == 0) { 
	                                /* an oscillator doesn't jump twice */ 
  	       count = 0;
               for (m = -1; m < 2; m++)
                 for (n = -1; n< 2; n++)
                   if ((i+m >=0) && (j+n>=0) && (i+m < N1) && (j+n < N2) && 
		       x[i+m][j+n] && (m != 0 || n != 0))
                      count++;

	       /* calculate S_ij */
	       Sij = 0.0;
               if (j-1 >= 0)
	          Sij += W[i][j][0]*x[i][j-1];
               if (i-1 >= 0 && j-1 >= 0)
	          Sij += W[i][j][1]*x[i-1][j-1];
               if (i-1 >= 0)
	          Sij += W[i][j][2]*x[i-1][j];
               if (i-1 >= 0 && j+1 < N2)
	          Sij += W[i][j][3]*x[i-1][j+1];
               if (j+1 < N2)
	          Sij += W[i][j][4]*x[i][j+1];
               if (i+1 < N1 && j+1 < N2) 
	          Sij += W[i][j][5]*x[i+1][j+1];
               if (i+1 < N1)
	          Sij += W[i][j][6]*x[i+1][j];
               if (i+1 < N1 && j-1 >= 0)
	          Sij += W[i][j][7]*x[i+1][j-1];

	       if (count > 0)
	         Sij = Sij/log(count+1);

  	       if (Sij > Wz*z) { 
	         x[i][j] = 1;
	         jump = 1;
	       }
	     }
       }
	
       /* record a new segment & erase leaders */
       for (i = 0; i < N1; i++)
         for (j = 0; j < N2; j++)
	   if (x[i][j]) {
	     result[i][j] = Num;
	     leader[i][j] = 0;
	   }
     }
   }
}

store()
{
   int i, j;
   FILE *out, *fopen();

   printf("The number of segments is %d\n", Num);
   
   out = fopen("Result.pgm", "w");
   fprintf(out,"%s\n%d %d\n%d\n","P2",N2,N1,Num);
   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++)
       fprintf(out, "%d ", result[i][j]);
   fclose(out);
}

main()
{
   int seed;

   inputs();
   connections();
   initiate();
   select();
   evaluate();
   store();
}
