       /*****************************************************************/
       /* This program is for 2D segmentation using Terman-Wang         */
       /* oscillators. The program implements the selective gating      */
       /* mechanism for LEGION networks. For a detailed description of  */
       /* the equations see D.L. Wang and D. Terman (1995): "Locally    */
       /* excitatory globally inhibitory oscillator networks,"          */
       /* IEEE Transactions on Neural Networks, vol. 6(1), 283-286.     */
       /* Implemented by D.L. Wang, 8/93.                               */
       /*****************************************************************/

#include "terman.h"

#define MAXVAL 32768.0 /* 2147483648.0, or 32768.0, machine dependable */
#define tnh(x) ((exp(x) - exp(-x))/(exp(x) + exp(-x)))
#define sgm(x,theta) (1.0/(1.0 + exp(-(x-theta)*50.0)))
#define H(x) 1/(1 + exp((-x)))       /* frequently used function */
#define h   0.2            /* delta t */

double eta;           /* value of the slow paramenter */
double rho;           /* controls the noise in the system         */
double Gamma;         /* the parameter for the inhibitory unit    */
double f(),g(),rnd(),noise(),s();
int bi();                 /* so we can change step size at will */

double x[N1][N2][2],y[N1][N2][2],ifx[2],ify[2],z[2],show[M][WIDTH2];
double W[N1][N2][4];   /* strength of the interactions for 4 directions */
double input[N1][N2];


inputs()
{
   int i, j;
   float stim = 0.2;

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++)
       input[i][j] = -0.02;

   for (i = 2; i < N1-1; i++) input[i][0] = stim;
   for (i = 2; i < N1-1; i++) input[i][3] = stim;
   for (i = 1; i < N1-1; i++) input[i][5] = stim;
   for (i = 1; i < N1-1; i++) input[i][8] = stim;
   for (i = 3; i < N1-3; i++) input[i][12] = stim;
   for (i = 1; i < N1-2; i++) input[i][16] = stim;
   for (i = 1; i < N1-2; i++) input[i][N2-1] = stim;

 /* for the first O */
   input[2][1] = stim; input[2][2] = stim;
   input[3][1] = stim; input[3][2] = stim;
   input[N1-3][1] = stim; input[N1-3][2] = stim;
   input[N1-2][1] = stim; input[N1-2][2] = stim;

 /* for H */
   input[9][6] = stim; input[9][7] = stim;
   input[10][6] = stim; input[10][7] = stim;

 /* for I */
   input[3][11] = stim; input[3][13] = stim;
   for (j = 10; j < 15; j++) input[16][j] = stim;

 /* for the second O */
   input[1][N2-3] = stim; input[1][N2-2] = stim;
   input[2][N2-3] = stim; input[2][N2-2] = stim;
   input[N1-4][N2-3] = stim; input[N1-4][N2-2] = stim;
   input[N1-3][N2-3] = stim; input[N1-3][N2-2] = stim;

  printf("inputs are \n");
  for (i = 0; i < N1; i++) {
    for (j = 0; j < N2; j++) {
      if (input[i][j] > stim/2)
	printf("%s ", "*");
      else
	printf("%s ", " ");
    }
    printf("\n");
  }
}

connections()
{
   int i, j, count;
   float baseline = 0.05, Wtotal = 6.0;

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++) {
       W[i][j][0] = W[i][j][1] = W[i][j][2] = W[i][j][3] = 0.0;
     }

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++) 
       if (input[i][j] > baseline) {   /* oscillator(i,j) has input */
         count = 0;
         if (i-1 >= 0 && input[i-1][j] > baseline) count++;
	 if (i+1 < N1 && input[i+1][j] > baseline) count++;
	 if (j-1 >= 0 && input[i][j-1] > baseline) count++;
	 if (j+1 < N2 && input[i][j+1] > baseline) count++;
	 if (count) {
           if (i-1 >= 0 && input[i-1][j] > baseline) W[i][j][1] = Wtotal/count;
	   if (i+1 < N1 && input[i+1][j] > baseline) W[i][j][3] = Wtotal/count;
	   if (j-1 >= 0 && input[i][j-1] > baseline) W[i][j][0] = Wtotal/count;
	   if (j+1 < N2 && input[i][j+1] > baseline) W[i][j][2] = Wtotal/count;

	 }
       }
  printf("weights are \n");
  printf("W[0][0] are %f %f %f %f\n", W[0][0][0], W[0][0][1], W[0][0][2],
	                    W[0][0][3]);
  printf("W[1][0] are %f %f %f %f\n", W[1][0][0], W[1][0][1], W[1][0][2],
	                    W[1][0][3]);
  printf("W[2][0] are %f %f %f %f\n", W[2][0][0], W[2][0][1], W[2][0][2],
	                    W[2][0][3]);
  printf("W[2][1] are %f %f %f %f\n", W[2][1][0], W[2][1][1], W[2][1][2],
	                    W[2][1][3]);
  printf("\n");
}

runge_kutta(x1,y1,x2,y2,input)
   double x1, y1, *x2, *y2, input;
{
   double kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4;

   kx1 = h*f(x1,y1,input); 
   ky1 = h*g(x1,y1);

   kx2 = h*f(x1+kx1/2.0, y1+ky1/2.0, input);
   ky2 = h*g(x1+kx1/2.0, y1+ky1/2.0); 

   kx3 = h*f(x1+kx2/2.0, y1+ky2/2.0, input);
   ky3 = h*g(x1+kx2/2.0, y1+ky2/2.0);

   kx4 = h*f(x1+kx3, y1+ky3, input);
   ky4 = h*g(x1+kx3, y1+ky3);

   *x2 = x1 + (kx1 + 2.0*kx2 + 2.0*kx3 + kx4)/6.0;   
   *y2 = y1 + (ky1 + 2.0*ky2 + 2.0*ky3 + ky4)/6.0; 
}

evaluate()
{
   float ki1, ki2, ki3, ki4, xinput[N1][N2], zinput;
   float temp, thetae = -0.5;
   int i,j,t,hop, flag=1;
   FILE *outfile, *fopen();

   for (i = 0; i < N1; i++)   /* randomize the initial phases */
     for (j = 0; j < N2; j++) {
        x[i][j][0] = 2.0*rnd();
        y[i][j][0] = Gamma*(rnd() + 1.0);   
    }
    z[0] = 0.0;
     
    for ( t = 0; t < WIDTH; t++) {
      for (i = 0; i < N1; i++)
	for (j = 0; j < N2; j++) { /* compute net input */
          xinput[i][j] = input[i][j] + noise() - 2.0 * sgm(z[0],0.1);  
	  if (i-1 >= 0)
	     xinput[i][j] += W[i][j][1] * sgm(x[i-1][j][0], thetae);
	  if (i+1 < N1)
	     xinput[i][j] += W[i][j][3] * sgm(x[i+1][j][0], thetae);
	  if (j-1 >= 0)
	     xinput[i][j] += W[i][j][0] * sgm(x[i][j-1][0], thetae);
	  if (j+1 < N2)
	     xinput[i][j] += W[i][j][2] * sgm(x[i][j+1][0], thetae);
        }

      for (i = 0; i < N1; i++)
	for (j = 0; j < N2; j++) /* differentiate */
          runge_kutta(x[i][j][0],y[i][j][0],&x[i][j][1],&y[i][j][1],
		                                              xinput[i][j]);

/* compute global inhibition */
      temp = 0.0;
      for (i = 0; i < N1; i++)
	for (j = 0; j < N2; j++) /* differentiate */
	  temp += sgm(x[i][j][0],0.1);  

        zinput = bi(temp, 0.1);  /* should be sgm */
        ki1 = h*s(z[0], zinput); 
        ki2 = h*s(z[0]+ki1/2.0, zinput);
        ki3 = h*s(z[0]+ki2/2.0, zinput);
        ki4 = h*s(z[0]+ki3, zinput);
        z[1] = z[0] + (ki1 + 2.0*ki2 + 2.0*ki3 + ki4)/6.0;

        hop = WIDTH/WIDTH2;

    /* store for display in 1d case*/
/*
        if (t%hop == 0) {   
	    for (j = 0; j < N2; j++) {
               show[j][t/hop] = x[0][j][0] + 2.0;       
               show[N2+j][t/hop] = x[16][j][0] + 2.0;      
	    }
          show[M-1][t/hop] = z[0];
        }
*/

    /* store for display in 2d case*/

        if (t%hop == 0) {   
          for (i = 0; i < N1; i++)
	    for (j = 0; j < N2; j++)
	      if (i*N2+j < M-1)
                show[i*N2+j][t/hop] = x[i][j][0] + 2.0;       
          show[M-1][t/hop] = z[0];
        }


      for (i = 0; i < N1; i++)
	for (j = 0; j < N2; j++) { /* replace OLD by NEW */
          x[i][j][0] = x[i][j][1];
          y[i][j][0] = y[i][j][1];
	}
      z[0] = z[1];
    }
}

store()
{
   int i,j,t;
   float max_v = 0.0, min_v = 0.0, baseline = 0.05;
   FILE *outfile,*out0, *out1, *out2, *out3, *out4, *fopen();

   outfile = fopen("Output", "w");

   out0 = fopen("Snapshot0", "w");
   out1 = fopen("Snapshot1", "w");
   out2 = fopen("Snapshot2", "w");
   out3 = fopen("Snapshot3", "w");
   out4 = fopen("Snapshot4", "w");

   for (i = 0; i < M-1; i++) 
      for (t = 0; t < WIDTH2; t++)  {
        if (show[i][t] > max_v)
          max_v = show[i][t];
	if (show[i][t] < min_v)
          min_v = show[i][t];
      }
	
printf("Maximum is %f,   Minimum is %f\n", max_v, min_v);

   for (i = 0; i < M; i++)
     for (t = 0; t < WIDTH2; t++)
       fprintf(outfile, "%f ", (show[i][t]-min_v)/(max_v-min_v));


   for (i = 0; i < M-1; i++) {
     fprintf(out0, "%f ", (show[i][0]-min_v)/(max_v-min_v));
     fprintf(out1, "%f ", (show[i][385]-min_v)/(max_v-min_v));
     fprintf(out2, "%f ", (show[i][408]-min_v)/(max_v-min_v));
     fprintf(out3, "%f ", (show[i][430]-min_v)/(max_v-min_v));
     fprintf(out4, "%f ", (show[i][450]-min_v)/(max_v-min_v));
   }

   out0 = fopen("image_O1", "w");
   out1 = fopen("image_H", "w");
   out2 = fopen("image_I", "w");
   out3 = fopen("image_O2", "w");
   out4 = fopen("inhibitor", "w");

   for (i = 0; i < N1; i++) /* document the first O */
     for (j = 0; j < 4; j++)
       if (input[i][j] > baseline)
         for (t = 0; t < WIDTH2; t++)
	   fprintf(out0, "%f ", (show[i*N2+j][t]-min_v)/(max_v-min_v));

   for (i = 0; i < N1; i++) /* document H */
     for (j = 5; j < 9; j++)
       if (input[i][j] > baseline)
         for (t = 0; t < WIDTH2; t++)
	   fprintf(out1, "%f ", (show[i*N2+j][t]-min_v)/(max_v-min_v));

   for (i = 0; i < N1; i++)  /* document I */
     for (j = 10; j < 15; j++)
       if (input[i][j] > baseline)
         for (t = 0; t < WIDTH2; t++)
	   fprintf(out2, "%f ", (show[i*N2+j][t]-min_v)/(max_v-min_v));

   for (i = 0; i < N1; i++)  /* document the second O */
     for (j = 16; j < N2; j++)
       if (input[i][j] > baseline)
         for (t = 0; t < WIDTH2; t++)
	   fprintf(out3, "%f ", (show[i*N2+j][t]-min_v)/(max_v-min_v));

   for (t = 0; t < WIDTH2; t++)
     fprintf(out4, "%f ", (show[M-1][t]-min_v)/(max_v-min_v));

   fclose(outfile);
}

main()
{
    int i,j,seed,n;

    printf("Input a seed now\n");
    scanf("%d",&seed); 
    srand(seed); /* 213 for dynamic normalization; 212121 for the other */

    eta = 0.02;
    rho = 0.02;
    Gamma = 6.0;

    inputs();
    connections();
    evaluate();
    store();
}
 
double noise()
{
   double x1,x2;
   x1 = 0.5+rnd()/2.0;
   x2 = 0.5+rnd()/2.0;
   return(rho*(4*x1*(1-x2)-0.5));
}

double rnd()
{
   return(1 - 2*(rand()/(MAXVAL - 1.0)));
}

double f(x,y,input)
double x,y,input;
{
  return((3*x - x*x*x + 2 - y + input));
}

double g(x,y)
double x,y;
{  
   return(eta*(Gamma*(1 + tnh(x/0.1)) - y));
}
  
int bi(x, threshold)
float x, threshold;
{
    if (x > threshold) return (1);
    else return (0);
}

double s(x,input)
double x,input;
{
   return (3.0*(input - x));
}
