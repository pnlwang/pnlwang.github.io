       /*****************************************************************/
       /* This program is for 2D segmentation using oscillators         */
       /* based on Terman-Wang equations.  These equations are solved   */
       /* in the sigular limit, as described by Linsay P.S. and         */
       /* Wang D.L. (1998): "Fast numerical integration of relaxation   */
       /* oscillator networks based on singular limit solutions," IEEE  */
       /* Transactions on Neural Networks, vol. 9, 523-532.             */
       /* Written by DeLiang Wang, 7/96, at OSU.                        */
       /*****************************************************************/

#include "singular.h"

#define MAXVAL 32768.0
               /* 2147483648.0 or 32768.0,machine (HP or SUN) dependable*/
#define MAX 32767 /* 2147483647 or 32767, machine dependable */

float epsilon = 0.02;           /* value of the slow parameter */
float Gamma = 6.5;         /* CAREFULL: this value offset Wz = 1.5 */
float alpha = 5.0;
float theta = 0.3;
float lambda = 0.1;        /* rate of developing potential */
float time = 0.0;          /* global time */
float rnd(), getx(), getxx();
int z[2], Hp[N1][N2];        /* argument to the Heaviside for potential */

int image1[N1][N2], image2[N1][N2], image3[N1][N2];
float y[N1][N2],cubic[N1][N2],p[N1][N2],show[M][STEPS];
float W[N1][N2][4];   /* strength of the interactions for 4 directions */
float input[N1][N2];  /* overall input an oscillator receives */
int branch[N1][N2];
     /* branch: keeps track of whether an oscillator is on LB or RB */

inputs()
{
   FILE *fopen(), *finput;
   int i, j, temp;
   float stim = 0.2;

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++)
       input[i][j] = baseline;

/* The input is "Sun-Mountain-Tree" */
   finput = fopen("Stim", "r");
   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++) {
       fscanf(finput, "%d", &temp);
       if (temp)
         input[i][j] = stim;
     }

 /* add noise */
   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++) 
       if (rand()/(MAXVAL - 1.0) < 0.2)  /* 20% noise */
           input[i][j] = stim;


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
   float Wtotal = 8.0;

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++) {
       W[i][j][0] = W[i][j][1] = W[i][j][2] = W[i][j][3] = 0.0;
     }

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++) 
       if (input[i][j] > 0.0) {   /* oscillator(i,j) has input */
         count = 0;
         if (i-1 >= 0 && input[i-1][j] > 0.0) count++;
	 if (i+1 < N1 && input[i+1][j] > 0.0) count++;
	 if (j-1 >= 0 && input[i][j-1] > 0.0) count++;
	 if (j+1 < N2 && input[i][j+1] > 0.0) count++;
	 if (count) {
           if (i-1 >= 0 && input[i-1][j] > 0.0) W[i][j][1] = Wtotal/count;
	   if (i+1 < N1 && input[i+1][j] > 0.0) W[i][j][3] = Wtotal/count;
	   if (j-1 >= 0 && input[i][j-1] > 0.0) W[i][j][0] = Wtotal/count;
	   if (j+1 < N2 && input[i][j+1] > 0.0) W[i][j][2] = Wtotal/count;

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

initiate()
{
   int i,j,seed;  /* labels for setting up display intervals */
   float angle;

   printf("Input another seed now\n"); /* 2342 for potential; 111 no potent.*/
   scanf("%d",&seed); 
   srand(seed); 

   for (i = 0; i < N1; i++)  /* place the initial phases on the left branch */
     for (j = 0; j < N2; j++) {
       if (input[i][j] > 0.0)
         y[i][j] = 2.0*Gamma*rnd()+input[i][j];
       branch[i][j] = 0;
       /*
       p[i][j] = 0.0;
       Hp[i][j] = 0;
       */
       cubic[i][j] = input[i][j];
     }

   for (i = 0; i < N1; i++)  /* place the initial phases on the left branch */
     for (j = 0; j < N2; j++)
       if (input[i][j] < 0.0)
	 y[i][j] = 2.0*Gamma*rnd();

   z[0] = 0; 

 /* initial display values for 1D */
/*
   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++)
       show[N2*i+j][0] = getxx(y[i][j]-cubic[i][j], 0);
*/

 /* initial display values for 2D */
/*
   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++)
       show[N2*i+j][0] = getxx(y[i][j]-cubic[i][j], 0);

   show[M-1][0] = 0.0;
*/
}

findtmin(tmin, mu, ilead, jlead) 
  float *tmin, *mu;
  int *ilead, *jlead;
     /* Find shortest time-to-knee when no oscillator is ready to jump */
{
   int i, j, il, jl, ir, jr;
   float temp, minl = MAXVAL, minr = MAXVAL, yknee;

   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++) {
       if (branch[i][j]) {          /* on RB */
	 yknee = (y[i][j]-2.0*Gamma)/(cubic[i][j]+4.0-2.0*Gamma);
	 if (yknee < minr && yknee > 0.0)  /* exclude stable fixed point */
	   {minr = yknee; ir = i; jr = j;}
       }
       else {                       /* on LB */
	 yknee = y[i][j]/cubic[i][j];
	 if (yknee < minl && yknee > 0.0)  /* exclude stable fixed point */
  	   {minl = yknee; il = i; jl = j;}
       }
     }

   if (minl < 1.0)
     {*tmin = 0.0; *mu = 1.0; *ilead = il; *jlead = jl;}
   else if (minr < 1.0)
     {*tmin = 0.0; *mu = 1.0; *ilead = ir; *jlead = jr;}
   else {
     *tmin = log(minl);
     *mu = minl;
     *ilead = il;
     *jlead = jl;
     temp = log(minr);
     if (temp < *tmin) {       
       *tmin = temp;
       *mu = minr;
       *ilead = ir; *jlead = jr;
     } 
   }
}

jump()      /* activity at the jump (up or down) time */
{       
   int i, j, jumped = 1, Tij = 2, flag = 1;  /* Tij: permanent weight */

   while (jumped) {
      jumped = 0;
      z[1] = z[0];
      for (i = 0; i < N1; i++)
	for (j = 0; j < N2; j++) { /* compute net input */
	  if (time < 7.0 || p[i][j] > theta || flag) { 
	                                      /* thrshld for exp(-alpha*t) */
  	    if (z[0] > 0)
              cubic[i][j] = input[i][j] - 1.5;  /* global inhibition */
	    else
	      cubic[i][j] = input[i][j];
          } else {
	    if (z[0] > 0)
	      cubic[i][j] = baseline - 1.5;
	    else
	      cubic[i][j] = baseline;
	  }
          if (i-1 >= 0 && branch[i-1][j])
	     cubic[i][j] += W[i][j][1];
	  if (i+1 < N1 && branch[i+1][j])
	     cubic[i][j] += W[i][j][3];
  	  if (j-1 >= 0 && branch[i][j-1])
	     cubic[i][j] += W[i][j][0];
	  if (j+1 < N2 && branch[i][j+1])
	     cubic[i][j] += W[i][j][2];

	  if (!branch[i][j] && y[i][j]<cubic[i][j]+error) {  /* jump up */
		branch[i][j] = 1;
		jumped = 1;
		z[1]++;
	  } else if (branch[i][j] && y[i][j]>cubic[i][j]+4.0-error) {
	                                                     /* jump down */
		branch[i][j] = 0;
		jumped = 1;
		z[1]--;
	  }
	}
      z[0] = z[1];
   }

   for (i = 0; i < N1; i++)  /* compute Hp[i][j] */
     for (j = 0; j < N2; j++) { 
       Hp[i][j] = 0;
       if (i-1 >= 0 && branch[i-1][j])
         Hp[i][j] += Tij;
       if (i+1 < N1 && branch[i+1][j])
	 Hp[i][j] += Tij;
       if (j-1 >= 0 && branch[i][j-1])
	 Hp[i][j] += Tij;
       if (j+1 < N2 && branch[i][j+1])
	 Hp[i][j] += Tij;

       if (Hp[i][j] > 7)
	 Hp[i][j] = 1;
       else 
	 Hp[i][j] = 0;
     }

}

evaluate()
{
   float intv, mu, current, length, xval, yval, angle, temp, cat;
   int ilead, jlead, num; /* leading oscillator, and no of jumpers */
   int i,j,step,bef,aft;  /* labels for setting up display intervals */

   bef = 0;
   length = delta*(STEPS-1);      /* reduce one step to facilitate display */
   while (time < length) {
     findtmin(&intv, &mu, &ilead, &jlead);

     aft = (time+intv)/delta;
     
    /* store for display in 1d case*/
/*
     for (step = bef; step < aft && step < STEPS-1; step++) {
       for (i = 0; i < N1; i++)
         for (j = 0; j < N2; j++) {
           temp = -((1.0+step)*delta - time);
	   if (branch[i][j])
	     yval = (y[i][j] - 2*Gamma)*exp(temp) + 2*Gamma - cubic[i][j];
	   else
	     yval = y[i][j]*exp(temp) - cubic[i][j];

           show[N2*i+j][step+1] = getxx(yval,branch[i][j]);
         }

       if (z[0] > 0)
         show[M-1][step+1] = 1.0-2.5;
       else
         show[M-1][step+1] = 0.0-2.5;
     }
*/

    /* store for display in 2d case*/

     for (step = bef; step < aft && step < STEPS-1; step++) {
       for (i = 0; i < N1; i++)
         for (j = 0; j < N2; j++) {
           temp = -((1.0+step)*delta - time);
	   if (branch[i][j])
	     yval = (y[i][j] - 2*Gamma)*exp(temp) + 2*Gamma - cubic[i][j];
	   else
	     yval = y[i][j]*exp(temp) - cubic[i][j];

           show[N2*i+j][step+1] = getxx(yval,branch[i][j]);
         }

       if (z[0] > 0)
         show[M-1][step+1] = 1.0-2.5;
       else
         show[M-1][step+1] = 0.0-2.5;

    /* store into images for later scanning for output display */
/*
	if (step == 255)
	  for (i = 0; i < N1; i++)
	    for (j = 0; j < N2; j++) 
		image1[i][j] = branch[i][j];
	if (step == 275)
	  for (i = 0; i < N1; i++)
	    for (j = 0; j < N2; j++) 
		image2[i][j] = branch[i][j];
	if (step == 295)
	  for (i = 0; i < N1; i++)
	    for (j = 0; j < N2; j++) 
		image3[i][j] = branch[i][j];
*/
     }


     for (i = 0; i < N1; i++) /* advance oscillators */
       for (j = 0; j < N2; j++) {
	 if (branch[i][j])
	   y[i][j] = (y[i][j]-2.0*Gamma)/mu + 2.0*Gamma;
	 else
	   y[i][j] = y[i][j]/mu;
       }


     for (i = 0; i < N1; i++) /* update oscillator potential */
       for (j = 0; j < N2; j++)
	 if (Hp[i][j] && input[i][j] > 0.0)
	   p[i][j] = 1.0-(1.0-p[i][j])/mu;


     if (branch[ilead][jlead]) {
       branch[ilead][jlead] = 0;
       z[0]--;
     } else {
       branch[ilead][jlead] = 1;
       z[0]++;
     }

     jump();            /* fast time activity */

     time += intv;      
     bef = aft;
   }
}

store()
{
   int i,j,t,count;
   float max_v = 0.0, min_v = 0.0;
   FILE *outfile,*out0, *out1, *out2, *out3, *out4, *fopen();

   outfile = fopen("Output", "w");

   out0 = fopen("Snapshot0", "w");
   out1 = fopen("Snapshot1", "w");
   out2 = fopen("Snapshot2", "w");
   out3 = fopen("Snapshot3", "w");
   out4 = fopen("Snapshot4", "w");

   for (i = 0; i < M-1; i++) 
      for (t = 0; t < STEPS; t++)  {
        if (show[i][t] > max_v)
          max_v = show[i][t];
        if (show[i][t] < min_v)
          min_v = show[i][t];
      }

printf("Maximum is %f,   Minimum is %f\n", max_v, min_v);

   for (i = 0; i < M; i++)
     for (t = 0; t < STEPS; t++)
       fprintf(outfile, "%f ", (show[i][t]-min_v)/(max_v-min_v));

   for (i = 0; i < M-1; i++) {
     fprintf(out0, "%f ", (show[i][0]-min_v)/(max_v-min_v));
     fprintf(out1, "%f ", (show[i][266]-min_v)/(max_v-min_v)); /* or 248 */
     fprintf(out2, "%f ", (show[i][286]-min_v)/(max_v-min_v)); /* or 263 */
     fprintf(out3, "%f ", (show[i][299]-min_v)/(max_v-min_v)); /* or 277 */
     /*     fprintf(out4, "%f ", (show[i][291]-min_v)/(max_v-min_v));*/
   }

   /*
   out0 = fopen("Sun", "w");
   out1 = fopen("Tree", "w");
   out2 = fopen("Mountain", "w");
   out3 = fopen("image_rest", "w");
   out4 = fopen("inhibitor", "w");

   count = 0;
   for (i = 0; i < N1; i++) * document Sun *
     for (j = 0; j < N2; j++)
       if (image1[i][j]) {
         count++;
         if (count%2 == 1)
           for (t = 0; t < STEPS; t++)
             fprintf(out0, "%f ", (show[i*N2+j][t]-min_v)/(max_v-min_v));
       }
   printf("Count for Sun == %d\n", (1+count)/2);

   count = 0;
   for (i = 0; i < N1; i++) * document Tree *
     for (j = 0; j < N2; j++)
       if (image2[i][j]) {
         count++;
         if (count%2 == 1) * reduce fig size *
           for (t = 0; t < STEPS; t++)
             fprintf(out1, "%f ", (show[i*N2+j][t]-min_v)/(max_v-min_v));
       }
   printf("Count for Tree == %d\n", (1+count)/2);

   count = 0;
   for (i = 0; i < N1; i++) * document Mountain *
     for (j = 0; j < N2; j++)
       if (image3[i][j]) {
         count++;
         if (count%3 == 1)
           for (t = 0; t < STEPS; t++)
             fprintf(out2, "%f ", (show[i*N2+j][t]-min_v)/(max_v-min_v));
       }
   printf("Count for Mountain == %d\n", (2+count)/3);

   count = 0;
   for (i = 0; i < N1; i++) * document "rest" *
     for (j = 0; j < N2; j++)
       if (input[i][j] > 0.0 && image1[i][j] == 0 && image2[i][j] == 0 &&
           image3[i][j] == 0) {
         count++;
         if (count%2 == 1)
           for (t = 0; t < STEPS; t++)
             fprintf(out3, "%f ", (show[i*N2+j][t]-min_v)/(max_v-min_v));
       }
   printf("Count for Background == %d\n", (1+count)/2);

   for (t = 0; t < STEPS; t++)
     fprintf(out4, "%f ", (show[M-1][t]-min_v)/(max_v-min_v));
   */

   fclose(outfile);
}

main()
{
    int seed;

    printf("Input a seed now\n");
    scanf("%d",&seed); 
    srand(seed); /* 213 gives a reasonable pattern for "Sun-Mountain-Tree" */

    inputs();
    connections();
    initiate();
    evaluate();

    store();

}
   
float rnd() /* returning a value in (0, 1) */
{
   return(rand()/(MAXVAL - 1.0));
}

float getxx(y,br) /* br: branch */
  float y;
  int br;
{
  float temp, cat;

  /* solving cubic equation */
  cat = -0.5*(y-2.0);
  if (y > 0.0 && y < 4.0) {
    temp = acos(cat);
    if (br)
      return(2*cos(temp/3.0));
    else
      return(2*cos(temp/3.0+0.66*3.1416));
  } else {
    temp = 0.5*sqrt(y*y-4.0*y);
    if (cat - temp > 0.0)
      return(pow(cat+temp,0.333) + pow(cat-temp,0.333));
    else if (cat + temp > 0.0)
      return(pow(cat+temp,0.333) - pow(-cat+temp,0.333));
    else
      return(-pow(-cat-temp,0.333) - pow(-cat+temp,0.333));
  }
}  

float getx(y,br) /* piecewise linear version */
  float y;
  int br;
{
  if (y > 0.0 && y < 4.0) {
    if (br)
      return(-0.25*y+2.0);
    else
      return(-0.25*y-1.0);
  } else {
    if (y > 0.0)
      return(-0.25*y-1.0);
    else
      return(-0.25*y+2.0);
  }
}  


