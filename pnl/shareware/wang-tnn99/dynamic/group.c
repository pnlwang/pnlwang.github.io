       /*****************************************************************/
       /* This program is for grouping in speech segregation using         */
       /* Terman-Wang oscillators. These equations are solved           */
       /* in the sigular limit (Linsay-Wang), and the potential is      */
       /* included here. Real signals.  by DeLiang Wang, 8/97           */
       /*****************************************************************/

#include "group.h"

#define MAXVAL 32768.0
               /* 2147483648.0 or 32768.0,machine (HP or SUN) dependable*/
#define MAX 32767 /* 2147483647 or 32767, machine dependable */

float Gamma = 6.5;         /* CAREFUL: this value offset Wz = 1.5 */
float alpha = 5.0;
float theta = -0.5;
float lambda = 0.1;        /* rate of developing potential */
float time = 0.0;          /* global time */
float rnd(), getx(), getxx(), lateral();
int Hp[N1][N2];        /* argument to the Heaviside for potential */
int Num = 0, longseg;           /* number of major phase updates */

float y[N1][N2],cubic[N1][N2],p[N1][N2],show[M][STEPS];
float shows[Ms][STEPS], showp[Mp][STEPS];
float W[N1][N2][N1];   /* strength of the interactions */
float input[N1][N2];  /* overall input an oscillator receives */
float We = 0.5, Wi = -0.5; /* We: total excit weight; C: count */
int branch[N1][N2], ext[N1][N2], ext2[N1][N2], segment[SEG][2];
int speech[N1][N2], phone[N1][N2]; /* for storing and display */
int chance[SEG], jumpup[SEG];
     /* branch: keeps track of whether an oscillator is on LB or RB */
     /* ext2: segment information */

inputs()
{
   FILE *fopen(), *finput, *fout;
   int i, j, k, l, maxy, temp, temp1, temp2, label[N1];
   int width, height, c;
   char magic[25];
   float stim = 0.2;

   for (i = 0; i < N1; i++) 
     for (j = 0; j < N2; j++)
       input[i][j] = -10.0;

   finput = fopen("v8n6.pitchmask.pgm", "r");
   fscanf(finput,"%s\n%d %d\n%d\n",magic,&width,&height,&c);
   for(j=0;j<N2;j++) 
     for(i=0;i<N1;i++)
        fscanf(finput, "%d",&ext[i][j]); 

   finput = fopen("v8n6-seg.pgm", "r");
   fscanf(finput,"%s\n%d %d\n%d\n",magic,&width,&height,&c);
   for(i=0;i<N1;i++)
     for(j=0;j<N2;j++)
        fscanf(finput, "%d",&ext2[N1-1-i][j]); 

   /* for storing and display purposes */
   finput = fopen("Snap1", "r");
   for(i=0;i<N1;i++)
     for(j=0;j<N2;j++)
        fscanf(finput, "%d",&speech[i][j]); 
   finput = fopen("Snap2", "r");
   for(i=0;i<N1;i++)
     for(j=0;j<N2;j++)
        fscanf(finput, "%d",&phone[i][j]); 

   /* masking by segments */
   for(i=0;i<N1;i++)
     for(j=0;j<N2;j++) 
       if (ext2[i][j] == 0)
	 ext[i][j] = 0;

   /*
   fout = fopen("phone1", "w");
   for(i=0;i<N1;i++)
     for(j=0;j<N2;j++)
        fprintf(fout, "%d ",ext[i][j]); 
   fclose(fout);

   fout = fopen("phone2", "w");
   for(i=0;i<N1;i++)
     for(j=0;j<N2;j++)
        fprintf(fout, "%d ",ext2[i][j]); 
   fclose(fout);
   */

   /* identify longest segment and clip accordingly */
   for (i = 0; i < SEG; i++) {
     segment[i][0] = N2;
     segment[i][1] = 0;
   }

   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++) {
       temp = ext2[i][j];
       if (segment[temp][0] > j)
	 segment[temp][0] = j;
       if (segment[temp][1] < j)
	 segment[temp][1] = j;
     }

   maxy = 0;
   for (i = 1; i < SEG; i++) /* excluding the background */
     if (segment[i][1] - segment[i][0] + 1 > maxy) {
       maxy = segment[i][1] - segment[i][0] + 1;
       longseg = i;
     }
   printf("longest segment is %d, and has columns %d\n", longseg, maxy);

   for (j = 0; j < segment[longseg][0]; j++)
     for (i = 0; i < N1; i++)
       ext[i][j] = 0;
   for (j = segment[longseg][1] + 1; j < N2; j++)
     for (i = 0; i < N1; i++)
       ext[i][j] = 0;

   /* recompute segment lengths */
   k = segment[longseg][0]; l = segment[longseg][1];
   for (i = 0; i < SEG; i++) {
     segment[i][0] = N2;
     segment[i][1] = 0;
   }

   for (i = 0; i < N1; i++)
     for (j = k; j <= l; j++) {
       temp = ext2[i][j];
       if (segment[temp][0] > j)
	 segment[temp][0] = j;
       if (segment[temp][1] < j)
	 segment[temp][1] = j;
     }

   for (i = 1; i < SEG; i++)
     if (segment[i][1] - segment[i][0] < 2) { /* threshold out <3-frame seg */
       for(k=0;k<N1;k++)                      /* caused by clipping         */ 
         for(l=0;l<N2;l++)
           if (ext2[k][l] == i)
	     ext[k][l] = 0;
     }

   for(i=0;i<N1;i++)
     for(j=0;j<N2;j++)
       if (ext[i][j] > 0)
	 input[i][j] = stim;


   /* smoothen within each column of each segment */
   for (j = segment[longseg][0]; j <= segment[longseg][1]; j++) {
     for (i = 0; i < N1; i++)
       label[i] = 0;

     for (i = 0; i < N1; i++)
       if (label[i] == 0 && ext[i][j] > 0) {
	 temp = ext2[i][j]; temp1 = temp2 = 0;
	 for (k = 0; k < N1; k++)
	   if (ext2[k][j] == temp) {
	     if (ext[k][j] == 1)
	       temp1++;
	     else if (ext[k][j] == 2)
	       temp2++;
	   }
	 for (k = 0; k < N1; k++)
	   if (ext2[k][j] == temp && ext[k][j] > 0) {
	     label[k] = 1;
	     if (temp1 > temp2)
	       ext[k][j] = 1;
             else
	       ext[k][j] = 2;
	   }
       }
   }


   fout = fopen("phone3", "w");
   for(i=0;i<N1;i++)
     for(j=0;j<N2;j++)
        fprintf(fout, "%d ",ext[i][j]); 
   fclose(fout);

/* high "potential" */
   input[8][N2/2] = 0.21; /* 20 for v3n6; 8 for v8n6 */

   printf("The segment number with high potential is %d\n", ext2[8][N2/2]);
}

connections()
{
   int i, j, k; /* e: excitatory; i: inhibitory; l: lateral */

   for (i = 0; i < N1; i++) {
     for (j = segment[longseg][0]; j <= segment[longseg][1]; j++) {
       for (k = 0; k < N1; k++) /* initialize connections */
	 W[i][j][k] = 0.0;
       if (ext[i][j] == 1) {
         for (k = 0; k < N1; k++) 
	   if (k != i && ext2[k][j] != ext2[i][j] ) {
	     if (ext[k][j] == 1)
	       W[i][j][k] = We;
             else if (ext[k][j] == 2)
	       W[i][j][k] = Wi;
	   }
       }
       else if (ext[i][j] == 2)
         for (k = 0; k < N1; k++) 
	   if (k != i && ext2[k][j] != ext2[i][j] ) {
  	     if (ext[k][j] == 2)
	       W[i][j][k] = We;
             else if (ext[k][j] == 1)
	       W[i][j][k] = Wi;
	   }
     }
   }

   printf("weights are \n");
   printf("weights are \n");
   printf("W[0][0] are %f %f\n", W[0][0][1], W[0][0][2]);
   printf("W[1][4] are %f %f\n", W[1][4][0], W[1][4][2]);
   printf("W[2][6] are %f %f\n", W[2][6][0], W[2][6][1]);

   printf("\n");
}

initiate()
{
   int i,j,seed,counts,countp;  /* labels for setting up display intervals */
   float angle, temp;

   printf("Input a seed now\n"); /* 2342 for potential; 111 no potent.*/
   scanf("%d",&seed); 
   srand(seed); 

/* randomize the initial phases */
   temp = 2.0*Gamma*rnd();
   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++) {  /* Initial phases */
       if (ext[i][j] > 0)
	 y[i][j] = temp + input[i][j];
       else
	 y[i][j] = 2.0*Gamma*rnd();

       branch[i][j] = 0;

       /*
       p[i][j] = 0.0;

       Hp[i][j] = 0;
       */

       cubic[i][j] = input[i][j];
     }

   for (i = 0; i < SEG; i++) {
     chance[i] = 1;
     jumpup[i] = 0;
   }

 /* initial display values for 1D */
   for (i = 0; i < N1; i = i+2)
       show[i/2][0] = getxx(y[i][COL]-cubic[i][COL], 0);

 /* initial display values for two storage arrays */
   counts = countp = 0; /* for collecting display data for two streams */
   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++) {
	if (speech[i][j] && (counts%20 == 0))
          shows[counts/20][0] = getxx(y[i][j]-cubic[i][j],0);
	else if (phone[i][j] && (countp%20 == 0))
          showp[countp/20][0] = getxx(y[i][j]-cubic[i][j],0);

	if (speech[i][j])
	  counts++;
	else if (phone[i][j])
	  countp++;
     }

 /* initial display values for 2D */
   /*
   for (i = 0; i < N1; i++)
     for (j = 0; j < N2; j++)
       show[N2*i+j][0] = getxx(y[i][j]-cubic[i][j], 0);
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
     for (j = segment[longseg][0]; j <= segment[longseg][1]; j++) {
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
   int i, j, k, l, count, jumped = 1, Tij = 2, flag = 1, cat[N1][N2];
                                               /* Tij: permanent weight */
   int temp, curseg;
   float We1 = 4.0, We2 = 0.1;

   while (jumped) {
      jumped = 0;
      for (i = 0; i < N1; i++)
	for (j = 0; j < N2; j++)
	  cat[i][j] = 0;

      for (i = 0; i < N1; i++)
	for (j = segment[longseg][0]; j <= segment[longseg][1]; j++) 
	if (ext[i][j] > 0 && ext2[i][j] > 0) { 
                                                      /* compute net input */
	  if (p[i][j] > theta)
	                                      /* thrshld for exp(-alpha*t) */
	    cubic[i][j] = input[i][j];
          else
	    cubic[i][j] = baseline;

	  /* update cubics within each segment */
	  count = 0;
	  curseg = ext2[i][j];
	  for (k = segment[curseg][0]; k <= segment[curseg][1]; k++) {
	    temp = 0;
	    for (l = 0; l < N1; l++)
	      if (ext2[l][k] == ext2[i][j] && branch[l][k])
		temp = 1;
	    if (temp == 1)
	      count++;
	  }

	  temp = ext2[i][j];
          if (ext[i][j]>0 && count>(1+segment[temp][1]-segment[temp][0])/2)  
	                           /* strict majarity rule: >; less >= (1+.) */
	    cubic[i][j] += We1;
	  else if (ext[i][j] > 0 
		   && segment[temp][1] == segment[temp][0]+count*2-1) { 
	                           /* break symmetry */
	    if (chance[curseg])
	      cubic[i][j] += We2;
	    else
	      cubic[i][j] += We1;
	  }
	  else if (ext[i][j] > 0 && count > 0)
	    cubic[i][j] += We2;

	  /* update cubics within each column */
	  if (branch[i][j] == 0) 
	    cubic[i][j] += lateral(i,j);

	  if (!branch[i][j] && y[i][j]<cubic[i][j]+error) {  /* jump up */
  	        cat[i][j] = 1;
		jumped = 1;
		cubic[i][j] -= lateral(i,j);
		jumpup[curseg] = 1;
	  } else if (branch[i][j] && y[i][j]>cubic[i][j]+4.0-error) {
	                                                     /* jump down */
	        cat[i][j] = -1;
		jumped = 1;
		cubic[i][j] += lateral(i,j);
	  }

      }

      for (i = 0; i < N1; i++)
	for (j = 0; j < N2; j++) { /* update jumps */
	  if (cat[i][j] > 0)
	    branch[i][j] = 1;
	  else if (cat[i][j] < 0)
	    branch[i][j] = 0;
	}
   }

   for (i = 1; i < SEG; i++)
     if (jumpup[i]) {
       temp = 0;

       for (k = segment[i][0]; k <= segment[i][1]; k++) 
	 for (l = 0; l < N1; l++)
	   if (ext2[l][k] == i && branch[l][k])
	     temp = 1;

       if (temp == 0) {
	 chance[i] = 1 - chance[i];
	 jumpup[i] = 0;
       }
     }

/* compute Hp[i][j] */
   /*
   for (i = 0; i < N1; i++)  
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
   */
}

evaluate()
{
   float intv, mu, current, length, xval, yval, angle, temp, cat;
   int ilead, jlead, num; /* leading oscillator, and no of jumpers */
   int counts, countp;    /* counts for speech and phone respectively */
   int i,j,step,bef,aft;  /* labels for setting up display intervals */
   FILE *outf, *fopen();

   outf = fopen("Out", "w");
   bef = 0;
   length = delta*(STEPS-1);      /* reduce one step to facilitate display */
   while (time < length) {
     findtmin(&intv, &mu, &ilead, &jlead);

     aft = (time+intv)/delta;
     
    /* store for display in 1d case*/

     for (step = bef; step < aft && step < STEPS-1; step++) {
       for (i = 0; i < N1; i = i+2) {
           temp = -((1.0+step)*delta - time);
	   if (branch[i][COL])
	     yval = (y[i][COL] - 2*Gamma)*exp(temp) + 2*Gamma - cubic[i][COL];
	   else
	     yval = y[i][COL]*exp(temp) - cubic[i][COL];

           show[i/2][step+1] = getxx(yval,branch[i][COL]);
       }
       counts = countp = 0; /* for collecting display data for two streams */
       for (i = 0; i < N1; i++)
         for (j = 0; j < N2; j++) {
	   if (speech[i][j] && (counts%20 == 0)) {
             temp = -((1.0+step)*delta - time);
	     if (branch[i][j])
	       yval = (y[i][j] - 2*Gamma)*exp(temp) + 2*Gamma - cubic[i][j];
	     else
	       yval = y[i][j]*exp(temp) - cubic[i][j];

             shows[counts/20][step+1] = getxx(yval,branch[i][j]);
	   } else if (phone[i][j] && (countp%20 == 0)) {
             temp = -((1.0+step)*delta - time);
	     if (branch[i][j])
	       yval = (y[i][j] - 2*Gamma)*exp(temp) + 2*Gamma - cubic[i][j];
	     else
	       yval = y[i][j]*exp(temp) - cubic[i][j];

             showp[countp/20][step+1] = getxx(yval,branch[i][j]);
	   }

	   if (speech[i][j])
	     counts++;
	   else if (phone[i][j])
	     countp++;
         }
     }


    /* store for display in 2d case*/
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
     }
*/

     for (i = 0; i < N1; i++) /* advance oscillators */
       for (j = 0; j < N2; j++) {
	 if (branch[i][j])
	   y[i][j] = (y[i][j]-2.0*Gamma)/mu + 2.0*Gamma;
	 else
	   y[i][j] = y[i][j]/mu;
       }

/* update oscillator potential */
     /*
     for (i = 0; i < N1; i++) 
       for (j = 0; j < N2; j++)
	 if (Hp[i][j] && input[i][j] > 0.0)
	   p[i][j] = 1.0-(1.0-p[i][j])/mu;

     */

     if (branch[ilead][jlead]) {
       branch[ilead][jlead] = 0;
       /*       priming(jlead);*/
       cubic[ilead][jlead] += lateral(ilead,jlead);
     }
     else {
       branch[ilead][jlead] = 1;
       /*       priming(jlead);*/
       cubic[ilead][jlead] -= lateral(ilead,jlead);
     }

     jump();            /* fast time activity */

     Num++;
     printf("Num = %d, time = %f, intv = %f\n", Num, time, intv);

     for(i=0;i<N1;i++)
       for(j=0;j<N2;j++)
         fprintf(outf, "%d ",branch[i][j]); 
   
     time += intv;      
     bef = aft;
   }
   fclose(outf);
}

store()
{
   int i,j,t,count;
   float max_v = 0.0, min_v = 0.0;
   FILE *outfile,*out,*fopen();

   outfile = fopen("Output", "w");
   out = fopen("NUMBER", "w");
   fprintf(out, "%d ", Num);
   printf("The number of the segments/cycles is %d\n", Num);

   /*
   for (i = 0; i < M; i++) 
      for (t = 0; t < STEPS; t++)  {
        if (show[i][t] > max_v)
          max_v = show[i][t];
        if (show[i][t] < min_v)
          min_v = show[i][t];
      }
      */
   for (i = 0; i < Ms; i++) 
      for (t = 0; t < STEPS; t++)  {
        if (shows[i][t] > max_v)
          max_v = shows[i][t];
        if (shows[i][t] < min_v)
          min_v = shows[i][t];
      }
   for (i = 0; i < Mp; i++) 
      for (t = 0; t < STEPS; t++)  {
        if (showp[i][t] > max_v)
          max_v = showp[i][t];
        if (showp[i][t] < min_v)
          min_v = showp[i][t];
      }

printf("Maximum is %f,   Minimum is %f\n", max_v, min_v);

   for (i = 0; i < M; i++)
     for (t = 0; t < STEPS; t++)
       fprintf(outfile, "%f ", (show[i][t]-min_v)/(max_v-min_v));
   fclose(outfile); 

   /*
   outfile = fopen("Speech", "w");
   for (i = 0; i < Ms; i++)
     for (t = 0; t < STEPS; t++)
       fprintf(outfile, "%f ", (shows[i][t]-min_v)/(max_v-min_v));
   fclose(outfile); 

   outfile = fopen("Phone", "w");
   for (i = 0; i < Mp; i++)
     for (t = 0; t < STEPS; t++)
       fprintf(outfile, "%f ", (showp[i][t]-min_v)/(max_v-min_v));
   fclose(outfile); 
   */
   fclose(out);
}

main()
{
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

float lateral(i,j) /* determine lateral connections */
  int i, j;
{
  int k;
  float temp1 = 0.0, temp2 = 0.0, base = 0.05;

  for (k = 0; k < N1; k++)
    if (k != i && ext2[k][j] != ext2[i][j] && branch[k][j] && ext[i][j] > 0) {
      if (W[i][j][k] > base)
	temp1 = We;
      else if (W[i][j][k] < -base)
	temp2 = Wi;
    }
  if (temp2 < -base)
    return(Wi);
  else if (temp1 > base)
    return(We);
  else return(0.0);
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


