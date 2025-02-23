/***************************************************************/
/* This C program implements the feature-preserving smoothing  */
/* (FPS) algorithm by K. Chen, D.L. Wang and X. Liu.           */
/* The algorithm was published in IEEE Trans. on Neural        */
/* Networks, vol. 11, No. 5, pp. 1106-1123, in an article      */
/* entitled "Weight adaptation and oscillatory correlation     */
/* for image segmentation". The algorithm performs selective   */
/* smoothing in a 3X3 window by incorporating variance.        */
/* Originally written by Ke Chen, Nov. 1997.                   */
/* Modified by DeLiang Wang, Feb. 2001.                        */
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N1  230        /* used to be 175: Height of an image */
#define N2  240        /* used to be 175: Width */
#define Maxpixel 255   /* brightest possible pixel intensity */

   /* following four are adjustable parameters */
#define Kappa 5.0   /* In Eqs. 11 and 15 */
#define SCALE 9.0 /* "s" in Eqs. 10 & 15 */
#define theta_sig 0.5 /* in the definition of function Phi */
#define R 9 /* Window size: (2*R+1)*(2*R+1) */

#define T 40    /* iteration number */
#define MOMENT 40 /* save partial result every MOMENT iterations */

float imbuf[N1][N2], weight[N1][N2], var_mt[N1][N2];
  /* image buffer, local discontinuity weights, and variance-based weights */
void inputs(); /* input original image */
float adapt(); /* calculate second line of Eq. 15, major step in smoothing */
float Dx(); /* D_H */
float Dy(); /* D_V */
float Dplus(); /* D_C */
float Dminus(); /* D_D */
float absd(); /* absolute value for "float" */  
void cal_weight(); /* calculate local discontinuity and related weights */
void cal_varcoef();
float variance(); /* calculate mean and variance */
float dmax(); /* maximum variance: sigma_2{max} */
float dmin(); /* minimum variance: sigma_2{min} */
int mod(); /* find modulo */

void inputs(stim)  /* input a 2D image of the size N1xN2 in .pgm format */
int stim[N1][N2];
{
   int i, j;
   int width, height, c;
   char magic[25], bigbuf[500];

   FILE *IN;

   IN = fopen("images/noise.pgm", "r");
   fscanf(IN,"%s\n", magic);
   while (1) {
     fscanf(IN, "%s", bigbuf);
     if (bigbuf[0] != '#') {
        sscanf(bigbuf, "%d", &width);
        break;
        }
     else
        while (getc(IN) != '\n');
     }
   fscanf(IN,"%d\n%d\n",&height,&c);
   fprintf(stderr,"[Type %s %dx%d colours %d]\n",magic,width,height,c);

   for(i=0;i<N1;i++)
     for(j=0;j<N2;j++)
       stim[i][j] = getc(IN); 
       /* fscanf(IN, "%d ", &stim[i][j]); */
   printf("Stimuli are %d %d %d %d\n",stim[0][0],stim[14][85],
              stim[10][10], stim[17][17]);

   fclose(IN);
}


void cal_weight() /* calculate local discontinuity measure */
{
  int i, j;
  float gx, gy, gplus, gminus, tmp; 

  for (i=0; i<N1; i++)
    for (j=0; j<N2; j++) {
      gx = Dx(i,j);
      gy = Dy(i,j);
      gplus = Dplus(i,j);
      gminus = Dminus(i,j);
      tmp = (absd(gx) + absd(gy) + absd(gplus) + absd(gminus))/4.0;
      weight[i][j] = exp(-tmp/SCALE);
    }
}

int mod(x,y) /* find modulo of x wrt y */
int x,y;
{
  int rtn;

  rtn = x;

  while (rtn >= y)
    rtn -= y;
 
  return rtn;
}


float adapt(i, j) /* calculate second line of Eq. 15 */
int i, j;
{
    int m, n;
    float mbuf, sum, rtn, test1;
    
    mbuf = 0.0; sum = 0.0;
    for (m = -1; m < 2; m++) /* summation over a 3x3 window */
      for (n = -1; n< 2; n++) 
        if ((i+m >=0) && (j+n>=0) && (i+m < N1) && (j+n < N2)) {
  mbuf += (imbuf[i+m][j+n] - imbuf[i][j])*weight[i+m][j+n]*var_mt[i+m][j+n];
             sum += weight[i+m][j+n] * var_mt[i+m][j+n]; 
        }
        
   if (sum > 0.0)
      rtn = mbuf / sum;
   else
      rtn = 0.0;

   return rtn;
}

float Dy(x,y)
int x, y;
{
   int y1, y2;
   float rtn;

   if (y == 0) {
      y1 = 0;
      y2 = 1;
    }
   else if (y == (N2-1)) {
      y1 = N2 - 2;  
      y2 = N2 - 1;
   }
   else {
      y1 = y - 1;
      y2 = y + 1;
   }

   rtn = imbuf[x][y2] - imbuf[x][y1];

   return rtn;
}

float Dx(x,y)
int x, y;
{
   int x1, x2;
   float rtn;

   if (x == 0) {
      x1 = 0;
      x2 = 1;
   } 
   else if (x == (N1-1)) {
      x1 = N1-2;
      x2 = N1-1;
   }
   else {
      x1 = x-1;
      x2 = x+1;
   }

   rtn = imbuf[x2][y] - imbuf[x1][y];

   return rtn;
}

float Dplus(x,y)
int x, y;
{
   int x1, y1, x2, y2;
   float rtn;

   if (x == 0 ) 
      x1 = 0;
   else
      x1 = x - 1;

   if (x == (N1-1))
      x2 = N1 - 1;
   else
      x2 = x + 1;

   if (y == 0)
      y1 = 0;
   else
      y1 = y - 1;

   if (y == (N2-1))
      y2 = N2 - 1;
   else
      y2 = y + 1;

  
  rtn = imbuf[x2][y2] - imbuf[x1][y1];

  return rtn;
}

float Dminus(x,y)
int x, y;
{
  int x1, y1, x2, y2;
  float rtn;

  if (x == 0)
     x1 = 0;
  else
     x1 = x - 1;
     
  if (y == (N2-1))
     y1 = N2 - 1;
  else
     y1 = y + 1;

  if (x == (N1-1))
     x2 = N1 - 1;
  else
     x2 = x + 1;


  if (y == 0)
     y2 = 0;
  else
     y2 = y - 1;

  rtn = imbuf[x2][y2] - imbuf[x1][y1];
  return rtn;
}


float absd(x) /* absolute value of "float" type */
float x;
{
  float rtn;

  if (x >= 0)
     rtn = x;
  else
     rtn = -x;

  return rtn;
}

void cal_varcoef()
{
  int i, j;
  float mdif, minvar, maxvar, tmp, varbuf[N1][N2];
  /* varbuf: variance buffer for the whole image */
  FILE *contxt;

  contxt = fopen("Context.pgm", "w");  /* contextual discontinuity map */
  fprintf(contxt,"%s\n%d %d\n%d\n","P2",N2,N1,255);

  for (i=0; i<N1; i++)
    for (j=0; j<N2; j++)
      varbuf[i][j] = variance(i,j);

  minvar = dmin(varbuf);
  maxvar = dmax(varbuf);
  mdif = maxvar - minvar;

  for (i=0; i<N1; i++)
    for (j=0; j<N2; j++) {
      tmp = (varbuf[i][j] - minvar)/mdif; 
/* sigma-hat: see Eq. 7. Value btw 0 and 1 */

      if (tmp > theta_sig)
         fprintf(contxt, "%d ", (int)(tmp*Maxpixel));
      else
         fprintf(contxt, "%d ", 0);

      if (tmp > theta_sig)
         var_mt[i][j] = exp(-Kappa*tmp);
      else
	 var_mt[i][j] = 1.0;
    }
  fclose(contxt);
}


float variance(i, j)  /* calculate mean and variance in an image window */
int i, j;
{
    int m, n;
    float rtn, ctr, mbuf, vbuf; 
    /* ctr: counter; mbuf: mean buffer; vbuf: variance buffer */

    ctr = 0.0; mbuf = 0.0; vbuf = 0.0;
    for (m = -R; m < R + 1; m++)
      for (n = -R; n< R + 1; n++)
        if ((i+m >=0) && (j+n>=0) && (i+m < N1) && (j+n < N2)) {
            mbuf += imbuf[i+m][j+n];
            vbuf += imbuf[i+m][j+n] * imbuf[i+m][j+n];
            ctr++;
        }

    if (ctr > 1.0) {
       mbuf = mbuf/ctr;
       rtn = vbuf/ctr - mbuf*mbuf;
    }
    else {
      printf("there is an error in moment.\n");
      exit(0);
      }

    return rtn;
}


float dmax(mt) /* maximum spatial variance from variance matrix: mt */
float mt[N1][N2];
{
   int i, j;
   float rtn = -1.0;

   for (i=0; i<N1; i++)
     for (j=0; j<N2; j++)
       if (mt[i][j] > rtn)  rtn = mt[i][j];

   return rtn;
}


float dmin(mt) /* minimum spatial variance from variance matrix: mt */
float mt[N1][N2];
{
   int i, j;
   float rtn;

   rtn = mt[0][0];
   for (i=0; i<N1; i++)
     for (j=0; j<N2; j++)
       if (mt[i][j] < rtn)  rtn = mt[i][j];

   return rtn;
}


main()
{
     int i, j, t, imagedata[N1][N2];
     float outbuf[N1][N2];
     char  dname1[100], dname2[100], numb1[20];
     FILE *ftr;

     printf("Give the filename to save the smoothed image: ");
     scanf("%s",numb1);
     strcpy(dname1,"./result/");
     strcat(dname1,numb1);

     inputs(imagedata);

     for (i=0; i<N1; i++)
       for (j=0; j<N2; j++) 
           imbuf[i][j] = (float) imagedata[i][j];

     cal_varcoef();

     for (t=0; t <= T; t++) {

       cal_weight();

       for (i=0; i<N1; i++) /* calculate Eq. 15 */
         for (j=0; j<N2; j++) 
           outbuf[i][j] = imbuf[i][j] + var_mt[i][j]*adapt(i,j);

       for (i=0; i<N1; i++)
	 for (j=0; j<N2; j++) 
	   imbuf[i][j] = outbuf[i][j];
       
       printf("iteration %d is done ...\n",t+1);

       /* the following stores smoothed images to .pgm files */     
      if ((mod(t,MOMENT) == 0) && (t >0)) {
 	sprintf(dname2,"%s_%d.pgm",dname1,t);
     
        ftr = fopen(dname2,"w");
        fprintf(ftr,"%s\n%d %d\n%d\n","P2",N2,N1,Maxpixel);
        for (i=0; i<N1; i++)
          for (j=0; j<N2; j++)
	    fprintf(ftr, "%d ", (int)(outbuf[i][j])); 

        fclose(ftr);
      }
    }
}
