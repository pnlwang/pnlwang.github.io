/***************************************************************************/
/* Part of speech segregation program for the 1999 IEEE Trans. Neural Net. */
/* paper by D.L. Wang and G.J. Brown.                                      */
/* This part performs resynthesis based on segregated streams (binary      */
/* masks). It takes a segregated stream (a binary mask) and the original   */
/* mixture and produces the corresponding segregated signal.               */
/***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* booleans */

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* auditory filterbank constants */

#define MAX_CHANNEL         150                                 /* maxmimum number of filters */
#define BW_CORRECTION       1.019                               /* ERB bandwidth correction 4th order */
#define SAMPLING_FREQUENCY  16000                               /* Hz */
#define MAX_WINDOW          320                               /* use a window of 20 ms */
#define OFFSET              160                               /* offset 10 ms */
#define MAX_SIGNAL          100000
#define MAX_FRAME           700
#define HWIN                MAX_WINDOW/2.0

/* frequency scale definitions from Moore and Glasberg 1990 */

#define erb(f) (24.7*(4.37e-3*(f)+1.0))
#define hzToERBrate(f) (21.4*log10(4.37e-3*(f)+1.0))
#define ERBrateToHz(f) ((pow(10.0,((f)/21.4))-1.0)/4.37e-3)
#define sqr(x) ((x)*(x))

typedef struct {
  double cf, bw, criticalRate, z, gain, expCoeff;
  double p0, p1, p2, p3, p4;
  double q0, q1, q2, q3, q4;
  double u0, u1;
  double v0, v1;
} channel;

/* function prototypes */

void help(void);
/* UNIX help */

void blip(void);
/* write dots on stderr */

void initChannels(int lowerCF, int upperCF, int numChannels);
/* initialise filterbank channels */

double updateCochlea(channel *c, float sigval, int tim);
/* process one sample of the input though the cochlea */

int msToSamples(float ms);
/* converts time in ms to samples at srate */

/* global variables */

channel cochlea[MAX_CHANNEL];
float t, dt, twoPi, twoPiT;
float mask[MAX_FRAME][MAX_CHANNEL];

/* function declarations */

void help(void)
{
  fprintf(stderr,"-l int  lowest filter centre frequency (Hz) (500)\n");
  fprintf(stderr,"-u int  highest filter centre frequency (Hz) (2000)\n");
  fprintf(stderr,"-n int  number of channels (32)\n");
  fprintf(stderr,"-a string name of left input file\n");
  fprintf(stderr,"-b string name of right input file\n");
  fprintf(stderr,"-d float buffer decay time in ms (20.0)\n");
  fprintf(stderr,"-v bool verbose output (FALSE)\n");
}

int readSignal(float s[])
{
  int sample=0;
  int i;

  for (i=0; i<MAX_SIGNAL; i++)
    s[i]=0.0;
  while (scanf("%f",&s[sample])!=EOF)
    sample++;


  return sample;
}

void blip(void)
{
  static int count=0;

  fprintf(stderr,".");
  count+=1;
  if (count>32) count=0;
}

float DBtoAmplitude(float dB)
{
  return pow(10.0,(dB/20.0));
}

void initChannels(int lowerCF, int upperCF, int numChannels)
{
  float lowerERB, upperERB, spaceERB;
  channel c;
  int chan;

  dt = 1.0/(float)SAMPLING_FREQUENCY;
  twoPi = 2.0*M_PI;
  twoPiT = 2.0*M_PI*dt;

  lowerERB = hzToERBrate(lowerCF);
  upperERB = hzToERBrate(upperCF);

  if (numChannels > 1)
    spaceERB = (upperERB-lowerERB)/(numChannels-1);
  else
    spaceERB = 0.0;

  for (chan=0; chan<numChannels; chan++) {
    c.criticalRate = lowerERB+chan*spaceERB;
    c.cf = ERBrateToHz(c.criticalRate);
    c.bw = erb(c.cf)*BW_CORRECTION;
    c.z = exp(-twoPiT*c.bw);
    c.expCoeff = c.cf*twoPiT;
    c.gain = sqr(sqr(2*M_PI*c.bw*dt))/3.0;
    c.p0 = 0.0; c.p1 = 0.0; c.p2 = 0.0; c.p3 = 0.0; c.p4 = 0.0;
    c.q0 = 0.0; c.q1 = 0.0; c.q2 = 0.0; c.q3 = 0.0; c.q4 = 0.0;
    c.u0 = 0.0; c.u1 = 0.0;
    c.v0 = 0.0; c.v1 = 0.0;
    cochlea[chan]=c;
  }
}

double updateCochlea(channel *c, float sigval, int tim)
{
  double zz, bm;

  zz = c->z;
  c->p0 = sigval*cos(c->expCoeff*tim)+zz*(4*c->p1-zz*(6*c->p2-zz*(4*c->p3-zz*c->p4)));
  c->q0 =-sigval*sin(c->expCoeff*tim)+zz*(4*c->q1-zz*(6*c->q2-zz*(4*c->q3-zz*c->q4)));
  c->u0 = zz*(c->p1+zz*(4*c->p2+zz*c->p3));
  c->v0 = zz*(c->q1+zz*(4*c->q2+zz*c->q3));
  bm = (c->u0*cos(c->expCoeff*tim)-c->v0*sin(c->expCoeff*tim))*c->gain;

  /* filter coefficients */

  c->p4 = c->p3; c->p3 = c->p2; c->p2 = c->p1; c->p1 = c->p0;
  c->q4 = c->q3; c->q3 = c->q2; c->q2 = c->q1; c->q1 = c->q0;
  c->u1 = c->u0; c->v1 = c->v0;

  return(bm);
}

int msToSamples(float ms)
{
  return (int)((float)SAMPLING_FREQUENCY*ms/1000.0);
}

void readMask(char fname[], int maxFrame, int maxChan)
{
  int frame, chan;
  FILE *ifp;

  fprintf(stderr,"reading mask from file %s...",fname);
  ifp=fopen(fname,"r");
  if (ifp==NULL) {
    fprintf(stderr,"Cannot open file %s\n",fname);
    exit(0);
  }
  for (frame=0; frame<maxFrame; frame++)
  {
    for (chan=0; chan<maxChan; chan++)
      fscanf(ifp,"%f",&mask[frame][chan]);
        fscanf(ifp,"\n");
    }
  fprintf(stderr,"ok.\n");
}

/*------------------------------------------------------*/
/* Main program */
/*------------------------------------------------------*/

int main (int argc, char **argv)
{
        int opt;
        extern char *optarg;
        int ok = TRUE;
        int numChannels=128;
        int lowerCF=80;
        int upperCF=5000;
        float signal[MAX_SIGNAL];
        float resynth[MAX_SIGNAL];
        float nerve[MAX_SIGNAL];
        float nerve2[MAX_SIGNAL];
        float w[MAX_SIGNAL];
        int numSamples;
        int numFrames;
        int chan;
        int frame, i;
        char fname[50];

        while ((opt=getopt(argc,argv,"Hl:u:n:f:")) != EOF)
                switch(opt) {
                case 'H': help(); ok=FALSE; break;
                case 'l': lowerCF = atoi(optarg); break;
                case 'u': upperCF = atoi(optarg); break;
                case 'n': numChannels = atoi(optarg); break;
                case 'f': strcpy(fname,optarg); break;
                case '?': ok = FALSE;
                }

        if (ok==FALSE) exit(0);

        initChannels(lowerCF,upperCF,numChannels);


        /* read the signal and mask */

        numSamples=readSignal(signal);
        fprintf(stderr,"read %d samples.\n",numSamples);
        numFrames=(int)(numSamples)/OFFSET;
        numSamples=MAX_WINDOW+OFFSET*(numFrames-1);
        fprintf(stderr,"max frames=%d\n",numFrames);
        readMask(fname,numFrames,numChannels-1);

        for (i=0; i<numSamples; i++) resynth[i]=0.0;

        /* do the filtering */

        for (chan=0; chan<numChannels-1; chan++) {
          blip();
          for (i=0; i<MAX_SIGNAL; i++) w[i]=0.0;
          for (frame=-1; frame<numFrames-1; frame++) {
            for (i=0; i<MAX_WINDOW/2; i++)
            {if ((frame*OFFSET+i)>=0)
              w[frame*OFFSET+i]=w[frame*OFFSET+i]+mask[frame][chan]*0.5*(1.0+cos(i*M_PI/(HWIN)+M_PI));}
            for (i=MAX_WINDOW/2; i<MAX_WINDOW; i++)
            {if ((frame*OFFSET+i)>=0)
              w[frame*OFFSET+i]=w[frame*OFFSET+i]+mask[frame][chan]*0.5*(1.0+cos((i-HWIN)*M_PI/(HWIN)));}
          }
          for (i=0; i<numSamples; i++)
            nerve[numSamples-i-1]=(float)updateCochlea(&cochlea[chan],signal[i],i);
          for (i=0; i<numSamples; i++)
            nerve2[numSamples-i-1]=(float)updateCochlea(&cochlea[chan],nerve[i],i);
          for (i=0; i<numSamples; i++)
            resynth[i]+=w[i]*nerve2[i];
        }

        /* write */
	for (i=0; i<OFFSET; i++)  resynth[i]=0.0;
        for (i=0; i<numSamples; i++) printf("%f\n",resynth[i]);

}
