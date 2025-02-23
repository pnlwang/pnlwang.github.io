/***************************************************************************/
/* Speech segregation program for the 1999 IEEE Trans. Neural Net. paper   */
/* by D.L. Wang and G.J. Brown.                                            */
/* The segments produced by this program have been discarded since an      */
/* older version was used. The two useful outputs are:                     */
/* (1) The file "PITCHMASK"; (2) The file "ENERGY". The latter is a        */
/* map (frequency-time) representation of the Meddis hair cell response.   */
/***************************************************************************/

#include <stdio.h>
#include <math.h>

#define MAX_CHANNEL 128
#define MAX_DELAY 200
#define MIN_DELAY 32
#define MAX_FRAME 150 /* 323 or 150 */
#define MIN_FRAME 0
#define MAX_SEGMENT 300

/* threshold WAS 0.95 */

#define THRESHOLD 0.95
#define ENERGY_THRESHOLD 50.0
#define GROUP_THRESHOLD 0.985
#define TRUE 1
#define FALSE 0
#define sqr(x) ((x)*(x))

struct _point {
  int frame;
  int chan;
  struct _point *next;
} _point;

typedef struct _point point;

struct _group {
  int first;
  int last;
  int size;
  short group_id;
  point *head;
  short oneCount[MAX_FRAME];
  short twoCount[MAX_FRAME];
} _group;

typedef struct _group group;

float acg[MAX_DELAY][MAX_CHANNEL];
float norm_acg[MAX_DELAY][MAX_CHANNEL];
float correlation[MAX_FRAME][MAX_CHANNEL];
float summary[MAX_DELAY];
short channels[MAX_FRAME][MAX_CHANNEL];
short energy[MAX_FRAME][MAX_CHANNEL];   /* newly added */
short mask[MAX_FRAME][MAX_CHANNEL];
group segment[MAX_SEGMENT];
short groupsize[MAX_SEGMENT];
short sourcemask[MAX_FRAME][MAX_CHANNEL];

int count_one[MAX_FRAME];
int count_two[MAX_FRAME];
int group_first;
int group_last;

void clearSourceMask()
{
  int frame, chan;

  for (frame=0; frame<MAX_FRAME; frame++)
    for (chan=0; chan<MAX_CHANNEL; chan++)
      sourcemask[frame][chan]=0;
}

void initGroup(group *g)
{
  int frame;

  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++) {
    g->oneCount[frame]=0;
    g->twoCount[frame]=0;
  }
  g->first=MAX_FRAME;
  g->last=MIN_FRAME;
  g->size=0;
  g->head=NULL;
  g->group_id=0;
}

void zeroCounts()
{
  int frame;

  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++) {
    count_one[frame]=0;
    count_two[frame]=0;
  }
}

void initChannels()
{
  int frame, chan;

  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++)
    for (chan=0; chan<MAX_CHANNEL; chan++) {
      channels[frame][chan]=2;
      mask[frame][chan]=0;
    }
}

void readCorrelogram(char *filename)
{
   int delay, chan, dummy;
   FILE *infile;

   infile = fopen(filename, "r");
   if (infile==NULL) {
       fprintf(stderr,"error opening file %s\n",filename);
       exit(1);
   }
   fscanf(infile,"%d %d",&dummy,&dummy);
   for (chan=0; chan<MAX_CHANNEL; chan++)
       for (delay=0; delay<MAX_DELAY; delay++)
           fscanf(infile,"%f",&acg[delay][chan]);
   
   fclose(infile);
}

void makeSummary(int frame)
{
  int delay, chan;

  for (delay=0; delay<MAX_DELAY; delay++)
    summary[delay]=0.0;
  for (delay=0; delay<MAX_DELAY; delay++) 
    for (chan=0; chan<MAX_CHANNEL; chan++)
      if (channels[frame][chan]==2)
	summary[delay]+=acg[delay][chan];
}

int findSummaryMax()
{
  int delay, maxDelay;
  float maxVal;

  maxVal=summary[MIN_DELAY];
  maxDelay=MIN_DELAY;
  for (delay=MIN_DELAY+1; delay<MAX_DELAY; delay++)
    if (summary[delay]>maxVal) {
      maxVal=summary[delay];
      maxDelay=delay;
    }
  return maxDelay;
}

void selectChannels(int frame, int maxd)
{
  int chan, delay;
  float maxval;

  for (chan=0; chan<MAX_CHANNEL; chan++) 
    if (channels[frame][chan]==2) {
      maxval=acg[MIN_DELAY][chan];
      for (delay=MIN_DELAY+1; delay<MAX_DELAY; delay++)
	if (acg[delay][chan]>maxval) maxval=acg[delay][chan];
      if ((acg[maxd][chan]/maxval)>THRESHOLD)
	channels[frame][chan]=1;
    }
}

void writeChannels()
{
  int chan, frame;
  FILE *ofp;

  ofp=fopen("PITCHMASK","w");
  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++) {
    for (chan=0; chan<MAX_CHANNEL-1; chan++) 
      fprintf(ofp,"%1d ",channels[frame][chan]*mask[frame][chan]);
    fprintf(ofp,"\n");
  }
  fclose(ofp);
}

void writeMask()
{
  int chan, frame;

  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++) {
    for (chan=0; chan<MAX_CHANNEL-1; chan++)
      printf("%d ",mask[frame][chan]);
    printf("\n");
  }
}

void weightSummary()
{
  int delay;

  for (delay=0; delay<MAX_DELAY; delay++)
    summary[delay]*=(1.0-0.3*delay/(MAX_DELAY-1.0));
}

void writeSourceMaskToFile(char fname[])
{
  int chan, frame;
  FILE *ofp;

  ofp=fopen(fname,"w");
  if (ofp==NULL) {
    fprintf(stderr,"cannot open file %s\n",fname);
    exit(0);
  }
  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++) {
    for (chan=0; chan<MAX_CHANNEL-1; chan++)
      fprintf(ofp,"%d ",sourcemask[frame][chan]);
    fprintf(ofp,"\n");
  }
  fclose(ofp);
}

void writeSummary()
{
  int delay;

  for (delay=0; delay<MAX_DELAY; delay++)
    printf("%1.2f ",summary[delay]);
  printf("\n");
}

void normalise()
{
  float mean, std;
  int delay, chan;

  for (chan=0; chan<MAX_CHANNEL; chan++) {
    mean=0.0;
    std=0.0;
    for (delay=0; delay<MAX_DELAY; delay++)
      mean+=acg[delay][chan];
    mean=mean/(float)MAX_DELAY;
    for (delay=0; delay<MAX_DELAY; delay++)
      std+=sqr(acg[delay][chan]-mean);
    std=sqrt(std/(float)MAX_DELAY);
    for (delay=0; delay<MAX_DELAY; delay++)
      norm_acg[delay][chan]=(acg[delay][chan]-mean)/std;
  }
}

void doCorrelation(int frame)
{
  float xy;
  int chan, delay;

  for (chan=0; chan<MAX_CHANNEL-1; chan++) {
    xy=0.0;
    for (delay=0; delay<MAX_DELAY; delay++) {
      xy+=norm_acg[delay][chan]*norm_acg[delay][chan+1];
    }
    correlation[frame][chan]=xy/(float)MAX_DELAY;
    if ((correlation[frame][chan]>GROUP_THRESHOLD) && (acg[0][chan]>ENERGY_THRESHOLD)) {
      mask[frame][chan]=1;
    }
  }
}

void addPoint(group *g, int f, int c)
{
  point *p;

  //fprintf(stderr,"adding (%d,%d)\n",f,c);
  p=(point *)malloc(sizeof(point));
  p->next=g->head;
  p->frame=f;
  p->chan=c;
  g->size++;
  g->head=p;
  if (f < g->first) 
    g->first=f;
  if (f > g->last) 
    g->last=f;
}

void traverse(point *p) 
{
  if (p!=NULL) {
    printf("%d %d\n",p->frame,p->chan);
    traverse(p->next);
  }
}

void addToMask(point *p)
{
  if (p!=NULL) {
    sourcemask[p->frame][p->chan]=1;
    addToMask(p->next);
  }
}

void traverseToFile(point *p, FILE *ofp)
{
  if (p!=NULL) {
    fprintf(ofp,"%d %d\n",p->frame,p->chan);
    traverseToFile(p->next,ofp);
  }
}

void writeGroup(int n)
{
  char fname[30];
  FILE *ofp;
  int frame;
  
  sprintf(fname,"SEGMENT.%03d",n);
  ofp=fopen(fname,"w");
  if (ofp==NULL) {
    fprintf(stderr,"ARGH!!! CANNOT OPEN FILE\n");
    exit(0);
  }
  fprintf(ofp,"first=%d last=%d size=%d\n",segment[n].first,segment[n].last,segment[n].size);
  fprintf(ofp,"Count for group 1:\n");
  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++)
    fprintf(ofp,"%d ",segment[n].oneCount[frame]);
  fprintf(ofp,"\n");
  fprintf(ofp,"Count for group 2:\n");
  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++)
    fprintf(ofp,"%d ",segment[n].twoCount[frame]);
  fprintf(ofp,"\n");
  traverseToFile(segment[n].head,ofp);
  fclose(ofp);
}

void searchGroup(group *g, int frame, int chan)
{
  if (mask[frame][chan]==1) {
    addPoint(g,frame,chan);
    mask[frame][chan]=0;
    if (channels[frame][chan]==1) g->oneCount[frame]++;
    if (channels[frame][chan]==2) g->twoCount[frame]++;
    if (chan<MAX_CHANNEL-1) 
      searchGroup(g,frame,chan+1);
    if (chan>0)
      searchGroup(g,frame,chan-1);
    if (frame>0)
      searchGroup(g,frame-1,chan);
    if (frame<MAX_FRAME-1)
      searchGroup(g,frame+1,chan);
  }
}

int add_pitch_info(int n)
{
  int frame;

  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++) {
    count_one[frame]+=segment[n].oneCount[frame];
    count_two[frame]+=segment[n].twoCount[frame];
  }
}

int matchPitch(int n, int frame) 
/* matches the periodicity information in segment "n" with the current group
at time frame "frame" */
{
  int test1, test2;

  test1=(segment[n].oneCount[frame] >= segment[n].twoCount[frame]);
  test2=(count_one[frame] >= count_two[frame]);
  return (test1 && test2) | (!test1 && !test2);
}

int min(int x1, int x2)
{
  return (x1<x2) ? x1 : x2;
}

int max(int x1, int x2)
{
  return (x1>x2) ? x1 : x2;
}

int matchSegment(int n)
{
  int st, fn, cnt, frame;
 
  cnt=0;
  st=max(segment[n].first,group_first);
  fn=min(segment[n].last,group_last);
  for (frame=st; frame<=fn; frame++)
    cnt+=matchPitch(n,frame);
  //fprintf(stderr,"TARGET=%d ST=%d FN=%d  CNT=%d\n",n,st,fn,cnt);
  if (((float)cnt/((float)fn-(float)st))>0.5)
    return TRUE;
  else
    return FALSE;
}

int overlap(int n)
/* does segment n overlap with group? */
{
  int s1, s2, f1, f2;

  s1=group_first;
  f1=group_last;
  s2=segment[n].first;
  f2=segment[n].last;

  return ((s1>=s2)&&(f1<=f2)) | ((s2>=s1)&&(f2<=f1)) | ((s2>s1)&&(f1<f2)&&(f1>s2)) | ((s1>s2)&&(f2<f1)&&(f2>s1));
}

main()
{
  int maxDelay;
  int frame, chan;
  char fname[50];
  short pg[MAX_DELAY];
  int delay;
  int numSegment;
  group newSegment;
  int seed, s, s1, MaxSize, groupNum;
  FILE *ofp1, *ofp2;

  int maxval = 0, temp;
  
  initChannels();

  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++) {
    sprintf(fname,"ACG.%03d",frame);
    fprintf(stderr,"%s\n",fname);
    readCorrelogram(fname);
    for (chan = 0; chan < MAX_CHANNEL; chan++)
      energy[frame][chan] = (short) acg[0][chan];
    normalise();
    doCorrelation(frame);
    makeSummary(frame); 
    maxDelay=findSummaryMax();
    selectChannels(frame,maxDelay);
  }

  writeChannels();

  ofp1 = fopen("ENERGY", "w");
  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++) 
    for (chan = 0; chan < MAX_CHANNEL; chan++) {
      if (energy[frame][chan] > maxval)
	maxval = energy[frame][chan];
        fprintf(ofp1, "%d ", energy[frame][chan]);
    }
  printf("maximum energy is %d\n", maxval);

  /*
  numSegment=0;
  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++) 
    for (chan=0; chan<MAX_CHANNEL; chan++) {
      initGroup(&segment[numSegment]);
      searchGroup(&segment[numSegment],frame,chan);
      if (segment[numSegment].size>2) {
	writeGroup(numSegment);
	numSegment++;
	if (numSegment==MAX_SEGMENT) {
	  fprintf(stderr,"TOO MANY SEGMENTS! TERMINATING THE SEARCH");
	  exit(0);
	}
      }
    }

  fprintf(stderr,"Made the elements ... doing the search now.\n");
  */
  /* find the longest segment to seed the search */

  /*
  zeroCounts();
  seed=0;
  MaxSize=0;
  for (s=0; s<numSegment; s++) {
    if ((segment[s].last-segment[s].first) > MaxSize) {
      MaxSize=(segment[s].last-segment[s].first);
      seed=s;
    }
  }
  fprintf(stderr,"searching from segment %d\n",seed);
  segment[seed].group_id=1;
  add_pitch_info(seed);
  group_first=segment[seed].first;
  group_last=segment[seed].last;
  for (s1=0; s1<numSegment; s1++)
    for (s=0; s<numSegment; s++) {
      if ((segment[s].group_id==0) && overlap(s) && matchSegment(s)) {
	add_pitch_info(s);
	segment[s].group_id=1;
	if (segment[s].first<group_first) 
	  group_first=segment[s].first;
	if (segment[s].last>group_last) 
	  group_last=segment[s].last;
      }
    }
  */
  /* now find the next group */
  /*
  zeroCounts();
  seed=0;
  MaxSize=0;
  for (s=0; s<numSegment; s++) {
    if ((segment[s].group_id==0) && ((segment[s].last-segment[s].first) > MaxSize)) {
      MaxSize=segment[s].last-segment[s].first;
      seed=s;
    }
  }
  fprintf(stderr,"searching from segment %d\n",seed);
  segment[seed].group_id=2;
  add_pitch_info(seed);
  group_first=segment[seed].first;
  group_last=segment[seed].last;
  for (s1=0; s1<numSegment; s1++)
    for (s=0; s<numSegment; s++) {
      if ((segment[s].group_id==0) && overlap(s) && matchSegment(s)) {
	add_pitch_info(s);
	segment[s].group_id=2;
	if (segment[s].first<group_first) 
	  group_first=segment[s].first;
	if (segment[s].last>group_last) 
	group_last=segment[s].last;
      }
    }
  */
  /* make the groups */
  /*
  clearSourceMask();
  for (s=0; s<numSegment; s++) {
    if (segment[s].group_id==1)
      addToMask(segment[s].head);
  }
  writeSourceMaskToFile("ONE.BROWN");

  clearSourceMask();
  for (s=0; s<numSegment; s++) {
    if (segment[s].group_id==2)
      addToMask(segment[s].head);
  }
  writeSourceMaskToFile("TWO.BROWN");
  */
}
