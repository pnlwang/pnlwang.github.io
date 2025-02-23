/***************************************************************************/
/* Part of speech segregation program for the 1999 IEEE Trans. Neural Net. */
/* paper by D.L. Wang and G.J. Brown.                                      */
/* It takes the results of segmentation produced by "trim.c",  and         */
/* performs grouping. The results of grouping are output to two files:     */
/* STREAM1 and STREAM2                                                   */
/***************************************************************************/

#include <stdio.h>
#include <math.h>

#define MAX_CHANNEL 127
#define MAX_DELAY 200
#define MIN_DELAY 32
#define MAX_FRAME 300
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

int numChannels, numFrames;

short channels[MAX_FRAME][MAX_CHANNEL];
short mask[MAX_FRAME][MAX_CHANNEL];
group segment[MAX_SEGMENT];
short groupsize[MAX_SEGMENT];
short sourcemask[MAX_FRAME][MAX_CHANNEL];

void clearSourceMask()
{
  int frame, chan;

  for (frame=0; frame<MAX_FRAME; frame++)
    for (chan=0; chan<MAX_CHANNEL; chan++)
      sourcemask[frame][chan]=0;
}

void initGroup(group *g, int numFrames)
{
  int frame;

  for (frame=0; frame<numFrames; frame++) {
    g->oneCount[frame]=0;
    g->twoCount[frame]=0;
  }
  g->first=numFrames;
  g->last=0;
  g->size=0;
  g->head=NULL;
  g->group_id=0;
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

void writeMask()
{
  int chan, frame;

  for (frame=MIN_FRAME; frame<MAX_FRAME; frame++) {
    for (chan=0; chan<MAX_CHANNEL-1; chan++)
      printf("%d ",mask[frame][chan]);
    printf("\n");
  }
}

void writeSourceMaskToFile(char fname[],int numFrames, int numChannels)
{
  int chan, frame;
  FILE *ofp;

  ofp=fopen(fname,"w");
  if (ofp==NULL) {
    fprintf(stderr,"cannot open file %s\n",fname);
    exit(0);
  }
  //fprintf(ofp,"%d %d\n",numFrames,numChannels);
  for (frame=0; frame<numFrames; frame++) {
    for (chan=0; chan<numChannels; chan++)
      fprintf(ofp,"%d ",sourcemask[frame][chan]);
    fprintf(ofp,"\n");
  }
  fclose(ofp);
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

void writeGroup(int n, int numFrames)
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
  for (frame=0; frame<numFrames; frame++)
    fprintf(ofp,"%d ",segment[n].oneCount[frame]);
  fprintf(ofp,"\n");
  fprintf(ofp,"Count for group 2:\n");
  for (frame=0; frame<numFrames; frame++)
    fprintf(ofp,"%d ",segment[n].twoCount[frame]);
  fprintf(ofp,"\n");
  traverseToFile(segment[n].head,ofp);
  fclose(ofp);
}

void searchGroup(group *g, int frame, int chan, int numFrames, int numChannels)
{
  if (channels[frame][chan]!=0) {
    addPoint(g,frame,chan);
    if (channels[frame][chan]==1) g->oneCount[frame]++;
    if (channels[frame][chan]==2) g->twoCount[frame]++;
    channels[frame][chan]=0;
    if (chan<numChannels-1) 
      searchGroup(g,frame,chan+1,numFrames,numChannels);
    if (chan>0)
      searchGroup(g,frame,chan-1,numFrames,numChannels);
    if (frame>0)
      searchGroup(g,frame-1,chan,numFrames,numChannels);
    if (frame<numFrames-1)
      searchGroup(g,frame+1,chan,numFrames,numChannels);
  }
}

int matchPitch(int n1, int n2, int frame) 
/* matches the periodicity information in segments "n1" and "n2"
at time frame "frame" */
{
  int test1, test2;

  test1=(segment[n1].oneCount[frame] > segment[n1].twoCount[frame]);
  test2=(segment[n2].oneCount[frame] > segment[n2].twoCount[frame]);
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

int matchSegment(int n1, int n2)
{
  int st, fn, cnt, frame;
 
  cnt=0;
  st=max(segment[n1].first,segment[n2].first);
  fn=min(segment[n1].last,segment[n2].last);
  for (frame=st; frame<=fn; frame++)
    cnt+=matchPitch(n1,n2,frame);
  //fprintf(stderr,"SEED=%d TARGET=%d ST=%d FN=%d  CNT=%d\n",n2,n1,st,fn,cnt);
  if (((float)cnt/((float)fn-(float)st+1.0))>0.5)
    return TRUE;
  else
    return FALSE;
}

int overlap(int n1, int n2)
/* do segments n1 and n2 overlap? */
{
  int s1, s2, f1, f2;

  s1=segment[n1].first;
  f1=segment[n1].last;
  s2=segment[n2].first;
  f2=segment[n2].last;

  return ((s1>=s2)&&(f1<=f2)) | ((s2>=s1)&&(f2<=f1)) | ((s2>s1)&&(f1<f2)&&(f1>s2)) | ((s1>s2)&&(f2<f1)&&(f2>s1));
}

main(int argc, char **argv)
{
	int opt; 
	extern char *optarg;
	int ok = TRUE;
	char path[100];
	char fname[150];
	int frame, chan, dummy;
	short pg[MAX_DELAY];
	int numSegment;
	group newSegment;
	int seed, s, MaxSize, groupNum;
	FILE *ofp1, *ofp2, *ifp;

	while ((opt=getopt(argc,argv,"n:")) != EOF)
		switch(opt) {
	   case 'n': numFrames=atoi(optarg); break;
		case '?': ok = FALSE;
		}

	if (ok==FALSE) exit(0);
  
	numChannels=127;

	for (frame=0; frame<numFrames; frame++)
		for (chan=0; chan<numChannels; chan++) {
			scanf("%d",&dummy);
			channels[frame][chan]=(short)dummy;
			}
	fprintf(stderr,"read %d frames and %d channels...\n",numFrames,numChannels);

	numSegment=0;
	for (frame=0; frame<numFrames; frame++) 
		for (chan=0; chan<numChannels; chan++) {
			initGroup(&segment[numSegment],numFrames);
      	searchGroup(&segment[numSegment],frame,chan,numFrames,numChannels);
      	if ((segment[numSegment].last-segment[numSegment].first)>=2) {
				numSegment++;
				if (numSegment==MAX_SEGMENT) {
					fprintf(stderr,"TOO MANY SEGMENTS! TERMINATING THE SEARCH");
					exit(0);
					}
				}
    		}

	fprintf(stderr,"Made the elements ... doing the search now.\n");

	/* find the biggest segment to seed the search */

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
	for (s=0; s<numSegment; s++) {
		if ((segment[s].group_id==0) && overlap(s,seed) && matchSegment(s,seed)) {
			segment[s].group_id=1;
			}
		}


	/* now find the next group */
/*
	seed=0;
	MaxSize=0;
	for (s=0; s<numSegment; s++) {
		if ((segment[s].group_id==0) && ((segment[s].last-segment[s].first) > MaxSize)) {
			MaxSize=(segment[s].last-segment[s].first);
			seed=s;
			}
		}
	fprintf(stderr,"searching from segment %d\n",seed);
	segment[seed].group_id=2;
	for (s=0; s<numSegment; s++) {
		if ((segment[s].group_id==0) && overlap(s,seed) && matchSegment(s,seed)) {
			segment[s].group_id=2;
			}
		}
*/

	/* make the groups */
 
	clearSourceMask();
	for (s=0; s<numSegment; s++) {
		if (segment[s].group_id==1)
			addToMask(segment[s].head);
		}
	sprintf(fname,"STREAM1");
	writeSourceMaskToFile(fname,numFrames,numChannels);

	clearSourceMask();
	for (s=0; s<numSegment; s++) {
		if (segment[s].group_id==0)
			addToMask(segment[s].head);
		}
	sprintf(fname,"STREAM2");
	writeSourceMaskToFile(fname,numFrames,numChannels);

}






